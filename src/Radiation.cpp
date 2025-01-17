/*
 * Radiation.cpp
 *
 *  Created on: Jan 15, 2025
 *      Author: dmarce1
 */

#include "Constants.hpp"
#include "Numbers.hpp"
#include "Real.hpp"
#include "Vector.hpp"

#include <algorithm>
#include <functional>

using namespace Math;

template<typename T, int N>
Vector<T, N> findRootNelderMead(std::function<Vector<T, N>(Vector<T, N> const&)> const&, Vector<T, N> const&);

template<typename T, int N>
Vector<T, N> findRootNelderMead(std::function<Vector<T, N>(Vector<T, N> const&)> const &testFunction, Vector<T, N> const &xGuess) {
	using simplex_type = std::pair<T, Vector<T, N>>;
	static constexpr T zero = T(0), one = T(1);
	static constexpr T alpha = T(1), rho = T(0.5), sigma = T(0.5), gamma = T(2), tolerance(1e-6);
	Vector<simplex_type, N + 1> xSimplex;
	T const iN = one / T(N);
	T const dX = T(1.0e-3) * vectorMagnitude(xGuess);
	for (int k = 0; k < N; k++) {
		xSimplex[k].second = xGuess;
		xSimplex[k].second[k] += dX;
		xSimplex[k].first = testFunction(xSimplex[k].second);
	}
	xSimplex[N].second = xGuess;
	xSimplex[N].first = testFunction(xSimplex[N].second);
	while (true) {
		bool shrink;
		simplex_type reflectionPoint;
		simplex_type contractionPoint;
		Real error = zero;
		std::sort(xSimplex.begin(), xSimplex.end(), [](simplex_type const &a, simplex_type const &b) {
			return abs(a.first) < abs(b.first);
		});
		for (int k = 0; k <= N; k++) {
			error += nSquared(xSimplex[k].first);
		}
		error = sqrt(error * iN);
		if (error <= tolerance) {
			return xSimplex[0].second;
		}
		Vector<T, N> xCentroid = zeroVector<T, N>();
		for (int k = 0; k < N; k++) {
			xCentroid += xSimplex[k].second;
		}
		xCentroid *= iN;
		reflectionPoint.second = xCentroid + alpha * (xCentroid - xSimplex[N].second);
		reflectionPoint.first = testFunction(reflectionPoint.second);
		shrink = false;
		if (abs(reflectionPoint.first) <= abs(xSimplex[0].first)) {
			simplex_type expansionPoint;
			expansionPoint.second = (one - gamma) * xCentroid + gamma * reflectionPoint.second;
			expansionPoint.first = testFunction(expansionPoint.second);
			if (abs(expansionPoint.first) < abs(reflectionPoint.first)) {
				xSimplex[N] = expansionPoint;
			} else {
				xSimplex[N] = reflectionPoint;
			}
		} else if (abs(reflectionPoint.first) <= abs(xSimplex[N - 1].first)) {
			xSimplex[N] = reflectionPoint;
		} else if (abs(reflectionPoint.first) <= abs(xSimplex[N].first)) {
			contractionPoint.second = (one - rho) * xCentroid + rho * reflectionPoint.second;
			contractionPoint.first = testFunction(contractionPoint.second);
			if (abs(contractionPoint.first) <= abs(reflectionPoint.first)) {
				xSimplex[N] = contractionPoint;
			} else {
				shrink = true;
			}
		} else {
			contractionPoint.second = (one - rho) * xCentroid + rho * xSimplex[N].second;
			contractionPoint.first = testFunction(contractionPoint.second);
			if (abs(contractionPoint.first) <= abs(xSimplex[N].first)) {
				xSimplex[N] = contractionPoint;
			} else {
				shrink = true;
			}
		}
		if (shrink) {
			for (int k = 1; k <= N; k++) {
				xSimplex[k].second = (one - sigma) * xSimplex[0].second + sigma * xSimplex[k].second;
				xSimplex[k].first = testFunction(xSimplex[k].second);
			}
		}
	}
}

static constexpr int XDIM = 0;
static constexpr int YDIM = 1;
static constexpr int ZDIM = 2;
static constexpr int NF = 4;

using ColumnVector = Vector<Real, NDIM>;
using RowVector = Matrix<Real, 1, NDIM>;

auto testImplicitRadiation(Real dt, ColumnVector F, Real Er, ColumnVector F0, Real Er0, ColumnVector Beta0, Real Eg0, Real rho, Real mu, Real kappa, Real Chi,
		Real gamma) {
	static constexpr Constants<Real> pc;
	static constexpr Real zero(0), third(1.0 / 3.0), half(0.5), one(1), two(2), three(3), four(4), twelve(12);
	static SquareMatrix<Real, NDIM> const I = identityMatrix<Real, NDIM>();
	auto const cHat = pc.c;
	ColumnVector const Beta = Beta0 + (F0 - F) / (rho * pc.c);
	Real const Eg = Eg0 + Er0 - Er;
	Real const F2 = vectorNorm(F);
	Real const absF = sqrt(F2);
	Real const f = absF / Er;
	ColumnVector const N = F / absF;
	SquareMatrix<Real, NDIM> NN = N * matrixTranspose(N);
	Real const Xi = three - two * sqrt(four - f * three * f);
	auto const D = half * ((Xi + one) * third * I + (one - Xi) * NN);
	RowVector const BetaT = matrixTranspose(Beta);
	auto const P = Er * D;
	Real const T = Eg * (mu * pc.m) / ((gamma - one) * pc.kB * rho);
	Real const T2 = nSquared(T);
	Real const T4 = nSquared(T2);
	Real const gk = kappa * (Er - pc.aR * T4 - two * BetaT * F);
	ColumnVector const Gx = Chi * (F - (Er * I + P) * Beta);
	ColumnVector const Gk = gk * Beta;
	Real const gx = BetaT * Gx;
	Real const hr = Er - Er0 + dt * cHat * (gk + gx);
	ColumnVector const Hr = F - F0 + dt * (Gk + Gx);
	Real const T3 = T * T2;
	Real const dXi_df = twelve * f / (three - Xi);
	Real const df_dEr = -f / Er;
	Real const dXi_dEr = dXi_df * df_dEr;
	SquareMatrix<Real, NDIM> const dD_dEr = half * (third * dXi_dEr * I - dXi_dEr * NN);
	SquareMatrix<Real, NDIM> const dP_dEr = dD_dEr * Er + D;
	Real const dT_dEr = -T / Eg;
	Real const dgk_dEr = kappa * (one - four * pc.aR * T3 * dT_dEr);
	ColumnVector const dGx_dEr = -Chi * (I + dP_dEr) * Beta;
	ColumnVector const dGk_dEr = dgk_dEr * Beta;
	Real const dgx_dEr = BetaT * dGx_dEr;
	Real const dhr_dEr = one + dt * cHat * (dgk_dEr + dgx_dEr);
	ColumnVector const dHr_dEr = dt * cHat * (dGk_dEr + dGx_dEr);
	ColumnVector const df_dF = N / Er;
	SquareMatrix<Real, NDIM> const dN_dF = (I - NN) / absF;
	Vector<SquareMatrix<Real, NDIM>, NDIM> dNN_dF;
	for (int k = 0; k < NDIM; k++) {
		for (int n = 0; n < NDIM; n++) {
			for (int m = 0; m < NDIM; m++) {
				dNN_dF[k][n, m] = N[k] * dN_dF[n, m] + dN_dF[k, n] * N[m];
			}
		}
	}
	ColumnVector const dXi_dF = dXi_df * df_dF;
	Vector<SquareMatrix<Real, NDIM>, NDIM> dD_dF;
	for (int k = 0; k < NDIM; k++) {
		dD_dF[k] = half * (dXi_dF[k] * (third * I - NN) + (one - Xi) * dNN_dF[k]);
	}
	Vector<SquareMatrix<Real, NDIM>, NDIM> const dP_dF = Er * dD_dF;
	ColumnVector const dgk_dF = -two * kappa * (Beta0 + F0 - two * F / (rho * pc.c));
	//Chi * (F - (Er * I + P) * Beta);
	SquareMatrix<Real, NDIM> const dGk_dF = dgk_dF * BetaT - gk * I / (rho * pc.c);
	//Chi * (F - (Er * I + P) * Beta)
	SquareMatrix<Real, NDIM> dGx_dF;
	for (int k = 0; k < NDIM; k++) {
		for (int n = 0; n < NDIM; n++) {
			dGx_dF[n, k] = zero;
			for (int m = 0; m < NDIM; m++) {
				dGx_dF[n, k] += Chi * (I[n, k] - dP_dF[k][n, m] * Beta[m] + (Er * I[n, m] + P[n, m]) * I[k, m] / (rho * pc.c));
			}
		}
	}
	ColumnVector const dgx_dF = dGx_dF * Beta - (I / (rho * pc.c)) * Gx;
	ColumnVector const dhr_dF = dt * cHat * (dgk_dF + dgx_dF);
	SquareMatrix<Real, NDIM> const dHr_dF = I + dt * cHat * (dGk_dF + dGx_dF);
	std::pair<Vector<Real, NF>, SquareMatrix<Real, NF>> rc;
	Vector<Real, NDIM + 1> &F4 = rc.first;
	SquareMatrix<Real, NDIM + 1> &dF4 = rc.second;
	for (int k = 0; k < NDIM; k++) {
		F4[k] = Hr[k];
		for (int n = 0; n < NDIM; n++) {
			dF4[k, n] = dHr_dF[k, n];
		}
		dF4[NDIM, k] = dhr_dF[k];
		dF4[k, NDIM] = dHr_dEr[k];
	}
	F4[NDIM] = hr;
	dF4[NDIM, NDIM] = dhr_dEr;
	return rc;
}

