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
Vector<T, N> findRootNelderMead(std::function<Vector<T, N>(Vector<T, N> const&)> const &testFunction,
		Vector<T, N> const &xGuess) {
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

Vector<Real, NF> testImplicitRadiation(Real dt, ColumnVector F, Real Erad, ColumnVector Beta, Real Egas,
		ColumnVector F0, Real Erad0, ColumnVector Beta0, Real Egas0, Real rho, Real mu, Real kappa, Real Chi,
		Real gamma) {
	static constexpr Constants<Real> pc;
	static constexpr Real third(1.0 / 3.0), half(0.5), one(1), two(2), three(3), four(4);
	static auto const diffusionTensor = third * identityMatrix<Real, NDIM>();
	Vector<Real, NDIM + 1> resid;
	Real const f = vectorMagnitude(F) / Erad;
	ColumnVector const unitF = F / (f * Erad);
	RowVector const unitFT = matrixTranspose(unitF);
	RowVector const BetaT = matrixTranspose(Beta);
	Real const Xi = three - two * sqrt(four - f * three * f);
	auto const streamingTensor = unitF * unitFT;
	auto const dEdd = (Xi + one) * diffusionTensor;
	auto const sEdd = (Xi - one) * streamingTensor;
	auto const P = half * Erad * (dEdd + sEdd);
	Real const T4 = integerPower<Real>((mu * pc.m * Egas) / ((gamma - one) * pc.kB * rho), 4);
	Real const gk = kappa * (Erad - pc.aR * T4 - two * BetaT * F);
	ColumnVector const Gk = gk * Beta;
	ColumnVector const Gx = Chi * (F - Erad * Beta - P * Beta);
	Real const gx = BetaT * Gx;
	Real const r = Erad - Erad0 + dt * (gk + gx);
	ColumnVector const R = F - F0 + dt * (Gk + Gx);
	std::copy(R.begin(), R.end(), resid.begin());
	resid[NDIM] = r;
	return resid;
}

