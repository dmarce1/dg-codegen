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
using Tensor = SquareMatrix<Real, NDIM>;

struct testImplicitRadiation {
	testImplicitRadiation(Real dt_, Real Er0_, ColumnVector F0_, Real Eg0_, ColumnVector Beta0_, Real rho_, Real mu_, Real kappa_, Real Chi_, Real gamma_) {
		/*	static constexpr Constants<Real> pc;
		 mAMU = pc.m;
		 kB = pc.kB;
		 aR = pc.aR;
		 c = pc.c;
		 */
		mAMU = kB = aR = c = Real(1);
		dt = dt_;
		Er0 = Er0_;
		Eg0 = Eg0_;
		rho = rho_;
		mu = mu_;
		kappa = kappa_;
		Chi = Chi_;
		gamma = gamma_;
		F0 = F0_;
		Beta0 = Beta0_;
	}
	auto operator()(Vector<Real, NF> U) const {
		auto Er = U[NDIM];
		ColumnVector F;
		for (int n = 0; n < NDIM; n++) {
			F[n] = U[n];
		}
		static constexpr Real zero(0), third(1.0 / 3.0), half(0.5), one(1), two(2), three(3), four(4), twelve(12);
		static SquareMatrix<Real, NDIM> const I = identityMatrix<Real, NDIM>();
		static auto const cHat = c;
		ColumnVector Beta;
		for (int n = 0; n < NDIM; n++) {
			Beta[n] = Beta0[n] + (F0[n] - F[n]) / (rho * c);
		}
		Real const Eg = Eg0 + Er0 - Er;
		Real const F2 = vectorNorm(F);
		Real const absF = sqrt(F2);
		Real const T = Eg * (mu * mAMU) / ((gamma - one) * kB * rho);
		Real const T2 = nSquared(T);
		Real const T4 = nSquared(T2);
		Tensor dBeta_dF;
		for (int n = 0; n < NDIM; n++) {
			for (int m = 0; m < NDIM; m++) {
				dBeta_dF[n, m] = -I[n, m] / (rho * c);
			}
		}
		ColumnVector N;
		for (int n = 0; n < NDIM; n++) {
			N[n] = F[n] / (absF + Real::tiny());
		}
		Real const T3 = T * T2;
		Real const dT_dEr = -T / Eg;
		Tensor dN_dF;
		for (int n = 0; n < NDIM; n++) {
			for (int m = 0; m < NDIM; m++) {
				dN_dF[n, m] = (I[n, m] - N[n] * N[m]) / (absF + Real::tiny());
			}
		}
		ColumnVector df_dF;
		Real const f = absF / Er;
		Real const df_dEr = -f / Er;
		for (int n = 0; n < NDIM; n++) {
			df_dF[n] = N[n] / Er;
		}
		Real const Xi = three - two * sqrt(four - f * three * f);
		Real const dXi_df = twelve * f / (three - Xi);
		Real const dXi_dEr = dXi_df * df_dEr;
		ColumnVector dXi_dF;
		for (int n = 0; n < NDIM; n++) {
			dXi_dF[n] = dXi_df * df_dF[n];
		}
		Tensor D;
		for (int n = 0; n < NDIM; n++) {
			for (int m = 0; m < NDIM; m++) {
				D[n, m] = half * ((Xi + one) * third * I[n, m] + (one - Xi) * N[n] * N[m]);
			}
		}
		Tensor dD_dEr;
		for (int n = 0; n < NDIM; n++) {
			for (int m = 0; m < NDIM; m++) {
				dD_dEr[n, m] = half * dXi_dEr * (third * I[n, m] - N[n] * N[m]);
			}
		}
		Vector<Tensor, NDIM> dD_dF;
		for (int l = 0; l < NDIM; l++) {
			for (int n = 0; n < NDIM; n++) {
				for (int m = 0; m < NDIM; m++) {
					dD_dF[l][n, m] = half * (dXi_dF[l] * (third * I[n, m] - N[n] * N[m]) + (one - Xi) * (dN_dF[n, l] * N[m] + N[n] * dN_dF[m, l]));
				}
			}
		}
		Tensor P;
		for (int n = 0; n < NDIM; n++) {
			for (int m = 0; m < NDIM; m++) {
				P[n, m] = Er * D[n, m];
			}
		}
		Tensor dP_dEr;
		for (int n = 0; n < NDIM; n++) {
			for (int m = 0; m < NDIM; m++) {
				dP_dEr[n, m] = Er * dD_dEr[n, m] + D[n, m];
			}
		}
		Vector<Tensor, NDIM> dP_dF;
		for (int n = 0; n < NDIM; n++) {
			for (int m = 0; m < NDIM; m++) {
				for (int l = 0; l < NDIM; l++) {
					dP_dF[l][n, m] = Er * dD_dF[l][n, m];
				}
			}
		}
		Real gk = kappa * (Er - aR * T4);
		for (int k = 0; k < NDIM; k++) {
			gk -= two * kappa * Beta[k] * F[k];
		}
		Real const dgk_dEr = kappa * (one - four * aR * T3 * dT_dEr);
		ColumnVector dgk_dF;
		for (int n = 0; n < NDIM; n++) {
			dgk_dF[n] = -two * kappa * Beta[n];
			for (int k = 0; k < NDIM; k++) {
				dgk_dF[n] -= two * kappa * F[k] * dBeta_dF[k, n];
			}
		}
		ColumnVector Gx;
		for (int n = 0; n < NDIM; n++) {
			Gx[n] = Chi * (F[n] - Beta[n] * Er);
			for (int k = 0; k < NDIM; k++) {
				Gx[n] -= Chi * P[n, k] * Beta[k];
			}
		}
		ColumnVector dGx_dEr;
		for (int n = 0; n < NDIM; n++) {
			dGx_dEr[n] = -Chi * Beta[n];
			for (int k = 0; k < NDIM; k++) {
				dGx_dEr[n] -= Chi * dP_dEr[n, k] * Beta[k];
			}
		}
		Tensor dGx_dF;
		for (int n = 0; n < NDIM; n++) {
			for (int m = 0; m < NDIM; m++) {
				dGx_dF[n, m] = Chi * (I[n, m] - dBeta_dF[n, m] * Er);
				for (int k = 0; k < NDIM; k++) {
					dGx_dF[n, m] -= Chi * dP_dF[m][n, k] * Beta[k];
					dGx_dF[n, m] -= Chi * P[n, k] * dBeta_dF[k, m];
				}
			}
		}
		Real gx = zero;
		for (int k = 0; k < NDIM; k++) {
			gx += Gx[k] * Beta[k];
		}
		Real dgx_dEr = zero;
		for (int k = 0; k < NDIM; k++) {
			dgx_dEr += dGx_dEr[k] * Beta[k];
		}
		ColumnVector dgx_dF;
		for (int n = 0; n < NDIM; n++) {
			dgx_dF[n] = zero;
			for (int k = 0; k < NDIM; k++) {
				dgx_dF[n] = dGx_dF[k, n] * Beta[k] + dBeta_dF[k, n] * Gx[k];
			}
		}
		ColumnVector Gk;
		for (int n = 0; n < NDIM; n++) {
			Gk[n] = gk * Beta[n];
		}
		ColumnVector dGk_dEr;
		for (int n = 0; n < NDIM; n++) {
			dGk_dEr[n] = dgk_dEr * Beta[n];
		}
		Tensor dGk_dF;
		for (int n = 0; n < NDIM; n++) {
			for (int m = 0; m < NDIM; m++) {
				dGk_dF[n, m] = Beta[n] * dgk_dF[m] + gk * dBeta_dF[n, m];
			}
		}
		Real const hr = Er - Er0 + dt * cHat * (gk + gx);
		Real const dhr_dEr = one + dt * cHat * (dgk_dEr + dgx_dEr);
		ColumnVector dhr_dF;
		for (int n = 0; n < NDIM; n++) {
			dhr_dF[n] = dt * cHat * (dgk_dF[n] + dgx_dF[n]);
		}
		ColumnVector Hr;
		for (int n = 0; n < NDIM; n++) {
			Hr[n] = F[n] - F0[n] + dt * (Gk[n] + Gx[n]);
		}
		ColumnVector dHr_dEr;
		for (int n = 0; n < NDIM; n++) {
			dHr_dEr[n] = dt * cHat * (dGk_dEr[n] + dGx_dEr[n]);
		}
		Tensor dHr_dF;
		for (int n = 0; n < NDIM; n++) {
			for (int m = 0; m < NDIM; m++) {
				dHr_dF[n, m] = I[n, m] + dt * cHat * (dGk_dF[n, m] + dGx_dF[n, m]);
			}
		}

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
private:
	Real dt;
	Real Er0;
	Real Eg0;
	Real rho;
	Real mu;
	Real kappa;
	Real Chi;
	Real gamma;
	ColumnVector F0;
	ColumnVector Beta0;
	Real mAMU;
	Real kB;
	Real aR;
	Real c;
};

void solveImplicitRadiation(Real &Er, ColumnVector &F, Real &Eg, ColumnVector &Mg, Real rho, Real mu, Real kappa, Real Chi, Real gamma, Real dt) {
	/*	static constexpr Constants<Real> pc;
	 * 	Real c;
	 *
	 mAMU = pc.m;
	 kB = pc.kB;
	 aR = pc.aR;
	 c = pc.c;
	 */
	Real const c = Real(1);
	auto const Beta0 = Mg / (rho * c);
	auto const Eg0 = Eg;
	testImplicitRadiation test(dt, Er, F, Eg0, Beta0, rho, mu, kappa, Chi, gamma);
	Vector<Real, NF> x;
	for (int n = 0; n < NDIM; n++) {
		x[n] = F[n] / c;
	}
	x[NDIM] = Er;
	Real toler = Real(1e-9);
	Real err;
	do {
		auto const f_and_dfdx = test(x);
		auto const &f = f_and_dfdx.first;
		auto const &dfdx = f_and_dfdx.second;
		auto const inv_dfdx = matrixInverse(dfdx);
		auto const dx = -inv_dfdx * f;
		err = vectorMagnitude(dx);
		x += Real(0.5) * dx;
	} while (err > toler);
	for (int n = 0; n < NDIM; n++) {
		Mg[n] = rho * c * Beta0[n] + F[n] - x[n];
		F[n] = c * x[n];
	}
	Eg = Eg0 + Er - x[NDIM];
	Er = x[NDIM];
}

void testRadiation() {
	Real a(.354);
	ColumnVector F = Real(Real(0.0) * a);
	ColumnVector Mg = Real(-Real(0.0) * a);
	F[0] = Real(0.99) * a;
	Real Er = Real(a);
	Real Eg = Real(2 * a);
	Real rho = Real(1);
	Real mu = Real(1);
	Real kappa = Real(1e6);
	Real Chi = Real(1e6);
	Real gamma = Real(5.0 / 3.0);
	Real dt = Real(1e-9);
	printf("Er0 = %e Fx = %e Fy = %e Fz = %e \n", Er, F[0], F[1], F[2]);
	printf("Eg0 = %e Mg = %e Mg = %e Mg = %e \n", Eg, Mg[0], Mg[1], Mg[2]);
	solveImplicitRadiation(Er, F, Eg, Mg, rho, mu, kappa, Chi, gamma, dt);
	printf("Er = %e Fx = %e Fy = %e Fz = %e \n", Er, F[0], F[1], F[2]);
	printf("Eg = %e Mg = %e Mg = %e Mg = %e \n", Eg, Mg[0], Mg[1], Mg[2]);

}
