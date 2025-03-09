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
		static SquareMatrix<Real, NDIM> const I = SquareMatrix<Real, NDIM>::identity();
		static auto const cHat = c;
		ColumnVector Beta;
		for (int n = 0; n < NDIM; n++) {
			Beta[n] = Beta0[n] + (F0[n] - F[n]) / (rho * c);
		}
		Real Ek = zero;
		Vector<Real, NDIM> dEk_dF;
		for (int n = 0; n < NDIM; n++) {
			dEk_dF[n] = -c * Beta[n];
		}
		Real const Eg = Eg0 + Er0 - Er;
		Real const dEg_dEr = -one;
		Real const eps = Eg - Ek;
		Vector<Real, NDIM> const deps_dF = -dEk_dF;
		Real const deps_dEr = -one;
		Real const F2 = vectorNorm(F);
		Real const absF = sqrt(F2);
		Real const iCv = (mu * mAMU) / ((gamma - one) * kB * rho);
		Real const T = eps * iCv;
		Real const dT_dEr = -iCv;
		Vector<Real, NDIM> const dT_dF = deps_dF * iCv;
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
			dgk_dF[n] = -kappa * (two * Beta[n] + three * aR * T3 * dT_dF[n]);
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

struct RadiationState: public Math::Vector<Real, 4> {
	static constexpr int NF = 4;
	using base_type = Math::Vector<Real, NF>;
	RadiationState() :
			F((Math::Vector<Real, NDIM>&) base_type::operator[](0)), E(base_type::operator[](NDIM)) {
	}
	Vector<Real, NDIM> &F;
	Real &E;
};

std::pair<Real, Real> eigenvalues(Real f, Real mu) {
	std::pair<Real, Real> rc;
	Real constexpr third = Real(1. / 3.), one = Real(1), two = Real(2), three = Real(3), four = Real(4);
	Real const f2 = f * f;
	Real const mu2 = mu * mu;
	Real const a = sqrt(four - three * f2);
	Real const b = a * a - a;
	Real const c = two - f2 - a;
	Real const d = two * (third * b + mu2 * c);
	Real const e = one / a;
	Real const g = mu * f;
	rc.first = e * (g - d);
	rc.second = e * (g + d);
	return rc;
}

struct FluxReturn {
	RadiationState flux;
	Real signalSpeed;
};

FluxReturn radiationFlux(RadiationState const &uL, RadiationState const &uR, int dim) {
	Real constexpr zero(0), sixth = Real(1. / 6.), half = Real(1. / 2.), one(1), two = Real(2), three = Real(3), four = Real(4);
	Real const tiny = Real::tiny();
	/*	static constexpr Constants<Real> pc;
	 mAMU = pc.m;
	 kB = pc.kB;
	 aR = pc.aR;
	 c = pc.c;
	 */
	Real const c = one;
	Real const cinv = one / c;
	Real const c2inv = cinv * cinv;
	FluxReturn rc;
	auto &flux = rc.flux;
	auto &signalSpeed = rc.signalSpeed;
	Math::Vector<Real, NDIM> pL;
	Math::Vector<Real, NDIM> pR;
	Math::Vector<Real, NDIM> FrL;
	Math::Vector<Real, NDIM> FrR;
	for (int k = 0; k < NDIM; k++) {
		FrL[k] = uL[k] * c2inv;
		FrR[k] = uR[k] * c2inv;
	}
	Real const ErL = uL[NDIM];
	Real const ErR = uR[NDIM];
	Real const F2L = Math::vectorNorm(FrL);
	Real const F2R = Math::vectorNorm(FrR);
	Real const FmagL = sqrt(F2L + tiny);
	Real const FmagR = sqrt(F2R + tiny);
	Real const fL = FmagL / ErL;
	Real const fR = FmagR / ErR;
	Real const FinvL = one / FmagL;
	Real const FinvR = one / FmagR;
	auto const nL = FrL * FinvL;
	auto const nR = FrR * FinvR;
	Real const muL = nL[dim];
	Real const muR = nR[dim];
	auto const lambdaL = eigenvalues(fL, muL);
	auto const lambdaR = eigenvalues(fR, muR);
	Real const sL = min(lambdaR.first, lambdaL.first);
	Real const sR = max(lambdaR.second, lambdaL.second);
	signalSpeed = max(abs(sR), abs(sL));
	Real const XiL = three - two * sqrt(four - fL * three * fL);
	Real const XiR = three - two * sqrt(four - fR * three * fR);
	for (int k = 0; k < NDIM; k++) {
		pL[k] = half * ErL * (one - XiL) * nL[k] * nL[dim];
		pR[k] = half * ErR * (one - XiR) * nR[k] * nR[dim];
	}
	pL[dim] += sixth * ErL * (XiL + one);
	pR[dim] += sixth * ErR * (XiR + one);
	if (sL >= zero) {
		for (int dim = 0; dim < NDIM; dim++) {
			flux[dim] = pL[dim];
		}
		flux[NDIM] = FrL[dim];
	} else if (/*sL < zero &&*/zero < sR) {
		Real const aInv = one / (sR - sL);
		Real const aL = sR * aInv;
		Real const aR = sL * aInv;
		Real const aRL = sR * sL * aInv;
		for (int dim = 0; dim < NDIM; dim++) {
			flux[dim] = aL * pL[dim] - aR * pR[dim] + aRL * (FrR[dim] - FrL[dim]);
		}
		flux[NDIM] = aL * FrL[dim] - aR * FrR[dim] + aRL * (ErR - ErL);
	} else /*if (sR >= zero)*/{
		for (int dim = 0; dim < NDIM; dim++) {
			flux[dim] = pR[dim];
		}
		flux[NDIM] = FrR[dim];
	}
	return rc;
}

struct RadiationGrid {
	static Real constexpr zero = Real(0), one = Real(1);
	static constexpr int NF = 4;
	static constexpr int BW = 1;
	using flux_array_type = std::array<std::vector<std::valarray<Real>>, NDIM>;
	RadiationGrid(int N_ = 64) {
		N = N_ + 2 * BW;
		N3 = N * N * N;
		dN[0] = 1;
		for (int dim = 0; dim < NDIM - 1; dim++) {
			dN[dim + 1] = dN[dim] * N;
		}
		bottomCorner = 0;
		for (int dim = 0; dim < NDIM; dim++) {
			bottomCorner += dN[dim];
		}
		topCorner = N3 - bottomCorner;
		for (int f = 0; f < NF; f++) {
			Ur[f] = std::valarray<Real>(N3);
			Ug[f] = std::valarray<Real>(N3);
		}
		rho = std::valarray<Real>(N3);
		dx = one / Real(N - 2 * BW);
	}
	Real computeFlux() {
		Real maxSignalSpeed = zero;
		for (int dim = 0; dim < NDIM; dim++) {
			for (int n = bottomCorner; n < N3; n++) {
				RadiationState uR, uL;
				F[dim] = std::vector<std::valarray<Real>>(NF, std::valarray<Real>(N3));
				for (int k = 0; k <= NDIM; k++) {
					uR[k] = Ur[k][n];
					uL[k] = Ur[k][n - dN[dim]];
				}
				auto const fluxReturn = radiationFlux(uL, uR, dim);
				for (int k = 0; k <= NDIM; k++) {
					F[dim][k][n] = fluxReturn.flux[k];
				}
				maxSignalSpeed = max(maxSignalSpeed, fluxReturn.signalSpeed);
			}
		}
		return maxSignalSpeed;
	}
	void implicitSolve(Real dt) {
		Real constexpr gamma = Real(5.0 / 3.0);
		Real const c = one;
		Real const cinv = one / c;
		ColumnVector Fr;
		ColumnVector Mg;
		for (int n = bottomCorner; n < topCorner; n++) {
			Real Er = Ur[NDIM][n];
			Real Eg = Ug[NDIM][n];
			for (int dim = 0; dim < NDIM; dim++) {
				Fr[dim] = Ur[dim][n] * cinv;
				Mg[dim] = Ug[dim][n];
			}
			Real const mu = Real(1);
			Real const kappa = Real(1e6);
			Real const Chi = Real(1e6);
			solveImplicitRadiation(Er, Fr, Eg, Mg, rho[n], mu, kappa, Chi, gamma, dt);
			Ur[NDIM][n] = Er;
			Ug[NDIM][n] = Eg;
			for (int dim = 0; dim < NDIM; dim++) {
				Ur[dim][n] = c * Fr[dim];
				Ug[dim][n] = Mg[dim];
			}
		}

	}
	void explicitSolve(Real beta, Real dt) {
		const Real alpha = one - beta;
		const Real lambda = dt / dx;
		for (int f = 0; f < NF; f++) {
			for (int n = bottomCorner; n < topCorner; n++) {
				Real dU = zero;
				for (int dim = 0; dim < NDIM; dim++) {
					dU -= lambda * (F[dim][f][n + dN[dim]] - F[dim][f][n]);
				}
				auto const U1 = Ur[f][n] + dU;
				Ur[f][n] = alpha * Ur0[f][n] + beta * U1;
			}
		}
	}
	Real prepareStep() {
		Ur0 = Ur;
		return dx / computeFlux();
	}
	void doStep(Real dt) {
		explicitSolve(one, dt);
		implicitSolve(dt);
	}
	void finishStep() {
		std::array<int, NDIM> n;
		auto &nx = n[0];
		auto &ny = n[1];
		auto &nz = n[2];
		for (nx = 0; nx < BW; nx++) {
			for (ny = BW; ny < N - BW; ny++) {
				for (nz = BW; nz < N - BW; nz++) {
					for (int dim = 0; dim < NDIM; dim++) {
						std::swap(n[0], n[dim]);
						int iLd = 0;
						for (int k = 0; k < NDIM; k++) {
							iLd += dN[k] * n[k];
						}
						int iLs = iLd, iHs = iLd, iHd = iLd;
						iLs += (BW - dN[dim]) * n[dim];
						iHs += (N - dN[dim] - BW) * n[dim];
						for (int f = 0; f < NF; f++) {
							Ur[f][iLd] = Ur[f][iLs];
							Ur[f][iHd] = Ur[f][iHs];
						}
						Ur[dim][iLd] = min(zero, Ur[dim][iLd]);
						Ur[dim][iHd] = max(zero, Ur[dim][iHd]);
						std::swap(n[0], n[dim]);
					}
				}
			}
		}
	}
private:
	std::array<std::valarray<Real>, NF> Ur0;
	std::array<std::valarray<Real>, NF> Ur;
	std::array<std::valarray<Real>, NF> Ug;
	std::valarray<Real> rho;
	flux_array_type F;
	int N;
	int N3;
	int bottomCorner;
	int topCorner;
	Real dx;
	Math::Vector<int, NDIM> dN;
};

RadiationGrid driver(RadiationGrid grid, Real tMax) {
	static constexpr Real zero = Real(0), cfl = Real(0.4);
	Real t = zero;
	while (t < tMax) {
		Real const dt = cfl * grid.prepareStep();
		grid.doStep(dt);
		grid.finishStep();
		t += dt;
	}
	return grid;
}
void testRadiation() {
	RadiationGrid test;
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
