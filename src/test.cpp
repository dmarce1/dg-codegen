/*
 * test.cpp
 *
 *  Created on: Dec 18, 2024
 *      Author: dmarce1
 */

#include "Real.hpp"
#include "HLLC.hpp"
#include "Basis.hpp"
#include <vector>

constexpr int N = 1024;
constexpr int BW = 2;
constexpr int P = 3;

static const Real dx = Real(1) / Real(N - 1);

using state_type = HydroState<Real, 1>;

void compute2() {
	HLLCSolver<Real, 1> riemann;
	std::vector<state_type> U(N);
	std::vector<state_type> U0(N);
	std::vector<state_type> UL(N);
	std::vector<state_type> UR(N);
	std::vector<state_type> F(N);
	std::vector<state_type> S1(N);
	std::vector<state_type> S2(N);

	Real tMax = Real(.1);
	for (int n = 0; n < N; n++) {
		state_type u;
		Real vel1 = Real(0);
		Real vel2 = Real(0);
		Real rho1 = Real(1);
		Real rho2 = Real(.1);
		Real ene1 = Real(2.5);
		Real ene2 = Real(.125);
		if (n < N / 2) {
			u[D_i] = rho1;
			u[S_i] = rho1 * vel1;
			u[E_i] = ene1 + rho1 * vel1 * vel1 * Real(0.5);
		} else {
			u[D_i] = rho2;
			u[S_i] = rho2 * vel2;
			u[E_i] = ene2 + rho2 * vel2 * vel2 * Real(0.5);
		}
		U[n] = std::move(u);
	}
	int n = 0;
	Real dt;
	const std::array<std::array<Real, 3>, 3> beta0 = { { { Real(1) }, { Real(1), Real(0.5) }, { Real(1), Real(0.25),
			Real(2.0 / 3.0) } } };
	const std::array beta = P == 1 ? beta0[0] : (P == 2 ? beta0[1] : beta0[2]);
	constexpr int rkOrder = std::min(P, 3);
	Real t = Real(0);
	while (t < tMax) {
		Real lambdaMax = Real(0);
		Real lambda;
		U0 = U;
		for (int rk = 0; rk < rkOrder; rk++) {
			auto const dN = N - 2 * BW;
			for (int n = 0; n < BW; n++) {
				U[n] = U[n + dN];
				U[N - n - 1] = U[N - n - 1 - dN];
			}
			for (int n = BW - 2; n < N - BW + 2; n++) {
				Vector<Real, Nfield> dW;
				auto const R = primToEigenvectors<Real, 1>(conToPrim<Real, 1>(U[n]));
				auto const invR = Math::matrixInverse(R);
				auto const dWm = toVector(invR * toColumnVector(U[n] - U[n - 1]));
				auto const dWp = toVector(invR * toColumnVector(U[n + 1] - U[n]));
				S1[n] = superBee(dWp, dWm);
			}
			for (int n = BW - 1; n < N - BW + 1; n++) {
				auto const sp = S1[n + 1];
				auto const s0 = S1[n];
				auto const sm = S1[n - 1];
				S2[n] = superBee(sp - s0, s0 - sm);
				for (int f = 0; f < Nfield; f++) {
					S2[n][f] = sign(S2[n][f]) * min(abs(S2[n][f]), abs(Real(2) * S1[n][f]));
				}
			}
			for (int n = BW - 1; n < N - BW + 1; n++) {
				Vector<Real, Nfield> dW;
				auto const R = primToEigenvectors<Real, 1>(conToPrim<Real, 1>(U[n]));
				S1[n] = toVector(R * toColumnVector(S1[n]));
				S2[n] = toVector(R * toColumnVector(S2[n]));
			}
			for (int n = BW - 1; n < N - BW + 1; n++) {
				Real const c2(0.5);
				Real const c6(1.0 / 6.0);
				UL[n + 1] = U[n] + c2 * S1[n] + c6 * S2[n];
				UR[n] = UL[n + 1] - S1[n];
			}
			for (int n = BW; n < N - BW + 1; n++) {
				auto const rc = riemann(UL[n], UR[n]);
				F[n] = rc.flux;
				lambdaMax = std::max(lambdaMax, rc.signalSpeed);
			}
			if (rk == 0) {
				dt = Real(0.4) * dx / lambdaMax;
				lambda = dt / dx;
				printf("%i %e %e %e\n", n, t, dt, lambdaMax);
			}
			for (int n = BW; n < N - BW; n++) {
				U[n] -= lambda * (F[n + 1] - F[n]);
			}
			Real const w1 = beta[rk];
			Real const w0 = Real(1) - w1;
			for (int n = BW; n < N - BW; n++) {
				U[n] = w1 * U[n] + w0 * U0[n];
			}
		}
		t += dt;
		n++;
	}
	std::string const fname = "test.dat";
	FILE *fp = fopen(fname.c_str(), "wt");
	for (int n = 0; n < N; n++) {
		fprintf(fp, "\t%i \t%e %e %e \n", n, U[n][0], U[n][1], U[n][2]);
	}
	fclose(fp);
}

void compute() {
	HLLCSolver<Real, 1> riemann;
	std::vector<Vector<state_type, P>> U(N);
	std::vector<Vector<state_type, P>> U0(N);
	std::vector<state_type> F(N);
	std::vector<Vector<state_type, P>> dUdt(N);
	std::vector<state_type> dUL(N);
	std::vector<state_type> dUR(N);
	constexpr int Q = P + 1;
	auto const quadRules = gaussLobattoRules<Q>();
	Real t = Real(0);
	Real tMax = Real(.1);
	for (int n = 0; n < N; n++) {
		Vector<state_type, P> u;
		Real vel1 = Real(1.5);
		Real vel2 = Real(-1.5);
		Real rho1 = Real(1);
		Real rho2 = Real(.1);
		Real ene1 = Real(2.5);
		Real ene2 = Real(.125);
		if (n < N / 2) {
			u[0][D_i] = rho1;
			u[0][S_i] = rho1 * vel1;
			u[0][E_i] = ene1 + rho1 * vel1 * vel1 * Real(0.5);
		} else {
			u[0][D_i] = rho2;
			u[0][S_i] = rho2 * vel2;
			u[0][E_i] = ene2 + rho2 * vel2 * vel2 * Real(0.5);
		}
		for (int p = 1; p < P; p++) {
			for (int f = 0; f < Nfield; f++) {
				u[p][f] = Real(0);
			}
		}
		U[n] = std::move(u);
	}
	int n = 0;
	Basis<P> const basis(dx);
	Real dt;
	Real const one = Real(1.0);
	Real const half = Real(0.5);
	Real const zero = Real(0.0);
	const std::array<std::array<Real, 3>, 3> beta0 = { { { Real(1) }, { Real(1), Real(0.5) }, { Real(1), Real(0.25),
			Real(2.0 / 3.0) } } };
	const std::array beta = ((P == 1) ? beta0[0] : ((P == 2) ? beta0[1] : beta0[2]));
	constexpr int rkOrder = std::min(P, 3);
	auto const an = basis.coeffs();
	while (t < tMax) {
		U0 = U;
		for (int rk = 0; rk < rkOrder; rk++) {
			for (int n = 0; n < BW; n++) {
				U[n] = U[n + N - 2 * BW];
				U[N - n - 1] = U[2 * BW - n - 1];
			}
			int nlimit = 0;
			int ntot = 0;
			for (int n = 1; n < N - 1; n++) {
				dUL[n] = state_type(zero);
				dUR[n] = state_type(zero);
				for (int l = 1; l < P; l++) {
					const auto ul = U[n][l];
					const auto al = an[l];
					dUL[n] += al * ul * basis.function(l, +half * dx);
					dUR[n] -= al * ul * basis.function(l, -half * dx);
				}
			}
			for (int n = 1; n < N - 1; n++) {
				auto const delP = U[n + 1][0] - U[n][0];
				auto const delM = U[n][0] - U[n - 1][0];
				for (int f = 0; f < Nfield; f++) {
					ntot++;
					auto const del = minmod(delP[f], delM[f]);
					Real const dL = minmod(dUL[n][f], del);
					Real const dR = minmod(dUR[n][f], del);
					if ((dL != dUL[n][f]) || (dR != dUR[n][f])) {
						nlimit++;
						for (int l = 1; l < std::min(P, 3); l++) {
							auto const cp = one / (an[l] * basis.function(l, +half * dx));
							auto const cn = one / (an[l] * basis.function(l, -half * dx));
							U[n][l][f] = half * (cp * dL - cn * dR);
						}
						for (int l = 3; l < P; l++) {
							U[n][l][f] = zero;
						}
						dUL[n][f] = zero;
						dUR[n][f] = zero;
						for (int l = 1; l < std::min(P, 3); l++) {
							const auto ul = U[n][l][f];
							const auto al = an[l];
							dUL[n][f] += al * ul * basis.function(l, +half * dx);
							dUR[n][f] -= al * ul * basis.function(l, -half * dx);
						}
					}
				}
			}
			Real lambdaMax = Real(0);
			for (int n = 2; n < N - 1; n++) {
				for (int f = 0; f < Nfield; f++) {
					auto const rc = riemann(U[n - 1][0] + dUL[n - 1], U[n][0] - dUR[n]);
					F[n] = rc.flux;
					lambdaMax = std::max(lambdaMax, rc.signalSpeed);
				}
			}
			if (rk == 0) {
				dt = Real(0.4) * dx / (Real(2 * P - 1) * lambdaMax);
				printf("%i %e %e %e %i %i\n", n, t, dt, lambdaMax, nlimit, ntot);
			}
			for (int n = 2; n < N - 2; n++) {
				Vector<state_type, Q> f;
				for (int q = 1; q < Q - 1; q++) {
					auto const x = half * dx * quadRules.x[q];
					state_type u = U[n][0];
					for (int l = 1; l < P; l++) {
						const auto ul = U[n][l];
						const auto al = an[l];
						const auto vl = basis.function(l, x);
						u += al * ul * vl;
					}
					f[q] = primToFlux<Real, 1>(conToPrim<Real, 1>(u));
				}
				for (int l = 0; l < P; l++) {
					auto const fp = F[n + 1];
					auto const fm = F[n];
					Real const Pp = basis.function(l, +half * dx);
					Real const dp = basis.derivative(l, +half * dx);
					Real const Pm = basis.function(l, -half * dx);
					Real const dm = basis.derivative(l, -half * dx);
					Real const wt = half * dx * quadRules.w[0];
					state_type volInt(zero);
					state_type surfInt(zero);
					surfInt += Pp * fp;
					surfInt -= Pm * fm;
					volInt += wt * dp * fp;
					volInt += wt * dm * fm;
					for (int q = 1; q < Q - 1; q++) {
						Real const x = half * dx * quadRules.x[q];
						Real const wt = half * dx * quadRules.w[q];
						Real const dpdx = basis.derivative(l, x);
						volInt += wt * dpdx * f[q];
					}
					dUdt[n][l] = (volInt - surfInt) * pow(dx, -l - 1);
				}
			}
			Real const rkWt1 = beta[rk];
			Real const rkWt0 = one - rkWt1;
			for (int n = BW; n < N - BW; n++) {
				for (int l = 0; l < P; l++) {
					U[n][l] += rkWt1 * dt * dUdt[n][l] + rkWt0 * (U0[n][l] - U[n][l]);
				}
			}
		}
		t += dt;
		n++;
	}
	{
		std::string fname = "test." + std::to_string(P);
		FILE *fp = fopen(fname.c_str(), "wt");
		for (int n = BW; n < N - BW; n++) {
			fprintf(fp, "\t%i \t%e %e %e %e %e %e \n", (n - BW), U[n][0][0], U[n][1][0], U[n][0][1], U[n][1][1],
					U[n][0][2], U[n][1][2]);
		}
		fclose(fp);
	}
}

