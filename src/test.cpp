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
constexpr int P = 4;
static const Real dx = Real(1) / Real(N - 1);

using state_type = HydroState<Real, 1>;

std::vector<Vector<state_type, P>> U(N);
std::vector<Vector<state_type, P>> U0(N);
std::vector<state_type> F(N);
std::vector<state_type> UL(N);
std::vector<state_type> UR(N);

void compute() {
	HLLCSolver<Real, 1> riemann;
	std::vector<state_type> U(N);
	std::vector<state_type> U0(N);
	std::vector<state_type> UL(N);
	std::vector<state_type> UR(N);
	std::vector<state_type> F(N);
	std::vector<state_type> S(N);

	Real tMax = Real(.1);
	for (int n = 0; n < N; n++) {
		state_type u;
		Real vel1 = Real(5);
		Real vel2 = Real(-5);
		Real rho1 = Real(1);
		Real rho2 = Real(1);
		Real ene1 = Real(1);
		Real ene2 = Real(1);
		if (n < N / 4 || n >= 3 * N / 4) {
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
//	const std::array beta = { Real(1) };
	const std::array beta = { Real(1), Real(0.5) };
//	const std::array beta = { Real(1), Real(0.25), Real(2.0 / 3.0) };
	constexpr int rkOrder = beta.size();
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
			for (int n = 1; n < N - 1; n++) {
				Vector<Real, Nfield> dW;
				auto const R = primToEigenvectors<Real, 1>(conToPrim<Real, 1>(U[n]));
				auto const invR = Math::matrixInverse(R);
				auto const dWm = toVector(invR * toColumnVector(U[n] - U[n - 1]));
				auto const dWp = toVector(invR * toColumnVector(U[n + 1] - U[n]));
				for (int f = 0; f < Nfield; f++) {
					dW[f] = vanLeer(dWp[f], dWm[f]);
				}
				S[n] = toVector(R * toColumnVector(dW));
			}
			for (int n = 1; n < N - 1; n++) {
				UR[n] = U[n] - Real(0.5) * S[n];
				UL[n + 1] = U[n] + Real(0.5) * S[n];
			}
			for (int n = 2; n < N - 1; n++) {
				auto const rc = riemann(UL[n], UR[n]);
				F[n] = rc.flux;
				lambdaMax = std::max(lambdaMax, rc.signalSpeed);
			}
			if (rk == 0) {
				dt = Real(0.4) * dx / lambdaMax;
				lambda = dt / dx;
				printf("%i %e %e %e\n", n, t, dt, lambdaMax);
			}
			for (int n = 2; n < N - 2; n++) {
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

void compute2() {
	HLLCSolver<Real, 1> riemann;
	//RoeSolver<Real, 1> riemann;
	Real t = Real(0);
	Real tMax = Real(.01);
	for (int n = 0; n < N; n++) {
		Vector<state_type, P> u;
		Real vel1 = Real(2);
		Real vel2 = Real(-2);
		Real rho1 = Real(1);
		Real rho2 = Real(1);
		Real ene1 = Real(1);
		Real ene2 = Real(1);
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
	for (int n = 0; n < N; n++) {
		auto &u = U[n][1];
		for (int f = 0; f < Nfield; f++) {
			auto const vp = U[n + 1][0][f];
			auto const v0 = U[n][0][f];
			auto const vm = U[n - 1][0][f];
			u[f] = minmod_theta(vp - v0, v0 - vm, Real(1));
		}
	}
	int n = 0;
	Basis<P> const basis;
	Real dt;
	//const std::array beta = { Real(1) };
	//	const std::array beta = { Real(1), Real(0.5) };
	const std::array beta = { Real(1), Real(0.25), Real(2.0 / 3.0) };
	constexpr int rkOrder = beta.size();
	while (t < tMax) {
		Real lambdaMax = Real(0);
		U0 = U;
		for (int rk = 0; rk < rkOrder; rk++) {
			for (int n = 0; n < BW; n++) {
				U[n] = U[n + N - 2 * BW];
				U[N - n - 1] = U[2 * BW - n - 1];
			}
			for (int n = 1; n < N; n++) {
				UL[n] = UR[n] = Real(0);
				for (int p = 0; p < P; p++) {
					auto const lN = U[n - 1][p];
					auto const rN = U[n][p];
					UL[n] += lN;
					UR[n] += rN * Real(nOnePow(p));
				}
			}
			int nlimits = 0;
			int ntotal = 0;
			for (int n = 1; n < N - 1; n++) {
				for (int f = 0; f < Nfield; f++) {
					auto const vp = U[n + 1][0][f];
					auto const v0 = U[n][0][f];
					auto const vm = U[n - 1][0][f];
					auto const ul = UL[n + 1][f];
					auto const ur = UR[n][f];
					auto slp = minmod(vp - v0, v0 - vm);
					bool doLimiter = false;
					doLimiter = doLimiter || ((ul - v0) * slp < Real(0));
					doLimiter = doLimiter || (abs(ul - v0) > abs(slp));
					doLimiter = doLimiter || ((v0 - ur) * slp < Real(0));
					doLimiter = doLimiter || (abs(v0 - ur) > abs(slp));
					if (doLimiter) {
						for (int p = 2; p < P; p++) {
							U[n][p][f] = Real(0);
						}
						slp = vanLeer(vp - v0, v0 - vm);
						//	slp = Real(0.5) * minmod_theta(vp - v0, v0 - vm, Real(2));
						slp = minmod(slp, U[n][1][f]);
						U[n][1][f] = slp;
						UL[n + 1][f] = v0 + slp;
						UR[n][f] = v0 - slp;
						nlimits++;
					}
					ntotal++;
				}
			}
			for (int n = 2; n < N - 1; n++) {
				auto const rc = riemann(UL[n], UR[n]);
				F[n] = rc.flux;
				lambdaMax = std::max(lambdaMax, rc.signalSpeed);
			}
			//printf( "%e\n", lambdaMax);
			if (rk == 0) {
				dt = dx / (Real(2 * P + 1) * lambdaMax);
				printf("%i %e %e %e %i %i \n", n, t, dt, lambdaMax, nlimits, ntotal);
			}
			for (int n = 2; n < N - 2; n++) {
				std::vector<state_type> up(P);
				std::vector<state_type> fp(P);
				for (int q = 0; q < P; q++) {
					auto const x = basis.qPoint(q);
					for (int f = 0; f < Nfield; f++) {
						up[q][f] = Real(0);
						for (int p = 0; p < P; p++) {
							const auto aN = U[n][p][f];
							up[q][f] += aN * basis.function(p, x);
						}
					}
				}
				for (int q = 0; q < P; q++) {
					fp[q] = primToFlux<Real, 1>(conToPrim<Real, 1>(up[q]));
				}
				for (int p = 0; p < P; p++) {
					Real const sgnP = basis.function(p, +Real(0.5));
					Real const sgnM = basis.function(p, -Real(0.5));
					for (int f = 0; f < Nfield; f++) {
						Real part1, part2;
						part1 = -(sgnP * F[n + 1][f] - sgnM * F[n][f]);
						part2 = Real(0);
						for (int q = 0; q < P; q++) {
							auto const x = basis.qPoint(q);
							auto const wt = basis.qWeight(q);
							auto const dgdx = basis.derivative(p, x);
							part2 += wt * dgdx * fp[q][f];
						}
						U[n][p][f] += (part1 + part2) * (dt / dx) * basis.massInv(p);
					}
				}
			}
			for (int n = BW; n < N - BW; n++) {
				Real const w1 = beta[rk];
				Real const w0 = Real(1) - w1;
				for (int p = 0; p < P; p++) {
					U[n][p] = w1 * U[n][p] + w0 * U0[n][p];
				}
			}
		}
		t += dt;
		n++;
	}
	FILE *fp = fopen("test.dat", "wt");
	for (int n = 0; n < N; n++) {
		fprintf(fp, "\t%i \t%e %e %e %e %e %e \n", n, U[n][0][0], U[n][1][0], U[n][0][1], U[n][1][1], U[n][0][2],
				U[n][1][2]);
	}
	fclose(fp);
}

