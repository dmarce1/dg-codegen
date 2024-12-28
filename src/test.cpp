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

constexpr int N = 256;
constexpr int BW = 2;
constexpr int P = 6;
constexpr int Nfield = ConservedState<Real, 1>::Nfield;
static const Real dx = Real(1) / Real(N - 1);

using state_type = ConservedState<Real, 1>;

std::vector<Vector<ConservedState<Real, 1>, P>> U(N);
std::vector<Vector<ConservedState<Real, 1>, P>> U0(N);
std::vector<ConservedState<Real, 1>> F(N);
std::vector<Vector<Real, Nfield>> UL(N);
std::vector<Vector<Real, Nfield>> UR(N);

void compute() {
	HLLCSolver<Real, 1> riemann;
	RoeSolver<Real, 1> riemann2;
	Real t = Real(0);
	Real tMax = Real(.12);
	for (int n = 0; n < N; n++) {
		Vector<ConservedState<Real, 1>, P> u;
		Real vel = Real(0.0);
		Real ene1 = Real(2.5);
		Real ene2 = Real(.125);
		Real rho1 = Real(1);
		Real rho2 = Real(.1);
		if (n < N / 2) {
			u[0].momentumDensity[0] = rho1 * vel;
			u[0].massDensity = rho1;
			u[0].energyDensity = ene1 + rho1 * vel * vel * Real(0.5);
		} else {
			u[0].momentumDensity[0] = rho2 * vel;
			u[0].massDensity = rho2;
			u[0].energyDensity = ene2 + rho2 * vel * vel * Real(0.5);
		}
		for (int p = 1; p < P; p++) {
			(Vector<Real, Nfield>&) u[p] = Real(0);
		}
		U[n] = std::move(u);
	}
	int n = 0;
	Basis<P> const basis;
	Real dt;
	//const std::array beta = { Real(1) };
	//const std::array beta = { Real(1), Real(0.5) };
	const std::array beta = { Real(1), Real(0.25), Real(2.0 / 3.0) };
	constexpr int rkOrder = beta.size();
	while (t < tMax) {
		Real lambdaMax = Real(0);
		U0 = U;
		for (int rk = 0; rk < rkOrder; rk++) {
			for (int n = 0; n < BW; n++) {
				U[n] = U[BW];
				U[N - n - 1] = U[N - BW - 1];
			}
			for (int n = 1; n < N; n++) {
				UL[n] = UR[n] = Real(0);
				for (int p = 0; p < P; p++) {
					auto const lN = U[n - 1][p].stateVector;
					auto const rN = U[n][p].stateVector;
					UL[n] += lN;
					UR[n] += rN * Real(nOnePow(p));
				}
			}
			int nlimits = 0;
			int ntotal = 0;
			for (int n = 1; n < N - 1; n++) {
				for (int f = 0; f < Nfield; f++) {
					auto const vp = U[n + 1][0].stateVector[f];
					auto const v0 = U[n][0].stateVector[f];
					auto const vm = U[n - 1][0].stateVector[f];
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
							U[n][p].stateVector[f] = Real(0);
						}
						slp = Real(0.5) * vanLeer(vp - v0, v0 - vm);
						slp = minmod(slp, U[n][1].stateVector[f]);
						U[n][1].stateVector[f] = slp;
						UL[n + 1][f] = v0 + slp;
						UR[n][f] = v0 - slp;
						nlimits++;
					}
					ntotal++;
				}
			}
			for (int n = 2; n < N - 1; n++) {
				auto const rc = riemann2(UL[n], UR[n]);
				auto const rc2 = riemann(UL[n], UR[n]);
				F[n] = rc.flux;
				lambdaMax = std::max(lambdaMax, rc.soundSpeed);
			}
			if (rk == 0) {
				dt = dx / (Real(2 * P + 1) * lambdaMax);
				printf("%i %e %e %e %i %i \n", n, t, dt, lambdaMax, nlimits, ntotal);
			}
			for (int n = 2; n < N - 2; n++) {
				std::vector<Vector<Real, Nfield>> up(P);
				std::vector<Vector<Real, Nfield>> fp(P);
				for (int q = 0; q < P; q++) {
					auto const x = basis.qPoint(q);
					for (int f = 0; f < Nfield; f++) {
						up[q][f] = Real(0);
						for (int p = 0; p < P; p++) {
							const auto aN = U[n][p].stateVector[f];
							up[q][f] += aN * basis.function(p, x);
						}
					}
				}
				for (int q = 0; q < P; q++) {
					fp[q] = PrimitiveState<Real, 1>(ConservedState<Real, 1>(up[q])).getFlux();
				}
				for (int p = 0; p < P; p++) {
					Real const sgnP = basis.function(p, +Real(0.5));
					Real const sgnM = basis.function(p, -Real(0.5));
					for (int f = 0; f < Nfield; f++) {
						Real part1, part2;
						part1 = -(sgnP * F[n + 1].stateVector[f] - sgnM * F[n].stateVector[f]);
						part2 = Real(0);
						for (int q = 0; q < P; q++) {
							auto const x = basis.qPoint(q);
							auto const wt = basis.qWeight(q);
							auto const dgdx = basis.derivative(p, x);
							part2 += wt * dgdx * fp[q][f];
						}
						U[n][p].stateVector[f] += (part1 + part2) * (dt / dx) * basis.massInv(p);
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
		fprintf(fp, "\t%i \t%e %e %e %e %e %e \n", n, U[n][0].stateVector[0], U[n][1].stateVector[0],
				U[n][0].stateVector[1], U[n][1].stateVector[1], U[n][0].stateVector[2], U[n][1].stateVector[2]);
	}
	fclose(fp);
}

