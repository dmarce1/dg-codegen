/******************************************************************************
 Copyright (C) 2024  Dominic C. Marcello
 *******************************************************************************/

#include <hpx/hpx_init.hpp>
#include "HLLC.hpp"
#include "Integrate.hpp"
#include "LegendreP.hpp"
#include "Polynomial.hpp"
#include "Real.hpp"
#include "Vector.hpp"
#include "Complex.hpp"
#include <unordered_map>
#include <numeric>
#include <stack>
#include <type_traits>
#include <debug/vector>
#include <span>
using namespace std;

#include "Complex.hpp"
#include "Real.hpp"
#include "HLLC.hpp"
#include <include/Quadrature.hpp>
#include "Sod.hpp"
#include <vector>

Real legendreP(int n, Real x) {
	Real const one(1.0);
	if (n == 0) {
		return one;
	} else if (n == 1) {
		return x;
	} else {
		Real Pnm1, Pn;
		Pnm1 = one;
		Pn = x;
		for (int l = 1; l < n; l++) {
			Real const a = Real(2 * l + 1) / Real(l + 1);
			Real const b = Real(l) / Real(l + 1);
			Real const Pnp1 = a * Pn * x - b * Pnm1;
			Pnm1 = Pn;
			Pn = Pnp1;
		}
		return Pn;
	}
}

Real dLegendrePdX(int n, Real x) {
	Real const zero(0), one(1), three(3);
	if (n == 0) {
		return zero;
	} else if (n == 1) {
		return one;
	} else {
		Real Pnm1, Pn;
		Real dPnm1dX, dPndX;
		Pnm1 = one;
		Pn = x;
		dPndX = one;
		for (int l = 1; l < n - 1; l++) {
			Real const a = Real(l + 1);
			Real const ainv = one / a;
			Real const b = Real(2 * l + 1) * ainv;
			Real const c = Real(l) * ainv;
			Real const Pnp1 = b * Pn * x - c * Pnm1;
			Real const dPnp1dX = a * Pn + x * dPndX;
			Pnm1 = Pn;
			Pn = Pnp1;
			dPndX = dPnp1dX;
		}
		return Real(n) * Pn + x * dPndX;
	}
}

template<int M, int P>
void compute() {
	using state_type = HydroState<Real, 1>;
	constexpr int BW = 2;
	constexpr int N = M + 2 * BW;
	Real const dx = Real(1) / Real(M);
	Real const half = Real(0.5);
	Real const quarter = Real(0.25);
	Real const third = Real(1.0 / 3.0);
	Real const sixth = Real(1.0 / 6.0);
	Real const ninth = Real(1.0 / 9.0);
	Real const zero = Real(0.0);
	Real const one = Real(1.0);
	Real const two = Real(2.0);
	Real const three = Real(3.0);
	const std::array<std::array<Real, 4>, 4> alpha1 = { { { one } } };
	const std::array<std::array<Real, 4>, 4> alpha2 = { { { one, zero }, { half, half } } };
	const std::array<std::array<Real, 4>, 4> alpha3 = { { { one, zero, zero }, { three * quarter, quarter, zero }, {
			third, zero, two * third } } };
	const std::array<std::array<Real, 4>, 4> alpha4 = { { { one, zero, zero, zero }, { half, half, zero, zero }, {
			ninth, two * ninth, two * third, zero }, { zero, third, third, third } } };
	const std::array<std::array<Real, 4>, 4> beta1 = { { { one } } };
	const std::array<std::array<Real, 4>, 4> beta2 = { { { one, zero }, { zero, half } } };
	const std::array<std::array<Real, 4>, 4> beta3 = { { { one, zero, zero }, { zero, quarter, zero }, { zero, zero, two
			* third } } };
	const std::array<std::array<Real, 4>, 4> beta4 = { { { half, zero, zero, zero }, { -quarter, half, zero, zero }, {
			-ninth, -third, one, zero }, { zero, sixth, zero, sixth } } };
	const std::array<std::array<Real, 4>, 4> alpha = (
			(P == 1) ? alpha1 : ((P == 2) ? alpha2 : ((P == 3) ? alpha3 : alpha4)));
	const std::array<std::array<Real, 4>, 4> beta = ((P == 1) ? beta1 : ((P == 2) ? beta2 : ((P == 3) ? beta3 : beta4)));
	const std::array cfl = { one, one, one, two * third };
	constexpr int rkOrder = std::min(P, 4);
	HLLCSolver<Real, 1> riemann;
	std::vector<std::vector<Vector<state_type, P>>> Uk(rkOrder + 1, std::vector<Vector<state_type, P>>(N));
	std::vector<std::vector<Vector<state_type, P>>> dUdt(rkOrder, std::vector<Vector<state_type, P>>(N));
	std::vector<state_type> dUL(N);
	std::vector<state_type> dUR(N);
	std::vector<state_type> F(N);
	Real const eps(0.1);
	constexpr int Q = P + 1;
	auto const quadRules = gaussLobattoRules<Q>();
	Real t = Real(0);
	Real tMax = Real(.125);
	for (int n = 0; n < N; n++) {
		std::vector<Vector<state_type, P>> &U = Uk[rkOrder];
		HydroState<Real, 1> v;
		if (n > N / 2) {
			v[D_i] = Real(sod_init.rhor);
			v[S_i] = zero;
			v[E_i] = Real(sod_init.pr / sod_init.rhor / (sod_init.gamma - 1.0));
		} else {
			v[D_i] = Real(sod_init.rhol);
			v[S_i] = zero;
			v[E_i] = Real(sod_init.pl / sod_init.rhol / (sod_init.gamma - 1.0));
		}
		U[n][0] = primToCon<Real, 1>(v);
		for (int l = 1; l < P; l++) {
			U[n][l] = zero;
		}
		/*		Vector<state_type, P> u;
		 Real const two(2);
		 Real const x0 = Real(n + 0.5) * dx;
		 for (int l = 0; l < P; l++) {
		 auto const u1 = quadRules.integrate(std::function<Real(Real)>([eps, x0, one, l, two, dx](Real x) {
		 return (one + eps * cos(x + x0)) * legendreP(l, two * x / dx);
		 }), dx) / dx;
		 auto const u2 = quadRules.integrate(std::function<Real(Real)>([eps, x0, one, half, l, two, dx](Real x) {
		 return (half + half * (one + eps * cos(x + x0))) * legendreP(l, two * x / dx);
		 }), dx) / dx;
		 U[n][l][D_i] = u1 * Real(2 * l + 1);
		 U[n][l][S_i] = u1 * Real(2 * l + 1);
		 U[n][l][E_i] = u2 * Real(2 * l + 1);
		 }*/
	}
	int n = 0;
	Real dt;
	printf("Starting...\n");
	const auto basisV = [dx](int l, Real x) {
		Real const two(2);
		return legendreP(l, two * x / dx);
	};
	const auto dBasisVdX = [dx](int l, Real x) {
		Real const two(2);
		return dLegendrePdX(l, two * x / dx) * two / dx;
	};
	for (int rk = 0; rk < rkOrder; rk++) {
		for (int k = 0; k <= rk; k++) {
			printf("%i % i %e %e\n", rk, k, double(alpha[rk][k]), double(beta[rk][k]));
		}
	}
//	return;
	int tn = 0;
	const auto output_check = [&]() {
		if (10.0 * t / tMax >= tn) {
			std::string name = "data." + std::to_string(M) + "." + std::to_string(P) + "." + std::to_string(tn)
					+ ".dat";
			FILE *fp = fopen(name.c_str(), "wt");
			std::vector<Vector<state_type, P>> &U = Uk[rkOrder];
			for (int n = BW; n < N - BW; n++) {
				Real const x = (Real(n - N / 2) + half) / Real(N);
				sod_state_t sod;
				exact_sod(&sod, &sod_init, x, t, dx);
				HydroState<Real, 1> vA;
				vA[rho_i] = Real(sod.rho);
				vA[vel_i] = Real(sod.v);
				vA[eps_i] = Real(sod.p / sod.rho / (sod_init.gamma - 1.0));
				auto const v = conToPrim<Real, 1>(U[n][0]);
				fprintf(fp, "%e ", x);
				for (int f = 0; f < Nfield; f++) {
					fprintf(fp, "%e ", v[f]);
					fprintf(fp, "%e ", vA[f]);
					fprintf(fp, "%e ", (v[f] - vA[f]) / (half * (Real(1.0e-100) + v[f] + vA[f])));
				}
				fprintf(fp, "\n");
			}
			fclose(fp);
			tn++;
		}
	};
	while (t < tMax) {
		output_check();
		constexpr int Q = P + 2;
		auto const quadRules = gaussLobattoRules<Q>();
		Uk[0] = Uk[rkOrder];
		for (int rk = 0; rk < rkOrder; rk++) {
			std::vector<Vector<state_type, P>> &U = Uk[rk];
			for (int n = 0; n < BW; n++) {
				U[n] = U[n + N - 2 * BW];
				U[N - n - 1] = U[2 * BW - n - 1];
			}
			auto Uc = U;
			int nlim = 0;
			int ntot = 0;
			for (int n = 0; n < N; n++) {
				auto const v = conToPrim<Real, 1>(Uc[n][0]);
				auto const Rinv = primToEigenvectors<Real, 1>(v);
				auto const R = matrixInverse(Rinv);
				for (int p = 1; p < P; p++) {
					Uc[n][p] = toVector(R * toColumnVector(Uc[n][p]));
				}
			}
			for (int f = 0; f < Nfield; f++) {
				std::vector<bool> limitDone(N, false);
				for (int p = P - 1; p > 0; p--) {
					int const pm1 = p - 1;
					for (int n = 1; n < N - 1; n++) {
						if (!limitDone[n]) {
							Real const norm = one / Real(2 * pm1 + 1);
							auto const cL = Uc[n][pm1][f] - Uc[n - 1][pm1][f];
							auto const cR = Uc[n + 1][pm1][f] - Uc[n][pm1][f];
							auto ulim = minmod(Uc[n][p][f], minmod(cL, cR) * norm);
							if (ulim != Uc[n][p][f]) {
								Uc[n][p][f] = ulim;
								nlim++;
							} else {
								limitDone[n] = true;
							}
						}
						ntot++;
					}
				}
			}
			for (int n = 0; n < N; n++) {
				auto const v = conToPrim<Real, 1>(U[n][0]);
				auto const R = primToEigenvectors<Real, 1>(v);
				for (int p = 1; p < P; p++) {
					U[n][p] = toVector(R * toColumnVector(Uc[n][p]));
				}
			}
			for (int n = 1; n < N - 1; n++) {
				dUL[n] = state_type(zero);
				dUR[n] = state_type(zero);
				for (int l = 1; l < P; l++) {
					const auto ul = U[n][l];
					dUL[n] += ul * basisV(l, +half * dx);
					dUR[n] -= ul * basisV(l, -half * dx);
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
				dt = Real(0.8) * cfl[rkOrder - 1] * dx / (Real(2 * P - 1) * lambdaMax);
				dt = std::min(dt, tMax - t);
				printf("%i %e %e %e %i %i\n", n, t, dt, lambdaMax, nlim, ntot);
			}
			for (int n = 2; n < N - 2; n++) {
				Vector<state_type, Q> f;
				int const bw = quadRules.withEndpoints ? 1 : 0;
				for (int q = bw; q < Q - bw; q++) {
					auto const x = half * dx * quadRules.x[q];
					state_type u = U[n][0];
					for (int l = 1; l < P; l++) {
						const auto ul = U[n][l];
						const auto vl = basisV(l, x);
						u += ul * vl;
					}
					f[q] = primToFlux<Real, 1>(conToPrim<Real, 1>(u));
				}
				for (int l = 0; l < P; l++) {
					auto const fp = F[n + 1];
					auto const fm = F[n];
					Real const Pp = basisV(l, +half * dx);
					Real const Pm = basisV(l, -half * dx);
					state_type volInt(zero);
					state_type surfInt(zero);
					surfInt += Pp * fp;
					surfInt -= Pm * fm;
					if (quadRules.withEndpoints) {
						f[0] = fm;
						f[Q - 1] = fp;
					}
					for (int q = 0; q < Q; q++) {
						Real const x = half * dx * quadRules.x[q];
						Real const wt = half * dx * quadRules.w[q];
						Real const dpdx = dBasisVdX(l, x);
						volInt += wt * dpdx * f[q];
					}
					dUdt[rk][n][l] = (volInt - surfInt) * (Real(2 * l + 1) / dx);
				}
			}
			int const i = rk;
			for (int n = BW; n < N - BW; n++) {
				for (int l = 0; l < P; l++) {
					Uk[i + 1][n][l] = zero;
					for (int k = 0; k <= i; k++) {
						Uk[i + 1][n][l] += alpha[i][k] * Uk[k][n][l] + beta[i][k] * dt * dUdt[k][n][l];
					}
				}
			}
		}
		t += dt;
		n++;
		output_check();
	}
	{
		/*		double l1 = 0.0;
		 for (int n = BW; n < N - BW; n++) {
		 std::vector<Vector<state_type, P>> &U = Uk[rkOrder];
		 Real const x0 = Real(n + 0.5) * dx;
		 Real const D = quadRules.integrate(std::function<Real(Real)>([eps, x0, one, dx, t](Real x) {
		 return (one + eps * cos(x + x0 - t));
		 }), dx) / dx;
		 double const err = abs(U[n][0][0] - D);
		 l1 += err * dx;
		 }
		 FILE *fp = fopen("err.dat", "at");
		 fprintf(fp, "%4i %4i %12.4e\n", M, P, l1);
		 fclose(fp);
		 */
	}
}

Math::Polynomial<Real> lagrangePolynomial(int k, int j, std::function<Real(int)> xn = nullptr) {
	if (xn == nullptr) {
		xn = [k](int n) {
			return Real(2 * n - k) / Real(k + 2);
		};
	}
	Math::Polynomial<Real> p;
	p[0] = Real(1);
	Real xj = xn(j);
	for (int m = 0; m <= k; m++) {
		if (m != j) {
			Real xm = xn(m);
			Real c0 = Real(1) / (xj - xm);
			Math::Polynomial<Real> q;
			q[0] = -xm * c0;
			q[1] = c0;
			p *= q;
		}
	}
	return p;
}

int hpx_main(int argc, char *argv[]) {
	compute<256, 1>();
	compute<256, 2>();
	compute<256, 3>();
	compute<256, 4>();
	return hpx::local::finalize();
}

int main(int argc, char *argv[]) {
	auto rc = hpx::init(argc, argv);
	return rc;
}
