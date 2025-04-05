/******************************************************************************
 Copyright (C) 2024  Dominic C. Marcello
 *******************************************************************************/
#include "ValArray.hpp"
#include <optional>
#include "Metric.hpp"
#include "Timer.hpp"
#include <hpx/hpx_init.hpp>
#include "Tensor.hpp"
#include "Interpolate.hpp"
#include "IndexTuple.hpp"
#include "Hydrodynamics.hpp"
#include "LegendreP.hpp"
#include "Polynomial.hpp"
#include "Permutation.hpp"
#include "Real.hpp"
#include "HydroGrid.hpp"
#include "Options.hpp"
#include "Radiation.hpp"
#include "Relativity.hpp"
#include "SparseMatrix.hpp"
#include "Root.hpp"
#include "Vector.hpp"
#include <unordered_map>
#include <numeric>
#include <stack>
#include <type_traits>
#include <debug/vector>
#include <span>
using namespace std;
#include "TriangularArray.hpp"
#include "Real.hpp"
#include "Interval.hpp"
#include "Limiters.hpp"
#include "Quadrature.hpp"
#include "Sod.hpp"
#include "Zobrist.hpp"
#include <random>
#include <vector>
#include <fstream>

template<int M, int rho>
void compute() {
	using namespace Hydrodynamics;
	using state_type = ConservedState<Real, 1>;
	using prim_type = PrimitiveState<Real, 1>;
	constexpr int BW = 2;
	constexpr int N = M + 2 * BW;
	constexpr int NF = state_type::NFields;
	Real const dx = Real(1) / Real(M);
	Real const half = Real(0.5);
	Real const zero = Real(0.0);
	Real const one = Real(1.0);
	std::vector<state_type> U(N);
	std::vector<state_type> S(N);
	std::vector<state_type> U0(N);
	std::vector<state_type> dUdt(N);
	std::vector<state_type> UL(N);
	std::vector<state_type> UR(N);
	std::vector<state_type> F(N);
	Real tMax = Real(.01);
	Real t = Real(0);
	Real dt;
	for (int n = 0; n < N; n++) {
		state_type v;
		if (n > N / 2) {
			v.D = Real(sod_init.rhor);
			v.S = zero;
			v.E = Real(sod_init.pr / (sod_init.gamma - 1.0));
			v.tau = Real(energy2Entropy(sod_init.rhor, sod_init.pr / (sod_init.gamma - 1.0)));
		} else {
			v.D = Real(sod_init.rhol);
			v.S = zero;
			v.E = Real(sod_init.pl / (sod_init.gamma - 1.0));
			v.tau = Real(energy2Entropy(sod_init.rhol, sod_init.pl / (sod_init.gamma - 1.0)));
		}
		U[n] = state_type(v);
	}
	int n = 0;
	printf("Starting...\n");
	while (t < tMax) {
		U0 = U;
		for (int n = BW - 1; n < BW + 1; n++) {
			const auto maxE = max(max(U0[n].E, U0[n + 1].E), U0[n - 1].E);
			U[n].dualEnergyUpdate(maxE);
		}
		U = U0;
		for (int n = BW; n < N - BW; n++) {

		}
		for (int rk = 0; rk < BW; rk++) {
			for (int n = 0; n < BW; n++) {
				U[n] = U[BW];
				U[N - n - 1] = U[N - BW - 1];
			}
			for (int n = 1; n < N - 1; n++) {
				state_type const up = U[n + 1];
				state_type const u0 = U[n];
				state_type const um = U[n - 1];
				for (int f = 0; f < NF; f++) {
					S[n][f] = minmod(u0[f] - um[f], up[f] - u0[f]);
				}
			}
			for (int n = 1; n < N - 1; n++) {
				UL[n + 1] = U[n] + half * S[n];
				UR[n] = U[n] + half * S[n];
			}
			Real lambdaMax = Real(0);
			for (int n = BW; n < N - BW + 1; n++) {
				auto const rc = riemannSolver(UL[n], UR[n]);
				F[n] = rc.flux;
				lambdaMax = std::max(lambdaMax, rc.signalSpeed);
			}
			if (rk == 0) {
				dt = Real(0.4) * dx / lambdaMax;
				dt = std::min(dt, tMax - t);
				printf("%i %e %e %e\n", n, t, dt, lambdaMax);
			}
			for (int n = BW; n < N - BW; n++) {
				dUdt[n] = -(F[n + 1] - F[n]) * (one / dx);
			}
			if (rk == 0) {
				for (int n = BW; n < N - BW; n++) {
					U[n] += dt * dUdt[n];
				}
			} else {
				for (int n = BW; n < N - BW; n++) {
					U[n] = half * (U0[n] + U[n] + dt * dUdt[n]);
				}
			}
		}
		t += dt;
		n++;
	}
	{
		FILE *fp = fopen("test.txt", "wt");
		double l1 = 0.0;
		constexpr sod_init_t sod_init = { 1.0, 0.125, 1.0, 0.1, 5.0 / 3.0 };
		for (int n = BW; n < N - BW; n++) {
			sod_state_t sod;
			prim_type V(U[n]);
			double const x = double((n - N / 2) + 0.5) * double(dx);
			exact_sod(&sod, &sod_init, x, t, dx);
			double const err = abs(double(U[n].D) - sod.rho);
			l1 += err * dx;
			fprintf(fp, "%e %e %e %e %e\n", x, U[n].D, sod.rho, U[n].tau, energy2Entropy(V.rho, V.rho * V.eps));
		}
		fclose(fp);
		printf("Error = %12.4e\n", l1);
	}
}

void testIntegers();
void testRadiation();

template<typename T>
Polynomial<T> legendrePolynomial(int n) {
	Polynomial<T> Pn;
	if (n == 0) {
		Pn[0] = Real(1);
	} else {
		if (n % 2 == 0) {
			Pn[0] = Real(1);
			Pn[1] = Real(0);
		} else {
			Pn[0] = Real(0);
			Pn[1] = Real(1);
		}
		for (int k = 0; k <= n - 2; k++) {
			Pn[k + 2] = -Real((n - k) * (n + k + 1)) / Real((k + 1) * (k + 2)) * Pn[k];
		}
		Real norm = Real(0);
		for (int k = 0; k <= n; k++) {
			norm += Pn[k];
		}
		norm = Real(1) / norm;
		for (int k = 0; k <= n; k++) {
			Pn[k] *= norm;
		}
	}
	return Pn;
}

#include <iostream>
#include <vector>
#include <array>
#include <cmath>
#include <cassert>

//-----------------------------------------------------------------
// Utility functions
//-----------------------------------------------------------------

#include <iostream>
#include <vector>
#include <array>
#include <cmath>
#include <cassert>

template<typename T>
auto genFunctor(std::string const &name, std::vector<Polynomial<T>> const &P, int Order, int shift = 0, int Dimension = 3) {
	static constexpr T zero(0), half(0.5), one(1), two(2);
	std::ostringstream code;
	code << std::scientific;
	code.precision(16);
	code << "template<typename T> \n";
	code << "struct " << name << "<T, " << Order << ", " << Dimension << "> {\n";
	code << "\tinline T operator()(T const& x) const {\n";
	code << "\t\tusing namespace std;\n";
	code << "\t\tT y;\n";
	int NumPieces;
	const int n = Order;
	const bool isOdd = (n % 2 == 1);
	int m;
	if (isOdd) {
		m = (n - 1) / 2;
		NumPieces = m + 1;
	} else {
		m = n / 2;
		NumPieces = m;
	}
	if (isOdd) {
		code << "\t\tswitch( floor(x * T(" << T(0.5 * n) << ") + T(" << T(0.5 * n - m) << ")) ) {\n";
	} else {
		code << "\t\tswitch( floor(x) ) {\n";
	}
	for (int p = 0; p < NumPieces; p++) {
		bool first = true;
		code << "\t\tcase " << p << ":\n";
		for (int i = n - 1; i >= 0; i--) {
			T const c = P[p][i];
			if (first) {
				if (c > T(0)) {
					code << "\t\t\ty = +T(" << c << ");\n";
					first = false;
				} else if (c < T(0)) {
					code << "\t\t\ty = -T(" << -c << ");\n";
					first = false;
				}
			} else {
				if (c > T(0)) {
					code << "\t\t\ty = fma(y, x, +T(" << c << "));\n";
				} else if (c < T(0)) {
					code << "\t\t\ty = fma(y, x, -T(" << -c << "));\n";
				} else {
					code << "\t\t\ty *= x;\n";
				}
			}
		}
		code << "\t\t\tbreak;\n";
	}
	code << "\t\tdefault:\n";
	code << "\t\t\ty = T(0);\n";
	code << "\t\t}\n";
	if (shift < 0) {
		code << "\t\tT const xinv = T(1) / x;\n";
		while (shift) {
			code << "\t\ty *= xinv;\n";
			shift++;
		}
	}
	code << "\t\treturn y;\n";
	code << "\t}\n";
	code << "};\n";
	return code.str();
}

auto genBSplineFunctor(int Order, int Dimension = 3) {
	using T = Real;
	static constexpr T zero(0), half(0.5), one(1), two(2);
	int NumPieces;
	std::vector<T> boundaries;
	std::vector<Polynomial<T>> rho;
	std::vector<Polynomial<T>> Menc;
	std::vector<Polynomial<T>> drho_dr;
	const int n = Order;
	const bool isOdd = (n % 2 == 1);
	int m;
	if (isOdd) {
		m = (n - 1) / 2;
		NumPieces = m + 1;
	} else {
		m = n / 2;
		NumPieces = m;
	}
	boundaries.resize(NumPieces + 1);
	if (isOdd) {
		boundaries[0] = zero;
		boundaries[1] = one / T(n);
		for (int r = 1; r < NumPieces; ++r) {
			boundaries[r + 1] = boundaries[1] + T(r) * (two / T(n));
		}
	} else {
		boundaries.resize(NumPieces + 1);
		for (int i = 0; i <= NumPieces; ++i) {
			boundaries[i] = T(i) * (two / T(n));
		}
	}
	rho.resize(NumPieces);
	for (int r = 0; r < NumPieces; ++r) {
		int k_full = m + r;
		for (int i = 0; i < n; ++i) {
			T sum_j = zero;
			for (int j = 0; j <= k_full; ++j) {
				T term = integerPower<T>(one - two * T(j) / T(n), n - 1 - i);
				sum_j += negativeOne2Power<T>(j) * nChooseK<T>(n, j) * term;
			}
			rho[r][i] = (nChooseK<T>(n - 1, i) * sum_j) / nFactorial<T>(n - 1);
		}
	}
	Polynomial<T> dV;
	using namespace std;
	dV[Dimension - 1] = two * pow(T(M_PI), T(Dimension) * half) / tgamma(T(Dimension) * half);
	T norm = zero;
	for (int r = 0; r < NumPieces; ++r) {
		norm += polynomialIntegrate(rho[r] * dV, boundaries[r], boundaries[r + 1]);
	}
	norm = one / norm;
	for (int r = 0; r < NumPieces; ++r) {
		rho[r] *= norm;
	}
	T Mtot = zero;
	for (int r = 0; r < NumPieces; ++r) {
		for (int i = 0; i < n; ++i) {
			if (abs(rho[r][i]) < T(1e-10)) {
				rho[r][i] = T(0);
			}
		}
		auto const M = polynomialAntiDerivative(rho[r] * dV);
		T const M0 = M(boundaries[r]);
		T const M1 = M(boundaries[r + 1]);
		Mtot += M1 - M0;
		Menc.push_back(M + Polynomial<T>(Mtot - M0));
		drho_dr.push_back(polynomialDerivative(rho[r]));
	}
	std::ostringstream code;
	code << "\n";
	code << genFunctor<T>("BSplineDensity", rho, Order, Dimension);
	code << "\n";
	code << genFunctor<T>("BSplineDerivative", drho_dr, Order, Dimension);
	code << "\n";
	code << genFunctor<T>("BSplineEnclosedMass", Menc, Order, -2, Dimension);
	code << "\n";
	return code.str();
}

void test();
void testStrings();
void testEinstein();

#include <cuda/std/mdspan>
#include <unordered_set>

template<size_t D, size_t R>
struct TensorIndexing {
	using poly_t = double;
	using index_t = IndexTuple<D, R>;
	TensorIndexing(std::vector<Permutation<R>> const &permutations, std::vector<bool> anti = std::vector<bool>()) {
		std::unordered_set<Permutation<R>> antisymmetric;

		for (size_t i = 0; i < anti.size(); i++) {
			auto const &p = permutations[i];
			if (anti[i] && (p.parity() < 0)) {
				antisymmetric.insert(p);
			}
		}
		auto const isAntisymmetric = [&antisymmetric](Permutation<R> const &P) {
			return antisymmetric.find(P) != antisymmetric.end();
		};

		auto const isSymmetric = [&antisymmetric](Permutation<R> const &P) {
			return antisymmetric.find(P) == antisymmetric.end();
		};
		auto const inOrder = [&permutations, &isSymmetric](index_t const &indices) {
			auto const &I = indices;
			for (const auto &p : permutations) {
				auto const J = p.apply(I);
				if (J > I) {
					return false;
				}
			}
			return true;
		};
		for (size_t i = 0; i < R; i++) {
			std::array<std::array<int, D>, D> forwardDifference;
			size_t jLast = -1;
			size_t k = 0;
			size_t count = 0;
			index_t indices = index_t::begin();
			for (; indices != index_t::end(); indices++) {
				if (indices[i] != jLast) {
					jLast = indices[i];
					forwardDifference[0][k++] = count;
					if (k % D == 0) {
						break;
					}
				}
				if (inOrder(indices)) {
					count++;
				}
			}
			for (size_t j = 1; j < D; j++) {
				for (size_t k = 0; k < D - j; k++) {
					forwardDifference[j][k] = forwardDifference[j - 1][k + 1] - forwardDifference[j - 1][k];
				}
			}
			for (size_t k = 0; k < D; k++) {
				for (size_t j = 0; j < D - k; j++) {
					printf("%i ", forwardDifference[j][k]);
				}
				printf("\n");
			}
			printf("\n");
			double q = 1.0;
			Polynomial<double> xn;
			xn[0] = 1.0;
			for (size_t j = 0; j < D; j++) {
				polynomials[i] += q * xn * forwardDifference[j][0];
				q /= (j + 1);
				Polynomial<double> xm;
				xm[0] = -double(j);
				xm[1] = 1.0;
				xn = xn * xm;
			}
			std::cout << "rank = " << std::to_string(i) << " " << toString(polynomials[i]) << "\n";
		}

		for (index_t indices = index_t::begin(); indices != index_t::end(); indices++) {
			if (inOrder(indices)) {
				std::cout << toString(indices) << " ";
				size_t sum = 0;
				for (size_t r = 0; r < R; r++) {
					size_t const term = std::round(polynomials[r](indices[r]));
					std::cout << " + " << std::to_string(term) << " ";
					sum += term;
				}
				std::cout << " = " << sum << "\n";
			}
		}

	}
private:

	std::array<Polynomial<double>, R> polynomials;
};

int hpx_main(int argc, char *argv[]) {
	printf("Reading options...\n");
	processOptions(argc, argv);
	printf("Done.\n");
	using namespace Tensors;
	constexpr size_t D = 3;
	constexpr size_t R = 2;
	using symmetry_type = std::pair<int, Permutation<R>>;
	static constexpr Index<'i'> i;
	static constexpr Index<'j'> j;
	static constexpr Index<'k'> k;
	Tensor<double, D, R, symmetry_type( { 1, { 1, 2 } }), symmetry_type( { 1, { 2, 1 } })> A;
	Tensor<double, D, R, symmetry_type( { 1, { 1, 2 } }), symmetry_type( { 1, { 2, 1 } })> B;
	Tensor<double, D, R, symmetry_type( { 1, { 1, 2 } }), symmetry_type( { 1, { 2, 1 } })> C;
	C(i, k) = A(i, j) * B(j, k);
	return hpx::local::finalize();
}

namespace Math {
Real normalDistribution() {
	static std::default_random_engine gen(42);
	static std::normal_distribution d { 0.0, 1.0 };
	return Real(d(gen));
}

template<int N>
Vector<Real, N> leastSquares(std::vector<Real> y, std::vector<Real> x) {
	using namespace Math;
	Real const zero(0), one(1), two(2);
	int const M = y.size();
	std::vector<Vector<Real, N>> F(M);
	SquareMatrix<Real, N> J;
	Vector<Real, N> beta;
	for (int m = 0; m < M; m++) {
		Real const xm = x[m];
		Real xn;
		for (int n = 0; n < N; n++) {
			xn = (n == 0) ? one : (xn * xm);
			F[m][n] = -xn;
		}
	}
	for (int n = 0; n < N; n++) {
		for (int m = 0; m < N; m++) {
			J[n, m] = zero;
			for (int k = 0; k < M; k++) {
				J[n, m] += F[k][m] * F[k][n];
			}
		}
	}
	auto const iJ = matrixInverse(J);
	for (int n = 0; n < N; n++) {
		beta[n] = zero;
	}
	for (int m = 0; m < M; m++) {
		beta -= y[m] * toVector(iJ * toColumnVector(F[m]));
	}
	Real r2 = zero;
	for (int m = 0; m < M; m++) {
		Real const xm = x[m];
		Real xn;
		Real r = y[m];
		for (int n = 0; n < N; n++) {
			xn = (n == 0) ? one : (xn * xm);
			r -= beta[n] * xn;
		}
		r2 += abs(r);
	}
	printf("%e\n\n", r2 / M);

	return beta;
}

Real derivative(std::function<Real(Real)> const &f, Real x) {
	static Real const eps = pow(Real(std::numeric_limits<double>::epsilon()), Real(1.0 / 6.0));
	static Real const mindx = Real(std::numeric_limits<double>::min());
	Real const third(1.0 / 3.0), half(0.5), zero(0), one(1), four(4);
	Real const dx = max(eps * abs(x), mindx);
	Real const dxinv = one / dx;
	Real const dfdx1 = half * (f(x + dx) - f(x - dx)) * dxinv;
	Real const dfdx2 = (f(x + half * dx) - f(x - half * dx)) * dxinv;
	return (four * dfdx2 - dfdx1) * third;
}
}

int main(int argc, char *argv[]) {
	std::vector<std::string> cfg = { "hpx.commandline.allow_unknown=1" };
	cfg.push_back("hpx.stacks.small_size=524288");
	hpx::init_params init_params;
	init_params.cfg = std::move(cfg);
	auto rc = hpx::init(argc, argv, init_params);
	return rc;
}
