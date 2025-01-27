/*
 * Grid.hpp
 *
 *  Created on: Jan 26, 2025
 *      Author: dmarce1
 */

#ifndef INCLUDE_GRID_HPP_
#define INCLUDE_GRID_HPP_

#include "TriangularArray.hpp"
#include "LegendreP.hpp"
#include "Options.hpp"
#include "Quadrature.hpp"

#include <vector>

template<typename T, int D, int P, int Q>
struct Grid {
	static T constexpr zero = T(0), half = T(0.5), one = T(1), two = T(2);
	static int constexpr Q3 = Math::integerPower(Q, D);
	static int constexpr P3 = []() {
		int p3 = 1;
		int n = P;
		for (int d = 0; d < D; d++) {
			p3 *= n + d;
			p3 /= d + 1;
		}
		return p3;
	}();
	Grid(int N_) {
		N = N_;
		N3 = Math::integerPower<int>(N, D);
		for (int p3 = 0; p3 < P3; p3++) {
			U[p3].resize(N3);
		}
	}
	void transform(int dir) {
		static constexpr T i2d = Math::integerPower(half, D);
		if (dir < 0) {
			for (int q3 = 0; q3 < Q3; q3++) {
				V[q3] = std::vector<T>(N3, zero);
			}
			for (int q3 = 0; q3 < Q3; q3++) {
				auto const x = quadratureRules.x[q3];
				auto const pn = legendreBasis<T, D, P>(x);
				for (int n = 0; n < N3; n++) {
					for (int p3 = 0; p3 < P3; p3++) {
						V[q3][n] += pn[p3] * U[p3][n];
					}
				}
			}
			U = decltype(U)();
		} else if (dir > 0) {
			for (int p3 = 0; p3 < P3; p3++) {
				U[p3] = std::vector<T>(N3, zero);
			}
			for (int q3 = 0; q3 < Q3; q3++) {
				auto const x = quadratureRules.x[q3];
				auto const w = quadratureRules.w[q3];
				auto const pn = legendreBasis<T, D, P>(x);
				for (int n = 0; n < N3; n++) {
					for (int p3 = 0; p3 < P3; p3++) {
						U[q3][n] += i2d * pn[p3] * U[p3][n] * pNorm[p3];
					}
				}
			}
			V = decltype(V)();
		}
	}
private:
	int N;
	int N3;
	Math::Vector<std::vector<T>, P3> U;
	Math::Vector<std::vector<T>, Q3> V;
	static Quadrature::MultivariateQuadratureRules<T, Q, D> const quadratureRules;
	static TriangularArray<T, D, P> const pNorm;
};

template<typename T, int D, int P, int Q>
Quadrature::MultivariateQuadratureRules<T, Q, D> const Grid<T, D, P, Q>::quadratureRules(Quadrature::gaussLobattoRules<T, Q>());

template<typename T, int D, int P, int Q>
TriangularArray<T, D, P> const Grid<T, D, P, Q>::pNorm = []() {
	TriangularArray<T, D, P> norm;
	for (TriangularIndices<D, P> p3; p3 != p3.end(); p3++) {
		norm[p3] = one;
		for (int dim = 0; dim < D; dim++) {
			norm[p3] *= two * T(p3.indexAt(dim)) + one;
		}
	}
	return norm;
}();

#endif /* INCLUDE_GRID_HPP_ */
