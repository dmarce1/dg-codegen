/******************************************************************************
 Copyright (C) 2024  Dominic C. Marcello
 *******************************************************************************/

#ifndef INCLUDE_QUADRATURE_HPP_
#define INCLUDE_QUADRATURE_HPP_

#include <cmath>
#include <vector>
#include <unordered_map>
#include "Util.hpp"
#include "Matrix.hpp"
#include "Polynomial.hpp"

template<typename T>
struct QuadraturePoint {
	T position;
	T weight;
};

template<typename T>
std::vector<QuadraturePoint<T>> const& gaussLegendreQuadraturePoints(int l) {
	using namespace std;
	static const T π = T(4) * atan(T(1));
	static const T zero = 0;
	static const T one = 1;
	static const T two = 2;
	static const T half = one / two;
	static unordered_map<int, std::shared_ptr<vector<QuadraturePoint<T>>>> memory;
	if (memory.find(l) == memory.end()) {
		vector<QuadraturePoint<T>> results;
		results.reserve(l);
		for (int pointIndex = 0; pointIndex < l; ++pointIndex) {
			if ((l % 2 == 0) && (2 * pointIndex + 1 == l)) {
				T const w = two / sqr(assoc_legendre(l, 1, one));
				results.push_back(QuadraturePoint<T>( { zero, T(w) }));
				continue;
			}
			T θ = π * (T(pointIndex) + half) / T(l);
			T θ0;
			int cnt = 0;
			do {
				θ0 = θ;
				θ += legendre(l, cos(θ)) / assoc_legendre(l, 1, cos(θ));
				cnt++;
				if (cnt > 32) {
					printf("Failed to converge - %lli\n", ulpDistance(cos(θ0), cos(θ)));
					abort();
				}
			} while (abs(ulpDistance(cos(θ0), cos(θ))) > 64);
			T const x = cos(θ);
			T const w = two / sqr(assoc_legendre(l, 1, cos(θ)));
			results.push_back(QuadraturePoint<T>( { T(x), T(w) }));
		}
		memory[l] = std::make_shared<std::vector<QuadraturePoint<T>>>(std::move(results));
	}
	return *(memory[l]);
}

template<typename T>
std::vector<QuadraturePoint<T>> const& gaussLaguerreQuadraturePoints(int n) {
//	n=4;
//	α=1;
	using namespace std;
	static const T zero = 0;
	static const T one = 1;
	static const T two = 2;
	static const T four = 4;
	static const T half = one / two;
	static unordered_map<int, shared_ptr<vector<QuadraturePoint<T>>>> memory;
	if (memory.find(n) == memory.end()) {
		std::vector<QuadraturePoint<T>> results;
		std::vector<T> xbnd;
		for (int m = 2;; m *= 2) {
			xbnd.resize(1);
			xbnd[0] = zero;
			double dx = (four * n + two) / (n * m);
			double x = zero;
			double Ln = laguerre(n, zero);
			int count = 0;
			for (int k = 0; k < m * n; k++) {
				x += dx;
				double Lnp1 = laguerre(n, x);
				if (copysign(one, Lnp1) * copysign(one, Ln) <= zero) {
					count++;
					xbnd.push_back(x);
				}
				Ln = Lnp1;
			}
			if (n == count) {
				break;
			}
		}
		for (int k = 0; k < n; k++) {
			T xmin = xbnd[k];
			T xmax = xbnd[k + 1];
			T x, w;
			do {
				x = half * (xmax + xmin);
				if (copysign(one, laguerre(n, x)) * copysign(one, laguerre(n, xmin)) < zero) {
					xmax = x;
				} else {
					xmin = x;
				}
			} while (xmax > nextafter(xmin, xmax));
			w = x / (sqr((n + 1) * laguerre(n + 1, x)));
//			if(n==4) printf( "%e %e\n", x,w);
			results.push_back( { x, w });
		}
//		if(n==4) abort();
		memory[n] = make_shared<vector<QuadraturePoint<T>>>(move(results));
	}
	return *(memory[n]);
}

#endif
