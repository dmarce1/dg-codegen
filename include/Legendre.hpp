#pragma once

#include <array>
#include <limits>

template<typename T, int N, int M>
constexpr  auto dMLegendrePdXm(T const &x) {
	std::array<std::array<T, N>, M + 1> P { };
	P[0][0] = T(1);
	if constexpr (N >= 2) {
		P[0][1] = x;
		for (int n = 1; n < N - 1; ++n) {
			T const a = T(2 * n + 1) / T(n + 1);
			T const b = T(n) / T(n + 1);
			P[0][n + 1] = a * x * P[0][n] - b * P[0][n - 1];
		}
	}
	for (int m = 1; m <= M; ++m) {
		P[m][0] = T(0);
		if constexpr (N >= 2) {
			P[m][1] = m == 1 ? T(1) : T(0); // Derivatives of x
		}
		for (int n = 1; n < N - 1; ++n) {
			T const a = T(2 * n + 1) / T(n + 1);
			T const b = T(n) / T(n + 1);
			P[m][n + 1] = a * (x * P[m][n] + T(m) * P[m - 1][n]) - b * P[m][n - 1];
		}
	}
	return P;
}

template<typename T, int N>
constexpr  auto legendreP(T const &x) {
	auto const a = dMLegendrePdXm<T, N, 0>(x);
	return a[0];
}

