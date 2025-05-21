#pragma once

#include <array>
#include <limits>

template<typename T, int N, int M>
constexpr auto dMLegendrePdXm(T const &x) {
	std::array<std::array<T, N>, M + 1> P;
	static constexpr auto nInv = []() {
		std::array<T, N> nInv;
		nInv[0] = std::numeric_limits<T>::infinity();
		for (int n = 1; n < N; n++) {
			nInv[n] = T(1) / T(n);
		}
		return nInv;
	}();
	if constexpr (N >= 1) {
		P[0][0] = T(1);
	}
	if constexpr (N >= 2) {
		P[0][1] = x;
	}
	if constexpr (N >= 3) {
		for (int n = 1; n < N - 1; n++) {
			P[0][n + 1] = (T(2 * n + 1) * P[0][n] * x - T(n) * P[0][n - 1]) * nInv[n + 1];
		}
	}
	for (int m = 1; m <= M; m++) {
		if constexpr (N >= 1) {
			P[m][0] = T(0);
		}
		if constexpr (N >= 2) {
			for (int n = 0; n < N - 1; n++) {
				P[m][n + 1] = T(n + m) * P[m - 1][n] + x * P[m][n];
			}
		}
	}
	return P;
}

template<typename T, int N>
constexpr auto legendreP(T const &x) {
	auto const a = dMLegendrePdXm<T, N, 0>(x);
	return a[0];
}

