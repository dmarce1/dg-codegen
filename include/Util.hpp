#pragma once

#include <array>

#include "Vector.hpp"

void enableFPE();
void disableFPE();

template<int N, typename T>
inline constexpr auto repeat(T const &value) {
	std::array<T, N> a;
	for (int n = 0; n < N; n++) {
		a[n] = value;
	}
	return a;
}

template<int N, typename T>
inline constexpr auto insert(T const &value, int i, std::array<T, N - 1> const &A0) {
	std::array<T, N> A1;
	std::copy(A0.begin(), &A0[i], A1.begin());
	A1[i] = value;
	std::copy(&A0[i], A0.end(), &A1[i + 1]);
	return A1;
}

template<typename T>
inline constexpr T ipow(T x, int n) {
	static constexpr T one = T(1);
	if (n >= 0) {
		T xm = x;
		T xn = one;
		while (n) {
			if (n & 1) {
				xn *= xm;
			}
			xm *= xm;
			n >>= 1;
		}
		return xn;
	} else {
		return one / ipow(x, -n);
	}
}

template<typename T>
inline constexpr T binco(T n, T k) {
	static constexpr T one = T(1);
	T num = one;
	T den = one;
	for (int i = 1; i <= k; i++) {
		num *= T(n + 1 - i);
		den *= T(i);
	}
	return num / den;
}

template<typename T>
inline constexpr T factorial(int n) {
	static constexpr T one = T(1);
	if (n <= 1) {
		return one;
	} else {
		return T(n) * nFactorial < T > (n - 1);
	}
}

template<typename T>
inline constexpr T sqr(T r) {
	return r * r;
}

template<typename T>
inline constexpr T sign(T number) {
	static constexpr T zero = T(0);
	static constexpr T one = T(1);
	if (number > zero) {
		return +one;
	} else if (number < zero) {
		return -one;
	} else {
		return zero;
	}
}

inline constexpr int nonepow(int k) {
	return 1 - 2 * (k & 1);
}

template<typename T, int D>
inline constexpr std::array<T, D> unit(int d) {
	std::array<T, D> u;
	u.fill(T(0));
	u[d] = T(1);
	return u;
}

