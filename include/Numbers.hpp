/*
 * Numbers.hpp
 *
 *  Created on: Dec 12, 2024
 *      Author: dmarce1
 */

#ifndef INCLUDE_NUMBERS_HPP_
#define INCLUDE_NUMBERS_HPP_

#include <type_traits>

#include "Vector.hpp"

namespace Math {

using Integer = int;
using Index = unsigned int;

constexpr auto sign(auto number) {
	using value_type = decltype(number);
	if (number > value_type(0)) {
		return value_type(1);
	}
	if (number < value_type(0)) {
		return -value_type(1);
	}
	return value_type(0);
}

constexpr Integer nOnePow(Integer k) {
	return 1 - 2 * (k & 1);
}

constexpr auto sqr(auto r) {
	return r * r;
}

template<typename RType, typename IType>
constexpr typename std::enable_if<std::is_integral_v<IType>, RType>::type Pow(RType const &x, IType n) {
	RType xm = x;
	RType xn = RType(1);
	while (n) {
		if (n & 1) {
			xn *= xm;
		}
		xm *= xm;
		n >>= 1;
	}
	return xn;
}

constexpr int factorial(int n) {
	if (n == 0) {
		return 1;
	} else {
		return n * factorial(n - 1);
	}
}

inline auto binomial(auto n, auto k) {
	using T = std::remove_reference_t<decltype(n)>;
	T result = T(1);
	for (int j = n - k + 1; j <= n; j++) {
		result *= T(j);
	}
	for (int j = 2; j <= k; j++) {
		result /= T(j);
	}
	return result;
}

inline auto minmod(auto const &a, auto const &b) {
	using type = std::remove_reference_t<decltype(a)>;
	type const half(0.5);
	return half * (sign(a) + sign(b)) * min(abs(a), abs(b));
}

inline auto maxmod(auto const &a, auto const &b) {
	using type = std::remove_reference_t<decltype(a)>;
	type const half(0.5);
	return half * (sign(a) + sign(b)) * max(abs(a), abs(b));
}

inline auto minmodTheta(auto const &a, auto const &b, auto const &theta) {
	using type = std::remove_reference_t<decltype(a)>;
	type const half(0.5);
	const auto ab = half * (a + b);
	return minmod(theta * minmod(a, b), ab);
}

inline auto superBee(auto const &a, auto const &b) {
	using T = std::remove_reference_t<decltype(a)>;
	return maxmod(minmod(T(2) * a, b), minmod(a, T(2) * b));
}

inline auto vanLeer(auto const &a, auto const &b) {
	using T = std::remove_reference_t<decltype(a)>;
	T const tiny(1e-20);
	T const ab = a * b;
	return (ab + abs(ab)) / (b + a + tiny);
}

template<typename T, int N>
inline Vector<T, N> vanLeer(Vector<T, N> const &a, Vector<T, N> const &b) {
	Vector<T, N> c;
	for (int n = 0; n < N; n++) {
		c[n] = vanLeer(a[n], b[n]);
	}
	return c;
}

template<typename T, int N>
inline Vector<T, N> superBee(Vector<T, N> const &a, Vector<T, N> const &b) {
	Vector<T, N> c;
	for (int n = 0; n < N; n++) {
		c[n] = superBee(a[n], b[n]);
	}
	return c;
}

template<typename T, int N>
inline Vector<T, N> minmod(Vector<T, N> const &a, Vector<T, N> const &b) {
	Vector<T, N> c;
	for (int n = 0; n < N; n++) {
		c[n] = minmod(a[n], b[n]);
	}
	return c;
}

template<typename T, int N>
inline Vector<T, N> minmodTheta(Vector<T, N> const &a, Vector<T, N> const &b, T const &theta) {
	Vector<T, N> c;
	T const half(0.5);
	for (int n = 0; n < N; n++) {
		c[n] = minmodTheta(a[n], b[n], theta);
	}
	return c;
}

}
#endif /* INCLUDE_NUMBERS_HPP_ */
