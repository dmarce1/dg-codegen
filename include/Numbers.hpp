/*
 * Numbers.hpp
 *
 *  Created on: Dec 12, 2024
 *      Author: dmarce1
 */

#ifndef INCLUDE_NUMBERS_HPP_
#define INCLUDE_NUMBERS_HPP_

#include <type_traits>

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
constexpr typename std::enable_if<std::is_integral_v<RType>, IType>::type Pow(RType const &x, IType const &n) {
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

template<typename Type>
Type binomial(int n, int k) {
	Type result = Type(1);
	for (int j = n - k + 1; j <= n; j++) {
		result *= Type(j);
	}
	for (int j = 2; j <= k; j++) {
		result /= Type(j);
	}
	return result;
}

inline auto minmod(auto const &a, auto const &b) {
	using type = std::remove_reference_t<decltype(a)>;
	type const half(0.5);
	return half * (sign(a) + sign(b)) * min(abs(a), abs(b));
}

inline auto vanLeer(auto const &a, auto const &b) {
	using type = std::remove_reference_t<decltype(a)>;
	auto const ab = a * b;
	if( ab > type(0)) {
		return 	ab  / (a + b);
	} else {
		return type(0);
	}
}

inline auto minmod_theta(auto const &a, auto const &b, auto const &theta) {
	using type = std::remove_reference_t<decltype(a)>;
	type const half(0.5);
	const auto ab = half * (a + b);
	return minmod(theta * minmod(a, b), ab);
}

}
#endif /* INCLUDE_NUMBERS_HPP_ */
