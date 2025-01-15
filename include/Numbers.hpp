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

template<typename T, typename U>
constexpr T integerPower(T x, U n) {
	static constexpr T one = T(1);
	static constexpr U iOne = U(1);
	T xm = x;
	T xn = one;
	while (n) {
		if (n & iOne) {
			xn *= xm;
		}
		xm *= xm;
		n >>= iOne;
	}
	return xn;
}

template<typename T>
constexpr T kroneckerDelta(T n, T m) {
	static constexpr T zero = T(0);
	static constexpr T one = T(1);
	if (n == m) {
		return one;
	} else {
		return zero;
	}
}

template<typename T>
constexpr T negativeOne2Power(T k) {
	static constexpr T one = T(1);
	if( k & 1 ) {
		return -one;
	} else {
		return +one;
	}
}

template<typename T>
constexpr T nChooseK(T n, T k) {
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
constexpr T nFactorial(T n) {
	static constexpr T one = T(1);
	if (n <= 1) {
		return one;
	} else {
		return n * nFactorial(n - 1);
	}
}

template<typename T>
constexpr T nSign(T number) {
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

template<typename T>
constexpr T nSquared(T r) {
	return r * r;
}


}
#endif /* INCLUDE_NUMBERS_HPP_ */
