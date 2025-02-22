/*
 * Numbers.hpp
 *
 *  Created on: Dec 12, 2024
 *      Author: dmarce1
 */

#ifndef INCLUDE_NUMBERS_HPP_
#define INCLUDE_NUMBERS_HPP_

#include <type_traits>
#include <vector>

#include "Vector.hpp"

namespace Math {

template<typename T>
constexpr T integerPower(T x, int n) {
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
		return one / integerPower(x, -n);
	}
}

template<typename T>
T invInteger(size_t n) {
	static constexpr T one(1);
	static thread_local std::vector<T> memory(1, -1);
	while (n >= memory.size()) {
		memory.push_back(one / T(memory.size()));
	}
	return memory[n];
}

template<typename T>
constexpr T kroneckerDelta(int n, int m) {
	static constexpr T zero = T(0);
	static constexpr T one = T(1);
	if (n == m) {
		return one;
	} else {
		return zero;
	}
}

template<typename T>
constexpr T negativeOne2Power(int k) {
	static constexpr T one = T(1);
	if (k & 1) {
		return -one;
	} else {
		return +one;
	}
}

template<typename T>
constexpr T nChooseK(int n, int k) {
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
constexpr T nFactorial(int n) {
	static constexpr T one = T(1);
	if (n <= 1) {
		return one;
	} else {
		return T(n) * nFactorial<T>(n - 1);
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

constexpr auto periodicModulus(auto a, auto b) {
	return a - b * (a / b - int(a < 0));
}

}
#endif /* INCLUDE_NUMBERS_HPP_ */
