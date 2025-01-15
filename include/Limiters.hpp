/*
 * Limiters.hpp
 *
 *  Created on: Jan 14, 2025
 *      Author: dmarce1
 */

#ifndef INCLUDE_LIMITERS_HPP_
#define INCLUDE_LIMITERS_HPP_

#include "Numbers.hpp"

template<typename T>
T minmod(T a, T b) {
	using namespace Math;
	static constexpr T half = T(0.5);
	return half * (nSign(a) + nSign(b)) * min(abs(a), abs(b));
}

template<typename T>
T monotonizedCentral(T a, T b) {
	static constexpr T half = T(0.5), two = T(2.0);
	return minmod(two * minmod(a, b), half * (a + b));
}

template<typename T>
T superBee(T a, T b) {
	using namespace Math;
	static constexpr T half = T(0.5), two = T(2.0);
	T const c = minmod(two * a, b);
	T const d = minmod(two * b, a);
	return half * (nSign(c) + nSign(d)) * max(abs(c), abs(d));
}

template<typename T>
T vanLeer(T a, T b) {
	static constexpr T zero = T(0.0), two = T(2.0);
	T const num = a * b;
	if (num > zero) {
		T const den = a + b;
		return two * num / den;
	} else {
		return zero;
	}
}

template<typename T>
T ospreLimiter(T a, T b) {
	using namespace Math;
	static constexpr T zero = T(0), two = T(2), three = T(3);
	T const ab = a * b;
	if (ab > zero) {
		T const a2 = nSquared(a);
		T const b2 = nSquared(b);
		T const num = three * (a2 * b + a * b2);
		T const den = two * (a2 + ab + b2);
		return num / den;
	} else {
		return zero;
	}
}

template<typename T>
T symmetricUMIST(T a, T b) {
	static constexpr T two = T(2.0), oneQuarter = T(0.25), threeQuarter = T(0.75);
	T const cLim1 = two * a;
	T const cLim2 = threeQuarter * a + oneQuarter * b;
	T const cLim3 = oneQuarter * a + threeQuarter * b;
	T const cLim4 = two * b;
	return minmod(minmod(cLim1, cLim2), minmod(cLim3, cLim4));
}

#endif /* INCLUDE_LIMITERS_HPP_ */
