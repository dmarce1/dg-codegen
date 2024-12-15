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
	if( number > value_type(0)) {
		return value_type(1);
	}
	if( number < value_type(0)) {
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
constexpr std::enable_if<std::is_integral_v<RType>, IType>::type Pow(RType const &x, IType const &n) {
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

}

#endif /* INCLUDE_NUMBERS_HPP_ */
