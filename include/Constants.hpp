/*
 * Constants.hpp
 *
 *  Created on: Jan 12, 2025
 *      Author: dmarce1
 */

#ifndef INCLUDE_CONSTANTS_HPP_
#define INCLUDE_CONSTANTS_HPP_

#include "Real.hpp"

template<typename T>
struct Constants {
	constexpr Constants() = default;
	static constexpr T c = T(2.99792458e10);
	static constexpr T G = T(6.67259e-8);
	static constexpr T m = T(1.6605402e-24);
	static constexpr T kB = T(1.380658e-16);
	static constexpr T h = T(6.6260755e-27);
	static constexpr T aR = T(7.5646e-15);
	static constexpr T sigma = T(5.67051e-5);
	static constexpr T Msol = T(1.99e33);
	static constexpr T Rsol = T(6.96e10);
	static constexpr T Lsol = T(3.9e33);
	static constexpr T AU = T(1.496e13);
	static constexpr T pc = T(3.086e18);
	static constexpr T ly = T(9.463e17);
};

#endif /* INCLUDE_CONSTANTS_HPP_ */
