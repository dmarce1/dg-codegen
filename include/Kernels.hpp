/*
 * Kernels.hpp
 *
 *  Created on: Feb 22, 2025
 *      Author: dmarce1
 */

#ifndef INCLUDE_KERNELS_HPP_
#define INCLUDE_KERNELS_HPP_

#include "Numbers.hpp"
#include "Utilities.hpp"

template<typename R>
auto kernelTSC(R x) {
	x = abs(x);
	if (x < R(0.5)) {
		return R(0.75) - x * x;
	} else if (x < R(1.5)) {
		return (x * (x * R(0.5) - R(1.5)) + R(1.125));
	} else {
		return R(0.0);
	}
}

#endif /* INCLUDE_KERNELS_HPP_ */
