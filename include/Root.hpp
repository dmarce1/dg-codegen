/*
 * Root.hpp
 *
 *  Created on: Jan 8, 2025
 *      Author: dmarce1
 */

#ifndef INCLUDE_ROOT_HPP_
#define INCLUDE_ROOT_HPP_

#include "Real.hpp"

#include <functional>

inline Real findRootNewtonRhapson(std::function<Real(Real)> const &F, std::function<Real(Real)> const &dFdx,
		Real const x0 = Real(0), Real const toler = Real(1e-10)) {
	Real x = x0;
	Real f;
	do {
		f = F(x);
		Real const dfdx = dFdx(x);
		Real const dx = -f / dfdx;
		x += dx;
	} while (abs(f) > toler);
	return x;
}

#endif /* INCLUDE_ROOT_HPP_ */
