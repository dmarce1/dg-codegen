/*
 * Integrate.hpp
 *
 *  Created on: Dec 16, 2024
 *      Author: dmarce1
 */

#ifndef INCLUDE_INTEGRATE_HPP_
#define INCLUDE_INTEGRATE_HPP_

#include <cmath>
#include <functional>
#include <unordered_map>

#include "Numbers.hpp"

namespace Math {

template<typename T>
T integrate(std::function<T(T)> const &f, T const &a, T const &b) {
	constexpr T toler = 1000.0 * std::numeric_limits<T>::epsilon();
	static std::vector<T> Rnp1(1);
	static std::vector<T> Rn(1);
	auto const tIntegrate = [f, a, b](int nPow, T dH) {
		T thisSum = 0.0;
		if (nPow == 0) {
			thisSum = (f(b) + f(a)) * dH;
		} else {
			for (int k = 1; k <= nPow; k++) {
				T const x = a + (2 * k - 1) * dH;
				thisSum += f(x) * dH;
			}
		}
		return thisSum;
	};
	Rn[0] = tIntegrate(0, T(0.5) * T(b - a));
	int nPow = 1;
	int n = 1;
	T dH = T(0.5) * T(b - a);
	T err;
	do {
		while (int(n * n + n) >= 2 * int(Rnp1.size())) {
			Rnp1.resize(2 * Rnp1.size());
			Rn.resize(Rnp1.size());
		}
		Rnp1[0] = T(0.5) * Rn[0] + tIntegrate(nPow, dH);
		int wt = 4;
		for (int m = 1; m <= n; m++) {
			Rnp1[m] = (T(wt) * Rnp1[m - 1] - Rn[m - 1]) / T(wt - 1);
			wt <<= 2;
		}
		err = abs(Rnp1[n] / Rn[n - 1] - 1.0);
		std::swap(Rn, Rnp1);
		nPow <<= 1;
		n++;
		dH *= T(0.5);
	} while (err > toler);
	return Rn[n];
}

}
#endif /* INCLUDE_INTEGRATE_HPP_ */
