/*
 * LegendreP.hpp
 *
 *  Created on: Dec 15, 2024
 *      Author: dmarce1
 */

#ifndef INCLUDE_LEGENDREP_HPP_
#define INCLUDE_LEGENDREP_HPP_

#include "Real.hpp"
#include "Polynomial.hpp"

Real legendreP(int n, Real x) {
	Real const one(1.0);
	if (n == 0) {
		return one;
	} else if (n == 1) {
		return x;
	} else {
		Real Pnm1, Pn;
		Pnm1 = one;
		Pn = x;
		for (int l = 1; l < n; l++) {
			Real const a = Real(2 * l + 1) / Real(l + 1);
			Real const b = Real(l) / Real(l + 1);
			Real const Pnp1 = a * Pn * x - b * Pnm1;
			Pnm1 = Pn;
			Pn = Pnp1;
		}
		return Pn;
	}
}

Real dLegendrePdX(int n, Real x) {
	Real const zero(0), one(1), three(3);
	if (n == 0) {
		return zero;
	} else if (n == 1) {
		return one;
	} else {
		Real Pnm1, Pn;
		Real dPnm1dX, dPndX;
		Pnm1 = one;
		Pn = x;
		dPndX = one;
		for (int l = 1; l < n - 1; l++) {
			Real const a = Real(l + 1);
			Real const ainv = one / a;
			Real const b = Real(2 * l + 1) * ainv;
			Real const c = Real(l) * ainv;
			Real const Pnp1 = b * Pn * x - c * Pnm1;
			Real const dPnp1dX = a * Pn + x * dPndX;
			Pnm1 = Pn;
			Pn = Pnp1;
			dPndX = dPnp1dX;
		}
		return Real(n) * Pn + x * dPndX;
	}
}

Real d2LegendrePdX2(int n, Real x) {
	Real const zero(0), one(1), three(3);
	if (n <= 1) {
		return zero;
	} else if (n == 2) {
		return three;
	} else {
		Real Pnm1, Pn;
		Real dPnm1dX, dPndX;
		Real d2Pnm1dX2, d2PndX2;
		Pnm1 = one;
		Pn = x;
		dPndX = one;
		d2PndX2 = three;
		for (int l = 1; l < n - 2; l++) {
			Real const a = Real(l + 1);
			Real const ainv = one / a;
			Real const b = Real(2 * l + 1) * ainv;
			Real const c = Real(l) * ainv;
			Real const d = Real(l + 2);
			Real const Pnp1 = b * Pn * x - c * Pnm1;
			Real const dPnp1dX = a * Pn + x * dPndX;
			Real const d2Pnp1dX2 = d * dPndX + x * d2PndX2;
			Pnm1 = Pn;
			Pn = Pnp1;
			dPndX = dPnp1dX;
			d2PndX2 = d2Pnp1dX2;
		}
		Real const dPnp1dX = Real(n - 1) * Pn + x * dPndX;
		d2PndX2 = Real(n) * dPndX + x * d2PndX2;
		dPndX = dPnp1dX;
		return Real(n + 1) * dPndX + x * d2PndX2;
	}
}

#endif /* INCLUDE_LEGENDREP_HPP_ */
