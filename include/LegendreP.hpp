/*
 * LegendreP.hpp
 *
 *  Created on: Dec 15, 2024
 *      Author: dmarce1
 */

#ifndef INCLUDE_LEGENDREP_HPP_
#define INCLUDE_LEGENDREP_HPP_

#include <hpx/hpx.hpp>
#include "Real.hpp"
#include "Polynomial.hpp"

namespace Legendre {

template<int Order>
struct Polynomials {
	Polynomials() {
		Math::Polynomial<Real> one, x;
		one[0] = 1.0_Real;
		x[1] = 1.0_Real;
		auto pnm1 = one;
		auto pn = x;
		legendreP.push_back(one);
		legendreP.push_back(x);
		for (int n = 1; n <= Order; n++) {
			auto const aCons = Math::Polynomial<Real>(Real(2 * n + 1) / Real(n + 1));
			auto const bCons = Math::Polynomial<Real>(Real(n) / Real(n + 1));
			auto const pnp1 = aCons * pn * x - bCons * pnm1;
			pnm1 = std::move(pn);
			pn = std::move(pnp1);
			legendreP.push_back(pn);
		}

		dLegendreP.push_back(Math::polynomialDerivative(legendreP[0]));
		for (int n = 1; n <= Order; n++) {
			dLegendreP.push_back(Math::polynomialDerivative(legendreP[n]));
		}
		auto const roots = Math::polynomialFindAllRoots(legendreP[Order]);
		for (auto const &pt : roots) {
			quadPos.push_back(pt.real());
		}
		std::sort(quadPos.begin(), quadPos.end());
		auto const dPdx = dLegendreP[Order];
		for (int p = 0; p < Order; p++) {
			auto const x = quadPos[p];
			auto const x2 = sqr(x);
			quadWt.push_back(Real(2) / (sqr(dPdx(x)) * (Real(1) - x2)));
		}
	}
	Real valueAt(int n, Real x) const {
		return legendreP[n](x);
	}
	Real derivativeAt(int n, Real x) const {
		return dLegendreP[n](x);
	}
	Real point(int n) const {
		return quadPos[n];
	}
	Real weight(int n) const {
		return quadWt[n];
	}
	Math::Polynomial<Real> operator()( int n) const {
		return legendreP[n];
	}
private:
	std::vector<Math::Polynomial<Real>> legendreP;
	std::vector<Math::Polynomial<Real>> dLegendreP;
	std::vector<Real> quadPos;
	std::vector<Real> quadWt;
}
;

}

#endif /* INCLUDE_LEGENDREP_HPP_ */
