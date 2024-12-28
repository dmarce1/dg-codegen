/*
 * Quadrature.hpp
 *
 *  Created on: Dec 27, 2024
 *      Author: dmarce1
 */

#ifndef INCLUDE_BASIS_HPP_
#define INCLUDE_BASIS_HPP_

#include "LegendreP.hpp"
#include "Numbers.hpp"
#include "Real.hpp"

#include <functional>

template<int P>
struct Basis {
	Basis() {
		Legendre::Polynomials<P> polys;
		for (int n = 0; n < P; n++) {
			Math::Polynomial<Real> X;
			X[1] = Real(2);
			Pn[n] = polys(n)(X);
			dPndx[n] = Math::polynomialDerivative(Pn[n]);
			auto const Pn2 = Pn[n] * Pn[n];
			M[n] = Math::polynomialIntegrate(Pn2, -Real(0.5), +Real(0.5));
			pt[n] = Real(0.5) * polys.point(n);
			wt[n] = Real(0.5) * polys.weight(n);
		}
	}
	Real function(int n, Real x) const {
		return Pn[n](x);
	}
	Polynomial<Real> function(int n) const {
		return Pn[n];
	}
	Real derivative(int n, Real x) const {
		return dPndx[n](x);
	}
	Real boundaryValue(int n, int sgn) const {
		return Pn[n](Real(sgn));
	}
	Real qPoint(int n) const {
		return pt[n];
	}
	Real qWeight(int n) const {
		return wt[n];
	}
	Real mass(int n) const {
		return M[n];
	}
	Real massInv(int n) const {
		return Real(1) / M[n];
	}
private:
	Vector<Math::Polynomial<Real>, P> Pn;
	Vector<Math::Polynomial<Real>, P> dPndx;
	Vector<Real, P> M;
	Vector<Real, P> pt;
	Vector<Real, P> wt;
};

#endif /* INCLUDE_BASIS_HPP_ */
