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
#include "Vector.hpp"

#include <functional>

template<int N>
struct QuadratureRules {
	Math::Vector<Real, N> x;
	Math::Vector<Real, N> w;
};

template<int N>
QuadratureRules<N> gaussLobattoRules() {
	Real const one(1);
	Real const two(2);
	QuadratureRules<N> rules;
	Legendre::Polynomials<N - 1> legendreP;
	auto const Pn = legendreP(N - 1);
	auto const dPndx = polynomialDerivative(Pn);
	auto const xRoots = Math::polynomialFindAllRoots(dPndx);
	rules.x[0] = -one;
	rules.w[0] = one;
	for (int n = 0; n < N - 2; n++) {
		Real const xi = xRoots[n].real();
		rules.x[n + 1] = xi;
		rules.w[n + 1] = one / sqr(Pn(xi));
	}
	rules.x[N - 1] = one;
	rules.w[N - 1] = one;
	rules.w *= two / Real(N * (N - 1));
	return rules;
}

template<int O>
struct Basis {
	Basis(Real dx_) :
			dx(dx_) {
		Real const two(2.0);
		Real const one(1.0);
		Real const half(0.5);
		Legendre::Polynomials<O> legendreP;
		Math::Polynomial<Real> x;
		x[1] = Real(1);
		for (int n = 0; n < O; n++) {
			Real const norm = pow(dx, n) / Real(Math::binomial(2 * n, n));
			Pn[n] = legendreP(n);
			Pn[n] = norm * Pn[n](two * x / dx);
			dPndx[n] = Math::polynomialDerivative(Pn[n]);
			a[n] = pow(dx, n + 1) / polynomialIntegrate(Pn[n] * Pn[n], -half * dx, half * dx);
			pt[n] = half * dx * legendreP.point(n);
			wt[n] = half * dx * legendreP.weight(n);
		}
	}
	Real function(int n, Real x) const {
		return Pn[n](x);
	}
	Real derivative(int n, Real x) const {
		return dPndx[n](x);
	}
	Real coeff(int n) const {
		return a[n];
	}
	Real qPoint(int n) const {
		return pt[n];
	}
	Real qWeight(int n) const {
		return wt[n];
	}
	Polynomial<Real> function(int n) const {
		return Pn[n];
	}
	Vector<Real, O> const& coeffs() const {
		return a;
	}
private:
	Vector<Math::Polynomial<Real>, O> Pn;
	Vector<Math::Polynomial<Real>, O> dPndx;
	Vector<Real, O> a;
	Vector<Real, O> pt;
	Vector<Real, O> wt;
	Real dx;
};

#endif /* INCLUDE_BASIS_HPP_ */
