#ifndef INCLUDE_QUADRATURE_HPP_
#define INCLUDE_QUADRATURE_HPP_

#include "LegendreP.hpp"
#include "Numbers.hpp"
#include "Real.hpp"
#include "Vector.hpp"

#include <functional>

template<int N>
struct QuadratureRules {
	Math::Vector<Real, N> x;
	Math::Vector<Real, N> w;
	bool withEndpoints;
	Real integrate(std::function<Real(Real)> const &f, Real dx) const {
		Real I = Real(0);
		for (int n = 0; n < N; n++) {
			I += f(Real(0.5) * x[n] * dx) * w[n] * Real(0.5) * dx;
		}
		return I;
	}
	std::string toString() const {
		std::string str;
		for (int n = 0; n < N; n++) {
			str += std::string("(") + to_string(x[n]);
			str += std::string(", ") + to_string(w[n]) + std::string(")");
			if (n != N - 1) {
				str += std::string(", ");
			}
		}
		return str;
	}
};

inline Polynomial<Real> openLagrangePolynomial(int n, int j, Real a, Real b) {
	Real const one(1);
	Real const half(0.5);
	Real const zero(0);
	Polynomial < Real > p(one);
	Polynomial < Real > X;
	X[1] = one;
	Real h = (b - a) / Real(n);
	Real xj = (Real(j) + half) * h;
	for (int m = 0; m <= n; m++) {
		if (m != j) {
			Real xm = Real(m) * h + a;
			Polynomial<Real> const X0(xm);
			Polynomial<Real> const dX = X - X0;
			auto const inv = one / (xj - xm);
			p = inv * p * dX;
		}
	}
	return p;
}

inline Polynomial<Real> closedLagrangePolynomial(int n, int j, Real a, Real b) {
	Real const one(1);
	Real const zero(0);
	Polynomial < Real > p(one);
	Polynomial < Real > X;
	X[1] = one;
	Real h = (b - a) / Real(n - 1);
	Real xj = Real(j) * h + a;
	for (int m = 0; m < n; m++) {
		if (m != j) {
			Real xm = Real(m) * h + a;
			Polynomial<Real> const X0(xm);
			Polynomial<Real> const dX = X - X0;
			auto const inv = one / (xj - xm);
			p = inv * p * dX;
		}
	}
	return p;
}

template<int N>
QuadratureRules<N> closedNewtonCotesRules(Real a = Real(-1), Real b = Real(1)) {
	Real const half(0.5);
	Real const one(1);
	Real const two(2);
	Real const h = (b - a) / Real(N - 1);
	QuadratureRules<N> rules;
	for (int i = 0; i < N; i++) {
		rules.x[i] = a + Real(i) * h;
		auto const pl = closedLagrangePolynomial(N, i, a, b);
		rules.w[i] = polynomialIntegrate(pl, a, b);
	}
	rules.withEndpoints = true;
	return rules;
}

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
	rules.withEndpoints = true;
	return rules;
}

template<int N>
QuadratureRules<N> gaussLegendreRules() {
	Real const half(0.5), one(1.0), two(2.0);
	QuadratureRules<N> rules;
	Math::Polynomial<Real> Pnm1, x;
	x[1] = one;
	Pnm1[0] = one;
	auto Pn = x;
	for (int n = 1; n < N; n++) {
		auto const aCons = Math::Polynomial<Real>(Real(2 * n + 1) / Real(n + 1));
		auto const bCons = Math::Polynomial<Real>(Real(n) / Real(n + 1));
		auto const Pnp1 = aCons * Pn * x - bCons * Pnm1;
		Pnm1 = std::move(Pn);
		Pn = std::move(Pnp1);
	}
	auto const dPn_dx = Math::polynomialDerivative(Pn);
	auto const roots = Math::polynomialFindAllRoots(Pn);
	for (int n = 0; n < N; n++) {
		rules.x[n] = roots[n].real();
	}
	std::sort(rules.x.begin(), rules.x.end());
	for (int p = 0; p < N; p++) {
		auto const x = rules.x[p];
		auto const x2 = sqr(x);
		rules.w[p] = two / (sqr(dPn_dx(x)) * (one - x2));
	}
	rules.withEndpoints = false;
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
			Real const norm = sqrt(Real(2 * n + 1));
			Pn[n] = legendreP(n);
			Pn[n] = norm * Pn[n](two * x / dx);
			dPndx[n] = Math::polynomialDerivative(Pn[n]);
			a[n] = dx / polynomialIntegrate(Pn[n] * Pn[n], -half * dx, half * dx);
			pt[n] = half * dx * legendreP.point(n);
			wt[n] = half * dx * legendreP.weight(n);
		}
	}
	std::vector<Real> transform(std::function<Real(Real)> const &func) const {
		std::vector<Real> u(O);
		for (int l = 0; l < O; l++) {
			u[l] = Real(0);
			for (int q = 0; q < O; q++) {
				u[l] += wt[q] * Pn[l](pt[q]) * func(pt[q]);
			}
		}
		return u;
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

#endif
