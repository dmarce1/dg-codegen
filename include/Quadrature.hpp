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

template<int N>
Math::Vector<Polynomial<Real>, N> lagrangePolynomials(Math::Vector<Real, N> const &X) {
	Real const one(1);
	Real const zero(0);
	Polynomial < Real > x1;
	x1[1] = one;
	Math::Vector<Polynomial<Real>, N> polys;
	for (int j = 0; j < N; j++) {
		Polynomial < Real > p(one);
		Real xj = X[j];
		for (int m = 0; m < N; m++) {
			if (m != j) {
				Real const xm = X[m];
				Polynomial<Real> const x0(xm);
				Polynomial<Real> const dx = x1 - x0;
				auto const inv = one / (xj - xm);
				p = inv * p * dx;
			}
		}
		polys[j] = p;
	}
	return polys;
}

template<int N>
QuadratureRules<N> closedNewtonCotesRules(Real a = Real(-1), Real b = Real(1)) {
	Real const h = (b - a) / Real(N - 1);
	QuadratureRules<N> rules;
	for (int i = 0; i < N; i++) {
		rules.x[i] = a + Real(i) * h;
	}
	auto const L = lagrangePolynomials<N>(rules.x);
	for (int i = 0; i < N; i++) {
		rules.w[i] = polynomialIntegrate(L[i], a, b);
	}
	rules.withEndpoints = true;
	return rules;
}

template<int N>
QuadratureRules<N> openNewtonCotesRules(Real a = Real(-1), Real b = Real(1)) {
	Real const h = (b - a) / Real(N + 1);
	QuadratureRules<N> rules;
	for (int i = 0; i < N; i++) {
		rules.x[i] = a + Real(i + 1) * h;
	}
	auto const L = lagrangePolynomials<N>(rules.x);
	for (int i = 0; i < N; i++) {
		rules.w[i] = polynomialIntegrate(L[i], a, b);
	}
	rules.withEndpoints = false;
	return rules;
}

template<int N>
QuadratureRules<N> gaussLegendreRules() {
	Real const half(0.5), one(1.0), two(2.0);
	Real const toler((1 << (2 + std::ilogb(N) / 3)) * std::numeric_limits<double>::epsilon());
	QuadratureRules<N> rules;
	Real const dx = one / Real(N);
	Real const pi = Real(M_PI);
	for (int l = 0; l < N; l++) {
		Real dTheta;
		Real theta = pi * (one - (Real(l) + half) * dx);
		do {
			Real const cosTheta = cos(theta);
			Real const sinTheta = sin(theta);
			Real const Pn = legendreP(N, cosTheta);
			Real const dPnDx = dLegendrePdX(N, cosTheta);
			Real const iDen = one / (sinTheta * dPnDx);
			dTheta = Pn * iDen;
			theta += dTheta;
		} while (abs(dTheta) > toler);
		Real const x = cos(theta);
		rules.x[l] = x;
		rules.w[l] = two / (sqr(dLegendrePdX(N, x)) * (one - sqr(x)));
	}
	rules.withEndpoints = false;
	return rules;
}

template<int P>
QuadratureRules<P> clenshawCurtisRules() {
	constexpr int N = P - 1;
	Real const zero(0), half(0.5), one(1), two(2), pi(M_PI);
	Real const Ninv = one / Real(N);
	Real const piNinv = pi * Ninv;
	QuadratureRules<P> rules;
	Vector < Real, N / 2 + 1 > d;
	d[0] = one;
	for (int k = 1; k <= N / 2; k++) {
		d[k] = two / Real(1 - sqr(2 * k));
	}
	for (int n = 0; n <= N / 2; n++) {
		Real const c = (Real(n != 0) + one);
		Real w = zero;
		for (int k = 0; k <= N / 2; k++) {
			w += c * d[k] * cos(Real(2 * n * k) * piNinv) * Ninv;
		}
		rules.w[n] = w;
	}
	for (int k = 0; k < N - k; k++) {
		rules.w[N - k] = rules.w[k];
	}
	for (int k = 0; k <= N; k++) {
		rules.x[k] = cos(pi * Real(k) / Real(N));
	}
	return rules;
}

template<int N>
QuadratureRules<N> gaussLobattoRules() {
	Real const half(0.5), one(1.0), two(2.0);
	Real const toler((3 << (1 + (std::ilogb(N)) / 3)) * std::numeric_limits<double>::epsilon());
	QuadratureRules<N> rules;
	Real const dx = one / Real(N);
	Real const pi = Real(M_PI);
	for (int l = 1; l < N - 1; l++) {
		Real dTheta;
		Real theta = pi * (one - (Real(l) + half) * dx);
		do {
			Real const cosTheta = cos(theta);
			Real const sinTheta = sin(theta);
			Real const cotTheta = cosTheta / sinTheta;
			Real const dPnDx = dLegendrePdX(N - 1, cosTheta);
			Real const d2PnDx2 = d2LegendrePdX2(N - 1, cosTheta);
			Real const iDen = one / (sinTheta * d2PnDx2 - dPnDx * cotTheta);
			dTheta = dPnDx * iDen;
			theta += dTheta;
		} while (abs(dTheta) > toler);
		Real const x = cos(theta);
		rules.x[l] = x;
		rules.w[l] = two / (Real(N * (N - 1)) * sqr(legendreP(N - 1, x)));
	}
	rules.x[0] = -one;
	rules.x[N - 1] = one;
	rules.w[0] = rules.w[N - 1] = two / (Real(N * (N - 1)));
	rules.withEndpoints = true;
	return rules;
}

#endif
