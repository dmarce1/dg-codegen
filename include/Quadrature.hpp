#ifndef INCLUDE_QUADRATURE_HPP_
#define INCLUDE_QUADRATURE_HPP_

#include "LegendreP.hpp"
#include "Numbers.hpp"
#include "Polynomial.hpp"
#include "Vector.hpp"

#include <functional>

namespace Quadrature {

using namespace Math;

template<typename T, int N>
struct QuadratureRules {
	Vector<T, N> x;
	Vector<T, N> w;
	bool withEndpoints;
	T integrate(std::function<T(T)> const &f, T dx) const {
		T I = T(0);
		for (int n = 0; n < N; n++) {
			I += f(T(0.5) * x[n] * dx) * w[n] * T(0.5) * dx;
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

template<typename T, int N>
Vector<Polynomial<T>, N> lagrangePolynomials(Vector<T, N> const &X) {
	T const one(1);
	T const zero(0);
	Polynomial<T> x1;
	x1[1] = one;
	Vector<Polynomial<T>, N> polys;
	for (int j = 0; j < N; j++) {
		Polynomial<T> p(one);
		T xj = X[j];
		for (int m = 0; m < N; m++) {
			if (m != j) {
				T const xm = X[m];
				Polynomial<T> const x0(xm);
				Polynomial<T> const dx = x1 - x0;
				auto const inv = one / (xj - xm);
				p = inv * p * dx;
			}
		}
		polys[j] = p;
	}
	return polys;
}

template<typename T, int N>
QuadratureRules<T, N> closedNewtonCotesRules(T a = T(-1), T b = T(1)) {
	T const h = (b - a) / T(N - 1);
	QuadratureRules<T, N> rules;
	for (int i = 0; i < N; i++) {
		rules.x[i] = a + T(i) * h;
	}
	auto const L = lagrangePolynomials<N>(rules.x);
	for (int i = 0; i < N; i++) {
		rules.w[i] = polynomialIntegrate(L[i], a, b);
	}
	rules.withEndpoints = true;
	return rules;
}

template<typename T, int N>
QuadratureRules<T, N> openNewtonCotesRules(T a = T(-1), T b = T(1)) {
	T const h = (b - a) / T(N + 1);
	QuadratureRules<T, N> rules;
	for (int i = 0; i < N; i++) {
		rules.x[i] = a + T(i + 1) * h;
	}
	auto const L = lagrangePolynomials<T, N>(rules.x);
	for (int i = 0; i < N; i++) {
		rules.w[i] = polynomialIntegrate(L[i], a, b);
	}
	rules.withEndpoints = false;
	return rules;
}

template<typename T, int N>
QuadratureRules<T, N> gaussLegendreRules() {
	T const half(0.5), one(1.0), two(2.0);
	T const toler((1 << (2 + std::ilogb(N) / 3)) * std::numeric_limits<double>::epsilon());
	QuadratureRules<T, N> rules;
	T const dx = one / T(N);
	T const pi = T(M_PI);
	for (int l = 0; l < N; l++) {
		T dTheta;
		T theta = pi * (one - (T(l) + half) * dx);
		do {
			T const cosTheta = cos(theta);
			T const sinTheta = sin(theta);
			T const Pn = legendreP(N, cosTheta);
			T const dPnDx = dLegendrePdX(N, cosTheta);
			T const iDen = one / (sinTheta * dPnDx);
			dTheta = Pn * iDen;
			theta += dTheta;
		} while (abs(dTheta) > toler);
		T const x = cos(theta);
		rules.x[l] = x;
		rules.w[l] = two / (nSquared(dLegendrePdX(N, x)) * (one - nSquared(x)));
	}
	rules.withEndpoints = false;
	return rules;
}

template<typename T, int P>
QuadratureRules<T, P> clenshawCurtisRules() {
	constexpr int N = P - 1;
	T const zero(0), half(0.5), one(1), two(2), pi(M_PI);
	T const Ninv = one / T(N);
	T const piNinv = pi * Ninv;
	QuadratureRules<T, P> rules;
	Vector<T, N / 2 + 1> d;
	d[0] = one;
	for (int k = 1; k <= N / 2; k++) {
		d[k] = two / T(1 - nSquared(2 * k));
	}
	for (int n = 0; n <= N / 2; n++) {
		T const c = (T(n != 0) + one);
		T w = zero;
		for (int k = 0; k <= N / 2; k++) {
			w += c * d[k] * cos(T(2 * n * k) * piNinv) * Ninv;
		}
		rules.w[n] = w;
	}
	for (int k = 0; k < N - k; k++) {
		rules.w[N - k] = rules.w[k];
	}
	for (int k = 0; k <= N; k++) {
		rules.x[k] = cos(pi * T(k) / T(N));
	}
	return rules;
}

template<typename T, int N>
QuadratureRules<T, N> gaussLobattoRules() {
	T const half(0.5), one(1.0), two(2.0);
	T const toler((3 << (1 + (std::ilogb(N)) / 3)) * std::numeric_limits<double>::epsilon());
	QuadratureRules<T, N> rules;
	T const dx = one / T(N);
	T const pi = T(M_PI);
	for (int l = 1; l < N - 1; l++) {
		T dTheta;
		T theta = pi * (one - (T(l) + half) * dx);
		do {
			T const cosTheta = cos(theta);
			T const sinTheta = sin(theta);
			T const cotTheta = cosTheta / sinTheta;
			T const dPnDx = dLegendrePdX(N - 1, cosTheta);
			T const d2PnDx2 = d2LegendrePdX2(N - 1, cosTheta);
			T const iDen = one / (sinTheta * d2PnDx2 - dPnDx * cotTheta);
			dTheta = dPnDx * iDen;
			theta += dTheta;
		} while (abs(dTheta) > toler);
		T const x = cos(theta);
		rules.x[l] = x;
		rules.w[l] = two / (T(N * (N - 1)) * nSquared(legendreP(N - 1, x)));
	}
	rules.x[0] = -one;
	rules.x[N - 1] = one;
	rules.w[0] = rules.w[N - 1] = two / (T(N * (N - 1)));
	rules.withEndpoints = true;
	return rules;
}

template<typename T, int N, int D>
struct MultivariateQuadratureRules {
	static constexpr int N3 = Math::integerPower(N, D);
	Vector<Vector<T, D>, N3> x;
	Vector<T, N3> w;
	bool withEndpoints;
	static constexpr int size() {
		return N3;
	}
	MultivariateQuadratureRules(QuadratureRules<T, N> const &rules1d) {
		static constexpr T one = T(1);
		withEndpoints = rules1d.withEndpoints;
		for (int index = 0; index < N3; index++) {
			Vector<int, D> indices;
			int n = index;
			for (int dim = 0; dim < D; dim++) {
				auto const qr = std::div(n, N);
				indices[dim] = qr.rem;
				n = qr.quot;
			}
			for (int dim = 0; dim < D; dim++) {
				x[index][dim] = rules1d.x[indices[dim]];
			}
			T weight = one;
			for (int dim = 0; dim < D; dim++) {
				weight *= rules1d.w[indices[dim]];
			}
			w[index] = weight;
		}
	}
};

}

#endif
