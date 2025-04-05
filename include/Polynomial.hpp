/*
 * Polynomial.hpp
 *
 *  Created on: Dec 11, 2024
 *      Author: dmarce1
 */

#ifndef INCLUDE_POLYNOMIAL_HPP_
#define INCLUDE_POLYNOMIAL_HPP_

#include <numeric>
#include <map>
#include <cmath>
#include <math.h>
#include <string>
#include <vector>

#include "Complex.hpp"
#include "Numbers.hpp"
#include "Vector.hpp"

namespace Math {

template<typename Type>
struct Polynomial: public std::vector<Type> {
	using base_type = std::vector<Type>;
	Polynomial(Type d) {
		*this = std::vector<Type>(1, Type(d));
	}
	Polynomial() = default;
	~Polynomial() = default;
	Polynomial(base_type const &other) :
			base_type(other) {
	}
	Polynomial(Polynomial const&) = default;
	Polynomial& operator=(Polynomial const&) = default;
	int degree() const {
		return base_type::size() - 1;
	}
	Type operator()(Type const &x) const {
		Polynomial const &a = *this;
		Type f = Type(a[degree()]);
		for (int n = degree() - 1; n >= 0; n--) {
			f = x * f + a[n];
		}
		return f;
	}
	Polynomial<Type> operator()(Polynomial<Type> const &x) const {
		Polynomial const &a = *this;
		Polynomial<Type> f = Type(a[degree()]);
		for (int n = degree() - 1; n >= 0; n--) {
			f = x * f + a[n];
		}
		return f;
	}
	Complex operator()(Complex const &x) const {
		Polynomial const &a = *this;
		Complex f = Complex(Type(a[degree()]));
		for (int n = degree() - 1; n >= 0; n--) {
			f = x * f + a[n];
		}
		return f;
	}
	Type& operator[](int i) {
		if (int(base_type::size()) <= i) {
			base_type::resize(i + 1, Type(0));
		}
		return base_type::operator[](i);
	}
	Type const operator[](int i) const {
		if (int(base_type::size()) <= i) {
			return Type(0);
		} else {
			return base_type::operator[](i);
		}
	}
	Polynomial operator-() const {
		Polynomial A = *this;
		for (int n = 0; n <= degree(); n++) {
			A[n] = -A[n];
		}
		return A;
	}
	Polynomial operator+(Polynomial const &B) const {
		Polynomial A = *this;
		for (int l = 0; l <= B.degree(); l++) {
			A[l] += B[l];
		}
		return A;
	}
	Polynomial operator-(auto const &B) const {
		Polynomial A = *this;
		for (int l = 0; l <= B.degree(); l++) {
			A[l] -= B[l];
		}
		return A;
	}
	Polynomial operator*(Polynomial<Type> const &B) const {
		Polynomial const &A = *this;
		Polynomial C;
		int const M = A.degree();
		int const N = B.degree();
		C.resize(N + M + 1);
		for (int n = 0; n <= N + M; n++) {
			C[n] = Type(0);
			for (int m = std::max(0, n - N); m <= std::min(n, M); m++) {
				C[n] += A[m] * B[n - m];
			}
		}
		return C;
	}
	Polynomial& operator+=(auto const &other) {
		*this = *this + other;
		return *this;
	}
	Polynomial& operator-=(auto const &other) {
		*this = *this - other;
		return *this;
	}
	Polynomial& operator*=(Polynomial<Type> const &A) {
		*this = *this * A;
		return *this;
	}
	Polynomial& operator*=(Type const &a) {
		*this = a * *this;
		return *this;
	}
	Polynomial& operator/=(Type const &a) {
		*this = a / *this;
		return *this;
	}
};

template<typename Type>
std::string toString(Polynomial<Type> const &P) {
	using std::to_string;
	std::string str;
	str += "(";
	for (int n = 0; n < int(P.size()); n++) {
		str += to_string(P[n]);
		if (n + 1 < int(P.size())) {
			str += ", ";
		}
	}
	str += ")";
	return str;
}

template<typename Type>
Polynomial<Type> operator*(Type const &a, Polynomial<Type> const &B) {
	Polynomial<Type> C;
	for (int n = 0; n <= B.degree(); n++) {
		C[n] = a * B[n];
	}
	return C;
}

template<typename Type>
Polynomial<Type> operator*(Polynomial<Type> const &B, Type const &a) {
	return a * B;
}

template<typename Type>
Polynomial<Type> operator/(Polynomial<Type> const &B, Type const &a) {
	return B * (Real(1) / a);
}

template<typename Type>
auto polynomialDivision(Polynomial<Type> const &D, Polynomial<Type> const &I) {
	int const N = D.degree();
	int const M = I.degree();
	Polynomial<Type> Q;
	Polynomial<Type> R = D;
	for (int n = N; n >= M; n--) {
		Type const q = R[n] / I[M];
		R[n] = Type(0);
		for (int m = 1; m <= M; m++) {
			R[n - m] -= q * I[M - m];
		}
		Q[n - M] = q;
	}
	struct div_t {
		Polynomial<Type> quot;
		Polynomial<Type> rem;
	};
	return div_t { Q, R };
}

template<typename Type>
Polynomial<Type> polynomialDerivative(Polynomial<Type> F) {
	static constexpr Type zero = Type(0);
	auto const deg = F.degree();
	for (int n = 0; n < deg; n++) {
		F[n] = Type(n + 1) * F[n + 1];
		F[n + 1] = zero;
	}
	return F;
}

template<typename Type>
Polynomial<Type> polynomialAntiDerivative(Polynomial<Type> const &dfdx) {
	Polynomial<Type> f;
	for (int n = 0; n <= dfdx.degree(); n++) {
		f[n + 1] = dfdx[n] / Type(n + 1);
	}
	return f;
}

template<typename Type>
Polynomial<Type> polynomialDiscreteAntiDerivative(Polynomial<Type> const &P, size_t N) {
	Polynomial<Type> Q;
	for (int p = 0; p <= P.degree(); p++) {
		Q[p + 1] += P[p] / Type(p + 1);
		for (int j = 0; j < p; j++) {
			Q[p - j] -= P[p] * Math::nChooseK<size_t>(p, j) * std::riemann_zeta(-j);
		}
	}
	return Q;
}

template<typename Type>
Polynomial<Polynomial<Type>> polynomialAntiDerivative(Polynomial<Polynomial<Type>> const &dfdx) {
	Polynomial<Type> f;
	for (int n = 0; n <= dfdx.degree(); n++) {
		f[n + 1] = dfdx[n] / Type(n + 1);
	}
	for (int n = 1; n <= dfdx.degree() + 1; n++) {
		f[n] = polynomialAntiDerivative(f[n]);
	}
	return f;
}

template<typename Type>
Type polynomialIntegrate(Polynomial<Type> const &f, Type const &a, Type const &b) {
	auto const F = polynomialAntiDerivative(f);
	return F(b) - F(a);
}

template<typename Type>
Type polynomialFindRoot(Polynomial<Type> const &F, Type guess = 0.5) {
	constexpr int maxIters = 64;
	constexpr Type tolerance = Type(100) * std::numeric_limits < Type > ::epsilon();
	auto const dFdx = polynomialDerivative(F);
	auto const normInv = abs(Type(1) / F.back());
	Type x = guess;
	for (int iterCount = 0; iterCount < maxIters; iterCount++) {
		auto const f = F(x);
		auto const dfdx = dFdx(x);
		auto const dx = -f / dfdx;
		x += dx;
		auto const error = abs(f * normInv);
		if (error < tolerance) {
			return x;
		}
	}
	printf("Root solver failed to converge after %i iterations\n", maxIters);
	throw;
}

template<typename Type>
Complex polynomialFindComplexRoot(Polynomial<Type> const &F, Complex guess = Complex(0.5, 0.1)) { constexpr
	int maxIters = 1000;
	const Type tolerance = Type(1000 * std::numeric_limits<double>::epsilon());
	auto const dFdx = polynomialDerivative(F);
	auto const normInv = abs(Type(1) / F.back());
	Complex x = guess;
	for (int iterCount = 0; iterCount < maxIters; iterCount++) {
		Complex const f = F(x);
		Complex const dfdx = dFdx(x);
		Complex dx;
		dx = -f / dfdx;
		x += dx;
		Type const error = abs(f) * normInv;
		if (error < tolerance) {
			return x;
		}
	}
	printf("Complex root solver failed to converge after %i iterations\n", maxIters);
	throw;
}

template<class T>
T randOne() {
	static const T one = T(1);
	static const T two = T(2);
	static const T norm = one / (T(RAND_MAX) + one);
	return norm * (two * T(rand()) + one) - one;
}

template<typename Type>
std::vector<Complex> polynomialFindAllRoots(Polynomial<Type> const &F) {
	if (F.degree() == 0) {
		return std::vector<Complex>();
	} else if (F.degree() == 1) {
		std::vector<Complex> root(1, Complex(-F[0] / F[1]));
		return root;
	} else {
		static const Type tinyEps = Type(1000.0 * std::numeric_limits<double>::epsilon());
		std::vector<Complex> roots;
		Type upperBound = Type(0);
		Type const scaleInv = abs(Type(1) / F.back());
		for (int n = 0; n < F.degree(); n++) {
			Type const toPower = Type(1) / Type(F.degree() - n);
			Type const thisBound = pow(abs(F[n]) * scaleInv, toPower);
			upperBound = std::max(upperBound, thisBound);
		}
		upperBound *= Type(2);
		Type const lowerBound = Type(1) / upperBound;
		Type const theta = Type(M_PI) * randOne<Type>();
		Type const radius = abs(randOne<Type>()) * (upperBound - lowerBound) + lowerBound;
		Complex z = radius * exp(1.0_Imag * theta);
		z = polynomialFindComplexRoot(F, z);
		roots.push_back(z);
		Polynomial<Type> xmz;
		xmz[1] = Type(1);
		xmz[0] = -z.real();
		if (abs(z.imaginary()) < tinyEps) {
			auto const divResult = polynomialDivision(F, xmz);
			auto const nextRoots = polynomialFindAllRoots(divResult.quot);
			for (auto const &r : nextRoots) {
				roots.push_back(r);
			}
		} else {
			Polynomial<Type> xmz2;
			xmz2[0] = norm(z);
			xmz2[1] = Type(-2) * z.real();
			xmz2[2] = Type(1);
			z = conj(z);
			roots.push_back(z);
			auto const divResult = polynomialDivision(F, xmz2);
			auto const nextRoots = polynomialFindAllRoots(divResult.quot);
			for (auto const &r : nextRoots) {
				roots.push_back(r);
			}
		}
		return roots;
	}
}

}

#endif /* INCLUDE_POLYNOMIAL_HPP_ */
