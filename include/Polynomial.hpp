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

#include "Complex.hpp"
#include "Vector.hpp"

namespace Math {

template<typename Type, int Degree = 0>
struct Polynomial: public std::array<Type, Degree + 1> {
	using base_type = std::array<Type, Degree + 1>;
	Polynomial(Type init = Type(0)) {
		(*this)[0] = init;
		for (int n = 1; n <= Degree; n++) {
			(*this)[n] = Type(0);
		}
	}
	Type operator()(Type const &x) const {
		Polynomial const &a = *this;
		Type f = Type(a[Degree]);
		for (int n = Degree - 1; n >= 0; n--) {
			f = x * f + a[n];
		}
		return f;
	}
	Complex<Type> operator()(Complex<Type> const &x) const {
		Polynomial const &a = *this;
		Complex<Type> f = Type(a[Degree]);
		for (int n = Degree - 1; n >= 0; n--) {
			f = x * f + a[n];
		}
		return f;
	}
	template<int M>
	Polynomial<Type, std::max(M, Degree)> operator+(Polynomial<Type, M> const &B) const {
		constexpr int N = Degree;
		constexpr int MaxNM = std::max(N, M);
		constexpr int MinNM = std::min(N, M);
		Polynomial const &A = *this;
		Polynomial<Type, MaxNM> C;
		for (int l = 0; l <= MinNM; l++) {
			C[l] = A[l] + B[l];
		}
		for (int n = MinNM + 1; n <= N; n++) {
			C[n] = A[n];
		}
		for (int m = MinNM + 1; m <= M; m++) {
			C[m] = B[m];
		}
		return C;
	}
	Polynomial operator-() const {
		Polynomial A = *this;
		for (int n = 0; n <= Degree; n++) {
			A[n] = -A[n];
		}
		return A;
	}
	template<int M>
	auto operator*(Polynomial<Type, M> const &B) const {
		constexpr int N = Degree;
		Polynomial const &A = *this;
		Polynomial<Type, N + M> C;
		for (int n = 0; n <= N + M; n++) {
			for (int m = std::max(0, n - N); m <= std::min(n, M); m++) {
				C[n] += A[m] * B[n - m];
			}
		}
		return C;
	}
	Polynomial& operator*=(Polynomial<Type, 0> const &B) {
		Polynomial A = *this;
		for (int n = 0; n <= Degree; n++) {
			A[n] *= B[0];
		}
		return *this;
	}
	auto operator-(auto const &other) const {
		return *this + -other;
	}
	auto& operator+=(auto const &other) {
		*this = *this + other;
		return *this;
	}
	auto& operator-=(auto const &other) {
		*this = *this - other;
		return *this;
	}
	std::string toString() const {
		std::string str;
		str += "(";
		for (int n = 0; n < int(base_type::size()); n++) {
			str += std::to_string((*this)[n]);
			if (n + 1 < int(base_type::size())) {
				str += ", ";
			}
		}
		str += ")";
		return str;
	}
};

template<typename Type, int N, int M>
auto polynomialDivision(Polynomial<Type, N> D, Polynomial<Type, M> I) {
	Polynomial<Type, N - M> Q;
	Polynomial<Type, M - 1> R;
	for (int n = N; n >= M; n--) {
		Type const q = D[n] / I[M];
		D[n] -= D[n];
		for (int m = 1; m <= M; m++) {
			D[n - m] -= q * I[M - m];
		}
		Q[n - M] = q;
	}
	for (int n = 0; n < M; n++) {
		R[n] = D[n];
	}
	struct div_t {
		Polynomial<Type, N - M> quot;
		Polynomial<Type, M - 1> rem;
	};
	return div_t { Q, R };
}

template<typename Type, int N>
Polynomial<Type, N - 1> polynomialDerivative(Polynomial<Type, N> const &F) {
	Polynomial<Type, N - 1> dfdx;
	for (int n = 0; n < N; n++) {
		dfdx[n] = Type(n + 1) * F[n + 1];
	}
	return dfdx;
}

template<typename Type, int Degree>
Type polynomialFindRoot(Polynomial<Type, Degree> const &F, Type guess = 0.5) {
	constexpr int maxIters = 64;
	constexpr Type tolerance = Type(100) * std::numeric_limits<Type>::epsilon();
	auto const dFdx = polynomialDerivative(F);
	auto const normInv = abs(Type(1) / F[Degree]);
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

template<typename Type, int Degree>
Complex<Type> polynomialFindComplexRoot(Polynomial<Type, Degree> const &F,
		Complex<Type> guess = Complex<Type>(0.5, 0.1)) { constexpr
	int maxIters = 64;
	constexpr Type tolerance = Type(100) * std::numeric_limits<Type>::epsilon();
	auto const dFdx = polynomialDerivative(F);
//auto const d2Fdx2 = polynomialDerivative(dFdx);
	auto const normInv = std::abs(Type(1) / F[Degree]);
	Complex<Type> x = guess;
	for (int iterCount = 0; iterCount < maxIters; iterCount++) {
		Complex<Type> const f = F(x);
		//	Complex<Type> const d2fdx2 = d2Fdx2(x);
		Complex<Type> const dfdx = dFdx(x);
		Complex<Type> dx;
		dx = -f / dfdx;
		//if (abs(d2fdx2 * dx) > abs(dfdx)) {
		//	dx /= abs(dx);
		//	dx *= abs(f / (dx * d2fdx2));
		//}
		x += dx;
		Type const error = abs(f) * normInv;
		if (error < tolerance) {
			return x;
		}
	}
	printf("Complex root solver failed to converge after %i iterations\n", maxIters);
	throw;
}

template<class T = double>
T randOne() {
	static constexpr T one = 1;
	static constexpr T two = 2;
	static constexpr T norm = one / (T(RAND_MAX) + one);
	return norm * (two * T(rand()) + one) - one;
}

template<typename Type, int Degree>
std::array<Complex<Type>, Degree> polynomialFindAllRoots(Polynomial<Type, Degree> const &F) {
	if constexpr (Degree == 0) {
		return std::array<Complex<Type>, 0>();
	} else if constexpr (Degree == 1) {
		std::array<Complex<Type>, 1> root;
		root[0] = -F[0] / F[1];
		return root;
	} else {
		static constexpr Type tinyEps = Type(100) * std::numeric_limits<Type>::epsilon();
		std::array<Complex<Type>, Degree> roots;
		Type upperBound = Type(0);
		Type const scaleInv = fabs(Type(1) / F[Degree]);
		for (int n = 0; n < Degree; n++) {
			Type const toPower = Type(1) / Type(Degree - n);
			Type const thisBound = pow(fabs(F[n]) * scaleInv, toPower);
			upperBound = std::max(upperBound, thisBound);
		}
		upperBound *= Type(2);
		Type const lowerBound = Type(1) / upperBound;
		Type const theta = M_PI * randOne();
		Type const radius = fabs(randOne()) * (upperBound - lowerBound) + lowerBound;
		Complex<Type> z = radius * exp(1.0_Imag * theta);
		z = polynomialFindComplexRoot(F, z);
		roots.back() = z;
		Polynomial<Type, 1> xmz;
		xmz[0] = -z.real();
		xmz[1] = Type(1);
		if (std::abs(z.imaginary()) < tinyEps) {
			auto const divResult = polynomialDivision(F, xmz);
			auto const nextArray = polynomialFindAllRoots(divResult.quot);
			for (int n = 0; n < Degree - 1; n++) {
				roots[n] = nextArray[n];
			}
		} else {
			Polynomial<Type, 2> xmz2;
			xmz2[0] = norm(z);
			xmz2[1] = Type(-2) * z.real();
			xmz2[2] = Type(1);
			z = conj(z);
			roots[Degree - 2] = conj(roots[Degree - 1]);
			auto const divResult = polynomialDivision(F, xmz2);
			auto const nextArray = polynomialFindAllRoots(divResult.quot);
			for (int n = 0; n < Degree - 2; n++) {
				roots[n] = nextArray[n];
			}
		}
		return roots;
	}
}

}

#endif /* INCLUDE_POLYNOMIAL_HPP_ */
