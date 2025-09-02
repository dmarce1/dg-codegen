/******************************************************************************
 Copyright (C) 2024  Dominic C. Marcello
 *******************************************************************************/

#ifndef INCLUDE_PADEAPPROXIMANT_HPP_
#define INCLUDE_PADEAPPROXIMANT_HPP_

#include "Matrix.hpp"
#include <functional>

template<int N>
struct PadeApproximant {
	PadeApproximant(std::function<long double(long double, int)> const &f, long double x0 = 0.0) {
		std::array<long double, 2 * N> c;
		std::array<long double, N> r;
		SquareMatrix<long double, N> A;
		for (int n = 0; n < 2 * N; n++) {
			c[n] = f(x0, n) / factorial(n);
		}
		for (int n = 0; n < N; n++) {
			r[n] = -c[N + n];
			for (int m = 0; m < N; m++) {
				A(n, m) = c[n - m + N - 1];
			}
		}
		b_[0] = static_cast<long double>(1);
		auto Ainv = A;
		matrixInverseAndDeterminant(Ainv);
		for (int n = 0; n < N; n++) {
			b_[n + 1] = 0.0;
			for (int m = 0; m < N; m++) {
				b_[n + 1] += Ainv(n, m) * r[m];
			}
		}
		for (int n = 0; n < N; n++) {
			a_[n] = c[n];
			for (int m = 1; m <= n; m++) {
				a_[n] += c[n - m] * b_[m];
			}
		}
	}
	double operator()(double x) const {
		long double y = 0.0;
		long double z = 0.0;
		for (int n = N - 1; n >= 0; n--) {
			y = y * x + a_[n];
		}
		for (int n = N; n >= 0; n--) {
			z = z * x + b_[n];
		}
		return y / z;
	}
private:
	std::array<long double, N + 1> b_;
	std::array<long double, N> a_;
};

#endif /* INCLUDE_PADEAPPROXIMANT_HPP_ */
