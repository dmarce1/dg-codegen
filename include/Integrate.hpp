/******************************************************************************
 Copyright (C) 2024  Dominic C. Marcello
 *******************************************************************************/

#ifndef INCLUDE_INTEGRATE_HPP_
#define INCLUDE_INTEGRATE_HPP_

#include <functional>

#include "Quadrature.hpp"
#include "Util.hpp"

template<typename T>
struct KahnSum {
	KahnSum() = default;
	KahnSum(T const &a) {
		sum_ = a;
		c_ = T(0);
	}
	KahnSum& operator+=(T const &x) {
		T const y = x - c_;
		T const t = sum_ + y;
		c_ = (t - sum_) - y;
		sum_ = t;
		return *this;
	}
	KahnSum& operator-=(T const &x) {
		*this += -x;
		return *this;
	}
	KahnSum operator+(T const &x) const {
		KahnSum a = *this;
		a += x;
		return a;
	}
	KahnSum operator-(T const &x) const {
		KahnSum a = *this;
		a -= x;
		return a;
	}
	operator T() const {
		return sum_;
	}
	bool operator==(KahnSum const &other) const {
		return (sum_ == other.sum_) && (c_ == other.c_);
	}
	bool operator<(KahnSum const &other) const {
		if (sum_ < other.sum_) {
			return true;
		}
		return (c_ > other.c_);
	}
	bool operator!=(KahnSum const &other) const {
		return !(*this == other);
	}
	bool operator>=(KahnSum const &other) const {
		return !(*this < other);
	}
	friend bool operator>(KahnSum const &a, KahnSum const &b) {
		return b < a;
	}
	friend bool operator<=(KahnSum const &a, KahnSum const &b) {
		return b >= a;
	}
private:
	T sum_ { 0 };
	T c_ { 0 };
};

#define frac(n, d) (T(n)/T(d))

//template<typename T, int nq = 100>
//T integrate(std::function<T(T const&)> const &f, T a, T b, int nh) {
//	static auto const quad = gaussLegendreQuadraturePoints<T>(nq);
//	T constexpr zero = T(0);
//	T constexpr one = T(1);
//	T constexpr two = T(2);
//	T constexpr half = one / two;
//	T constexpr infinity = std::numeric_limits<T>::infinity();
//	T sgn = one;
//	std::function<T(T const&)> g;
//	if (b < a) {
//		std::swap(a, b);
//		sgn = -sgn;
//	}
//	if (a == -infinity) {
//		if (b == +infinity) {
//			g = [f, b](T const &t) {
//				T const y = one / sqrt(one - t);
//				T const x = t * y;
//				return y * sqr(y) * f(x);
//			};
//			a = -one;
//			b = +one;
//		} else {
//			g = [f, b](T const &t) {
//				T const y = one / (one - t);
//				T const x = b + (t - b) * y;
//				return sqr(y) * f(x);
//			};
//		}
//		a = -one;
//	} else {
//		if (b == +infinity) {
//			g = [f, a](T const &t) {
//				T const y = one / (one - t);
//				T const x = a + (t - a) * y;
//				return sqr(y) * f(x);
//			};
//			b = +one;
//		} else {
//			g = f;
//		}
//	}
//	KahnSum<T> F = zero;
//	T const dn = (b - a) / T(nh);
//	T const dq = half * dn;
//	for (int n = 0; n < nh; n++) {
//		for (int q = 0; q < nq; q++) {
//			T const x = quad[q].position;
//			T const w = quad[q].weight;
//			T const t = a + (n + half * (x + one)) * dn;
//			F += g(t) * w * dq;
//		}
//	}
//	return sgn * T(F);
//}

template<typename T>
struct IntegrationSegment {
	T sum;
	T a;
	T b;
	int nq;
};

constexpr int Q0 = 16;

template<typename T, int N = 16, bool first = true>
T integrate(std::function<T(T const&)> const &g, T a, T b, T I0 = 0) {
	using std::numeric_limits;
	using std::sqrt;
	using std::min;
	using std::abs;
	constexpr int Nq = std::min(1024, 1 << (std::ilogb(sqrt(N) - 1) + 1));
	constexpr int Nh = N / Nq;
	auto const &quad = gaussLegendreQuadraturePoints<T>(Nq);
	static T const ε = numeric_limits<T>::epsilon();
	static T const smallest = numeric_limits<T>::min();
	static T const zero = T(0);
	static T const one = T(1);
	static T const two = T(2);
	static T const half = one / two;
	static T const infinity = numeric_limits<T>::infinity();
	std::function<T(T const&)> f;
	T sgn = one;
	if constexpr (first) {
		if (b < a) {
			std::swap(a, b);
			sgn = -sgn;
		}
		if (a == -infinity) {
			if (b == +infinity) {
				f = [g, b](T const &t) {
					T const y = one / sqrt(one - t);
					T const x = t * y;
					return y * sqr(y) * g(x);
				};
				a = -one;
				b = +one;
			} else {
				f = [g, b](T const &t) {
					T const y = one / (one - t);
					T const x = b + (t - b) * y;
					return sqr(y) * g(x);
				};
			}
			a = -one;
		} else {
			if (b == +infinity) {
				f = [g, a](T const &t) {
					T const y = one / (one - t);
					T const x = a + (t - a) * y;
					return sqr(y) * g(x);
				};
				b = +one;
			} else {
				f = g;
			}
		}
	} else {
		f = g;
	}
	T const dx = (b - a) / T(Nh);
	KahnSum<T> I = zero;
	for (int i = 0; i < Nh; i++) {
		T const x0 = a + T(i) * dx;
		for (int q = 0; q < Nq; q++) {
			T const x = quad[q].position;
			T const w = quad[q].weight;
			T const t = x0 + half * (x + one) * dx;
			I += half * dx * w * f(t);
		}
	}
	if constexpr (first) {
		return integrate<T, 2 * N, false>(f, a, b, T(I) * sgn);
	} else {
		T const err = abs(T(I - I0) / T(I + smallest));
		//printf("%i %e\n", N, err);
		if (err >= 128*ε) {
			if constexpr (N < std::numeric_limits<int>::max() / 2) {
				return integrate<T, 2 * N, false>(f, a, b, I);
			} else {
				printf("!!!!!!\n");
				abort();
			}
		}
		return I;
	}
}

template<typename T, int N = 1024>
std::complex<T> complexIntegrate(std::function<std::complex<T>(T const&)> const &f, T a, T b) {
	using std::numeric_limits;
	using std::sqrt;
	using std::min;
	using std::abs;
	constexpr int Nq = std::min(1024, 1 << (std::ilogb(sqrt(N) - 1) + 1));
	constexpr int Nh = N / Nq;
	auto const &quad = gaussLegendreQuadraturePoints<T>(Nq);
	static T const ε = T(8) * numeric_limits<T>::epsilon();
	static T const smallest = numeric_limits<T>::min();
	static T const zero = T(0);
	static T const one = T(1);
	static T const two = T(2);
	static T const half = one / two;
	static T const infinity = numeric_limits<T>::infinity();
	T sgn = one;
	T const dx = (b - a) / T(Nh);
	KahnSum<std::complex<T>> I = std::complex<T>(zero, zero);
	for (int i = 0; i < Nh; i++) {
		T const x0 = a + T(i) * dx;
		for (int q = 0; q < Nq; q++) {
			T const x = quad[q].position;
			T const w = quad[q].weight;
			T const t = x0 + half * (x + one) * dx;
			I += half * dx * w * f(t);
		}
	}
	return std::complex<T>(I);
}


template<typename T, typename F>
T rhombergIntegrate(F f, T a, T b, int Nlev = 32) {
	constexpr int kb = 1;
	std::vector<KahnSum<T>> r0, r1;
	r0.resize((1 << Nlev) + 1);
	r1.resize((1 << Nlev) + 1);
	T h = b - a;
	r1[0] = frac(1, 2) * h * (f(a) + f(b));
	int ke = 1;
	for (int n = 1; n < Nlev; n++) {
		std::swap(r0, r1);
		r1.resize(n + 1);
		h = (b - a) / (1 << (n + 1));
		ke <<= 1;
		r1[0] = frac(1, 2) * r0[0];
		for (int k = kb; k <= ke; k++) {
			T const x = a + (2 * k - 1) * h;
			r1[0] += h * f(x);
		}
		for (int m = 1; m <= n; m++) {
			r1[m] = r1[m - 1] + (r1[m - 1] - r0[m - 1]) / sqr(T((1 << m) - 1));
		}
		long long err = ulpDistance(T(r0.back()), T(r1.back()));
		if (err <= 4) {
			break;
		}
	}
	return r1.back();
}

#endif /* INCLUDE_INTEGRATE_HPP_ */
