/******************************************************************************
 Copyright (C) 2024  Dominic C. Marcello
 *******************************************************************************/

#ifndef INCLUDE_DOUBLEREAL_HPP_
#define INCLUDE_DOUBLEREAL_HPP_

#include <cmath>

template<typename T>
struct DoubleReal {
	T h_;
	T l_;
	static inline auto fast2Sum(T a, T b) {
		auto volatile const s = a + b;
		auto volatile const z = s - a;
		auto volatile const t = s - a;
		return DoubleReal(s, t);;
	}
	static inline auto _2Sum(T a, T b) {
		auto volatile const s = a + b;
		auto volatile const α = s - b;
		auto volatile const β = s - α;
		auto volatile const δ = a - α;
		auto volatile const Δ = b - β;
		auto volatile const t = δ + Δ;
		return DoubleReal(s, t);
	}
	static inline auto _2Prod(T a, T b) {
		auto volatile const π = a * b;
		auto volatile const ρ = std::fma(a, b, -π);
		return DoubleReal(π, ρ);
	}
	static inline auto _2Inv(DoubleReal const &x) {
		auto volatile const y = T(1) / x.h_;
		return y * (T(1) - y * x);
	}
public:
	auto& operator=(DoubleReal const &other) {
		h_ = other.h_;
		l_ = other.l_;
		return *this;
	}
	auto& operator=(T other) {
		h_ = other;
		l_ = T(0);
		return *this;
	}
	friend auto operator-(DoubleReal x) {
		x.h_ = -x.h_;
		x.l_ = -x.l_;
		return x;
	}
	friend auto operator+(DoubleReal const &x, DoubleReal const &y) {
		auto const s = _2Sum(x.h_, y.h_);
		auto const t = _2Sum(x.l_, y.l_);
		auto const c = s.l_ + t.h_;
		auto const v = fast2Sum(s.h_, c);
		auto const w = t.l_ + v.l_;
		auto const z = fast2Sum(v.h_, w);
		return z;
	}
	friend auto operator-(DoubleReal const &x, DoubleReal const &y) {
		return x + -y;
	}
	friend auto operator*(DoubleReal const &x, DoubleReal const &y) {
		auto const c = _2Prod(x.h_, y.h_);
		auto const volatile t1 = x.h_ * y.l_;
		auto const volatile t2 = x.l_ * y.h_;
		auto const volatile c1 = t1 + t2;
		auto const volatile c2 = c.l_ + c1;
		auto const z = fast2Sum(c.h_, c2);
		return z;
	}
	friend auto operator/(DoubleReal const &x, DoubleReal const &y) {
		return x * _2Inv(y);
	}
	friend auto operator+(DoubleReal const &x, T y) {
		auto const s = _2Sum(x.h_, y);
		auto const v = x.l_ + s.l_;
		auto const z = fast2Sum(s.h_, v);
		return z;

	}
	friend auto operator-(DoubleReal const &x, T y) {
		return x + -y;

	}
	friend auto operator*(DoubleReal const &x, T y) {
		auto const s = _2Prod(x.h_, y);
		auto const volatile t = x.l_ * y;
		auto const volatile c = s.l_ + t;
		auto const z = fast2Sum(s.h_, c);
		return z;
	}
	friend auto operator/(DoubleReal const &x, T y) {
		return x * (T(1) / y);
	}
	friend auto operator+(T x, DoubleReal const &y) {
		return y + x;
	}
	friend auto operator-(T x, DoubleReal const &y) {
		return -y + x;
	}
	friend auto operator*(T x, DoubleReal const &y) {
		return y * x;
	}
	friend auto operator/(T x, DoubleReal const &y) {
		return x * _2Inv(y);
	}
	friend auto exp(DoubleReal const &y) {
		auto const s = std::exp(y.h_);
		auto const t = y.l_ * s;
		return _2Sum(s, t);
	}
	friend auto sqrt(DoubleReal const &y) {
		auto const s = std::sqrt(y.h_);
		auto const t = T(0.5) * y.l_ / s;
		return _2Sum(s, t);
	}
	friend auto pow(DoubleReal const &y, T n) {
		auto const s = std::pow(y.h_, n);
		auto const t = n * s * y.l_ / y.h_;
		return _2Sum(s, t);
	}
};

#endif /* INCLUDE_DOUBLEREAL_HPP_ */
