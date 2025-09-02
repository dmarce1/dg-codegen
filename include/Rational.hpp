/******************************************************************************
 Copyright (C) 2024  Dominic C. Marcello
 *******************************************************************************/

#ifndef INCLUDE_RATIONAL_HPP_
#define INCLUDE_RATIONAL_HPP_

#include "MultiPrecision.hpp"

struct Rational {
	void normalize() {
		auto const d = std::copysign(std::gcd(n_, d_), d_);
		n_ /= d;
		d_ /= d;
	}
	Rational(intmax_t n, intmax_t d = 1) :
			n_(n), d_(d) {
	}
	Rational& operator=(intmax_t integer) {
		n_ = integer;
		d_ = 1;
		return *this;
	}
	friend Rational operator-(Rational const &b) {
		Rational a;
		a.n_ = -b.n_;
		a.d_ = b.d_;
		return a;
	}
	friend Rational operator+(Rational const &b, Rational const &c) {
		Rational a;
		a.n_ = b.d_ * c.n_ + b.n_ * c.d_;
		a.d_ = b.d_ * c.d_;
		a.normalize();
		return a;
	}
	friend Rational operator-(Rational const &a, Rational const &b) {
		return a + -b;
	}
	friend Rational operator*(Rational const &b, Rational const &c) {
		Rational a;
		a.n_ = b.n_ * c.n_;
		a.d_ = b.d_ * c.d_;
		a.normalize();
		return a;
	}
	friend Rational operator/(Rational const &b, Rational const &c) {
		Rational a;
		a.n_ = b.n_ * c.d_;
		a.d_ = b.d_ * c.n_;
		a.normalize();
		return a;
	}
	template<typename T>
	operator T() const {
		return T(n_) / T(d_);
	}
private:
	intmax_t n_ = 0;
	intmax_t d_ = 1;
};

#endif /* INCLUDE_RATIONAL_HPP_ */
