///******************************************************************************
// Copyright (C) 2024  Dominic C. Marcello
// *******************************************************************************/
//
//#ifndef INCLUDE_MULTIPRECISION_HPP_
//#define INCLUDE_MULTIPRECISION_HPP_
//
//#include <cmath>
//#include <cstdint>
//#include <numeric>
//#include <gmp.h>
//#include <mpfr.h>
//
//#define MULTI_PRECISION_FUNCTION1( name )                  \
//	friend MultiPrecision name (MultiPrecision const& a) { \
//		MultiPrecision b;                                  \
//		mpfr_##name (b.value_, a.value_, roundingMode);    \
//		return b;                                          \
//	}
//
//#define MULTI_PRECISION_FUNCTION2( name )                                    \
//	friend MultiPrecision name (MultiPrecision const& a, MultiPrecision const& b) { \
//		MultiPrecision c;                                                            \
//		mpfr_##name (c.value_, a.value_, b.value_, roundingMode);                   \
//		return c;                                                                    \
//	}
//
//template<int precision, mpfr_rnd_t roundingMode = MPFR_RNDN>
//struct MultiPrecision {
//	MultiPrecision() {
//		mpfr_init2(value_, precision);
//	}
//	MultiPrecision(long double other) {
//		mpfr_init2(value_, precision);
//		*this = other;
//	}
//	MultiPrecision(MultiPrecision const &other) {
//		mpfr_init2(value_, precision);
//		*this = other;
//	}
//	MultiPrecision& operator=(MultiPrecision const &other) {
//		mpfr_set(value_, other.value_, roundingMode);
//		return *this;
//	}
//	MultiPrecision& operator=(long double real) {
//		mpfr_set_ld(value_, real, roundingMode);
//		return *this;
//	}
//	MultiPrecision operator-() const {
//		MultiPrecision b;
//		mpfr_neg(b.value_, value_, roundingMode);
//		return b;
//	}
//	MultiPrecision operator+() const {
//		return *this;
//	}
//	MultiPrecision operator+(MultiPrecision const &b) const {
//		MultiPrecision c;
//		mpfr_add(c.value_, value_, b.value_, roundingMode);
//		return c;
//	}
//	MultiPrecision operator-(MultiPrecision const &a) const {
//		return *this + -a;
//	}
//	MultiPrecision operator*(MultiPrecision const &b) const {
//		MultiPrecision c;
//		mpfr_mul(c.value_, value_, b.value_, roundingMode);
//		return c;
//	}
//	MultiPrecision operator/(MultiPrecision const &b) const {
//		MultiPrecision c;
//		mpfr_div(c.value_, value_, b.value_, roundingMode);
//		return c;
//	}
//	MultiPrecision& operator+=(MultiPrecision const &a) {
//		*this = *this + a;
//		return *this;
//	}
//	MultiPrecision& operator-=(MultiPrecision const &a) {
//		*this = *this - a;
//		return *this;
//	}
//	MultiPrecision& operator*=(MultiPrecision const &a) {
//		*this = *this * a;
//		return *this;
//	}
//	MultiPrecision& operator/=(MultiPrecision const &a) {
//		*this = *this / a;
//		return *this;
//	}
//	operator double() const {
//		return mpfr_get_d(value_, roundingMode);
//	}
//	friend bool operator==(MultiPrecision const &a, MultiPrecision const &b) {
//		return bool(mpfr_cmp(a.value_, b.value_) == 0);
//	}
//	friend bool operator!=(MultiPrecision const &a, MultiPrecision const &b) {
//		return bool(mpfr_cmp(a.value_, b.value_) != 0);
//	}
//	friend bool operator>(MultiPrecision const &a, MultiPrecision const &b) {
//		return bool(mpfr_cmp(a.value_, b.value_) > 0);
//	}
//	friend bool operator<(MultiPrecision const &a, MultiPrecision const &b) {
//		return bool(mpfr_cmp(a.value_, b.value_) < 0);
//	}
//	friend bool operator<=(MultiPrecision const &a, MultiPrecision const &b) {
//		return !(b < a);
//	}
//	friend bool operator>=(MultiPrecision const &a, MultiPrecision const &b) {
//		return !(b > a);
//	}
//	friend MultiPrecision nexttoward(MultiPrecision const &a, MultiPrecision const &b) {
//		MultiPrecision c = a;
//		mpfr_nexttoward(c.value_, b.value_);
//		return c;
//	}
//	friend MultiPrecision legendre(int n, MultiPrecision const &x) {
//		MultiPrecision p2, p1, p0;
//		p0 = MultiPrecision(0);
//		p1 = MultiPrecision(1);
//		for (int l = 0; l < n; l++) {
//			p2 = (MultiPrecision(2 * l + 1) * x * p1 - MultiPrecision(l) * p0) / MultiPrecision(l + 1);
//			p0 = p1;
//			p1 = p2;
//		}
//		return p1;
//	}
//	MULTI_PRECISION_FUNCTION1(abs)
//	MULTI_PRECISION_FUNCTION1(cos)
//	MULTI_PRECISION_FUNCTION1(sin)
//	MULTI_PRECISION_FUNCTION1(tan)
//	MULTI_PRECISION_FUNCTION1(cosh)
//	MULTI_PRECISION_FUNCTION1(sinh)
//	MULTI_PRECISION_FUNCTION1(tanh)
//	MULTI_PRECISION_FUNCTION1(acos)
//	MULTI_PRECISION_FUNCTION1(asin)
//	MULTI_PRECISION_FUNCTION1(atan)
//	MULTI_PRECISION_FUNCTION1(acosh)
//	MULTI_PRECISION_FUNCTION1(asinh)
//	MULTI_PRECISION_FUNCTION1(atanh)
//	MULTI_PRECISION_FUNCTION1(exp)
//	MULTI_PRECISION_FUNCTION1(log)
//	MULTI_PRECISION_FUNCTION1(sqrt)
//	MULTI_PRECISION_FUNCTION2(copysign)
//	MULTI_PRECISION_FUNCTION2(max)
//	MULTI_PRECISION_FUNCTION2(min)
//	MULTI_PRECISION_FUNCTION2(pow)
//	template<int, int>
//	friend struct std::numeric_limits;
//private:
//	mpfr_t value_;
//};
//
//namespace std {
//template<int prec, mpfr_rnd_t rndMode>
//struct numeric_limits<MultiPrecision<prec, rndMode>> {
//	static constexpr bool is_specialized = true;
//	static constexpr bool is_signed = true;
//	static constexpr bool is_integer = false;
//	static constexpr bool is_exact = false;
//	static constexpr int radix = 2;
//
//	static constexpr int digits = prec; // bits of precision
//
//	static MultiPrecision<prec, rndMode> min() {
//		MultiPrecision<prec, rndMode> m;
//		mpfr_set_ui_2exp(m.value_, 1, mpfr_get_emin(), rndMode);
//		return m;
//	}
//	static MultiPrecision<prec, rndMode> max() {
//		MultiPrecision<prec, rndMode> m;
//		mpfr_set_inf(m.value_, 1);
//		mpfr_nextbelow(m.value_);
//		return m;
//	}
//	static MultiPrecision<prec, rndMode> epsilon() {
//		MultiPrecision<prec, rndMode> a(1);
//		MultiPrecision<prec, rndMode> b(a);
//		MultiPrecision<prec, rndMode> eps;
//		mpfr_nextabove(b.value_);
//		mpfr_sub(b.value_, b.value_, a.value_, rndMode);
//		return b;
//	}
//	static MultiPrecision<prec, rndMode> infinity() {
//		MultiPrecision<prec, rndMode> m;
//		mpfr_set_inf(m.value_, 1);
//		return m;
//	}
//	static MultiPrecision<prec, rndMode> quiet_NaN() {
//		MultiPrecision<prec, rndMode> m;
//		mpfr_set_nan(m.value_);
//		return m;
//	}
//};
//}
//
//template<typename T>
//T legendreP(int n, T const &x) {
//	T p2, p1, p0;
//	p0 = 0;
//	p1 = 1;
//	for (int l = 0; l < n; l++) {
//		p2 = ((2 * l + 1) * x * p1 - l * p0) / (l + 1);
//		p0 = p1;
//		p1 = p2;
//	}
//	return p1;
//}
//
//#endif /* INCLUDE_MULTIPRECISION_HPP_ */
