/******************************************************************************
 Copyright (C) 2024  Dominic C. Marcello
 *******************************************************************************/

#ifndef INCLUDE_HALFINTEGER_HPP_
#define INCLUDE_HALFINTEGER_HPP_

#include <limits>
#include <cmath>

template<typename Integer = int>
struct HalfInteger {
	friend constexpr HalfInteger operator+(HalfInteger I) {
		return I;
	}
	friend constexpr HalfInteger operator-(HalfInteger J) {
		HalfInteger I;
		I.i_ = -J.i_;
		return I;
	}
	friend constexpr HalfInteger operator+(HalfInteger J, HalfInteger K) {
		HalfInteger I;
		I.i_ = J.i_ + K.i_;
		return I;
	}
	friend constexpr HalfInteger operator+(Integer j, HalfInteger K) {
		HalfInteger I;
		I.i_ = 2 * j + K.i_;
		return I;
	}
	friend constexpr HalfInteger operator+(HalfInteger J, Integer k) {
		return k + J;
	}
	friend constexpr HalfInteger operator-(HalfInteger J, HalfInteger K) {
		HalfInteger I;
		I.i_ = J.i_ - K.i_;
		return I;
	}
	friend constexpr HalfInteger operator-(Integer j, HalfInteger K) {
		HalfInteger I;
		I.i_ = 2 * j - K.i_;
		return I;
	}
	friend constexpr HalfInteger operator-(HalfInteger J, Integer k) {
		HalfInteger I;
		I.i_ = J.i_ - 2 * k;
		return I;
	}
	friend constexpr HalfInteger operator*(Integer j, HalfInteger K) {
		HalfInteger I;
		I.i_ = j * K.i_;
		return I;
	}
	friend constexpr HalfInteger operator*(HalfInteger J, Integer k) {
		return k * J;
	}
	friend constexpr double operator*(HalfInteger J, double k) {
		return k * J.i_ * 0.5;
	}
	friend constexpr double operator*(double j, HalfInteger K) {
		return j * K.i_ * 0.5;
	}
	friend constexpr HalfInteger operator/(Integer j, HalfInteger K) {
		HalfInteger I;
		I.i_ = (4 * j) / K.i_;
		return I;
	}
	friend constexpr HalfInteger operator/(HalfInteger J, Integer k) {
		HalfInteger I;
		I.i_ = J.i_ / k * 0.5;
		return I;
	}
	friend constexpr double operator/(HalfInteger J, double k) {
		return k / J.i_ * 2.0;
	}
	friend constexpr HalfInteger operator/(double j, HalfInteger K) {
		return j / K.i_;
	}
	constexpr HalfInteger& operator+=(HalfInteger I) {
		*this = *this + I;
		return *this;
	}
	constexpr HalfInteger& operator+=(int i) {
		*this = *this + i;
		return *this;
	}
	constexpr HalfInteger& operator-=(HalfInteger I) {
		*this = *this - I;
		return *this;
	}
	constexpr HalfInteger& operator-=(Integer i) {
		*this = *this - i;
		return *this;
	}
	constexpr HalfInteger& operator*=(HalfInteger I) {
		*this = *this * I;
		return *this;
	}
	constexpr HalfInteger& operator*=(Integer i) {
		*this = *this * i;
		return *this;
	}
	constexpr HalfInteger& operator/=(HalfInteger I) {
		*this = *this / I;
		return *this;
	}
	constexpr HalfInteger& operator/=(Integer i) {
		*this = *this / i;
		return *this;
	}
	constexpr HalfInteger& operator++() {
		i_ += 2;
		return *this;
	}
	constexpr HalfInteger operator++(int) {
		auto const before = *this;
		++(*this);
		return before;
	}
	constexpr operator double() const {
		return i_ / 2.0;
	}
	constexpr bool operator==(HalfInteger other) const {
		return i_ == other.i_;
	}
	constexpr bool operator!=(HalfInteger other) const {
		return i_ != other.i_;
	}
	constexpr bool operator>(HalfInteger other) const {
		return i_ > other.i_;
	}
	constexpr bool operator<(HalfInteger other) const {
		return i_ < other.i_;
	}
	constexpr bool operator>=(HalfInteger other) const {
		return i_ >= other.i_;
	}
	constexpr bool operator<=(HalfInteger other) const {
		return i_ <= other.i_;
	}
	constexpr bool isInteger() const {
		return !(abs(i_) & 1);
	}
	constexpr Integer toInteger() const {
		return i_ / 2;
	}
	template<typename T>
	friend constexpr T pow(T z, HalfInteger n) {
		if (n < 0) {
			return 1.0 / pow(z, -n);
		} else if (n > 0) {
			if (n.i_ & 1) {
				return sqrt(ipow(z, n.i_));
			}
			return ipow(z, n.i_ >> 1);
		}
		return T(1);
	}
	friend constexpr double tgamma(HalfInteger z) {
		constexpr double pi = M_PI;
		constexpr double infinity = std::numeric_limits<double>::infinity();
		if (z < 0) {
			if (((-z.i_) & 1) == 0) {
				return infinity;
			}
			z += 1;
			return pi / (sin(pi * z) * tgamma(z));
		} else if (z > 0) {
			if (z.i_ == 1) {
				return sqrt(pi);
			}
			if (z.i_ == 2) {
				return 1.0;
			}
			z -= 1;
			return z * tgamma(z);
		}
		return infinity;
	}
	constexpr HalfInteger() = default;
	constexpr HalfInteger(HalfInteger const &I) = default;
	constexpr HalfInteger(HalfInteger &&I) = default;
	constexpr HalfInteger& operator=(HalfInteger const &I) = default;
	constexpr HalfInteger& operator=(HalfInteger &&I) = default;
	static constexpr HalfInteger zero() {
		HalfInteger I;
		I.i_ = 0;
		return I;
	}
	static constexpr HalfInteger half() {
		HalfInteger I;
		I.i_ = 1;
		return I;
	}
	static constexpr HalfInteger one() {
		HalfInteger I;
		I.i_ = 1;
		return I;
	}
private:
	Integer i_;
};

#endif /* INCLUDE_HALFINTEGER_HPP_ */
