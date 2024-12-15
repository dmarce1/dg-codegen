/*
 * Real.hpp
 *
 *  Created on: Dec 11, 2024
 *      Author: dmarce1
 */

#ifndef INCLUDE_REAL_HPP_
#define INCLUDE_REAL_HPP_

#include <limits>
#include <type_traits>

namespace Math {

template<typename >
struct Real;

template<typename Type = double>
struct Real {
	constexpr Real() {
		if constexpr (CHECK_REALS) {
			value = std::numeric_limits<Type>::signaling_NaN();
		}
	}
	constexpr Real(Type _value) :
			value(_value) {
	}
	Real& operator=(Real const &a) {
		value = a.value;
		return *this;
	}
	Real operator+() const {
		return *this;
	}
	Real operator-() const {
		return Real(-value);
	}
	Real operator+(Real const &a) const {
		return Real(value + a.value);
	}
	Real operator-(Real const &a) const {
		return Real(value - a.value);
	}
	Real operator*(Real const &a) const {
		return Real(value * a.value);
	}
	Real operator/(Real const &a) const {
		return Real(value / a.value);
	}
	Real& operator+=(Real const &a) {
		*this = *this + a;
		return *this;
	}
	Real& operator-=(Real const &a) {
		*this = *this - a;
		return *this;
	}
	Real& operator*=(Real const &a) {
		*this = *this * a;
		return *this;
	}
	Real& operator/=(Real const &a) {
		*this = *this / a;
		return *this;
	}
	bool operator==(Real const &a) const {
		return value == a.value;
	}
	bool operator!=(Real const &a) const {
		return value != a.value;
	}
	bool operator<=(Real const &a) const {
		return value <= a.value;
	}
	bool operator>=(Real const &a) const {
		return value >= a.value;
	}
	bool operator<(Real const &a) const {
		return value < a.value;
	}
	bool operator>(Real const &a) const {
		return value > a.value;
	}
private:
	Type value;
	void debug_check(Real a = Real(0)) const {
		if constexpr (CHECK_REALS) {
			if (!isfinite(a.value) || !isfinite(value)) {
				assert(false);
			}
		}
	}
};

}

#endif /* INCLUDE_REAL_HPP_ */
