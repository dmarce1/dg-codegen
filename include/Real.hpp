/*
 * Real.hpp
 *
 *  Created on: Dec 11, 2024
 *      Author: dmarce1
 */

#ifndef INCLUDE_REAL_HPP_
#define INCLUDE_REAL_HPP_

#include <cassert>
#include <cmath>
#include <iostream>
#include <limits>
#include <stacktrace>
#include <stdexcept>
#include <type_traits>

namespace Math {

struct Real {
	using Type = double;
	constexpr Real() {
		if constexpr (CHECK_REALS) {
			value = std::numeric_limits<Type>::signaling_NaN();
		}
	}
	explicit Real(double a) {
		value = Type(a);
		debug_check(*this);
	}
	Real& operator=(Real const &a) {
		value = a.value;
		debug_check(a);
		debug_check(*this);
		return *this;
	}
	Real operator+() const {
		debug_check(*this);
		return *this;
	}
	Real operator-() const {
		debug_check(*this);
		return Real(-value);
	}
	Real operator+(Real const &a) const {
		debug_check(a);
		debug_check(*this);
		return Real(value + a.value);
	}
	Real operator-(Real const &a) const {
		debug_check(a);
		debug_check(*this);
		return Real(value - a.value);
	}
	Real operator*(Real const &a) const {
		debug_check(a);
		debug_check(*this);
		return Real(value * a.value);
	}
	Real operator/(Real const &a) const {
		Real result;
		zero_check(a);
		debug_check(a);
		debug_check(*this);
		result.value = value / a.value;
		return result;
	}
	Real& operator+=(Real const &a) {
		debug_check(a);
		debug_check(*this);
		*this = *this + a;
		return *this;
	}
	Real& operator-=(Real const &a) {
		debug_check(a);
		debug_check(*this);
		*this = *this - a;
		return *this;
	}
	Real& operator*=(Real const &a) {
		debug_check(a);
		debug_check(*this);
		*this = *this * a;
		return *this;
	}
	Real& operator/=(Real const &a) {
		zero_check(a);
		debug_check(a);
		*this = *this / a;
		return *this;
	}
	bool operator==(Real const &a) const {
		debug_check(a);
		debug_check(*this);
		return value == a.value;
	}
	bool operator!=(Real const &a) const {
		debug_check(a);
		debug_check(*this);
		return value != a.value;
	}
	bool operator<=(Real const &a) const {
		debug_check(a);
		debug_check(*this);
		return value <= a.value;
	}
	bool operator>=(Real const &a) const {
		debug_check(a);
		debug_check(*this);
		return value >= a.value;
	}
	bool operator<(Real const &a) const {
		debug_check(a);
		debug_check(*this);
		return value < a.value;
	}
	bool operator>(Real const &a) const {
		debug_check(a);
		debug_check(*this);
		return value > a.value;
	}
	static Real zero() {
		Real z;
		z.value = Type(0);
		return z;
	}
	friend Real abs(Real a) {
		debug_check(a);
		a.value = std::fabs(a.value);
		return a;
	}
	friend Real exp(Real a) {
		debug_check(a);
		a.value = std::exp(a.value);
		return a;
	}
	friend Real cos(Real a) {
		debug_check(a);
		a.value = std::cos(a.value);
		return a;
	}
	friend Real sin(Real a) {
		debug_check(a);
		a.value = std::sin(a.value);
		return a;
	}
	friend Real sqrt(Real a) {
		nonneg_check(a);
		debug_check(a);
		a.value = std::sqrt(a.value);
		return a;
	}
	friend Real pow(Real a, Real b) {
		nonneg_check(a);
		debug_check(a);
		a.value = std::pow(a.value, b.value);
		return a;
	}
	friend Real max(Real a, Real b) {
		debug_check(a);
		debug_check(b);
		a.value = std::max(a.value, b.value);
		return a;
	}
	friend Real min(Real a, Real b) {
		debug_check(a);
		debug_check(b);
		a.value = std::min(a.value, b.value);
		return a;
	}
	friend std::string to_string(Real r) {
		return std::to_string(r.value);
	}
private:
	Type value;
	static void nonneg_check(Real a) {
		if constexpr (CHECK_REALS) {
			if (a.value < 0.0) {
				std::string errorString = "FATAL ERROR: Illegal operation on negative number.\n";
				errorString += "Stack trace:\n";
				errorString += std::to_string(std::stacktrace::current());
				std::cout << errorString;
				abort();
//				throw std::invalid_argument(errorString);
			}
		}
	}
	static void zero_check(Real a) {
		if constexpr (CHECK_REALS) {
			if (a.value == 0.0) {
				std::string errorString = "FATAL ERROR: Divide by zero\n";
				errorString += "Stack trace:\n";
				errorString += std::to_string(std::stacktrace::current());
				std::cout << errorString;
				abort();
//				throw std::invalid_argument(errorString);
			}
		}
	}
	static void debug_check(Real a) {
		if constexpr (CHECK_REALS) {
			if (!std::isfinite(a.value)) {
				std::string errorString = "FATAL ERROR: Operation on NaN\n";
				errorString += "Stack trace:\n";
				errorString += std::to_string(std::stacktrace::current());
				std::cout << errorString;
				abort();
//				throw std::invalid_argument(errorString);
			}
		}
	}
};

}

constexpr Math::Real operator""_Real(long double a) {
	return Math::Real(a);
}

#endif /* INCLUDE_REAL_HPP_ */
