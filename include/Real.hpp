#pragma once

#ifndef NDEBUG

#include <cassert>
#include <cmath>
#include <iostream>
#include <limits>
#include <numbers>
#include <stacktrace>
#include <stdexcept>
#include <type_traits>

struct Real {
	using Type = double;
	inline constexpr Real() {
	}
	inline constexpr explicit Real(double a) {
		value = Type(a);
	}
	inline constexpr Real& operator=(Real const &a) {
		value = a.value;
		return *this;
	}
	inline constexpr operator Type() const {
		return value;
	}
	inline constexpr Real operator+() const {
		debug_check(*this);
		return *this;
	}
	inline constexpr Real operator-() const {
		debug_check(*this);
		return Real(-value);
	}
	inline constexpr Real operator+(Real const &a) const {
		debug_check(a);
		debug_check(*this);
		return Real(value + a.value);
	}
	inline constexpr Real operator-(Real const &a) const {
		debug_check(a);
		debug_check(*this);
		return Real(value - a.value);
	}
	inline constexpr Real operator*(Real const &a) const {
		debug_check(a);
		debug_check(*this);
		return Real(value * a.value);
	}
	inline constexpr Real operator/(Real const &a) const {
		Real result;
		zero_check(a);
		debug_check(a);
		debug_check(*this);
		result.value = value / a.value;
		return result;
	}
	inline constexpr Real& operator+=(Real const &a) {
		debug_check(a);
		debug_check(*this);
		*this = *this + a;
		return *this;
	}
	inline constexpr Real& operator-=(Real const &a) {
		debug_check(a);
		debug_check(*this);
		*this = *this - a;
		return *this;
	}
	inline constexpr Real& operator*=(Real const &a) {
		debug_check(a);
		debug_check(*this);
		*this = *this * a;
		return *this;
	}
	inline constexpr Real& operator/=(Real const &a) {
		zero_check(a);
		debug_check(a);
		*this = *this / a;
		return *this;
	}
	inline constexpr bool operator==(Real const &a) const {
		debug_check(a);
		debug_check(*this);
		return value == a.value;
	}
	inline constexpr bool operator!=(Real const &a) const {
		debug_check(a);
		debug_check(*this);
		return value != a.value;
	}
	inline constexpr bool operator<=(Real const &a) const {
		debug_check(a);
		debug_check(*this);
		return value <= a.value;
	}
	inline constexpr bool operator>=(Real const &a) const {
		debug_check(a);
		debug_check(*this);
		return value >= a.value;
	}
	inline constexpr bool operator<(Real const &a) const {
		debug_check(a);
		debug_check(*this);
		return value < a.value;
	}
	inline constexpr bool operator>(Real const &a) const {
		debug_check(a);
		debug_check(*this);
		return value > a.value;
	}
	static inline constexpr Real zero() {
		Real z;
		z.value = Type(0);
		return z;
	}
	static inline constexpr Real tiny() {
		Real z;
		z.value = Type(std::numeric_limits<double>::min());
		return z;
	}
	friend inline constexpr Real abs(Real a) {
		debug_check(a);
		a.value = std::fabs(a.value);
		return a;
	}
	friend inline constexpr Real copysign(Real a, Real b) {
		debug_check(a);
		debug_check(b);
		a.value = std::copysign(a.value, b.value);
		return a;
	}
	friend inline constexpr Real expm1(Real a) {
		debug_check(a);
		a.value = std::expm1(a.value);
		return a;
	}
	friend inline constexpr Real exp(Real a) {
		debug_check(a);
		a.value = std::exp(a.value);
		return a;
	}
	friend inline constexpr Real log(Real a) {
		debug_check(a);
		a.value = std::log(a.value);
		return a;
	}
	friend inline constexpr Real cos(Real a) {
		debug_check(a);
		a.value = std::cos(a.value);
		return a;
	}
	friend inline constexpr Real sin(Real a) {
		debug_check(a);
		a.value = std::sin(a.value);
		return a;
	}
	friend inline constexpr Real acos(Real a) {
		debug_check(a);
		range_check(-Real(1), a, Real(1));
		a.value = std::acos(a.value);
		return a;
	}
	friend inline constexpr Real asin(Real a) {
		debug_check(a);
		range_check(-Real(1), a, Real(1));
		a.value = std::asin(a.value);
		return a;
	}
	friend inline constexpr auto lround(Real a) {
		debug_check(a);
		return std::lround(a.value);
	}
	friend inline constexpr Real tgamma(Real a) {
		nonneg_check(a);
		debug_check(a);
		a.value = std::tgamma(a.value);
		return a;
	}
	friend inline constexpr Real sqrt(Real a) {
		nonneg_check(a);
		debug_check(a);
		a.value = std::sqrt(a.value);
		return a;
	}
	friend inline constexpr Real pow(Real a, Real b) {
		nonneg_check(a);
		debug_check(a);
		debug_check(b);
		a.value = std::pow(a.value, b.value);
		return a;
	}
	friend inline constexpr Real max(Real a, Real b) {
		debug_check(a);
		debug_check(b);
		a.value = std::max(a.value, b.value);
		return a;
	}
	friend inline constexpr Real min(Real a, Real b) {
		debug_check(a);
		debug_check(b);
		a.value = std::min(a.value, b.value);
		return a;
	}
	friend std::string to_string(Real r) {
		return std::to_string(r.value);
	}
	Type value = std::numeric_limits<Type>::signaling_NaN();
	static inline constexpr void nonneg_check(Real a) {
		if (a.value < 0.0) {
			std::string errorString = "FATAL ERROR: Illegal operation on negative number.\n";
			errorString += "Stack trace:\n";
			errorString += std::to_string(std::stacktrace::current());
			std::cout << errorString;
			assert(false);
			throw std::invalid_argument(errorString);
		}
	}
	static inline constexpr void zero_check(Real a) {
		if (a.value == 0.0) {
			std::string errorString = "FATAL ERROR: Divide by zero\n";
			errorString += "Stack trace:\n";
			errorString += std::to_string(std::stacktrace::current());
			std::cout << errorString;
			assert(false);
			throw std::invalid_argument(errorString);
		}
	}
	static inline constexpr void debug_check(Real a) {
		if (!std::isfinite(a.value)) {
			std::string errorString = "FATAL ERROR: Operation on NaN\n";
			errorString += "Stack trace:\n";
			errorString += std::to_string(std::stacktrace::current());
			std::cout << errorString;
			assert(false);
			throw std::invalid_argument(errorString);
		}
	}
	static inline constexpr void range_check(Real a, Real b, Real c) {
		if ((b < a) || (b > c)) {
			std::string errorString = "FATAL ERROR: Range violation\n";
			errorString += "Stack trace:\n";
			errorString += std::to_string(std::stacktrace::current());
			std::cout << errorString;
			assert(false);
		throw std::invalid_argument(errorString);
		}
	}
};

namespace std {
namespace numbers {
template<>
constexpr Real pi_v<Real> = Real(pi_v<typename Real::Type>);
}
}

#else

using Real = double;

#endif
