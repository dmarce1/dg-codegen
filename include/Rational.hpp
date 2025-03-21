/*
 * Rational.hpp
 *
 *  Created on: Mar 9, 2025
 *      Author: dmarce1
 */

#ifndef INCLUDE_RATIONAL_HPP_
#define INCLUDE_RATIONAL_HPP_

#include <numeric>

template<typename INT>
struct Rational {
	Rational operator+(Rational const &b) const {
		Rational a;
		a.num = num * b.den + b.num * den;
		a.den = den * b.den;
		a.normalize();
		return a;
	}
	Rational operator-(Rational const &b) const {
		Rational a;
		a.num = num * b.den - b.num * den;
		a.den = den * b.den;
		a.normalize();
		return a;
	}
	Rational operator*(Rational const &b) const {
		Rational a;
		a.num = num * b.num;
		a.den = den * b.den;
		a.normalize();
		return a;
	}
	Rational operator/(Rational const &b) const {
		Rational a;
		a.num = num * b.den;
		a.den = den * b.num;
		a.normalize();
		return a;
	}
	Rational(INT n = 0, INT d = 1) :
			num(n), den(d) {
		normalize();
	}
	Rational(Rational const&) = default;
	Rational(Rational&&) = default;
	Rational& operator=(Rational const&) = default;
	Rational& operator=(Rational&&) = default;
	explicit operator double() const {
		return double(num) / double(den);
	}
	bool operator==(Rational const &other) const {
		return (num == other.num) && (den == other.den);
	}
	bool operator<(Rational const &other) const {
		return num * other.den < den * other.num;
	}
	bool operator!=(Rational const &other) const {
		return !operator==(other);
	}
	bool operator>(Rational const &other) const {
		return other < *this;
	}
	bool operator>=(Rational const &other) const {
		return !(*this < other);
	}
	bool operator<=(Rational const &other) const {
		return !(*this > other);
	}
	std::string toString() const {
		if( den != 1 ) {
			return std::to_string(num) + ".0/" + std::to_string(den) + ".0";
		} else {
			return std::to_string(num);
		}
	}
private:
	void normalize() {
		INT const div = std::gcd(num, den);
		if (div > 1) {
			num /= div;
			den /= div;
		}
	}
	INT num;
	INT den;
};

#endif /* INCLUDE_RATIONAL_HPP_ */
