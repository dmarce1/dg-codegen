/*
 * Interval.hpp
 *
 *  Created on: Jan 13, 2025
 *      Author: dmarce1
 */

#ifndef INCLUDE_INTERVAL_HPP_
#define INCLUDE_INTERVAL_HPP_

#include "Vector.hpp"

template<typename Type, int Ndim>
struct Interval;

template<typename Type, int Ndim>
Interval<Type, Ndim> nullInterval();

template<typename Type, int Ndim>
struct Interval {
	static constexpr Type one = Type(1);
	static constexpr Type two = Type(2);
	constexpr Interval() {
	}
	Interval(Interval const &other) {
		*this = other;
	}
	Interval(Interval &&other) {
		*this = std::move(other);
	}
	Interval(std::array<Type, Ndim> const &b, std::array<Type, Ndim> const &e) {
		_begin = b;
		_end = e;
	}
	Interval& operator=(Interval const &other) {
		_begin = other._begin;
		_end = other._end;
		return *this;
	}
	Interval& operator=(Interval &&other) {
		_begin = std::move(other._begin);
		_end = std::move(other._end);
		return *this;
	}
	Interval operator+() const {
		return *this;
	}
	Interval operator-() const {
		Interval I;
		I.begin = -begin;
		I.end = -end;
		return I;
	}
	Interval operator+(Math::Vector<Type, Ndim> const &v) const {
		Interval I;
		I.begin = begin + v;
		I.end = end + v;
		return I;
	}
	Interval operator-(Math::Vector<Type, Ndim> const &v) const {
		Interval I;
		I.begin = begin - v;
		I.end = end - v;
		return I;
	}
	Interval operator*(Type a) const {
		Interval I;
		for (int k = 0; k < Ndim; k++) {
			I.begin(k) = a * begin(k);
			I.end(k) = a * end(k);
		}
		return I;
	}
	Interval operator/(Type a) const {
		return *this * (one / a);
	}
	Interval& operator+=(Type a) {
		*this = *this + a;
		return *this;
	}
	Interval& operator-=(Type a) {
		*this = *this - a;
		return *this;
	}
	Interval& operator*=(Type a) {
		*this = *this * a;
		return *this;
	}
	Interval& operator/=(Type a) {
		*this = *this / a;
		return *this;
	}
	bool operator==(Interval const &B) const {
		Interval const &A = *this;
		for (int k = 0; k < Ndim; k++) {
			if (A.begin(k) != B.begin(k)) {
				return false;
			}
			if (A.end(k) != B.end(k)) {
				return false;
			}
		}
		return true;
	}
	bool operator!=(Interval const &B) const {
		return !operator==(B);
	}
	Type& begin(int k = 0) {
		return _begin[k];
	}
	Type& end(int k = 0) {
		return _end[k];
	}
	Type begin(int k = 0) const {
		return _begin[k];
	}
	Type end(int k = 0) const {
		return _end[k];
	}
	Type span(int k) const {
		return end(k) - begin(k);
	}
	Type volume() const {
		Type v = one;
		for (int k = 0; k < Ndim; k++) {
			v *= span(k);
		}
		return v;
	}
	int longest() const {
		int n;
		Type L = -1;
		for (int k = 0; k < Ndim; k++) {
			Type const l = span(k);
			if (L < l) {
				L = l;
				n = k;
			}
		}
		return n;
	}
	std::pair<Interval, Interval> split(int k = -1) const {
		std::pair<Interval, Interval> rc;
		if (k >= 0) {
			Type const mid = (begin(k) + end(k)) / two;
			rc.first = *this;
			rc.second = *this;
			rc.first.end(k) = mid;
			rc.second.begin(k) = mid;
		} else {
			rc = split(longest());
		}
		return rc;
	}
	bool intersects(Interval const &B) const {
		Interval const &A = *this;
		for (int k = 0; k < Ndim; k++) {
			Type const Ibegin = max(A.begin(k), B.begin(k));
			Type const Iend = min(A.end(k), B.end(k));
			if (Iend < Ibegin) {
				return false;
			}
		}
		return true;
	}
	Interval intersection(Interval const &B) const {
		Interval const &A = *this;
		Interval I;
		for (int k = 0; k < Ndim; k++) {
			I.begin(k) = max(A.begin(k), B.begin(k));
			I.end(k) = min(A.end(k), B.end(k));
			if (I.end(k) < I.begin(k)) {
				return nullInterval<Type, Ndim>();
			}
		}
		return I;
	}
	friend Interval operator*(Type a, Interval const &B) {
		return B * a;
	}
private:
	Math::Vector<Type, Ndim> _begin;
	Math::Vector<Type, Ndim> _end;
};

template<typename Type, int Ndim>
Interval<Type, Ndim> nullInterval() {
	static constexpr Type zero = Type(0);
	Interval<Type, Ndim> I;
	for (int k = 0; k < Ndim; k++) {
		I.begin(k) = zero;
		I.end(k) = zero;
	}
	return I;
}

template<typename Type, int Ndim>
Interval<Type, Ndim> unitInterval() {
	static constexpr Type zero = Type(0);
	static constexpr Type one = Type(1);
	Interval<Type, Ndim> I;
	for (int k = 0; k < Ndim; k++) {
		I.begin(k) = zero;
		I.end(k) = one;
	}
	return I;
}

#endif /* INCLUDE_INTERVAL_HPP_ */
