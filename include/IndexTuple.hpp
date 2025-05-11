/*
 * IndexTuple.hpp
 *
 *  Created on: Mar 27, 2025
 *      Author: dmarce1
 */

#ifndef INCLUDE_INDEXTUPLE_HPP_
#define INCLUDE_INDEXTUPLE_HPP_

#include <array>
#include <cassert>
#include <cstddef>
#include <iostream>

#include "Numbers.hpp"

template<typename T, int P>
std::ostream& operator<<(std::ostream &os, const std::array<T, P> &arr) {
	os << "(";
	for (int i = 0; i < P; ++i) {
		os << arr[i];
		if (i + 1 < P) {
			os << ", ";
		}
	}
	os << ")";
	return os;
}

template<int R>
constexpr bool isNonIncreasing(std::array<int, R> const &indices) {
	for (int i = 1; i < R; ++i) {
		if (indices[i] > indices[i - 1]) {
			std::cout << indices << std::endl;
			return false;
		}
	}
	return true;
}

template<int D, int R>
struct TriangularIndex: public std::array<int, R> {
	using value_type = int;
	constexpr TriangularIndex() = default;
	constexpr TriangularIndex(std::initializer_list<int> init) {
		std::copy(init.begin(), init.end(), indices.begin());
	}
	constexpr TriangularIndex(std::array<int, R> const &init) {
		std::copy(init.begin(), init.end(), indices.begin());
	}
	constexpr TriangularIndex& operator++() {
		int idx = R - 1;
		while (true) {
			if (indices[idx] >= (!idx ? (D - 1) : std::min(indices[idx - 1], D - 1))) {
				if (idx == 0) {
					std::fill(indices.begin(), indices.end(), D);
					return *this;
				} else {
					indices[idx] = 0;
					idx--;
					continue;
				}
			} else {
				indices[idx]++;
				return *this;
			}
		}
	}
	constexpr int& operator[](int k) {
		return indices[k];
	}
	constexpr TriangularIndex operator++(int) const {
		auto const temp = *this;
		const_cast<TriangularIndex*>(this)->operator++();
		return temp;
	}
	constexpr int const& operator[](int k) const {
		return indices[k];
	}
	constexpr int flatIndex() const {
		int k = 0;
		for (int i = 0; i < R; i++) {
			int const ci = indices[i];
			k += Math::binCo<int>(ci + R - i - 1, R - i);
		}
		return k;
	}
	static constexpr int elementCount() {
		return Math::binCo<int>(D + R - 1, R);
	}
	constexpr int size() const {
		return R;
	}
	static constexpr TriangularIndex begin() {
		TriangularIndex<D, R> b;
		std::fill(b.indices.begin(), b.indices.end(), 0);
		return b;
	}
	static constexpr TriangularIndex end() {
		TriangularIndex<D, R> e;
		std::fill(e.indices.begin(), e.indices.end(), D);
		return e;
	}
	constexpr operator std::array<int, R>() const {
		return *this;
	}
	constexpr bool operator==(TriangularIndex const &other) const {
		for (int k = 0; k < R; k++) {
			if (indices[k] != other[k]) {
				return false;
			}
		}
		return true;
	}
	constexpr bool operator<(TriangularIndex const &other) const {
		for (int k = 0; k < R; k++) {
			if (indices[k] < other[k]) {
				return true;
			} else if (indices[k] > other[k]) {
				return false;
			}
		}
		return false;
	}
	constexpr bool operator!=(TriangularIndex const &other) const {
		return !(*this == other);
	}
	constexpr bool operator<=(TriangularIndex const &other) const {
		return !(*this < other);
	}
	constexpr bool operator>(TriangularIndex const &other) const {
		return other < *this;
	}
	constexpr bool operator>=(TriangularIndex const &other) const {
		return !(other < *this);
	}
	explicit constexpr operator int() const {
		return flatIndex();
	}
private:
	std::array<int, R> indices;
};

template<int D, int R>
struct IndexTuple {
	using value_type = int;
	constexpr IndexTuple(std::initializer_list<int> init) {
		std::copy(init.begin(), init.end(), indices.begin());
	}
	constexpr IndexTuple(std::array<int, R> const &init) {
		std::copy(init.begin(), init.end(), indices.begin());
	}
	constexpr IndexTuple& operator++() {
		int k = R;
		while (++indices[--k] == D) {
			indices[k] = 0;
			if (!k) {
				*this = end();
				break;
			}
		}
		return *this;
	}
	constexpr int& operator[](int k) {
		return indices[k];
	}
	constexpr IndexTuple operator++(int) const {
		auto const temp = *this;
		const_cast<IndexTuple*>(this)->operator++();
		return temp;
	}
	constexpr int const& operator[](int k) const {
		return indices[k];
	}
	constexpr int flatIndex() const {
		int iFlat = indices[0];
		for (int i = 1; i < R; i++) {
			iFlat = D * iFlat + indices[i];
		}
		return iFlat;
	}
	static constexpr int elementCount() {
		int count = 1;
		int placeValue = D;
		int i = R;
		while (i) {
			if (i & 1) {
				count *= placeValue;
			}
			placeValue *= placeValue;
			i >>= 1;
		}
		return count;
	}
	constexpr int size() const {
		return R;
	}
	constexpr IndexTuple() {
	}
	static constexpr IndexTuple begin() {
		IndexTuple<D, R> b;
		std::fill(b.indices.begin(), b.indices.end(), 0);
		return b;
	}
	static constexpr IndexTuple end() {
		IndexTuple e = begin();
		e.indices[0] = D;
		return e;
	}
	operator std::array<int, R>() const {
		return indices;
	}
	constexpr bool operator==(IndexTuple const &other) const {
		for (int k = 0; k < R; k++) {
			if (indices[k] != other[k]) {
				return false;
			}
		}
		return true;
	}
	constexpr bool operator<(IndexTuple const &other) const {
		for (int k = 0; k < R; k++) {
			if (indices[k] < other[k]) {
				return true;
			} else if (indices[k] > other[k]) {
				return false;
			}
		}
		return false;
	}
	constexpr bool operator!=(IndexTuple const &other) const {
		return !(*this == other);
	}
	constexpr bool operator<=(IndexTuple const &other) const {
		return !(*this < other);
	}
	constexpr bool operator>(IndexTuple const &other) const {
		return other < *this;
	}
	constexpr bool operator>=(IndexTuple const &other) const {
		return !(other < *this);
	}
private:
	std::array<int, R> indices;
};

#endif /* INCLUDE_INDEXTUPLE_HPP_ */
