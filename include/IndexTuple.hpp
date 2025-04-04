/*
 * IndexTuple.hpp
 *
 *  Created on: Mar 27, 2025
 *      Author: dmarce1
 */

#ifndef INCLUDE_INDEXTUPLE_HPP_
#define INCLUDE_INDEXTUPLE_HPP_

#include <array>
#include <cstddef>
#include <iostream>

template<size_t D, size_t R>
struct IndexTuple {
	using base_type = std::array<size_t, R>;
	using value_type = size_t;
	constexpr IndexTuple(std::initializer_list<size_t> init) {
		std::copy(init.begin(), init.end(), indices.begin());
	}
	constexpr IndexTuple& operator++() {
		size_t k = R;
		while (++indices[--k] == D) {
			indices[k] = 0;
			if (!k) {
				*this = end();
				break;
			}
		}
		return *this;
	}
	constexpr size_t& operator[](size_t k) {
		return indices[k];
	}
	constexpr IndexTuple operator++(int) const {
		auto const temp = *this;
		const_cast<IndexTuple*>(this)->operator++();
		return temp;
	}
	constexpr size_t const& operator[](size_t k) const {
		return indices[k];
	}
	constexpr size_t flatIndex() const {
		size_t iFlat = indices[0];
		for (size_t i = 1; i < R; i++) {
			iFlat = D * iFlat + indices[i];
		}
		return iFlat;
	}
	constexpr bool operator==(IndexTuple const &other) const {
		for (size_t k = 0; k < R; k++) {
			if (indices[k] != other[k]) {
				return false;
			}
		}
		return true;
	}
	constexpr bool operator<(IndexTuple const &other) const {
		for (size_t k = 0; k < R; k++) {
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
	constexpr 	size_t elementCount() const {
		size_t count = 1;
		size_t placeValue = D;
		size_t i = R;
		while (i) {
			if (i & 1) {
				count *= placeValue;
			}
			placeValue *= placeValue;
			i >>= 1;
		}
		return count;
	}
	constexpr size_t size() const {
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
private:
	std::array<size_t, R> indices;
};

#endif /* INCLUDE_INDEXTUPLE_HPP_ */
