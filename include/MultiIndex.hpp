#pragma once

#include "Range.hpp"
#include <cassert>


template<auto RANGE, auto InnerRange = RANGE>
struct MultiIndex {
	static_assert(InnerRange.ndim == RANGE.ndim);
	static constexpr int D = RANGE.ndim;
	MultiIndex() = default;
	MultiIndex(MultiIndex const&) = default;
	MultiIndex(MultiIndex&&) = default;
	MultiIndex& operator=(MultiIndex const&) = default;
	MultiIndex& operator=(MultiIndex&&) = default;
	constexpr MultiIndex(std::array<int, D> const &I) :
			I_(I) {
	}
	constexpr MultiIndex& operator++() {
		int dim = D - 1;
		while (++I_[dim] == InnerRange.end[dim]) {
			if (dim == 0) {
				I_ = InnerRange.end;
				return *this;
			}
			I_[dim] = InnerRange.begin[dim];
			dim--;
		}
		return *this;
	}
	constexpr MultiIndex operator++(int) const {
		auto const rc = *this;
		const_cast<MultiIndex*>(this)->operator++();
		return rc;
	}
	constexpr operator int() const {
		if constexpr (D) {
			static constexpr auto span = rangeSpan(RANGE);
			int i = I_[0] - RANGE.begin[0];
			for (int d = 1; d < D; d++) {
				i = i * span[d] + I_[d] - RANGE.begin[d];
			}
			return i;
		} else {
			return 0;
		}
	}
	constexpr int operator[](int i) const {
		return I_[i];
	}
	constexpr int& operator[](int i) {
		return I_[i];
	}
	constexpr operator std::array<int, D>() const {
		return I_;
	}
	static constexpr auto begin() {
		return MultiIndex(InnerRange.begin);
	}
	static constexpr auto end() {
		return MultiIndex(InnerRange.end);
	}
	static constexpr int size() {
		return rangeVolume(InnerRange);
	}
	friend std::ostream& operator<<(std::ostream &os, MultiIndex const &I) {
		os << "(";
		for (int d = 0; d < D; d++) {
			os << std::to_string(I.I_[d]);
			if (d + 1 < D) {
				os << ", ";
			}
		}
		os << ")";
		return os;
	}

private:
	std::array<int, D> I_;
};

template<auto RANGE, auto InnerRange = RANGE>
static constexpr auto createMultiIndexMap() {
	constexpr int N = rangeVolume(InnerRange);
	using index_type = MultiIndex<RANGE, InnerRange>;
	std::array<index_type, N> map;
	index_type I = index_type::begin();
	for (int i = 0; i < N; i++) {
		map[i] = I;
		I++;
	}
	return map;

}

template<auto RANGE, typename InputIt, typename OutputIt>
constexpr void reverseMultiIndexData(InputIt srcBegin, InputIt srcEnd, OutputIt dstBegin) {
	using index_type = MultiIndex<RANGE>;
	static constexpr int D = index_type::D;
	static constexpr int N = index_type::size();
	assert(std::distance(srcBegin, srcEnd) == N);
	constexpr auto indexMap = createMultiIndexMap<RANGE>();
	for (int linear = 0; linear < N; ++linear) {
		auto const &value = *(srcBegin + linear);
		index_type idx = indexMap[linear];
		std::array<int, D> rev { };
		for (int d = 0; d < D; ++d) {
			rev[d] = idx[D - 1 - d];
		}
		index_type revIdx { rev };
		*(dstBegin + static_cast<int>(revIdx)) = value;
	}
}
