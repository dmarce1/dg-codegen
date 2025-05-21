#pragma once

#include "Range.hpp"

template<auto OuterRange, auto InnerRange = OuterRange>
struct MultiIndex {
	static_assert(InnerRange.ndim == OuterRange.ndim);
	static constexpr int D = OuterRange.ndim;
	MultiIndex() = default;
	MultiIndex(MultiIndex const&) = default;
	MultiIndex(MultiIndex&&) = default;
	MultiIndex& operator=(MultiIndex const&) = default;
	MultiIndex& operator=(MultiIndex&&) = default;
	MultiIndex(std::array<int, D> const &I) :
			I_(I) {
	}
	MultiIndex& operator++() {
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
	MultiIndex operator++(int) const {
		auto const rc = *this;
		const_cast<MultiIndex*>(this)->operator++();
		return rc;
	}
	constexpr operator int() const {
		if constexpr (D) {
			static constexpr auto span = rangeSpan(OuterRange);
			int i = I_[0] - OuterRange.begin[0];
			for (int d = 1; d < D; d++) {
				i = i * span[d] + I_[d] - OuterRange.begin[d];
			}
			return i;
		} else {
			return 0;
		}
	}
	constexpr int operator[](int i) const {
		return I_[i];
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

template<auto OuterRange, auto InnerRange = OuterRange>
static constexpr auto createMultiIndexMap() {
	constexpr int D = OuterRange.D;
	constexpr int N = rangeVolume(InnerRange);
	using index_type = MultiIndex<OuterRange, InnerRange>;
	std::array<index_type, N> map;
	index_type I = index_type::begin();
	for (int i = 0; i < N; i++) {
		map[i] = I;
		I++;
	}
	return map;

}
