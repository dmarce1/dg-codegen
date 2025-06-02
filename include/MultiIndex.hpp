#pragma once

#include "ContainerArithmetic.hpp"
#include "Range.hpp"
#include <cassert>

template<auto firstArgument, auto secondArgument = firstArgument>
struct MultiIndex {
	static constexpr auto initializeDimensionCount() {
		if constexpr (std::is_integral<decltype(firstArgument)>::value) {
			return secondArgument;
		} else {
			return firstArgument.dimensionCount;
		}
	}
	static constexpr int dimensionCount = initializeDimensionCount();
	static constexpr auto initializeInteriorRange() {
		if constexpr (std::is_integral<decltype(firstArgument)>::value) {
			return Range<int, dimensionCount> {
					repeat<dimensionCount>(0),
					repeat<dimensionCount>(firstArgument) };
		} else {
			return secondArgument;
		}
	}
	static constexpr auto initializeExteriorRange() {
		if constexpr (std::is_integral<decltype(firstArgument)>::value) {
			return Range<int, dimensionCount> {
					repeat<dimensionCount>(0),
					repeat<dimensionCount>(firstArgument) };
		} else {
			return firstArgument;
		}
	}
	static constexpr Range<int, dimensionCount> interiorRange = initializeInteriorRange();
	static constexpr Range<int, dimensionCount> exteriorRange = initializeExteriorRange();
	MultiIndex() = default;
	MultiIndex(MultiIndex const&) = default;
	MultiIndex(MultiIndex&&) = default;
	MultiIndex& operator=(MultiIndex const&) = default;
	MultiIndex& operator=(MultiIndex&&) = default;
	constexpr MultiIndex(std::array<int, dimensionCount> const &other) :
			indexValues_(other) {
	}
	constexpr MultiIndex& operator++() {
		int dim = dimensionCount - 1;
		while (++indexValues_[dim] == interiorRange.end[dim]) {
			if (dim == 0) {
				indexValues_ = interiorRange.end;
				return *this;
			}
			indexValues_[dim] = interiorRange.begin[dim];
			dim--;
		}
		return *this;
	}
	static constexpr std::array<int, dimensionCount> strides() {
		std::array<int, dimensionCount> strides;
		int stride = 1;
		for( int dimensionIndex = 0; dimensionIndex < dimensionCount; dimensionIndex++) {
			strides[dimensionIndex] = stride;
			stride *= exteriorRange.end[dimensionIndex] - exteriorRange.begin[dimensionIndex];
		}
		return strides;
	}
	constexpr MultiIndex operator++(int) const {
		auto const rc = *this;
		const_cast<MultiIndex*>(this)->operator++();
		return rc;
	}
	constexpr operator int() const {
		if constexpr (dimensionCount) {
			static constexpr auto span = rangeSpan(exteriorRange);
			int i = indexValues_[0] - exteriorRange.begin[0];
			for (int d = 1; d < dimensionCount; d++) {
				i = i * span[d] + indexValues_[d] - exteriorRange.begin[d];
			}
			return i;
		} else {
			return 0;
		}
	}
	constexpr int operator[](int i) const {
		return indexValues_[i];
	}
	constexpr int& operator[](int i) {
		return indexValues_[i];
	}
	constexpr operator std::array<int, dimensionCount>() const {
		return indexValues_;
	}
	static constexpr auto begin() {
		return MultiIndex(interiorRange.begin);
	}
	static constexpr auto end() {
		return MultiIndex(interiorRange.end);
	}
	static constexpr int count() {
		return rangeVolume(interiorRange);
	}
	static constexpr int size() {
		return dimensionCount;
	}
	friend std::ostream& operator<<(std::ostream &os, MultiIndex const &I) {
		os << "(";
		for (int d = 0; d < dimensionCount; d++) {
			os << std::to_string(I.indexValues_[d]);
			if (d + 1 < dimensionCount) {
				os << ", ";
			}
		}
		os << ")";
		return os;
	}
private:
	std::array<int, dimensionCount> indexValues_;
};

template<auto outerRange, auto innerRange>
struct CanDoArithmetic<MultiIndex<outerRange, innerRange>> {
	static constexpr bool value = true;
};

template<auto outerRange, auto innerRange = outerRange>
constexpr auto createMultiIndexMap() {
	constexpr int N = rangeVolume(innerRange);
	using index_type = MultiIndex<outerRange, innerRange>;
	std::array<index_type, N> map;
	index_type I = index_type::begin();
	for (int i = 0; i < N; i++) {
		map[i] = I;
		I++;
	}
	return map;
}

template<auto outerRange, typename InputIt, typename OutputIt>
constexpr void reverseMultiIndexData(InputIt srcBegin, InputIt srcEnd, OutputIt dstBegin) {
	using index_type = MultiIndex<outerRange>;
	static constexpr int D = index_type::dimensionCount;
	static constexpr int N = index_type::count();
	assert(std::distance(srcBegin, srcEnd) == N);
	constexpr auto indexMap = createMultiIndexMap<outerRange>();
	for (int linear = 0; linear < N; ++linear) {
		auto const &value = *(srcBegin + linear);
		index_type idx = indexMap[linear];
		std::array<int, D> rev { };
		for (int d = 0; d < D; ++d) {
			rev[d] = idx[D - 1 - d];
		}
		index_type revIdx {
				rev };
		*(dstBegin + static_cast<int>(revIdx)) = value;
	}
}
