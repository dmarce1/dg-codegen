#pragma once

#include <array>
#include <numeric>

#include "Util.hpp"

template<typename Type, int dimensionCount_>
struct Range {
	static constexpr int dimensionCount = dimensionCount_;
	std::array<Type, dimensionCount> begin;
	std::array<Type, dimensionCount> end;
};

template<typename Type, int dimensionCount>
inline constexpr auto rangeSpan(Range<Type, dimensionCount> const &range) {
	std::array<Type, dimensionCount> span;
	for (int d = 0; d < dimensionCount; d++) {
		span[d] = std::max(range.end[d] - range.begin[d], Type(0));
	}
	return span;
}

template<typename Type, int dimensionCount>
inline constexpr Type rangeVolume(Range<Type, dimensionCount> const &range) {
	auto const span = rangeSpan(range);
	return std::accumulate(span.begin(), span.end(), Type(1), std::multiplies<Type>());
}

template<typename TypeA, typename TypeB, int dimensionCountA, auto dimensionCountB>
inline constexpr bool rangeWithin(Range<TypeA, dimensionCountA> const &range,
		std::array<TypeB, dimensionCountB> const &x) {
	static_assert(dimensionCountA == dimensionCountB);
	for (int d = 0; d < dimensionCountA; d++) {
		if ((x[d] < range.begin[d]) || (x[d] >= range.end[d])) {
			return false;
		}
	}
	return true;
}
