#pragma once

#include <array>
#include <numeric>

#include "Util.hpp"

template<typename T, int D>
struct Range {
	static constexpr int ndim = D;
	std::array<T, D> begin;
	std::array<T, D> end;
};

template<typename T, int D>
inline constexpr auto rangeSpan(Range<T, D> const &range) {
	std::array<T, D> span;
	for (int d = 0; d < D; d++) {
		span[d] = std::max(range.end[d] - range.begin[d], T(0));
	}
	return span;
}

template<typename T, int D>
inline constexpr T rangeVolume(Range<T, D> const &range) {
	auto const span = rangeSpan(range);
	return std::accumulate(span.begin(), span.end(), T(1), std::multiplies<T>());
}
