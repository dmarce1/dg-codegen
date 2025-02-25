/*
 * Einstein.cpp
 *
 *  Created on: Feb 25, 2025
 *      Author: dmarce1
 */

#include <algorithm>
#include <array>
#include <numeric>
#include <tuple>
#include <stdio.h>

template<int RANK, int DIM>
static constexpr int indices2Index(std::array<int, RANK> indices) {
	return std::accumulate(indices.rbegin(), indices.rend(), 0, [](int i, int j) {
		return DIM * i + j;
	});
}

template<int RANK, int DIM>
static constexpr std::array<int, RANK> index2Indices(int index) {
	std::array<int, RANK> indices;
	for (int i = 0; i < RANK; i++) {
		indices[i] = index % DIM;
		index /= DIM;
	}
	return indices;
}

constexpr int pow(int x, int n) {
	int xm = x;
	int xn = 1;
	while (n > 0) {
		if (n & 1) {
			xn *= xm;
		}
		xm *= xm;
		n >>= 1;
	}
	return xn;
}

template<typename ...ARGS>
struct FirstArgument;

template<>
struct FirstArgument<> {
	using type = void;
};

template<typename FIRST, typename ...ARGS>
struct FirstArgument<FIRST, ARGS...> {
	using type = FIRST;
	static constexpr FIRST operator()(FIRST const &first, ARGS ...args) {
		return first;
	}
};

template<typename F, typename FIRST, typename ...ARGS>
static constexpr auto skipFirst(F const &func, FIRST f1, ARGS ...args) {
	return func(std::forward<ARGS>(args)...);
}

template<int RANK, typename ...ARGS>
static constexpr auto sortIndices(std::array<int, RANK> const &indices, ARGS ...args) {
	FirstArgument<ARGS...> first;
	using array_type = FirstArgument<ARGS...>::type;
	if constexpr (std::is_void<array_type>::value) {
		return indices;
	} else {
		auto const symmetrySet = first(std::forward<ARGS>(args)...);
		auto sortedSymmetries = symmetrySet;
		auto sortedIndices = indices;
		std::sort(sortedSymmetries.begin(), sortedSymmetries.end(), [indices](int a, int b) {
			return indices[a] > indices[b];
		});
		for (size_t i = 0; i < sortedSymmetries.size(); i++) {
			sortedIndices[symmetrySet[i]] = indices[sortedSymmetries[i]];
		}
		if constexpr (sizeof...(ARGS)) {
			auto const function = [sortedIndices]<typename...ARGS2>(ARGS2 ...args) {
				return sortIndices<RANK>(sortedIndices, std::forward<ARGS2>(args)...);
			};
			return skipFirst(function, std::forward<ARGS>(args)...);
		} else {
			return sortedIndices;
		}
	}
}

template<int RANK, int DIM, typename ...ARGS>
static constexpr auto indexMap(ARGS ...args) {
	static constexpr int size = pow(DIM, RANK);
	std::array<int, size> map;
	std::array<bool, size> filled;
	std::fill(filled.begin(), filled.end(), false);
	int nextIndex = 0;
	for (int thisIndex = 0; thisIndex < size; thisIndex++) {
		auto const theseIndices = index2Indices<RANK, DIM>(thisIndex);
		auto const commonIndices = sortIndices<RANK>(theseIndices, std::forward<ARGS>(args)...);
		int const commonIndex = indices2Index<RANK, DIM>(commonIndices);
		printf( "%i %i\n", thisIndex, commonIndex);
		if (!filled[commonIndex]) {
			map[commonIndex] = nextIndex++;
			filled[thisIndex] = true;
		}
		map[thisIndex] = map[commonIndex];
	}
	printf( "\n");
	return map;
}

void testEinstein() {
	constexpr int RANK = 4;
	constexpr int DIM = 3;
	auto I = indexMap<RANK, DIM>(std::array { int(0), int(1) }, std::array { int(2), int(3) });
}
