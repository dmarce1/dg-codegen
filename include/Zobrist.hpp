/*
 * Zobrist.hpp
 *
 *  Created on: Jan 15, 2025
 *      Author: dmarce1
 */

#ifndef INCLUDE_ZOBRIST_HPP_
#define INCLUDE_ZOBRIST_HPP_

#include <limits>
#include <random>
#include <set>

template<class Type = int>
struct ZobristGenerator {
	ZobristGenerator(Type len, Type off = 0) :
			offset(off), length(len), zKeys(len) {
		static constexpr size_t keyMax = std::numeric_limits<size_t>::max() - 1;
		std::set<size_t> used;
		std::mt19937 pseudoRandomGenerator(1234);
		std::uniform_int_distribution<size_t> randomInteger(0, keyMax);
		int i = length;
		size_t key;
		while (i > 0) {
			i--;
			do {
				key = randomInteger(pseudoRandomGenerator);
			} while (!key || (used.find(key) != used.end()));
			zKeys[i] = key;
			used.insert(key);
		}
	}
	auto operator()(auto begin, auto end) const {
		size_t key;
		key ^= key;
		for (auto iterator = begin; iterator != end; iterator++) {
			key ^= zKeys[*iterator];
		}
		return key;
	}
private:
	Type offset;
	Type length;
	std::vector<size_t> zKeys;
};

#endif /* INCLUDE_ZOBRIST_HPP_ */
