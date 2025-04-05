/*
 * Permutation.hpp
 *
 *  Created on: Mar 26, 2025
 *      Author: dmarce1
 */

#ifndef INCLUDE_PERMUTATION_HPP_
#define INCLUDE_PERMUTATION_HPP_

#include <array>
#include <algorithm>
#include <concepts>
#include <cstddef>
#include <numeric>
#include <string>
#include <type_traits>
#include <vector>

#include "IO.hpp"

template<typename T>
concept RandomAccess = requires(T a, size_t i) {
	{	a[i]}-> std::convertible_to<typename T::value_type>;
	{	a.size()}-> std::convertible_to<size_t>;
};

constexpr size_t factorial(size_t n) {
	if (n) {
		return n * factorial(n - 1);
	} else {
		return 1;
	}
}

template<typename Container>
bool lexCompare(Container const &A, Container const &B) {
	size_t const len = std::min(A.size(), B.size());
	for (size_t i = 0; i < len; i++) {
		if (A[i] < B[i]) {
			return true;
		} else if (A[i] > B[i]) {
			return false;
		}
	}
	return A.size() < B.size();
}

template<size_t R>
struct Cycle {
	Cycle() {
		subcycles = std::vector<std::vector<size_t>>();
		for (size_t r1 = 0; r1 < subcycles.size(); r1++) {
			std::vector<size_t> cycle(1, r1 + 1);
			subcycles.push_back(cycle);
		}
	}
	Cycle(std::initializer_list<std::initializer_list<size_t>> list) {
		std::array<bool, R> has = { 0 };
		for (size_t j = 0; j < list.size(); j++) {
			for (auto i : list) {
				std::vector<size_t> cycle;
				for (auto j : i) {
					has[j] = true;
					cycle.push_back(j);
				}
				std::sort(cycle.begin(), cycle.end());
				subcycles.push_back(cycle);
			}
			for (size_t j = 0; j < R; j++) {
				if (!has[j]) {
					std::vector<size_t> cycle(1, j + 1);
					printf("R %i j = %i\n", R, j);
					subcycles.push_back(cycle);
				}
			}
		}
	}
	bool operator<(const Cycle &other) const {
		auto const &A = subcycles;
		auto const &B = other.subcycles;
		size_t const len = std::min(A.size(), B.size());
		for (size_t i = 0; i < len; i++) {
			if (lexCompare(A[i], B[i])) {
				return true;
			} else if (lexCompare(B[i], A[i])) {
				return false;
			}
		}
		return A.size() < B.size();
	}
	bool hasMember(size_t slot) const {
		for (auto const &cycle : subcycles) {
			for (auto const n : cycle) {
				if (slot == n) {
					return true;
				}
			}
		}
		return false;
	}
	size_t subCycleCount() const {
		return subcycles.size();
	}
	std::vector<size_t> const& getSubcycle(size_t i = 0) {
		return subcycles[i];
	}
	auto begin() const {
		return subcycles.begin();
	}
	auto end() const {
		return subcycles.end();
	}
private:
	std::vector<std::vector<size_t>> subcycles;
}
;

template<size_t N>
struct Permutation: public std::array<size_t, N> {
	using base_type = std::array<size_t, N>;
	constexpr static Permutation identity = []() {
		Permutation P;
		std::iota(P.begin(), P.end(), 1);
		return P;
	}();
	constexpr static Permutation reversing = []() {
		Permutation p = identity;
		std::reverse(p.begin(), p.end());
		return p;
	}();
	constexpr Permutation() {
	}
	constexpr Permutation(Permutation const &other) :
			base_type { other.size() } {
		reinterpret_cast<base_type&>(*this) = other;
	}
	constexpr Permutation& operator=(Permutation const &other) {
		reinterpret_cast<base_type&>(*this) = other;
		return *this;
	}
	template<typename IT>
	Permutation(IT itBegin, IT itEnd) {
		size_t i = 0;
		for (auto it = itBegin; it != itEnd; it++) {
			base_type::operator[](i++) = *it;
		}
	}
	constexpr Permutation(std::initializer_list<size_t> const &list) {
		std::copy(list.begin(), list.end(), this->begin());
	}
	Permutation inverse() const {
		base_type const &P = *this;
		Permutation I;
		for (size_t n = 0; n < N; n++) {
			I[P[n] - 1] = n + 1;
		}
		return I;
	}
	template<RandomAccess Container>
	constexpr Container apply(Container const &B) const {
		base_type const &P = *this;
		Container A = B;
		for (size_t n = 0; n < N; n++) {
			A[P[n] - 1] = B[n];
		}
		return A;
	}
	Permutation operator*(Permutation const &B) const {
		base_type const &P = *this;
		Permutation A;
		for (size_t n = 0; n < N; n++) {
			A[P[n] - 1] = B[n];
		}
		return A;
	}
	Permutation operator/(Permutation const &B) const {
		return operator*(B.inverse());
	}
	Permutation& operator*=(Permutation const &B) {
		*this = *this * B;
		return *this;
	}
	Permutation& operator/=(Permutation const &B) {
		*this = *this / B;
		return *this;
	}
	size_t cycleLength() const {
		Permutation Q = identity;
		size_t count = 0;
		do {
			Q = *this * Q;
			count++;
		} while (Q != identity);
		return count;
	}
	size_t inversionCount() const {
		base_type const &P = *this;
		size_t count = 0;
		for (size_t n = 0; n < N; n++) {
			for (size_t m = n + 1; m < N; m++) {
				if (P[n] > P[m]) {
					count++;
				}
			}
		}
		return count;
	}
	std::vector<Permutation> generateSubgroup() const {
		std::vector<Permutation> G;
		Permutation Q = identity;
		do {
			G.emplace_back(Q);
			Q = apply(Q);
		} while (Q != identity);
		return G;
	}
	size_t genHashKey() const {
		static constexpr auto zobristKeys = []() {
			std::array<size_t, N * N> keys;
			size_t key = 42;
			for (auto &thisKey : keys) {
				key = 69069 * key + 1;
				thisKey = key;
			}
			return keys;
		}();
		size_t key = 0;
		for (size_t n = 0; n < N; n++) {
			key ^= zobristKeys[n * N + (*this)[n] - 1];
		}
		return key;
	}
	Permutation next() const {
		Permutation nextP = *this;
		if (!std::next_permutation(nextP.begin(), nextP.end())) {
			nextP = identity;
		}
		return nextP;
	}
	Permutation previous() const {
		Permutation prevP = *this;
		if (!std::prev_permutation(prevP.begin(), prevP.end())) {
			prevP = reversing;
		}
		return prevP;
	}
	int parity() const {
		return (inversionCount() & 1) ? -1 : +1;
	}
	template<size_t M>
	Permutation<N - M> reduceRank(std::integral_constant<size_t, M> = std::integral_constant<size_t, M>()) {
		Permutation<N> Q = *this;
		for (size_t m = 0; m < M; m++) {
			size_t const q = Q[m];
			std::transform(Q.begin(), Q.end(), Q.begin(), [q](size_t value) {
				if (value > q) {
					value--;
				}
				return value;
			});
		}
		return Permutation<N - M>(Q.begin() + M, Q.end());
	}
	std::vector<Cycle<N>> cycles() const {
		std::vector<Cycle<N>> C;
		std::vector<bool> visited(N, false);
		Cycle<N> thisCycle;
		for (size_t r = 0; r < N; r++) {
			if (!visited[r]) {
				visited[r] = true;
				std::vector<size_t> cycle(1, r + 1);
				size_t next = (*this)[r];
				while (next != r + 1) {
					cycle.push_back(next);
					visited[next - 1] = true;
					next = (*this)[next - 1];
				}
				std::rotate(cycle.begin(), std::max_element(cycle.begin(), cycle.end()), cycle.end());
				thisCycle.addSubcycle(cycle);
			}
			C.push_back(thisCycle);
		}
		std::sort(C.begin(), C.end());
		return C;
	}
};

template<size_t N>
struct FactorialNumber {
	FactorialNumber() = default;
	FactorialNumber(FactorialNumber const &other) {
		other = digits.other;
	}
	FactorialNumber(size_t number) {
		digits[0] = 0;
		for (size_t n = 1; n < N; n++) {
			digits[n] = number % (n + 1);
			number /= (n + 1);
		}
	}
	FactorialNumber(Permutation<N> P) {
		size_t index;
		std::array<size_t, N> theseDigits = digits;
		for (int n = N - 1; n >= 0; n--) {
			digits[n] = P[n];
			for (int m = n - 1; m >= 0; m--) {
				if (P[m] > P[n]) {
					P[m]--;
				}
			}
		}
	}
	operator Permutation<N>() const {
		Permutation<N> P;
		std::array<size_t, N> theseDigits = digits;
		for (int n = N - 1; n >= 0; n--) {
			P[n] = theseDigits[n];
			for (int m = n + 1; m < N; m++) {
				if (theseDigits[m] > P[n]) {
					theseDigits--;
				}
			}
		}
		return P;
	}
	operator size_t() const {
		size_t number = 0;
		size_t base = 1;
		for (size_t n = 1; n < N; n++) {
			base *= (n + 1);
			number += digits[n] * base;
		}
		return number;
	}
private:
	std::array<size_t, N> digits;
};

namespace std {

template<size_t N>
struct hash<Permutation<N>> {
	size_t operator()(Permutation<N> const &P) const {
		return P.genHashKey();
	}
};

}
#endif /* INCLUDE_PERMUTATION_HPP_ */
