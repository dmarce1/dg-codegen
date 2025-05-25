#pragma once

#include "Util.hpp"

template<int O, int D>
struct TriIndex {
	TriIndex() = default;
	TriIndex(TriIndex const&) = default;
	TriIndex(TriIndex&&) = default;
	TriIndex& operator=(TriIndex const&) = default;
	TriIndex& operator=(TriIndex&&) = default;
	TriIndex(std::array<int, D> const &I) :
			I_(I), deg_(std::accumulate(I_.begin(), I_.end(), 0)) {
	}
	TriIndex& operator++() {
		while (true) {
			for (int i = 1; i < D; i++) {
				if (I_[i] > 0) {
					I_[i]--;
					int sum = 1;
					for (int j = i - 1; j >= 0; j--) {
						sum += I_[j];
						I_[j] = 0;
					}
					I_[i - 1] = sum;
					return *this;
				}
			}
			deg_++;
			if (deg_ == O) {
				*this = end();
				return *this;
			}
			I_.fill(0);
			I_.back() = deg_;
			return *this;
		}
	}
	TriIndex operator++(int) const {
		auto const rc = *this;
		const_cast<TriIndex*>(this)->operator++();
		return rc;
	}
	TriIndex inc(int d) {
		TriIndex rc = *this;
		rc.I_[d]++;
		rc.deg_++;
		return rc;
	}
	TriIndex dec(int d) {
		TriIndex rc = *this;
		rc.I_[d]--;
		rc.deg_--;
		return rc;
	}
	constexpr operator int() const {
		int i = 0, j = 0;
		for (int d = 0; d < D; d++) {
			j += I_[d];
			i += binco(j + d, d + 1);
		}
		return i;
	}
	constexpr int operator[](int i) const {
		return I_[i];
	}
	constexpr int degree() const {
		return deg_;
	}
	constexpr bool operator==(TriIndex const &other) const {
		return I_ == other.I_;
	}
	static constexpr auto begin(int deg) {
		TriIndex I;
		I.deg_ = deg;
		I.I_.fill(0);
		I.I_.back() = deg;
		return I;
	}
	static constexpr auto end(int deg) {
		if (deg + 1 == O) {
			return end();
		}
		return begin(deg + 1);
	}
	static auto begin() {
		return TriIndex(repeat<D>(0));
	}
	static auto end() {
		return TriIndex(repeat<D>(O));
	}
	static constexpr int count() {
		return binco(O + D - 1, D);
	}
	static constexpr auto size() {
		return D;
	}
private:
	std::array<int, D> I_;
	int deg_;
};

template<int D, int O, int Deg = -1>
static constexpr auto createTriIndexMap() {
	using index_type = TriIndex<D, O>;
	if constexpr (Deg < 0) {
		std::array<index_type, index_type::size()> map;
		index_type I = index_type::begin();
		for (int i = 0; i < index_type::size(); i++) {
			map[i] = I;
			I++;
		}
	}
	if constexpr (Deg == 0) {
		std::array<index_type, 1> map;
		map[0].fill(0);
		return map;
	}
	constexpr int size = index_type::size() - TriIndex<D, O - 1>::size();
	std::array<index_type, size> map;
	index_type I = index_type::begin(Deg);
	for (int i = 0; i < size; i++) {
		map[i] = I;
		I++;
	}
	return map;
}

template<int D, int O>
std::ostream& operator<<(std::ostream &os, TriIndex<D, O> const &I) {
	os << std::to_string(I.degree());
	os << ":(";
	for (int i = 0; i < D; i++) {
		os << std::to_string(I[i]);
		if (i + 1 < D) {
			os << ", ";
		}
	}
	os << ")";
	return os;
}
