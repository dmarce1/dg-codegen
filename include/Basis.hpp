#pragma once

#include <cmath>

#include "Legendre.hpp"
#include "TriIndex.hpp"

template<typename T, int O, int D>
struct Basis {
	using index_type = TriIndex<D, O>;
	constexpr auto operator()(std::array<T, D> const &x) const {
		std::array<T, size()> phi;
		phi.fill(T(1));
		for (int d = 0; d < D; d++) {
			auto const Pn = legendreP<T, O>(x[d]);
			for (auto P = index_type::begin(); P != index_type::end(); P++) {
				phi[P] *= Pn[P[d]];
			}
		}
		return phi;
	}
	constexpr auto gradient(int dim, std::array<T, D> const &x) const {
		std::array<T, size()> phi;
		phi.fill(T(1));
		for (int d = 0; d < D; d++) {
			auto const Pn = dMLegendrePdXm<T, O, 1>(x[d]);
			for (auto P = index_type::begin(); P != index_type::end(); P++) {
				phi[P] *= Pn[int(d == dim)][P[d]];
			}
		}
		return phi;
	}
	constexpr auto norm() const {
		std::array<T, size()> N;
		for (auto P = index_type::begin(); P != index_type::end(); P++) {
			T &n = N[P];
			n = T(1);
			for (int d = 0; d < D; d++) {
				n *= T(P[d]) + T(0.5);
			}
		}
		return N;
	}
	static constexpr int size() {
		return index_type::size();
	}
};

