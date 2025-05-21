#pragma once

#include <cmath>
#include <numbers>

#include "Legendre.hpp"
#include "MultiIndex.hpp"
#include "Util.hpp"

template<typename T, int N>
struct GaussLegendreQuadrature {
	using T_ = T;
	static constexpr int N_ = N;
	GaussLegendreQuadrature() {
		constexpr T pi = std::numbers::pi_v<T>;
		for (int n = 0; n < N; n++) {
			std::array<std::array<T, N + 1>, 2> Pn;
			T theta, x;
			T newTheta = pi * (T(1) - (T(n) + T(0.5)) / T(N));
			do {
				theta = newTheta;
				x = cos(theta);
				Pn = dMLegendrePdXm<T, N + 1, 1>(x);
				newTheta = theta + Pn[0].back() / (sin(theta) * Pn[1].back());
			} while (theta != newTheta);
			T const z = Pn[1].back();
			points_[n] = x;
			weights_[n] = T(2) / (z * z * (T(1) - x * x));
		}
	}
	constexpr T const& weight(int i) const {
		return weights_[i];
	}
	constexpr T const& point(int i) const {
		return points_[i];
	}
private:
	std::array<T, N> weights_;
	std::array<T, N> points_;
};

template<typename T, int N>
struct GaussLobattoQuadrature {
	using T_ = T;
	static constexpr int N_ = N;
	GaussLobattoQuadrature() {
		constexpr T pi = std::numbers::pi_v<T>;
		for (int n = 1; n < N - 1; n++) {
			T theta, x, y;
			T newTheta = pi * (T(1) - (T(n) + T(0.5)) / T(N));
			do {
				theta = newTheta;
				x = cos(theta);
				y = sin(theta);
				auto const dPnDx = legendreP<T, N, 1>(x);
				auto const d2PnDx2 = legendreP<T, N, 2>(x);
				newTheta = theta + dPnDx.back() * y / (y * y * d2PnDx2.back() - dPnDx.back() * x);
			} while (newTheta != theta);
			std::array<T, N> const z = legendreP<T, N>(x);
			points_[n] = x;
			weights_[n] = T(2) / (T(N * N - N) * z.back() * z.back());
		}
		points_[0] = -T(1);
		points_[N - 1] = +T(1);
		weights_[0] = weights_[N - 1] = T(2) / T(N * N - N);
	}
	constexpr T const& weight(int i) const {
		return weights_[i];
	}
	constexpr T const& point(int i) const {
		return points_[i];
	}
private:
	std::array<T, N> weights_;
	std::array<T, N> points_;
};

template<typename Q, int D>
struct Quadrature {
	using T = typename Q::T_;
	static constexpr int N = Q::N_;
	using index_type = MultiIndex<Range<int, D> {repeat<D>(0), repeat<D>(N)}>;
	Quadrature() {
		static Q const quad { };
		for (int dim = 0; dim < D; dim++) {
			for (auto I = index_type::begin(); I != index_type::end(); I++) {
				int const pi = I;
				weights_[pi] = T(1);
				for (int d = 0; d < D; d++) {
					weights_[pi] *= quad.weight(I[d]);
					points_[pi][d] = quad.point(I[d]);
				}
			}
		}
	}
	constexpr T const& weight(int i) const {
		return weights_[i];
	}
	constexpr auto const& point(int i) const {
		return points_[i];
	}
	static constexpr int size() {
		return index_type::size();
	}
private:
	std::array<T, size()> weights_;
	std::array<std::array<T, D>, size()> points_;
};

template<typename Q>
struct Quadrature<Q, 0> {
	using T = typename Q::T_;
	using index_type = MultiIndex<Range<int, 0> {repeat<0>(0), repeat<0>(1)}>;
	constexpr T const& weight(int i) const {
		return T(1);
	}
	constexpr auto point(int i) const {
		return std::array<T, 0> {};
	}
	static constexpr int size() {
		return 1;
	}
};

