#pragma once

#include <cmath>
#include <numbers>

#include "Legendre.hpp"
#include "MultiIndex.hpp"
#include "Util.hpp"

template<typename T, int N>
using quad_rc_type = std::tuple<std::array<T, N>, std::array<T, N>>;

template<typename T, int N>
constexpr auto GaussLegendreQuadrature() {
	constexpr T pi = std::numbers::pi_v<T>;
	std::array<T, N> weights;
	std::array<T, N> points;
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
		points[n] = x;
		weights[n] = T(2) / (z * z * (T(1) - x * x));
	}
	return quad_rc_type<T, N>(points, weights);
}

template<typename T, int N>
constexpr auto GaussLobattoQuadrature() {
	constexpr T pi = std::numbers::pi_v<T>;
	std::array<T, N> weights;
	std::array<T, N> points;
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
		points[n] = x;
		weights[n] = T(2) / (T(N * N - N) * z.back() * z.back());
	}
	points[0] = -T(1);
	points[N - 1] = +T(1);
	weights[0] = weights[N - 1] = T(2) / T(N * N - N);
	return quad_rc_type<T, N>(points, weights);
}

template<typename T, int N, int D, quad_rc_type<T, N> (*Q)() = &GaussLegendreQuadrature<T, N>>
struct Quadrature {
	using index_type = MultiIndex<Range<int, D> {repeat<D>(0), repeat<D>(N)}>;
	constexpr T const& weight(int i) const {
		return std::get<1>(points_)[i];
	}
	constexpr auto const& point(int i) const {
		return std::get<0>(points_)[i];
	}
	static constexpr int size() {
		return index_type::size();
	}
private:
	static constexpr auto initialize() {
		auto [points1, weights1] = Q();
		std::array<T, size()> weights;
		std::array<std::array<T, D>, size()> points;
		for (int dim = 0; dim < D; dim++) {
			for (auto I = index_type::begin(); I != index_type::end(); I++) {
				int const pi = I;
				weights[pi] = T(1);
				for (int d = 0; d < D; d++) {
					weights[pi] *= weights1[I[d]];
					points[pi][d] = points1[I[d]];
				}
			}
		}
		return std::tuple(points, weights);
	}
	static constexpr std::tuple<std::array<std::array<T, D>, size()>, std::array<T, size()>> points_ = initialize();
};

template<typename T, int N, quad_rc_type<T, N> (*Q)()>
struct Quadrature<T, N, 0, Q> {
	using index_type = MultiIndex<Range<int, 0> {repeat<0>(0), repeat<0>(1)}>;
	constexpr T const& weight(int i) const {
		return T(1);
	}
	constexpr auto point(int i) const {
		return std::array<T, 0> { };
	}
	static constexpr int size() {
		return 1;
	}
};

