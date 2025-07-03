#pragma once

#include "Matrix.hpp"

#include <array>

template<typename T, int S, T cflCondition, std::array<std::array<T, S>, S> A, std::array<T, S> B, std::array<T, S> C>
struct ButcherTable {
	static constexpr T a(int i, int j) {
		return A[i][j];
	}

	static constexpr T b(int k) {
		return B[k];
	}

	static constexpr T c(int k) {
		return C[k];
	}

	static constexpr int stageCount() {
		return S;
	}

	static constexpr T cfl() {
		return cflCondition;
	}
};

template<typename T>
using RK_1_1 = ButcherTable<T, 1, T(1),
	std::array<std::array<T, 1>, 1> {
		{
			std::array<T, 1> { {T(0)}}
		}
	},
	std::array<T, 1> {
		{T(1)}
	},
	std::array<T, 1> {
		{T(0)}
	}
>;

template<typename T>
using RK_2_2 = ButcherTable<T, 2, T(1),
	std::array<std::array<T, 2>, 2> {
		{
			std::array<T, 2> { {T(0), T(0)}},
			std::array<T, 2> { {T(1), T(0)}}
		}
	}, std::array<T, 2> {
		{T(0.5), T(0.5)}
	},
	std::array<T, 2> {
		{T(0), T(1)}
	}
>;

template<typename T>
using RK_3_3 = ButcherTable<T, 3, T(1),
	std::array<std::array<T, 3>, 3> {
		{
			std::array<T, 3> { {T(0), T(0), T(0)}},
			std::array<T, 3> { {T(1), T(0), T(0)}},
			std::array<T, 3> { {T(0.25), T(0.25), T(0)}}
		}
	}, std::array<T, 3> {
		{T(1./6.), T(1./6.), T(2./3.)}
	},
	std::array<T, 3> {
		{T(0), T(1), T(0.5)}
	}
>;

template<typename T>
using RK_10_4 = ButcherTable<T, 10, T(6.0),
		std::array<std::array<T, 10>, 10> {
			{
				std::array<T, 10> { {T(0),       T(0),       T(0),       T(0),       T(0),       T(0),      T(0),      T(0),      T(0),      T(0)}},
				std::array<T, 10> { {T(1)/T(6),  T(0),       T(0),       T(0),       T(0),       T(0),      T(0),      T(0),      T(0),      T(0)}},
				std::array<T, 10> { {T(1)/T(6),  T(1)/T(6),  T(0),       T(0),       T(0),       T(0),      T(0),      T(0),      T(0),      T(0)}},
				std::array<T, 10> { {T(1)/T(6),  T(1)/T(6),  T(1)/T(6),  T(0),       T(0),       T(0),      T(0),      T(0),      T(0),      T(0)}},
				std::array<T, 10> { {T(1)/T(6),  T(1)/T(6),  T(1)/T(6),  T(1)/T(6),  T(0),       T(0),      T(0),      T(0),      T(0),      T(0)}},
				std::array<T, 10> { {T(1)/T(15), T(1)/T(15), T(1)/T(15), T(1)/T(15), T(1)/T(15), T(0),      T(0),      T(0),      T(0),      T(0)}},
				std::array<T, 10> { {T(1)/T(15), T(1)/T(15), T(1)/T(15), T(1)/T(15), T(1)/T(15), T(1)/T(6), T(0),      T(0),      T(0),      T(0)}},
				std::array<T, 10> { {T(1)/T(15), T(1)/T(15), T(1)/T(15), T(1)/T(15), T(1)/T(15), T(1)/T(6), T(1)/T(6), T(0),      T(0),      T(0)}},
				std::array<T, 10> { {T(1)/T(15), T(1)/T(15), T(1)/T(15), T(1)/T(15), T(1)/T(15), T(1)/T(6), T(1)/T(6), T(1)/T(6), T(0),      T(0)}},
				std::array<T, 10> { {T(1)/T(15), T(1)/T(15), T(1)/T(15), T(1)/T(15), T(1)/T(15), T(1)/T(6), T(1)/T(6), T(1)/T(6), T(1)/T(6), T(0)}}}
			},
		std::array<T, 10> {
			{T(0.1), T(0.1), T(0.1), T(0.1), T(0.1), T(0.1), T(0.1), T(0.1), T(0.1), T(0.1)}
		},
		std::array<T, 10> {
			{T(0), T(1)/T(6), T(1)/T(3), T(1)/T(2), T(2)/T(3), T(1)/T(3), T(1)/T(2), T(2)/T(3), T(5)/T(6), T(1)}
		}
>;

template<typename T, int N>
struct RungeKutta;

template<typename T>
struct RungeKutta<T, 1> {
	using type = RK_1_1<T>;
};

template<typename T>
struct RungeKutta<T, 2> {
	using type = RK_2_2<T>;
};

template<typename T>
struct RungeKutta<T, 3> {
	using type = RK_3_3<T>;
};

template<typename T>
struct RungeKutta<T, 4> {
	using type = RK_10_4<T>;
};

