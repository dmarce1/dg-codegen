#pragma once

#include "Matrix.hpp"
#include "Polynomial.hpp"
#include "Vector.hpp"

template<typename T>
struct RK_1_1 {
	constexpr RK_1_1() :
			a_( { { T(0) } }), b_( { T(1) }), c_( { T(0) }) {
	}
	T a(int i, int j) const {
		return a_[i, j];
	}
	T b(int k) const {
		return b_[k];
	}
	T c(int k) const {
		return c_[k];
	}
	static constexpr int stageCount() {
		return 1;
	}
	static constexpr T cfl() {
		return T(1);
	}
private:
	Math::SquareMatrix<T, 1> a_;
	Math::Vector<T, 1> b_;
	Math::Vector<T, 1> c_;
};

template<typename T>
struct RK_2_2 {
	constexpr RK_2_2() :
			a_( { { T(0), T(0) }, { T(1), T(0) } }), b_( { T(0.5), T(0.5) }), c_( { T(0), T(1) }) {
	}
	T a(int i, int j) const {
		return a_[i, j];
	}
	T b(int k) const {
		return b_[k];
	}
	T c(int k) const {
		return c_[k];
	}
	static constexpr int stageCount() {
		return 2;
	}
	static constexpr T cfl() {
		return T(1);
	}
private:
	Math::SquareMatrix<T, 2> a_;
	Math::Vector<T, 2> b_;
	Math::Vector<T, 2> c_;
};

template<typename T>
struct RK_3_3 {
	constexpr RK_3_3() :
			a_( { { T(0), T(0), T(0) }, { T(1), T(0), T(0) }, { T(0.25), T(0.25), T(0) } }), b_( { T(1) / T(6), T(1) / T(6), T(2) / T(3) }), c_( { T(0), T(1),
					T(0.5) }) {
	}
	T a(int i, int j) const {
		return a_[i, j];
	}
	T b(int k) const {
		return b_[k];
	}
	T c(int k) const {
		return c_[k];
	}
	static constexpr int stageCount() {
		return 3;
	}
	static constexpr T cfl() {
		return T(1);
	}
private:
	Math::SquareMatrix<T, 3> a_;
	Math::Vector<T, 3> b_;
	Math::Vector<T, 3> c_;
};

template<typename T>
struct RK_10_4 {
	static constexpr T c0 = T(0);
	static constexpr T c1 = T(1) / T(15);
	static constexpr T c2 = T(1) / T(6);
	constexpr RK_10_4() :
			a_( { { c0, c0, c0, c0, c0, c0, c0, c0, c0, c0 }, { c2, c0, c0, c0, c0, c0, c0, c0, c0, c0 }, { c2, c2, c0, c0, c0, c0, c0, c0, c0, c0 }, { c2, c2,
					c2, c0, c0, c0, c0, c0, c0, c0 }, { c2, c2, c2, c2, c0, c0, c0, c0, c0, c0 }, { c1, c1, c1, c1, c1, c0, c0, c0, c0, c0 }, { c1, c1, c1, c1,
					c1, c2, c0, c0, c0, c0 }, { c1, c1, c1, c1, c1, c2, c2, c0, c0, c0 }, { c1, c1, c1, c1, c1, c2, c2, c2, c0, c0 }, { c1, c1, c1, c1, c1, c2,
					c2, c2, c2, c0 } }), b_( { T(.1), T(.1), T(.1), T(.1), T(.1), T(.1), T(.1), T(.1), T(.1), T(.1) }), c_(
					{ T(0), T(1) / T(6), T(1) / T(3), T(1) / T(2), T(2) / T(3), T(1) / T(3), T(1) / T(2), T(2) / T(3), T(5) / T(6), T(1) }) {
	}
	T a(int i, int j) const {
		return a_[i, j];
	}
	T b(int k) const {
		return b_[k];
	}
	T c(int k) const {
		return c_[k];
	}
	static constexpr int stageCount() {
		return 10;
	}
	static constexpr T cfl() {
		return T(6);
	}
private:
	Math::SquareMatrix<T, 10> a_;
	Math::Vector<T, 10> b_;
	Math::Vector<T, 10> c_;
};

template<typename T, int>
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

template<typename T, typename RK>
Math::Polynomial<T> rungeKuttaStabilityPolynomial() {
	static const RK rk { };
	static constexpr int S = rk.stageCount();
	Math::Polynomial<T> R;
	Math::Vector<Math::Polynomial<T>, S> w;
	Math::Polynomial<T> z;
	z[1] = T(1);
	for (int i = 0; i < S; i++) {
		for (int j = 0; j < i; j++) {
			w[i] += rk.a(i, j) * w[j];
		}
		w[i] *= z;
		w[i] += T(1);
	}
	for (int k = 0; k < S; k++) {
		R += rk.b(k) * w[k];
	}
	R *= z;
	R += T(1);
	return R;
}

template<typename T, typename RK>
T rungeKuttaStabilityRadius() {
	static const RK rk { };
	auto const R = rungeKuttaStabilityPolynomial<T, RK>();
	auto g = [R](T const &x) {
		return std::abs(R(-x)) - T(1);
	};
	T x1 = T(0);
	T x2 = T(2);
	while (g(x2) <= T(0)) {
		x2 *= T(2);
	}
	for (int iter = 0; iter < 64; ++iter) {
		T xm = (x1 + x2) * T(0.5);
		if (g(xm) <= T(0)) {
			x1 = xm;
		} else {
			x2 = xm;
		}
	}
	return x1;
}
