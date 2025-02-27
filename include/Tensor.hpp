/*
 * Tensor.hpp
 *
 *  Created on: Feb 17, 2025
 *      Author: dmarce1
 */

#ifndef INCLUDE_TENSOR_HPP_
#define INCLUDE_TENSOR_HPP_

#include "Vector.hpp"

using namespace Math;

template<typename Type, int Rank, int Ndim>
struct Tensor;


template<typename Type, int Rank, int Ndim>
struct Tensor {
	constexpr Tensor() = default;
	template<typename ...Args>
	auto& operator[](int i, Args ...args) {
		return values[i].operator[](args...);
	}
	template<typename ...Args>
	auto const& operator[](int i, Args ...args) const {
		return values[i].operator[](args...);
	}
	template<int Mdim>
	auto& operator[](Vector<int, Mdim> const &indices) {
		if constexpr (Mdim > 1) {
			return values[indices[0]][*reinterpret_cast<Vector<int, Mdim - 1> const*>(&indices[1])];
		} else {
			return values[indices[0]];
		}
	}
	template<int Mdim>
	auto const& operator[](Vector<int, Mdim> const &indices) const {
		if constexpr (Mdim > 1) {
			return values[indices[0]][*reinterpret_cast<Vector<int, Mdim - 1> const*>(&indices[1])];
		} else {
			return values[indices[0]];
		}
	}
	auto& operator[](int i) {
		return values[i];
	}
	auto const& operator[](int i) const {
		return values[i];
	}
	Tensor(Tensor const &A) :
			values(A.values) {
	}
	Tensor(Type a) {
		*this = a;
	}
	Tensor& operator=(Tensor const &A) {
		values = A.values;
		return *this;
	}
	Tensor& operator=(Type a) {
		for (int i = 0; i < Ndim; i++) {
			values[i] = a;
		}
		return *this;
	}
	Tensor& operator*=(Type a) {
		values *= a;
		return *this;
	}
	Tensor& operator+=(Tensor const &A) {
		values += A.values;
		return *this;
	}
	Tensor operator-=(Tensor const &A) {
		values -= A.values;
		return *this;
	}
	Tensor operator*(Type a) {
		Tensor B;
		B.values = values * a;
		return B;
	}
	friend Tensor operator*(Tensor<Type, 0, Ndim> a, Tensor const &B) {
		Tensor C;
		for (int d = 0; d < Ndim; d++) {
			C.values[d] = a * B.values[d];
		}
		return C;
	}
	Tensor operator+(Tensor const &A) const {
		Tensor B;
		B.values = values + A.values;
		return B;
	}
	Tensor operator-(Tensor const &A) const {
		Tensor B;
		for (int d = 0; d < Ndim; d++) {
			B.values[d] = values[d] - A.values[d];
		}
		return B;
	}
	Tensor operator+() const {
		return *this;
	}
	Tensor operator-() const {
		Tensor A;
		for (int d = 0; d < Ndim; d++) {
			A.values[d] = -values[d];
		}
		return A;
	}
	template<int Rank1>
	auto operator*(Tensor<Type, Rank1, Ndim> const &other) const {
		static constexpr int Rank2 = Rank + Rank1 - 2;
		Tensor<Type, Rank2, Ndim> C;
		if constexpr (Rank > 1) {
			for (int i = 0; i < Ndim; i++) {
				C[i] = values[i] * other;
			}
		} else {
			C = values[0] * other[0];
			for (int i = 1; i < Ndim; i++) {
				C += values[i] * other[i];
			}
		}
		return C;
	}
	template<int Rank1>
	auto operator^(Tensor<Type, Rank1, Ndim> const &A) const {
		static constexpr int Rank2 = Rank + Rank1;
		Tensor<Type, Rank2, Ndim> B;
		if constexpr (Rank > 1) {
			for (int i = 0; i < Ndim; i++) {
				B[i] = values[i] ^ A;
			}
		} else {
			for (int i = 0; i < Ndim; i++) {
				for (int j = 0; j < Ndim; j++) {
					B[i, j] = values[i] * A[j];
				}
			}
		}
		return B;
	}
	template<int K>
	auto subtensor(int k) const {
		static constexpr int Rank1 = Rank - 1;
		Tensor<Type, Rank1, Ndim> B;
		if constexpr (K > 0) {
			for (int i = 0; i < Ndim; i++) {
				B[i] = values[i].template subtensor<K - 1>(k);
			}
		} else {
			B = values[k];
		}
		return B;
	}

	template<int I0 = 0, int I1 = 1>
	Tensor<Type, Rank - 2, Ndim> trace() const {
		static constexpr int Rank1 = Rank - 2;
		Tensor<Type, Rank1, Ndim> B;
		if constexpr (Rank1 == 0) {
			B = values[0][0];
			for (int k = 1; k < Ndim; k++) {
				B += values[k][k];
			}
		} else if constexpr (I0 > I1) {
			B = trace<I1, I0>();
		} else if constexpr (I0 > 0) {
			for (int k = 0; k < Ndim; k++) {
				B[k] = values[k].template trace<I0 - 1, I1 - 1>();
			}
		} else {
			B = values[0].template subtensor<I1 - 1>(0);
			for (int k = 1; k < Ndim; k++) {
				B += values[k].template subtensor<I1 - 1>(k);
			}
		}
		return B;
	}

	template<int I0 = 0, int I1 = 1>
	Tensor transpose() const {
		Tensor A;
		if constexpr (Rank == 2) {
			for (int i = 0; i < Ndim; i++) {
				for (int k = 0; k < Ndim; k++) {
					A[i, k] = values[k][i];
				}
			}
		} else if constexpr (I0 > I1) {
			A = transpose<I1, I0>();
		} else if constexpr (I0 > 0) {
			for (int i = 0; i < Ndim; i++) {
				A[i] = values[i].template transpose<I0 - 1, I1 - 1>();
			}
		} else {
			for (int i = 0; i < Ndim; i++) {
				for (int k = 0; k < Ndim; k++) {
					A[i, k] = values[k].template subtensor<I1 - 1>(i);
				}
			}
		}
		return A;
	}
	template<int I0 = 0, int I1 = 1>
	auto symmetric() const {
		return (*this + transpose<I0, I1>());
	}
	template<int I0 = 0, int I1 = 1>
	auto antisymmetric() const {
		return (*this - transpose<I0, I1>());
	}

private:
	Vector<Tensor<Type, Rank - 1, Ndim>, Ndim> values;
};

template<typename Type, int Ndim>
struct Tensor<Type, 0, Ndim> {
	Tensor() = default;
	template<int Mdim>
	constexpr Tensor(Tensor<Type, 0, Mdim> const &v) :
			value(Type(v)) {
	}
	constexpr Tensor(Type const &v) :
			value(v) {
	}
	Type& operator[]() {
		return value;
	}
	auto const& operator[]() const {
		return value;
	}
	Tensor& operator=(Type a) {
		value = a;
		return *this;
	}
	Tensor& operator+=(Tensor a) {
		value += a.value;
		return *this;
	}
	Tensor& operator-=(Tensor a) {
		value -= a.value;
		return *this;
	}
	Tensor& operator*=(Tensor a) {
		value *= a.value;
		return *this;
	}
	Tensor& operator/=(Tensor a) {
		value /= a.value;
		return *this;
	}
	Tensor operator+() const {
		return Tensor(+value);
	}
	Tensor operator-() const {
		return Tensor(-value);
	}
	Tensor operator+(Tensor a) const {
		return Tensor(value + a.value);
	}
	Tensor operator-(Tensor a) const {
		return Tensor(value - a.value);
	}
	Tensor operator*(Tensor a) const {
		return Tensor(value * a.value);
	}
	Tensor operator/(Tensor a) const {
		return Tensor(value / a.value);
	}
	operator Type() const {
		return value;
	}
	friend Tensor operator+(Type a, Tensor b) {
		return Tensor(a + b.value);
	}
	friend Tensor operator-(Type a, Tensor b) {
		return Tensor(a - b.value);
	}
	friend Tensor operator*(Type a, Tensor b) {
		return Tensor(a * b.value);
	}
	friend Tensor operator/(Type a, Tensor b) {
		return Tensor(a / b.value);
	}
private:
	Type value;
};

template<typename Type>
struct KroneckerDelta {
	Type operator[](int i, int j) const {
		return Type(i == j);
	}
};

template<typename Type, int Ndim>
auto tensorCotensor(Tensor<Type, 2, Ndim> const &A, int I, int K) {
	Tensor<Type, 2, Ndim - 1> B;
	for (int i = 0; i < Ndim - 1; i++) {
		for (int k = 0; k < Ndim - 1; k++) {
			B[i, k] = A[i + int(i >= I), k + int(k >= K)];
		}
	}
	return B;
}

template<typename Type, int Ndim>
Type tensorDeterminant(Tensor<Type, 2, Ndim> const &A) {
	static constexpr Type zero = Type(0);
	Type det;
	if constexpr (Ndim == 1) {
		return Type(A[0, 0]);
	} else {
		det = zero;
		for (int i = 0; i < Ndim - 1; i++) {
			det += Type(negativeOne2Power<Type>(i) * A[0, i] * tensorDeterminant(tensorCotensor(A, 0, i)));
		}
	}
	return det;
}

template<typename Type, int Ndim>
auto tensorInverse(Tensor<Type, 2, Ndim> const &A) {
	static constexpr Type zero = Type(0), one = Type(1);
	Tensor<Type, 2, Ndim> Ainv;
	Vector<Type, Ndim> cofactor;
	Type det = zero;
	for (int i = 0; i < Ndim; i++) {
		for (int k = 0; k < Ndim; k++) {
			Ainv[i, k] = negativeOne2Power<Type>(i + k) * tensorDeterminant(tensorCotensor(A, k, i));
		}
	}
	for (int i = 0; i < Ndim; i++) {
		det += Type(A[0, i] * Ainv[i, 0]);
	}
	Ainv = (one / det) * Ainv;
	return Ainv;
}

#endif /* INCLUDE_TENSOR_HPP_ */
