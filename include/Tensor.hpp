/*
 * Tensor.hpp
 *
 *  Created on: Mar 2, 2025
 *      Author: dmarce1
 */

#ifndef INCLUDE_TENSOR_HPP_
#define INCLUDE_TENSOR_HPP_

#include <array>
#include <functional>

namespace Tensors {

template<char C>
struct Index {
};

template<int I1, int I2>
constexpr int Sym = ((1 << I1) | (1 << I2));

template<typename T, int D, int R, int ...S>
struct Tensor;

template<typename T, int D>
struct Tensor<T, D, 0> {
	using value_type = T;
	auto operator[]() const {
		return a;
	}
	auto operator()() const;
private:
	T a;
};

template<typename T, int D>
struct Tensor<T, D, 1> : public std::array<T, D> {
	using value_type = T;
	auto operator()(int i) const {
		return (*this)[i];
	}
	template<char C>
	auto operator()(Index<C>) const;
};

template<typename T, int D, int ...S>
struct Tensor<T, D, 2, S...> {
	using value_type = T;
	static constexpr bool isSymmetric = sizeof...(S);
	auto operator[](int i, int j) const {
		if constexpr (isSymmetric) {
			if (i < j) {
				std::swap(i, j);
			}
			return U[((i * (i + 1)) >> 1) + j];
		} else {
			return U[D * i + j];
		}
	}
	auto& operator[](int i, int j) {
		if constexpr (isSymmetric) {
			if (i < j) {
				std::swap(i, j);
			}
			return U[((i * (i + 1)) >> 1) + j];
		} else {
			return U[D * i + j];
		}
	}
	template<char C1, char C2>
	auto operator()(Index<C1>, Index<C2>) const;
	template<char C>
	auto operator()(int, Index<C>) const;
	template<char C>
	auto operator()(Index<C>, int) const;
	template<bool TR>
	struct access1d {
		access1d(Tensor const &t, int i) :
				tensor(t), index(i) {
		}
		auto operator[](int i) const {
			if constexpr (TR) {
				return tensor[i, index];
			} else {
				return tensor[index, i];
			}
		}
	private:
		Tensor const &tensor;
		int index;
	};
	Tensor<T, D - 1, 2> subTensor(int j, int k) const {
		Tensor<T, D - 1, 2> sub;
		for (int n = 0; n < j; n++) {
			for (int m = 0; m < k; m++) {
				sub[n, m] = (*this)[n, m];
			}
			for (int m = k + 1; m < D; m++) {
				sub[n, m - 1] = (*this)[n, m];
			}
		}
		for (int n = j + 1; n < D; n++) {
			for (int m = 0; m < k; m++) {
				sub[n - 1, m] = (*this)[n, m];
			}
			for (int m = k + 1; m < D; m++) {
				sub[n - 1, m - 1] = (*this)[n, m];
			}
		}
		return sub;
	}
	T determinant() const {
		if constexpr (D == 1) {
			return U[0];
		} else {
			T sum = T(0);
			T sgn = T(1);
			for (int i = 0; i < D; i++) {
				sum += sgn * subTensor(0, i).determinant();
				sgn = -sgn;
			}
			return sum;
		}
	}
	Tensor inverse() const {
		Tensor I;
		T sgn_n = T(1);
		T det = T(0);
		for (int n = 0; n < D; n++) {
			T sgn = sgn_n;
			for (int m = 0; m <= n; m++) {
				I[n, m] = sgn * subTensor(n, m).determinant();
				sgn = -sgn;
			}
			sgn_n = -sgn_n;
		}
		for (int n = 0; n < D; n++) {
			det += I[0, n];
		}
		T const dinv = T(1) / det;
		for (int n = 0; n < D; n++) {
			for (int m = 0; m < D; m++) {
				I[n, m] *= dinv;
			}
		}
		return I;
	}
private:
	static constexpr int Size = isSymmetric ? (D * (D + 1) / 2) : (D * D);
	std::array<T, Size> U;
};

template<typename T, int D>
struct Unit {
	using value_type = T;
	Unit(int j_) :
			j(j_) {
	}
	auto operator[](int i) const {
		return T(i == j);
	}
	template<char C>
	auto operator()(Index<C>) const;
private:
	int const j;
};

template<typename T, int D>
struct Delta {
	using value_type = T;
	auto operator[](int i, int j) const {
		return T(i == j);
	}
	template<char C1, char C2>
	auto operator()(Index<C1>, Index<C2>) const;
	template<char C>
	auto operator()(int i, Index<C>) const;
	template<char C>
	auto operator()(Index<C>, int i) const;
};

template<typename T, int D>
struct Tensor<T, D, 3> {
	using value_type = T;
	auto operator[](int i, int j, int k) const {
		return U[D * (D * i + j) + k];
	}
	auto& operator[](int i, int j, int k) {
		return U[D * (D * i + j) + k];
	}
	template<char C1, char C2, char C3>
	auto operator()(Index<C1>, Index<C2>, Index<C3>) const;
private:
	std::array<T, D * D * D> U;
};

template<typename T, int D>
struct Tensor<T, D, 3, Sym<0, 1>> {
	using value_type = T;
	auto operator[](int i, int j, int k) const {
		return U[i * D2 + ((j * (j + 1)) >> 1) + k];
	}
	auto& operator[](int i, int j, int k) {
		return U[i * D2 + ((j * (j + 1)) >> 1) + k];
	}
	template<char C1, char C2, char C3>
	auto operator()(Index<C1>, Index<C2>, Index<C3>) const;
	template<char C2, char C3>
	auto operator()(int i, Index<C2>, Index<C3>);
private:
	static constexpr int D2 = D * (D + 1) / 2;
	std::array<T, D * D2> U;
};

template<typename A>
class Tensor0Expr {
	A ptr;
public:
	using value_type = decltype(ptr.operator[]());
	Tensor0Expr(A a) :
			ptr(a) {
	}
	auto operator[]() const {
		return ptr.operator[]();
	}
	auto operator=(Tensor0Expr const &other) {
		ptr.operator[]() = other.operator[]();
		return *this;
	}
	operator value_type() const {
		return ptr.operator[]();
	}
};

template<typename A, int D, char I>
class Tensor1Expr {
	A ptr;
public:
	using value_type = decltype(ptr[0]);
	Tensor1Expr(A a) :
			ptr(a) {
	}
	auto operator[](int i) const {
		return ptr[i];
	}
	template<typename U>
	auto operator=(U const &other) {
		if constexpr (std::is_convertible<value_type, U>::value) {
			for (int i = 0; i < D; i++) {
				ptr[i] = other;
			}
		} else {
			for (int i = 0; i < D; i++) {
				ptr[i] = other[i];
			}
		}
		return *this;
	}
	template<typename U>
	auto operator+=(U const &other) {
		if constexpr (std::is_convertible<value_type, U>::value) {
			for (int i = 0; i < D; i++) {
				ptr[i] += other;
			}
		} else {
			for (int i = 0; i < D; i++) {
				ptr[i] += other[i];
			}
		}
		return *this;
	}
};

template<typename T>
struct isSymmetric2 {
	static constexpr bool value = false;
};

template<typename T, int D>
struct isSymmetric2<Tensor<T, D, 2, Sym<0, 1>>> {
	static constexpr bool value = true;
};

template<typename A, int D, char I, char J>
class Tensor2Expr {
	A ptr;
public:
	using value_type = decltype(ptr[0, 0]);
	Tensor2Expr(A a) :
			ptr(a) {
	}
	constexpr auto operator[](int i, int j) const {
		return ptr[i, j];
	}
	template<template<typename, int, char, char> class Expr, typename B>
	auto operator=(Expr<B, D, I, J> const &other) {
		if constexpr (isSymmetric2<A>::value) {
			for (int i = 0; i < D; i++) {
				for (int j = 0; j <= i; j++) {
					ptr[i, j] = other[i, j];
				}
			}
		} else {
			for (int i = 0; i < D; i++) {
				for (int j = 0; j < D; j++) {
					ptr[i, j] = other[i, j];
				}
			}
		}
		return *this;
	}
	template<typename B>
	auto operator=(B const &other) {
		if constexpr (isSymmetric2<A>::value) {
			for (int i = 0; i < D; i++) {
				for (int j = 0; j <= i; j++) {
					ptr[i, j] = other;
				}
			}
		} else {
			for (int i = 0; i < D; i++) {
				for (int j = 0; j < D; j++) {
					ptr[i, j] = other;
				}
			}
		}
		return *this;
	}
	template<template<typename, int, char, char> class Expr, typename B>
	auto operator+=(Expr<B, D, I, J> const &other) {
		if constexpr (isSymmetric2<A>::value) {
			for (int i = 0; i < D; i++) {
				for (int j = 0; j <= i; j++) {
					ptr[i, j] += other[i, j];
				}
			}
		} else {
			for (int i = 0; i < D; i++) {
				for (int j = 0; j < D; j++) {
					ptr[i, j] += other[i, j];
				}
			}
		}
		return *this;
	}
	template<template<typename, int, char, char> class Expr, typename B>
	auto operator-=(Expr<B, D, I, J> const &other) {
		if constexpr (isSymmetric2<A>::value) {
			for (int i = 0; i < D; i++) {
				for (int j = 0; j <= i; j++) {
					ptr[i, j] -= other[i, j];
				}
			}
		} else {
			for (int i = 0; i < D; i++) {
				for (int j = 0; j < D; j++) {
					ptr[i, j] -= other[i, j];
				}
			}
		}
		return *this;
	}
	template<template<typename, int, char, char> class Expr, typename B>
	auto operator-=(Expr<B, D, J, I> const &other) {
		if constexpr (isSymmetric2<A>::value) {
			for (int i = 0; i < D; i++) {
				for (int j = 0; j <= i; j++) {
					ptr[i, j] -= other[j, i];
				}
			}
		} else {
			for (int i = 0; i < D; i++) {
				for (int j = 0; j < D; j++) {
					ptr[i, j] -= other[j, i];
				}
			}
		}
		return *this;
	}
	template<template<typename, int, char, char> class Expr, typename B>
	auto operator+=(Expr<B, D, J, I> const &other) {
		if constexpr (isSymmetric2<A>::value) {
			for (int i = 0; i < D; i++) {
				for (int j = 0; j <= i; j++) {
					ptr[i, j] += other[j, i];
				}
			}
		} else {
			for (int i = 0; i < D; i++) {
				for (int j = 0; j < D; j++) {
					ptr[i, j] += other[j, i];
				}
			}
		}
		return *this;
	}
};

template<typename A, int D, char I, char J, char K>
class Tensor3Expr {
	A ptr;
public:
	using value_type = decltype(ptr[0,0,0]);
	Tensor3Expr(A a) :
			ptr(a) {
	}
	auto operator[](int i, int j, int k) const {
		return ptr[i, j, k];
	}
	template<typename B>
	auto operator=(B const &other) {
		for (int i = 0; i < D; i++) {
			for (int j = 0; j < D; j++) {
				for (int k = 0; k < D; k++) {
					ptr[i, j, k] = other;
				}
			}
		}
		return *this;
	}
	template<template<typename, int, char, char, char> class Expr, typename B>
	auto operator=(Expr<B, D, I, J, K> const &other) {
		for (int i = 0; i < D; i++) {
			for (int j = 0; j < D; j++) {
				for (int k = 0; k < D; k++) {
					ptr[i, j, k] = other[i, j, k];
				}
			}
		}
		return *this;
	}
	template<template<typename, int, char, char, char> class Expr, typename B>
	auto operator=(Expr<B, D, J, K, I> const &other) {
		for (int i = 0; i < D; i++) {
			for (int j = 0; j < D; j++) {
				for (int k = 0; k < D; k++) {
					ptr[i, j, k] = other[j, k, i];
				}
			}
		}
		return *this;
	}
	template<template<typename, int, char, char, char> class Expr, typename B>
	auto operator=(Expr<B, D, K, J, I> const &other) {
		for (int i = 0; i < D; i++) {
			for (int j = 0; j < D; j++) {
				for (int k = 0; k < D; k++) {
					ptr[i, j, k] = other[k, j, i];
				}
			}
		}
		return *this;
	}
};

template<typename T, int D>
auto Tensor<T, D, 0>::operator()() const {
	return Tensor0Expr<Tensor<T, D, 0>>(&a);
}

template<typename T, int D>
template<char C>
auto Tensor<T, D, 1>::operator()(Index<C>) const {
	return Tensor1Expr<Tensor<T, D, 1>, D, C>(*this);
}

template<typename T, int D, int ...S>
template<char C1, char C2>
auto Tensor<T, D, 2, S...>::operator()(Index<C1>, Index<C2>) const {
	return Tensor2Expr<Tensor<T, D, 2, S...>, D, C1, C2>(*this);
}

template<typename T, int D>
template<char C1, char C2>
auto Delta<T, D>::operator()(Index<C1>, Index<C2>) const {
	return Tensor2Expr<Delta<T, D>, D, C1, C2>(*this);
}

template<typename T, int D>
template<char C1, char C2, char C3>
auto Tensor<T, D, 3>::operator()(Index<C1>, Index<C2>, Index<C3>) const {
	return Tensor3Expr<Tensor<T, D, 3>, D, C1, C2, C3>(*this);
}

template<typename T, int D>
template<char C1, char C2, char C3>
auto Tensor<T, D, 3, Sym<0, 1>>::operator()(Index<C1>, Index<C2>, Index<C3>) const {
	return Tensor3Expr<Tensor<T, D, 3, Sym<0, 1>>, D, C1, C2, C3>(*this);
}

template<typename T, int D>
template<char C2, char C3>
auto Tensor<T, D, 3, Sym<0, 1>>::operator()(int i, Index<C2>, Index<C3>) {
	return Tensor2Expr<value_type*, D, C2, C3>(U.data() + i * D2);
}

template<typename T, int D, int ...S>
template<char C>
auto Tensor<T, D, 2, S...>::operator()(int i, Index<C>) const {
	Tensor<T, D, 2, S...>::access1d<false> a { *this, i };
	return Tensor1Expr<Tensor<T, D, 2, S...>::access1d<false>, D, C>(a);
}

template<typename T, int D, int ...S>
template<char C>
auto Tensor<T, D, 2, S...>::operator()(Index<C>, int i) const {
	Tensor<T, D, 2, S...>::access1d<true> a { *this, i };
	return Tensor1Expr<Tensor<T, D, 2, S...>::access1d<true>, D, C>(a);
}

template<typename T, int D>
template<char C>
auto Delta<T, D>::operator()(int i, Index<C>) const {
	Unit<T, D> unit(i);
	return Tensor1Expr<Unit<T, D>, D, C>(unit);
}

template<typename T, int D>
template<char C>
auto Delta<T, D>::operator()(Index<C>, int i) const {
	Unit<T, D> unit(i);
	return Tensor1Expr<Unit<T, D>, D, C>(unit);
}

template<typename T, typename U>
using sum_type = decltype(T() + U());

template<typename A>
struct TensorNegate {
	TensorNegate(A const &a) :
			ptrA(a) {
	}
	auto operator[](auto ... i) const {
		return -ptrA.operator[](i...);
	}
private:
	A ptrA;
};

template<typename A, typename B, int D, int T = 0>
struct Tensor3Add {

	Tensor3Add(A const &a, B const &b) :
			ptrA(a), ptrB(b) {
	}

	auto operator[](int i, int j, int k) const {
		if constexpr (T == 0x01) {
			return ptrA[i, j, k] + ptrB[i, k, j];
		} else if constexpr (T == 0x02) {
			return ptrA[i, j, k] + ptrB[k, j, i];
		} else if constexpr (T == 0x12) {
			return ptrA[i, j, k] + ptrB[j, i, k];
		} else {
			return ptrA[i, j, k] + ptrB[i, j, k];
		}
	}

private:
	A ptrA;
	B ptrB;
};

template<typename A, typename B, int D, typename O, int I0 = -1, int I1 = -1>
struct TensorBinaryOp {
	TensorBinaryOp(A const &a, B const &b, O const &op_) :
			ptrA(a), ptrB(b), op(op_) {
	}
	auto operator[](auto ...i) const {
		static constexpr int Rank = sizeof...(i);
		std::array<int, Rank> I;
		int j = 0;
		((I[j++] = i),...);
		if (I0 != -1) {
			std::swap(I[I0], I[I1]);
		}
		return op(ptrA.operator[](i...), std::apply([this](auto ...i) {
			return ptrB.operator[](i...);
		}, I));
	}
private:
	A ptrA;
	B ptrB;
	O op;
};

template<typename A, typename B>
struct Tensor0xN {

	Tensor0xN(A const &a_, B const &b) :
			a(a_), ptrB(b) {
	}

	auto operator[](auto ...j) const {
		return a * ptrB.operator[](j...);
	}

private:
	A a;
	B ptrB;
};

template<typename A, typename B>
struct Tensor1xN {

	Tensor1xN(A const &a, B const &b) :
			ptrA(a), ptrB(b) {
	}

	auto operator[](int i, auto ...j) const {
		return ptrA[i] * ptrB.operator[](j...);
	}

private:
	A ptrA;
	B ptrB;
};

template<typename A, typename B>
struct TensorNx1 {

	TensorNx1(A const &a, B const &b) :
			ptrA(a), ptrB(b) {
	}

	auto operator[](auto ...j, int i) const {
		return ptrA.operator[](j...) * ptrB[i];
	}

private:
	A ptrA;
	B ptrB;
};

template<typename A, typename B, int D>
struct Tensor1dot1 {

	Tensor1dot1(A const &a, B const &b) :
			ptrA(a), ptrB(b) {
	}

	auto operator[]() const {
		auto sum = ptrA[0] * ptrB[0];
		for (int d = 1; d < D; d++) {
			sum += ptrA[d] * ptrB[d];
		}
		return sum;
	}

private:
	A ptrA;
	B ptrB;
};

template<typename A, typename B, int D>
struct Tensor1dotN {

	Tensor1dotN(A const &a, B const &b) :
			ptrA(a), ptrB(b) {
	}

	auto operator[](auto ...i) const {
		auto sum = ptrA[0] * ptrB.operator[](0, i...);
		for (int d = 1; d < D; d++) {
			sum += ptrA[d] * ptrB.operator[](d, i...);
		}
		return sum;
	}

private:
	A ptrA;
	B ptrB;
};

template<typename A, typename B, int D>
struct Tensor2dotN {

	Tensor2dotN(A const &a, B const &b) :
			ptrA(a), ptrB(b) {
	}

	auto operator[](int i, auto ...j) const {
		auto sum = ptrA[i, 0] * ptrB.operator[](0, j...);
		for (int d = 1; d < D; d++) {
			sum += ptrA[i, d] * ptrB.operator[](d, j...);
		}
		return sum;
	}

private:
	A ptrA;
	B ptrB;
};

template<typename A, typename B, int D>
struct TensorNdot1 {

	TensorNdot1(A const &a, B const &b) :
			ptrA(a), ptrB(b) {
	}

	auto operator[](auto ...i) const {
		auto sum = ptrA.operator[](i..., 0) * ptrB[0];
		for (int d = 1; d < D; d++) {
			sum += ptrA.operator[](i..., d) * ptrB[d];
		}
		return sum;
	}

private:
	A ptrA;
	B ptrB;
};

template<typename A, typename B, int D>
struct TensorNdotdot2 {

	TensorNdotdot2(A const &a, B const &b) :
			ptrA(a), ptrB(b) {
	}
	auto operator[](auto ...k) const {
		auto sum = ptrA.operator[](k..., 0, 0) * ptrB[0, 0];
		for (int j = 1; j < D; j++) {
			sum += ptrA.operator[](k..., 0, j) * ptrB[0, j];
			for (int i = 0; i < D; i++) {
				sum += ptrA.operator[](k..., j, i) * ptrB[j, i];
			}
		}
		return sum;
	}

private:
	A ptrA;
	B ptrB;
};

template<typename A, typename B, int D>
struct Tensor2dotdot3 {

	Tensor2dotdot3(A const &a, B const &b) :
			ptrA(a), ptrB(b) {
	}
	auto operator[](auto k) const {
		auto sum = ptrA.operator[](0, 0) * ptrB[0, 0, k];
		for (int j = 1; j < D; j++) {
			sum += ptrA.operator[](0, j) * ptrB[j, 0, k];
			for (int i = 0; i < D; i++) {
				sum += ptrA.operator[](j, i) * ptrB[i, j, k];
			}
		}
		return sum;
	}

private:
	A ptrA;
	B ptrB;
};

template<typename A, typename B, int D>
struct Tensor3dotdotdot3 {

	Tensor3dotdotdot3(A const &a, B const &b) :
			ptrA(a), ptrB(b) {
	}
	auto operator[]() const {
		using type = decltype(ptrA.operator[](0, 0, 0) * ptrB[0, 0, 0]);
		type sum = type(0);
		for (int k = 0; k < D; k++) {
			for (int j = 0; j < D; j++) {
				for (int i = 0; i < D; i++) {
					sum += ptrA.operator[](k, j, i) * ptrB[k, j, i];
				}
			}
		}
		return sum;
	}

private:
	A ptrA;
	B ptrB;
};

template<typename A, typename B, int D, int TR>
struct Tensor3dotdot3 {

	Tensor3dotdot3(A const &a, B const &b) :
			ptrA(a), ptrB(b) {
	}
	auto operator[](auto k, int l) const {
		auto sum = ptrA.operator[](k, 0, 0) * ptrB[0, 0, l];
		if constexpr (TR == 0x02) {
			for (int j = 1; j < D; j++) {
				sum += ptrA.operator[](k, 0, j) * ptrB[j, 0, l];
				for (int i = 0; i < D; i++) {
					sum += ptrA.operator[](k, j, i) * ptrB[i, j, l];
				}
			}
		} else if constexpr (TR == -0x02) {
			for (int j = 1; j < D; j++) {
				sum += ptrA.operator[](k, 0, j) * ptrB[0, j, l];
				for (int i = 0; i < D; i++) {
					sum += ptrA.operator[](k, j, i) * ptrB[j, i, l];
				}
			}
		} else if constexpr (TR == 0x22) {
			for (int j = 1; j < D; j++) {
				sum += ptrA.operator[](0, j, k) * ptrB[j, 0, l];
				for (int i = 0; i < D; i++) {
					sum += ptrA.operator[](j, i, k) * ptrB[i, j, l];
				}
			}
		}
		return sum;
	}

private:
	A ptrA;
	B ptrB;
};

template<typename A, typename B, int D, int TR>
struct TensorNdot2 {

	TensorNdot2(A const &a, B const &b) :
			ptrA(a), ptrB(b) {
	}
	template<class ...Args>
	auto operator[](Args ...ks) const {
		static constexpr int N = sizeof...(ks) - 1;
		int index = 0;
		std::array<int, N> K;
		int i = -1;
		((((++index < N) ? K[index] : i) = ks),...);
		auto const fa = [this](auto ...k) {
			return ptrA.operator[](k...);
		};
		auto const fb = [this](int i, int j) {
			if constexpr (TR & 0x01) {
				std::swap(i, j);
			}
			return ptrB.operator[](i, j);
		};
		if constexpr (TR & 0x10) {
			auto sum = std::apply(fa, std::tuple_cat(K, std::tuple<int>(0))) * fb(0, i);
			for (int j = 1; j < D; j++) {
				sum += std::apply(fa, std::tuple_cat(K, std::tuple<int>(j))) * fb(j, i);
			}
			return sum;
		} else {
			auto sum = std::apply(fa, std::tuple_cat(std::tuple<int>(0), K)) * fb(0, i);
			for (int j = 1; j < D; j++) {
				sum += std::apply(fa, std::tuple_cat(std::tuple<int>(j), K)) * fb(j, i);
			}
			return sum;
		}
	}

private:
	A ptrA;
	B ptrB;
};

template<typename A, typename B, int D>
struct Tensor3dot2 {

	Tensor3dot2(A const &a, B const &b) :
			ptrA(a), ptrB(b) {
	}
	auto operator[](int i, int j, int k) const {
		auto sum = ptrA.operator[](i, 0, j) * ptrB.operator[](0, k);
		for (int l = 1; l < D; l++) {
			sum += ptrA.operator[](i, l, j) * ptrB.operator[](l, k);
		}
		return sum;
	}

private:
	A ptrA;
	B ptrB;
};

template<typename A, typename B, int D>
struct TensorNdotdot2Tr {

	TensorNdotdot2Tr(A const &a, B const &b) :
			ptrA(a), ptrB(b) {
	}

	auto operator[](auto ...k) const {
		auto sum = ptrA.operator[](0, k..., 0) * ptrB[0, 0];
		for (int j = 1; j < D; j++) {
			sum += ptrA.operator[](0, k..., j) * ptrB[0, j];
			for (int i = 0; i < D; i++) {
				sum += ptrA.operator[](j, k..., i) * ptrB[j, i];
			}
		}
		return sum;
	}

private:
	A ptrA;
	B ptrB;
};

template<typename A, typename B, int D, bool TR>
struct Tensor2dotdot2 {
	Tensor2dotdot2(A const &a, B const &b) :
			ptrA(a), ptrB(b) {
	}
	auto operator[]() const {
		auto sum = ptrA[0, 0] * ptrB[0, 0];
		if constexpr(TR) {
			for (int k = 1; k < D; k++) {
				sum += ptrA[k, 0] * ptrB[0, k];
				for (int n = 0; n < D; n++) {
					sum += ptrA[n, k] * ptrB[k, n];
				}
			}
			return sum;
		} else {
			for (int k = 1; k < D; k++) {
				sum += ptrA[k, 0] * ptrB[k, 0];
				for (int n = 0; n < D; n++) {
					sum += ptrA[n, k] * ptrB[n, k];
				}
			}
			return sum;
		}
	}

private:
	A ptrA;
	B ptrB;
};

template<template<typename, int, char...> typename Expr, typename A, int D, char...I>
auto operator-(Expr<A, D, I...> const &a) {
	TensorNegate<Expr<A, D, I...>> const neg(a);
	return Expr<decltype(neg), D, I...>(neg);
}

template<template<typename, int, char...> typename Expr, typename A, typename B, int D, char...I>
auto operator+(Expr<A, D, I...> const &a, Expr<B, D, I...> const &b) {
	using typeA = std::remove_reference<typename Expr<A, D, I...>::value_type>::type;
	using typeB = std::remove_reference<typename Expr<B, D, I...>::value_type>::type;
	using value_type = decltype(typeA() + typeB());
	static constexpr std::plus<value_type> o {};
	TensorBinaryOp<Expr<A, D, I...>, Expr<B, D, I...>, D, std::plus<value_type>> const sum(a, b, o);
	return Expr<decltype(sum), D, I...>(sum);
}

template<template<typename, int, char...> typename Expr, typename A, typename B, int D, char I, char J, char...K>
auto operator+(Expr<A, D, I, J, K...> const &a, Expr<B, D, J, I, K...> const &b) {
	using typeA = std::remove_reference<typename Expr<A, D, I, J, K...>::value_type>::type;
	using typeB = std::remove_reference<typename Expr<B, D, J, I, K...>::value_type>::type;
	using value_type = decltype(typeA() + typeB());
	static constexpr std::plus<value_type> o {};
	TensorBinaryOp<Expr<A, D, I, J, K...>, Expr<B, D, J, I, K...>, D, std::plus<value_type>, 0, 1> const sum(a, b, o);
	return Expr<decltype(sum), D, I, J, K...>(sum);
}

template<template<typename, int, char...> typename Expr, typename A, typename B, int D, char...I>
auto operator-(Expr<A, D, I...> const &a, Expr<B, D, I...> const &b) {
	using typeA = std::remove_reference<typename Expr<A, D, I...>::value_type>::type;
	using typeB = std::remove_reference<typename Expr<B, D, I...>::value_type>::type;
	using value_type = decltype(typeA() - typeB());
	static constexpr std::minus<value_type> o {};
	TensorBinaryOp<Expr<A, D, I...>, Expr<B, D, I...>, D, std::minus<value_type>> const sum(a, b, o);
	return Expr<decltype(sum), D, I...>(sum);
}

template<template<typename, int, char...> typename Expr, typename A, typename B, int D, char I, char J, char...K>
auto operator-(Expr<A, D, I, J, K...> const &a, Expr<B, D, J, I, K...> const &b) {
	using typeA = std::remove_reference<typename Expr<A, D, I, J, K...>::value_type>::type;
	using typeB = std::remove_reference<typename Expr<B, D, J, I, K...>::value_type>::type;
	using value_type = decltype(typeA() - typeB());
	static constexpr std::minus<value_type> o {};
	TensorBinaryOp<Expr<A, D, I, J, K...>, Expr<B, D, J, I, K...>, D, std::minus<value_type>, 0, 1> const sum(a, b, o);
	return Expr<decltype(sum), D, I, J, K...>(sum);
}

template<template<typename, int, char...> typename Expr, typename A, typename B, int D, char I, char J, char K>
auto operator-(Expr<A, D, J, K, I> const &a, Expr<B, D, I, K, J> const &b) {
	using typeA = std::remove_reference<typename Expr<A, D, J, K, I>::value_type>::type;
	using typeB = std::remove_reference<typename Expr<B, D, I, K, J>::value_type>::type;
	using value_type = decltype(typeA() - typeB());
	static constexpr std::minus<value_type> o {};
	TensorBinaryOp<Expr<A, D, J, K, I>, Expr<B, D, I, K, J>, D, std::minus<value_type>, 0, 2> const sum(a, b, o);
	return Expr<decltype(sum), D, J, K, I>(sum);
}

template<typename A, typename B, int D, char I, char J, char K>
auto operator+(Tensor3Expr<A, D, I, J, K> const &a, Tensor3Expr<B, D, I, K, J> const &b) {
	Tensor3Add<Tensor3Expr<A, D, I, J, K>, Tensor3Expr<B, D, I, K, J>, D, 0x12> const sum(a, b);
	return Tensor3Expr<decltype(sum), D, I, J, K>(sum);
}

template<typename A, typename B, int D, char I, char J, char K>
auto operator+(Tensor3Expr<A, D, I, J, K> const &a, Tensor3Expr<B, D, K, J, I> const &b) {
	Tensor3Add<Tensor3Expr<A, D, I, J, K>, Tensor3Expr<B, D, I, K, J>, D, 0x02> const sum(a, b);
	return Tensor3Expr<decltype(sum), D, I, J, K>(sum);
}

template<typename A, typename B, int D, char I, char J, char K>
auto operator+(Tensor3Expr<A, D, I, J, K> const &a, Tensor3Expr<B, D, J, I, K> const &b) {
	Tensor3Add<Tensor3Expr<A, D, I, J, K>, Tensor3Expr<B, D, J, I, K>, D, 0x01> const sum(a, b);
	return Tensor3Expr<decltype(sum), D, I, J, K>(sum);
}

/* 0 ? */
template<template<typename, int, char...> typename Expr, typename B, int D, char... I>
auto operator*(auto a, Expr<B, D, I...> const &b) {
	Tensor0xN<decltype(a), Expr<B, D, I...>> product(a, b);
	return Expr<Tensor0xN<decltype(a), Expr<B, D, I...>>, D, I...>(product);
}

/* 1 0 */
template<typename A, typename B, int D, char I, char J>
auto operator*(Tensor1Expr<A, D, I> const &a, Tensor2Expr<B, D, J, I> const &b) {
	return b * a;
}

/* 1 1 */
template<typename A, typename B, int D, char I>
auto operator*(Tensor1Expr<A, D, I> const &a, Tensor1Expr<B, D, I> const &b) {
	using type1 = Tensor1Expr<A, D, I>;
	using type2 = Tensor1Expr<B, D, I>;
	using rctype = Tensor1dot1<type1, type2, D>;
	return Tensor0Expr<rctype>(rctype(a, b));
}

template<typename A, typename B, int D, char I, char J>
auto operator*(Tensor1Expr<A, D, I> const &a, Tensor1Expr<B, D, J> const &b) {
	Tensor1xN<Tensor1Expr<A, D, I>, Tensor1Expr<B, D, J>> product(a, b);
	return Tensor2Expr<Tensor1xN<Tensor1Expr<A, D, I>, Tensor1Expr<B, D, J>>, D, I, J>(product);
}

/* 1 2 */
template<typename A, typename B, int D, char I, char J>
auto operator*(Tensor1Expr<A, D, I> const &a, Tensor2Expr<B, D, I, J> const &b) {
	using type2 = Tensor2Expr<B, D, I, J>;
	using type1 = Tensor1Expr<A, D, I>;
	using rctype = Tensor1dotN<type1, type2, D>;
	return Tensor1Expr<rctype, D, J>(rctype(a, b));
}

template<typename A, typename B, int D, char I, char J, char K>
auto operator*(Tensor1Expr<A, D, I> const &a, Tensor2Expr<B, D, J, K> const &b) {
	using type2 = Tensor2Expr<B, D, J, K>;
	using type1 = Tensor1Expr<A, D, I>;
	using rctype = Tensor1xN<type1, type2>;
	return Tensor3Expr<rctype, D, I, J, K>(rctype(a, b));
}

/* 1 3 */
template<typename A, typename B, int D, char I, char J, char K>
auto operator*(Tensor1Expr<A, D, I> const &a, Tensor3Expr<B, D, I, J, K> const &b) {
	using type1 = Tensor1Expr<A, D, I>;
	using type2 = Tensor3Expr<B, D, I, J, K>;
	using rctype = Tensor1dotN<type1, type2, D>;
	return Tensor2Expr<rctype, D, J, K>(rctype(a, b));
}

template<typename A, typename B, int D, char I, char J, char K>
auto operator*(Tensor1Expr<A, D, I> const &a, Tensor3Expr<B, D, J, K, I> const &b) {
	return b * a;
}

/* 2 1 */
template<typename A, typename B, int D, char I, char J>
auto operator*(Tensor2Expr<B, D, I, J> const &a, Tensor1Expr<A, D, J> const &b) {
	using type2 = Tensor1Expr<A, D, J>;
	using type1 = Tensor2Expr<B, D, I, J>;
	using rctype = TensorNdot1<type1, type2, D>;
	return Tensor1Expr<rctype, D, I>(rctype(a, b));
}

template<typename A, typename B, int D, char I, char J>
auto operator*(Tensor2Expr<B, D, J, I> const &a, Tensor1Expr<A, D, J> const &b) {
	return b * a;
}

template<typename A, typename B, int D, char I, char J, char K>
auto operator*(Tensor2Expr<A, D, I, J> const &a, Tensor1Expr<B, D, K> const &b) {
	TensorNx1<Tensor2Expr<A, D, I, J>, Tensor1Expr<B, D, K>> product(a, b);
	return Tensor3Expr<TensorNx1<Tensor2Expr<A, D, I, J>, Tensor1Expr<B, D, K>>, D, I, J, K>(product);
}

/* 2 2 */

template<typename A, typename B, int D, char I, char J>
auto operator*(Tensor2Expr<A, D, I, J> const &a, Tensor2Expr<B, D, J, I> const &b) {
	using type1 = Tensor2Expr<A, D, I, J>;
	using type2 = Tensor2Expr<B, D, J, I>;
	Tensor2dotdot2<type1, type2, D, true> dot(a, b);
	return Tensor0Expr<Tensor2dotdot2<type1, type2, D, true>>(dot);
}

template<typename A, typename B, int D, char I, char J>
auto operator*(Tensor2Expr<A, D, I, J> const &a, Tensor2Expr<B, D, I, J> const &b) {
	using type1 = Tensor2Expr<A, D, I, J>;
	using type2 = Tensor2Expr<B, D, I, J>;
	Tensor2dotdot2<type1, type2, D, false> dot(a, b);
	return Tensor0Expr<Tensor2dotdot2<type1, type2, D, false>>(dot);
}

template<typename A, typename B, int D, char I, char J, char K>
auto operator*(Tensor2Expr<B, D, J, I> const &a, Tensor2Expr<A, D, J, K> const &b) {
	using type1 = Tensor2Expr<B, D, J, I>;
	using type2 = Tensor2Expr<A, D, J, K>;
	using rctype = TensorNdot2<type1, type2, D, 0x00>;
	return Tensor2Expr<rctype, D, I, K>(rctype(a, b));
}

template<typename A, typename B, int D, char I, char J, char K>
auto operator*(Tensor2Expr<B, D, I, J> const &a, Tensor2Expr<A, D, J, K> const &b) {
	using type1 = Tensor2Expr<B, D, I, J>;
	using type2 = Tensor2Expr<A, D, J, K>;
	using rctype = TensorNdot2<type1, type2, D, 0x10>;
	return Tensor2Expr<rctype, D, I, K>(rctype(a, b));
}

template<typename A, typename B, int D, char I, char J, char K>
auto operator*(Tensor2Expr<B, D, J, I> const &a, Tensor2Expr<A, D, K, J> const &b) {
	using type1 = Tensor2Expr<B, D, J, I>;
	using type2 = Tensor2Expr<A, D, K, J>;
	using rctype = TensorNdot2<type1, type2, D, 0x01>;
	return Tensor2Expr<rctype, D, I, K>(rctype(a, b));
}

template<typename A, typename B, int D, char I, char J, char K>
auto operator*(Tensor2Expr<B, D, I, J> const &a, Tensor2Expr<A, D, K, J> const &b) {
	using type1 = Tensor2Expr<B, D, I, J>;
	using type2 = Tensor2Expr<A, D, K, J>;
	using rctype = TensorNdot2<type1, type2, D, 0x11>;
	return Tensor2Expr<rctype, D, I, K>(rctype(a, b));
}

/* 2 3 */
template<typename A, typename B, int D, char I, char J, char K, char L>
auto operator*(Tensor2Expr<A, D, I, J> const &a, Tensor3Expr<B, D, J, K, L> const &b) {
	using type1 = Tensor2Expr<A, D, I, J>;
	using type2 = Tensor3Expr<B, D, J, K, L>;
	using rctype = Tensor2dotN<type1, type2, D>;
	return Tensor3Expr<rctype, D, I, K, L>(rctype(a, b));
}

template<typename A, typename B, int D, char I, char J, char K, char L>
auto operator*(Tensor2Expr<A, D, I, J> const &a, Tensor3Expr<B, D, K, L, J> const &b) {
	using type1 = Tensor2Expr<A, D, I, J>;
	using type2 = Tensor3Expr<B, D, K, L, J>;
	using rctype = Tensor2dotN<type1, type2, D>;
	return Tensor3Expr<rctype, D, I, K, L>(rctype(a, b));
}

template<typename A, typename B, int D, char I, char J, char K>
auto operator*(Tensor2Expr<A, D, I, J> const &a, Tensor3Expr<B, D, J, I, K> const &b) {
	using type1 = Tensor2Expr<A, D, I, J>;
	using type2 = Tensor3Expr<B, D, J, I, K>;
	using rctype = Tensor2dotdot3<type1, type2, D>;
	return Tensor1Expr<rctype, D, K>(rctype(a, b));
}

/* 3 1 */
template<typename A, typename B, int D, char I, char J, char K>
auto operator*(Tensor3Expr<B, D, I, J, K> const &a, Tensor1Expr<A, D, K> const &b) {
	using type1 = Tensor3Expr<B, D, I, J, K>;
	using type2 = Tensor1Expr<A, D, K>;
	using rctype = TensorNdot1<type1, type2, D>;
	return Tensor2Expr<rctype, D, I, J>(rctype(a, b));
}

/* 3 2 */
template<typename A, typename B, int D, char I, char J, char K>
auto operator*(Tensor3Expr<B, D, I, J, K> const &a, Tensor2Expr<A, D, J, K> const &b) {
	using type1 = Tensor3Expr<B, D, I, J, K>;
	using type2 = Tensor2Expr<A, D, J, K>;
	using rctype = TensorNdotdot2<type1, type2, D>;
	return Tensor1Expr<rctype, D, I>(rctype(a, b));
}

template<typename A, typename B, int D, char I, char J, char K>
auto operator*(Tensor3Expr<B, D, J, I, K> const &a, Tensor2Expr<A, D, J, K> const &b) {
	using type1 = Tensor3Expr<B, D, J, I, K>;
	using type2 = Tensor2Expr<A, D, J, K>;
	using rctype = TensorNdotdot2Tr<type1, type2, D>;
	return Tensor1Expr<rctype, D, I>(rctype(a, b));
}

template<typename A, typename B, int D, char I, char J, char K, char L>
auto operator*(Tensor3Expr<B, D, I, J, K> const &a, Tensor2Expr<A, D, K, L> const &b) {
	using type1 = Tensor3Expr<B, D, I, J, K>;
	using type2 = Tensor2Expr<A, D, K, L>;
	using rctype = TensorNdot2<type1, type2, D, 0x10>;
	return Tensor3Expr<rctype, D, I, J, L>(rctype(a, b));
}

template<typename A, typename B, int D, char I, char J, char K, char L>
auto operator*(Tensor3Expr<B, D, I, J, K> const &a, Tensor2Expr<A, D, J, L> const &b) {
	using type1 = Tensor3Expr<B, D, I, J, K>;
	using type2 = Tensor2Expr<A, D, J, L>;
	using rctype = Tensor3dot2<type1, type2, D>;
	return Tensor3Expr<rctype, D, I, K, L>(rctype(a, b));
}

/* 3 3 */
template<typename A, typename B, int D, char I, char J, char K, char L>
auto operator*(Tensor3Expr<B, D, I, J, K> const &a, Tensor3Expr<A, D, J, K, L> const &b) {
	using type1 = Tensor3Expr<B, D, I, J, K>;
	using type2 = Tensor3Expr<A, D, J, K, L>;
	using rctype = Tensor3dotdot3<type1, type2, D, 0x02>;
	return Tensor2Expr<rctype, D, I, L>(rctype(a, b));
}

template<typename A, typename B, int D, char I, char J, char K, char L>
auto operator*(Tensor3Expr<B, D, J, K, I> const &a, Tensor3Expr<A, D, K, J, L> const &b) {
	using type1 = Tensor3Expr<B, D, J, K, I>;
	using type2 = Tensor3Expr<A, D, K, J, L>;
	using rctype = Tensor3dotdot3<type1, type2, D, 0x22>;
	return Tensor2Expr<rctype, D, I, L>(rctype(a, b));
}

template<typename A, typename B, int D, char I, char J, char K, char L>
auto operator*(Tensor3Expr<B, D, I, J, K> const &a, Tensor3Expr<A, D, K, J, L> const &b) {
	using type1 = Tensor3Expr<B, D, I, J, K>;
	using type2 = Tensor3Expr<A, D, K, J, L>;
	using rctype = Tensor3dotdot3<type1, type2, D, -0x02>;
	return Tensor2Expr<rctype, D, I, L>(rctype(a, b));
}

template<typename A, typename B, int D, char I, char J, char K>
auto operator*(Tensor3Expr<B, D, I, J, K> const &a, Tensor3Expr<A, D, I, J, K> const &b) {
	using type2 = Tensor3Expr<A, D, I, J, K>;
	using type1 = Tensor3Expr<B, D, I, J, K>;
	using rctype = Tensor3dotdotdot3<type1, type2, D>;
	return Tensor0Expr<rctype>(rctype(a, b));
}

/* ? 0 */
template<template<typename, int, char...> typename Expr, typename B, int D, char... I>
auto operator*(Expr<B, D, I...> const &b, auto a) {
	Tensor0xN<decltype(a), Expr<B, D, I...>> product(a, b);
	return Expr<Tensor0xN<decltype(a), Expr<B, D, I...>>, D, I...>(product);
}

}

#endif
