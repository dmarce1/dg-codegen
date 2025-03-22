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
#include <optional>

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
	auto operator()() const {
		return a;
	}
	auto getHandle() const {
		return [this](int i) -> T& {
			return a;
		};
	}
private:
	T a;
};

template<typename T, int D>
struct Tensor<T, D, 1> {
	using value_type = T;
	auto operator()(int i) const {
		return U[i];
	}
	auto& operator()(int i) {
		return U[i];
	}
	template<char C>
	auto operator()(Index<C>);
	template<char C>
	auto operator()(Index<C>) const;
	auto getHandle() {
		return [this](int i) -> T& {
			return U[i];
		};
	}
	auto getHandle() const {
		return [this](int i) -> T const& {
			return U[i];
		};
	}
private:
	std::array<T, D> U;
};

template<typename T, int D, int ...S>
struct Tensor<T, D, 2, S...> {
	using value_type = T;
	static constexpr bool isSymmetric = sizeof...(S);
	static constexpr int computeIndex(int i, int j) {
		if constexpr (isSymmetric) {
			if (i < j) {
				std::swap(i, j);
			}
			return ((i * (i + 1)) >> 1) + j;
		} else {
			return D * i + j;
		}
	}
	auto const& operator()(int i, int j) const {
		return U[computeIndex(i, j)];
	}
	auto& operator()(int i, int j) {
		return U[computeIndex(i, j)];
	}
	template<char A, char B> auto operator()(Index<A>, Index<B>);
	template<char A> auto operator()(int, Index<A>);
	template<char A> auto operator()(Index<A>, int);
	template<char A, char B> auto operator()(Index<A>, Index<B>) const;
	template<char A> auto operator()(int, Index<A>) const;
	template<char A> auto operator()(Index<A>, int) const;
	auto getHandle() {
		return [this](int i, int j) -> T& {
			return this->operator()(i, j);
		};
	}
	std::function<T& (int)> getHandle(std::optional<int> I, std::optional<int> J) {
		if (bool(I)) {
			return [this, I](int j) -> T& {
				return this->operator()(I.value(), j);
			};
		} else {
			return [this, J](int i) -> T& {
				return this->operator()(i, J.value());
			};
		}
	}
	auto getHandle() const {
		return const_cast<Tensor&>(*this).getHandle();
	}
	std::function<T const& (int)> getHandle(std::optional<int> I, std::optional<int> J) const {
		return const_cast<Tensor&>(*this).getHandle(I, J);
	}
	Tensor<T, D - 1, 2> subTensor(int j, int k) const {
		Tensor<T, D - 1, 2> sub;
		for (int n = 0; n < j; n++) {
			for (int m = 0; m < k; m++) {
				sub(n, m) = (*this)(n, m);
			}
			for (int m = k + 1; m < D; m++) {
				sub(n, m - 1) = (*this)(n, m);
			}
		}
		for (int n = j + 1; n < D; n++) {
			for (int m = 0; m < k; m++) {
				sub(n - 1, m) = (*this)(n, m);
			}
			for (int m = k + 1; m < D; m++) {
				sub(n - 1, m - 1) = (*this)(n, m);
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
				I(n, m) = sgn * subTensor(n, m).determinant();
				sgn = -sgn;
			}
			sgn_n = -sgn_n;
		}
		for (int n = 0; n < D; n++) {
			det += I(0, n);
		}
		T const dinv = T(1) / det;
		for (int n = 0; n < D; n++) {
			for (int m = 0; m < D; m++) {
				I(n, m) *= dinv;
			}
		}
		return I;
	}
private:
	static constexpr int Size = isSymmetric ? (D * (D + 1) / 2) : (D * D);
	std::array<T, Size> U;
}
;

template<typename T, int D>
struct Unit {
	using value_type = T;
	Unit(int j_) :
			j(j_) {
	}
	auto operator()(int i) const {
		return T(i == j);
	}
	template<char C>
	auto operator()(Index<C>) const;
	auto getHandle() const {
		return [this](int i) {
			return T(i == j);
		};
	}
private:
	int const j;
};

template<typename T, int D>
struct Delta {
	using value_type = T;
	auto operator()(int i, int j) const {
		return T(i == j);
	}
	template<char C1, char C2>
	auto operator()(Index<C1>, Index<C2>) const;
	template<char C>
	auto operator()(int i, Index<C>) const;
	template<char C>
	auto operator()(Index<C>, int i) const;
	auto getHandle() const {
		return [this](int i, int j) {
			return T(i == j);
		};
	}
};

template<typename T, int D>
struct Tensor<T, D, 3> {
	using value_type = T;
	auto const& operator()(int i, int j, int k) const {
		return U[D * (D * i + j) + k];
	}
	auto& operator()(int i, int j, int k) {
		return U[D * (D * i + j) + k];
	}
	template<char C1, char C2, char C3>
	auto operator()(Index<C1>, Index<C2>, Index<C3>);
	template<char C1, char C2, char C3>
	auto operator()(Index<C1>, Index<C2>, Index<C3>) const;
	auto getHandle() {
		return [this](int i, int j, int k) -> const T& {
			return this->operator()(i, j, k);
		};
	}
	auto getHandle() const {
		return const_cast<Tensor&>(*this).getHandle();
	}
private:
	std::array<T, D * D * D> U;
};

template<typename T, int D>
struct Tensor<T, D, 3, Sym<0, 1>> {
	using value_type = T;
	auto operator()(int i, int j, int k) const {
		return U[i * D2 + ((j * (j + 1)) >> 1) + k];
	}
	auto& operator()(int i, int j, int k) {
		return U[i * D2 + ((j * (j + 1)) >> 1) + k];
	}
	template<char C1, char C2, char C3> auto operator()(Index<C1>, Index<C2>, Index<C3>);
	template<char C1, char C2, char C3> auto operator()(Index<C1>, Index<C2>, Index<C3>) const;
	template<char C2, char C3> auto operator()(int i, Index<C2>, Index<C3>);
	auto getHandle() {
		return [this](int i, int j, int k) -> T& {
			return this->operator()(i, j, k);
		};
	}
	std::function<T const& (int, int)> getHandle(std::optional<int> I, std::optional<int> J, std::optional<int> K) const {
		if (bool(I)) {
			return [this, I](int j, int k) -> T const& {
				return this->operator()(I.value(), j, k);
			};
		} else if (bool(J)) {
			return [this, J](int i, int k) -> T const& {
				return this->operator()(i, J.value(), k);
			};
		} else {
			return [this, K](int i, int j) -> T const& {
				return this->operator()(i, j, K.value());
			};
		}
	}
	auto getHandle() const {
		return const_cast<Tensor&>(*this).getHandle();
	}
	std::function<T& (int, int)> getHandle(std::optional<int> I, std::optional<int> J, std::optional<int> K) {
		return const_cast<Tensor&>(*this).getHandle(I, J, K);
	}
private:
	static constexpr int D2 = D * (D + 1) / 2;
	std::array<T, D * D2> U;
};

template<typename A>
class Tensor0Expr {
	A ptr;
public:
	using value_type = decltype(ptr.operator()());
	Tensor0Expr(A a) :
			ptr(a) {
	}
	auto operator()() const {
		return ptr.operator()();
	}
	auto operator=(Tensor0Expr const &other) {
		ptr.operator()() = other.operator()();
		return *this;
	}
	operator value_type() const {
		return ptr.operator()();
	}
};

template<typename A, int D, char I>
class Tensor1Expr {
	A ptr;
public:
	using value_type = decltype(ptr(0));
	Tensor1Expr(A a) :
			ptr(a) {
	}
	auto operator()(int i) const {
		return ptr(i);
	}
	auto& operator()(int i) {
		return ptr(i);
	}
	template<typename U>
	auto operator=(U const &other) {
		if constexpr (std::is_convertible<value_type, U>::value) {
			for (int i = 0; i < D; i++) {
				ptr(i) = other;
			}
		} else {
			for (int i = 0; i < D; i++) {
				ptr(i) = other(i);
			}
		}
		return *this;
	}
	template<typename U>
	auto operator+=(U const &other) {
		if constexpr (std::is_convertible<value_type, U>::value) {
			for (int i = 0; i < D; i++) {
				ptr(i) += other;
			}
		} else {
			for (int i = 0; i < D; i++) {
				ptr(i) += other(i);
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
	using value_type = decltype(ptr(0, 0));
	Tensor2Expr(A a) :
			ptr(a) {
	}
	constexpr auto operator()(int i, int j) const {
		return ptr(i, j);
	}
	template<template<typename, int, char, char> class Expr, typename B>
	auto operator=(Expr<B, D, I, J> const &other) {
		if constexpr (isSymmetric2<A>::value) {
			for (int i = 0; i < D; i++) {
				for (int j = 0; j <= i; j++) {
					ptr(i, j) = other(i, j);
				}
			}
		} else {
			for (int i = 0; i < D; i++) {
				for (int j = 0; j < D; j++) {
					ptr(i, j) = other(i, j);
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
					ptr(i, j) = other;
				}
			}
		} else {
			for (int i = 0; i < D; i++) {
				for (int j = 0; j < D; j++) {
					ptr(i, j) = other;
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
					ptr(i, j) += other(i, j);
				}
			}
		} else {
			for (int i = 0; i < D; i++) {
				for (int j = 0; j < D; j++) {
					ptr(i, j) += other(i, j);
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
					ptr(i, j) -= other(i, j);
				}
			}
		} else {
			for (int i = 0; i < D; i++) {
				for (int j = 0; j < D; j++) {
					ptr(i, j) -= other(i, j);
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
					ptr(i, j) -= other(j, i);
				}
			}
		} else {
			for (int i = 0; i < D; i++) {
				for (int j = 0; j < D; j++) {
					ptr(i, j) -= other(j, i);
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
					ptr(i, j) += other(j, i);
				}
			}
		} else {
			for (int i = 0; i < D; i++) {
				for (int j = 0; j < D; j++) {
					ptr(i, j) += other(j, i);
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
	using value_type = decltype(ptr(0,0,0));
	Tensor3Expr(A a) :
			ptr(a) {
	}
	auto operator()(int i, int j, int k) const {
		return ptr(i, j, k);
	}
	template<typename B>
	auto operator=(B const &other) {
		for (int i = 0; i < D; i++) {
			for (int j = 0; j < D; j++) {
				for (int k = 0; k < D; k++) {
					ptr(i, j, k) = other;
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
					ptr(i, j, k) = other(i, j, k);
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
					ptr(i, j, k) = other(j, k, i);
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
					ptr(i, j, k) = other(k, j, i);
				}
			}
		}
		return *this;
	}
};

template<typename T, int D>
template<char C>
auto Tensor<T, D, 1>::operator()(Index<C>) {
	return Tensor1Expr<decltype(getHandle()), D, C>(getHandle());
}

template<typename T, int D>
template<char C>
auto Tensor<T, D, 1>::operator()(Index<C>) const {
	return Tensor1Expr<decltype(getHandle()), D, C>(getHandle());
}

template<typename T, int D, int ...S>
template<char C1, char C2>
auto Tensor<T, D, 2, S...>::operator()(Index<C1>, Index<C2>) {
	return Tensor2Expr<decltype(getHandle()), D, C1, C2>(getHandle());
}

template<typename T, int D, int ...S>
template<char A>
auto Tensor<T, D, 2, S...>::operator()(int i, Index<A>) {
	return Tensor1Expr<decltype(getHandle(i, std::nullopt)), D, A>(getHandle(i, std::nullopt));
}

template<typename T, int D, int ...S>
template<char A>
auto Tensor<T, D, 2, S...>::operator()(Index<A>, int i) {
	return Tensor1Expr<decltype(getHandle(std::nullopt, i)), D, A>(getHandle(std::nullopt, i));
}

template<typename T, int D, int ...S>
template<char C1, char C2>
auto Tensor<T, D, 2, S...>::operator()(Index<C1>, Index<C2>) const {
	return Tensor2Expr<decltype(getHandle()), D, C1, C2>(getHandle());
}

template<typename T, int D, int ...S>
template<char A>
auto Tensor<T, D, 2, S...>::operator()(int i, Index<A>) const {
	return Tensor1Expr<decltype(getHandle(i, std::nullopt)), D, A>(getHandle(i, std::nullopt));
}

template<typename T, int D, int ...S>
template<char A>
auto Tensor<T, D, 2, S...>::operator()(Index<A>, int i) const {
	return Tensor1Expr<decltype(getHandle(std::nullopt, i)), D, A>(getHandle(std::nullopt, i));
}

template<typename T, int D>
template<char C1, char C2>
auto Delta<T, D>::operator()(Index<C1>, Index<C2>) const {
	return Tensor2Expr<Delta<T, D>, D, C1, C2>(*this);
}

template<typename T, int D>
template<char C1, char C2, char C3>
auto Tensor<T, D, 3>::operator()(Index<C1>, Index<C2>, Index<C3>) {
	return Tensor3Expr<decltype(getHandle()), D, C1, C2, C3>(getHandle());
}

template<typename T, int D>
template<char C1, char C2, char C3>
auto Tensor<T, D, 3>::operator()(Index<C1>, Index<C2>, Index<C3>) const {
	return Tensor3Expr<decltype(getHandle()), D, C1, C2, C3>(getHandle());
}

template<typename T, int D>
template<char C1, char C2, char C3>
auto Tensor<T, D, 3, Sym<0, 1>>::operator()(Index<C1>, Index<C2>, Index<C3>) {
	return Tensor3Expr<decltype(getHandle()), D, C1, C2, C3>(getHandle());
}

template<typename T, int D>
template<char C1, char C2, char C3>
auto Tensor<T, D, 3, Sym<0, 1>>::operator()(Index<C1>, Index<C2>, Index<C3>) const {
	return Tensor3Expr<decltype(getHandle()), D, C1, C2, C3>(getHandle());
}

template<typename T, int D>
template<char C2, char C3>
auto Tensor<T, D, 3, Sym<0, 1>>::operator()(int i, Index<C2>, Index<C3>) {
	return Tensor2Expr<decltype(getHandle(i, std::nullopt, std::nullopt)), D, C2, C3>(getHandle(i, std::nullopt, std::nullopt));
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

template<template<typename, int, char...> typename Expr, typename T, int D, char...I>
auto operator-(Expr<T, D, I...> const &A) {
	auto const f = [A](auto... is) {
		return -A(is...);
	};
	return Expr<decltype(f), D, I...>(f);
}

template<typename T, typename U, int D, char I>
auto operator+(Tensor1Expr<T, D, I> const &A, Tensor1Expr<U, D, I> const &B) {
	auto const f = [A, B](int i) {
		return A(i) + B(i);
	};
	return Tensor1Expr<decltype(f), D, I>(f);
}

template<typename T, typename U, int D, char I>
auto operator-(Tensor1Expr<T, D, I> const &A, Tensor1Expr<U, D, I> const &B) {
	auto const f = [A, B](int i) {
		return A(i) - B(i);
	};
	return Tensor1Expr<decltype(f), D, I>(f);
}

template<typename T, typename U, int D, char I, char J, char K, char L>
auto operator+(Tensor2Expr<T, D, I, J> const &A, Tensor2Expr<U, D, K, L> const &B) {
	if constexpr (K == I && L == J) {
		auto const f = [A, B](int i, int j) {
			return A(i, j) + B(i, j);
		};
		return Tensor2Expr<decltype(f), D, I, J>(f);
	} else if constexpr (K == J && L == I) {
		auto const f = [A, B](int i, int j) {
			return A(i, j) + B(j, i);
		};
		return Tensor2Expr<decltype(f), D, I, J>(f);
	}
}

template<typename T, typename U, int D, char I, char J, char K, char L>
auto operator-(Tensor2Expr<T, D, I, J> const &A, Tensor2Expr<U, D, K, L> const &B) {
	if constexpr (K == I && L == J) {
		auto const f = [A, B](int i, int j) {
			return A(i, j) - B(i, j);
		};
		return Tensor2Expr<decltype(f), D, I, J>(f);
	} else if constexpr (K == J && L == I) {
		auto const f = [A, B](int i, int j) {
			return A(i, j) - B(j, i);
		};
		return Tensor2Expr<decltype(f), D, I, J>(f);
	}
}

template<typename T, typename U, int D, char I, char J, char K, char L, char M, char N>
auto operator+(Tensor3Expr<T, D, I, J, K> const &A, Tensor3Expr<U, D, L, M, N> const &B) {
	if constexpr (L == I && M == J && N == K) {
		auto const f = [A, B](int i, int j, int k) {
			return A(i, j, k) + B(i, j, k);
		};
		return Tensor3Expr<decltype(f), D, I, J, K>(f);
	} else if constexpr (L == J && M == K && N == I) {
		auto const f = [A, B](int i, int j, int k) {
			return A(i, j, k) + B(j, k, i);
		};
		return Tensor3Expr<decltype(f), D, I, J, K>(f);
	} else if constexpr (L == K && M == I && N == J) {
		auto const f = [A, B](int i, int j, int k) {
			return A(i, j, k) + B(k, i, j);
		};
		return Tensor3Expr<decltype(f), D, I, J, K>(f);
	} else if constexpr (L == J && M == I && N == K) {
		auto const f = [A, B](int i, int j, int k) {
			return A(i, j, k) + B(j, i, k);
		};
		return Tensor3Expr<decltype(f), D, I, J, K>(f);
	} else if constexpr (L == I && M == K && N == J) {
		auto const f = [A, B](int i, int j, int k) {
			return A(i, j, k) + B(i, k, j);
		};
		return Tensor3Expr<decltype(f), D, I, J, K>(f);
	} else if constexpr (L == K && M == J && N == I) {
		auto const f = [A, B](int i, int j, int k) {
			return A(i, j, k) + B(k, j, i);
		};
		return Tensor3Expr<decltype(f), D, I, J, K>(f);
	}
}

template<typename T, typename U, int D, char I, char J, char K, char L, char M, char N>
auto operator-(Tensor3Expr<T, D, I, J, K> const &A, Tensor3Expr<U, D, L, M, N> const &B) {
	if constexpr (L == I && M == J && N == K) {
		auto const f = [A, B](int i, int j, int k) {
			return A(i, j, k) - B(i, j, k);
		};
		return Tensor3Expr<decltype(f), D, I, J, K>(f);
	} else if constexpr (L == J && M == K && N == I) {
		auto const f = [A, B](int i, int j, int k) {
			return A(i, j, k) - B(j, k, i);
		};
		return Tensor3Expr<decltype(f), D, I, J, K>(f);
	} else if constexpr (L == K && M == I && N == J) {
		auto const f = [A, B](int i, int j, int k) {
			return A(i, j, k) - B(k, i, j);
		};
		return Tensor3Expr<decltype(f), D, I, J, K>(f);
	} else if constexpr (L == J && M == I && N == K) {
		auto const f = [A, B](int i, int j, int k) {
			return A(i, j, k) - B(j, i, k);
		};
		return Tensor3Expr<decltype(f), D, I, J, K>(f);
	} else if constexpr (L == I && M == K && N == J) {
		auto const f = [A, B](int i, int j, int k) {
			return A(i, j, k) - B(i, k, j);
		};
		return Tensor3Expr<decltype(f), D, I, J, K>(f);
	} else if constexpr (L == K && M == J && N == I) {
		auto const f = [A, B](int i, int j, int k) {
			return A(i, j, k) - B(k, j, i);
		};
		return Tensor3Expr<decltype(f), D, I, J, K>(f);
	}
}

/* 0 ? */
template<template<typename, int, char...> typename Expr, typename T, int D, char... I>
auto operator*(auto a, Expr<T, D, I...> const &B) {
	auto const f = [a, B](auto...is) {
		return a * B(is...);
	};
	return Expr<decltype(f), D, I...>(f);
}

template<template<typename, int, char...> typename Expr, typename T, int D, char... I>
auto operator*(Expr<T, D, I...> const &B, auto a) {
	return a * B;
}

/* 1 1 */

template<typename T, typename U, int D, char I, char J>
auto operator*(Tensor1Expr<T, D, I> const &A, Tensor1Expr<U, D, J> const &B) {
	if constexpr (I == J) {
		auto const f = [A, B]() {
			using type = decltype(A(0) * B(0));
			type sum = type(0);
			for (int j = 0; j < D; j++) {
				sum += A(j) * B(j);
			}
			return sum;
		};
		return Tensor0Expr<decltype(f)>(f);
	} else {
		auto const f = [A, B](int i, int j) {
			return A(i) * B(j);
		};
		return Tensor2Expr<decltype(f), D, I, J>(f);
	}
}

/* 1 2 */

template<typename T, typename U, int D, char I, char J, char K>
auto operator*(Tensor1Expr<T, D, I> const &A, Tensor2Expr<U, D, J, K> const &B) {
	if constexpr (I == J) {
		auto const f = [A, B](int k) {
			using type = decltype(A(0) * B(0, 0));
			type sum = type(0);
			for (int j = 0; j < D; j++) {
				sum += A(j) * B(j, k);
			}
			return sum;
		};
		return Tensor1Expr<decltype(f), D, K>(f);
	} else if constexpr (I == K) {
		auto const f = [A, B](int j) {
			using type = decltype(A(0) * B(0, 0));
			type sum = type(0);
			for (int k = 0; k < D; k++) {
				sum += A(k) * B(j, k);
			}
			return sum;
		};
		return Tensor1Expr<decltype(f), D, J>(f);
	} else {
		auto const f = [A, B](int i, int j, int k) {
			return A(i) * B(j, k);
		};
		return Tensor3Expr<decltype(f), D, I, J, K>(f);
	}
}

/* 1 3 */

template<typename T, typename U, int D, char I, char J, char K, char L>
auto operator*(Tensor1Expr<T, D, I> const &A, Tensor3Expr<U, D, J, K, L> const &B) {
	if constexpr (I == J) {
		auto const f = [A, B](int k, int l) {
			using type = decltype(A(0) * B(0, 0, 0));
			type sum = type(0);
			for (int j = 0; j < D; j++) {
				sum += A(j) * B(j, k, l);
			}
			return sum;
		};
		return Tensor2Expr<decltype(f), D, K, L>(f);
	} else if constexpr (I == K) {
		auto const f = [A, B](int j, int l) {
			using type = decltype(A(0) * B(0, 0, 0));
			type sum = type(0);
			for (int k = 0; k < D; k++) {
				sum += A(k) * B(j, k, l);
			}
			return sum;
		};
		return Tensor2Expr<decltype(f), D, J, L>(f);
	} else if constexpr (I == L) {
		auto const f = [A, B](int j, int k) {
			using type = decltype(A(0) * B(0, 0, 0));
			type sum = type(0);
			for (int l = 0; l < D; l++) {
				sum += A(l) * B(j, k, l);
			}
			return sum;
		};
		return Tensor2Expr<decltype(f), D, J, K>(f);
	}
}

/* 2 1 */
template<typename T, typename U, int D, char I, char J, char K>
auto operator*(Tensor2Expr<T, D, I, J> const &A, Tensor1Expr<U, D, K> const &B) {
	if constexpr (J == K) {
		auto const f = [A, B](int i) {
			using type = decltype(A(0, 0) * B(0));
			type sum = type(0);
			for (int j = 0; j < D; j++) {
				sum += A(i, j) * B(j);
			}
			return sum;
		};
		return Tensor1Expr<decltype(f), D, I>(f);
	} else if constexpr (I == K) {
		auto const f = [A, B](int j) {
			using type = decltype(A(0, 0) * B(0));
			type sum = type(0);
			for (int i = 0; i < D; i++) {
				sum += A(i, j) * B(i);
			}
			return sum;
		};
		return Tensor1Expr<decltype(f), D, J>(f);
	} else {
		auto const f = [A, B](int i, int j, int k) {
			return A(i, j) * B(k);
		};
		return Tensor3Expr<decltype(f), D, I, J, K>(f);
	}
}

/* 2 2 */
template<typename T, typename U, int D, char I, char J, char K, char L>
auto operator*(Tensor2Expr<T, D, I, J> const &A, Tensor2Expr<U, D, K, L> const &B) {
	if constexpr (I == K && J == L) {
		auto const f = [A, B]() {
			using type = decltype(A(0, 0) * B(0, 0));
			type sum = type(0);
			for (int i = 0; i < D; i++) {
				for (int j = 0; j < D; j++) {
					sum += A(j, i) * B(i, j);
				}
			}
			return sum;
		};
		return Tensor0Expr<decltype(f)>(f);
	} else if constexpr (I == L && J == K) {
		auto const f = [A, B]() {
			using type = decltype(A(0, 0) * B(0, 0));
			type sum = type(0);
			for (int i = 0; i < D; i++) {
				for (int j = 0; j < D; j++) {
					sum += A(j, i) * B(j, i);
				}
			}
			return sum;
		};
		return Tensor0Expr<decltype(f)>(f);
	} else if constexpr (J == L) {
		auto const f = [A, B](int i, int k) {
			using type = decltype(A(0, 0) * B(0, 0));
			type sum = type(0);
			for (int j = 0; j < D; j++) {
				sum += A(i, j) * B(k, j);
			}
			return sum;
		};
		return Tensor2Expr<decltype(f), D, I, K>(f);
	} else if constexpr (I == L) {
		auto const f = [A, B](int i, int k) {
			using type = decltype(A(0, 0) * B(0, 0));
			type sum = type(0);
			for (int j = 0; j < D; j++) {
				sum += A(j, i) * B(k, j);
			}
			return sum;
		};
		return Tensor2Expr<decltype(f), D, J, K>(f);
	} else if constexpr (J == K) {
		auto const f = [A, B](int i, int k) {
			using type = decltype(A(0, 0) * B(0, 0));
			type sum = type(0);
			for (int j = 0; j < D; j++) {
				sum += A(i, j) * B(j, k);
			}
			return sum;
		};
		return Tensor2Expr<decltype(f), D, I, L>(f);
	} else if constexpr (I == K) {
		auto const f = [A, B](int i, int k) {
			using type = decltype(A(0, 0) * B(0, 0));
			type sum = type(0);
			for (int j = 0; j < D; j++) {
				sum += A(j, i) * B(j, k);
			}
			return sum;
		};
		return Tensor2Expr<decltype(f), D, J, L>(f);
	}
}

/* 2 3 */
template<typename T, typename U, int D, char I, char J, char K, char L, char M>
auto operator*(Tensor2Expr<T, D, I, J> const &A, Tensor3Expr<U, D, K, L, M> const &B) {
	if constexpr (I == L && J == K) {
		auto const f = [A, B](int k) {
			using type = decltype(A(0, 0) * B(0, 0, 0));
			type sum = type(0);
			for (int i = 0; i < D; i++) {
				for (int j = 0; j < D; j++) {
					sum += A(i, j) * B(j, i, k);
				}
			}
			return sum;
		};
		return Tensor1Expr<decltype(f), D, M>(f);
	} else if constexpr (K == J) {
		auto const f = [A, B](int i, int k, int l) {
			using type = decltype(A(0, 0) * B(0, 0, 0));
			type sum = type(0);
			for (int j = 0; j < D; j++) {
				sum += A(i, j) * B(j, k, l);
			}
			return sum;
		};
		return Tensor3Expr<decltype(f), D, I, L, M>(f);
	} else if constexpr (M == J) {
		auto const f = [A, B](int i, int k, int l) {
			using type = decltype(A(0, 0) * B(0, 0, 0));
			type sum = type(0);
			for (int j = 0; j < D; j++) {
				sum += A(i, j) * B(k, l, j);
			}
			return sum;
		};
		return Tensor3Expr<decltype(f), D, I, K, L>(f);
	}
}

/* 3 1 */
template<typename T, typename U, int D, char I, char J, char K>
auto operator*(Tensor3Expr<T, D, I, J, K> const &A, Tensor1Expr<U, D, K> const &B) {
	auto const f = [A, B](int i, int j) {
		using type = decltype(A(0, 0, 0) * B(0));
		type sum = type(0);
		for (int k = 0; k < D; k++) {
			sum += A(i, j, k) * B(k);
		}
		return sum;
	};
	return Tensor2Expr<decltype(f), D, I, J>(f);
}

/* 3 2 */
template<typename T, typename U, int D, char I, char J, char K, char L, char M>
auto operator*(Tensor3Expr<T, D, I, J, K> const &A, Tensor2Expr<U, D, L, M> const &B) {
	if constexpr (L == J && M == K) {
		auto const f = [A, B](int i) {
			using type = decltype(A(0, 0, 0) * B(0, 0));
			type sum = type(0);
			for (int j = 0; j < D; j++) {
				for (int k = 0; k < D; k++) {
					sum += A(i, j, k) * B(j, k);
				}
			}
			return sum;
		};
		return Tensor1Expr<decltype(f), D, I>(f);
	} else if constexpr (L == I && M == K) {
		auto const f = [A, B](int j) {
			using type = decltype(A(0, 0, 0) * B(0, 0));
			type sum = type(0);
			for (int i = 0; i < D; i++) {
				for (int k = 0; k < D; k++) {
					sum += A(i, j, k) * B(i, k);
				}
			}
			return sum;
		};
		return Tensor1Expr<decltype(f), D, J>(f);
	} else if constexpr (L == J) {
		auto const f = [A, B](int i, int k, int l) {
			using type = decltype(A(0, 0, 0) * B(0, 0));
			type sum = type(0);
			for (int j = 0; j < D; j++) {
				sum += A(i, j, k) * B(j, l);
			}
			return sum;
		};
		return Tensor3Expr<decltype(f), D, I, K, M>(f);
	} else if constexpr (L == K) {
		auto const f = [A, B](int i, int j, int l) {
			using type = decltype(A(0, 0, 0) * B(0, 0));
			type sum = type(0);
			for (int k = 0; k < D; k++) {
				sum += A(i, j, k) * B(k, l);
			}
			return sum;
		};
		return Tensor3Expr<decltype(f), D, I, J, M>(f);
	}
}

/* 3 3 */
template<typename T, typename U, int D, char I, char J, char K, char L, char M, char N>
auto operator*(Tensor3Expr<T, D, I, J, K> const &A, Tensor3Expr<U, D, L, M, N> const &B) {
	if constexpr (I == L && J == M && K == N) {
		auto const f = [A, B]() {
			using type = decltype(A(0, 0, 0) * B(0, 0, 0));
			type sum = type(0);
			for (int k = 0; k < D; k++) {
				for (int j = 0; j < D; j++) {
					for (int i = 0; i < D; i++) {
						sum += A(i, j, k) * B(i, j, k);
					}
				}
			}
			return sum;
		};
		return Tensor0Expr<decltype(f)>(f);
	} else if constexpr (L == K && M == J) {
		auto const f = [A, B](int i, int l) {
			using type = decltype(A(0, 0, 0) * B(0, 0, 0));
			type sum = type(0);
			for (int k = 0; k < D; k++) {
				for (int j = 0; j < D; j++) {
					sum += A(i, j, k) * B(k, j, l);
				}
			}
			return sum;
		};
		return Tensor2Expr<decltype(f), D, I, N>(f);
	} else if constexpr (L == J && M == I) {
		auto const f = [A, B](int k, int l) {
			using type = decltype(A(0, 0, 0) * B(0, 0, 0));
			type sum = type(0);
			for (int i = 0; i < D; i++) {
				for (int j = 0; j < D; j++) {
					sum += A(i, j, k) * B(j, i, l);
				}
			}
			return sum;
		};
		return Tensor2Expr<decltype(f), D, K, N>(f);
	} else if constexpr (L == J && M == K) {
		auto const f = [A, B](int i, int l) {
			using type = decltype(A(0, 0, 0) * B(0, 0, 0));
			type sum = type(0);
			for (int k = 0; k < D; k++) {
				for (int j = 0; j < D; j++) {
					sum += A(i, j, k) * B(j, k, l);
				}
			}
			return sum;
		};
		return Tensor2Expr<decltype(f), D, I, N>(f);
	}
}

}

#endif
