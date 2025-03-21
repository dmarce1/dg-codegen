/*
 * Tensor.hpp
 *
 *  Created on: Mar 2, 2025
 *      Author: dmarce1
 */

#ifndef INCLUDE_TENSOR_HPP_
#define INCLUDE_TENSOR_HPP_

namespace Tensors {

template<char C>
struct Index {
};

template<typename T>
struct Tensor0 {
	using value_type = T;
	auto operator[]() const {
		return a;
	}
	auto operator()() const;
private:
	T a;
};

template<typename T, int D>
struct Tensor1: public std::array<T, D> {
	using value_type = T;
	template<char C>
	auto operator()(Index<C>) const;
};

template<typename T, int D>
struct Tensor2 {
	using value_type = T;
	auto operator[](int i, int j) const {
		return U[D * i + j];
	}
	template<char C1, char C2>
	auto operator()(Index<C1>, Index<C2>) const;
private:
	std::array<T, D * D> U;
};

template<typename T, int D>
struct SymmetricTensor2 {
	using value_type = T;
	auto operator[](int i, int j) const {
		if (i < j) {
			std::swap(i, j);
		}
		return U[((i * (i + 1)) >> 1) + j];
	}
	template<char C1, char C2>
	auto operator()(Index<C1>, Index<C2>) const;
private:
	static constexpr int Size = D * (D + 1) / 2;
	std::array<T, Size> U;
};

template<typename T, int D>
struct Tensor3 {
	using value_type = T;
	auto operator[](int i, int j, int k) const {
		return U[D * (D * i + j) + k];
	}
	template<char C1, char C2, char C3>
	auto operator()(Index<C1>, Index<C2>, Index<C3>) const;
private:
	std::array<T, D * D * D> U;
};

template<typename A>
struct Tensor0Expr {
	using value_type = A::value_type;
	Tensor0Expr(A a) :
			ptr(a) {
	}
	auto operator[]() const {
		return ptr[0];
	}
	operator value_type() const {
		return ptr[0];
	}
private:
	A ptr;
};

template<typename A, int D, char I>
struct Tensor1Expr {
	using value_type = A::value_type;
	Tensor1Expr(A a) :
			ptr(a) {
	}
	auto operator[](int i) const {
		return ptr[i];
	}
private:
	A ptr;
};

template<typename A, int D, char I, char J>
struct Tensor2Expr {
	using value_type = A::value_type;
	Tensor2Expr(A a) :
			ptr(a) {
	}
	auto operator[](int i, int j) const {
		return ptr[i, j];
	}
private:
	A ptr;
};

template<typename A, int D, char I, char J, char K>
struct Tensor3Expr {
	using value_type = A::value_type;
	Tensor3Expr(A a) :
			ptr(a) {
	}
	auto operator[](int i, int j, int k) const {
		return ptr[i, j, k];
	}
private:
	A ptr;
};

template<typename T>
auto Tensor0<T>::operator()() const {
	return Tensor0Expr<Tensor0<T>>(&a);
}

template<typename T, int D>
template<char C>
auto Tensor1<T, D>::operator()(Index<C>) const {
	return Tensor1Expr<Tensor1<T, D>, D, C>(*this);
}

template<typename T, int D>
template<char C1, char C2>
auto Tensor2<T, D>::operator()(Index<C1>, Index<C2>) const {
	return Tensor2Expr<Tensor2<T, D>, D, C1, C2>(*this);
}

template<typename T, int D>
template<char C1, char C2>
auto SymmetricTensor2<T, D>::operator()(Index<C1>, Index<C2>) const {
	return Tensor2Expr<SymmetricTensor2<T, D>, D, C1, C2>(*this);
}

template<typename T, int D>
template<char C1, char C2, char C3>
auto Tensor3<T, D>::operator()(Index<C1>, Index<C2>, Index<C3>) const {
	return Tensor3Expr<Tensor3<T, D>, D, C1, C2, C3>(*this);
}

template<typename T, typename U>
using product_type = decltype(T()*U());

template<typename A, typename B>
struct Tensor1xN {
	using value_type = product_type<typename A::value_type, typename B::value_type>;

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
	using value_type = product_type<typename A::value_type, typename B::value_type>;

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
	using value_type = product_type<typename A::value_type, typename B::value_type>;

	Tensor1dot1(A const &a, B const &b) :
			ptrA(a), ptrB(b) {
	}

	auto operator[]() const {
		value_type sum = ptrA[0] * ptrB[0];
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
	using value_type = product_type<typename A::value_type, typename B::value_type>;

	Tensor1dotN(A const &a, B const &b) :
			ptrA(a), ptrB(b) {
	}

	auto operator[](auto ...i) const {
		value_type sum = ptrA[0] * ptrB.operator[](0, i...);
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
struct TensorNdot1 {
	using value_type = product_type<typename A::value_type, typename B::value_type>;

	TensorNdot1(A const &a, B const &b) :
			ptrA(a), ptrB(b) {
	}

	auto operator[](auto ...i) const {
		value_type sum = ptrA.operator[](i..., 0) * ptrB[0];
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
struct Tensor2dot2 {
	using value_type = product_type<typename A::value_type, typename B::value_type>;

	Tensor2dot2(A const &a, B const &b) :
			ptrA(a), ptrB(b) {
	}

	auto operator[]() const {
		value_type sum = value_type(0);
		for (int k = 0; k < D; k++) {
			for (int n = 0; n < D; n++) {
				sum += ptrA[k, n] * ptrB[k, n];
			}
		}
		return sum;
	}

private:
	A ptrA;
	B ptrB;
};

template<typename A, typename B, int D, char I>
auto operator*(Tensor1Expr<A, D, I> const &a, Tensor1Expr<B, D, I> const &b) {
	using type2 = Tensor1Expr<A, D, I>;
	using type1 = Tensor1Expr<B, D, I>;
	using rctype = Tensor1dot1<type1, type2, D>;
	return Tensor0Expr<rctype>(rctype(a, b));
}

template<typename A, typename B, int D, char I, char J>
auto operator*(Tensor1Expr<A, D, I> const &a, Tensor1Expr<B, D, J> const &b) {
	Tensor1xN<Tensor1Expr<A, D, I>, Tensor1Expr<B, D, J>> product(a, b);
	return Tensor2Expr<Tensor1xN<Tensor1Expr<A, D, I>, Tensor1Expr<B, D, J>>, D, I, J>(product);
}

template<typename A, typename B, int D, char I, char J>
auto operator*(Tensor1Expr<A, D, I> const &a, Tensor2Expr<B, D, I, J> const &b) {
	using type2 = Tensor2Expr<B, D, I, J>;
	using type1 = Tensor1Expr<A, D, I>;
	using rctype = Tensor1dotN<type1, type2, D>;
	return Tensor2Expr<rctype, D, I, J>(rctype(a, b));
}

template<typename A, typename B, int D, char I, char J>
auto operator*(Tensor2Expr<B, D, I, J> const &a, Tensor1Expr<A, D, J> const &b) {
	using type2 = Tensor1Expr<A, D, J>;
	using type1 = Tensor2Expr<B, D, I, J>;
	using rctype = TensorNdot1<type1, type2, D>;
	return Tensor2Expr<rctype, D, I, J>(rctype(a, b));
}

template<typename A, typename B, int D, char I, char J>
auto operator*(Tensor1Expr<A, D, I> const &a, Tensor2Expr<B, D, J, I> const &b) {
	return b * a;
}

template<typename A, typename B, int D, char I, char J>
auto operator*(Tensor2Expr<B, D, J, I> const &a, Tensor1Expr<A, D, J> const &b) {
	return b * a;
}

template<typename A, typename B, int D, char I, char J>
auto operator*(Tensor2Expr<B, D, J, I> const &a, Tensor2Expr<B, D, J, I> const &b) {
	using type2 = Tensor2Expr<A, D, I, J>;
	using type1 = Tensor2Expr<B, D, I, J>;
	using rctype = Tensor2dot2<type1, type2, D>;
	return Tensor0Expr<rctype>(rctype(a, b));
}

template<typename A, typename B, int D, char I, char J, char K>
auto operator*(Tensor1Expr<A, D, I> const &a, Tensor3Expr<B, D, I, J, K> const &b) {
	using type1 = Tensor1Expr<A, D, I>;
	using type2 = Tensor3Expr<B, D, I, J, K>;
	using rctype = Tensor1dotN<type1, type2, D>;
	return Tensor2Expr<rctype, D, J, K>(rctype(a, b));
}

template<typename A, typename B, int D, char I, char J, char K>
auto operator*(Tensor3Expr<B, D, I, J, K> const &a, Tensor1Expr<A, D, K> const &b) {
	using type1 = Tensor3Expr<B, D, I, J, K>;
	using type2 = Tensor1Expr<A, D, K>;
	using rctype = TensorNdot1<type1, type2, D>;
	return Tensor2Expr<rctype, D, I, J>(rctype(a, b));
}

}

#endif
