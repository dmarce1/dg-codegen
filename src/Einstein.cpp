/*
 * Einstein.cpp
 *
 *  Created on: Feb 25, 2025
 *      Author: dmarce1
 */

#include "Vector.hpp"

constexpr int pow(int x, int n) {
	int y = 1;
	for (int m = 0; m < n; m++) {
		y *= x;
	}
	return y;
}

enum class symmetry_t : int {
	asymmetric, symmetric, antisymmetric
};

template<int R, int D>
auto symmetricMapping() {
	std::array<int, R> indices;
	std::array<int, pow(D, R)> map;
	std::fill(indices.begin(), indices.end(), 0);
	std::fill(map.begin(), map.end(), -1);
	auto const size = pow(D, R);
	int next = 0;
	for (int l = 0; l < size; l++) {
		auto common = indices;
		std::sort(common.begin(), common.end());
		int k = 0;
		for (int d = 0; d < R; d++) {
			k = D * k + common[d];
		}
		if (map[k] < 0) {
			map[k] = next++;
		}
		map[l] = map[k];
		int dim = 0;
		while ((++indices[dim] == D) && (dim < R)) {
			indices[dim++] = 0;
		}
	}
	std::fill(indices.begin(), indices.end(), 0);
	for (int l = 0; l < size; l++) {
		int dim = 0;
		while ((++indices[dim] == D) && (dim < R)) {
			indices[dim++] = 0;
		}
	}
}

template<typename T, int R, int D, symmetry_t SYMMETRY = symmetry_t::asymmetric>
class Tensor;

template<typename DERIVED, typename T, int S, int R, int D>
struct TensorBase: public Math::Vector<T, S> {
	using REAL = T;
	static constexpr int SIZE = S;
	static constexpr int RANK = R;
	static constexpr int DIMS = D;
	using base_type = Math::Vector<REAL, SIZE>;

	REAL const& operator[](auto ...i) const {
		return static_cast<DERIVED*>(this)->access(i...);
	}

	REAL& operator[](auto ...i) {
		return static_cast<DERIVED*>(this)->assign(i...);
	}

	template<int I0 = 0, int I1 = 1>
	Tensor<T, R - 2, D> trace() const {
		Tensor<T, R - 2, D> result;
		static constexpr int J0 = std::min(I0, I1);
		static constexpr int J1 = std::max(I0, I1);
		static constexpr int N1 = pow(D, J0);
		static constexpr int N2 = pow(D, J1 - J0);
		static constexpr int N3 = pow(D, R - J1);
		static constexpr int DK = N3 * (D * N2 + 1);
		for (int n1 = 0; n1 < N1; n1++) {
			for (int n2 = 0; n2 < N2; n2++) {
				for (int n3 = 0; n3 < N3; n3++) {
					int const index1 = n3 + N3 * (n2 + N2 * n1);
					int const index2 = n3 + N3 * D * (n2 + N2 * D * n1);
					auto &res = (Math::Vector<T, pow(D, R - 2)>&) result[index1];
					auto const &src = (Math::Vector<T, S>&) (*this);
					res = REAL(0);
					for (int k = 0; k < D; k++) {
						res += src[index2 + DK * k];
					}
				}
			}
		}
		return result;
	}
};

template<typename T, int R, int D>
class TensorDelta: public TensorBase<TensorDelta<T, R, D>, T, 0, R, D> {
	T const& access(int i, int j) const {
		return (i == j);
	}
};

template<typename T, int R, int D, symmetry_t SYMMETRY>
class Tensor: public TensorBase<Tensor<T, R, D, SYMMETRY>, T, Tensor<T, R, D, SYMMETRY>::size(), R, D> {
	static constexpr int size() {
		if constexpr (SYMMETRY == symmetry_t::asymmetric) {
			int n = 1;
			int d = 1;
			for (int r = 0; r < R; r++) {
				n *= D + r;
				d *= r + 1;
			}
			return n / d;
		} else if constexpr (SYMMETRY == symmetry_t::antisymmetric) {
			int n = 1;
			int d = 1;
			for (int r = 0; r < R; r++) {
				n *= D - r;
				d *= r + 1;
			}
			return n / d;
		} else if constexpr (SYMMETRY == symmetry_t::symmetric) {
			return pow(D, R);
		}
	}
	;
	using base_type = TensorBase<Tensor, T, size(), R, D>;
	base_type::base_type &data;
	Tensor() :
			data(*this) {
	}
	int findIndex(int sum, int i0, auto ...i1) const {
		if constexpr (SYMMETRY == symmetry_t::asymmetric) {
			auto constexpr static map = symmetricMapping<R, D>();
			return findIndex(D * sum + i0, i1...);
		} else if constexpr (SYMMETRY == symmetry_t::symmetric) {
			return findIndex(D * sum + i0, i1...);
		}
	}
	T const& access(auto ...i) const {
		return data[findIndex(0, i...)];
	}
	T& assign(auto ...i) {
		return data[findIndex(0, i...)];
	}
	friend base_type;
};

template<typename TENSOR, int I1 = 0, int I2 = 1>
struct TensorTranspose: public TensorBase<TensorTranspose<TENSOR, I1, I2>, typename TENSOR::REAL, 0, TENSOR::RANK, TENSOR::DIMS> {
	using REAL = typename TENSOR::REAL;
	static constexpr int RANK = TENSOR::RANK;
	static constexpr int DIMS = TENSOR::DIMS;
	TensorTranspose(TENSOR const &t1) :
			T1(t1) {
	}
	REAL const& access(auto ...index) const {
		std::array<int, RANK> pack;
		int i = 0;
		((pack[i++] = index),...);
		std::swap(pack[I1], pack[I2]);
		return std::apply([this](auto ...i) {
			return T1(i...);
		}, pack);
	}
private:
	TENSOR const &T1;

};

template<typename TENSOR1, typename TENSOR2>
struct TensorProduct: public TensorBase<TensorProduct<TENSOR1, TENSOR2>, typename TENSOR1::REAL, 0, (TENSOR1::RANK + TENSOR2::RANK), TENSOR1::DIMS> {
	static constexpr int RANK1 = TENSOR1::RANK;
	static constexpr int RANK2 = TENSOR2::RANK;
	static constexpr int RANK = RANK1 + RANK2;
	static constexpr int DIMS = TENSOR1::DIMS;
	using REAL = typename TENSOR1::REAL;
	TensorProduct(TENSOR1 const &t1, TENSOR2 const &t2) :
			T1(t1), T2(t2) {
	}
	REAL const& access(auto ...index) const {
		std::array<int, RANK1> pack1;
		std::array<int, RANK2> pack2;
		int i = 0;
		((((i < RANK1) ? &pack1[0] : &pack2[-RANK1])[i++] = index),...);
		return std::apply([this](auto ...i) {
			return T1(i...);
		}, pack1) * std::apply([this](auto ...i) {
			return T2(i...);
		}, pack2);
	}
private:
	TENSOR1 const &T1;
	TENSOR2 const &T2;
};

void testEinstein() {
	symmetricMapping<2, 3>();
}
