/*
 * Tensor.hpp
 *
 *  Created on: Mar 2, 2025
 *      Author: dmarce1
 */

#ifndef INCLUDE_TENSOR_HPP_
#define INCLUDE_TENSOR_HPP_
#include <algorithm>
#include <array>
#include <cstdio>
#include <cctype>
#include <functional>

#include "Vector.hpp"

template<int N>
struct Factorial {
	static constexpr int value = N * Factorial<N - 1>::value;
};

template<int X, int N>
struct Pow {
	static constexpr int value = X * Pow<X, N - 1>::value;
};

template<int X>
struct Pow<X, 0> {
	static constexpr int value = 1;
};

template<int X, int N>
struct PowArray {
	static constexpr std::array<int, N + 1> values = []() {
		std::array<int, N + 1> values;
		auto const A = PowArray<X, N - 1>::values;
		values[0] = 1;
		for (int n = 0; n < N; n++) {
			values[n + 1] = X * A[n];
		}
		return values;
	}();
};

template<int X>
struct PowArray<X, 0> {
	static constexpr std::array<int, 1> values { 1 };
};

template<>
struct Factorial<0> {
	static constexpr int value = 1;
};

template<int N, int K>
struct Binomial {
	static constexpr int value = Factorial<N>::value / (Factorial<N - K>::value * Factorial<K>::value);
};

template<int DimCount, int Rank, int ...Symmetries>
struct TensorIndexMap {
	static constexpr auto createMap() {
		static constexpr int MapSize = Pow<DimCount, Rank>::value;
		std::array<int, MapSize> imap;
		std::array<int, Rank> symmetries, I;
		int i = 0;
		((symmetries[i++] = Symmetries),...);
		while (i < Rank) {
			symmetries[i++] = 0;
		}
		bool done = false;
		std::fill(I.begin(), I.end(), 0);
		std::fill(imap.begin(), imap.end(), -1);
		int TensorSize = 0;
		while (!done) {
			int i1 = I[0];
			for (int r = 1; r < Rank; r++) {
				i1 = DimCount * i1 + I[r];
			}
			auto J = I;
			bool zeroFlag = false;
			for (int i = 0; i < Rank; i++) {
				if (symmetries[i]) {
					for (int j = i + 1; j < Rank; j++) {
						if (symmetries[j] == symmetries[i]) {
							if ((symmetries[i] < 0) && J[i] == J[j]) {
								zeroFlag = true;
								break;
							} else if (J[i] < J[j]) {
								std::swap(J[i], J[j]);
							}
						}
					}
					if (zeroFlag) {
						break;
					}
				}
			}
			if (!zeroFlag) {
				int j1 = J[0];
				for (int r = 1; r < Rank; r++) {
					j1 = DimCount * j1 + J[r];
				}
				if (imap[j1] < 0) {
					imap[j1] = TensorSize++;
				}
				imap[i1] = imap[j1];
			}
			int r = Rank - 1;
			while (++I[r] == DimCount) {
				I[r--] = 0;
				if (r < 0) {
					done = true;
					break;
				}
			}
		}
		return imap;
	}
	static constexpr auto value = createMap();
};

template<typename T>
struct TensorBase {
	auto& operator()(auto ...I) {
		return static_cast<T*>(this)->assign(I...);
	}
	auto const& operator()(auto ...I) const {
		return static_cast<T*>(this)->access(I...);
	}
	template<typename C>
	auto& operator[](C I) {
		return static_cast<T*>(this)->assign(I);
	}
	template<typename C>
	auto const& operator[](C I) const {
		return static_cast<T*>(this)->access(I);
	}
};

template<typename T, int D, int R, int ...S>
struct Tensor: public TensorBase<Tensor<T, D, R, S...>> {
	static constexpr int Rank = R;
	static constexpr int Ndims = D;
	using base_type = TensorBase<Tensor>;
	using value_type = T;
	T& assign(auto ...I) {
		int i = 0;
		((i = D * i + I),...);
		return value[iMap[i]];
	}
	template<typename C>
	T& assign(C const &I) {
		int i = 0;
		for (auto j : I) {
			i = D * i + j;
		}
		return value[iMap[i]];
	}
	T const& access(auto ...I) const {
		int i = 0;
		((i = D * i + I),...);
		return value[iMap[i]];
	}
	template<typename C>
	T const& access(C const &I) const {
		int i = 0;
		for (auto j : I) {
			i = D * i + j;
		}
		return value[iMap[i]];
	}
private:
	static constexpr auto iMap = TensorIndexMap<D, R, S...>::value;
	static constexpr int Size = *std::max_element(iMap.begin(), iMap.end()) + 1;
	std::array<T, Size> value;
};

template<typename T, int I0 = 0, int I1 = 1>
struct TensorTranspose: public TensorBase<TensorTranspose<T, I0, I1>> {
	using tensor_type = T;
	using value_type = T::value_type;
	using base_type = TensorBase<TensorTranspose<T, I0, I1>>;
	static constexpr int Rank = T::Rank;
	static constexpr int Ndims = T::Ndims;
	TensorTranspose(tensor_type &A_) :
			A(A_) {
	}
	value_type& assign(auto ...K) {
		std::array<int, sizeof...(K)> J;
		int i = 0;
		((J[i++] = K),...);
		return A.assign(J);
	}
	template<typename C>
	value_type& assign(C K) {
		int i = 0;
		std::swap(K[I0], K[I1]);
		return A.assign(K);
	}
	value_type const& access(auto ...K) const {
		std::array<int, sizeof...(K)> J;
		int i = 0;
		((J[i++] = K),...);
		return A.assign(J);
	}
	template<typename C>
	value_type const& access(C K) const {
		int i = 0;
		std::swap(K[I0], K[I1]);
		return A.assign(K);
	}
private:
	tensor_type &A;
};

template<typename T, int I0 = 0, int I1 = 1>
TensorTranspose<T, I0, I1> tensorTranspose(T &A) {
	return TensorTranspose<T, I0, I1>(A);
}


#endif /* INCLUDE_TENSOR_HPP_ */
