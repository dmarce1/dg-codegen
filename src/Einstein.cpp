/*
 * Einstein.cpp
 *
 *  Created on: Feb 25, 2025
 *      Author: dmarce1
 */

#include <algorithm>
#include <array>
#include <numeric>
#include <tuple>
#include <stdio.h>
#include <tuple>
#include <string>
#include <utility>

#include <utility>
#include <type_traits>

template<int I1, int ...I2>
struct Symmetry {
	constexpr auto operator()(std::integer_sequence<int, I1, I2...> = std::integer_sequence<int, I1, I2...>()) {
		static constexpr int size = 1 + sizeof...(I2);
		std::array<int, size> indices;
		indices[0] = I1;
		if constexpr (size > 1) {
			Symmetry<I2...> nextFunctor;
			auto const theRest = nextFunctor();
			static_assert(theRest.size() == size - 1);
			for (size_t i = 0; i < theRest.size(); i++) {
				indices[i + 1] = theRest[i];
			}
		}
		return indices;
	}
	template<int J, int K>
	using contracted_type = Symmetry< I1 - int(I1>J) - int(I1>K), (I2 - int(I2>J) - int(I2>K))... >;
};

template<int I, typename Seq>
struct contains;

template<int I, int Head, int ... Tail>
struct contains<I, Symmetry<Head, Tail...> > {
	static constexpr bool value = ((I == Head) ? true : contains<I, Symmetry<Tail...>>::value);
};

template<int I>
struct contains<I, void> {
	static constexpr bool value = false;
};

template<int I, int ... Is>
struct symmetryContract {
	using type = Symmetry< ((Is > I) ? (Is - 1) : Is)... >;
};

template<int RANK, int DIM>
static constexpr int indices2Index(std::array<int, RANK> indices) {
	return std::accumulate(indices.rbegin(), indices.rend(), 0, [](int i, int j) {
		return DIM * i + j;
	});
}

template<int RANK, int DIM>
static constexpr std::array<int, RANK> index2Indices(int index) {
	std::array<int, RANK> indices;
	for (int i = 0; i < RANK; i++) {
		indices[i] = index % DIM;
		index /= DIM;
	}
	return indices;
}

constexpr int pow(int x, int n) {
	int xm = x;
	int xn = 1;
	while (n > 0) {
		if (n & 1) {
			xn *= xm;
		}
		xm *= xm;
		n >>= 1;
	}
	return xn;
}

template<typename ...ARGS>
struct FirstArgument;

template<>
struct FirstArgument<> {
	using type = void;
};

template<typename ... ARGS>
constexpr std::array<int, sizeof...(ARGS)> ints2Array(ARGS ... args) {
	static_assert((std::is_integral_v<ARGS> && ...), "All arguments must be integers.");
	return {static_cast<int>(args)...};
}

template<typename FIRST, typename ...ARGS>
struct FirstArgument<FIRST, ARGS...> {
	using type = FIRST;
};

template<int RANK, typename FIRST = void, typename ...ARGS>
static constexpr auto sortIndices(std::array<int, RANK> const &indices) {
	if constexpr (std::is_void<FIRST>::value) {
		return indices;
	} else {
		FIRST firstSymmetry;
		auto const symmetrySet = firstSymmetry();
		auto sortedSymmetries = symmetrySet;
		auto sortedIndices = indices;
		std::sort(sortedSymmetries.begin(), sortedSymmetries.end(), [indices](int a, int b) {
			return indices[a] > indices[b];
		});
		for (size_t i = 0; i < sortedSymmetries.size(); i++) {
			sortedIndices[symmetrySet[i]] = indices[sortedSymmetries[i]];
		}
		if constexpr (sizeof...(ARGS)) {
			return sortIndices<RANK, ARGS...>(sortedIndices);
		} else {
			return sortedIndices;
		}
	}
}

template<int RANK, int DIM, typename ...ARGS>
static constexpr auto createIndexMap() {
	static constexpr int size = pow(DIM, RANK);
	std::array<int, size> map;
	std::array<bool, size> filled;
	std::fill(filled.begin(), filled.end(), false);
	int nextIndex = 0;
	for (int thisIndex = 0; thisIndex < size; thisIndex++) {
		auto const theseIndices = index2Indices<RANK, DIM>(thisIndex);
		auto const commonIndices = sortIndices<RANK, ARGS...>(theseIndices);
		int const commonIndex = indices2Index<RANK, DIM>(commonIndices);
		if (!filled[commonIndex]) {
			map[commonIndex] = nextIndex++;
			filled[thisIndex] = true;
		}
		map[thisIndex] = map[commonIndex];
	}
	return map;
}

template<int I1, int I2, typename Head, typename ... Tail>
struct SymmetryContraction {
	using head_type = std::conditional<
	contains<I1, Head, Tail...>::value || contains<I2, Head, Tail...>::value,
	std::tuple<>,
	std::tuple<typename Head::contracted_type<I1, I2>>
	>::type;
	using tail_type = SymmetryContraction<I1, I2, Tail...>::type;
	using type = std::remove_reference<decltype(std::tuple_cat(head_type(), tail_type()))>::type;
};

template<typename DERIVED, int RANK, int DIM, typename ...SYMMETRIES>
struct TensorBase {

	template<typename ...ARGS>
	constexpr auto const& operator[](ARGS &&... args) const {
		return derived.access_impl(std::forward<ARGS>(args)...);
	}

	template<typename ...ARGS>
	constexpr auto& operator[](ARGS &&... args) {
		return derived.assign_impl(std::forward<ARGS>(args)...);
	}

private:
	DERIVED &derived = *static_cast<DERIVED*>(this);
};

template<typename, int, int, typename ...>
struct contracted_type;

//template<int I1, int I2, int RANK, int DIM>
//static constexpr auto contractionList() {
//	static constexpr int SIZE = pow(DIM, RANK - 2);
//	std::array<std::array<int, DIM>, SIZE> lists;
//	int next_index = 0;
//	for (int k = 0; k < RANK * SIZE; k++) {
//		auto const theseIndices = index2Indices<RANK, DIM>(k);
//		if (theseIndices[I1] = theseIndices[I2]) {
//			for (int l = 0; l < DIM; l++) {
//				lists[next_index++][theseIndices[I1]] = k;
//			}
//		}
//	}
//	return lists;
//}

template<typename TYPE, int RANK, int DIM, typename ...SYMMETRIES>
struct Tensor: public TensorBase<Tensor<TYPE, RANK, DIM, SYMMETRIES...>, RANK, DIM, SYMMETRIES...> {
	using value_type = TYPE;
	static constexpr int rank = RANK;
	static constexpr int ndim = DIM;

	template<int I1 = 0, int I2 = 1>
	constexpr auto contract() const {
		using sym_tuple = SymmetryContraction<I1, I2, SYMMETRIES...>::type;
		using return_type = contracted_type<TYPE, RANK, DIM, sym_tuple>::type;
		return_type result;

		return result;
	}

	constexpr auto const& access_impl(std::array<int, RANK> const &indices) const {
		return values[indexMap[indices2Index<rank, ndim>(indices)]];
	}

	constexpr auto& access_impl(std::array<int, RANK> const &indices) {
		return values[indexMap[indices2Index<rank, ndim>(indices)]];
	}

	template<typename ...ARGS>
	constexpr auto const& access_impl(int I, ARGS &&... args) const {
		return access_impl(ints2Array(I, std::forward<ARGS>(args)...));
	}

	template<typename ...ARGS>
	auto& assign_impl(int I, ARGS &&... args) {
		return access_impl(ints2Array(I, std::forward<ARGS>(args)...));
	}

private:
	static constexpr auto indexMap = createIndexMap<rank, ndim, SYMMETRIES...>();
	static constexpr int size = *std::max_element(indexMap.begin(), indexMap.end()) + 1;
	std::array<TYPE, size> values;
};

template<typename TYPE, int RANK, int DIM>
struct Tensor<TYPE, RANK, DIM> : public TensorBase<Tensor<TYPE, RANK, DIM>, RANK, DIM> {
	using value_type = TYPE;
	static constexpr int rank = RANK;
	static constexpr int ndim = DIM;

	template<int I1 = 0, int I2 = 1>
	constexpr auto contract() const {
		if constexpr (I1 > I2) {
			return contract<I2, I1>();
		} else {
			Tensor<TYPE, RANK - 2, DIM> result;
			std::array<int, RANK> Is;
			std::fill(Is.begin(), Is.end(), 0);
			constexpr int Imin = std::min(I1, I2);
			constexpr int Imax = std::max(I1, I2);
			constexpr int N1 = pow(DIM, Imin);
			constexpr int N2 = pow(DIM, Imax - Imin);
			constexpr int N3 = pow(DIM, RANK - Imax);
			for (int n3 = 0; n3 < N3; n3++) {
				for (int n2 = 0; n2 < N2; n2++) {
					for (int n1 = 0; n1 < N1; n1++) {
						TYPE sum = TYPE(0);
						int const index1 = N1 * (N2 * n3 + n2) + n1;
						int const index2 = N1 * DIM * (N2 * DIM * n3 + n2) + n1;
						int const index3 = N1 * (DIM * N2 + 1);
						for (int k = 0; k < DIM; k++) {
						}
					}
				}
			}
			return result;
		}
	}

	constexpr auto const& access_impl(std::array<int, RANK> const &indices) const {
		return values[indices2Index<rank, ndim>(indices)];
	}

	constexpr auto access_impl(std::array<int, RANK> const &indices) {
		return values[indices2Index<rank, ndim>(indices)];
	}

	template<typename ...ARGS>
	constexpr auto const& access_impl(int I, ARGS &&... args) const {
		return access_impl(ints2Array(I, std::forward<ARGS>(args)...));
	}

	template<typename ...ARGS>
	constexpr auto& assign_impl(int I, ARGS &&... args) {
		return access_impl(ints2Array(I, std::forward<ARGS>(args)...));
	}

private:
	static constexpr int size = pow(ndim, rank);
	std::array<TYPE, size> values;
};

template<typename TYPE, int RANK, int DIM, typename ...SYMMETRIES>
struct contracted_type<TYPE, RANK, DIM, std::tuple<SYMMETRIES...>> {
	using type = Tensor<TYPE, RANK - 2, DIM, SYMMETRIES...>;
};

void testEinstein() {
	constexpr int RANK = 4;
	constexpr int DIM = 3;
	Tensor<double, RANK, DIM, Symmetry<1, 2>> Ta;
	Ta[2, 3, 4, 5];
}
