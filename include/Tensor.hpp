#pragma once

#include <array>
#include <bitset>
#include <functional>
#include <utility>
#include <type_traits>

#include "IndexTuple.hpp"
#include "Numbers.hpp"
#include "Permutation.hpp"

namespace Tensors {

template<char C>
struct Index {
	static constexpr char value = C;
	inline constexpr operator char() const {
		return C;
	}
};

template<typename T>
struct IsIndex {
	static constexpr bool value = false;
};

template<char C>
struct IsIndex<Index<C>> {
	static constexpr bool value = true;
	inline operator char() const {
		return value;
	}
};

struct map_index_t;

template<size_t N>
using SymmetryPermutation = std::pair<int, Permutation<N>>;

template<typename T, typename F, size_t D, size_t R, size_t N, std::array<SymmetryPermutation<R>, N> S, char ...I>
struct TensorExpression;

template<typename T, size_t D, size_t R, SymmetryPermutation<R> ...S>
struct Tensor;

struct map_index_t {
	signed char sign;
	unsigned index;
};

template<typename U, size_t ... I>
inline consteval auto getTupleTail(const U &tup, std::index_sequence<I...>) {
	return std::make_tuple(std::get<I + 1>(tup)...);
}

template<typename ... Args>
inline consteval auto tupleTail(const std::tuple<Args...> &tup) {
	return getTupleTail(tup, std::make_index_sequence<sizeof...(Args) - 1> {});
}

template<typename ... Args>
inline consteval std::tuple<> createIndexTuple(std::tuple<Args...> tup) {
	return std::tuple<Args...>();
}

template<char C, typename ...Args>
inline consteval auto createIndexTuple(std::tuple<Index<C>, Args...> tup, size_t i, auto ...J) {
	std::tuple<size_t> const thisTuple(i);
	auto const nextTuple = createIndexTuple(tupleTail(tup), J...);
	return std::tuple_cat(thisTuple, nextTuple);
}

template<typename ...Args>
inline consteval auto createIndexTuple(std::tuple<int, Args...> tup, auto ...is) {
	std::tuple<size_t> const thisTuple(std::get < 0 > (tup));
	auto const nextTuple = createIndexTuple(tupleTail(tup), is...);
	return std::tuple_cat(thisTuple, nextTuple);
}

template<typename ...Args>
inline consteval std::tuple<> extractChars(Args&&...) {
	return std::tuple<>();
}

template<typename ...Args>
inline consteval auto extractChars(int i, Args &&...args) {
	return extractChars(std::forward<Args>(args)...);
}

template<char C, typename ...Args>
inline consteval auto extractChars(Index<C>, Args &&...args) {
	return std::tuple_cat(std::tuple<Index<C>>(), extractChars(std::forward<Args>(args)...));
}

template<typename T, size_t R, Permutation<R> P>
inline consteval auto permuteList(T const &list) {
	auto permutedList = list;
	for (size_t i = 0; i < list.size(); i++) {
		for (size_t j = 0; j < R; j++) {
			permutedList[i][j] = list[i][P[j] - 1];
		}
	}
	return permutedList;
}

template<size_t N, std::array<char, N> A, std::array<char, N> B>
inline consteval Permutation<N> computePermutation() {
	Permutation<N> P;
	for (size_t i = 0; i < N; i++) {
		size_t j = 0;
		while (A[i] != B[j]) {
			j++;
		}
		P[i] = j + 1;
	}
	return P;
}

template<size_t R, size_t N, std::array<SymmetryPermutation<R>, N> P, typename ...Args>
inline consteval size_t prunePermutationsCount(std::tuple<Args...>) {
	static constexpr size_t M = (std::is_integral<Args>::value + ...);
	std::array < size_t, M > indices;
	size_t i = 0;
	size_t k = 0;
	((std::is_integral<Args>::value ? (indices[i++] = k++) : k++),...);
	std::sort(indices.begin(), indices.end());
	size_t count = 0;
	for (i = 0; i < N; i++) {
		size_t flag = 1;
		for (k = 0; k < M; k++) {
			size_t const j = indices[k];
			if (P[i].second[j] != j + 1) {
				flag = 0;
				break;
			}
		}
		count += flag;
	}
	return count;
}

template<size_t R, size_t N, std::array<SymmetryPermutation<R>, N> P, typename ...Args>
inline consteval auto prunePermutations(std::tuple<Args...>) {
	static constexpr size_t M = (std::is_integral<Args>::value + ...);
	static constexpr size_t L = prunePermutationsCount<N, P, Args...>(std::tuple<Args...>());
	std::array<int, M + 2> indices;
	size_t i = 1;
	size_t k = 0;
	((std::is_integral<Args>::value ? (indices[i++] = k++) : k++),...);
	indices[0] = -1;
	indices[M + 1] = R;
	std::sort(indices.begin(), indices.end());
	std::array<SymmetryPermutation<R - M>, L> Q;
	size_t nextIndex = 0;
	for (size_t n = 0; n < N; n++) {
		bool flag = true;
		for (k = 1; k <= M; k++) {
			size_t const j = indices[k];
			if (P[n].second[j] != j + 1) {
				flag = false;
				break;
			}
		}
		if (flag) {
			Permutation<R - M> thisPerm;
			k = 0;
			for (size_t m = 0; m <= M; m++) {
				for (int j = indices[m] + 1; j < indices[m + 1]; j++) {
					thisPerm[k++] = P[n].second[j];
				}
			}
			for (size_t m = 0; m < R - M; m++) {
				bool flag;
				do {
					flag = true;
					for (size_t n = 0; n < R - M; n++) {
						if (thisPerm[n] == m + 1) {
							flag = false;
							break;
						}
					}
					if (flag) {
						for (size_t n = 0; n < R - M; n++) {
							if (thisPerm[n] > m + 1) {
								thisPerm[n]--;
							}
						}
					}
				}while (flag);
			}
			SymmetryPermutation<R - M> q {P[n].first, thisPerm};
			Q[nextIndex].first = q.first;
			Q[nextIndex].second = q.second;
			nextIndex++;
		}
	}

	return Q;
}

template<size_t D, size_t R, size_t N, std::array<SymmetryPermutation<R>, N> S>
inline consteval size_t computeUniqueIndicesCount() {
	using i_type = IndexTuple<D, R>;
	static constexpr size_t Size = Math::integerPower(D, R);
	std::bitset < Size > haveVisited;
	size_t count = 0;
	for (i_type tensorIndices = i_type::begin(); tensorIndices != i_type::end(); tensorIndices++) {
		size_t const tensorIndex = tensorIndices.flatIndex();
		if (!haveVisited[tensorIndex]) {
			count++;
			haveVisited[tensorIndex] = true;
			for (size_t n = 0; n < N; n++) {
				auto const permutation = S[n].second;
				auto permutedIndices = permutation.apply(tensorIndices);
				do {
					haveVisited[permutedIndices.flatIndex()] = true;
					permutedIndices = permutation.apply(permutedIndices);
				}while (permutedIndices != tensorIndices);
			}
		}
	}
	return count;
}

template<size_t D, size_t R, size_t N, std::array<SymmetryPermutation<R>, N> S>
inline consteval auto computeUniqueIndices() {
	using i_type = IndexTuple<D, R>;
	static constexpr size_t Size = Math::integerPower(D, R);
	static constexpr size_t Count = computeUniqueIndicesCount<D, R, N, S>();
	std::array<i_type, Count> uniqeIndices;
	std::bitset < Size > haveVisited;
	size_t count = 0;
	for (i_type tensorIndices = i_type::begin(); tensorIndices != i_type::end(); tensorIndices++) {
		size_t const tensorIndex = tensorIndices.flatIndex();
		if (!haveVisited[tensorIndex]) {
			haveVisited[tensorIndex] = true;
			uniqeIndices[count] = tensorIndices;
			count++;
			for (size_t n = 0; n < N; n++) {
				auto const permutation = S[n].second;
				auto permutedIndices = permutation.apply(tensorIndices);
				do {
					haveVisited[permutedIndices.flatIndex()] = true;
					permutedIndices = permutation.apply(permutedIndices);
				}while (permutedIndices != tensorIndices);
			}
		}
	}
	return uniqeIndices;
}

template<typename T, size_t D, size_t R, SymmetryPermutation<R> ...S>
inline consteval auto createIndexMap() {
	using i_type = IndexTuple<D, R>;
	static constexpr size_t Size = Math::integerPower(D, R);
	static constexpr size_t N = sizeof...(S);
	static constexpr std::array<SymmetryPermutation<R>, N> permutations = {S...};
	std::pair<std::array<T, Size>, std::array<size_t, Size> > rc;
	auto &signMap = rc.first;
	auto &indexMap = rc.second;
	std::bitset < Size > haveVisited;
	size_t elementCount = 0;
	for (i_type tensorIndices = i_type::begin(); tensorIndices != i_type::end(); tensorIndices++) {
		size_t const tensorIndex = tensorIndices.flatIndex();
		if (!haveVisited[tensorIndex]) {
			haveVisited[tensorIndex] = true;
			unsigned const commonMapIndex = elementCount++;
			indexMap[tensorIndex] = commonMapIndex;
			signMap[tensorIndex] = T(1);
			for (size_t n = 0; n < N; n++) {
				auto const &symmetrySign = permutations[n].first;
				auto const &permutation = permutations[n].second;
				auto permutedIndices = permutation.apply(tensorIndices);
				do {
					size_t const mapIndex = permutedIndices.flatIndex();
					haveVisited[mapIndex] = true;
					signed char const sign = ((symmetrySign < 0) && (permutation.parity() < 0)) ? ((permutedIndices == tensorIndices) ? 0 : -1) : +1;
					indexMap[mapIndex] = commonMapIndex;
					signMap[mapIndex] = T(sign);
					permutedIndices = permutation.apply(permutedIndices);
				}while (permutedIndices != tensorIndices);
			}
		}
	}
	return rc;
}

template<size_t L, typename T, size_t N, std::array<std::pair<size_t, size_t>, N> I, size_t J = 0, size_t K = 0>
inline consteval auto pruneCharTuple() {
	if constexpr (K == std::tuple_size < T > ::value) {
		return std::tuple<> {};
	} else if constexpr ((L == 0) && (I[std::min(J, N - 1)].first == K)) {
		return pruneCharTuple<L, T, N, I, J + 1, K + 1>();
	} else if constexpr ((L == 1) && (I[std::min(J, N - 1)].second == K)) {
		return pruneCharTuple<L, T, N, I, J + 1, K + 1>();
	} else {
		static constexpr T tup {};
		static constexpr char C = std::get < K > (tup).value;
		static constexpr Index<C> index {};
		return std::tuple_cat(std::tuple < Index < C >> (index), pruneCharTuple<L, T, N, I, J, K + 1>());
	}
}

template<size_t D, char ...I>
struct CharPackTraits {

	template<char ...J>
	inline static consteval size_t contractionCount() {
		static constexpr std::array<char, sizeof...(J)> jArray = {J...};
		size_t count = 0;
		for (size_t i = 0; i < sizeof...(I); i++) {
			for (size_t j = 0; j < sizeof...(J); j++) {
				count += size_t(iArray[i] == jArray[j]);
			}
		}
		return count;
	}

	template<char ...J>
	inline static consteval auto contractionPairs() {
		static constexpr std::array<char, sizeof...(J)> jArray = {J...};
		static constexpr size_t M = contractionCount<J...>();
		std::array<std::pair<size_t, size_t>, M> pairs;
		size_t count = 0;
		for (size_t i = 0; i < sizeof...(I); i++) {
			for (size_t j = 0; j < sizeof...(J); j++) {
				if (iArray[i] == jArray[j]) {
					pairs[count++] = std::pair<size_t, size_t>(i, j);
				}
			}
		}
		return pairs;
	}

	template <size_t R3>
	inline static consteval auto normalizePermutation(const Permutation<R3>& input) {
		Permutation<R3> sorted = input;
		std::sort(sorted.begin(), sorted.end());
		Permutation<R3> result {};
		for (std::size_t i = 0; i < R3; ++i) {
			auto it = std::find(sorted.begin(), sorted.end(), input[i]);
			result[i] = size_t(std::distance(sorted.begin(), it)) + 1;
		}
		return result;
	}

	template<size_t R1, size_t N1, std::array<SymmetryPermutation<R1>, N1> P1, size_t R2, size_t N2, std::array<SymmetryPermutation<R2>, N2> P2, char ...J>
	inline static consteval auto concatenatePermutationsHelper() {
		static constexpr size_t M = contractionCount<J...>();
		static constexpr size_t R3 = R1 + R2 - 2 * M;
		static constexpr size_t N3 = N1 * N2;
		static constexpr auto pairs = contractionPairs<J...>();
		std::array<SymmetryPermutation<R3>, N3> P3;
		std::array<size_t, R1> shift1 = {0};
		std::array<size_t, R2> shift2 = {0};
		for(auto const& pair : pairs) {
			for(size_t r1 = pair.first; r1 < R1; r1++) {
				shift1[r1]++;
			}
			for(size_t r2 = pair.second; r2 < R2; r2++) {
				shift2[r2]++;
			}
		}
		for( size_t n1 = 0; n1 < N1; n1++) {
			for( size_t n2 = 0; n2 < N2; n2++) {
				SymmetryPermutation<R3> p3;
				for(size_t r1 = 0; r1 < R1 - M; r1++) {
					p3.second[r1] = P1[n1].second[r1 + shift1[r1]];
				}
				for(size_t r2 = 0; r2 < R2 - M; r2++) {
					p3.second[R1 - M + r2] = P2[n2].second[r2 + shift2[r2]] + R1;
				}
				p3.second = normalizePermutation(p3.second);
				p3.first = P1[n1].first * P2[n2].first;
				P3[n1 * N2 + n2] = p3;
			}
		}
		std::sort(P3.begin(), P3.end());
		size_t const count = std::unique(P3.begin(), P3.end()) - P3.begin();
		return std::pair<size_t, std::array<SymmetryPermutation<R3>, N3>>(count, P3);
	}

	template<size_t R1, size_t N1, std::array<SymmetryPermutation<R1>, N1> P1, size_t R2, size_t N2, std::array<SymmetryPermutation<R2>, N2> P2, char ...J>
	inline static consteval auto concatenatePermutations() {
		static constexpr auto P3 = concatenatePermutationsHelper<R1, N1, P1, R2, N2, P2, J...>();
		static constexpr size_t N3 = P3.first;
		static constexpr size_t M = contractionCount<J...>();
		static constexpr size_t R3 = R1 + R2 - 2 * M;
		std::array<SymmetryPermutation<R3>, N3> P4;
		std::copy(P3.second.begin(), P3.second.begin() + N3, P4.begin());
		return P4;
	}

	template<char ...J>
	inline static consteval auto contractionGuide() {
		using namespace Math;
		static constexpr auto conPairs = contractionPairs<J...>();
		static constexpr size_t M = conPairs.size();
		static constexpr size_t D3M = integerPower(D, M);
		static constexpr size_t R1a = sizeof...(I);
		static constexpr size_t R2a = sizeof...(J);
		static constexpr size_t R1b = R1a - M;
		static constexpr size_t R2b = R2a - M;
		static constexpr size_t R3 = R1b + R2b;
		static constexpr size_t N3 = integerPower(D, R3);
		using i1a_type = IndexTuple<D, R1a>;
		using i2a_type = IndexTuple<D, R2a>;
		using i1b_type = IndexTuple<D, R1b>;
		using i2b_type = IndexTuple<D, R2b>;
		using i3_type = IndexTuple<D, R3>;
		using i4_type = IndexTuple<D, M>;
		using pair_type = std::pair<i1a_type, i2a_type>;
		using factors_type = std::array<pair_type, D3M>;
		std::array<factors_type, N3> U;
		for (i3_type I3 = i3_type::begin(); I3 != i3_type::end(); I3++) {
			i1a_type J1;
			i2a_type J2;
			i1b_type I1;
			i2b_type I2;
			for (size_t r1b = 0; r1b < R1b; r1b++) {
				I1[r1b] = I3[r1b];
			}
			for (size_t r2b = 0; r2b < R2b; r2b++) {
				I2[r2b] = I3[R1b + r2b];
			}
			size_t n = 0;
			for (size_t m = 0; m <= M; m++) {
				auto const nEnd = (m == M) ? R1b : conPairs[m].first;
				for (; n < nEnd; n++) {
					J1[n + m] = I1[n];
				}
			}
			n = 0;
			for (size_t m = 0; m <= M; m++) {
				auto const nEnd = (m == M) ? R2b : conPairs[m].second;
				for (; n < nEnd; n++) {
					J2[n + m] = I2[n];
				}
			}
			auto &u3 = U[I3.flatIndex()];
			for (i4_type I4 = i4_type::begin(); I4 != i4_type::end(); I4++) {
				for (size_t m = 0; m < M; m++) {
					auto const p = conPairs[m];
					J1[p.first] = J2[p.second] = I4[m];
				}
				u3[I4.flatIndex()] = pair_type {J1, J2};
			}
		}
		return U;
	}

private:
	static constexpr std::array<char, sizeof...(I)> iArray = { I... };
};

template<typename T, typename F, size_t D, size_t R, size_t N, std::array<SymmetryPermutation<R>, N> S, char ...I>
struct TensorExpression {
	using index_type = IndexTuple<D, R>;
	using value_type = T;
	template<typename F1, size_t N1, std::array<SymmetryPermutation<R>, N1> S1, char ...I1>
	using assignment_type = std::function<T&(T &, T )>;
	static constexpr auto assignment = [](T &x, T y) {
		return (x = y);
	};
	static constexpr auto plus = [](T &x, T y) {
		return (x += y);
	};
	static constexpr auto minus = [](T &x, T y) {
		return (x -= y);
	};
	inline constexpr TensorExpression(F const &h) :
			handle(h) {
	}
	inline constexpr auto operator()(auto ...is) const {
		return handle(is...);
	}
	inline constexpr auto& operator()(auto ...is) {
		return handle(is...);
	}
	inline constexpr auto operator()(IndexTuple<D, R> const &iii) const {
		return std::apply([this](auto ...is) {
			return this->operator()(is...);
		}, std::array<size_t, R>(iii));
	}
	inline constexpr auto& operator()(IndexTuple<D, R> const &iii) {
		return std::apply([this](auto ...is) {
			return this->operator()(is...);
		}, std::array<size_t, R>(iii));
	}

	template<typename F2, size_t R2, size_t N2, std::array<SymmetryPermutation<R2>, N2> S2, char ...C>
	inline constexpr auto createTensorExpression(F2 const &handle, std::tuple<Index<C> ...>) const {
		return TensorExpression<T, F2, D, R2, N2, S2, C...>(handle);
	}

	template<typename F1, size_t R1, size_t N1, std::array<SymmetryPermutation<R1>, N1> S1, char ...J>
	inline constexpr auto operator*(TensorExpression<T, F1, D, R1, N1, S1, J...> other) const {
		static constexpr size_t R2 = R + R1 - 2 * charPack.template contractionCount<J...>();
		auto const thisHandle = [this, other](auto ...is) {
			static constexpr auto U = charPack.template contractionGuide<J...>();
			IndexTuple<D, R2> const indices(std::array<size_t, R2> { is... });
			auto const &u = U[indices.flatIndex()];
			T sum = handle(u[0].first) * other(u[0].second);
			for (size_t m = 1; m < u.size(); m++) {
				sum += handle(u[m].first) * other(u[m].second);
			}
			return sum;
		};
		static constexpr auto pairs = charPack.template contractionPairs<J...>();
		static constexpr size_t N2 = pairs.size();
		static constexpr auto tup1 = pruneCharTuple<0, std::tuple<Index<I> ...>, N2, pairs>();
		static constexpr auto tup2 = pruneCharTuple<1, std::tuple<Index<J> ...>, N2, pairs>();
		static constexpr auto charTuple = std::tuple_cat(tup1, tup2);
		static constexpr auto S2 = charPack.template concatenatePermutations<R, N, S, R1, N1, S1, J...>();
		//		static constexpr std::array<SymmetryPermutation<R2>, 1> S2 = { SymmetryPermutation<R2>( { 1, Permutation<R2>::identity }) };
		return createTensorExpression<decltype(thisHandle), R2, S2.size(), S2>(thisHandle, charTuple);
	}
	template<typename AS, typename F1, size_t N1, std::array<SymmetryPermutation<R>, N1> S1, char ...I1>
	inline constexpr auto& assignOperation(AS assign, TensorExpression<T, F1, D, R, N1, S1, I1...> const &other) {
		static constexpr std::array<char, R> A = { I... };
		static constexpr std::array<char, R> B = { I1... };
		static constexpr auto P = computePermutation<R, A, B>();
		static constexpr auto indicesList = computeUniqueIndices<D, R, N1, S1>();
		static constexpr auto permutedIndicesList = permuteList<decltype(indicesList), R, P>(indicesList);
		for (size_t i = 0; i < indicesList.size(); i++) {
			T &ref = handle(indicesList[i]);
			assign(ref, other(permutedIndicesList[i]));
		}
		return *this;
	}
	template<typename F1, size_t N1, std::array<SymmetryPermutation<R>, N1> S1, char ...I1>
	inline constexpr auto& operator=(TensorExpression<T, F1, D, R, N1, S1, I1...> const &other) {
		return assignOperation(assignment, other);
	}
	template<typename F1, size_t N1, std::array<SymmetryPermutation<R>, N1> S1, char ...I1>
	inline constexpr auto& operator+=(TensorExpression<T, F1, D, R, N1, S1, I1...> const &other) {
		return assignOperation(plus, other);
	}
	template<typename F1, size_t N1, std::array<SymmetryPermutation<R>, N1> S1, char ...I1>
	inline constexpr auto& operator-=(TensorExpression<T, F1, D, R, N1, S1, I1...> const &other) {
		return assignOperation(minus, other);
	}
private:
	static constexpr CharPackTraits<D, I...> charPack { };
	F handle;
};

template<typename T, size_t D, size_t R, bool hasAntisymmetry, SymmetryPermutation<R> ...S>
struct TensorAccess {
	inline constexpr operator T() const {
		return sign * A.data[index];
	}
	inline constexpr auto& operator=(T const &other) const {
		A.data[index] = sign * other;
		return A;
	}
	inline constexpr TensorAccess(Tensor<T, D, R, S...> &a, T const &b, size_t const &c) :
			A(a), sign(b), index(c) {
	}
private:
	Tensor<T, D, R, S...> &A;
	T const &sign;
	size_t const &index;
};

template<typename T, size_t D, size_t R, SymmetryPermutation<R> ...S>
struct TensorAccess<T, D, R, false, S...> {
	inline constexpr operator T() const {
		return A.data[index];
	}
	inline constexpr operator T&() {
		return A.data[index];
	}
	inline constexpr auto& operator=(T const &other) const {
		A.data[index] = other;
		return A;
	}
	inline constexpr TensorAccess(Tensor<T, D, R, S...> &a, size_t const &c) :
			A(a), index(c) {
	}
private:
	Tensor<T, D, R, S...> &A;
	size_t const &index;
};

template<typename F, typename T, size_t D, size_t R, size_t N, std::array<SymmetryPermutation<R>, N> S, char ...C>
inline constexpr auto createTensorExpression(F const &handle, std::tuple<Index<C> ...>) {
	return TensorExpression<F, T, D, R, N, S, C...>(handle);
}

template<typename T, size_t D, size_t R, SymmetryPermutation<R> ...S>
struct Tensor {
	using value_type = T;
	using index_type = IndexTuple<D, R>;

	inline constexpr auto operator()(IndexTuple<D, R> const &iii) {
		return std::apply([this](auto ...is) {
			return this->operator()(is...);
		}, std::array<size_t, R>(iii));
	}
	inline constexpr auto operator()(IndexTuple<D, R> const &iii) const {
		return const_cast<Tensor&>(*this).operator()(iii);
	}
	inline constexpr auto operator()(auto ...is) {
		if constexpr (hasAntisymmetry) {
			size_t const flatIndex = index_type( { size_t(is)... }).flatIndex();
			return TensorAccess<T, D, R, true, S...>(*this, signMap[flatIndex], indexMap[flatIndex]);
		} else {
			return TensorAccess<T, D, R, false, S...>(*this, indexMap[index_type( { size_t(is)... }).flatIndex()]);
		}
	}
	inline constexpr T operator()(auto ...is) const {
		if constexpr (hasAntisymmetry) {
			size_t const flatIndex = index_type( { is... }).flatIndex();
			return signMap[flatIndex] * data[indexMap[flatIndex]];
		} else {
			return data[indexMap[index_type( { is... }).flatIndex()]];
		}
	}

	template<char ...I>
	inline constexpr auto operator()(Index<I> ...) const {
		return const_char<Tensor const&>(*this).operator()(Index<I>() ...);
	}

	template<char ...I>
	inline constexpr auto operator()(Index<I> ...) {
		auto const handle = [this](auto ...is) {
			return this->operator()(is...);
		};
		return TensorExpression<T, decltype(handle), D, R, sizeof...(S), Symmetries, I...>(handle);
	}

	template<typename ...Args, std::enable_if<((std::is_integral<Args>::value || ...) && (IsIndex<Args>::value || ...)), int>::type = 0>
	inline constexpr auto operator()(Args &&...args) {
		std::tuple<Args...> const argTuple(std::forward<Args>(args)...);
		auto const handle = [this, argTuple](auto ...is) {
			auto const intTuple = createIndexTuple(argTuple, is...);
			return std::apply([this](auto ...is) {
				return this->operator()(is...);
			}, intTuple);
		};
		static constexpr auto prunedSymmetries = prunePermutations<sizeof...(S), Symmetries, Args...>(std::tuple<Args...>());
		auto const chars = extractChars(args...);
		return createTensorExpression<decltype(handle), std::tuple_size<decltype(chars)>::value, prunedSymmetries.size(), prunedSymmetries>(handle, chars);
	}

	template<typename ...Args, std::enable_if<((std::is_integral<Args>::value || ...) && (IsIndex<Args>::value || ...)), int>::type = 0>
	inline constexpr auto operator()(Args &&...args) const {
		return const_cast<Tensor&>(*this)(std::forward<Args>(args)...);
	}

	inline constexpr std::string toString() const {
		std::string str;
		for (index_type I = index_type::begin(); I != index_type::end(); I++) {
			str += "(";
			for (size_t r = 0; r < R; r++) {
				str += std::to_string(I[r]);
				if (r + 1 < R) {
					str += ", ";
				}
			}
			str += ") ";
			str += std::to_string((*this)(I));
			str += "\n";
		}
		str += "\n";
		return str;
	}

private:
	static constexpr std::array<SymmetryPermutation<R>, sizeof...(S)> Symmetries = { S... };
	static constexpr size_t Size = Math::integerPower(D, R);
	static constexpr auto signIndexMap = createIndexMap<T, D, R, S...>();
	static constexpr auto signMap = signIndexMap.first;
	static constexpr auto indexMap = signIndexMap.second;
	static constexpr bool hasAntisymmetry = bool(T(0) == std::accumulate(signMap.begin(), signMap.end(), T(1), [](T const &a, T const &b) {
		return a * b;
	}));
	friend class TensorAccess<T, D, R, hasAntisymmetry, S...> ;
	std::array<T, Size> data;
};

}

