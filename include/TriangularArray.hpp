/*
 * TriangularArray.hpp
 *
 *  Created on: Jan 14, 2025
 *      Author: dmarce1
 */

#ifndef INCLUDE_TRIANGULARARRAY_HPP_
#define INCLUDE_TRIANGULARARRAY_HPP_

#include <unordered_map>

#include "Vector.hpp"

template<int, int>
struct TriangularIndices;

template<int DimCnt, int MomCnt>
struct TriangularIndices {
	using indices_type = Math::Vector<int, DimCnt>;
	using map_type = std::unordered_map<indices_type, int, Math::vectorHashKey<int, DimCnt, MomCnt>>;
	TriangularIndices() {
		reset();
	}
	TriangularIndices(TriangularIndices const&) = default;
	TriangularIndices(TriangularIndices&&) = default;
	TriangularIndices(indices_type const &indices) {
		set(indices);
	}
	TriangularIndices& operator=(TriangularIndices const&) = default;
	TriangularIndices& operator=(TriangularIndices&&) = default;
	TriangularIndices& operator++() {
		index++;
		return *this;
	}
	TriangularIndices& operator--() {
		index--;
		return *this;
	}
	TriangularIndices operator++(int) {
		TriangularIndices rc(*this);
		operator++();
		return rc;
	}
	TriangularIndices operator--(int) {
		TriangularIndices rc(*this);
		operator--();
		return rc;
	}
	operator int() const {
		return index;
	}
	bool begin() const {
		return bool(index == 0);
	}
	int degree(int k) const {
		return vectorSum(index2indices[index]);
	}
	bool end() const {
		return bool(index == size());
	}
	indices_type getIndices() const {
		return index2indices[index];
	}
	int indexAt(int k) const {
		return index2indices[index][k];
	}
	void reset() {
		index = 0;
	}
	void set(indices_type const &indices) {
		index = indices2index[indices];
	}
	static constexpr size_t size() {
		static constexpr int Size = []() {
			int num = 1, den = 1;
			for (int k = 0; k < DimCnt; k++) {
				num *= k + MomCnt;
				den *= k + 1;
			}
			return num / den;
		}();
		return Size;
	}
private:
	int index;
	static std::array<Math::Vector<int, DimCnt>, TriangularIndices<DimCnt, MomCnt>::size()> const index2indices;
	static map_type const indices2index;
};

template<typename T, int D, int N>
using TriangularArray = Math::Vector<T, TriangularIndices<D, N>::size()>;

template<int DimCnt, int MomCnt>
std::array<Math::Vector<int, DimCnt>, TriangularIndices<DimCnt, MomCnt>::size()> const TriangularIndices<DimCnt, MomCnt>::index2indices = []() {
	using namespace Math;
	std::array<Vector<int, DimCnt>, size()> indicesList;
	int currentIndex = 0;
	for (int pOrder = 0; pOrder < MomCnt; pOrder++) {
		Vector<int, DimCnt> vIndices = zeroVector<int, DimCnt>();
		bool continueLoop;
		do {
			int k = DimCnt - 1;
			vIndices[k]++;
			while ((k > 0) && (vIndices[k] > vIndices[k - 1])) {
				vIndices[k--] = 0;
				vIndices[k]++;
			}
			if (vIndices[0] >= MomCnt) {
				continueLoop = false;
			} else {
				if (pOrder == vectorSum(vIndices)) {
					indicesList[currentIndex++] = vIndices;
				}
				continueLoop = true;
			}
		} while (continueLoop);
	}
	return indicesList;
}();

template<int DimCnt, int MomCnt>
TriangularIndices<DimCnt, MomCnt>::map_type indices2index = [](
		std::array<Math::Vector<int, DimCnt>, TriangularIndices<DimCnt, MomCnt>::Size> const &index2indices) {
	static constexpr size_t maxValue = size_t(1 << 16);
	static constexpr size_t Size = TriangularIndices<DimCnt, MomCnt>::Size;
	using namespace Math;
	using hash_type = vectorHashKey<int, DimCnt, maxValue>;
	using value_type = Vector<int, DimCnt>;
	std::unordered_map<value_type, int, hash_type> indexMap;
	for (size_t i = 0; i != Size; ++i) {
		indexMap[index2indices[i]] = i;
	}
	return indexMap;
}(TriangularIndices<DimCnt, MomCnt>::index2indices);

#endif
