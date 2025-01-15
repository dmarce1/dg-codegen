/*
 * TriangularArray.hpp
 *
 *  Created on: Jan 14, 2025
 *      Author: dmarce1
 */

#ifndef INCLUDE_TRIANGULARARRAY_HPP_
#define INCLUDE_TRIANGULARARRAY_HPP_

#include "Vector.hpp"

constexpr int triangularArraySize(int Ndim, int M) {
	int num = 1, den = 1;
	for (int k = 0; k < Ndim; k++) {
		num *= M + k;
		den *= k + 1;
	}
	return num / den;
}

template<typename Type, int Ndim, int M>
struct TriangularArray: public Math::Vector<Type, triangularArraySize(Ndim, M)> {
	using base_type = Math::Vector<Type, triangularArraySize(Ndim, M)>;
	static constexpr size_t size() {
		return triangularArraySize(Ndim, M);
	}
	constexpr TriangularArray() :
			dataVector(reinterpret_cast<base_type&>(*this)) {
	}
	TriangularArray(TriangularArray const &other) :
			dataVector(reinterpret_cast<base_type&>(*this)) {
		*this = other;
	}
	TriangularArray(TriangularArray &&other) :
			dataVector(reinterpret_cast<base_type&>(*this)) {
		*this = std::move(other);
	}
	TriangularArray& operator=(TriangularArray const &other) {
		dataVector = other.dataVector;
		return *this;
	}
	TriangularArray& operator=(TriangularArray &&other) {
		dataVector = std::move(other.dataVector);
		return *this;
	}
	Type operator[](Math::Vector<int, Ndim> const &arrayIndices) const {
		return dataVector[computeIndex(arrayIndices)];
	}
	Type& operator[](Math::Vector<int, Ndim> const &arrayIndices) {
		return dataVector[computeIndex(arrayIndices)];
	}
private:
	static constexpr int computeIndex(Math::Vector<int, Ndim> const &arrayIndices) {
		int returnIndex = 0;
		for (int k = 0; k < Ndim; k++) {
			int num = 1, den = 1, reducedIndex = 0;
			for (int n = 0; n < Ndim - k; n++) {
				reducedIndex += arrayIndices[n];
			}
			for (int n = 0; n < Ndim - k; n++) {
				num *= reducedIndex + n;
				den *= n + 1;
			}
			returnIndex += num / den;
		}
		return returnIndex;
	}
	base_type &dataVector;
};

#endif
