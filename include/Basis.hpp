#pragma once

#include <cassert>
#include <cmath>

#include "Legendre.hpp"
#include "Matrix.hpp"
#include "TriIndex.hpp"

template<int basisOrder, int dimensionCount>
using BasisIndexType = MultiIndex<basisOrder, dimensionCount>;

template<typename Type_, int basisOrder_, int dimensionCount_>
struct Basis {
	static constexpr int basisOrder = basisOrder_;
	static constexpr int dimensionCount = dimensionCount_;
	using Type = Type_;
	using IndexType = BasisIndexType<basisOrder, dimensionCount>;
	constexpr auto operator()(std::array<Type, dimensionCount> const &position) const {
		std::array<Type, size()> basis;
		basis.fill(Type(1));
		for (int dimension = 0; dimension < dimensionCount; dimension++) {
			auto const legendrePolynomial = legendreP<Type, basisOrder>(position[dimension]);
			for (auto basisIndex = IndexType::begin(); basisIndex != IndexType::end(); basisIndex++) {
					basis[basisIndex] *= legendrePolynomial[basisIndex[dimension]];
			}
		}
		return basis;
	}
	constexpr auto gradient(int gradientDirection, std::array<Type, dimensionCount> const &position) const {
		std::array<Type, size()> basisGradient;
		basisGradient.fill(Type(1));
		for (int dimension = 0; dimension < dimensionCount; dimension++) {
			auto const legendrePolynomialWithDerivative = dMLegendrePdXm<Type, basisOrder, 1>(position[dimension]);
			for (auto basisIndex = IndexType::begin(); basisIndex != IndexType::end(); basisIndex++) {
				basisGradient[basisIndex] *= legendrePolynomialWithDerivative[int(dimension == gradientDirection)][basisIndex[dimension]];
			}
		}
		return basisGradient;
	}
	constexpr auto massMatrix() const {
		DiagonalMatrix<Type, size()> massMatrix;
		for (auto basisMultiIndex = IndexType::begin(); basisMultiIndex != IndexType::end(); basisMultiIndex++) {
			int const basisFlatIndex = basisMultiIndex;
			Type mass = Type(1);
			for (int dimension = 0; dimension < dimensionCount; dimension++) {
				mass /= Type(0.5) * Type(2 * basisMultiIndex[dimension] + 1);
			}
			massMatrix(basisFlatIndex, basisFlatIndex) = mass;
		}
		return massMatrix;
	}
	static constexpr int size() {
		return IndexType::count();
	}
};

