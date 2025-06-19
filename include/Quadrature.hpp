#pragma once
#ifdef ABCD
#include <cmath>
#include <numbers>

#include "Legendre.hpp"
#include "MultiIndex.hpp"
#include "Util.hpp"

template<typename Type, int oneDimensionalPointCount>
using QuadratureReturnType = std::tuple<std::array<Type, oneDimensionalPointCount>, std::array<Type, oneDimensionalPointCount>>;

template<typename Type, int oneDimensionalPointCount>
constexpr auto GaussLegendreQuadrature() {
	constexpr Type pi = std::numbers::pi_v<Type>;
	std::array<Type, oneDimensionalPointCount> weights;
	std::array<Type, oneDimensionalPointCount> points;
	for (int pointIndex = 0; pointIndex < oneDimensionalPointCount; pointIndex++) {
		std::array<std::array<Type, oneDimensionalPointCount>, 2> legendrePolynomial;
		Type theta, position;
		Type newTheta = pi * (Type(1) - (Type(pointIndex) + Type(0.5)) / Type(oneDimensionalPointCount));
		do {
			theta = newTheta;
			position = cos(theta);
			legendrePolynomial = dMLegendrePdXm<Type, oneDimensionalPointCount, 1>(position);
			newTheta = theta + legendrePolynomial[0].back() / (sin(theta) * legendrePolynomial[1].back());
		} while (newTheta != std::nextafter(theta, newTheta));
		points[pointIndex] = position;
		weights[pointIndex] = Type(2) / (sqr(legendrePolynomial[1].back()) * (Type(1) - position * position));
	}
	return QuadratureReturnType<Type, oneDimensionalPointCount>(points, weights);
}

template<typename Type, int oneDimensionalPointCount, int dimensionCount>
struct Quadrature {
	using IndexType = MultiIndex<Range<int, dimensionCount> {repeat<dimensionCount>(0), repeat<dimensionCount>(oneDimensionalPointCount)}>;
	constexpr Type const& weight(int i) const {
		return std::get<1>(pointsAndWeights_)[i];
	}
	constexpr auto const& point(int i) const {
		return std::get<0>(pointsAndWeights_)[i];
	}
	static constexpr int size() {
		return IndexType::count();
	}
private:
	static constexpr auto initialize() {
		auto [oneDimensionalPoints, oneDimensionalWeights] = GaussLegendreQuadrature<Type, oneDimensionalPointCount>();
		std::array<Type, size()> multiDimensionalWeights;
		std::array<std::array<Type, dimensionCount>, size()> multiDimensionalPoints;
		for (int dimension = 0; dimension < dimensionCount; dimension++) {
			for (auto multiIndex = IndexType::begin(); multiIndex != IndexType::end(); multiIndex++) {
				int const pi = multiIndex;
				multiDimensionalWeights[pi] = Type(1);
				for (int d = 0; d < dimensionCount; d++) {
					multiDimensionalWeights[pi] *= oneDimensionalWeights[multiIndex[d]];
					multiDimensionalPoints[pi][d] = oneDimensionalPoints[multiIndex[d]];
				}
			}
		}
		return std::tuple(multiDimensionalPoints, multiDimensionalWeights);
	}
	static constexpr std::tuple<std::array<std::array<Type, dimensionCount>, size()>, std::array<Type, size()>> pointsAndWeights_ = initialize();
};

template<typename Type, int oneDimensionalPointCount>
struct Quadrature<Type, oneDimensionalPointCount, 0> {
	using IndexType = MultiIndex<Range<int, 0> {repeat<0>(0), repeat<0>(1)}>;
	constexpr Type const& weight(int i) const {
		return Type(1);
	}
	constexpr auto point(int i) const {
		return std::array<Type, 0> { };
	}
	static constexpr int size() {
		return 1;
	}
};

template<typename Basis>
constexpr auto fourierLegendreTransform() {
	using Type = Basis::Type;
	constexpr Basis basis { };
	constexpr int pointCount = basis.size();
	constexpr int dimensionCount = basis.dimensionCount;
	constexpr int oneDimensionalPointCount = basis.basisOrder;
	constexpr auto massMatrix = basis.massMatrix();
	constexpr Quadrature<Type, oneDimensionalPointCount, dimensionCount> quadrature { };
	SquareMatrix<Type, pointCount> transformMatrix;
	for (int quadratureIndex = 0; quadratureIndex < pointCount; quadratureIndex++) {
		auto const basisValues = basis(quadrature.point(quadratureIndex));
		for (int basisIndex = 0; basisIndex < pointCount; basisIndex++) {
			auto const mass = massMatrix(basisIndex, basisIndex);
			auto const weight = quadrature.weight(quadratureIndex);
			transformMatrix(basisIndex, quadratureIndex) = mass * weight * basisValues[basisIndex];
		}
	}
	return transformMatrix;
}

template<typename Basis>
constexpr auto inverseFourierLegendreTransform() {
	using Type = Basis::Type;
	constexpr Basis basis { };
	constexpr int pointCount = basis.size();
	constexpr int dimensionCount = basis.dimensionCount;
	constexpr int oneDimensionalPointCount = basis.basisOrder;
	constexpr Quadrature<Type, oneDimensionalPointCount, dimensionCount> quadrature { };
	SquareMatrix<Type, pointCount> inverseTransformMatrix;
	for (int quadratureIndex = 0; quadratureIndex < pointCount; quadratureIndex++) {
		auto const basisValues = basis(quadrature.point(quadratureIndex));
		for (int basisIndex = 0; basisIndex < pointCount; basisIndex++) {
			inverseTransformMatrix(quadratureIndex, basisIndex) = basisValues[basisIndex];
		}
	}
	return inverseTransformMatrix;
}
#endif
