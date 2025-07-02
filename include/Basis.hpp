#pragma once

#include <cassert>
#include <cmath>

#include "Legendre.hpp"
#include "Matrix.hpp"
#include "TriIndex.hpp"
#include "Util.hpp"

//
//template<int basisOrder, int dimensionCount>
//using BasisIndexType = MultiIndex<basisOrder, dimensionCount>;

//template<int basisOrder, int dimensiontCount>
//using BasisIndexType = TriIndex<basisOrder, dimensiontCount>;
//
//template<typename Type_, int basisOrder_, int dimensionCount_>
//struct Basis {
//	static constexpr int basisOrder = basisOrder_;
//	static constexpr int dimensionCount = dimensionCount_;
//	using Type = Type_;
//	using IndexType = BasisIndexType<basisOrder, dimensionCount>;
//	constexpr auto operator()(std::array<Type, dimensionCount> const &position) const {
//		std::array<Type, size()> basis;
//		basis.fill(Type(1));
//		for (int dimension = 0; dimension < dimensionCount; dimension++) {
//			auto const legendrePolynomial = legendreP<Type, basisOrder>(position[dimension]);
//			for (auto basisIndex = IndexType::begin(); basisIndex != IndexType::end(); basisIndex++) {
//				basis[basisIndex] *= legendrePolynomial[basisIndex[dimension]];
//			}
//		}
//		return basis;
//	}
//	constexpr auto gradient(int gradientDirection, std::array<Type, dimensionCount> const &position) const {
//		std::array<Type, size()> basisGradient;
//		basisGradient.fill(Type(1));
//		for (int dimension = 0; dimension < dimensionCount; dimension++) {
//			auto const legendrePolynomialWithDerivative = dMLegendrePdXm<Type, basisOrder, 1>(position[dimension]);
//			for (auto basisIndex = IndexType::begin(); basisIndex != IndexType::end(); basisIndex++) {
//				basisGradient[basisIndex] *= legendrePolynomialWithDerivative[int(dimension == gradientDirection)][basisIndex[dimension]];
//			}
//		}
//		return basisGradient;
//	}
//	constexpr auto massMatrix() const {
//		DiagonalMatrix<Type, size()> massMatrix;
//		for (auto basisMultiIndex = IndexType::begin(); basisMultiIndex != IndexType::end(); basisMultiIndex++) {
//			int const basisFlatIndex = basisMultiIndex;
//			Type mass = Type(1);
//			for (int dimension = 0; dimension < dimensionCount; dimension++) {
//				mass /= Type(0.5) * Type(2 * basisMultiIndex[dimension] + 1);
//			}
//			massMatrix(basisFlatIndex, basisFlatIndex) = mass;
//		}
//		return massMatrix;
//	}
//	static constexpr int size() {
//		return IndexType::count();
//	}
//};
//
//template<typename Type>
//struct QuadraturePoint {
//	Type position;
//	Type weight;
//};
//
//template<typename Type>
//constexpr Type legendrePolynomial(int n, Type const x) {
//	if (n > 1) {
//		return (Type(2 * n - 1) * x * legendrePolynomial(n - 1, x) - Type(n - 1) * legendrePolynomial(n - 2, x)) / Type(n);
//	} else if (n == 1) {
//		return x;
//	} else {
//		return Type(1);
//	}
//}
//
//template<typename Type>
//constexpr Type legendrePolynomialDerivative(int n, Type const x) {
//	using namespace Math;
//	if (abs(x) != Type(1)) {
//		return Type(n + 1) * (legendrePolynomial(n + 1, x) - x * legendrePolynomial(n, x)) / (x * x - Type(1));
//	} else {
//		Type sign = (((n & 1) == 1) || (x > Type(0))) ? Type(+1) : Type(-1);
//		return sign * Type((n * (n + 1)) >> 1);
//	}
//}
//
//template<typename Type, int basisOrder>
//constexpr QuadraturePoint<Type> gaussLegendreQuadraturePoint(int pointIndex) {
//	constexpr Type pi = std::numbers::pi_v<Type>;
//	QuadraturePoint<Type> point;
//	Type theta, basisDerivative;
//	Type newTheta = pi * (Type(1) - (Type(pointIndex) + Type(0.5)) / Type(basisOrder));
//	do {
//		theta = newTheta;
//		point.position = cos(theta);
//		basisDerivative = legendrePolynomialDerivative(basisOrder, point.position);
//		newTheta = theta + legendrePolynomial(basisOrder, point.position) / (sin(theta) * basisDerivative);
//	} while (newTheta != std::nextafter(theta, newTheta));
//	point.weight = Type(2) / (sqr(basisDerivative) * (Type(1) - sqr(point.position)));
//	return point;
//}
//
//template<typename Type, int basisOrder>
//constexpr QuadraturePoint<Type> gaussLegendreQuadraturePoints() {
//	std::array<QuadraturePoint<Type>, basisOrder> points;
//	for (int basisIndex = 0; basisIndex < basisOrder; basisIndex++) {
//		points[basisIndex] = gaussLegendreQuadraturePoint<Type, basisOrder>(basisIndex);
//	}
//}
//
//enum class TransformDirection : int {
//	forward, backward
//};
//
////template<typename Type, int basisOrder, TransformDirection transformDirection = TransformDirection::forward>
////constexpr auto fourierLegendreTransform() {
////	if constexpr (transformDirection == TransformDirection::forward) {
////		return matrixInverse(fourierLegendreTransform<Type, basisOrder, TransformDirection::backward>());
////	} else {
////		using namespace Math;
////		SquareMatrix<Type, basisOrder> transform;
////		for (int basisIndex = 0; basisIndex < basisOrder; basisIndex++) {
////			for (int quadratureIndex = 0; quadratureIndex < basisOrder; quadratureIndex++) {
////				auto const [quadraturePosition, quadratureWeight] = gaussLegendreQuadraturePoint<Type, basisOrder>(quadratureIndex);
////				transform(basisIndex, quadratureIndex) = legendrePolynomial(basisIndex, quadraturePosition);
////			}
////		}
////		return transform;
////	}
////}
////
////template<typename Type, int basisOrder, int stride>
////auto fourierLegendreTransform1d(Type *v) {
////	constexpr auto A = matrixLUDecompose(fourierLegendreTransform<Type, basisOrder, TransformDirection::forward>());
////	for (int row = 0; row < basisOrder; row += stride) {
////		v[row] *= A(row, row);
////		for (int column = row + 1; column < basisOrder; column++) {
////			v[row] += A(row, column) * v[column];
////		}
////	}
////	for (int row = basisOrder - stride; row > 0; row -= stride) {
////		for (int column = 0; column < row; column++) {
////			v[row] += A(row, column) * v[column];
////		}
////	}
////}
////
////template<typename Type, int basisOrder, int dimensionCount>
////auto fourierLegendreTransform(Type *x) {
////	constexpr int size = ipow(basisOrder, dimensionCount);
////	constexpr int stride = ipow(basisOrder, dimensionCount - 1);
////	if constexpr (dimensionCount > 1) {
////		for (int offset = 0; offset < size; offset += stride) {
////			fourierLegendreTransform<Type, basisOrder, dimensionCount - 1>(x + offset);
////		}
////	}
////	for (int offset = 0; offset < stride; offset++) {
////		fourierLegendreTransform1d<Type, basisOrder, stride>(x + offset);
////	}
////}
//
//#include <array>
//
//
////// --- truncated transform ---------------------------------------------------
////template<typename T, int O, int D, TransformDirection dir = TransformDirection::forward>
////constexpr auto fourierLegendreTransformTruncated(std::array<T, MultiIndex<O, D>::count()> input) {
////	using Indices = MultiIndex<O, D>;
////	constexpr int N = O;
////	constexpr int Size = Indices::count();     // = N^D
////	constexpr auto strides = Indices::strides();  // precomputed
////	constexpr auto T1d = oneDTransformMatrix<T, O, dir>();
////
////	// 1) do the first D-1 dimensions in-place on a full-size buffer
////	std::array<T, Size> bufA = input;
////	std::array<T, Size> bufB { };
////
////	auto *inBuf = &bufA;
////	auto *outBuf = &bufB;
////
////	for (int d = 0; d < D - 1; ++d) {
////		int const s = strides[d];
////		int const blockSz = s * N;
////
////		for (int base = 0; base < Size; base += blockSz) {
////			for (int low = 0; low < s; ++low) {
////				for (int orig = 0; orig < N; ++orig) {
////					int const outFlat = base + orig * s + low;
////					T acc = T(0);
////					for (int q = 0; q < N; ++q) {
////						acc += T1d(orig, q) * (*inBuf)[base + q * s + low];
////					}
////					(*outBuf)[outFlat] = acc;
////				}
////			}
////		}
////		std::swap(inBuf, outBuf);
////	}
////
////	// 2) Now for the final dimension (d = D-1), only do the M valid flats
////	constexpr auto truncList = makeTruncList<O, D>();
////	constexpr int M = truncList.size();
////	std::array<T, M> result { };
////
////	int const dLast = D - 1;
////	int const s = strides[dLast];
////
////	for (auto [flat, outSlot] : truncList) {
////		// sum along dimension dLast for this flat
////		int const orig = (flat / s) % N;
////		int const low = flat % s;
////		T acc = T(0);
////
////		// mat-vec length N on just this line
////		for (int q = 0; q < N; ++q) {
////			int const inFlat = low + q * s + (flat - (orig * s + low));
////			// equivalently: blockBase = flat - orig*s - low
////			acc += T1d(orig, q) * (*inBuf)[inFlat];
////		}
////		result[outSlot] = acc;
////	}
////
////	return result;  // length = binomial(D+O-1, D)
////}
