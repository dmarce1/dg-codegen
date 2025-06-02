//#pragma once
//
//#include <array>
//#include <cmath>
//#include <cstddef>
//#include <numbers>
//#include <utility>
//
//template<typename Type>
//struct SparseMatrixEntry {
//	int row;
//	int column;
//	Type value;
//};
//
//template<typename Type, int rowCount, int columnCount, SparseMatrixEntry<Type> ... entries>
//struct SparseMatrix {
//	constexpr auto operator*(std::array<Type, columnCount> const &vector) noexcept {
//		std::array<Type, rowCount> result { };
//		(...,( result[ entries.row ] += entries.value * vector[entries.column]));
//		return result;
//	}
//};
//
//template<typename Type>
//constexpr Type legendrePolynomialDerivative(int n, Type const x) {
//	if (std::abs(x) != Type(1)) {
//		return Type(n + 1) * (std::legendre(n + 1, x) - x * std::legendre(n, x)) / (x * x - Type(1));
//	} else {
//		Type sign = (((n & 1) == 1) || (x > Type(0))) ? Type(+1) : Type(-1);
//		return sign * Type((n * (n + 1)) >> 1);
//	}
//}
//
//template<typename Type, int order, int flatIndex, int rowOffset, int columnOffset, int stride>
//constexpr SparseMatrixEntry<Type> legendreTransformEntry1D() {
//	SparseMatrixEntry<Type> result;
//	constexpr Type pi = std::numbers::pi_v<Type>;
//	constexpr int basisIndex = flatIndex / order;
//	constexpr int quadratureIndex = flatIndex % order;
//	Type theta, basisDerivative, position, weight;
//	Type newTheta = pi * (Type(1) - (Type(quadratureIndex) + Type(1) / Type(2)) / Type(order));
//	do {
//		theta = newTheta;
//		position = std::cos(theta);
//		basisDerivative = legendrePolynomialDerivative(order, position);
//		newTheta = theta + std::legendre(order, position) / (std::sin(theta) * basisDerivative);
//	} while (newTheta != std::nextafter(theta, newTheta));
//	weight = Type(2) / (basisDerivative * basisDerivative * (Type(1) - position * position));
//	result.value =
//	result.row = rowOffset + stride * basisIndex;
//	result.column = columnOffset + stride * quadratureIndex;
//	return result;
//}
