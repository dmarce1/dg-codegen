/*
 * Matrix.hpp
 *
 *  Created on: Dec 12, 2024
 *      Author: dmarce1
 */

#ifndef INCLUDE_MATRIX_HPP_
#define INCLUDE_MATRIX_HPP_

#include <array>
#include <type_traits>

#include "Numbers.hpp"
#include "Utilities.hpp"

namespace Math {

template<typename, int, int>
struct Matrix;

template<typename Type, int Ndim>
using SquareMatrix = Matrix<Type, Ndim, Ndim>;

template<typename Type, int Ndim>
using ColumnVector = Matrix<Type, Ndim, 1>;

template<typename Type, int Ndim>
using RowVector = Matrix<Type, 1, Ndim>;

template<typename Type, int RowCount, int ColumnCount>
struct Matrix {
	USE_STANDARD_DEFAULTS(Matrix)
	USE_STANDARD_ARITHMETIC(Matrix, Type)
	;/**/
	static constexpr std::size_t size() {
		return ColumnCount * RowCount;
	}
	Matrix(std::initializer_list<Type> const &initList) {
		std::copy(initList.begin(), initList.end(), values.begin());
	}
	Matrix(Type const &init) {
		for (int n = 0; n != size(); n++) {
			values[n] = init;
		}
	}
	Type& operator[](int n, int m) {
		return values[n * ColumnCount + m];
	}
	Type operator[](int n, int m) const {
		return values[n * ColumnCount + m];
	}
	template<int OtherColumnCount>
	auto operator*(Matrix<Type, ColumnCount, OtherColumnCount> const &C) const {
		Matrix const &B = *this;
		Matrix<Type, RowCount, OtherColumnCount> A;
		for (int n = 0; n < RowCount; n++) {
			for (int l = 0; l < OtherColumnCount; l++) {
				A[n, l] = Type(0);
				for (int m = 0; m < ColumnCount; m++) {
					A[n, l] += B[n, m] * C[m, l];
				}
			}
		}
		return A;
	}
private:
	std::array<Type, size()> values;
	Type& operator[](int n) {
		return values[n];
	}
	Type operator[](int n) const {
		return values[n];
	}
};

template<typename Type>
struct Matrix<Type, 1, 1> {
	USE_STANDARD_DEFAULTS(Matrix)
	USE_STANDARD_ARITHMETIC(Matrix, Type)
	;/**/
	static constexpr std::size_t size() {
		return 1;
	}
	Type& operator[](int, int) {
		return value;
	}
	Type operator[](int, int) const {
		return value;
	}
	operator Type() const {
		return value;
	}
	operator Type&() {
		return value;
	}
	template<int RowCount, int ColumnCount>
	auto operator*(Matrix<Type, RowCount, ColumnCount> const &C) const {
		return C * value;
	}
private:
	Type value;
};

template<typename Type, int Ndim>
SquareMatrix<Type, Ndim> operator*=(SquareMatrix<Type, Ndim> &A, SquareMatrix<Type, Ndim> const &C) {
	SquareMatrix<Type, Ndim> const B = A;
	for (int n = 0; n < Ndim; n++) {
		for (int l = 0; l < Ndim; l++) {
			A[n, l] = Type(0);
			for (int m = 0; m < Ndim; m++) {
				A[n, l] += B[n, m] * C[m, l];
			}
		}
	}
	return A;
}

template<typename Type, int RowCount, int ColumnCount>
auto matrixRow(Matrix<Type, RowCount, ColumnCount> const &A, int r) {
	RowVector<Type, ColumnCount> row;
	for (int c = 0; c < ColumnCount; c++) {
		row[r, 0] = A[r, c];
	}
	return row;
}

template<typename Type, int RowCount, int ColumnCount>
auto matrixColumn(Matrix<Type, RowCount, ColumnCount> const &A, int c) {
	ColumnVector<Type, RowCount> column;
	for (int r = 0; r < RowCount; r++) {
		column[0, c] = A[r, c];
	}
	return column;
}

template<typename Type, int RowCount, int ColumnCount>
auto matrixTranspose(Matrix<Type, RowCount, ColumnCount> const &B) {
	Matrix<Type, ColumnCount, RowCount> A;
	for (int r = 0; r < RowCount; r++) {
		for (int c = 0; c < ColumnCount; c++) {
			A[c, r] = B[r, c];
		}
	}
	return A;
}

template<typename Type, int Ndim>
SquareMatrix<Type, Ndim - 1> subMatrix(SquareMatrix<Type, Ndim> const &A, int row, int column) {
	SquareMatrix<Type, Ndim - 1> subMatrix;
	for (int r = 0; r < row; r++) {
		for (int c = 0; c < column; c++) {
			subMatrix[r, c] = A[r, c];
		}
		for (int c = column; c < Ndim - 1; c++) {
			subMatrix[r, c] = A[r, c + 1];
		}
	}
	for (int r = row; r < Ndim - 1; r++) {
		for (int c = 0; c < column; c++) {
			subMatrix[r, c] = A[r + 1, c];
		}
		for (int c = column; c < Ndim - 1; c++) {
			subMatrix[r, c] = A[r + 1, c + 1];
		}
	}
	return subMatrix;
}

template<typename Type, int Ndim>
Type matrixDeterminant(SquareMatrix<Type, Ndim> const &A) {
	if constexpr (Ndim == 1) {
		return A[0, 0];
	} else {
		Type result = Type(0);
		for (int c = 0; c < Ndim; c++) {
			result += nOnePow(c) * A[0, c] * matrixDeterminant(subMatrix(A, 0, c));
		}
		return result;
	}
}

template<typename Type, int Ndim>
Type matrixTrace(SquareMatrix<Type, Ndim> const &A) {
	Type result = Type(0);
	for (int c = 0; c < Ndim; c++) {
		result += A[c, c];
	}
	return result;
}

template<typename Type, int Ndim>
SquareMatrix<Type, Ndim> matrixCofactor(SquareMatrix<Type, Ndim> const &A) {
	SquareMatrix<Type, Ndim> aCofactor;
	for (int r = 0; r < Ndim; r++) {
		for (int c = 0; c < Ndim; c++) {
			aCofactor[r, c] = nOnePow(r + c) * matrixDeterminant(subMatrix(A, r, c));
		}
	}
	return aCofactor;
}

template<typename Type, int Ndim>
SquareMatrix<Type, Ndim> matrixAdjoint(SquareMatrix<Type, Ndim> const &A) {
	return matrixTranspose(matrixCofactor(A));
}

template<typename Type, int Ndim>
SquareMatrix<Type, Ndim> matrixInverse(SquareMatrix<Type, Ndim> const &A) {
	return matrixAdjoint(A) * (Type(1) / matrixDeterminant(A));
}

template<typename Type, int P, int Q, int M, int N>
Matrix<Type, P * M, Q * N> matrixKroneckerProduct(Matrix<Type, P, Q> const &A, Matrix<Type, M, N> const &B) {
	Matrix<Type, P * M, Q * N> C;
	for (int p = 0; p < P; p++) {
		auto const pQ = Q * p;
		for (int q = 0; q < Q; q++) {
			for (int m = 0; m < M; m++) {
				auto const mN = N * m;
				for (int n = 0; n < N; n++) {
					C[pQ + q, mN + n] = A[p, q] * B[n, m];
				}
			}
		}
	}
	return C;
}

template<typename Type, int Ndim>
SquareMatrix<Type, Ndim> identityMatrix() {
	SquareMatrix<Type, Ndim> identity;
	for (int n = 0; n < Ndim; n++) {
		for (int m = 0; m < Ndim; m++) {
			identity[n, m] = Type(n == m);
		}
	}
	return identity;
}

}
#endif /* INCLUDE_MATRIX_HPP_ */
