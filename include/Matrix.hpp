/*
 * Matrix.hpp
 *
 *  Created on: Dec 12, 2024
 *      Author: dmarce1
 */

#ifndef INCLUDE_MATRIX_HPP_
#define INCLUDE_MATRIX_HPP_

#include <algorithm>
#include <array>
#include <iostream>
#include <string>
#include <type_traits>

#include "Numbers.hpp"
#include "Utilities.hpp"
#include "Vector.hpp"

namespace Math {

template<typename, int, int>
struct Matrix;

template<typename Type, int Ndim>
using SquareMatrix = Matrix<Type, Ndim, Ndim>;

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
	auto begin() const {
		return values.begin();
	}
	auto end() const {
		return values.end();
	}
	auto begin() {
		return values.begin();
	}
	auto end() {
		return values.end();
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
private:
	Type value;
};

template<typename T, int N, int M, int L>
Matrix<T, N, L> operator*(Matrix<T, N, M> const &A, Matrix<T, M, L> const &B) {
	Matrix<T, N, L> C;
	for (int n = 0; n < N; n++) {
		for (int l = 0; l < L; l++) {
			C[n, l] = T(0);
			for (int m = 0; m < M; m++) {
				C[n, l] += A[n, m] * B[m, l];
			}
		}
	}
	return C;
}

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
	Vector<Type, ColumnCount> row;
	for (int c = 0; c < ColumnCount; c++) {
		row[c] = A[r, c];
	}
	return row;
}

template<typename Type, int RowCount, int ColumnCount>
auto matrixColumn(Matrix<Type, RowCount, ColumnCount> const &A, int c) {
	Vector<Type, RowCount> column;
	for (int r = 0; r < RowCount; r++) {
		column[r] = A[r, c];
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
		Type sum = Type(0);
		Type result = Type(0);
		Type sgn = Type(1);
		for (int c = 0; c < Ndim; c++) {
			result += sgn * A[0, c] * matrixDeterminant(subMatrix(A, 0, c));
			sgn = -sgn;
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
			auto const sgn = ((r + c) & 1) ? -Type(1) : Type(1);
			aCofactor[r, c] = sgn * matrixDeterminant(subMatrix(A, r, c));
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

template<typename T, int N>
void matrixSwapRow(SquareMatrix<T, N> &A, int row1, int row2) {
	for (int col = 0; col < N; col++) {
		std::swap(A[row1, col], A[row2, col]);
	}
}

template<typename T, int N>
void matrixScaleRow(SquareMatrix<T, N> &A, int row, T scale) {
	for (int col = 0; col < N; col++) {
		A[row, col] *= scale;
	}
}

template<typename T, int N>
void matrixSubtractRow(SquareMatrix<T, N> &A, int row1, int row2) {
	for (int col = 0; col < N; col++) {
		A[row1, col] -= A[row2, col];
	}
}

template<typename T, int N>
void matrixSubtractRow(SquareMatrix<T, N> &A, int row1, int row2, T coefficient) {
	for (int col = 0; col < N; col++) {
		A[row1, col] -= coefficient * A[row2, col];
	}
}

template<typename T, int N>
struct matrixReduction_t {
	T det;
};

template<typename T, int N>
matrixReduction_t<T, N> matrixRowReduction(SquareMatrix<T, N> &A, SquareMatrix<T, N> &D) {
	D = identityMatrix<T, N>();
	T det = T(1);
	for (int l = 0; l < N; l++) {
		T scale = A[l, l];
		T scaleInv = T(1) / scale;
		det *= scale;
		matrixScaleRow(A, l, scaleInv);
		std::cout << to_string(A) << "\n\n";
		matrixScaleRow(D, l, scaleInv);
		for (int n = l + 1; n < N; n++) {
			scale = A[n, l];
			if (scale != T(0)) {
				scaleInv = T(1) / scale;
				det *= scale;
				matrixScaleRow(A, n, scaleInv);
				std::cout << to_string(A) << "\n\n";
				matrixScaleRow(D, n, scaleInv);
				matrixSubtractRow(A, n, l);
				std::cout << to_string(A) << "\n\n";
				matrixSubtractRow(D, n, l);
			}
		}
	}
	for (int l = 0; l < N; l++) {
		det *= A[l, l];
	}
	for (int l = 0; l < N - 1; l++) {
		for (int n = l + 1; n < N; n++) {
			T scale = A[l, n];
			matrixSubtractRow(A, l, n, scale);
			matrixSubtractRow(D, l, n, scale);
		}
	}
	matrixReduction_t<T, N> rc;
	rc.det = det;
	return rc;
}

template<typename Type, int RowCount, int ColumnCount>
std::string to_string(Matrix<Type, RowCount, ColumnCount> const &M) {
	using std::to_string;
	int maxLength = 0;
	std::string result;
	for (int row = 0; row < RowCount; row++) {
		for (int column = 0; column < ColumnCount; column++) {
			auto const string = to_string(M[row, column]);
			maxLength = std::max(maxLength, int(string.size()));
		}
	}
	int const charsPerRow = ColumnCount * (maxLength + 3);
	std::string hLine;
	for (int n = 0; n < charsPerRow; n++) {
		hLine += "-";
	}
	hLine += "-\n";
	result = hLine;
	for (int row = 0; row < RowCount; row++) {
		for (int column = 0; column < ColumnCount; column++) {
			auto string = to_string(M[row, column]);
			while (int(string.size()) < maxLength) {
				string = std::string(" ") + string;
			}
			result += "| ";
			result += string;
			result += " ";
		}
		result += "|\n";
		result += hLine;
	}
	return result;
}

template<typename Type, int N>
Matrix<Type, 1, N> toColumnVector(Vector<Type, N> const &vec) {
	Matrix<Type, N, 1> C;
	for (int r = 0; r < N; r++) {
		C[r, 0] = vec[r];
	}
	return C;
}

template<typename Type, int N>
Matrix<Type, 1, N> toRowVector(Vector<Type, N> const &vec) {
	Matrix<Type, 1, N> R;
	for (int c = 0; c < N; c++) {
		R[0, c] = vec[c];
	}
	return R;
}

template<typename Type, int N>
Vector<Type, N> toVector(Matrix<Type, 1, N> const &R) {
	Vector<Type, N> vec;
	for (int c = 0; c < N; c++) {
		vec[c] = R[0, c];
	}
	return vec;
}

template<typename Type, int N>
Vector<Type, N> toVector(Matrix<Type, N, 1> const &C) {
	Vector<Type, N> vec;
	for (int r = 0; r < N; r++) {
		vec[r] = C[r, 0];
	}
	return vec;
}

template<typename Type, int Ndim>
SquareMatrix<Type, Ndim> vectorTensorProduct(Vector<Type, Ndim> const &A, Vector<Type, Ndim> const &B) {
	return A * matrixTranspose(B);
}

}
#endif /* INCLUDE_MATRIX_HPP_ */
