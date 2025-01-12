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
#include <stacktrace>
#include <type_traits>

#include "Utilities.hpp"

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
	Matrix(std::array<std::array<Type, ColumnCount>, RowCount> const &initList) {
		for (int n = 0; n < RowCount; n++) {
			for (int m = 0; m < ColumnCount; m++) {
				operator[](n, m) = initList[n][m];
			}
		}
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
};

template<typename Type, int RowCount>
struct Matrix<Type, RowCount, 1> {
	USE_STANDARD_DEFAULTS(Matrix)
	USE_STANDARD_ARITHMETIC(Matrix, Type)
	;/**/
	static constexpr std::size_t size() {
		return RowCount;
	}
	Matrix(std::array<std::array<Type, 1>, RowCount> const &initList) {
		for (int n = 0; n < RowCount; n++) {
			operator[](n) = initList[n][0];
		}
	}
	Matrix(std::array<Type, RowCount> const &initList) {
		for (int n = 0; n < RowCount; n++) {
			operator[](n) = initList[n];
		}
	}
	Matrix(Type const &init) {
		for (int n = 0; n != size(); n++) {
			values[n] = init;
		}
	}
	Type& operator[](int n, int m) {
		assert(m == 0);
		return operator[](n);
	}
	Type operator[](int n, int m) const {
		assert(m == 0);
		return operator[](n);
	}
	Type operator[](int n) const {
		return values[n];
	}
	Type& operator[](int n) {
		return values[n];
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
	std::array<Type, RowCount> values;
};

template<typename Type>
struct Matrix<Type, 1, 1> {
	USE_STANDARD_DEFAULTS(Matrix)
	USE_STANDARD_ARITHMETIC(Matrix, Type)
	;/**/
	Matrix(Type const &init) {
		value = init;
	}
	static constexpr std::size_t size() {
		return 1;
	}
	Type& operator[](int, int) {
		return value;
	}
	Type operator[](int, int) const {
		return value;
	}
	Type operator[](int n) const {
		return value;
	}
	Type& operator[](int n) {
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

template<typename Type, int Ndim>
using Vector = Matrix<Type, Ndim, 1>;

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

template<typename Type, int Ndim>
SquareMatrix<Type, Ndim> zeroMatrix() {
	SquareMatrix<Type, Ndim> Z;
	Type const zero = Type(0);
	for (int n = 0; n < Ndim; n++) {
		for (int m = 0; m < Ndim; m++) {
			Z[n, m] = zero;
		}
	}
	return Z;
}

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
	Matrix<Type, 1, ColumnCount> row;
	for (int c = 0; c < ColumnCount; c++) {
		row[0, c] = A[r, c];
	}
	return row;
}

template<typename Type, int RowCount, int ColumnCount>
auto matrixColumn(Matrix<Type, RowCount, ColumnCount> const &A, int c) {
	Matrix<Type, RowCount, 1> column;
	for (int r = 0; r < RowCount; r++) {
		column[r, 0] = A[r, c];
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

template<typename T, int R>
T matrixInverse(SquareMatrix<T, R> &A) {
	T const zero(0);
	T const one(1);
	T matrixDeterminant = one;
	auto D = identityMatrix<T, R>();
	for (int i = 0; i < R; ++i) {
		T pivot = A[i, i];
		if (pivot == zero) {
			int k = i;
			do {
				k++;
				if (k >= R) {
					return zero;
				}
				pivot = A[k, i];
			} while (pivot == zero);
			bool isSingular = true;
			for (int j = 0; j < R; ++j) {
				if (A[i, j] != zero) {
					isSingular = false;
				}
				std::swap(A[i, j], A[k, j]);
				std::swap(D[i, j], D[k, j]);
			}
			matrixDeterminant = -matrixDeterminant;
			if (isSingular) {
				return zero;
			}
		}
		T const iPivot = one / pivot;
		for (int j = 0; j < R; ++j) {
			A[i, j] *= iPivot;
			D[i, j] *= iPivot;
		}
		matrixDeterminant *= pivot;
		for (int k = 0; k < R; ++k) {
			if (k != i) {
				T const factor = A[k, i];
				for (int j = 0; j < R; ++j) {
					A[k, j] -= factor * A[i, j];
					D[k, j] -= factor * D[i, j];
				}
			}
		}
	}
	return matrixDeterminant;
}

template<typename Type, int Ndim>
Type matrixDeterminant(SquareMatrix<Type, Ndim> A) {
	if constexpr (Ndim == 1) {
		return A[0, 0];
	} else {
		return matrixInverse(A);
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

template<typename T, int N>
void matrixQRDecomposition(SquareMatrix<T, N> const &A, SquareMatrix<T, N> &Q, SquareMatrix<T, N> &R) {
	T const one(1);
	R = A;
	Q = identityMatrix<T, N>();
	for (int j = 0; j < N; j++) {
		for (int i = j + 1; i < N; i++) {
			SquareMatrix<T, N> G = identityMatrix<T, N>();
			T const x = R[j, j];
			T const y = R[i, j];
			T const hinv = one / sqrt(sqr(x) + sqr(y));
			T const c = x * hinv;
			T const s = -y * hinv;
			G[j, j] = G[i, i] = c;
			G[i, j] = s;
			G[j, i] = -s;
			R = G * R;
			Q = G * Q;
		}
	}
	Q = matrixTranspose(Q);
}

template<typename T, int R>
SquareMatrix<T, R> matrixLUDecomposition(SquareMatrix<T, R> const &A, SquareMatrix<T, R> &L, SquareMatrix<T, R> &U) {
	const T zero(0);
	SquareMatrix<T, R> P = identityMatrix<T, R>();
	L = P;
	U = A;
	for (int i = 0; i < R; ++i) {
		int k = i;
		while (U[i, i] == zero) {
			k++;
			for (int j = 0; j < R; ++j) {
				std::swap(U[i, j], U[k, j]);
				std::swap(P[i, j], P[k, j]);
			}
		}
		for (int k = i + 1; k < R; ++k) {
			L[k, i] = U[k, i] / U[i, i];
			for (int j = 0; j < R; ++j) {
				U[k, j] -= L[k, i] * U[i, j];
			}
		}
	}
	return matrixTranspose(P);
}

template<typename Type, int RowCount, int ColumnCount>
std::string toString(Matrix<Type, RowCount, ColumnCount> const &M) {
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

}
#endif /* INCLUDE_MATRIX_HPP_ */
