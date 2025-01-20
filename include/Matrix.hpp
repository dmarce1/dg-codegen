/*
 * Matrix.hpp
 *
 *  Created on: Dec 12, 2024
 *      Author: dmarce1
 */

#ifndef INCLUDE_MATRIX_HPP_
#define INCLUDE_MATRIX_HPP_

#include "ForwardDeclarations.hpp"
#include "Utilities.hpp"

#include <algorithm>
#include <array>
#include <cassert>
#include <iostream>
#include <string>
#include <stacktrace>
#include <type_traits>

namespace Math {

template<typename Type, int RowCount, int ColumnCount, typename Container = std::array<Type, ColumnCount * RowCount>>
struct Matrix;

template<typename Type, int Ndim, typename Container = std::array<Type, Ndim * Ndim>>
using SquareMatrix = Matrix<Type, Ndim, Ndim, Container>;

template<typename Type, int RowCount, int ColumnCount, typename Container>
struct Matrix {
	static constexpr std::size_t size() {
		return ColumnCount * RowCount;
	}
	Matrix() :
			values(), createContainer(values) {
	}
	Matrix(std::array<std::array<Type, ColumnCount>, RowCount> const &initList) :
			values(), createContainer(values) {
		for (int n = 0; n < RowCount; n++) {
			for (int m = 0; m < ColumnCount; m++) {
				operator[](n, m) = initList[n][m];
			}
		}
	}
	Matrix(Type const &init) :
			values(), createContainer(values) {
		for (int n = 0; n != size(); n++) {
			values[n] = init;
		}
	}
	Matrix(Matrix const &other) :
			values(), createContainer(values) {
		values = other.values;
	}
	Matrix(Matrix &&other) :
			values(), createContainer(values) {
		values = std::move(other.values);
	}
	Matrix& operator=(Matrix const &other) {
		values = other.values;
		return *this;
	}
	Matrix& operator=(Matrix &&other) {
		values = std::move(other.values);
		return *this;
	}
	Type& operator[](int n, int m) {
		return values[n * ColumnCount + m];
	}
	Type operator[](int n, int m) const {
		return values[n * ColumnCount + m];
	}
	Matrix& operator+=(Matrix const &A) {
		*this = *this + A;
		return *this;
	}
	Matrix& operator-=(Matrix const &A) {
		*this = *this - A;
		return *this;
	}
	Matrix& operator*=(Type const &a) {
		*this = *this * a;
		return *this;
	}
	Matrix& operator/=(Type const &a) {
		*this = *this / a;
		return *this;
	}
	Matrix operator*(Type const &a) const {
		Matrix B;
		for (int k = 0; k < int(size()); k++) {
			B.values[k] = a * values[k];
		}
		return B;
	}
	Matrix operator/(Type const &a) const {
		static constexpr Type one = Type(1);
		Matrix B;
		Type const aInv = one / a;
		for (int k = 0; k < int(size()); k++) {
			B.values[k] = aInv * values[k];
		}
		return B;
	}
	Matrix operator+() const {
		return *this;
	}
	Matrix operator-() const {
		Matrix B;
		for (int k = 0; k < int(size()); k++) {
			B.values[k] = -values[k];
		}
		return B;
	}
	Matrix operator+(Matrix const &A) const {
		Matrix B;
		for (int k = 0; k < int(size()); k++) {
			B.values[k] = values[k] + A.values[k];
		}
		return B;
	}
	Matrix operator-(Matrix const &A) const {
		Matrix B;
		for (int k = 0; k < int(size()); k++) {
			B.values[k] = values[k] + A.values[k];
		}
		return B;
	}
	Matrix operator==(Matrix const &A) const {
		return bool(this->values == A.values);
	}
	Matrix operator!=(Matrix const &A) const {
		return !operator==(A);
	}
	friend Matrix operator*(Type const &a, Matrix const &B) {
		return B * a;
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
	Container values;
	ContainerResizer<Container, size()> const createContainer;
};

template<typename Type, int RowCount, typename Container>
struct Matrix<Type, RowCount, 1, Container> {
	static constexpr std::size_t size() {
		return RowCount;
	}
	Matrix() :
			values(), createContainer(values) {
	}
	Matrix(Type const &initValue) :
			values(), createContainer(values) {
		std::fill(begin(), end(), initValue);
	}
	Matrix(std::initializer_list<Type> const &initList) :
			values(), createContainer(values) {
		std::copy(initList.begin(), initList.end(), begin());
	}
	Matrix(Matrix const &other) :
			values(), createContainer(values) {
		values = other.values;
	}
	Matrix(Matrix &&other) :
			values(), createContainer(values) {
		values = std::move(other.values);
	}
	Matrix& operator=(Matrix const &other) {
		values = other.values;
		return *this;
	}
	Matrix& operator=(Matrix &&other) {
		values = std::move(other.values);
		return *this;
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
	Matrix& operator+=(Matrix const &A) {
		*this = *this + A;
		return *this;
	}
	Matrix& operator-=(Matrix const &A) {
		*this = *this - A;
		return *this;
	}
	Matrix& operator*=(Type const &a) {
		*this = *this * a;
		return *this;
	}
	Matrix& operator/=(Type const &a) {
		*this = *this / a;
		return *this;
	}
	Matrix operator*(Type const &a) const {
		Matrix B;
		for (int k = 0; k < int(size()); k++) {
			B.values[k] = a * values[k];
		}
		return B;
	}
	Matrix operator/(Type const &a) const {
		static constexpr Type one = Type(1);
		Matrix B;
		Type const aInv = one / a;
		for (int k = 0; k < int(size()); k++) {
			B.values[k] = aInv * values[k];
		}
		return B;
	}
	Matrix operator+() const {
		return *this;
	}
	Matrix operator-() const {
		Matrix B;
		for (int k = 0; k < int(size()); k++) {
			B.values[k] = -values[k];
		}
		return B;
	}
	Matrix operator+(Matrix const &A) const {
		Matrix B;
		for (int k = 0; k < int(size()); k++) {
			B.values[k] = values[k] + A.values[k];
		}
		return B;
	}
	Matrix operator-(Matrix const &A) const {
		Matrix B;
		for (int k = 0; k < int(size()); k++) {
			B.values[k] = values[k] + A.values[k];
		}
		return B;
	}
	Matrix operator==(Matrix const &A) const {
		return bool(this->values == A.values);
	}
	Matrix operator!=(Matrix const &A) const {
		return !operator==(A);
	}
	friend Matrix operator*(Type const &a, Matrix const &B) {
		return B * a;
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
	Container values;
	ContainerResizer<Container, RowCount> const createContainer;
};

template<typename Type, typename Container>
struct Matrix<Type, 1, 1, Container> {
	Matrix() {
	}
	Matrix(std::initializer_list<Type> const &initList) :
			value(*(initList.begin())) {
	}
	Matrix(Type const &init) :
			value(init) {
	}
	Matrix(Matrix const &other) :
			value(other.value) {
	}
	Matrix(Matrix &&other) :
			value(std::move(other.value)) {
	}
	Matrix& operator=(Matrix const &other) {
		value = other.value;
		return *this;
	}
	Matrix& operator=(Matrix &&other) {
		value = std::move(other.value);
		return *this;
	}
	static constexpr size_t size() {
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
	Matrix& operator+=(Matrix const &A) {
		*this = *this + A;
		return *this;
	}
	Matrix& operator-=(Matrix const &A) {
		*this = *this - A;
		return *this;
	}
	Matrix& operator*=(Type const &a) {
		*this = *this * a;
		return *this;
	}
	Matrix& operator/=(Type const &a) {
		*this = *this / a;
		return *this;
	}
	Matrix operator*(Type const &a) const {
		Matrix B;
		B.value = a * value;
		return B;
	}
	Matrix operator/(Type const &a) const {
		static constexpr Type one = Type(1);
		Matrix B;
		Type const aInv = one / a;
		for (int k = 0; k < int(size()); k++) {
			B.value = aInv * value;
		}
		return B;
	}
	Matrix operator+() const {
		return *this;
	}
	Matrix operator-() const {
		Matrix B;
		for (int k = 0; k < int(size()); k++) {
			B.value = -value;
		}
		return B;
	}
	Matrix operator+(Matrix const &A) const {
		Matrix B;
		for (int k = 0; k < int(size()); k++) {
			B.value = value + A.value;
		}
		return B;
	}
	Matrix operator-(Matrix const &A) const {
		Matrix B;
		for (int k = 0; k < int(size()); k++) {
			B.value = value + A.value;
		}
		return B;
	}
	Matrix operator==(Matrix const &A) const {
		return bool(this->values == A.values);
	}
	Matrix operator!=(Matrix const &A) const {
		return !operator==(A);
	}
	friend Matrix operator*(Type const &a, Matrix const &B) {
		return B * a;
	}
private:
	Type value;
};

template<typename Type, int Ndim, typename Container = std::array<Type, Ndim * Ndim> >
constexpr SquareMatrix<Type, Ndim, Container> identityMatrix() {
	SquareMatrix<Type, Ndim> identity;
	for (int n = 0; n < Ndim; n++) {
		for (int m = 0; m < Ndim; m++) {
			identity[n, m] = kroneckerDelta<Type>(n, m);
		}
	}
	return identity;
}

template<typename Type, int Ndim, typename Container = std::array<Type, Ndim * Ndim> >
SquareMatrix<Type, Ndim, Container> zeroMatrix() {
	SquareMatrix<Type, Ndim, Container> Z;
	Type const zero = Type(0);
	for (int n = 0; n < Ndim; n++) {
		for (int m = 0; m < Ndim; m++) {
			Z[n, m] = zero;
		}
	}
	return Z;
}

template<typename T, int N, int M, int L, typename Container2, typename Container3>
Matrix<T, N, L> operator*(Matrix<T, N, M, Container2> const &A, Matrix<T, M, L, Container3> const &B) {
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

template<typename Type, int Ndim, typename Container>
SquareMatrix<Type, Ndim, Container> operator*=(SquareMatrix<Type, Ndim, Container> &A,
		SquareMatrix<Type, Ndim, Container> const &C) {
	SquareMatrix<Type, Ndim, Container> const B = A;
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

template<typename Type, int RowCount, int ColumnCount, typename Container>
auto matrixRow(Matrix<Type, RowCount, ColumnCount, Container> const &A, int r) {
	Matrix<Type, 1, ColumnCount, Container> row;
	for (int c = 0; c < ColumnCount; c++) {
		row[0, c] = A[r, c];
	}
	return row;
}

template<typename Type, int RowCount, int ColumnCount, typename Container>
auto matrixColumn(Matrix<Type, RowCount, ColumnCount, Container> const &A, int c) {
	Matrix<Type, RowCount, 1, Container> column;
	for (int r = 0; r < RowCount; r++) {
		column[r, 0] = A[r, c];
	}
	return column;
}

template<typename Type, int RowCount, int ColumnCount, typename Container>
Matrix<Type, ColumnCount, RowCount, Container> matrixTranspose(
		Matrix<Type, RowCount, ColumnCount, Container> const &B) {
	Matrix<Type, ColumnCount, RowCount, Container> A;
	for (int r = 0; r < RowCount; r++) {
		for (int c = 0; c < ColumnCount; c++) {
			A[c, r] = B[r, c];
		}
	}
	return A;
}

template<typename Type, int Ndim, typename Container>
SquareMatrix<Type, Ndim - 1, Container> subMatrix(SquareMatrix<Type, Ndim, Container> const &A, int row, int column) {
	SquareMatrix<Type, Ndim - 1, Container> subMatrix;
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

template<typename T, int R, typename Container>
T matrixAndDeterminateInverse(SquareMatrix<T, R, Container> &A) {
	T constexpr zero(0), one(1);
	T matrixDeterminant = one;
	auto D = identityMatrix<T, R, Container>();
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
	A = D;
	return matrixDeterminant;
}

template<typename T, int R, typename Container>
auto matrixInverse(SquareMatrix<T, R, Container> const A) {
	auto iA = A;
	matrixAndDeterminateInverse(iA);
	return iA;
}

template<typename Type, int Ndim, typename Container>
Type matrixDeterminant(SquareMatrix<Type, Ndim, Container> A) {
	if constexpr (Ndim == 1) {
		return A[0, 0];
	} else {
		return matrixInverseAndDeterminate(A);
	}
}

template<typename Type, int Ndim, typename Container>
Type matrixTrace(SquareMatrix<Type, Ndim, Container> const &A) {
	Type result = Type(0);
	for (int c = 0; c < Ndim; c++) {
		result += A[c, c];
	}
	return result;
}

template<typename Type, int Ndim, typename Container>
SquareMatrix<Type, Ndim, Container> matrixCofactor(SquareMatrix<Type, Ndim, Container> const &A) {
	SquareMatrix<Type, Ndim, Container> aCofactor;
	for (int r = 0; r < Ndim; r++) {
		for (int c = 0; c < Ndim; c++) {
			auto const sgn = ((r + c) & 1) ? -Type(1) : Type(1);
			aCofactor[r, c] = sgn * matrixDeterminant(subMatrix(A, r, c));
		}
	}
	return aCofactor;
}

template<typename Type, int Ndim, typename Container>
SquareMatrix<Type, Ndim, Container> matrixAdjoint(SquareMatrix<Type, Ndim, Container> const &A) {
	return matrixTranspose(matrixCofactor(A));
}

template<typename Type, int P, int Q, int M, int N, typename Container>
Matrix<Type, P * M, Q * N, Container> matrixKroneckerProduct(Matrix<Type, P, Q, Container> const &A,
		Matrix<Type, M, N, Container> const &B) {
	Matrix<Type, P * M, Q * N, Container> C;
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

template<typename T, int N, typename Container>
void matrixQRDecomposition(SquareMatrix<T, N, Container> const &A, SquareMatrix<T, N, Container> &Q,
		SquareMatrix<T, N, Container> &R) {
	T const one(1);
	R = A;
	Q = identityMatrix<T, N, Container>();
	for (int j = 0; j < N; j++) {
		for (int i = j + 1; i < N; i++) {
			SquareMatrix<T, N, Container> G = identityMatrix<T, N>();
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

template<typename T, int R, typename Container>
SquareMatrix<T, R, Container> matrixLUDecomposition(SquareMatrix<T, R, Container> const &A,
		SquareMatrix<T, R, Container> &L, SquareMatrix<T, R, Container> &U) {
	const T zero(0);
	SquareMatrix<T, R> P = identityMatrix<T, R, Container>();
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

template<typename Type, int RowCount, int ColumnCount, typename Container>
std::string toString(Matrix<Type, RowCount, ColumnCount, Container> const &M) {
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

#include "Numbers.hpp"

#endif /* INCLUDE_MATRIX_HPP_ */
