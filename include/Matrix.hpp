#pragma once

#include <algorithm>
#include <array>
#include <cassert>
#include <cmath>
#include <iostream>
#include <string>
#include <stacktrace>
#include <type_traits>

#include "ContainerArithmetic.hpp"

template<typename Type, int RowCount, int ColumnCount>
struct Matrix;

template<typename Type, int Ndim>
using SquareMatrix = Matrix<Type, Ndim, Ndim>;

template<typename Type, int Ndim>
struct DiagonalMatrix;

template<typename T, int D>
struct SymmetricMatrix;

template<typename T>
struct IsMatrix {
	static constexpr bool value = false;
};

template<typename T, int Nr, int Nc>
struct IsMatrix<Matrix<T, Nr, Nc>> {
	static constexpr bool value = true;
};

template<typename T, int N>
struct IsMatrix<SquareMatrix<T, N>> {
	static constexpr bool value = true;
};

template<typename T, int N>
struct IsMatrix<SymmetricMatrix<T, N>> {
	static constexpr bool value = true;
};

template<typename T, int N>
struct IsMatrix<DiagonalMatrix<T, N>> {
	static constexpr bool value = true;
};

template<typename T>
struct IsSquareMatrix {
	static constexpr bool value = false;
};

template<typename T, int N>
struct IsSquareMatrix<SquareMatrix<T, N>> {
	static constexpr bool value = true;
};

template<typename T, int N>
struct IsSquareMatrix<SymmetricMatrix<T, N>> {
	static constexpr bool value = true;
};

template<typename Type, int RowCount, int ColumnCount>
struct Matrix {
	static constexpr std::size_t size() {
		return ColumnCount * RowCount;
	}
	static constexpr std::size_t rowCount() {
		return RowCount;
	}
	static constexpr std::size_t columnCount() {
		return ColumnCount;
	}

	constexpr Matrix() {
	}

	constexpr Matrix(std::array<std::array<Type, ColumnCount>, RowCount> const &initList) :
			values { } {
		for (int n = 0; n < RowCount; n++) {
			for (int m = 0; m < ColumnCount; m++) {
				(*this)(n, m) = initList[n][m];
			}
		}
	}

	constexpr Matrix(std::initializer_list<Type> initList) {
		int i = 0;
		for (auto const &v : initList) {
			values[i++] = v;
		}
	}

	constexpr Matrix(std::initializer_list<std::initializer_list<Type>> init) {
		int i = 0;
		for (auto const &row : init) {
			for (auto const &val : row) {
				values[i++] = val;
			}
		}
	}

	constexpr Matrix(Type const &init) :
			values { } {
		for (std::size_t i = 0; i < size(); i++) {
			values[i] = init;
		}
	}

	constexpr Matrix(Matrix const &other) :
			values(other.values) {
	}
	constexpr Matrix(Matrix &&other) :
			values(std::move(other.values)) {
	}

	constexpr Matrix& operator=(Matrix const &other) {
		values = other.values;
		return *this;
	}
	constexpr Matrix& operator=(Matrix &&other) {
		values = std::move(other.values);
		return *this;
	}

	constexpr Type& operator()(int n, int m) {
		return values[n * ColumnCount + m];
	}
	constexpr Type const& operator()(int n, int m) const {
		return values[n * ColumnCount + m];
	}

	constexpr Matrix& operator+=(Matrix const &A) {
		*this = *this + A;
		return *this;
	}
	constexpr Matrix& operator-=(Matrix const &A) {
		*this = *this - A;
		return *this;
	}
	constexpr Matrix& operator*=(Type const &a) {
		*this = *this * a;
		return *this;
	}
	constexpr Matrix& operator/=(Type const &a) {
		*this = *this / a;
		return *this;
	}

	constexpr Matrix operator*(Type const &a) const {
		Matrix B;
		for (std::size_t k = 0; k < size(); k++) {
			B.values[k] = a * values[k];
		}
		return B;
	}

	constexpr Matrix operator/(Type const &a) const {
		static constexpr Type one = Type(1);
		Matrix B;
		Type const aInv = one / a;
		for (std::size_t k = 0; k < size(); k++) {
			B.values[k] = aInv * values[k];
		}
		return B;
	}

	constexpr Matrix operator+() const {
		return *this;
	}
	constexpr Matrix operator-() const {
		Matrix B;
		for (std::size_t k = 0; k < size(); k++) {
			B.values[k] = -values[k];
		}
		return B;
	}

	constexpr Matrix operator+(Matrix const &A) const {
		Matrix B;
		for (std::size_t k = 0; k < size(); k++) {
			B.values[k] = values[k] + A.values[k];
		}
		return B;
	}

	constexpr Matrix operator-(Matrix const &A) const {
		Matrix B;
		for (std::size_t k = 0; k < size(); k++) {
			B.values[k] = values[k] - A.values[k];
		}
		return B;
	}

	constexpr bool operator==(Matrix const &A) const {
		return values == A.values;
	}
	constexpr bool operator!=(Matrix const &A) const {
		return !(*this == A);
	}

	static constexpr std::enable_if_t<RowCount == ColumnCount, Matrix> identity() {
		Matrix I;
		for (int n = 0; n < RowCount; n++) {
			for (int m = 0; m < ColumnCount; m++) {
				I(n, m) = Type(n == m);
			}
		}
		return I;
	}

	static constexpr Matrix zero() {
		Matrix Z;
		std::fill(Z.begin(), Z.end(), Type(0));
		return Z;
	}

	friend Matrix operator*(Type const &a, Matrix const &B) {
		return B * a;
	}

	auto begin() {
		return values.begin();
	}
	auto end() {
		return values.end();
	}
	auto begin() const {
		return values.begin();
	}
	auto end() const {
		return values.end();
	}

private:
	std::array<Type, size()> values;
};

template<typename Type, int RowCount>
struct Matrix<Type, RowCount, 1> {
	static constexpr std::size_t size() {
		return RowCount;
	}
	static constexpr std::size_t rowCount() {
		return RowCount;
	}
	static constexpr std::size_t columnCount() {
		return 1;
	}

	Matrix() :
			values { } {
	}
	Matrix(Type const &init) :
			values { } {
		std::fill(begin(), end(), init);
	}
	Matrix(std::initializer_list<Type> initList) :
			values { } {
		std::copy(initList.begin(), initList.end(), begin());
	}
	Matrix(Matrix const &other) :
			values(other.values) {
	}
	Matrix(Matrix &&other) :
			values(std::move(other.values)) {
	}

	Matrix& operator=(Matrix const &other) {
		values = other.values;
		return *this;
	}
	Matrix& operator=(Matrix &&other) {
		values = std::move(other.values);
		return *this;
	}

	Type& operator()(int n, int m) {
		assert(m == 0);
		return (*this)(n);
	}
	Type const& operator()(int n, int m) const {
		assert(m == 0);
		return (*this)(n);
	}

	Type& operator()(int n) {
		return values[n];
	}
	Type const& operator()(int n) const {
		return values[n];
	}

	Matrix& operator+=(Matrix const &A) {
		return *this = *this + A;
	}
	Matrix& operator-=(Matrix const &A) {
		return *this = *this - A;
	}
	Matrix& operator*=(Type const &a) {
		return *this = *this * a;
	}
	Matrix& operator/=(Type const &a) {
		return *this = *this / a;
	}

	Matrix operator*(Type const &a) const {
		Matrix B;
		for (std::size_t k = 0; k < size(); k++)
			B.values[k] = a * values[k];
		return B;
	}
	Matrix operator/(Type const &a) const {
		static constexpr Type one = Type(1);
		Matrix B;
		Type const aInv = one / a;
		for (std::size_t k = 0; k < size(); k++)
			B.values[k] = aInv * values[k];
		return B;
	}

	Matrix operator+() const {
		return *this;
	}
	Matrix operator-() const {
		Matrix B;
		for (std::size_t k = 0; k < size(); k++)
			B.values[k] = -values[k];
		return B;
	}
	Matrix operator+(Matrix const &A) const {
		Matrix B;
		for (std::size_t k = 0; k < size(); k++)
			B.values[k] = values[k] + A.values[k];
		return B;
	}
	Matrix operator-(Matrix const &A) const {
		Matrix B;
		for (std::size_t k = 0; k < size(); k++)
			B.values[k] = values[k] - A.values[k];
		return B;
	}

	bool operator==(Matrix const &A) const {
		return values == A.values;
	}
	bool operator!=(Matrix const &A) const {
		return !(*this == A);
	}

	static constexpr Matrix zero() {
		Matrix Z;
		std::fill(Z.begin(), Z.end(), Type(0));
		return Z;
	}

	friend Matrix operator*(Type const &a, Matrix const &B) {
		return B * a;
	}

	auto begin() {
		return values.begin();
	}
	auto end() {
		return values.end();
	}
	auto begin() const {
		return values.begin();
	}
	auto end() const {
		return values.end();
	}

	auto* data() {
		return values.data();
	}
	auto const* data() const {
		return values.data();
	}

private:
	std::array<Type, size()> values;
};

template<typename Type>
struct Matrix<Type, 1, 1> {
	static constexpr std::size_t rowCount() {
		return 1;
	}
	static constexpr std::size_t columnCount() {
		return 1;
	}

	Matrix() = default;
	Matrix(std::initializer_list<Type> const &init) :
			value(*init.begin()) {
	}
	Matrix(Type const &v) :
			value(v) {
	}
	Matrix(Matrix const &o) :
			value(o.value) {
	}
	Matrix(Matrix &&o) :
			value(std::move(o.value)) {
	}

	Matrix& operator=(Matrix const &o) {
		value = o.value;
		return *this;
	}
	Matrix& operator=(Matrix &&o) {
		value = std::move(o.value);
		return *this;
	}

	static constexpr std::size_t size() {
		return 1;
	}

	Type& operator()(int, int) {
		return value;
	}
	Type const& operator()(int, int) const {
		return value;
	}

	Type& operator()(int) {
		return value;
	}
	Type const& operator()(int) const {
		return value;
	}

	operator Type() const {
		return value;
	}
	operator Type&() {
		return value;
	}

	Matrix& operator+=(Matrix const &A) {
		return *this = *this + A;
	}
	Matrix& operator-=(Matrix const &A) {
		return *this = *this - A;
	}
	Matrix& operator*=(Type const &a) {
		return *this = *this * a;
	}
	Matrix& operator/=(Type const &a) {
		return *this = *this / a;
	}

	Matrix operator*(Type const &a) const {
		Matrix B;
		B.value = a * value;
		return B;
	}
	Matrix operator/(Type const &a) const {
		static constexpr Type one = Type(1);
		Matrix B;
		B.value = (one / a) * value;
		return B;
	}

	Matrix operator+() const {
		return *this;
	}
	Matrix operator-() const {
		return Matrix(-value);
	}

	Matrix operator+(Matrix const &A) const {
		return Matrix(value + A.value);
	}
	Matrix operator-(Matrix const &A) const {
		return Matrix(value - A.value);
	}

	bool operator==(Matrix const &A) const {
		return value == A.value;
	}
	bool operator!=(Matrix const &A) const {
		return !(*this == A);
	}

	friend Matrix operator*(Type const &a, Matrix const &B) {
		return B * a;
	}

	static constexpr Matrix identity() {
		return Matrix(Type(1));
	}
	static constexpr Matrix zero() {
		return Matrix(Type(0));
	}

private:
	Type value;
};

template<typename T, int D>
SymmetricMatrix<T, D> matrixSymmetrize(SquareMatrix<T, D> const &A);

template<typename T, int D>
struct SymmetricMatrix: public std::array<T, ((D * D + D) >> 1)> {
	using base_type = std::array<T, ((D * D + D) >> 1)>;

	SymmetricMatrix() = default;
	SymmetricMatrix(SymmetricMatrix const&) = default;
	SymmetricMatrix(SymmetricMatrix&&) = default;
	SymmetricMatrix& operator=(SymmetricMatrix const&) = default;
	SymmetricMatrix& operator=(SymmetricMatrix&&) = default;

	static constexpr std::size_t rowCount() {
		return D;
	}
	static constexpr std::size_t columnCount() {
		return D;
	}

	SymmetricMatrix(T init) {
		std::fill(base_type::begin(), base_type::end(), init);
	}

	SymmetricMatrix(std::initializer_list<std::initializer_list<T>> init) {
		*this = matrixSymmetrize(SquareMatrix<T, D>(init));
	}

	T& operator()(int n, int k) {
		return base_type::operator[](index(n, k));
	}
	T const& operator()(int n, int k) const {
		return base_type::operator[](index(n, k));
	}

	operator SquareMatrix<T, D>() const {
		SquareMatrix<T, D> M;
		for (int p = 0; p < D; p++) {
			for (int q = p; q < D; q++) {
				auto v = base_type::operator[](((p * p + p) >> 1) + q);
				M(p, q) = M(q, p) = v;
			}
		}
		return M;
	}

	static constexpr SymmetricMatrix identity() {
		SymmetricMatrix I;
		for (int n = 0; n < D; n++) {
			for (int m = 0; m <= n; m++) {
				I(n, m) = T(n == m);
			}
		}
		return I;
	}

	static constexpr SymmetricMatrix zero() {
		SymmetricMatrix Z;
		std::fill(Z.begin(), Z.end(), T(0));
		return Z;
	}

	SymmetricMatrix& operator+=(SymmetricMatrix const &A) {
		*this = *this + A;
		return *this;
	}
	SymmetricMatrix& operator-=(SymmetricMatrix const &A) {
		*this = *this - A;
		return *this;
	}
	SymmetricMatrix& operator*=(T const &a) {
		*this = *this * a;
		return *this;
	}

	friend SymmetricMatrix operator+(SymmetricMatrix const &A, SymmetricMatrix const &B) {
		SymmetricMatrix C;
		for (std::size_t i = 0; i < A.size(); i++) {
			C.base_type::operator[](i) = A.base_type::operator[](i) + B.base_type::operator[](i);
		}
		return C;
	}

	friend SymmetricMatrix operator-(SymmetricMatrix const &A, SymmetricMatrix const &B) {
		SymmetricMatrix C;
		for (std::size_t i = 0; i < A.size(); i++) {
			C.base_type::operator[](i) = A.base_type::operator[](i) - B.base_type::operator[](i);
		}
		return C;
	}

	friend SymmetricMatrix operator*(SymmetricMatrix const &A, T b) {
		SymmetricMatrix C;
		for (std::size_t i = 0; i < A.size(); i++) {
			C.base_type::operator[](i) = A.base_type::operator[](i) * b;
		}
		return C;
	}

	friend SymmetricMatrix operator*(T a, SymmetricMatrix const &B) {
		return B * a;
	}

	friend SymmetricMatrix operator/(SymmetricMatrix const &A, T b) {
		return A * (T(1) / b);
	}

private:
	static constexpr int index(int n, int k) {
		int p = std::max(n, k), q = std::min(n, k);
		return ((p * p + p) >> 1) + q;
	}
};

template<typename T, int D>
SymmetricMatrix<T, D> matrixSymmetrize(SquareMatrix<T, D> const &A) {
	SymmetricMatrix<T, D> B;
	for (int n = 0; n < D; n++) {
		B(n, n) = A(n, n);
		for (int m = 0; m < n; m++) {
			B(n, m) = T(0.5) * (A(n, m) + A(m, n));
		}
	}
	return B;
}

template<typename T, int N, int M, int L>
Matrix<T, N, L> operator*(Matrix<T, N, M> const &A, Matrix<T, M, L> const &B) {
	Matrix<T, N, L> C;
	for (int n = 0; n < N; n++) {
		for (int l = 0; l < L; l++) {
			C(n, l) = T(0);
			for (int m = 0; m < M; m++) {
				C(n, l) += A(n, m) * B(m, l);
			}
		}
	}
	return C;
}

template<typename Type, int Ndim>
SquareMatrix<Type, Ndim> operator*=(SquareMatrix<Type, Ndim> &A, SquareMatrix<Type, Ndim> const &C) {
	SquareMatrix<Type, Ndim> B = A;
	for (int i = 0; i < Ndim; i++) {
		for (int j = 0; j < Ndim; j++) {
			A(i, j) = Type(0);
			for (int k = 0; k < Ndim; k++) {
				A(i, j) += B(i, k) * C(k, j);
			}
		}
	}
	return A;
}

template<typename Type, int R, int C>
auto matrixRow(Matrix<Type, R, C> const &A, int r) {
	Matrix<Type, 1, C> row;
	for (int c = 0; c < C; c++) {
		row(0, c) = A(r, c);
	}
	return row;
}

template<typename Type, int R, int C>
Matrix<Type, R, 1> matrixColumn(Matrix<Type, R, C> const &A, int c) {
	Matrix<Type, R, 1> col;
	for (int r = 0; r < R; r++) {
		col(r, 0) = A(r, c);
	}
	return col;
}

template<typename Type, int N>
Matrix<Type, N, 1> matrixColumn(SymmetricMatrix<Type, N> const &A, int c) {
	Matrix<Type, N, 1> col;
	for (int r = 0; r < N; r++) {
		col(r, 0) = A(r, c);
	}
	return col;
}

template<typename T, int N>
T matrixInverseAndDeterminant(SquareMatrix<T, N> &A) {
	T constexpr zero(0), one(1);
	T det = one;
	auto D = SquareMatrix<T, N>::identity();
	for (int i = 0; i < N; ++i) {
		T pivot = A(i, i);
		if (pivot == zero) {
			int k = i;
			do {
				k++;
				if (k >= N)
					std::abort();
				pivot = A(k, i);
			} while (pivot == zero);
			bool singular = true;
			for (int j = 0; j < N; ++j) {
				if (A(i, j) != zero)
					singular = false;
				std::swap(A(i, j), A(k, j));
				std::swap(D(i, j), D(k, j));
			}
			det = -det;
			if (singular)
				return zero;
		}
		T invP = one / pivot;
		for (int j = 0; j < N; ++j) {
			A(i, j) *= invP;
			D(i, j) *= invP;
		}
		det *= pivot;
		for (int r = 0; r < N; ++r) {
			if (r == i)
				continue;
			T f = A(r, i);
			for (int j = 0; j < N; ++j) {
				A(r, j) -= f * A(i, j);
				D(r, j) -= f * D(i, j);
			}
		}
	}
	A = D;
	return det;
}

template<typename Type, int Ndim>
Type matrixTrace(SquareMatrix<Type, Ndim> const &A) {
	Type sum = Type(0);
	for (int i = 0; i < Ndim; i++) {
		sum += A(i, i);
	}
	return sum;
}

template<typename Type, int R, int C>
Matrix<Type, C, R> matrixTranspose(Matrix<Type, R, C> const &B) {
	Matrix<Type, C, R> A;
	for (int r = 0; r < R; r++) {
		for (int c = 0; c < C; c++) {
			A(c, r) = B(r, c);
		}
	}
	return A;
}

template<typename Type, int Ndim>
SquareMatrix<Type, Ndim - 1> subMatrix(SquareMatrix<Type, Ndim> const &A, int row, int col) {
	SquareMatrix<Type, Ndim - 1> M;
	for (int r = 0; r < row; r++) {
		for (int c = 0; c < col; c++)
			M(r, c) = A(r, c);
		for (int c = col; c < Ndim - 1; c++)
			M(r, c) = A(r, c + 1);
	}
	for (int r = row; r < Ndim - 1; r++) {
		for (int c = 0; c < col; c++)
			M(r, c) = A(r + 1, c);
		for (int c = col; c < Ndim - 1; c++)
			M(r, c) = A(r + 1, c + 1);
	}
	return M;
}

template<typename T, int N>
T matrixDeterminant(SquareMatrix<T, N> const &A);

template<typename T, int N>
T matrixCofactor(SquareMatrix<T, N> const &A, int r, int c) {
	if constexpr (N > 1) {
		T sgn = ((r + c) & 1) ? -T(1) : T(1);
		return sgn * matrixDeterminant(subMatrix(A, r, c));
	} else {
		return A(0, 0);
	}
}

template<typename T, int N>
SquareMatrix<T, N> matrixCofactor(SquareMatrix<T, N> const &A) {
	SquareMatrix<T, N> C;
	for (int r = 0; r < N; r++) {
		for (int c = 0; c < N; c++) {
			C(r, c) = matrixCofactor(A, r, c);
		}
	}
	return C;
}

template<typename T, int N>
T matrixDeterminant(SquareMatrix<T, N> const &A) {
	if constexpr (N > 1) {
		T sum = T(0);
		for (int c = 0; c < N; c++) {
			sum += A(0, c) * matrixCofactor(A, 0, c);
		}
		return sum;
	} else {
		return A(0, 0);
	}
}

template<typename T, int N>
SquareMatrix<T, N> matrixAdjoint(SquareMatrix<T, N> const &A) {
	return matrixTranspose(matrixCofactor(A));
}

template<typename T, int N>
auto matrixInverse(SquareMatrix<T, N> const &A) {
	return matrixAdjoint(A) / matrixDeterminant(A);
}

template<typename T, int N>
auto matrixInverse(SymmetricMatrix<T, N> const &A) {
	return matrixAdjoint(A) / matrixDeterminant(A);
}

template<typename T, int N>
SymmetricMatrix<T, N - 1> subMatrix(SymmetricMatrix<T, N> const &A, int d) {
	SymmetricMatrix<T, N - 1> M;
	for (int i = 0; i < d; i++) {
		for (int j = 0; j <= i; j++) {
			M(i, j) = A(i, j);
		}
	}
	for (int i = d; i < N - 1; i++) {
		for (int j = 0; j < d; j++)
			M(i, j) = A(i + 1, j);
		for (int j = d; j <= i; j++)
			M(i, j) = A(i + 1, j + 1);
	}
	return M;
}

template<typename T, int N>
SymmetricMatrix<T, N> matrixCofactor(SymmetricMatrix<T, N> const &A) {
	SymmetricMatrix<T, N> C;
	for (int i = 0; i < N; i++) {
		for (int j = 0; j <= i; j++) {
			C(i, j) = matrixCofactor(A, i, j);
		}
	}
	return C;
}

template<typename T, int N>
T matrixDeterminant(SymmetricMatrix<T, N> const &A) {
	if constexpr (N > 1) {
		T sum = T(0);
		for (int c = 0; c < N; c++) {
			sum += A(0, c) * matrixCofactor(A, 0, c);
		}
		return sum;
	} else {
		return A(0, 0);
	}
}

template<typename T, int N>
SymmetricMatrix<T, N> matrixAdjoint(SymmetricMatrix<T, N> const &A) {
	return matrixCofactor(A);
}

template<typename T, int N, int M>
SquareMatrix<T, N> operator*(SymmetricMatrix<T, N> const &A, Matrix<T, N, M> const &B) {
	SquareMatrix<T, N> C;
	for (int i = 0; i < N; i++) {
		for (int j = 0; j < M; j++) {
			C(i, j) = T(0);
			for (int k = 0; k < N; k++) {
				C(i, j) += A(i, k) * B(k, j);
			}
		}
	}
	return C;
}

template<typename T, int N, int M>
SquareMatrix<T, N> operator*(Matrix<T, N, M> const &A, SymmetricMatrix<T, M> const &B) {
	SquareMatrix<T, N> C;
	for (int i = 0; i < N; i++) {
		for (int j = 0; j < M; j++) {
			C(i, j) = T(0);
			for (int k = 0; k < M; k++) {
				C(i, j) += A(i, k) * B(k, j);
			}
		}
	}
	return C;
}

template<typename T, int N>
SquareMatrix<T, N> operator*(SymmetricMatrix<T, N> const &A, SymmetricMatrix<T, N> const &B) {
	SquareMatrix<T, N> C;
	for (int i = 0; i < N; i++) {
		for (int j = 0; j < N; j++) {
			C(i, j) = T(0);
			for (int k = 0; k < N; k++) {
				C(i, j) += A(i, k) * B(k, j);
			}
		}
	}
	return C;
}

template<typename Type, int P, int Q, int M, int N>
Matrix<Type, P * M, Q * N> matrixKroneckerProduct(Matrix<Type, P, Q> const &A, Matrix<Type, M, N> const &B) {
	Matrix<Type, P * M, Q * N> C;
	for (int p = 0; p < P; p++) {
		int pQ = Q * p;
		for (int q = 0; q < Q; q++) {
			for (int m = 0; m < M; m++) {
				int mN = N * m;
				for (int n = 0; n < N; n++) {
					C(pQ + q, mN + n) = A(p, q) * B(n, m);
				}
			}
		}
	}
	return C;
}

template<typename T, int N>
void matrixQRDecomposition(SquareMatrix<T, N> const &A, SquareMatrix<T, N> &Q, SquareMatrix<T, N> &R) {
	T one = T(1);
	R = A;
	Q = SquareMatrix<T, N>::identity();
	for (int j = 0; j < N; j++) {
		for (int i = j + 1; i < N; i++) {
			SquareMatrix<T, N> G = SquareMatrix<T, N>::identity();
			T x = R(j, j), y = R(i, j);
			T inv = one / std::sqrt(x * x + y * y);
			T c = x * inv, s = -y * inv;
			G(j, j) = G(i, i) = c;
			G(i, j) = s;
			G(j, i) = -s;
			R = G * R;
			Q = G * Q;
		}
	}
	Q = matrixTranspose(Q);
}

template<typename T, int N>
void matrixLUDecompose(SquareMatrix<T, N> &A) {
	for (int n = 0; n < N - 1; n++) {
		T inv = T(1) / A(n, n);
		for (int k = n + 1; k < N; k++) {
			A(k, n) *= inv;
		}
		for (int j = n + 1; j < N; j++) {
			for (int k = n + 1; k < N; k++) {
				A(k, j) -= A(k, n) * A(n, j);
			}
		}
	}
}

template<typename T, int N, int M>
void matrixLUMultiply(SquareMatrix<T, N> const &LU, Matrix<T, N, M> &X) {
	for (int k = 0; k < N; k++) {
		for (int m = 0; m < M; m++) {
			X(k, m) *= LU(k, k);
		}
		for (int j = k + 1; j < N; j++) {
			for (int m = 0; m < M; m++) {
				X(k, m) += LU(k, j) * X(j, m);
			}
		}
	}
	for (int k = N - 1; k > 0; k--) {
		for (int j = 0; j < k; j++) {
			for (int m = 0; m < M; m++) {
				X(k, m) += LU(k, j) * X(j, m);
			}
		}
	}
}

template<typename T, int N>
void matrixLURecompose(SquareMatrix<T, N> &A) {
	for (int n = N - 2; n >= 0; n--) {
		for (int j = n + 1; j < N; j++) {
			for (int k = n + 1; k < N; k++) {
				A(k, j) += A(k, n) * A(n, j);
			}
		}
		T a = A(n, n);
		for (int k = n + 1; k < N; k++) {
			A(k, n) *= a;
		}
	}
}

template<typename T, int N>
struct DiagonalMatrix {
	constexpr T operator()(int c, int r) const {
		if (c == r) {
			return D[r];
		} else {
			return T(0);
		}
	}
	constexpr T& operator()(int c, int r) {
		assert(c == r);
		if (c != r) {
			throw std::runtime_error("Attempt to assign to diagnal matrix off diagonal\n");
		}
		return D[r];
	}
private:
	std::array<T, N> D;
};

template<typename T, int N>
constexpr T matrixDeterminant(DiagonalMatrix<T, N> const &A) {
	T det = T(1);
	for (int r = 0; r < N; r++) {
		det *= A(r, r);
	}
	return det;
}

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wmaybe-uninitialized"
template<typename T, int N>
constexpr DiagonalMatrix<T, N> matrixInverse(DiagonalMatrix<T, N> const &A) {
	DiagonalMatrix<T, N> Ainv;
	for (int r = 0; r < N; r++) {
		Ainv(r, r) = T(1) / A(r, r);
	}
	return Ainv;
}
#pragma GCC diagnostic pop

template<typename Type, int R, int C>
std::string toString(Matrix<Type, R, C> const &M) {
	auto formatEntry = [](Type const &value) {
		std::ostringstream oss;
		oss << std::scientific << std::setprecision(3) << value;
		return oss.str();
	};
	int maxLen = 0;
	for (int i = 0; i < R; ++i) {
		for (int j = 0; j < C; ++j) {
			std::string const entry = formatEntry(M(i, j));
			maxLen = std::max(maxLen, static_cast<int>(entry.size()));
		}
	}
	std::string line((C * (maxLen + 3) + 1), '-');
	line += "\n";
	std::string out = line;
	for (int i = 0; i < R; ++i) {
		for (int j = 0; j < C; ++j) {
			std::string cellStr = formatEntry(M(i, j));
			// pad on the left to align
			while (static_cast<int>(cellStr.size()) < maxLen) {
				cellStr = " " + cellStr;
			}
			out += "| " + cellStr + " ";
		}
		out += "|\n" + line;
	}

	return out;
}

template<typename Type, int R, int C>
std::string toMathematica(Matrix<Type, R, C> const &M) {
	using std::to_string;
	std::string out = "A =: {\n";
	for (int i = 0; i < R; i++) {
		out += "\t{";
		for (int j = 0; j < C; j++) {
			out += to_string(M(i, j));
			if (j + 1 < C)
				out += ",";
		}
		out += "}";
		if (i + 1 < R)
			out += ",";
		out += "\n";
	}
	out += "}";
	return out;
}

template<typename T, auto N, int M, typename Container>
inline constexpr auto operator*(Matrix<T, N, M> const &A, Container const &B) {
	Container C;
	for (int n = 0; n < N; n++) {
		C[n] = T(0);
		for (int m = 0; m < M; m++) {
			C[n] += A(n, m) * B[m];
		}
	}
	return C;
}

template<typename T, auto N, typename Container>
inline constexpr auto operator*(DiagonalMatrix<T, N> const &A, Container const &B) {
	Container C;
	for (int row = 0; row < N; row++) {
		int const &col = row;
		C[row] = A(row, col) * B[col];
	}
	return C;
}

