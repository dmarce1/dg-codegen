#pragma once

#include <algorithm>
#include <array>
#include <cassert>
#include <iostream>
#include <string>
#include <stacktrace>
#include <type_traits>

namespace Math {

template<typename Type, int RowCount, int ColumnCount>
struct Matrix;

template<typename Type, int Ndim>
using SquareMatrix = Matrix<Type, Ndim, Ndim>;

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
	Matrix() :
			values() {
	}
	Matrix(std::array<std::array<Type, ColumnCount>, RowCount> const &initList) :
			values() {
		for (int n = 0; n < RowCount; n++) {
			for (int m = 0; m < ColumnCount; m++) {
				operator[](n, m) = initList[n][m];
			}
		}
	}
	Matrix(std::initializer_list<Type> initList) {
		int i = 0;
		for (auto const &v : initList) {
			values[i++] = v;
		}
	}
	Matrix(std::initializer_list<std::initializer_list<Type>> init) {
		int i = 0;
		for (const auto &row : init) {
			for (const auto &val : row) {
				values[i++] = val;
			}
		}
	}
	Matrix(Type const &init) :
			values() {
		for (int n = 0; n != size(); n++) {
			values[n] = init;
		}
	}
	Matrix(Matrix const &other) :
			values() {
		values = other.values;
	}
	Matrix(Matrix &&other) :
			values() {
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
	Type const& operator[](int n, int m) const {
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
	static constexpr std::enable_if<RowCount == ColumnCount, Matrix>::type identity() {
		Matrix I;
		for (int n = 0; n < RowCount; n++) {
			for (int m = 0; m < RowCount; m++) {
				I[n, m] = Type(n == m);
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
			values() {
	}
	Matrix(Type const &initValue) :
			values() {
		std::fill(begin(), end(), initValue);
	}
	Matrix(std::initializer_list<Type> initList) :
			values() {
		std::copy(initList.begin(), initList.end(), begin());
	}
	Matrix(Matrix const &other) :
			values() {
		values = other.values;
	}
	Matrix(Matrix &&other) :
			values() {
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
	Type const& operator[](int n, int m) const {
		assert(m == 0);
		return operator[](n);
	}
	Type& operator[](int n) {
		return values[n];
	}
	Type const& operator[](int n) const {
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
	static constexpr Matrix zero() {
		Matrix Z;
		std::fill(Z.begin(), Z.end(), Type(0));
		return Z;
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
	Type const& operator[](int, int) const {
		return value;
	}
	Type& operator[](int n) {
		return value;
	}
	Type const& operator[](int n) const {
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
	auto begin() const {
		return &value;
	}
	auto end() const {
		return &value + 1;
	}
	auto begin() {
		return &value;
	}
	auto end() {
		return &value + 1;
	}
	friend Matrix operator*(Type const &a, Matrix const &B) {
		return B * a;
	}
	static constexpr Matrix identity() {
		return Matrix<Type, 1, 1> { Type(1) };
	}
	static constexpr Matrix zero() {
		return Matrix<Type, 1, 1> { Type(0) };
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
		*this = SymmetricMatrix(matrixSymmetrize(SquareMatrix<T, D>(init)));
	}
	T& operator[](int n, int k) {
		return ((base_type&) (*this))[index(n, k)];
	}
	T const& operator[](int n, int k) const {
		return ((base_type&) (*this))[index(n, k)];
	}
	operator SquareMatrix<T, D>() const {
		SquareMatrix<T, D> A;
		for (int p = 0; p < D; p++) {
			for (int q = p + 1; q < D; q++) {
				A[p, q] = A[q, p] = (*this).base_type::operator[](((p * p + p) >> 1) + q);
			}
		}
		return A;
	}
	static constexpr SymmetricMatrix identity() {
		SymmetricMatrix<T, D> identity;
		for (int n = 0; n < D; n++) {
			for (int m = 0; m <= n; m++) {
				identity[n, m] = T(n == m);
			}
		}
		return identity;
	}
	static constexpr SymmetricMatrix zero() {
		SymmetricMatrix Z;
		std::fill(Z.begin(), Z.end(), T(0));
		return Z;

	}
	SymmetricMatrix& operator+=(SymmetricMatrix const &A) {
		(*this) = (*this) + A;
		return *this;
	}
	SymmetricMatrix& operator-=(SymmetricMatrix const &A) {
		(*this) = (*this) - A;
		return *this;
	}
	SymmetricMatrix& operator*=(T const &a) {
		(*this) = (*this) * a;
		return *this;
	}
	friend SymmetricMatrix operator+(SymmetricMatrix const &A, SymmetricMatrix const &B) {
		SymmetricMatrix C;
		for (size_t n = 0; n < A.size(); n++) {
			C.base_type::operator[](n) = A.base_type::operator[](n) + B.base_type::operator[](n);
		}
		return C;
	}
	friend SymmetricMatrix operator-(SymmetricMatrix const &A, SymmetricMatrix const &B) {
		SymmetricMatrix C;
		for (size_t n = 0; n < A.size(); n++) {
			C.base_type::operator[](n) = A.base_type::operator[](n) - B.base_type::operator[](n);
		}
		return C;
	}
	friend SymmetricMatrix operator*(SymmetricMatrix const &A, T b) {
		SymmetricMatrix C;
		for (size_t n = 0; n < A.size(); n++) {
			C.base_type::operator[](n) = A.base_type::operator[](n) * b;
		}
		return C;
	}
	friend SymmetricMatrix operator*(T a, SymmetricMatrix const &B) {
		SymmetricMatrix C;
		for (size_t n = 0; n < B.size(); n++) {
			C.base_type::operator[](n) = B.base_type::operator[](n) * a;
		}
		return C;
	}
	friend SymmetricMatrix operator/(SymmetricMatrix const &A, T b) {
		return A * (T(1) / b);
	}
private:
	static constexpr int index(int n, int k) {
		int const p = std::max(n, k);
		int const q = std::min(n, k);
		int const i = ((p * p + p) >> 1) + q;
		return i;
	}
};

template<typename T, int D>
SymmetricMatrix<T, D> matrixSymmetrize(SquareMatrix<T, D> const &A) {
	SymmetricMatrix<T, D> B;
	for (int n = 0; n < D; n++) {
		B[n, n] = A[n, n];
		for (int m = 0; m < n; m++) {
			B[n, m] = T(0.5) * (A[n, m] + A[m, n]);
		}
	}
	return B;
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
Matrix<Type, RowCount, 1> matrixColumn(Matrix<Type, RowCount, ColumnCount> const &A, int c) {
	Matrix<Type, RowCount, 1> column;
	for (int r = 0; r < RowCount; r++) {
		column[r, 0] = A[r, c];
	}
	return column;
}

template<typename Type, int N>
Matrix<Type, N, 1> matrixColumn(SymmetricMatrix<Type, N> const &A, int c) {
	Matrix<Type, N, 1> column;
	for (int r = 0; r < N; r++) {
		column[r, 0] = A[r, c];
	}
	return column;
}

template<typename T, int N>
T matrixInverseAndDeterminant(SquareMatrix<T, N> &A) {
	T constexpr zero(0), one(1);
	T matrixDeterminant = one;
	auto D = SquareMatrix<T, N>::identity();
	for (int i = 0; i < N; ++i) {
		T pivot = A[i, i];
		if (pivot == zero) {
			int k = i;
			do {
				k++;
				if (k >= N) {
					printf("(%i >= %i)\n", k, N);
					abort();
					return zero;
				}
				pivot = A[k, i];
			} while (pivot == zero);
			bool isSingular = true;
			for (int j = 0; j < N; ++j) {
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
		for (int j = 0; j < N; ++j) {
			A[i, j] *= iPivot;
			D[i, j] *= iPivot;
		}
		matrixDeterminant *= pivot;
		for (int k = 0; k < N; ++k) {
			if (k != i) {
				T const factor = A[k, i];
				for (int j = 0; j < N; ++j) {
					A[k, j] -= factor * A[i, j];
					D[k, j] -= factor * D[i, j];
				}
			}
		}
	}
	A = D;
	return matrixDeterminant;
}

template<typename Type, int Ndim>
Type matrixTrace(SquareMatrix<Type, Ndim> const &A) {
	Type result = Type(0);
	for (int c = 0; c < Ndim; c++) {
		result += A[c, c];
	}
	return result;
}

template<typename Type, int RowCount, int ColumnCount>
Matrix<Type, ColumnCount, RowCount> matrixTranspose(Matrix<Type, RowCount, ColumnCount> const &B) {
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

template<typename T, int N>
SymmetricMatrix<T, N - 1> subMatrix(SymmetricMatrix<T, N> const &A, int d) {
	SymmetricMatrix<T, N - 1> subMatrix;
	for (int i = 0; i < d; i++) {
		for (int j = 0; j <= i; j++) {
			subMatrix[i, j] = A[i, j];
		}
	}
	for (int i = d; i < N - 1; i++) {
		for (int j = 0; j < d; j++) {
			subMatrix[i, j] = A[i + 1, j];
		}
		for (int j = d; j <= i; j++) {
			subMatrix[i, j] = A[i + 1, j + 1];
		}
	}
	return subMatrix;
}

template<typename T, int N>
T matrixCofactor(SquareMatrix<T, N> const &A, int r, int c) {
	if constexpr (N > 1) {
		auto const sgn = ((r + c) & 1) ? -T(1) : T(1);
		return sgn * matrixDeterminant(subMatrix(A, r, c));
	} else {
		return A[0, 0];
	}
}

template<typename T, int N>
T matrixCofactor(SymmetricMatrix<T, N> const &A, int r, int c) {
	if constexpr (N > 1) {
		if (r == c) {
			return matrixDeterminant(subMatrix(A, r));
		} else {
			auto const sgn = ((r + c) & 1) ? -T(1) : T(1);
			return sgn * matrixDeterminant(subMatrix(SquareMatrix<T, N>(A), r, c));
		}
	} else {
		return A[0, 0];
	}
}

template<typename T, int N>
SquareMatrix<T, N> matrixCofactor(SquareMatrix<T, N> const &A) {
	SquareMatrix<T, N> aCofactor;
	for (int r = 0; r < N; r++) {
		for (int c = 0; c < N; c++) {
			aCofactor[r, c] = matrixCofactor(A, r, c);
		}
	}
	return aCofactor;
}

template<typename T, int N>
SymmetricMatrix<T, N> matrixCofactor(SymmetricMatrix<T, N> const &A) {
	SymmetricMatrix<T, N> aCofactor;
	for (int i = 0; i < N; i++) {
		for (int j = 0; j <= i; j++) {
			aCofactor[i, j] = matrixCofactor(A, i, j);
		}
	}
	return aCofactor;
}

template<typename T, int N>
T matrixDeterminant(SymmetricMatrix<T, N> A) {
	if constexpr (N > 1) {
		T result = matrixCofactor(A, 0, 0);
		for (int c = 1; c < N; c++) {
			result += matrixCofactor(A, 0, c);
		}
		return result;
	} else {
		return A[0, 0];
	}
}

template<typename T, int N>
T matrixDeterminant(SquareMatrix<T, N> A) {
	if constexpr (N > 1) {
		T result = matrixCofactor(A, 0, 0);
		for (int c = 1; c < N; c++) {
			result += matrixCofactor(A, 0, c);
		}
		return result;
	} else {
		return A[0, 0];
	}
}

template<typename T, int N>
SquareMatrix<T, N> matrixAdjoint(SquareMatrix<T, N> const &A) {
	return matrixTranspose(matrixCofactor(A));
}

template<typename T, int N>
SymmetricMatrix<T, N> matrixAdjoint(SymmetricMatrix<T, N> const &A) {
	return matrixCofactor(A);
}

template<typename T, int N>
auto matrixInverse(SymmetricMatrix<T, N> const A) {
	return matrixAdjoint(A) / matrixDeterminant(A);
}
template<typename T, int N>
auto matrixInverse(SquareMatrix<T, N> const A) {
	return matrixAdjoint(A) / matrixDeterminant(A);
}

template<typename T, int N, int M>
SquareMatrix<T, N> operator*(SymmetricMatrix<T, N> const &A, Matrix<T, N, M> const &B) {
	SquareMatrix<T, N> C;
	for (int i = 0; i < N; i++) {
		for (int j = 0; j < M; j++) {
			C[i, j] = T(0);
			for (int m = 0; m < N; m++) {
				C[i, j] += A[i, m] * B[m, j];
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
			C[i, j] = T(0);
			for (int m = 0; m < M; m++) {
				C[i, j] += A[i, m] * B[m, j];
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
			C[i, j] = T(0);
			for (int m = 0; m < N; m++) {
				C[i, j] += A[i, m] * B[m, j];
			}
		}
	}
	return C;
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
	Q = SquareMatrix<T, N>::identity();
	for (int j = 0; j < N; j++) {
		for (int i = j + 1; i < N; i++) {
			SquareMatrix<T, N> G = SquareMatrix<T, N>::identity();
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

template<typename T, int N>
void matrixLUDecompose(SquareMatrix<T, N> &A) {
	for (int n = 0; n < N - 1; n++) {
		T const inv = T(1) / A[n, n];
		for (int k = n + 1; k < N; k++) {
			A[k, n] *= inv;
		}
		for (int j = n + 1; j < N; j++) {
			for (int k = n + 1; k < N; k++) {
				A[k, j] -= A[k, n] * A[n, j];
			}
		}
	}
}

template<typename T, int N, int M>
void matrixLUMultiply(SquareMatrix<T, N> const &LU, Matrix<T, N, M> &X) {
	Matrix<T, N, M> Y;
	for (int k = 0; k < N; k++) {
		for (int m = 0; m < M; m++) {
			X[k, m] *= LU[k, k];
		}
		for (int j = k + 1; j < N; j++) {
			for (int m = 0; m < M; m++) {
				X[k, m] += LU[k, j] * X[j, m];
			}
		}
	}
	for (int k = N - 1; k > 0; k--) {
		for (int j = 0; j < k; j++) {
			for (int m = 0; m < M; m++) {
				X[k, m] += LU[k, j] * X[j, m];
			}
		}
	}
}

template<typename T, int N>
void matrixLURecompose(SquareMatrix<T, N> &A) {
	for (int n = N - 2; n >= 0; n--) {
		for (int j = n + 1; j < N; j++) {
			for (int k = n + 1; k < N; k++) {
				A[k, j] += A[k, n] * A[n, j];
			}
		}
		T const a = A[n, n];
		for (int k = n + 1; k < N; k++) {
			A[k, n] *= a;
		}
	}
}

/***
 * 1   0   0 | u00 u01 u02    |u00     u01              u02
 * l10 1   0 |   0 u11 u12    |u00*l10 u01*l10+u11      u02*l10+u12
 * l20 l21 1 |   0   0 u22    |u00*l20 u01*l20+u11*l21  u02*l20+u12*l21+u22
 *
 *
 */
/*template<typename T, int R>
 SquareMatrix<T, R> matrixLUDecomposition(SquareMatrix<T, R> const &A, SquareMatrix<T, R> &L, SquareMatrix<T, R> &U) {
 const T zero(0);
 SquareMatrix<T, R> P = SquareMatrix<T, R>::identity();
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
 }*/

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

template<typename Type, int RowCount, int ColumnCount>
std::string toMathematica(Matrix<Type, RowCount, ColumnCount> const &M) {
	using std::to_string;
	std::string list;
	list += "A =: {\n";
	for (int row = 0; row < RowCount; row++) {
		list += "\t{";
		for (int column = 0; column < ColumnCount; column++) {
			list += to_string(M[row, column]);
			if (column != ColumnCount - 1) {
				list += ",";
			}
		}
		list += " }";
		if (row != RowCount - 1) {
			list += ",";
		}
	}
	list += "}";
	return list;
}

template<typename T>
SymmetricMatrix<T, 4> minkowskiMetric() {
	static constexpr T zero = T(0), one = T(1);
	using return_type = SymmetricMatrix<T, DIM4>;
	return_type nu = return_type(zero);
	nu[0, 0] = -one;
	nu[1, 1] = +one;
	nu[2, 2] = +one;
	nu[3, 3] = +one;
	return nu;
}

}
