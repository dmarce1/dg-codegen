/*
 * Interpolate.hpp
 *
 *  Created on: Jan 20, 2025
 *      Author: dmarce1
 */

#include "Matrix.hpp"
#include "TriangularArray.hpp"

template<typename T, int N>
constexpr Math::SquareMatrix<T, N> initCoefficients() {
	using namespace Math;
	static constexpr T one = T(1), two = T(2);
	SquareMatrix<T, N> A;
	for (int n = 0; n < N; n++) {
		T const x = T(2 * n - N + 1) / T(N - 1);
		printf("%e\n", x);
		T xm = one;
		for (int m = 0; m < N; m++) {
			A[n, m] = xm;
			xm *= x;
		}
	}
	A = matrixInverse(A);
	for (int n = 0; n < N; n++) {
		T const norm = nSquared(nFactorial<T>(n)) * integerPower(two, n) / nFactorial<T>(2 * n);
		for (int m = 0; m < N; m++) {
			A[n, m] *= norm;
		}
	}
	return A;
}
;

template<typename T, int N>
Math::Vector<T, N> derivativeWeights(int k = 1) {
	using namespace Math;
	static const SquareMatrix<T, N> Wk = initCoefficients<T, N>();
	static const SquareMatrix<T, N> wT = Math::matrixTranspose(Wk);
	Math::Vector<T, N> w;
	for (int r = 0; r < N; r++) {
		w[r] = wT[r, k];
	}
	//std::cout << toString(w);
	return w;
}

