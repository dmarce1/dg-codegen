/*
 * Interpolate.hpp
 *
 *  Created on: Jan 20, 2025
 *      Author: dmarce1
 */

#include "Vector.hpp"

#include <iostream>
#include <valarray>

template<typename T>
struct TricubicSpline {
	static int constexpr N = 64;
	Math::SquareMatrix<int, N> Ainv;
	TricubicSpline() {
		Math::SquareMatrix<int, N> A;
		int row = 0;
		for (int xyz = 0; xyz < 8; xyz++) {
			int const x = (xyz >> 2) & 1;
			int const y = (xyz >> 1) & 1;
			int const z = (xyz >> 0) & 1;
			int col = 0;
			for (int i = 0; i < 4; i++) {
				for (int j = 0; j < 4; j++) {
					for (int k = 0; k < 4; k++) {
						int xi = (i == 0) ? 1 : x;
						int yj = (j == 0) ? 1 : y;
						int zk = (k == 0) ? 1 : z;
						int xim1 = i * ((i == 0) ? 0 : ((i == 1) ? 1 : x));
						int yjm1 = j * ((j == 0) ? 0 : ((j == 1) ? 1 : y));
						int zkm1 = k * ((k == 0) ? 0 : ((k == 1) ? 1 : z));
						// f
						A[row + 0, col] = xi * yj * zk;
						A[row + 1, col] = xim1 * yj * zk;
						A[row + 2, col] = xi * yjm1 * zk;
						A[row + 3, col] = xi * yj * zkm1;
						A[row + 4, col] = xim1 * yjm1 * zk;
						A[row + 5, col] = xim1 * yj * zkm1;
						A[row + 6, col] = xi * yjm1 * zkm1;
						A[row + 7, col] = xim1 * yjm1 * zkm1;
						col++;

					}
				}
			}
			row += 8;
		}
		Ainv = matrixInverse(A);
		std::cout << toMathematica(Ainv) << "\n";

	}

};
