/******************************************************************************
 Copyright (C) 2024  Dominic C. Marcello
 *******************************************************************************/

#include <hpx/hpx_init.hpp>
#include "HLLC.hpp"
#include "Integrate.hpp"
#include "LegendreP.hpp"
#include "Polynomial.hpp"
#include "Real.hpp"
#include "Vector.hpp"
#include "Complex.hpp"
#include <unordered_map>
#include <numeric>
#include <stack>
#include <type_traits>
#include <debug/vector>
#include <span>
using namespace std;

Math::Polynomial<Real> lagrangePolynomial(int k, int j, std::function<Real(int)> xn = nullptr) {
	if (xn == nullptr) {
		xn = [k](int n) {
			return Real(2 * n - k) / Real(k + 2);
		};
	}
	Math::Polynomial<Real> p;
	p[0] = Real(1);
	Real xj = xn(j);
	for (int m = 0; m <= k; m++) {
		if (m != j) {
			Real xm = xn(m);
			Real c0 = Real(1) / (xj - xm);
			Math::Polynomial<Real> q;
			q[0] = -xm * c0;
			q[1] = c0;
			p *= q;
		}
	}
	return p;
}

void compute();

int hpx_main(int argc, char *argv[]) {
	compute();
	/*Matrix<Real, N, N> M;
	 Matrix<Real, N, N> D;
	 srand(42);
	 for (int n = 0; n < N; n++) {
	 for (int m = 0; m < N; m++) {
	 M[n, m] = Real(2 * (rand() & 1) - 1);
	 }
	 }
	 M[1, 0] = Real(0);
	 M[1, 2] = Real(0);
	 std::cout << to_string(M);
	 matrixRowReduction(M, D);
	 std::cout << to_string(M);
	 std::cout << to_string(D);*/
	return hpx::local::finalize();
}

int main(int argc, char *argv[]) {
	auto rc = hpx::init(argc, argv);
	return rc;
}
