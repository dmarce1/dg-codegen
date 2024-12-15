/******************************************************************************
 Copyright (C) 2024  Dominic C. Marcello
 *******************************************************************************/

#include <hpx/hpx_init.hpp>
#include "Matrix.hpp"
#include "Polynomial.hpp"
#include "Real.hpp"
#include "Vector.hpp"
#include "Complex.hpp"

#include <type_traits>
using namespace std;

int hpx_main(int argc, char *argv[]) {
	constexpr int N = 10;
	Math::Polynomial<double, N> P;
	P[0] = 1;
	P[1] = 1;
	P[2] = 1;
	P[3] = 1;
	P[4] = 1;
	P[5] = 1;
	P[6] = 1;
	P[7] = 1;
	P[8] = 1;
	P[9] = 1;
	P[10] = 1;
	std::cout << "P:   " << P.toString() << "\n";
	auto roots = Math::polynomialFindAllRoots(P);
	auto complex2str = [](Math::Complex<double> z) {
		std::string str;
		str += std::to_string(z.real());
		if (z.imaginary() > 0.0) {
			str += " + I" + std::to_string(abs(z.imaginary()));
		} else if (z.imaginary() < 0.0) {
			str += " - I" + std::to_string(abs(z.imaginary()));
		}
		return str;
	};
	for (int n = 0; n < N; n++) {
		printf("%i %s\n", n, complex2str(roots[n]).c_str());
	}
	return hpx::local::finalize();
}

int main(int argc, char *argv[]) {
	auto rc = hpx::init(argc, argv);
	return rc;
}
