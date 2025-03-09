/*
 * Einstein.cpp
 *
 *  Created on: Feb 25, 2025
 *      Author: dmarce1
 */

#include <algorithm>
#include <cstdio>
#include <cstring>
#include <iostream>
#include <string>
#include <vector>

template<typename T>
struct ArgConverter {
	static T convert(const T &arg) {
		return arg;
	}
};

template<>
struct ArgConverter<std::string> {
	static const char* convert(const std::string &arg) {
		return arg.c_str();
	}
};

template<typename ... Args>
int printf_wrapper(const std::string &format, Args &&... args) {
	return std::printf(format.c_str());
}

struct CodeGen {
	static constexpr char const *const tabString = "   ";
	CodeGen() {
		tabs = 0;
		nl = false;
	}
	template<typename ... Args>
	void print(std::string const &fmt, Args &&...args) {
		std::string line;
		char *buffer;
		asprintf(&buffer, fmt.c_str(), ArgConverter<typename std::decay<Args>::type>::convert(std::forward<Args>(args))...);

		for (int n = 0; n < std::strlen(buffer); n++) {
			code.push_back(buffer[n]);
			if (buffer[n] == '\n') {
				for (int t = 0; t < tabs; t++) {
					code += tabString;
				}
			}
		}
		free(buffer);
	}
	void indent() {
		tabs++;
	}
	void dedent() {
		tabs--;
	}
	std::string get() const {
		return code;
	}
private:
	int tabs;
	bool nl;
	std::string code;
};

struct TensorVar {
	TensorVar(std::string const &name_, std::string const &sig, bool symmetric_ = false) {
		name = name_;
		order = sig.size();
		symmetric = symmetric_;
		for (int i = 0; i < sig.size(); i++) {
			if (std::isupper(sig[i])) {
				contravar.push_back(true);
			} else {
				contravar.push_back(false);
			}
		}
	}
	std::string get(std::vector<int> I) const {
		std::string element = name;
		element += "_";
		if (symmetric) {
			int const i0 = I.size() - 1;
			int const i1 = I.size() - 2;
			if (I[i0] < I[i1]) {
				std::swap(I[i1], I[i0]);
			}
		}
		for (int i = 0; i < I.size(); i++) {
			if (contravar[i]) {
				element.push_back('X' + I[i]);
			} else {
				element.push_back('x' + I[i]);
			}
		}
		return element;
	}
	std::string get(auto ...I) const {
		std::vector<int> V;
		((V.push_back(I)),...);
		return get(V);
	}
	int getOrder() const {
		return order;
	}
	int isSymmetric() const {
		return symmetric;
	}
	const char* getName() const {
		return name.c_str();
	}
private:
	std::string name;
	int order;
	bool symmetric;
	std::vector<bool> contravar;
};

//void tensorContract(CodeGen &gen, TensorVar const &A, TensorVar const &B, TensorVar const &C, int bi, int ci, std::string prefix = "T const ") {
//	std::vector<int> aI(A.size(), 0);
//	std::vector<int> bI(B.size(), 0);
//	std::vector<int> cI(C.size(), 0);
//	bool done = false;
//	while (!done) {
//		if (!(A.isSymmetric() && aI[aI.size() - 2] < aI.back())) {
//			gen.print("%s %s = ", prefix.c_str(), A.get(aI));
//			for (int dim = 0; dim < NDIM; dim++) {
//				for( int i = 0; i < B.size(); i++) {
//					if( i < bi) {
//						B[i] = A[i];
//					}  else if( i > bi) {
//						B[i] = A[i - 1];
//					} else {
//						B[i] = dim;
//					}
//				}
//				for( int i = 0; i < C.size(); i++) {
//					if( i < Ci) {
//						C[i] = A[i + B.size() - 1];
//					}  else if( i > ci) {
//						C[i] = A[i + B.size() - 2];
//					} else {
//						C[i] = dim;
//					}
//				}
//				gen.print("%s * %s%s ", B.get(bI), C.get(cI), dim < NDIM - 1 ? " + " : ";\n");
//			}
//		}
//		int dim = 0;
//		while (++I[dim] == NDIM) {
//			I[dim++] = 0;
//			if (dim == A.size()) {
//				done = true;
//			}
//		}
//	}
//}

int main(int, char*[]) {
	TensorVar alpha0("alpha", "");
	TensorVar beta1("beta", "U");
	TensorVar gamma2("gamma", "ll", true);
	TensorVar A1("A", "l");
	TensorVar B2("B", "lU");
	TensorVar D3("D", "lll", true);
	TensorVar K2("K", "ll", true);
	TensorVar Theta0("Theta", "");
	TensorVar Z1("Z", "l");
	constexpr int NFIELDS = 50;
	TensorVar D3UU("D", "lUU", true);
	TensorVar K2U("K", "Ul");
	TensorVar alpha0inv("ahpla", "");
	TensorVar gamma0("gamma", "");
	TensorVar gamma0inv("igamma", "");
	TensorVar gamma2inv("gamma", "UU", true);
	TensorVar B0("B", "");
	TensorVar K0("K", "");
	TensorVar D1("D", "l");
	TensorVar E1("E", "l");
	TensorVar D1U("D", "U");
	TensorVar E1U("E", "U");
	TensorVar Z1U("Z", "U");
	TensorVar V1("V", "l");
	TensorVar Gamma3l("Gamma", "lll", true);
	TensorVar Gamma3U("Gamma", "Ull", true);
	TensorVar lambda3U("lambda", "Ull");
	TensorVar lambda3("lambda", "lll");
	TensorVar Q0("Q", "");
	TensorVar Q1("Q", "l");
	TensorVar Q2("Q", "ll");
	int I[6] = { 0, 1, 2, 0, 0, 1 };
	int J[6] = { 0, 1, 2, 1, 2, 2 };
	char const *plus[NDIM] = { " + ", " + ", "" };
	CodeGen gen;
	gen.print("\n#pragma once\n\n");
	gen.print("#include <array>\n\n");
	gen.print("#include <stdexcept>\n\n");
	gen.print("template<typename T>\n", NFIELDS);
	gen.print("struct Spacetime { \n");
	gen.indent();
	gen.print("\nSpacetime() :");
	gen.indent();
	gen.print("\n");
	int n = 0;
	gen.print("%s(U[%i]), ", alpha0.get(), n++);
	gen.print("\n");
	for (int i = 0; i < NDIM; i++) {
		gen.print("%s(U[%i]), ", beta1.get(i), n++);
	}
	gen.print("\n");
	for (int i = 0; i < NDIM; i++) {
		for (int j = i; j < NDIM; j++) {
			gen.print("%s(U[%i]), ", gamma2.get(i, j), n++);
		}
		gen.print("\n");
	}
	for (int i = 0; i < NDIM; i++) {
		for (int j = i; j < NDIM; j++) {
			gen.print("%s(U[%i]), ", K2.get(i, j), n++);
		}
		gen.print("\n");
	}
	for (int i = 0; i < NDIM; i++) {
		gen.print("%s(U[%i]), ", A1.get(i), n++);
	}
	gen.print("\n");
	for (int i = 0; i < NDIM; i++) {
		for (int j = 0; j < NDIM; j++) {
			gen.print("%s(U[%i]), ", B2.get(i, j), n++);
		}
		gen.print("\n");
	}
	for (int k = 0; k < NDIM; k++) {
		for (int i = 0; i < NDIM; i++) {
			for (int j = i; j < NDIM; j++) {
				gen.print("%s(U[%i]), ", D3.get(k, i, j), n++);
			}
			gen.print("\n");
		}
	}
	gen.print("%s(U[%i]), ", Theta0.get(), n++);
	gen.print("\n");
	for (int i = 0; i < NDIM - 1; i++) {
		gen.print("%s(U[%i]), ", Z1.get(i), n++);
	}
	gen.print("%s(U[%i])", Z1.get(NDIM - 1), n++);
	gen.dedent();
	gen.print("\n{}\n");
	gen.print("\n");
	gen.print("Spacetime flux(int dim) const {");
	gen.indent();
	gen.print("\nSpacetime F;\n");
	gen.print("\n");

	auto const compInvGamma = [&]() {
		for (int ij = 0; ij < 6; ij++) {
			int i = I[ij];
			int j = J[ij];
			int n1 = i == 0 ? 1 : 0;
			int n2 = i == 2 ? 1 : 2;
			int m1 = j == 0 ? 1 : 0;
			int m2 = j == 2 ? 1 : 2;
			if ((i + j) % 2 == 1) {
				std::swap(n1, n2);
				std::swap(m1, m2);
			}
			gen.print("T %s = ", gamma2inv.get(i, j));
			std::string arg1 = gamma2.get(n1, m1);
			std::string arg2 = gamma2.get(n2, m2);
			std::string arg3 = gamma2.get(n1, m2);
			std::string arg4 = gamma2.get(n2, m1);
			gen.print("%s * %s - %s * %s", arg1.c_str(), arg2.c_str(), arg3.c_str(), arg4.c_str());
			gen.print(";\n");
		}
		gen.print("\n");
		gen.print("T const  %s = ", gamma0.get());
		for (int i = 0; i < NDIM; i++) {
			gen.print("%s * %s%s", gamma2.get(0, i), gamma2inv.get(0, i), plus[i]);
		}

		gen.print(";\n\n");
		gen.print("T const  %s = T(1) / %s;\n", gamma0inv.get(), gamma0.get());

		gen.print("\n");
		for (int ij = 0; ij < 6; ij++) {
			int i = I[ij];
			int j = J[ij];
			gen.print("%s *= %s;\n", gamma2inv.get(i, j), gamma0inv.get());
		}

	};
	compInvGamma();
	gen.print("\n");
	for (int n = 0; n < NDIM; n++) {
		gen.print("T const  %s = ", D1.get(n));
		for (int k = 0; k < NDIM; k++) {
			gen.print("%s * %s + ", D3.get(n, k, k), gamma2inv.get(k, k));
		}
		gen.print("T(2) * (");
		for (int k = 0; k < NDIM; k++) {
			for (int m = k + 1; m < NDIM; m++) {
				gen.print("%s * %s", D3.get(n, k, m), gamma2inv.get(m, k));
				if (k + m < 2 * NDIM - 3) {
					gen.print(" + ");
				} else {
					gen.print(")");
				}
			}
		}
		gen.print(";\n");
	}

	gen.print("\n");
	for (int n = 0; n < NDIM; n++) {
		gen.print("T const  %s = ", E1.get(n));
		for (int k = 0; k < NDIM; k++) {
			for (int m = 0; m < NDIM; m++) {
				gen.print("%s * %s", D3.get(k, n, m), gamma2inv.get(m, k));
				if (k + m < 2 * NDIM - 2) {
					gen.print(" + ");
				}
			}
		}
		gen.print(";\n");
	}
	gen.print("\n");
	for (int k = 0; k < NDIM; k++) {
		gen.print("T const  %s = %s - %s - %s;\n", V1.get(k), D1.get(k), E1.get(k), Z1.get(k));
	}
	gen.print("\nT const  %s = ", K0.get());
	for (int n = 0; n < NDIM; n++) {
		gen.print("%s * %s + ", gamma2inv.get(n, n), K2.get(n, n));
	}
	gen.print("T(2) * (");
	for (int n = 0; n < NDIM; n++) {
		for (int k = n + 1; k < NDIM; k++) {
			gen.print("%s * %s", gamma2inv.get(n, k), K2.get(k, n));
			if (k + n < 2 * NDIM - 3) {
				gen.print(" + ");
			}
		}
	}
	gen.print(");\n");

	gen.print("\n");
	gen.print("T const  %s = T(1) / %s;\n", alpha0inv.get(), alpha0.get());

	gen.print("\n");
	gen.print("T const  %s = %s - T(2) * %s;\n", Q0.get(), K0.get(), Theta0.get());

	gen.print("\n");
	for (int k = 0; k < NDIM; k++) {
		gen.print("T const  %s = %s * (%s - %s + T(2) * %s);\n", Q1.get(k), alpha0.get(), A1.get(k), D1.get(k), V1.get(k));
	}

	gen.print("\n");
	for (int i = 0; i < NDIM; i++) {
		gen.print("T const  %s = %s - %s * %s;\n", Q2.get(i, i), K2.get(i, i), alpha0inv.get(), B2.get(i, i));
	}
	for (int i = 0; i < NDIM; i++) {
		for (int j = i + 1; j < NDIM; j++) {
			std::string tmp = B2.get(i, j);
			gen.print("T const  %s = %s - T(0.5) * %s * (%s + %s);\n", Q2.get(i, j), K2.get(i, j), alpha0inv.get(), tmp.c_str(), B2.get(j, i));
		}
	}
	gen.print("\n");
	for (int i = 0; i < NDIM; i++) {
		for (int j = 0; j < NDIM; j++) {
			gen.print("T const  %s = ", K2U.get(i, j));
			for (int m = 0; m < NDIM; m++) {
				gen.print("%s * %s", gamma2inv.get(i, m), K2.get(m, j));
				if (m < NDIM - 1) {
					gen.print(" + ");
				}
			}
			gen.print(";\n");
		}
	}
	gen.print("\n");
	gen.print("switch( dim ) {\n");
	for (int dim = 0; dim < NDIM; dim++) {
		gen.print("\n");
		gen.print("case %i: {", dim);
		gen.indent();
		gen.print("\n");
		gen.print("\n");
		for (int i = 0; i < NDIM; i++) {
			for (int j = i; j < NDIM; j++) {
				gen.print("T const  %s = ", D3UU.get(dim, i, j));
				for (int m = 0; m < NDIM; m++) {
					gen.print("%s * (", gamma2inv.get(i, m));
					for (int r = 0; r < NDIM; r++) {
						gen.print("%s * %s%s", D3.get(dim, m, r), gamma2inv.get(r, j), plus[r]);
					}
					gen.print(")%s", plus[m]);
				}
				gen.print(";\n");
			}
		}

		gen.print("\n");
		for (int k = 0; k < NDIM; k++) {
			for (int i = 0; i < NDIM; i++) {
				for (int j = i; j < NDIM; j++) {
					char c;
					if (i == k || j == k) {
						c = ' ';
					} else {
						c = '&';
					}
					gen.print("T const %c%s = %s", c, lambda3.get(k, i, j), D3.get(k, i, j));
					if (i == k) {
						if (j == k) {
							gen.print(" + %s + %s - T(2) * (%s + %s)", A1.get(j), D1.get(j), E1.get(j), Z1.get(j));
						} else {
							gen.print(" + T(0.5) * (%s + %s - T(2) * (%s + %s))", A1.get(j), D1.get(j), E1.get(j), Z1.get(j));
						}
					} else {
						if (j == k) {
							gen.print(" + T(0.5) * (%s + %s - T(2) * (%s + %s))", A1.get(i), D1.get(i), E1.get(i), Z1.get(i));
						}
					}
					gen.print(";\n");
				}
			}
		}
		for (int i = 0; i < NDIM; i++) {
			for (int j = i; j < NDIM; j++) {
				gen.print("T const  %s = ", lambda3U.get(dim, i, j));
				for (int m = 0; m < NDIM; m++) {
					gen.print("%s * %s%s", gamma2inv.get(dim, m), lambda3.get(m, i, j), plus[m]);
				}
				gen.print(";\n");
			}
		}
		//alpha
		gen.print("\n");
		gen.print("F.%s = ", alpha0.get());
		// beta
		for (int i = 0; i < NDIM; i++) {
			gen.print("F.%s = ", beta1.get(i));
		}
		for (int i = 0; i < NDIM; i++) {
			for (int j = i; j < NDIM; j++) {
				gen.print("F.%s = ", gamma2.get(i, j));
			}
		}
		gen.print("T(0);\n");

		// K
		gen.print("\n");
		for (int i = 0; i < NDIM; i++) {
			for (int j = i; j < NDIM; j++) {
				gen.print("F.%s = -%s * %s + %s * %s;\n", K2.get(i, j), beta1.get(dim), K2.get(i, j), alpha0.get(), lambda3U.get(dim, i, j));
			}
		}
		// A
		gen.print("\n");
		for (int k = 0; k < NDIM; k++) {
			gen.print("F.%s = -%s * %s", A1.get(k), beta1.get(dim), A1.get(k));
			if (k == dim) {
				gen.print(" + %s * %s", alpha0.get(), Q0.get());
			}
			gen.print(";\n");
		}
		// B
		gen.print("\n");
		for (int k = 0; k < NDIM; k++) {
			for (int i = 0; i < NDIM; i++) {
				gen.print("F.%s = -%s * %s", B2.get(k, i), beta1.get(dim), B2.get(k, i));
				if (k == dim) {
					gen.print(" + %s * (", alpha0.get());
					for (int m = 0; m < NDIM; m++) {
						gen.print("%s * %s%s", gamma2inv.get(i, m), Q1.get(m), plus[m]);
					}
					gen.print(");\n");
				} else {
					gen.print(";\n");
				}
			}
		}
		// D
		gen.print("\n");
		for (int k = 0; k < NDIM; k++) {
			for (int i = 0; i < NDIM; i++) {
				for (int j = i; j < NDIM; j++) {
					gen.print("F.%s = -%s * %s", D3.get(k, i, j), beta1.get(dim), D3.get(k, i, j));
					if (k == dim) {
						gen.print(" + %s * %s;\n", alpha0.get(), Q2.get(i, j));
					} else {
						gen.print(";\n");
					}
				}
			}
		}
		// Theta
		gen.print("\n");
		gen.print("F.%s = -%s * %s + %s * (", Theta0.get(), beta1.get(dim), Theta0.get(), alpha0.get());
		for (int k = 0; k < NDIM; k++) {
			gen.print("%s * %s%s", gamma2inv.get(dim, k), V1.get(k), plus[k]);
		}
		gen.print(");\n");
		// Z
		gen.print("\n");
		for (int k = 0; k < NDIM; k++) {
			if (k == dim) {
				gen.print("F.%s = -%s * %s - %s * (%s - %s + %s);\n", Z1.get(k), beta1.get(dim), Z1.get(k), alpha0.get(), K2U.get(dim, k), K0.get(),
						Theta0.get());
			} else {
				gen.print("F.%s = -%s * %s - %s * %s;\n", Z1.get(k), beta1.get(dim), Z1.get(k), alpha0.get(), K2U.get(dim, k));
			}
		}
		gen.print("\n");
		gen.print("break;");
		gen.dedent();
		gen.print("\n}\n");
	}
	gen.print("\ndefault:");
	gen.indent();
	gen.print("\nthrow std::invalid_argument( \"Index for Spacetime flux must be 0, 1, or 2.\" );");
	gen.dedent();
	gen.print("\n\n}\n");
	gen.print("\nreturn F;");
	gen.dedent();
	gen.print("\n}\n\n");

	gen.print("Spacetime source() const {");
	gen.indent();
	gen.print("\nSpacetime S;\n");

	gen.print("\n");
	compInvGamma();

	gen.print("\n");
	gen.print("T const  %s = ", B0.get());
	for (int k = 0; k < NDIM; k++) {
		gen.print("%s", B2.get(k, k));
		if (k + 1 < NDIM) {
			gen.print(" + ");
		}
	}
	gen.print(";\n");

	gen.print("\n");
	gen.print("T const  %s = ", K0.get());
	for (int n = 0; n < NDIM; n++) {
		gen.print("%s * %s + ", gamma2inv.get(n, n), K2.get(n, n));
	}
	gen.print("T(2) * (");
	for (int n = 0; n < NDIM; n++) {
		for (int k = n + 1; k < NDIM; k++) {
			gen.print("%s * %s", gamma2inv.get(n, k), K2.get(k, n));
			if (k + n < 2 * NDIM - 3) {
				gen.print(" + ");
			}
		}
	}
	gen.print(");\n");

	gen.print("\n");
	for (int n = 0; n < NDIM; n++) {
		gen.print("T const  %s = ", D1.get(n));
		for (int k = 0; k < NDIM; k++) {
			gen.print("%s * %s + ", D3.get(n, k, k), gamma2inv.get(k, k));
		}
		gen.print("T(2) * (");
		for (int k = 0; k < NDIM; k++) {
			for (int m = k + 1; m < NDIM; m++) {
				gen.print("%s * %s", D3.get(n, k, m), gamma2inv.get(m, k));
				if (k + m < 2 * NDIM - 3) {
					gen.print(" + ");
				} else {
					gen.print(")");
				}
			}
		}
		gen.print(";\n");
	}

	gen.print("\n");
	for (int n = 0; n < NDIM; n++) {
		gen.print("T const  %s = ", E1.get(n));
		for (int k = 0; k < NDIM; k++) {
			for (int m = 0; m < NDIM; m++) {
				gen.print("%s * %s", D3.get(k, n, m), gamma2inv.get(m, k));
				if (k + m < 2 * NDIM - 2) {
					gen.print(" + ");
				}
			}
		}
		gen.print(";\n");
	}

	gen.print("\n");
	for (int n = 0; n < NDIM; n++) {
		gen.print("T const  %s = ", D1U.get(n));
		for (int k = 0; k < NDIM; k++) {
			gen.print("%s * %s%s", D1.get(k), gamma2inv.get(k, n), plus[k]);
		}
		gen.print(";\n");
	}

	gen.print("\n");
	for (int n = 0; n < NDIM; n++) {
		gen.print("T const  %s = ", E1U.get(n));
		for (int k = 0; k < NDIM; k++) {
			gen.print("%s * %s%s", E1.get(k), gamma2inv.get(k, n), plus[k]);
		}
		gen.print(";\n");
	}

	gen.print("\n");
	for (int n = 0; n < NDIM; n++) {
		gen.print("T const  %s = ", Z1U.get(n));
		for (int k = 0; k < NDIM; k++) {
			gen.print("%s * %s%s", Z1.get(k), gamma2inv.get(k, n), plus[k]);
		}
		gen.print(";\n");
	}

	gen.print("\n");
	for (int i = 0; i < NDIM; i++) {
		for (int j = 0; j < NDIM; j++) {
			gen.print("T const  %s = ", K2U.get(i, j));
			for (int m = 0; m < NDIM; m++) {
				gen.print("%s * %s", gamma2inv.get(i, m), K2.get(m, j));
				if (m < NDIM - 1) {
					gen.print(" + ");
				}
			}
			gen.print(";\n");
		}
	}

	gen.print("\n");
	for (int k = 0; k < NDIM; k++) {
		for (int i = 0; i < NDIM; i++) {
			for (int j = i; j < NDIM; j++) {
				gen.print("T const  %s = ", D3UU.get(k, i, j));
				for (int m = 0; m < NDIM; m++) {
					gen.print("%s * (", gamma2inv.get(i, m));
					for (int r = 0; r < NDIM; r++) {
						gen.print("%s * %s%s", D3.get(k, m, r), gamma2inv.get(r, j), plus[r]);
					}
					gen.print(")%s", plus[m]);
				}
				gen.print(";\n");
			}
		}
	}

	gen.print("\n");
	for (int k = 0; k < NDIM; k++) {
		for (int ij = 0; ij < 6; ij++) {
			int i = I[ij];
			int j = J[ij];
			gen.print("T const %c%s = ", (i == k || j == k) ? '&' : ' ', Gamma3l.get(k, i, j));
			if (i == j && i == k) {
				gen.print("%s", D3.get(i, j, k));
			} else if (i == j) {
				gen.print("T(2) * %s", D3.get(i, i, k));
				gen.print(" - ");
				gen.print("%s", D3.get(k, i, j));
			} else if (i == k || j == k) {
				gen.print("%s", D3.get(j, i, k));
			} else {
				gen.print("%s", D3.get(i, j, k));
				gen.print(" + ");
				gen.print("%s", D3.get(j, i, k));
				gen.print(" - ");
				gen.print("%s", D3.get(k, i, j));
			}
			gen.print(";\n");
		}
	}
	gen.print("\n");
	for (int k = 0; k < NDIM; k++) {
		for (int ij = 0; ij < 6; ij++) {
			int i = I[ij];
			int j = J[ij];
			gen.print("T const  %s = T(0.5) * (", Gamma3U.get(k, i, j));
			for (int m = 0; m < NDIM; m++) {
				gen.print("%s * %s", gamma2inv.get(k, m), Gamma3l.get(m, i, j));
				if (m + 1 < NDIM) {
					gen.print(" + ");
				}
			}
			gen.print(");\n");
		}
	}

	gen.print("\n");
	for (int k = 0; k < NDIM; k++) {
		gen.print("T const  %s = %s - %s - %s;\n", V1.get(k), D1.get(k), E1.get(k), Z1.get(k));
	}

	gen.print("\n");
	gen.print("T const  %s = T(1) / %s;\n", alpha0inv.get(), alpha0.get());

	gen.print("\n");
	gen.print("T const  %s = %s - T(2) * %s;\n", Q0.get(), K0.get(), Theta0.get());

	gen.print("\n");
	for (int k = 0; k < NDIM; k++) {
		gen.print("T const  %s = %s * (%s - %s + T(2) * %s);\n", Q1.get(k), alpha0.get(), A1.get(k), D1.get(k), V1.get(k));
	}

	gen.print("\n");
	for (int i = 0; i < NDIM; i++) {
		gen.print("T const  %s = %s - %s * %s;\n", Q2.get(i, i), K2.get(i, i), alpha0inv.get(), B2.get(i, i));
	}
	for (int i = 0; i < NDIM; i++) {
		for (int j = i + 1; j < NDIM; j++) {
			std::string tmp = B2.get(i, j);
			gen.print("T const  %s = %s - T(0.5) * %s * (%s + %s);\n", Q2.get(i, j), K2.get(i, j), alpha0inv.get(), tmp.c_str(), B2.get(j, i));
		}
	}

	//alpha
	gen.print("\n");
	gen.print("S.%s = %s * (", alpha0.get(), alpha0.get());
	for (int k = 0; k < NDIM; k++) {
		gen.print("%s * %s", beta1.get(k), A1.get(k));
		if (k < NDIM - 1) {
			gen.print(" + ");
		}
	}
	gen.print(" - %s * %s);\n", alpha0.get(), Q0.get());

	// beta
	gen.print("\n");
	for (int i = 0; i < NDIM; i++) {
		gen.print("S.%s = ", beta1.get(i));
		for (int k = 0; k < NDIM; k++) {
			gen.print("%s * %s", beta1.get(k), B2.get(k, i));
			if (k < NDIM - 1) {
				gen.print(" + ");
			}
		}
		gen.print(" - %s * (", alpha0.get());
		for (int k = 0; k < NDIM; k++) {
			gen.print("%s * %s", gamma2inv.get(i, k), Q1.get(k));
			if (k < NDIM - 1) {
				gen.print(" + ");
			}
		}
		gen.print(");\n");
	}

	// gamma
	gen.print("\n");
	for (int i = 0; i < NDIM; i++) {
		gen.print("S.%s = T(2) * (", gamma2.get(i, i));
		for (int m = 0; m < NDIM; m++) {
			gen.print("%s * %s%s", beta1.get(m), D3.get(m, i, i), plus[m]);
		}
		gen.print(" - %s * %s + %s);\n", alpha0.get(), K2.get(i, i), B2.get(i, i));
	}
	for (int i = 0; i < NDIM; i++) {
		for (int j = i + 1; j < NDIM; j++) {
			gen.print("S.%s = T(2) * (", gamma2.get(i, j));
			for (int m = 0; m < NDIM; m++) {
				gen.print("%s * %s%s", beta1.get(m), D3.get(m, i, j), plus[m]);
			}
			std::string tmp = B2.get(i, j);
			gen.print(" - %s * %s) + %s + %s;\n", alpha0.get(), K2.get(i, j), tmp.c_str(), B2.get(j, i));
		}
	}

// K
	gen.print("\n");
	for (int ij = 0; ij < 6; ij++) {
		int i = I[ij];
		int j = J[ij];
		gen.print("S.%s = ", K2.get(i, j));
		gen.print("-%s * %s + ", B0.get(), K2.get(i, j));
		for (int k = 0; k < NDIM; k++) {
			gen.print("%s * %s + %s * %s", K2.get(i, k), B2.get(j, k), K2.get(j, k), B2.get(i, k));
			gen.print(" + ");
		}
		gen.print(" %s * (", alpha0.get());
		gen.indent();
		gen.print("\n");
		if (i == j) {
			gen.print("(%s - T(2) * %s) * %s", A1.get(i), Z1.get(i), D1.get(j));
		} else {
			gen.print("T(0.5) * ((%s - T(2) * %s) * %s + (%s - T(2) * %s) * %s)", A1.get(i), Z1.get(i), D1.get(j), A1.get(j), Z1.get(j), D1.get(i));
		}
		gen.print(" + ");
		for (int k = 0; k < NDIM; k++) {
			gen.print("(%s - T(2) * %s) * %s", D1.get(k), Z1.get(k), Gamma3U.get(k, i, j));
			if (k + 1 < NDIM) {
				gen.print(" + ");
			}
		}
		gen.print(" +\n");
		gen.print("-(");
		for (int m = 0; m < NDIM; m++) {
			for (int k = 0; k < NDIM; k++) {
				std::string arg1 = Gamma3U.get(k, m, j);
				gen.print("%s * %s", arg1.c_str(), Gamma3U.get(m, k, i));
				if (k + m + 2 < 2 * NDIM) {
					gen.print(" + ");
				} else {
					gen.print(")");
				}
			}
		}
		gen.print(" +\n");
		gen.print("-T(2) * (");
		for (int m = 0; m < NDIM; m++) {
			for (int k = 0; k < NDIM; k++) {
				gen.print("%s * %s", K2U.get(k, i), K2.get(k, j));
				if (k + m + 2 < 2 * NDIM) {
					gen.print(" + ");
				} else {
					gen.print(") +");
				}
			}
		}
		gen.print("\n%s * %s", Q0.get(), K2.get(i, j));
		gen.print(");");
		gen.dedent();
		gen.print("\n");
	}

	// A
	gen.print("\n");
	for (int k = 0; k < NDIM; k++) {
		gen.print("S.%s = ", A1.get(k));
		for (int l = 0; l < NDIM; l++) {
			gen.print("%s * %s%s", B2.get(k, l), A1.get(l), plus[l]);
		}
		gen.print(" - %s * %s;\n", B0.get(), A1.get(k));
	}

	// B
	gen.print("\n");
	for (int k = 0; k < NDIM; k++) {
		for (int i = 0; i < NDIM; i++) {
			gen.print("S.%s = ", B2.get(k, i));
			for (int l = 0; l < NDIM; l++) {
				std::string tmp = B2.get(k, l);
				gen.print("%s * %s%s", tmp.c_str(), B2.get(l, i), plus[l]);
			}
			gen.print(" - %s * %s;\n", B0.get(), B2.get(k, i));
		}
	}

	// D
	gen.print("\n");
	for (int k = 0; k < NDIM; k++) {
		for (int i = 0; i < NDIM; i++) {
			for (int j = i; j < NDIM; j++) {
				gen.print("S.%s = ", D3.get(k, i, j));
				for (int l = 0; l < NDIM; l++) {
					std::string tmp = B2.get(k, l);
					gen.print("%s * %s%s", tmp.c_str(), D3.get(l, i, j), plus[l]);
				}
				gen.print(" - %s * %s;\n", B0.get(), D3.get(k, i, j));
			}
		}
	}

	// Theta
	gen.print("\n");
	gen.print("S.%s = -%s * %s + T(0.5) * %s * (", Theta0.get(), B0.get(), Theta0.get(), alpha0.get());
	gen.print("T(2) * (");
	for (int k = 0; k < NDIM; k++) {
		gen.print("%s * (%s - %s - T(2) * %s)%s", A1.get(k), D1U.get(k), E1U.get(k), Z1U.get(k), plus[k]);
	}
	gen.print(") + ");
	gen.indent();
	for (int k = 0; k < NDIM; k++) {
		gen.print("\n");
		for (int r = 0; r < NDIM; r++) {
			for (int s = 0; s < NDIM; s++) {
				gen.print("%s * %s", D3UU.get(k, r, s), Gamma3U.get(k, r, s));
				gen.print("%s", plus[s]);
			}
			gen.print("%s", plus[r]);
		}
		gen.print("%s", plus[k]);
	}
	gen.print(" - (\n");
	for (int k = 0; k < NDIM; k++) {
		for (int r = 0; r < NDIM; r++) {
			gen.print("%s * %s + ", K2U.get(k, r), K2U.get(r, k));
		}
	}
	gen.print("\n");
	for (int k = 0; k < NDIM; k++) {
		gen.print("%s * (%s - T(2) * %s) + ", D1U.get(k), D1.get(k), Z1.get(k));
	}
	gen.print("%s * %s)", K0.get(), Q0.get());
	gen.dedent();
	gen.print("\n");
	gen.print(");\n");

	// Z
	gen.print("\n");
	for (int i = 0; i < NDIM; i++) {
		gen.print("S.%s = -%s * %s + ", Z1.get(i), B0.get(), Z1.get(i));
		for (int k = 0; k < NDIM; k++) {
			gen.print("%s * %s + ", B2.get(k, i), Z1.get(k));
		}
		gen.print("%s * (%s * %s + ", alpha0.get(), A1.get(i), Q0.get());
		for (int k = 0; k < NDIM; k++) {
			gen.print("(%s - %s - T(2) * %s) * %s%s", D1.get(k), A1.get(k), Z1.get(k), K2.get(k, i), plus[k]);
		}
		gen.print(" - ");
		gen.indent();
		gen.print("\n(");
		for (int k = 0; k < NDIM; k++) {
			for (int r = 0; r < NDIM; r++) {
				gen.print("%s * %s%s", K2.get(k, r), Gamma3U.get(r, k, i), plus[r]);
			}
			gen.print("%s", plus[k]);
		}
		gen.print(")");
		gen.dedent();
		gen.print("\n);\n");
	}

	gen.print("\n");
	gen.print("return S;");
	gen.dedent();
	gen.print("\n}");

	gen.print("\n");

	gen.dedent();
	gen.print("\nprivate:");
	gen.print("\n");
	gen.indent();
	gen.print("\n");
	gen.print("std::array<T, %i> U;\n", NFIELDS);
	gen.print("\n");
	gen.print("T& %s;\n", alpha0.get());
	for (int i = 0; i < NDIM; i++) {
		gen.print("T %s;\n", beta1.get(i));
	}
	gen.print("\n");
	for (int i = 0; i < NDIM; i++) {
		for (int j = i; j < NDIM; j++) {
			gen.print("T& %s;\n", gamma2.get(i, j));
		}
	}
	gen.print("\n");
	for (int i = 0; i < NDIM; i++) {
		for (int j = i; j < NDIM; j++) {
			gen.print("T& %s;\n", K2.get(i, j));
		}
	}
	gen.print("\n");
	for (int i = 0; i < NDIM; i++) {
		gen.print("T& %s;\n", A1.get(i));
	}
	gen.print("\n");
	for (int i = 0; i < NDIM; i++) {
		for (int j = 0; j < NDIM; j++) {
			gen.print("T& %s;\n", B2.get(i, j));
		}
	}
	gen.print("\n");
	for (int i = 0; i < NDIM; i++) {
		for (int j = 0; j < NDIM; j++) {
			for (int k = j; k < NDIM; k++) {
				gen.print("T& %s;\n", D3.get(i, j, k));
			}
		}
	}
	gen.print("\n");
	gen.print("T& %s;\n", Theta0.get());
	gen.print("\n");
	for (int i = 0; i < NDIM; i++) {
		gen.print("T& %s;\n", Z1.get(i));
	}
	gen.dedent();
	gen.print("\n};\n");
	std::cout << gen.get();
	FILE *fp = fopen("Z4.hpp", "wt");
	fprintf(fp, gen.get().c_str());
	fclose(fp);
	return 0;
}
