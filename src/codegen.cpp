#include "Indent.hpp"
#include "Util.hpp"

#include <algorithm>
#include <unordered_map>
#include <vector>
#include <iostream>

using Real = long double;

struct Constants {
	Constants() {
		setType("double");
	}
	std::string operator()(Real c) {
		std::string format;
		format = "static " + type_ + " const %s(%." + std::to_string(std::numeric_limits<double>::max_digits10) + "e);\n";
		auto it = names_.find(c);
		char *ptr;
		if (it == names_.end()) {
			std::string name = "C" + std::to_string(index_);
			names_.insert(std::make_pair(c, name));
			asprintf(&ptr, format.c_str(), name.c_str(), (double) c);
			code_ += ptr;
			free(ptr);
			it = names_.find(c);
			index_++;
		}
		return it->second;
	}
	std::string getCode() const {
		return code_;
	}
	void reset() {
		*this = Constants { };
	}
	void setType(std::string type) {
		type_ = type;
	}
private:
	std::unordered_map<double, std::string> names_;
	std::string code_;
	std::string type_;
	int index_ = 0;
};

static Constants getConstant;
using Matrix = std::vector<std::vector<Real>>;

Matrix createMatrix(int N, int M = -1) {
	if (M < 0) {
		M = N;
	}
	return std::vector<std::vector<Real>>(N, std::vector<Real>(M, 0.0));
}

Matrix identityMatrix(int N) {
	auto I = createMatrix(N);
	for (int n = 0; n < N; n++) {
		I[n][n] = 1.0;
	}
	return I;
}

Matrix matrixTranspose(Matrix const &A) {
	int rows = A.size(), cols = A[0].size();
	Matrix AT(cols, std::vector<Real>(rows));
	for (int i = 0; i < rows; ++i) {
		for (int j = 0; j < cols; ++j) {
			AT[j][i] = A[i][j];
		}
	}
	return AT;
}

Matrix matrixMultiply(Matrix const &A, Matrix const &B) {
	int rows = A.size(), cols = B[0].size(), inner = A[0].size();
	Matrix C(rows, std::vector<Real>(cols, 0.0));
	for (int i = 0; i < rows; ++i) {
		for (int k = 0; k < inner; ++k) {
			for (int j = 0; j < cols; ++j) {
				C[i][j] += A[i][k] * B[k][j];
			}
		}
	}
	return C;
}

Matrix matrixInverse(Matrix const &A) {
	int const N = A.size();
	for (auto const &row : A) {
		if (int(row.size()) != N) {
			throw std::runtime_error("Matrix is not square");
		}
	}
	Matrix I = identityMatrix(N);
	Matrix B = A;
	for (int i = 0; i < N; ++i) {
		//	std::cout << matrixToString(A) << "\n";
		int pivot = i;
		if (B[i][i] < std::sqrt(std::numeric_limits<double>::epsilon())) {
			int maxB = B[i][i];
			for (int row = i + 1; row < N; ++row) {
				if (std::abs(B[row][i]) > std::abs(maxB)) {
					maxB = B[row][i];
					pivot = row;
				}
			}
		}
		if (std::abs(B[pivot][i]) < std::sqrt(std::numeric_limits<Real>::epsilon())) {
			throw std::runtime_error("Matrix is singular or nearly singular");
		}
		if (pivot != i) {
			std::swap(B[i], B[pivot]);
			std::swap(I[i], I[pivot]);
		}
		Real diag = B[i][i];
		for (int j = 0; j < N; ++j) {
			B[i][j] /= diag;
			I[i][j] /= diag;
		}
		for (int k = 0; k < N; ++k) {
			if (k == i) {
				continue;
			}
			Real factor = B[k][i];
			for (int j = 0; j < N; ++j) {
				B[k][j] -= factor * B[i][j];
				I[k][j] -= factor * I[i][j];
			}
		}
	}
	return I;
}

Matrix leftInverse(Matrix const &R) {
	// (Rt * R)^-1 * Rt   R
	Matrix Rt = matrixTranspose(R);
	Matrix RtR = matrixMultiply(Rt, R);
	Matrix RtR_inv = matrixInverse(RtR);
	return matrixMultiply(RtR_inv, Rt);
}

Matrix rightInverse(Matrix const &L) {
	// L   Lt (L * Lt)^-1
	Matrix Lt = matrixTranspose(L);
	Matrix LLt = matrixMultiply(L, Lt);
	Matrix LLt_inv = matrixInverse(LLt);
	auto result = matrixMultiply(Lt, LLt_inv);
	return result;
}

Matrix kroneckerProduct(Matrix const &A, Matrix const &B) {
	int const NR = A.size();
	int const NC = A[0].size();
	int const MR = B.size();
	int const MC = B[0].size();
	int const LR = NR * MR;
	int const LC = NC * MC;
	auto C = createMatrix(LR, LC);
	for (int nr = 0; nr < NR; nr++) {
		for (int nc = 0; nc < NC; nc++) {
			for (int mr = 0; mr < MR; mr++) {
				for (int mc = 0; mc < MC; mc++) {
					C[nr * MR + mr][nc * MC + mc] = A[nr][nc] * B[mr][mc];
				}
			}
		}
	}
	return C;
}

std::string matrixToString(Matrix const &A) {
	char *ptr;
	std::string str;
	for (int n = 0; n < (int) A.size(); n++) {
		for (int m = 0; m < (int) A[n].size(); m++) {
			asprintf(&ptr, "%.0f ", (double) A[n][m]);
			str += ptr;
			free(ptr);
		}
		str += "\n";
	}
	str += "\n";
	return str;
}

enum class Quadrature : int {
	gaussLegendre, gaussLobatto
};

enum class TransformDirection : int {
	forward, backward
};

struct QuadraturePoint {
	Real position;
	Real weight;
};

Real legendreP(int n, Real x, int m = 0) {
	if (std::abs(x) != Real(1)) {
		return Real(nonepow(m)) * std::sqrt(Real(1) / ipow(Real(1) - x * x, m)) * std::assoc_legendre(n, m, x);
	} else {
		Real value = Real(1) / Real(2);
		for (int m = 1; m < n; m++) {
			value *= Real(n + m) / Real(2 * m * (n - m));
		}
		return Real((x > Real(0)) ? nonepow(n + m) : +1) * value;
	}
}

std::vector<QuadraturePoint> gaussQuadrature(int nodeCount, Quadrature quadrature) {
	constexpr Real pi = std::numbers::pi_v<Real>;
	std::vector<QuadraturePoint> results;
	switch (quadrature) {
	case Quadrature::gaussLegendre: {
		for (int pointIndex = 0; pointIndex < nodeCount; pointIndex++) {
			QuadraturePoint point;
			Real theta, rootEquationDerivative, rootEquation;
			Real newTheta = pi * (Real(1) - Real(2 * pointIndex + 1) / Real(2) / Real(nodeCount));
			do {
				theta = newTheta;
				point.position = std::cos(theta);
				rootEquation = legendreP(nodeCount, point.position);
				rootEquationDerivative = legendreP(nodeCount, point.position, 1) * std::sin(theta);
				newTheta = theta - rootEquation / rootEquationDerivative;
			} while (double(newTheta) != std::nextafter(double(theta), double(newTheta)));
			point.weight = Real(2) / (sqr(legendreP(nodeCount, point.position, 1)) * (Real(1) - sqr(point.position)));
			results.push_back(point);
		}
		break;
	}
	case Quadrature::gaussLobatto: {
		Real const baseWeight = Real(2) / (nodeCount * (nodeCount - 1));
		QuadraturePoint edgePoint = { -Real(1), baseWeight };
		results.push_back(edgePoint);
		for (int pointIndex = 1; pointIndex < nodeCount - 1; pointIndex++) {
			QuadraturePoint point;
			Real theta, rootEquationDerivative, rootEquation;
			Real newTheta = pi * (Real(1) - Real(2 * pointIndex + 1) / Real(2) / Real(nodeCount));
			do {
				theta = newTheta;
				point.position = std::cos(theta);
				rootEquation = legendreP(nodeCount - 1, point.position, 1);
				rootEquationDerivative = -legendreP(nodeCount - 1, point.position, 2) * std::sin(theta);
				newTheta = theta + rootEquation / rootEquationDerivative;
			} while (double(newTheta) != std::nextafter(double(theta), double(newTheta)));
			point.weight = baseWeight / sqr(legendreP(nodeCount - 1, point.position));
			results.push_back(point);
		}
		edgePoint.position = -edgePoint.position;
		results.push_back(edgePoint);
		break;
	}
	}
	return results;
}

//Matrix legendreTransform1d(int basisOrder, int quadratureLength, int surfaceSign) {
//	if (surfaceSign == 0) {
//		auto T = createMatrix(basisOrder, quadratureLength);
//		for (int n = 0; n < basisOrder; n++) {
//			for (int m = 0; m < quadratureLength; m++) {
//				auto qpt = gaussLegendreQuadrature(m, quadratureLength);
//				T[n][m] = qpt.weight * (n + 0.5) * std::legendre(n, qpt.position);
//			}
//		}
//		return T;
//	} else {
//		auto T = createMatrix(basisOrder, 1);
//		for (int n = 0; n < basisOrder; n++) {
//			T[n][0] = Real(surfaceSign) * std::pow(surfaceSign, n);
//			printf("%c1 ", T[n][0] > 0 ? '+' : '-');
//		}
//		printf("\n");
//		return T;
//	}
//}
//
//Matrix legendreTransformDerivative1d(int basisOrder, int quadratureLength) {
//	auto T = createMatrix(basisOrder, quadratureLength);
//	for (int n = 0; n < basisOrder; n++) {
//		for (int m = 0; m < quadratureLength; m++) {
//			auto qpt = gaussLegendreQuadrature(m, quadratureLength);
//			T[n][m] = qpt.weight * (n + 0.5) * dlegendre_dx(n, qpt.position);
//		}
//	}
//	return T;
//}
//
//
//constexpr long long factorial(int n) {
//	if (n) {
//		return (long long) (n) * factorial(n - 1);
//	} else {
//		return (long long) 1;
//	}
//}
//
//constexpr long long binco(int n, int k) {
//	return factorial(n) / (factorial(n - k) * factorial(k));
//}
//
//std::vector<int> computeVectorSizes(int basisOrder, int quadratureLength, int dim, int D) {
//	std::vector<int> sizes;
//	int size = 1;
//	int bDim = 0;
//	for (int d = 0; d < D; d++) {
//		if (d <= dim) {
//			size *= quadratureLength;
//			sizes.push_back(size);
//		} else {
//			size = size * (basisOrder + bDim) / (bDim + 1);
//			sizes.push_back(size);
//			bDim++;
//		}
//	}
//	return sizes;
//}
//// P = A^-1 x I x I x A
//Matrix legendreTransform(int basisOrder, int quadratureLength, int dim, int D, int surfaceDim, bool derivative) {
//	auto const vectorSizes = computeVectorSizes(basisOrder, quadratureLength, dim, D);
//	auto A = identityMatrix(1);
//	printf("%i\n", dim);
//	for (int d = D - 1; d >= 0; d--) {
//		Matrix B;
//		int const surfaceSign = int(d + 1 == std::abs(surfaceDim)) * (surfaceDim > 0 ? +1 : -1);
//		printf("   %i %i\n", d, surfaceSign);
//		if (d == dim) {
//			if (derivative) {
//				B = legendreTransformDerivative1d(basisOrder, quadratureLength);
//			} else if (surfaceSign == 0) {
//				B = legendreTransform1d(basisOrder, quadratureLength, 0);
//			} else {
//				B = legendreTransform1d(basisOrder, 1, surfaceSign);
//			}
//		} else if (d < dim) {
//			B = identityMatrix(quadratureLength);
//		} else {
//			B = identityMatrix(basisOrder);
//		}
//		A = kroneckerProduct(B, A);
//	}
//	std::vector<int> indices(D);
//	std::vector<bool> columnZeros(A[0].size());
//	std::vector<bool> rowZeros(A.size());
//	std::fill(columnZeros.begin(), columnZeros.end(), false);
//	std::fill(rowZeros.begin(), rowZeros.end(), false);
//	for (int n = 0; n < (int) A.size(); n++) {
//		int i = n;
//		for (int d = D - 1; d >= 0; d--) {
//			int stride = d >= dim ? basisOrder : quadratureLength;
//			indices[d] = i % stride;
//			i /= stride;
//		}
//		int sum = 0;
//		for (int d = 0; d < D; d++) {
//			if (d >= dim) {
//				sum += indices[d];
//			}
//		}
//		if (sum >= basisOrder) {
//			rowZeros[n] = true;
//		}
//	}
//	for (int m = 0; m < (int) A[0].size(); m++) {
//		int i = m;
//		for (int d = D - 1; d >= 0; d--) {
//			int stride = d > dim ? basisOrder : quadratureLength;
//			indices[d] = i % stride;
//			i /= stride;
//		}
//		int sum = 0;
//		for (int d = 0; d < D; d++) {
//			if (d > dim) {
//				sum += indices[d];
//			}
//		}
//		if (sum >= basisOrder) {
//			columnZeros[m] = true;
//		}
//	}
//	auto A0 = std::move(A);
//	A.resize(0);
//	auto it = rowZeros.begin();
//	for (auto &&a : A0) {
//		if (*it == false) {
//			A.push_back(std::move(a));
//		}
//		it++;
//	}
//	for (int n = 0; n < (int) A.size(); n++) {
//		auto A0 = std::move(A[n]);
//		A[n].resize(0);
//		auto it = columnZeros.begin();
//		for (auto &&a : A0) {
//			if (*it == false) {
//				A[n].push_back(std::move(a));
//			}
//			it++;
//		}
//	}
//	return A;
//}
//
//constexpr size_t getMaxSimdLaneWidth() {
//	return std::experimental::simd_abi::max_fixed_size<double> / sizeof(double);
//}
//
//
//using Matrix = std::vector<std::vector<Real>>;
//
//
//std::string matrixVectorProduct(std::vector<std::string> const &v, Matrix const &A, std::vector<std::string> const &x, Constants &getConstant) {
//	std::string code;
//	int const M = x.size();
//	int const N = A.size();
//	assert((int ) v.size() == N);
//	assert((int ) A[0].size() == M);
//	for (int n = 0; n < N; n++) {
//		if (v[n] == "") {
//			continue;
//		}
//		bool first = true;
//		char *ptr;
//		for (int m = 0; m < M; m++) {
//			if ((std::abs(A[n][m]) < std::numeric_limits<double>::epsilon()) || (x[m] == "")) {
//				continue;
//			}
//			if (double(A[n][m]) == double(1)) {
//				if (first) {
//					asprintf(&ptr, "%s = %s;\n", v[n].c_str(), x[m].c_str());
//				} else {
//					asprintf(&ptr, "%s += %s;\n", v[n].c_str(), x[m].c_str());
//				}
//			} else {
//				if (first) {
//					asprintf(&ptr, "%s = %s * %s;\n", v[n].c_str(), getConstant(A[n][m]).c_str(), x[m].c_str());
//				} else {
//					asprintf(&ptr, "%s = fma(%s, %s, %s);\n", v[n].c_str(), getConstant(A[n][m]).c_str(), x[m].c_str(), v[n].c_str());
//				}
//			}
//			code += indent;
//			code += ptr;
//			free(ptr);
//			first = false;
//		}
//		if (first) {
//			asprintf(&ptr, "%s = 0.0;\n", v[n].c_str());
//			code += std::string(indent) + ptr;
//			free(ptr);
//		}
//	}
//	return code;
//}
//
//Matrix leftInverse(Matrix const &R) {
//	// (Rt * R)^-1 * Rt   R
//	Matrix Rt = transpose(R);
//	Matrix RtR = multiply(Rt, R);
//	Matrix RtR_inv = inverse(RtR);
//	return multiply(RtR_inv, Rt);
//}
//
//Matrix rightInverse(Matrix const &L) {
//	// L   Lt (L * Lt)^-1
//	Matrix Lt = transpose(L);
//	Matrix LLt = multiply(L, Lt);
//	Matrix LLt_inv = inverse(LLt);
//	auto result = multiply(Lt, LLt_inv);
//	return result;
//}
//
//std::string generateLegendreTransformCode(std::string const &type, std::string header, int basisOrder, int quadratureLength, int D, int surfaceDim,
//		int derivativeDim = -1) {
//	std::vector<std::string> input;
//	std::vector<std::string> output;
//	std::string code;
//	Constants constants(type);
//	std::string arrays;
//	for (int dim = D - 1; dim >= 0; dim--) {
//		auto A = legendreTransform(basisOrder, quadratureLength, dim, D, surfaceDim, derivativeDim == dim);
//		input.resize(A[0].size());
//		output.resize(A.size());
//		arrays += std::string(indent) + "std::array<" + type + ", " + std::to_string(A.size()) + "> " + std::string(1, 'x' + dim) + ";\n";
//		for (int i = 0; i < (int) output.size(); i++) {
//			output[i] = std::string(1, 'x' + dim) + "[" + std::to_string(i) + "]";
//		}
//		for (int i = 0; i < (int) input.size(); i++) {
//			input[i] = std::string(1, 'x' + dim - 1) + "[" + std::to_string(i) + "]";
//		}
//		code += matrixVectorProduct(output, A, input, constants);
//	}
//	code += std::string(indent) + "return " + std::string(1, 'x') + ";\n";
//	return header + constants.getCode() + arrays + code;
//}
//
//std::string generateInverseLegendreTransformCode(std::string const &type, std::string header, int basisOrder, int quadratureLength, int D, int surfaceSign) {
//	std::vector<std::string> input;
//	std::vector<std::string> output;
//	std::string code;
//	Constants constants(type);
//	std::string arrays;
//	for (int dim = 0; dim < D; dim++) {
//		auto A = rightInverse(legendreTransform(basisOrder, quadratureLength, dim, D, 0, false));
//		input.resize(A[0].size());
//		output.resize(A.size());
//		arrays += std::string(indent) + "std::array<" + type + ", " + std::to_string(A.size()) + "> " + std::string(1, 'x' + dim) + ";\n";
//		for (int i = 0; i < (int) output.size(); i++) {
//			output[i] = std::string(1, 'x' + dim) + "[" + std::to_string(i) + "]";
//		}
//		for (int i = 0; i < (int) input.size(); i++) {
//			input[i] = std::string(1, 'x' + dim - 1) + "[" + std::to_string(i) + "]";
//		}
//		code += matrixVectorProduct(output, A, input, constants);
//	}
//	code += std::string(indent) + "return " + std::string(1, 'x' + D - 1) + ";\n";
//	return header + constants.getCode() + arrays + code;
//}
//
//std::string generatePair(std::string const &name, std::string const &type, int basisOrder, int quadratureLength, int D) {
//	int transformSize = binco(basisOrder + D - 1, D);
//	int inverseTransformSize = int(std::pow(quadratureLength, D));
//	std::string code;
//	std::string header;
//	header = "std::array<" + type + ", ";
//	header += std::to_string(transformSize);
//	header += "> " + name + "ForwardTransform(std::array<" + type + ", " + std::to_string(inverseTransformSize) + "> const& w) {\n";
//	indent++;
//	code = generateLegendreTransformCode(type, header, basisOrder, quadratureLength, D, 0);
//	indent--;
//	code += "}\n\n";
//	header = std::string(indent) + "std::array<" + type + ", ";
//	header += std::to_string(inverseTransformSize);
//	header += "> " + name + "BackwardTransform(std::array<" + type + ", " + std::to_string(transformSize) + "> const& w) {\n";
//	indent++;
//	code += generateInverseLegendreTransformCode(type, header, basisOrder, quadratureLength, D, 0);
//	indent--;
//	code += "}\n\n";
//	return code;
//}
//
//std::string generateSurface(std::string const &name, std::string const &type, int basisOrder, int quadratureLength, int D) {
//	int transformSize = binco(basisOrder + D - 1, D);
//	int inverseTransformSize = int(std::pow(quadratureLength, D));
//	std::string code;
//	std::string header;
//	for (int dir = +1; dir >= -1; dir -= 2) {
//		auto const transform = [dir](std::string A, std::string B, int a, int b, int c, int d) {
//			if (dir > 0) {
//				return generateLegendreTransformCode(A, B, a, b, c, d);
//			}
//			return generateInverseLegendreTransformCode(A, B, a, b, c, d);
//		};
//		code += "std::array<" + type + ", ";
//		code += std::to_string(transformSize);
//		code +=
//				"> " + std::string(dir > 0 ? "forward" : "backward") + "Transform(std::array<" + type + ", " + std::to_string(inverseTransformSize) + "> const& w, Face face = Face(-1)) {\n";
//		indent++;
//		code += std::string(indent) + "if(face == -1) {\n";
//		indent++;
//		code += transform(type, header, basisOrder, quadratureLength, D, 0);
//		indent--;
//		for (int dim = 0; dim < D; dim++) {
//			code += std::string(indent) + "} else if(face == " + std::to_string(2 * dim) + ") {\n";
//			indent++;
//			code += transform(type, header, basisOrder, quadratureLength, D, +(dim + 1));
//			indent--;
//			code += std::string(indent) + "} else ";
//			if (dim + 1 == D) {
//				code += "/*";
//			}
//			code += "if(face == " + std::to_string(2 * dim + 1) + ")";
//			if (dim + 1 == D) {
//				code += "*/";
//			}
//			code += +" {\n";
//			indent++;
//			code += transform(type, header, basisOrder, quadratureLength, D, -(dim + 1));
//			indent--;
//			if (dim + 1 == D) {
//				code += std::string(indent) + "}\n";
//			}
//		}
//		indent--;
//		code += std::string(indent) + "}\n\n";
//		std::swap(transformSize, inverseTransformSize);
//	}
//	return code;
//}
//
//std::string generateDerivatives(std::string const &name, std::string const &type, int basisOrder, int quadratureLength, int D) {
//	int transformSize = binco(basisOrder + D - 1, D);
//	int inverseTransformSize = int(std::pow(quadratureLength, D));
//	std::string code;
//	std::string header;
//	header = "std::array<" + type + ", ";
//	header += std::to_string(transformSize);
//	header += "> " + name + "Transform(int dim, std::array<" + type + ", " + std::to_string(inverseTransformSize) + "> const& w) {\n";
//	for (int dim = 0; dim < D; dim++) {
//		indent++;
//		if (dim == 0) {
//			code += std::string(indent);
//		}
//		if (dim + 1 == D) {
//			code += "/*";
//		}
//		code += "if(dim == " + std::to_string(dim) + ")";
//		if (dim + 1 == D) {
//			code += "*/";
//		}
//		code += " {\n";
//		indent++;
//		code += generateLegendreTransformCode(type, "", basisOrder, quadratureLength, D, 0, dim);
//		indent--;
//		code += std::string(indent) + "}";
//		code += (dim + 1 < D) ? " else " : "\n";
//		indent--;
//	}
//	code += "}\n\n";
//	return header + code;
//}
//
//std::vector<std::vector<signType>> makePermutationMatrix(Matrix const &B, int simdWidth) {
//	int const N = B.size();
//	int const M = B[0].size();
//	std::vector<signType> X(M);
//	std::vector<std::vector<signType>> A(N, std::vector<signType>(M));
//	std::vector<std::vector<signType>> P(M, std::vector<signType>(M));
//	std::iota(X.begin(), X.end(), signType(0));
//	for (int n = 0; n < M; n++) {
//		for (int m = 0; m < M; m++) {
//			A[n][m] = signType(std::abs(B[n][m]) >= std::numeric_limits<double>::epsilon());
//		}
//	}
//	return P;
//}

Matrix transformMatrix1D(TransformDirection transformDirection, int modeCount, int nodeCount, Quadrature quadratureType, int derivativeOrder = 0) {
	Matrix transform;
	if (transformDirection == TransformDirection::forward) {
		transform = createMatrix(modeCount, nodeCount);
		auto const quadratureRules = gaussQuadrature(nodeCount, quadratureType);
		for (int modeIndex = 0; modeIndex < modeCount; modeIndex++) {
			for (int nodeIndex = 0; nodeIndex < nodeCount; nodeIndex++) {
				Real const position = quadratureRules[nodeIndex].position;
				Real const weight = quadratureRules[nodeIndex].weight;
				Real const inverseMass = Real(2 * modeIndex + 1) / Real(2);
				transform[modeIndex][nodeIndex] = legendreP(modeIndex, position, derivativeOrder) * inverseMass * weight;
			}
		}
	} else {
		assert(derivativeOrder == 0);
		transform = createMatrix(nodeCount, modeCount);
		auto const quadratureRules = gaussQuadrature(nodeCount, quadratureType);
		for (int modeIndex = 0; modeIndex < modeCount; modeIndex++) {
			for (int nodeIndex = 0; nodeIndex < nodeCount; nodeIndex++) {
				Real const position = quadratureRules[nodeIndex].position;
				transform[nodeIndex][modeIndex] = legendreP(modeIndex, position);
			}
		}
	}
	return transform;
}

std::string matrixVectorProduct(std::vector<std::string> const &v, Matrix const &A, std::vector<std::string> const &x) {
	std::string code;
	int const M = x.size();
	int const N = A.size();
	assert((int ) v.size() == N);
	assert((int ) A[0].size() == M);
	for (int n = 0; n < N; n++) {
		if (v[n] == "") {
			continue;
		}
		bool first = true;
		char *ptr;
		for (int m = 0; m < M; m++) {
			if ((std::abs(A[n][m]) < std::numeric_limits<double>::epsilon()) || (x[m] == "")) {
				continue;
			}
			if (double(A[n][m]) == double(1)) {
				if (first) {
					asprintf(&ptr, "%s = %s;\n", v[n].c_str(), x[m].c_str());
				} else {
					asprintf(&ptr, "%s += %s;\n", v[n].c_str(), x[m].c_str());
				}
			} else {
				if (first) {
					asprintf(&ptr, "%s = %s * %s;\n", v[n].c_str(), getConstant(A[n][m]).c_str(), x[m].c_str());
				} else {
					asprintf(&ptr, "%s = fma(%s, %s, %s);\n", v[n].c_str(), getConstant(A[n][m]).c_str(), x[m].c_str(), v[n].c_str());
				}
			}
			code += ptr;
			free(ptr);
			first = false;
		}
		if (first) {
			asprintf(&ptr, "%s = 0.0;\n", v[n].c_str());
			code += ptr;
			free(ptr);
		}
	}
	return code;
}

std::vector<std::string> stringToLines(std::string text) {
	std::vector<std::string> lines;
	std::string line;
	for (int i = 0; i < (int) text.size(); i++) {
		auto const c = text[i];
		line.push_back(c);
		if (c == '\n') {
			lines.push_back(std::move(line));
		}
	}
	return lines;
}

struct Function {
	Function(std::string name, std::string returnName, std::string body) :
			name_(name), returnName_(returnName), body_(body) {
	}
	std::string getPrototype() const {
		std::string code;
		code += returnType_ + " " + name_ + "(";
		bool firstPass = true;
		for (auto parameter : parameters_) {
			if (!firstPass) {
				code += ", ";
			}
			code += parameter;
			firstPass = false;
		}
		code += ")";
		return code;
	}
	std::string toHpp() const {
		return getPrototype() + ";\n";
	}
	std::string toCpp() const {
		std::string code = getPrototype() + " {\n";
		auto text = std::string(getConstant.getCode()) + body_ + "return " + returnName_ + ";\n";
		auto lines = stringToLines(text);
		for (auto const &line : lines) {
			code += "   " + line;
		}
		code += "}\n\n";
		return code;
	}
	template<typename ...Args>
	void setParameters(Args ...args) {
		int i = 0;
		parameters_.resize(sizeof...(Args));
		((parameters_[i++] = args), ... );
	}
private:
	std::string name_;
	std::string returnName_;
	std::string body_;
	std::string returnType_ = "auto";
	std::vector<std::string> parameters_ = std::vector<std::string>(1, "auto");
};

struct TransformParameters {
	int dimensionCount = DIMENSION_COUNT;
	int normalDimension = 0;
	int modeCount = MODE_COUNT;
	int normalNodes = MODE_COUNT;
	int tangentialNodes = MODE_COUNT;
	Quadrature normalQuadrature = Quadrature::gaussLegendre;
	Quadrature tangentialQuadrature = Quadrature::gaussLegendre;
	TransformDirection transformDirection = TransformDirection::backward;
	int derivativeOrder = 0;
};

int flipDim(int dimensionIndex, int dimensionCount = DIMENSION_COUNT) {
	return dimensionCount - 1 - dimensionIndex;
}

Matrix permutationMatrix(int N, int modeCount, int triIndexCount) {
	std::vector<bool> ones(N, true);
	int M = 0;
	for (int i = 0; i < N; i++) {
		int deg = 0;
		int k = i;
		for (int j = 0; j < triIndexCount; j++) {
			deg += k % modeCount;
			k /= modeCount;
		}
		ones[i] = (deg < modeCount);
		if (ones[i]) {
			M++;
		}
	}
	Matrix P = createMatrix(M, N);
	int m = 0;
	for (int n = 0; n < N; n++) {
		if (ones[n]) {
			P[m++][n] = Real(1.0);
		}
	}
	return P;
}

Matrix transformMatrix(int transformDimension, TransformParameters const &parameters) {
	auto A = identityMatrix(1);
	Matrix B;
	std::vector<int> sizes;
	assert(parameters.transformDirection != TransformDirection::backward);
	for (int thisDimension = parameters.dimensionCount - 1; thisDimension >= 0; thisDimension--) {
		int nodeCount, derivativeOrder;
		Quadrature quadratureType;
		if (thisDimension == parameters.normalDimension) {
			nodeCount = parameters.normalNodes;
			quadratureType = parameters.normalQuadrature;
			derivativeOrder = parameters.derivativeOrder;
		} else {
			nodeCount = parameters.tangentialNodes;
			quadratureType = parameters.tangentialQuadrature;
			derivativeOrder = 0;
		}
		if (thisDimension < transformDimension) {
			A = kroneckerProduct(identityMatrix(nodeCount), A);
		} else if (thisDimension > transformDimension) {
			A = kroneckerProduct(identityMatrix(parameters.modeCount), A);
		} else/*if( thisDimension == transformDimension*/{
			A = kroneckerProduct(transformMatrix1D(TransformDirection::forward, parameters.modeCount, nodeCount, quadratureType, derivativeOrder), A);
		}
	}
	int const columnCount = A[0].size();
	int const rowCount = A.size();
	int const outputTriIndexCount = parameters.dimensionCount - transformDimension;
	int const inputTriIndexCount = outputTriIndexCount - 1;
	if (inputTriIndexCount >= 2) {
		auto permute = permutationMatrix(columnCount, parameters.modeCount, inputTriIndexCount);
		A = matrixMultiply(A, matrixTranspose(permute));
	}
	if (outputTriIndexCount >= 2) {
		auto permute = permutationMatrix(rowCount, parameters.modeCount, outputTriIndexCount);
		A = matrixMultiply(permute, A);
	}
	return A;
}

std::vector<std::string> generateVariableNames(std::string const &name, int count) {
	std::vector<std::string> names(count);
	for (int nameIndex = 0; nameIndex < count; nameIndex++) {
		names[nameIndex] = name + "[" + std::to_string(nameIndex) + "]";
	}
	return names;
}

std::vector<std::string> generateVariableNames(char name, int count) {
	return generateVariableNames(std::string(1, name), count);
}

std::string generateArgumentDeclaration(char name, int count, std::string type = "double") {
	return std::string("std::array<") + type + std::string(", ") + std::to_string(count) + std::string("> const& ") + std::string(1, name);
}

std::string generateVariableDeclaration(char name, int count, std::string type = "double") {
	return std::string("std::array<") + type + std::string(", ") + std::to_string(count) + std::string("> ") + std::string(1, name) + ";\n";
}

int main(int, char*[]) {
	constexpr int dimensionCount = DIMENSION_COUNT;
	std::string cppCode;
	std::string hppCode;
	hppCode += "#pragma once\n";
	hppCode += "\n";
	hppCode += "#include <array>\n";
	hppCode += "\n";
	cppCode += "\n";
	cppCode += "#include \"transforms.hpp\"\n";
	cppCode += "\n";
	cppCode += "#include <cmath>\n";
	cppCode += "\n";
	TransformParameters parameters;
	parameters.dimensionCount = DIMENSION_COUNT;
	parameters.modeCount = MODE_COUNT;
	parameters.normalNodes = MODE_COUNT + 2;
	parameters.tangentialNodes = MODE_COUNT;
	parameters.normalQuadrature = Quadrature::gaussLobatto;
	parameters.tangentialQuadrature = Quadrature::gaussLegendre;
	parameters.derivativeOrder = 0;
	for (parameters.normalDimension = 0; parameters.normalDimension < dimensionCount; parameters.normalDimension++) {
		std::vector<Matrix> factors;
		parameters.transformDirection = TransformDirection::forward;
		for (int transformIndex = 0; transformIndex < dimensionCount; transformIndex++) {
			auto A = transformMatrix(flipDim(transformIndex), parameters);
			factors.push_back(A);
		}
		for (int dir = +1; dir >= -1; dir -= 2) {
			std::string impls, decls, argument;
			for (int dimensionIndex = 0; dimensionIndex < dimensionCount; dimensionIndex++) {
				auto const &A = factors[dimensionIndex];
				int N = A.size();
				int M = A[0].size();
				std::vector<std::string> inputs;
				std::vector<std::string> outputs;
				inputs = generateVariableNames('u' + dimensionIndex, M);
				outputs = generateVariableNames('v' + dimensionIndex, N);
				if (argument == "") {
					argument = generateArgumentDeclaration('u', M);
				}
				impls += matrixVectorProduct(outputs, A, inputs);
				decls += generateVariableDeclaration('v' + dimensionIndex, N);
			}
			impls = decls + impls;
			std::string filename = (parameters.transformDirection == TransformDirection::forward) ? "forward" : "backward";
			filename += "Transform" + std::string(1, 'X' + parameters.normalDimension);
			auto F = Function(filename, std::string(1, 'u' + dimensionCount), impls);
			F.setParameters(argument);
			cppCode += F.toCpp();
			hppCode += F.toHpp();
			getConstant.reset();
			if (dir == +1) {
				for (auto &factor : factors) {
					factor = rightInverse(factor);
				}
				std::reverse(factors.begin(), factors.end());
				parameters.transformDirection = TransformDirection::backward;
			}
		}
	}
	toFile(hppCode, "./generated_source/transforms.hpp");
	toFile(cppCode, "./generated_source/transforms.cpp");
	return 0;
}
