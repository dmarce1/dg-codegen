#include "Definitions.hpp"
#include "Indent.hpp"
#include "Util.hpp"

#include <algorithm>
#include <iostream>
#include <functional>
#include <numeric>
#include <map>
#include <vector>
#include <utility>
constexpr auto tiny = std::numeric_limits<double>::epsilon();

template<typename T>
std::vector<int> sortWithPermutation(std::vector<T> &vec, std::function<bool(T const&, T const&)> const &less) {
	std::vector<std::pair<T, int>> paired;
	paired.reserve(vec.size());
	for (int i = 0; i < static_cast<int>(vec.size()); ++i) {
		paired.emplace_back(vec[i], i);
	}
	std::stable_sort(paired.begin(), paired.end(), [less](auto const &a, auto const &b) {
		return less(a.first, b.first);
	});
	std::vector<int> permutation(paired.size());
	for (int i = 0; i < int(permutation.size()); i++) {
		vec[i] = std::move(paired[i].first);
		permutation[i] = paired[i].second;
	}
	return permutation;
}

using Real = long double;

static Indent indent { };

struct Constants {
	Constants(std::string const &name = "C") :
			name_(name) {
		setType("double");
	}
	std::string operator()(Real constant, int simdWidth) {
		auto iterator = indexes_.find(constant);
		if (iterator == indexes_.end()) {
			indexes_.insert(std::make_pair(constant, nextIndex_++));
			iterator = indexes_.find(constant);
		}
		return name_ + std::to_string(simdWidth) + "[" + std::to_string(iterator->second) + std::string("]");
	}
//	std::string getHeader(std::string const &realTypename, int simdWidth) const {
//		std::ostringstream code;
//		code << std::string(indent) << "static const std::array<" << realTypename << ", " << std::to_string(indexes_.size()) << "> " << name_ << simdWidth
//				<< ";\n";
//		return code.str();
//	}
	std::string getCode(std::string const &realTypename, int simdWidth) const {
		std::ostringstream code;
		std::ostringstream line;
		code << std::string(indent) << "const std::array<" << realTypename << ", " << std::to_string(indexes_.size()) << "> " << name_ << simdWidth << " {";
		indent++;
		bool first = true;
		code << std::string(indent);
		int counter = 0;
		for (auto iterator = indexes_.begin(); iterator != indexes_.end(); iterator++) {
			if (!first) {
				line << ", ";
			}
			if (counter++ % 3 == 0) {
				code << line.str();
				line = std::ostringstream { };
				line << "\n" << std::string(indent);
			}
			line << realTypename << "{" << std::setprecision(std::numeric_limits<double>::max_digits10 - 1) << std::scientific << iterator->first << "}";
			first = false;
		}
		code << line.str();
		line.clear();
		indent--;
		code << "\n" << std::string(indent) + "};\n";
		return code.str();
	}
	void reset() {
		*this = Constants { };
	}
	void setType(std::string type) {
		type_ = type;
	}
private:
	std::map<double, int, std::less<double>> indexes_;
	std::string type_;
	std::string name_ = "C";
	int nextIndex_ = 0;
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
		I[n][n] = Real(1);
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
	assert(A[0].size() == B.size());
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
			asprintf(&ptr, "%.3f ", (double) A[n][m]);
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
		return std::sqrt(Real(1) / ipow(Real(1) - x * x, m)) * std::assoc_legendre(n, m, x);
	} else {
		Real value = Real(1);
		for (int k = 1; k < m; k++) {
			value *= Real(n + k) / Real(2 * k * (n - k));
		}
		return Real((x < Real(0)) ? nonepow(n + m) : +1) * value;
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
				newTheta = theta + rootEquation / rootEquationDerivative;
			} while (double(newTheta) != std::nextafter(double(theta), double(newTheta)));
			point.weight = Real(2) / (sqr(legendreP(nodeCount, point.position, 1)) * (Real(1) - sqr(point.position)));
			results.push_back(point);
		}
		break;
	}
	case Quadrature::gaussLobatto: {
		Real const baseWeight = Real(2) / Real(nodeCount * (nodeCount - 1));
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
				rootEquationDerivative = legendreP(nodeCount - 1, point.position, 2) * std::sin(theta);
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

Matrix transformMatrix1D(TransformDirection transformDirection, int modeCount, int nodeCount, Quadrature quadratureType, bool isDerivative = false) {
	Matrix transform;
	if (transformDirection == TransformDirection::forward) {
		transform = createMatrix(modeCount, nodeCount);
		auto const quadratureRules = gaussQuadrature(nodeCount, quadratureType);
		for (int modeIndex = 0; modeIndex < modeCount; modeIndex++) {
			for (int nodeIndex = 0; nodeIndex < nodeCount; nodeIndex++) {
				Real const position = quadratureRules[nodeIndex].position;
				Real const weight = quadratureRules[nodeIndex].weight;
				Real const inverseMass = Real(2 * modeIndex + 1) / Real(2);
				transform[modeIndex][nodeIndex] = weight * inverseMass * legendreP(modeIndex, position);
			}
		}
	} else {
		assert(!isDerivative);
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

std::string matrixVectorProduct(std::vector<std::string> const &v, Matrix const &A, std::vector<std::string> const &x, int simdWidth) {
	std::string code;
	int const M = x.size();
	int const N = A.size();
	if (v.size() != N) {
		printf("---- %i %i %i %i\n", v.size(), A.size(), A[0].size(), x.size());
	}
	assert((int ) v.size() == N);
	assert((int ) A[0].size() == M);

	for (int n = 0; n < N; n++) {
		if (v[n] == "") {
			continue;
		}
		bool first = true;
		char *ptr;
		for (int m = 0; m < M; m++) {
			Real const C = A[n][m];
			if (std::abs(C) < tiny) {
				continue;
			}
			std::string cons;
			if (std::abs(std::abs(C) - Real(1)) >= tiny) {
				cons = std::to_string(std::abs(C));
			}
			if (first) {
				if (std::abs(C - Real(+1)) < tiny) {
					asprintf(&ptr, "%s = %s;\n", v[n].c_str(), x[m].c_str());
				} else if (std::abs(C - Real(-1)) < tiny) {
					asprintf(&ptr, "%s = -%s;\n", v[n].c_str(), x[m].c_str());
				} else if (C < Real(0)) {
					assert(cons != "");
					asprintf(&ptr, "%s = -%s * %s;\n", v[n].c_str(), cons.c_str(), x[m].c_str());
				} else/*if (C > Real(0))*/{
					assert(cons != "");
					asprintf(&ptr, "%s = %s * %s;\n", v[n].c_str(), cons.c_str(), x[m].c_str());
				}
				first = false;
			} else {
				if (std::abs(C - Real(+1)) < tiny) {
					asprintf(&ptr, "%s += %s;\n", v[n].c_str(), x[m].c_str());
				} else if (std::abs(C - Real(-1)) < tiny) {
					asprintf(&ptr, "%s -= %s;\n", v[n].c_str(), x[m].c_str());
				} else if (C < Real(0)) {
					assert(cons != "");
					asprintf(&ptr, "%s -= %s * %s;\n", v[n].c_str(), cons.c_str(), x[m].c_str());
				} else/*if (C > Real(0))*/{
					assert(cons != "");
					asprintf(&ptr, "%s += %s * %s;\n", v[n].c_str(), cons.c_str(), x[m].c_str());
				}
			}
			code += indent + ptr;
			free(ptr);
		}
		if (first) {
			asprintf(&ptr, "%s  = %s(0);\n", v[n].c_str(), (std::string("Real") + std::to_string(simdWidth)).c_str());
			code += indent + ptr;
			free(ptr);
		}
	}
	return code;
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

std::vector<int> totalOrderPermute(int dimensionCount, int modeCount) {
	std::function<bool(std::vector<int> const&, std::vector<int> const&)> const less = [](std::vector<int> const &a, std::vector<int> const &b) {
		int aSum = std::accumulate(a.begin(), a.end(), 0);
		int bSum = std::accumulate(b.begin(), b.end(), 0);
		for (int i = 0; i < int(a.size()); i++) {
			if (aSum < bSum) {
				return true;
			} else if (aSum > bSum) {
				return false;
			}
			aSum = aSum - a[i];
			bSum = bSum - b[i];
		}
		return false;
	};
	int const flatSize = std::pow(modeCount, dimensionCount);
	std::vector<std::vector<int>> allowed;
	std::vector<int> indices(dimensionCount, 0);
	for (int index = 0; index < flatSize; index++) {
		int const deg = std::accumulate(indices.begin(), indices.end(), 0);
		if (deg < modeCount) {
			allowed.push_back(indices);
		}
		if (index + 1 < flatSize) {
			int dim = dimensionCount;
			while (++indices[--dim] == modeCount) {
				indices[dim] = 0;
			}
		}
	}
	return sortWithPermutation(allowed, less);
}

Matrix permutationToMatrix(std::vector<int> const &permutation) {
	auto P = createMatrix(permutation.size());
	int const N = P.size();
	for (int n = 0; n < N; n++) {
		P[n][permutation[n]] = Real(1);
	}
	return P;
}

Matrix transformMatrix(TransformDirection transformDirection, int transformDimension) {
	auto A = identityMatrix(1);
	Matrix B;
	std::vector<int> sizes;
	for (int thisDimension = dimensionCount - 1; thisDimension >= 0; thisDimension--) {
		int const nodeCount = modeCount + 1;
		if (thisDimension < transformDimension) {
			A = kroneckerProduct(identityMatrix(nodeCount), A);
		} else if (thisDimension > transformDimension) {
			A = kroneckerProduct(identityMatrix(modeCount), A);
		} else/*if( thisDimension == transformDimension*/{
			auto const transform = transformMatrix1D(transformDirection, modeCount, nodeCount, Quadrature::gaussLobatto);
			A = kroneckerProduct(transform, A);
		}
	}
	int const columnCount = A[0].size();
	int const rowCount = A.size();
	if (transformDirection == TransformDirection::forward) {
		auto const P1 = permutationMatrix(columnCount, modeCount, dimensionCount - transformDimension - 1);
		auto const P2 = permutationMatrix(rowCount, modeCount, dimensionCount - transformDimension);
		A = matrixMultiply(P2, matrixMultiply(A, matrixTranspose(P1)));
	} else {
		auto const P1 = permutationMatrix(columnCount, modeCount, dimensionCount - transformDimension);
		auto const P2 = permutationMatrix(rowCount, modeCount, dimensionCount - transformDimension - 1);
//		printf("%i x %i\n", P1.size(), P1[0].size());
//		printf("%i x %i\n", A.size(), A[0].size());
//		printf("%i x %i\n\n", P2.size(), P2[0].size());
		A = matrixMultiply(P2, matrixMultiply(A, matrixTranspose(P1)));
	}
//	if ((transformDimension == 0) && (dimensionCount > 1)) {
//		auto P = permutationToMatrix(totalOrderPermute(dimensionCount, modeCount));
//		A = matrixMultiply(P, A);
//	}
	return A;
}

std::vector<std::string> generateVariableNames(std::string const &name, int count, int offset = 0) {
	std::vector<std::string> names(count);
	for (int nameIndex = 0; nameIndex < count; nameIndex++) {
		names[nameIndex] = name + "[" + std::to_string(nameIndex + offset) + "]";
	}
	return names;
}

std::vector<std::string> generateVariableNames(char name, int count) {
	return generateVariableNames(std::string(1, name), count);
}

std::string generateArgumentDeclaration(char name, int count, std::string type) {
	return std::string("std::array<") + type + std::string(", ") + std::to_string(count) + std::string("> const& ") + std::string(1, name);
}

std::string generateVariableDeclaration(char name, int count, std::string type) {
	return indent + std::string("std::array<") + type + std::string(", ") + std::to_string(count) + std::string("> ") + std::string(1, name) + ";\n";
}

Matrix swapMatrix(int swapDimension, int M) {
	assert(swapDimension + 1 < dimensionCount);
	int const N = ipow(M, dimensionCount);
	auto P = createMatrix(N);
	int const stride = ipow(M, dimensionCount - 1 - swapDimension);
	for (int n = 0; n < N; n++) {
		int k = n + (1 - stride) * ((n / stride) % M - n % M);
		P[n][k] = Real(1);
	}
//	std::cout << matrixToString(P) << "\n";
	return P;
}

Matrix stiffnessMatrix(int modeCount) {
	auto S = createMatrix(modeCount);
	auto const rules = gaussQuadrature(modeCount, Quadrature::gaussLegendre);
	for (int n = 0; n < modeCount; n++) {
		for (int m = 0; m < modeCount; m++) {
			S[n][m] = Real(0);
			for (auto const &rule : rules) {
				Real const x = rule.position;
				Real const w = rule.weight;
				S[n][m] += w * legendreP(n, x, 1) * legendreP(m, x);
			}
		}
	}
	auto S3d = S;
	for (int dimension = dimensionCount - 2; dimension >= 0; dimension--) {
		S3d = kroneckerProduct(identityMatrix(modeCount), S3d);
	}
	return S3d;
}

Matrix massMatrix(int N, bool inverse) {
	auto M = createMatrix(N);
	auto const rules = gaussQuadrature(N + 1, Quadrature::gaussLegendre);
	for (int n = 0; n < N; n++) {
		for (int m = 0; m < N; m++) {
			M[n][m] = Real(0);
			for (auto const &rule : rules) {
				Real const x = rule.position;
				Real const w = rule.weight;
				M[n][m] += w * legendreP(n, x) * legendreP(m, x);
			}
		}
	}
	if (inverse) {
		for (int m = 0; m < M.size(); m++) {
			for (auto &value : M[m]) {
				if (std::abs(value) > tiny) {
					value = Real(1) / value;
				}
			}
		}
	}
	auto M3d = M;
	for (int dimension = dimensionCount - 2; dimension >= 0; dimension--) {
		M3d = kroneckerProduct(M3d, M);
	}
	return M3d;
}

int main(int, char*[]) {
	constexpr int dimensionCount = DIMENSION_COUNT;
	std::string cppCode;
	std::string hppCode;
	int nodeCount = modeCount + 1;
	std::unordered_map<int, std::string> realTypenames;
	for (int width = 1; width <= maximumSimdWidth; width *= 2) {
		realTypenames[width] = std::string("Real") + std::to_string(width);
	}
	hppCode += "#pragma once\n";
	hppCode += "\n";
	hppCode += "#include <array>\n";
	hppCode += "#include <experimental/simd>\n";
	hppCode += "\n";
	hppCode += "using Real1 = double;\n";
	for (int width = 2; width <= maximumSimdWidth; width *= 2) {
		char *ptr;
		asprintf(&ptr, "using Real%i = std::experimental::fixed_size_simd<double, %i>;\n", width, width);
		hppCode += ptr;
		free(ptr);
	}
	int const nodeCount3d = ipow(modeCount + 1, dimensionCount);
	int const modeCount3d = binco(modeCount + dimensionCount - 1, dimensionCount);
	for (int width = 1; width <= maximumSimdWidth; width *= 2) {
		char *ptr;
		asprintf(&ptr, "using ModalCoefficients%i = std::array<Real%i, %i>;\n", width, width, modeCount3d);
		hppCode += ptr;
		free(ptr);
	}
	for (int width = 1; width <= maximumSimdWidth; width *= 2) {
		char *ptr;
		asprintf(&ptr, "using NodalValues%i = std::array<Real%i, %i>;\n", width, width, nodeCount3d);
		hppCode += ptr;
		free(ptr);
	}
	hppCode += "\n";
	int simdWidth = 1;
	std::vector<Matrix> factors;
	std::vector<std::string> inputs, outputs;
	auto realTypename = realTypenames[simdWidth];
	std::string cppHeader1;
	std::string cppHeader2;
	cppHeader1 += "\n";
	cppHeader1 += "#include \"transforms.hpp\"\n";
	cppHeader1 += "\n";
	cppHeader1 += "#include <cmath>\n";
	cppHeader1 += "\n";
	hppCode += "\n";
	auto const genArray = [](std::string typeName, int count) {
		return std::string("std::array<") + typeName + ", " + std::to_string(count) + ">";
	};
	for (int simdWidth = 1; simdWidth <= maximumSimdWidth; simdWidth *= 2) {
		hppCode += genArray(realTypename, modeCount3d) + " legendreAnalyze(" + genArray(realTypename, nodeCount3d) + " const&);\n";
	}
	hppCode += "\n";
	for (int simdWidth = 1; simdWidth <= maximumSimdWidth; simdWidth *= 2) {
		hppCode += genArray(realTypename, nodeCount3d) + " legendreSynthesize(" + genArray(realTypename, modeCount3d) + " const&);\n";
	}
	hppCode += "\n";
	for (int simdWidth = 1; simdWidth <= maximumSimdWidth; simdWidth *= 2) {
		hppCode += genArray(realTypename, modeCount3d) + " applyInverseMass(" + genArray(realTypename, modeCount3d) + " const&);\n";
	}
	hppCode += "\n";
	for (int simdWidth = 1; simdWidth <= maximumSimdWidth; simdWidth *= 2) {
		hppCode += genArray(realTypename, modeCount3d) + " applyStiffness(" + genArray(realTypename, modeCount3d) + " const&);\n";
	}
	hppCode += "\n";
	std::string auxCode;
	for (int simdWidth = 1; simdWidth <= maximumSimdWidth; simdWidth *= 2) {
		auto realTypename = realTypenames[simdWidth];
		auto M = massMatrix(modeCount, true);
		auto permute = permutationMatrix(M[0].size(), modeCount, dimensionCount);
		M = matrixMultiply(M, matrixTranspose(permute));
		M = matrixMultiply(permute, M);
		auxCode += indent + genArray(realTypename, modeCount3d) + " applyInverseMass(" + genArray(realTypename, modeCount3d) + " const& input) {\n";
		indent++;
		inputs = generateVariableNames("input", M[0].size());
		outputs = generateVariableNames("output", M.size());
		auxCode += indent + genArray(realTypename, M.size()) + " output;\n";
		auxCode += matrixVectorProduct(outputs, M, inputs, simdWidth);
		auxCode += indent + "return output;\n";
		indent--;
		auxCode += indent + "}\n\n";

		auto S = stiffnessMatrix(modeCount);
		permute = permutationMatrix(S[0].size(), modeCount, dimensionCount);
		S = matrixMultiply(S, matrixTranspose(permute));
		S = matrixMultiply(permute, S);
		auxCode += indent + genArray(realTypename, modeCount3d) + " applyStiffness(" + genArray(realTypename, modeCount3d) + " const& input) {\n";
		indent++;
		inputs = generateVariableNames("input", S[0].size());
		outputs = generateVariableNames("output", S.size());
		auxCode += indent + genArray(realTypename, S.size()) + " output;\n";
		auxCode += matrixVectorProduct(outputs, S, inputs, simdWidth);
		auxCode += indent + "return output;\n";
		indent--;
		auxCode += indent + "}\n\n";
		for (int dim = 0; dim < dimensionCount - 1; dim++) {
			auto P = swapMatrix(dim, modeCount);
			P = matrixMultiply(permute, P);
			permute = permutationMatrix(P[0].size(), modeCount, dimensionCount);
			P = matrixMultiply(P, matrixTranspose(permute));
			auxCode +=
					indent + genArray(realTypename, modeCount3d) + " applySwap" + std::string(1, 'X' + dim) + "Z(" + genArray(realTypename, modeCount3d) + " const& input) {\n";
			indent++;
			inputs = generateVariableNames("input", P[0].size());
			outputs = generateVariableNames("output", P.size());
			auxCode += indent + genArray(realTypename, P.size()) + " output;\n";
			auxCode += matrixVectorProduct(outputs, P, inputs, simdWidth);
			auxCode += indent + "return output;\n";
			indent--;
			auxCode += indent + "}\n\n";
		}
		for (int dim = 0; dim < dimensionCount - 1; dim++) {
			auto P = swapMatrix(dim, nodeCount);
			auxCode +=
					indent + genArray(realTypename, nodeCount3d) + " applySwap" + std::string(1, 'X' + dim) + "Z(" + genArray(realTypename, nodeCount3d) + " const& input) {\n";
			indent++;
			inputs = generateVariableNames("input", P[0].size());
			outputs = generateVariableNames("output", P.size());
			auxCode += indent + genArray(realTypename, P.size()) + " output;\n";
			auxCode += matrixVectorProduct(outputs, P, inputs, simdWidth);
			auxCode += indent + "return output;\n";
			indent--;
			auxCode += indent + "}\n\n";
		}
	}
	cppCode += "\n\n" + auxCode;
	auxCode.clear();
	for (int simdWidth = 1; simdWidth <= maximumSimdWidth; simdWidth *= 2) {
		inputs = decltype(inputs)();
		outputs = decltype(outputs)();
		auto const realTypename = realTypenames[simdWidth];
		bool isAnalyze;
		isAnalyze = true;
		std::string functionName = "legendreAnalyze";
		int inCount = nodeCount3d;
		int outCount = modeCount3d;
SYNTHESIZE:
		factors.clear();
		if (!isAnalyze) {
			functionName = "legendreSynthesize";
		}
		outputs = generateVariableNames("input", inCount);
		cppCode += genArray(realTypename, outCount) + " " + functionName + "(" + genArray(realTypename, inCount) + " const& input) {\n";
		std::string deferredCode;
		indent++;
		std::vector<int> arraySizes;
		arraySizes = std::vector<int>(1, inCount);
		factors.clear();
		auto const direction = isAnalyze ? TransformDirection::forward : TransformDirection::backward;
		for (int transformDimension = 0; transformDimension < dimensionCount; transformDimension++) {
			factors.push_back(transformMatrix(direction, transformDimension));
		}
		if (isAnalyze) {
			std::reverse(factors.begin(), factors.end());
		}
		for (int transformDimension = 0; transformDimension < dimensionCount; transformDimension++) {
			auto const sz = factors[transformDimension].size();
			arraySizes.push_back(sz);
		}
		int bufferSize, currentSize;
		bufferSize = *(std::max_element(arraySizes.begin() + 1, arraySizes.end() - 1));
		currentSize = arraySizes[0] + arraySizes[1];
		for (int i = 0; i < dimensionCount - 2; i++) {
			currentSize -= arraySizes[i];
			currentSize += arraySizes[i + 2];
			bufferSize = std::max(bufferSize, currentSize);
		}
		cppCode += indent + genArray(realTypename, bufferSize) + " buffer;\n";
		int bufferOffset;
		bufferOffset = 0;
		for (int transformDimension = 0; transformDimension < dimensionCount; transformDimension++) {
			auto const &A = factors[transformDimension];
			inputs = std::move(outputs);
			if (transformDimension != (dimensionCount - 1)) {
				outputs = generateVariableNames("buffer", A.size(), bufferOffset);
				if (bufferOffset == 0) {
					bufferOffset = A.size();
				} else {
					bufferOffset = 0;
				}
			} else {
				outputs = generateVariableNames("output", outCount);
				cppCode += indent + genArray(realTypename, A.size()) + " output;\n";
			}
			deferredCode += matrixVectorProduct(outputs, A, inputs, simdWidth);
		}
		cppCode += deferredCode;
		deferredCode.clear();
		cppCode += indent + "return output;\n";
		indent--;
		cppCode += indent + "}\n\n";
		if (isAnalyze) {
			isAnalyze = false;
			std::swap(inCount, outCount);
			goto SYNTHESIZE;
		}
	}
	for (int simdWidth = 1; simdWidth <= maximumSimdWidth; simdWidth *= 2) {
		cppHeader2 += std::string(getConstant.getCode(realTypename, simdWidth));
		cppHeader2 += "\n";
	}
	cppCode = cppHeader1 + "\n" + cppHeader2 + "\n" + cppCode + "\n";
	hppCode += "\n";

	toFile(hppCode, "./generated_source/transforms.hpp");
	toFile(cppCode, "./generated_source/transforms.cpp");
	getConstant.reset();
	return 0;
}
