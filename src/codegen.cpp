#include "Definitions.hpp"
#include "Indent.hpp"
#include "Util.hpp"

#include <algorithm>
#include <iostream>
#include <functional>
#include <numeric>
#include <map>
#include <set>
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
		setType("T");
	}
	std::string operator()(Real constant) {
		auto iterator = indexes_.find(constant);
		if (iterator == indexes_.end()) {
			indexes_.insert(std::make_pair(constant, nextIndex_++));
			iterator = indexes_.find(constant);
			values_.push_back(constant);
		}
		return name_ + std::to_string(iterator->second);
	}
//	std::string getHeader(std::string const &realTypename, int simdWidth) const {
//		std::ostringstream code;
//		code << std::string(indent) << "static const std::array<" << realTypename << ", " << std::to_string(indexes_.size()) << "> " << name_ << simdWidth
//				<< ";\n";
//		return code.str();
//	}
	std::string getCode() const {
		std::ostringstream code;
		std::ostringstream line;
		for (int i = 0; i < values_.size(); i++) {
			code << indent + "static const T " + name_ + std::to_string(i) + " = ";
			code << "T(" << std::setprecision(std::numeric_limits<double>::max_digits10 - 1) << std::scientific << values_[i] << ");\n";
		}
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
	std::vector<double> values_;
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
			asprintf(&ptr, "%6.3f ", (double) A[n][m]);
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
		int const dOrder = isDerivative ? 1 : 0;
		transform = createMatrix(modeCount, nodeCount);
		auto const quadratureRules = gaussQuadrature(nodeCount, quadratureType);
		for (int modeIndex = 0; modeIndex < modeCount; modeIndex++) {
			for (int nodeIndex = 0; nodeIndex < nodeCount; nodeIndex++) {
				Real const position = quadratureRules[nodeIndex].position;
				Real const weight = quadratureRules[nodeIndex].weight;
				transform[modeIndex][nodeIndex] = weight /* * inverseMass*/* legendreP(modeIndex, position, dOrder);
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

std::string matrixVectorProduct(std::vector<std::string> const &v, Matrix const &A, std::vector<std::string> const &x) {
	std::string code;
	int const M = x.size();
	int const N = A.size();
	if (v.size() != N) {
		std::cout << matrixToString(A) << "\n";
		printf("---- %i %i %i %i\n", v.size(), A.size(), A[0].size(), x.size());
	}
	assert((int ) v.size() == N);
	assert((int ) A[0].size() == M);

	for (int n = 0; n < N; n++) {
		if (v[n] == "") {
			continue;
		}
		code += std::string(indent);
		struct Less {
			bool operator()(Real a, Real b) const {
				if (std::abs(a - b) < tiny) {
					return false;
				} else {
					return a < b;
				}
			}
		};
		std::map<Real, std::vector<std::pair<signed char, std::string>>, Less> terms;
		std::set<Real, Less> coefficients;
		for (int m = 0; m < M; m++) {
			Real const C = A[n][m];
			if (std::abs(C) < tiny) {
				continue;
			}
			if (C > Real(0)) {
				terms[std::abs(C)].push_back( { +1, x[m] });
			} else {
				terms[std::abs(C)].push_back( { -1, x[m] });
			}
			coefficients.insert(std::abs(C));
		}
		code += v[n] + " = ";
		if (coefficients.size()) {
			bool first1 = true;
			for (auto it = coefficients.begin(); it != coefficients.end(); it++) {
				Real const C = *it;
				std::cerr << std::to_string(C) << " ";
				auto const &theseTerms = terms[C];
				if (first1) {
					first1 = false;
				} else {
					code += " + ";
				}
				bool const isOne = std::abs(C - Real(1)) < tiny;
				if (theseTerms.size() > 1) {
					if (!isOne) {
						code += getConstant(C) + " * (";
					}
					bool first2 = true;
					for (auto const &term : terms[C]) {
						auto const s = term.first;
						auto const v = term.second;
						if (first2) {
							if (s < Real(0)) {
								code += "-";
							}
							first2 = false;
						} else {
							code += (s > 0) ? std::string(" + ") : std::string(" - ");
						}
						code += v;
					}
					if (!isOne) {
						code += ")";
					}
				} else {
					auto const s = theseTerms[0].first;
					auto const v = theseTerms[0].second;
					if (s < Real(0)) {
						code += "-";
					}
					if (!isOne) {
						code += getConstant(C) + " * ";
					}
					code += v;
				}
			}
			code += ";\n";
		} else {
			code += "T(0);\n";
		}
		std::cerr << "\n";
		std::string const str = "+ -";
		size_t i = code.find(str);
		while(i != std::string::npos) {
			code.replace(i, str.size(), "- ");
			i = code.find(str, i);
		}
	}
	return code;
}

Matrix permutationMatrix(int N, int modeCount, int triIndexCount) {
	std::vector<int> ones(N, -1);
	int M = 0;
	for (int targetDeg = 0; targetDeg < modeCount; targetDeg++) {
		for (int i = 0; i < N; i++) {
			int deg = 0;
			int k = i;
			for (int j = 0; j < triIndexCount; j++) {
				deg += k % modeCount;
				k /= modeCount;
			}
			if (deg == targetDeg) {
				ones[i] = M++;
			}
		}
	}
	Matrix P = createMatrix(M, N);
	for (int n = 0; n < N; n++) {
		if (ones[n] >= 0) {
			P[ones[n]][n] = Real(1.0);
		}
	}
	return P;
}

Matrix traceMatrix(int D, int modeCount, int traceFace, bool inverse = false) {
	int const traceDim = traceFace >> 1;
	int const traceSign = 2 * (traceFace & 1) - 1;
	std::vector<std::vector<int>> from, to;
	int bigSize = binco(modeCount + D - 1, D);
	int smallSize = binco(modeCount + D - 2, D - 1);
	for (int targetDeg = 0; targetDeg < modeCount; targetDeg++) {
		for (int i = 0; i < ipow(modeCount, D); i++) {
			int deg = 0;
			int k = i;
			std::vector<int> indices(D);
			for (int j = 0; j < D; j++) {
				indices[D - j - 1] = k % modeCount;
				deg += k % modeCount;
				k /= modeCount;
			}
			if (deg == targetDeg) {
				from.push_back(indices);
			}
		}
		for (int i = 0; i < ipow(modeCount, D - 1); i++) {
			int deg = 0;
			int k = i;
			std::vector<int> indices(D - 1);
			for (int j = 0; j < D - 1; j++) {
				indices[D - j - 2] = k % modeCount;
				deg += k % modeCount;
				k /= modeCount;
			}
			if (deg == targetDeg) {
				to.push_back(indices);
			}
		}
	}
	//printf("%i %i\n", to.size(), smallSize);
	assert(to.size() == smallSize);
	assert(from.size() == bigSize);
	Matrix P = createMatrix(smallSize, bigSize);
	for (int n = 0; n < smallSize; n++) {
		for (int m = 0; m < bigSize; m++) {
			bool flag = true;
			for (int d = 0; d < D - 1; d++) {
				if ((to[n][d] != from[m][d + int(d >= traceDim)])) {
					flag = false;
				}
			}
			if (flag) {
				int sign = (from[m][traceDim] & 1) ? traceSign : 1;
				P[n][m] = sign;
			}
		}
	}
	if (inverse) {
		return matrixTranspose(P);
	} else {
		return P;
	}
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

template<Quadrature quadratureType = Quadrature::gaussLegendre>
Matrix transformMatrix(int dimensionCount, int modeCount, TransformDirection transformDirection, int transformDimension, int derivDimension = -1) {
	auto A = identityMatrix(1);
	Matrix B;
	std::vector<int> sizes;
	for (int thisDimension = dimensionCount - 1; thisDimension >= 0; thisDimension--) {
		int const nodeCount = modeCount + int(quadratureType == Quadrature::gaussLobatto);
		if (thisDimension < transformDimension) {
			A = kroneckerProduct(identityMatrix(nodeCount), A);
		} else if (thisDimension > transformDimension) {
			A = kroneckerProduct(identityMatrix(modeCount), A);
		} else/*if( thisDimension == transformDimension*/{
			auto const transform = transformMatrix1D(transformDirection, modeCount, nodeCount, quadratureType, derivDimension == transformDimension);
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
		A = matrixMultiply(P2, matrixMultiply(A, matrixTranspose(P1)));
	}
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

Matrix stiffnessMatrix(int dimensionCount, int modeCount, int derivDimension) {
	auto S = createMatrix(modeCount);
	auto const rules = gaussQuadrature(modeCount + 1, Quadrature::gaussLegendre);
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
	auto S3d = identityMatrix(1);
	;
	for (int dimension = dimensionCount - 1; dimension >= 0; dimension--) {
		S3d = kroneckerProduct(dimension == derivDimension ? S : identityMatrix(modeCount), S3d);
	}
	auto const P1 = permutationMatrix(S3d[0].size(), modeCount, dimensionCount);
	S3d = matrixMultiply(P1, matrixMultiply(S3d, matrixTranspose(P1)));
	return S3d;
}

Matrix massMatrix(int dimensionCount, int N, bool inverse) {
	auto M = createMatrix(N);
	auto const rules = gaussQuadrature(N + 3, Quadrature::gaussLegendre);
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
	auto const P1 = permutationMatrix(M3d[0].size(), N, dimensionCount);
	M3d = matrixMultiply(P1, matrixMultiply(M3d, matrixTranspose(P1)));
	return M3d;
}
std::string tag(int D, int O) {
	return std::string("_") + std::to_string(D) + "D_O" + std::to_string(O);
}

std::string generate(int dimensionCount, int modeCount) {
	std::string hppCode;
	int nodeCount = modeCount;
	int const nodeCount3d = ipow(nodeCount, dimensionCount);
	int const modeCount3d = binco(modeCount + dimensionCount - 1, dimensionCount);
	std::vector<Matrix> factors;
	std::vector<std::string> inputs, outputs;
	auto const genArray = [](std::string typeName, int count) {
		return std::string("std::array<") + typeName + ", " + std::to_string(count) + ">";
	};
	hppCode += "\n";
	inputs = decltype(inputs)();
	outputs = decltype(outputs)();
	bool isAnalyze;
	isAnalyze = true;
	std::string functionTag = tag(dimensionCount, modeCount);
	std::string functionName = "dgAnalyze" + functionTag;
	int inCount = nodeCount3d;
	int outCount = modeCount3d;
SYNTHESIZE:
	factors.clear();
	if (!isAnalyze) {
		functionName = "dgSynthesize" + functionTag;
	}
	hppCode += "template<typename T>\n";
	outputs = generateVariableNames("input", inCount);
	if (isAnalyze) {
		hppCode += "std::array<T, triangleSize<" + std::to_string(dimensionCount) + ", " + std::to_string(modeCount) + ">> ";
		hppCode += functionName;
		hppCode += "(std::array<T, squareSize<" + std::to_string(dimensionCount) + ", " + std::to_string(modeCount) + ">> const& input) {\n";
	} else {
		hppCode += "std::array<T, squareSize<" + std::to_string(dimensionCount) + ", " + std::to_string(modeCount) + ">> ";
		hppCode += functionName;
		hppCode += "(std::array<T, triangleSize<" + std::to_string(dimensionCount) + ", " + std::to_string(modeCount) + ">> const& input) {\n";
	}
	std::string deferredCode;
	indent++;
	std::vector<int> arraySizes;
	arraySizes = std::vector<int>(1, inCount);
	factors.clear();
	auto direction = isAnalyze ? TransformDirection::forward : TransformDirection::backward;
	for (int transformDimension = 0; transformDimension < dimensionCount; transformDimension++) {
		factors.push_back(transformMatrix(dimensionCount, modeCount, direction, transformDimension));
//		std::cout << matrixToString(factors.back()) << std::endl;
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
	hppCode += indent + genArray("T", bufferSize) + " buffer;\n";
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
			hppCode += indent + genArray("T", A.size()) + " output;\n";
		}
		deferredCode += matrixVectorProduct(outputs, A, inputs);
	}
	hppCode += std::string(getConstant.getCode());
	getConstant.reset();
	hppCode += deferredCode;
	deferredCode.clear();
	hppCode += indent + "return output;\n";
	indent--;
	hppCode += indent + "}\n\n";
	if (isAnalyze) {
		isAnalyze = false;
		std::swap(inCount, outCount);
		goto SYNTHESIZE;
	}

	return hppCode;
}

std::string generateGaussLobattoSynthesize(int dimensionCount, int modeCount) {
	std::string hppCode;
	int nodeCount = modeCount + 1;
	int const nodeCount3d = ipow(nodeCount, dimensionCount);
	int const modeCount3d = binco(modeCount + dimensionCount - 1, dimensionCount);
	std::vector<Matrix> factors;
	std::vector<std::string> inputs, outputs;
	auto const genArray = [](std::string typeName, int count) {
		return std::string("std::array<") + typeName + ", " + std::to_string(count) + ">";
	};
	hppCode += "\n";
	inputs = decltype(inputs)();
	outputs = decltype(outputs)();
	std::string functionTag = tag(dimensionCount, modeCount);
	std::string functionName = "dgAnalyze" + functionTag;
	int inCount = modeCount3d;
	int outCount = nodeCount3d;
	factors.clear();
	functionName = "dgSynthesizeGaussLobatto" + functionTag;
	hppCode += "template<typename T>\n";
	outputs = generateVariableNames("input", inCount);
	hppCode += "std::array<T, squareSize<" + std::to_string(dimensionCount) + ", " + std::to_string(nodeCount) + ">> ";
	hppCode += functionName;
	hppCode += "(std::array<T, triangleSize<" + std::to_string(dimensionCount) + ", " + std::to_string(modeCount) + ">> const& input) {\n";
	std::string deferredCode;
	indent++;
	std::vector<int> arraySizes;
	arraySizes = std::vector<int>(1, inCount);
	factors.clear();
	for (int transformDimension = 0; transformDimension < dimensionCount; transformDimension++) {
		factors.push_back(transformMatrix<Quadrature::gaussLobatto>(dimensionCount, modeCount, TransformDirection::backward, transformDimension));
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
	hppCode += indent + genArray("T", bufferSize) + " buffer;\n";
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
			hppCode += indent + genArray("T", A.size()) + " output;\n";
		}
//		std::cout << functionTag << "\n";
		deferredCode += matrixVectorProduct(outputs, A, inputs);
	}
	hppCode += std::string(getConstant.getCode());
	getConstant.reset();
	hppCode += deferredCode;
	deferredCode.clear();
	hppCode += indent + "return output;\n";
	indent--;
	hppCode += indent + "}\n\n";
	return hppCode;
}

std::string genMassMatrix(int dimensionCount, int modeCount) {
	std::string hppCode, code1, code2;
	int const size = binco(modeCount + dimensionCount - 1, dimensionCount);
	std::string arrayType = "std::array<T, triangleSize<" + std::to_string(dimensionCount) + ", " + std::to_string(modeCount) + ">>";
	code1 += "template<typename T>\n";
	auto inputs = generateVariableNames("input", size);
	auto outputs = generateVariableNames("output", size);
	code1 += indent + arrayType + " ";
	code1 += "dgMassInverse" + tag(dimensionCount, modeCount);
	code1 += "(" + arrayType + " const& input) {\n";
	indent++;
	code2 += indent + arrayType + " output;\n";
	auto A = massMatrix(dimensionCount, modeCount, true);
	code2 += matrixVectorProduct(outputs, A, inputs);
	hppCode += code1 + std::string(getConstant.getCode()) + code2;
	getConstant.reset();
	hppCode += indent + "return output;\n";
	indent--;
	hppCode += "}\n";
	return hppCode;
}

std::string genStiffnessMatrix(int dimensionCount, int modeCount) {
	std::string hppCode, code1, code2;
	int const size = binco(modeCount + dimensionCount - 1, dimensionCount);
	std::string arrayType = "std::array<T, triangleSize<" + std::to_string(dimensionCount) + ", " + std::to_string(modeCount) + ">>";
	code1 += "\n";
	code1 += "template<typename T>\n";
	auto inputs = generateVariableNames("input", size);
	auto outputs = generateVariableNames("output", size);
	code1 += indent + arrayType + " ";
	code1 += "dgStiffness" + tag(dimensionCount, modeCount);
	code1 += "(int dimension, " + arrayType + " const& input) {\n";
	indent++;
	code2 += indent + arrayType + " output;\n";
	if (dimensionCount > 1) {
		code2 += std::string(indent);
	}
	for (int dim = 0; dim < dimensionCount; dim++) {
		if (dimensionCount > 1) {
			if (dim == dimensionCount - 1) {
				code2 += "/*";
			}
			code2 += "if(dimension == " + std::to_string(dim) + ")";
			if (dim == dimensionCount - 1) {
				code2 += "*/";
			}
			code2 += " {\n";
			indent++;
		}
		auto A = stiffnessMatrix(dimensionCount, modeCount, dim);
		code2 += matrixVectorProduct(outputs, A, inputs);
		if (dimensionCount > 1) {
			indent--;
			code2 += indent + "}";
			if (dim + 1 < dimensionCount) {
				code2 += " else ";
			} else {
				code2 += "\n";
			}
		}
	}
	hppCode += code1 + std::string(getConstant.getCode()) + code2;

	getConstant.reset();
	hppCode += indent + "return output;\n";
	indent--;
	hppCode += "}\n";
	return hppCode;
}
//Matrix traceMatrix(int D, int modeCount, int traceDim, bool inverse = false) {

std::string genTrace(int dimensionCount, int modeCount, bool inverse) {
	std::string hppCode, code1, code2;
	int size1 = binco(modeCount + dimensionCount - 1, dimensionCount);
	int size2 = binco(modeCount + dimensionCount - 2, dimensionCount - 1);
	if (inverse) {
		std::swap(size1, size2);
	}
	std::string arrayType2 = "std::array<T, triangleSize<" + std::to_string(dimensionCount - int(inverse)) + ", " + std::to_string(modeCount) + ">>";
	std::string arrayType1 = "std::array<T, triangleSize<" + std::to_string(dimensionCount + int(inverse) - 1) + ", " + std::to_string(modeCount) + ">>";
	code1 += "\ntemplate<typename T>\n";
	auto inputs = generateVariableNames("input", size1);
	auto outputs = generateVariableNames("output", size2);
	code1 += indent + arrayType1 + " ";
	code1 += "dgTrace" + std::string(inverse ? "Inverse" : "");
	code1 += tag(dimensionCount, modeCount);
	code1 += "(int face, " + arrayType2 + " const& input) {\n";
	indent++;
	code2 += indent + arrayType1 + " output;\n";
	code2 += std::string(indent);
	for (int face = 0; face < 2 * dimensionCount; face++) {
		if (face == 2 * dimensionCount - 1) {
			code2 += "/*";
		}
		code2 += "if(face == " + std::to_string(face) + ")";
		if (face == 2 * dimensionCount - 1) {
			code2 += "*/";
		}
		code2 += " {\n";
		indent++;
		auto A = traceMatrix(dimensionCount, modeCount, face, inverse);
		code2 += matrixVectorProduct(outputs, A, inputs);
		indent--;
		code2 += indent + "}";
		if (face + 1 < 2 * dimensionCount) {
			code2 += " else ";
		} else {
			code2 += "\n";
		}
	}
	hppCode += code1 + std::string(getConstant.getCode()) + code2;
	getConstant.reset();
	hppCode += indent + "return output;\n";
	indent--;
	hppCode += "}\n";
	return hppCode;
}

int main(int, char*[]) {
	std::string hppCode;
	hppCode += "#pragma once\n";
	hppCode += "\n";
	hppCode += "#include \"Util.hpp\"\n";
	hppCode += "#include <array>\n";
	hppCode += "#include <cmath>\n";
	hppCode += "\n";
	hppCode += "enum class Quadrature : int {\n";
	hppCode += "	gaussLegendre, gaussLobatto\n";
	hppCode += "};\n";
	hppCode += "\n";
	hppCode += indent + "template<int D, int O>\n";
	hppCode += indent + "constexpr int triangleSize = binco(O + D - 1, D);\n";
	hppCode += indent + "\n";
	hppCode += indent + "template<int D, int O>\n";
	hppCode += indent + "constexpr int squareSize = ipow(O, D);\n";
	hppCode += indent + "\n";
	for (int dim = 1; dim <= 3; dim++) {
		for (int order = 1; order <= 4; order++) {
			hppCode += generate(dim, order);
			hppCode += genMassMatrix(dim, order);
			hppCode += genStiffnessMatrix(dim, order);
			hppCode += genTrace(dim, order, false);
			hppCode += genTrace(dim, order, true);
//			hppCode += generateGaussLobattoSynthesize(dim, order);
		}
	}
	for (int iter = 0; iter <= 5; iter++) {
		hppCode += indent + "\n";
//		if (iter == 1) {
//			hppCode += indent + "template<typename T, int D, int O, Quadrature Q = Quadrature::gaussLegendre>\n";
//		} else {
		hppCode += indent + "template<typename T, int D, int O>\n";
//		}
		hppCode += std::string(indent);
		std::string fname, varname;
		if (iter == 0) {
			fname = "dgAnalyze";
			hppCode += "auto " + fname + "(std::array<T, squareSize<D, O>> const& input) {\n";
		} else if (iter == 1) {
			fname = "dgSynthesize";
			hppCode += "auto " + fname + "(std::array<T, triangleSize<D, O>> const& input) {\n";
		} else if (iter == 2) {
			fname = "dgMassInverse";
			hppCode += "auto " + fname + "(std::array<T, triangleSize<D, O>> const& input) {\n";
		} else if (iter == 3) {
			fname = "dgStiffness";
			varname = "dimension, ";
			hppCode += "auto " + fname + "(int dimension, std::array<T, triangleSize<D, O>> const& input) {\n";
		} else if (iter == 4) {
			varname = "face, ";
			fname = "dgTrace";
			hppCode += "auto " + fname + "(int face, std::array<T, triangleSize<D, O>> const& input) {\n";
		} else if (iter == 5) {
			varname = "face, ";
			fname = "dgTraceInverse";
			hppCode += "auto " + fname + "(int face, std::array<T, triangleSize<D - 1, O>> const& input) {\n";
		}
		indent++;
		for (int dim = 1; dim <= 3; dim++) {
			if (dim == 1) {
				hppCode += std::string(indent);
			}
			hppCode += "if constexpr(D == " + std::to_string(dim) + ") {\n";
			indent++;
			for (int order = 1; order <= 4; order++) {
				if (order == 1) {
					hppCode += std::string(indent);
				}
				hppCode += "if constexpr(O == " + std::to_string(order) + ") {\n";
				indent++;
				if (fname == "dgSynthesize") {
//					hppCode += indent + "if constexpr(Q == Quadrature::gaussLobatto) {\n";
//					indent++;
//					hppCode += indent + "return " + fname + "GaussLobatto" + tag(dim, order) + "(" + varname + "input);\n";
//					indent--;
//					hppCode += indent + "} else /*if constexpr(Q == Quadrature::gaussLegendre)*/ {\n";
//					indent++;
					hppCode += indent + "return " + fname + tag(dim, order) + "(" + varname + "input);\n";
//					indent--;
//					hppCode += indent + "}\n";
				} else {
					hppCode += indent + "return " + fname + tag(dim, order) + "(" + varname + "input);\n";
				}
				indent--;
				hppCode += indent + "}";
				if (order < 4) {
					hppCode += " else ";
				} else {
					hppCode += "\n";
				}
			}
			indent--;
			hppCode += indent + "}";
			if (dim < 3) {
				hppCode += " else ";
			} else {
				hppCode += "\n";
			}
		}
		indent--;
		hppCode += indent + "}\n";
	}
	hppCode += "\n";
	hppCode += indent + "template<typename T, int D, int O>\n";
	hppCode += "std::array<T, D> getQuadraturePoint(int flatIndex) {\n";
	indent++;
	for (int dim = 1; dim <= 3; dim++) {
		if (dim == 1) {
			hppCode += std::string(indent);
		}
		hppCode += "if constexpr(D == " + std::to_string(dim) + ") {\n";
		indent++;
		for (int order = 1; order <= 4; order++) {
			if (order == 1) {
				hppCode += std::string(indent);
			}
			auto Q = gaussQuadrature(order, Quadrature::gaussLegendre);
			hppCode += "if constexpr(O == " + std::to_string(order) + ") {\n";
			indent++;
			bool first = true;
			hppCode += indent + "constexpr std::array<std::array<T, D>, squareSize<D, O>> map{{\n";
			indent++;
			hppCode += std::string(indent);
			int cnt = 0;
			for (int iii = 0; iii < ipow(order, dim); iii++) {
				std::vector<int> I;
				int k = iii;
				for (int d = 0; d < dim; d++) {
					I.push_back(k % order);
					k /= order;
				}
				std::reverse(I.begin(), I.end());
				if (first) {
					first = false;
				} else {
					hppCode += ", ";
				}
				if (cnt == 1) {
					hppCode += "\n" + std::string(indent);
					cnt = 0;
				}
				hppCode += "{";
				for (int d = 0; d < dim; d++) {
					char *ptr;
					if (std::abs(Q[I[d]].position) < 1e-14) {
						Q[I[d]].position = 0.0;
					}
					asprintf(&ptr, "T(%24.17e)", double(Q[I[d]].position));
					hppCode += ptr;
					free(ptr);
					if (d + 1 < dim) {
						hppCode += ", ";
					}
				}
				hppCode += "}";
				cnt++;
			}
			indent--;
			hppCode += "\n" + std::string(indent) + "}};\n";
			hppCode += indent + "return map[flatIndex];\n";
			indent--;
			hppCode += indent + "}";
			if (order < 4) {
				hppCode += " else ";
			} else {
				hppCode += "\n";
			}
		}
		indent--;
		hppCode += indent + "}";
		if (dim < 3) {
			hppCode += " else ";
		} else {
			hppCode += "\n";
		}
	}
	indent--;
	hppCode += indent + "}\n";
	hppCode += "\n";

	hppCode += indent + "template<int D, int O>\n";
	hppCode += "std::array<int, D> flatToTriangular(int flatIndex) {\n";
	indent++;
	for (int dim = 1; dim <= 3; dim++) {
		if (dim == 1) {
			hppCode += std::string(indent);
		}
		hppCode += "if constexpr(D == " + std::to_string(dim) + ") {\n";
		indent++;
		for (int order = 1; order <= 4; order++) {
			if (order == 1) {
				hppCode += std::string(indent);
			}
			hppCode += "if constexpr(O == " + std::to_string(order) + ") {\n";
			indent++;
			bool first = true;
			hppCode += indent + "constexpr std::array<std::array<int, D>, triangleSize<D, O>> map{{\n";
			indent++;
			hppCode += std::string(indent);
			int cnt = 0;
			for (int targetDeg = 0; targetDeg < order; targetDeg++) {
				for (int i = 0; i < ipow(order, dim); i++) {
					int deg = 0;
					int k = i;
					std::vector<int> indices(dim);
					for (int j = 0; j < dim; j++) {
						indices[dim - j - 1] = k % order;
						deg += k % order;
						k /= order;
					}
					if (deg == targetDeg) {
						if (first) {
							first = false;
						} else {
							hppCode += ", ";
						}
						if (cnt == 4) {
							hppCode += "\n" + std::string(indent);
							cnt = 0;
						}
						hppCode += "{";
						for (int d = 0; d < dim; d++) {
							hppCode += std::to_string(indices[d]);
							if (d + 1 < dim) {
								hppCode += ", ";
							}
						}
						hppCode += "}";
						cnt++;
					}
				}
			}
			indent--;
			hppCode += "\n" + std::string(indent) + "}};\n";
			hppCode += indent + "return map[flatIndex];\n";
			indent--;
			hppCode += indent + "}";
			if (order < 4) {
				hppCode += " else ";
			} else {
				hppCode += "\n";
			}
		}
		indent--;
		hppCode += indent + "}";
		if (dim < 3) {
			hppCode += " else ";
		} else {
			hppCode += "\n";
		}
	}
	indent--;
	hppCode += indent + "}\n";
	hppCode += "\n";
	hppCode += "template<int D, int O>\n"
			"int triangularToFlat(std::array<int, D> ti) {\n"
			"\tstd::array<int, D> num, den;\n"
			"\tfor (int d = D - 1; d > 0; d--) {\n"
			"\t\tti[d - 1] += ti[d];\n"
			"\t}\n"
			"\tnum = ti;\n"
			"\tden.fill(1);\n"
			"\tfor (int d1 = 1; d1 < D; d1++) {\n"
			"\t\tfor (int d2 = 0; d2 < D - d1; d2++) {\n"
			"\t\t\tnum[d2] *= ti[d2] + d1;\n"
			"\t\t\tden[d2] *= d1 + 1;\n"
			"\t\t}\n"
			"\t}\n"
			"\tint flat = 0;\n"
			"\tfor (int d = 0; d < D; d++) {\n"
			"\t\tflat += num[d] / den[d];\n"
			"\t}\n"
			"\treturn flat;\n"
			"}\n"
			"";

	toFile(hppCode, "./generated_source/dgTransforms.hpp");
	getConstant.reset();
	return 0;
}

