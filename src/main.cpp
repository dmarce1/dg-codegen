#include <hpx/hpx_init.hpp>
#include "transforms.hpp"
#include "Options.hpp"
#include "MultiIndex.hpp"
#include "HyperGrid.hpp"
#include "EulerState.hpp"
#include "Real.hpp"
#include "RungeKutta.hpp"

template<typename Type, int basisOrder, TransformDirection transformDirection>
constexpr auto fourierLegendreTransform() {
	if constexpr (transformDirection == TransformDirection::forward) {
		return matrixInverse(fourierLegendreTransform<Type, basisOrder, TransformDirection::backward>());
	} else {
		using namespace Math;
		SquareMatrix<Type, basisOrder> transform;
		for (int basisIndex = 0; basisIndex < basisOrder; basisIndex++) {
			for (int quadratureIndex = 0; quadratureIndex < basisOrder; quadratureIndex++) {
				auto const [quadraturePosition, quadratureWeight] = gaussLegendreQuadraturePoint<Type, basisOrder>(quadratureIndex);
				transform(quadratureIndex, basisIndex) = legendrePolynomial(basisIndex, quadraturePosition);
			}
		}
		return transform;
	}
}

template<typename Type, int basisOrder>
auto fourierLegendreTransform1d(int offset, int stride) {
	constexpr auto A = matrixLUDecompose(fourierLegendreTransform<Type, basisOrder, TransformDirection::forward>());
	for (int row = 0; row < basisOrder; row++) {
//		v[row] *= A(row, row);
		printf("\tv[%i] *= %.17e;\n", offset + stride * row, A(row, row));
		for (int column = row + 1; column < basisOrder; column++) {
//			v[row] += A(row, column) * v[column];
			printf("\tv[%i] = std::fma(%.17e, v[%i], v[%i]);\n", offset + stride * row, A(row, column), offset + stride * column, offset + stride * row);
		}
	}
	for (int row = basisOrder - 1; row > 0; row--) {
		for (int column = 0; column < row; column++) {
//			v[row] += A(row, column) * v[column];
			printf("\tv[%i] = std::fma(%.17e, v[%i], v[%i]);\n", offset + stride * row, A(row, column), offset + stride * column, offset + stride * row);
		}
	}
}

template<typename Type, int basisOrder, int dimensionCount>
auto fourierLegendreTransform(int base = 0) {
	constexpr int size = ipow(basisOrder, dimensionCount);
	constexpr int stride = ipow(basisOrder, dimensionCount - 1);
	if constexpr (dimensionCount > 1) {
		for (int offset = 0; offset < size; offset += stride) {
			fourierLegendreTransform<Type, basisOrder, dimensionCount - 1>(base + offset);
		}
	}
	for (int offset = 0; offset < stride; offset++) {
		fourierLegendreTransform1d<Type, basisOrder>(base + offset, stride);
	}
}

template<class T, auto N>
std::string arrayToString(std::array<T, N> const &A) {
	std::ostringstream oss;
	for (auto const &a : A) {
		oss << a << " ";
	}
	return oss.str();
}

int hpx_main(int argc, char *argv[]) {
	enableFPE();
	printf("Reading options...\n");
	processOptions(argc, argv);
	constexpr int basisOrder = 2;
	constexpr int dimensionCount = 3;

//	NodalValues1 A;
//	ModalCoefficients1 B;
//	std::fill(A.begin(), A.end(), 0);
//	std::fill(B.begin(), B.end(), 0);
//	auto f = [](double x, double y, double z) {
//		return x * x * y * y;
//	};
//	std::array<double, 4> points = { -1, -0.447214, 0.447214, 1 };
//	std::array<double, 4> weights = { 1.0 / 6.0, 5.0 / 6.0, 5.0 / 6.0, 1.0 / 6.0 };
//	for (int n = 0; n < A.size(); n++) {
//		std::array<int, 3> I;
//		int k = n;
//		for (int d = 0; d < 3; d++) {
//			I[2 - d] = k % 4;
//			k /= 4;
//		}
//		A[n] = f(points[I[0]], points[I[1]], points[I[2]]);
//	}
//
//	B = legendreAnalyze(A);
//	//A = legendreSynthesize(B);
//	std::cout << arrayToString(A) << "\n" << arrayToString(B) << "\n";
//	using Indices = MultiIndex<basisOrder, dimensionCount>;
//	constexpr int size = Indices::count();
//	constexpr auto strides = Indices::strides();
//	for (int dimensionIndex = 0; dimensionIndex < dimensionCount; ++dimensionIndex) {
//		printf("dimensionIndex = %i\n", dimensionIndex);
//		int const stride = strides[dimensionIndex];
//		int const blockSize = basisOrder * stride;
//		for (int blockIndex = 0; blockIndex < size; blockIndex += blockSize) {
//			printf("\tblockIndex = %i\n", blockIndex);
//			for (int offsetIndex = 0; offsetIndex < stride; ++offsetIndex) {
//				printf("\t\toffsetIndex = %i\n", offsetIndex);
//				for (int basisIndex = 0; basisIndex < basisOrder; ++basisIndex) {
//					printf("\t\t\tbasisIndex = %i\n", basisIndex);
//					printf("\t\t\t\toutputIndex = %i\n", blockIndex + offsetIndex + stride * basisIndex);
//					for (int quadratureIndex = 0; quadratureIndex < basisOrder; ++quadratureIndex) {
//						//input[blockIndex + offsetIndex + stride * quadratureIndex];
//					}
//					//output[blockIndex + offsetIndex + stride * basisIndex];
//				}
//			}
//		}
//	}
//	using RK = typename RungeKutta<T, P>::type;
//	using S = EulerState<T, D>;
//	HyperGrid<S, N, P, RK> grid;
//	grid.initialize(initSodShockTube<T, D>);
//	grid.enforceBoundaryConditions();
//	T t = T(0);
//	T tmax = T(.15);
//	T dt;
//	RK const rk;
//	int iter = 0;
//	while (t < tmax) {
//		grid.output("X", iter, t);
//		std::cout << "i = " << std::to_string(iter);
//		std::cout << "  t = " << std::to_string(t);
//		dt = grid.beginStep();
//		std::cout << "  dt = " << dt << std::endl;
//		for (int s = 0; s < rk.stageCount(); s++) {
//			grid.subStep(dt, s);
//			grid.enforceBoundaryConditions();
//		}
//		grid.endStep();
//		grid.enforceBoundaryConditions();
//		iter++;
//		t += dt;
//	}
	return hpx::local::finalize();
}

int main(int argc, char *argv[]) {
#ifndef NDEBUG
	installFpeHandler();
#endif
	std::vector<std::string> cfg = { "hpx.commandline.allow_unknown=1" };
	cfg.push_back("hpx.stacks.small_size=524288");
	hpx::init_params init_params;
	init_params.cfg = std::move(cfg);
	auto rc = hpx::init(argc, argv, init_params);
	return rc;
}
