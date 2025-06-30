#include <hpx/hpx_init.hpp>
#include "ValArray.hpp"
#include "dgTransforms.hpp"
#include "Options.hpp"
#include "MultiIndex.hpp"
#include "HyperGrid.hpp"
#include "EulerState.hpp"
#include "Real.hpp"
#include "RungeKutta.hpp"

void testRadiation();

int hpx_main(int argc, char *argv[]) {
	printf("\nStarting\n");
	enableFPE();
	processOptions(argc, argv);
	printf("\nPrologue complete\n");
//	constexpr int D = 4;
//	constexpr int N = 3;
//	SquareMatrix<ValArray<double, N>, D> A;
//	for (int n = 0; n < D; n++) {
//		for (int m = 0; m < D; m++) {
//			A(n, m).randomize(-1.0, 1.0);
//		}
//	}
//	for (int n = 0; n < D; n++) {
//		for (int m = 0; m < D; m++) {
//			for (int l = 0; l < N; l++) {
//				printf("%e ", A(n, m)[l]);
//			}
//			printf("\n\n");
//		}
//	}
//	auto const B = matrixInverse(A);
//	auto const C = matrixInverse(B);
//	auto const E = A - C;
//	printf("\n\n");
//	for (int n = 0; n < D; n++) {
//		for (int m = 0; m < D; m++) {
//			for (int l = 0; l < N; l++) {
//				printf("%e ", E(n, m)[l]);
//			}
//			printf("\n\n");
//		}
//	}
//	return hpx::local::finalize();
	constexpr int P = 3;
	constexpr int D = 3;
	constexpr int N = 32;
	using T = Real;
	using RK = typename RungeKutta<T, P>::type;
	HyperGrid<T, D, N, P, RK, EulerState> grid;
	grid.initialize(initSodShockTube<T, D>);
	grid.enforceBoundaryConditions();
	grid.output("X", 0, Real(0.0));
	grid.applyLimiter();
	T t = T(0);
	T tmax = T(1);
	T dt;
	RK const rk;
	int iter = 0;
	while (t < tmax) {
		grid.output("X", iter, t);
		std::cout << "i = " << std::to_string(iter);
		std::cout << "  t = " << std::to_string(t);
		dt = grid.beginStep();
		std::cout << "  dt = " << dt << std::endl;
		for (int s = 0; s < rk.stageCount(); s++) {
			grid.subStep(dt, s);
		}
		grid.endStep();
		iter++;
		t += dt;
	}
	printf("\nStopping\n");
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
