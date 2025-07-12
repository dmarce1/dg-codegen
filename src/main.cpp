#include <algorithm>
#include <memory>
#include <numeric>
#include <type_traits>
#include <vector>
#include <hpx/hpx_init.hpp>
#include "Valarray.hpp"
#include "dgTransforms.hpp"
#include "Options.hpp"
#include "MultiIndex.hpp"
#include "HyperGrid.hpp"
#include "EulerState.hpp"
#include "Real.hpp"
#include "RungeKutta.hpp"

void testRadiation();

int hpx_main(int argc, char *argv[]) {
	Valarray<double> v1;
	Valarray<double> v2;
	printf("\nStarting\n");
	enableFPE();
	processOptions(argc, argv);
	printf("\nPrologue complete\n");
	constexpr int P = 3;
	constexpr int D = 2;
	constexpr int N = 128;
	using T = Real;
//	Valarray<T> a(1.0, N);
//	Valarray<T> b(1.0, N);
//	Valarray<T> d(1.0, N);
//	Valarray<T> e(1.0, N);
//	Valarray<T> c(1.0, N);
//	SquareMatrix<double, 4> A( { { 0.000000, -1.290994, 0.000000, 1.290994 }, { 1.000000, 0.000000, 0.000000, 0.000000 }, { 0.000000, 1.000000, 1.000000,
//			1.000000 }, { 0.000000, 2.500000, 0.000000, 2.500000 } });
//	auto B = matrixInverse(A);
//	auto tmp1 =  (a + b);
//	auto tmp2 = 2.0 * tmp1;
//	c = tmp2;
	using RK = RungeKutta<T, P>::type;
	RK const rk;
	HyperGrid<T, D, N, P, RK, EulerStateHLL> grid;
	grid.initialize(initSodShockTube<T, D>);
	grid.applyLimiter();
	grid.output("X", 0, Real(0.0));
	T t = T(0);
	T tmax = T(.125);
	T dt;
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
	cfg.push_back("hpx.stacks.small_size=1048576");
	hpx::init_params init_params;
	init_params.cfg = std::move(cfg);
	auto rc = hpx::init(argc, argv, init_params);
	return rc;
}
