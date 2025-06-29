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
	enableFPE();
//	testRadiation();
//	return hpx::local::finalize();
	printf("Reading options...\n");
	processOptions(argc, argv);
	constexpr int P = 3;
	constexpr int D = 3;
	constexpr int N = 128;
	using T = double;
	using RK = typename RungeKutta<T, P>::type;
	HyperGrid<T, D, N, P, RK, EulerState> grid;
	grid.initialize(initSodShockTube<T, D>);
	grid.enforceBoundaryConditions();
	grid.applyLimiter();
	T t = T(0);
	T tmax = T(.15);
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
			grid.enforceBoundaryConditions();
		}
		grid.endStep();
		grid.enforceBoundaryConditions();
		iter++;
		t += dt;
	}
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
