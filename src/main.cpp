#include <algorithm>
#include <memory>
#include <numeric>
#include <type_traits>
#include <vector>
#include <hpx/hpx_init.hpp>
#include "Octogrid.hpp"
#include "HyperSubgrid.hpp"
#include "dgTransforms.hpp"
#include "Options.hpp"
#include "MultiIndex.hpp"
#include "EulerState.hpp"
#include "Real.hpp"
#include "RungeKutta.hpp"

void testRadiation();
constexpr int P = 3;
constexpr int D = 2;
constexpr int N = 64;
using T = double;
using RK = RungeKutta<T, P>::type;

using SubgridType = HyperSubgrid<T, D, N, P, RK, EulerStateHLLC>;
using OctogridServerType = OctogridServer<SubgridType>;
using OctogridType = hpx::components::component<OctogridServerType>;

HPX_REGISTER_COMPONENT_MODULE();
HPX_REGISTER_COMPONENT (OctogridType);

using GetFaceChildrenAction = typename OctogridServerType::refineAction;
using RefineAction = typename OctogridServerType::getFaceChildrenAction;
using SetAction = typename OctogridServerType::setAction;

HPX_REGISTER_ACTION(GetFaceChildrenAction, getFaceChildrenAction);
HPX_REGISTER_ACTION(RefineAction, refineAction);
HPX_REGISTER_ACTION(SetAction, setAction);

int hpx_main(int argc, char *argv[]) {
	printf("\nStarting\n");
	enableFPE();
	processOptions(argc, argv);
	printf("\nPrologue complete\n");
	RK const rk;
	SubgridType grid;
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
