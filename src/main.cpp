#include <hpx/hpx_init.hpp>
#include "Options.hpp"
#include "MultiIndex.hpp"
#include "HyperGrid.hpp"
#include "EulerState.hpp"
#include "RungeKutta.hpp"

int hpx_main(int argc, char *argv[]) {
	enableFPE();
	printf("Reading options...\n");
	processOptions(argc, argv);
	using T = double;
	constexpr int P = 3;
	constexpr int D = 2;
	constexpr int N = 16;
	using RK = typename RungeKutta<T, P>::type;
	using S = EulerState<double, D>;
	HyperGrid<S, N, P, RK> grid;
	grid.initialize(initSodShockTube<T, D>);
	grid.enforceBoundaryConditions();
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
			grid.enforceBoundaryConditions();
			grid.subStep(dt, s);
		}
		grid.enforceBoundaryConditions();
		grid.endStep();
		iter++;
	}
//	constexpr int D = 2;
//	constexpr int P = 1;
//	constexpr int Nc = 1 << D;
//	using vector_t = Basis<double, D, P>::vector_type<D>;
//	using function_t = std::function<double(vector_t const&)>;
//	function_t fC = [](vector_t const &X) {
//		return X[0] * X[0] * X[0] * X[0];
//	};
//	std::array<function_t, Nc> fF;
//	for (Octant<D> o = Octant<D>::begin(); o != Octant<D>::end(); o++) {
//		fF[o] = [fC, o](vector_t const &X) {
//			auto x = X;
//			for (int d = 0; d < D; d++) {
//				x[d] = 2.0 * x[d] - (2 * int(o[d]) - 1);
//			}
//			return fC(x);
//		};
//	}
//	using rk_type = typename RungeKutta<double, P>::type;
//	HyperGrid<double, D, 16, P, EulerState<double, D>, rk_type> grid;
//	Math::Vector<double, D> origin = 0;
//	grid.initialize(origin, initSodShockTube);
//	grid.output("X");
//	std::cout << "dt = " << std::to_string(grid.step()) << std::endl;
//	grid.output("X");
//	std::cout << std::endl;
	return hpx::local::finalize();
}

int main(int argc, char *argv[]) {
	std::vector<std::string> cfg = { "hpx.commandline.allow_unknown=1" };
	cfg.push_back("hpx.stacks.small_size=524288");
	hpx::init_params init_params;
	init_params.cfg = std::move(cfg);
	auto rc = hpx::init(argc, argv, init_params);
	return rc;
}
