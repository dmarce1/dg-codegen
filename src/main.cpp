/******************************************************************************
 Copyright (C) 2024  Dominic C. Marcello
 *******************************************************************************/
#include <hpx/hpx_init.hpp>
#include "Hdf5.hpp"
#include "EulersState.hpp"
#include "HyperGrid.hpp"
#include "Options.hpp"

int hpx_main(int argc, char *argv[]) {
	printf("Reading options...\n");
	processOptions(argc, argv);
	constexpr int D = 2;
	constexpr int N = 3;
	constexpr int Nc = 1 << D;
	using vector_t = Basis<double, D, N>::vector_type<D>;
	using function_t = std::function<double(vector_t const&)>;
	function_t fC = [](vector_t const &X) {
		return X[0] * X[0] * X[0] * X[0];
	};
	std::array<function_t, Nc> fF;
	for (Octant<D> o = Octant<D>::begin(); o != Octant<D>::end(); o++) {
		fF[o] = [fC, o](vector_t const &X) {
			auto x = X;
			for (int d = 0; d < D; d++) {
				x[d] = 2.0 * x[d] - (2 * int(o[d]) - 1);
			}
			return fC(x);
		};
	}
	HyperGrid<double, D, 16, N, EulerState<double, D>> grid;
	Math::Vector<double, D> origin = 0;
	grid.initialize(origin, initSodShockTube);
	grid.output("X");
	std::cout << std::endl;
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
