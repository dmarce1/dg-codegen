/*
 * HydrodynamicsGrid.hpp
 *
 *  Created on: Jan 15, 2025
 *      Author: dmarce1
 */

#ifndef INCLUDE_HYDROGRID_HPP_
#define INCLUDE_HYDROGRID_HPP_

#include "Hydrodynamics.hpp"
#include "MultiArray.hpp"
#include "TriangularArray.hpp"
#include "Options.hpp"

template<typename, int, int>
struct HydroGrid;

using namespace Hydrodynamics;
using namespace Math;

template<typename Type>
struct GridOptions {
	int gridLength;
	int spatialOrder;
	Type gridScale;
	GridOptions() {
		auto const &opts = getOptions();
		gridLength = opts.gridLength;
		gridScale = Type(opts.gridScale);
	}
};

template<typename Type, int DimCount, int MomCount>
struct HydroGrid {
	using array_type = MultiArray<Type, DimCount>;
	using state_type = ConservedState<Type, DimCount>;
	using state_array_type = ConservedState<array_type, DimCount>;
	static constexpr int NF = state_type::NFields;
	HydroGrid() :
			gridOptions(), P3(TriangularIndices<DimCount, MomCount>::Size), U(P3) {
		int const nCells = gridOptions.gridLength;
		for (int p = 0; p < P3; p++) {
			for (int f = 0; f < NF; f++) {
				U[p][f] = array_type(nCells);
			}
		}
	}
private:
	GridOptions<Type> gridOptions;
	int const P3;
	std::vector<state_array_type> U;
};

#endif /* INCLUDE_HYDROGRID_HPP_ */
