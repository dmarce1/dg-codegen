/*
 * HydrodynamicsGrid.hpp
 *
 *  Created on: Jan 15, 2025
 *      Author: dmarce1
 */

#ifndef INCLUDE_HYDROGRID_HPP_
#define INCLUDE_HYDROGRID_HPP_

#include "Hydrodynamics.hpp"
#include "Limiters.hpp"
#include "TriangularArray.hpp"
#include "Options.hpp"

#include <valarray>

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

template<typename T, int D, int P>
struct HydroGrid {
	using triindices_t = TriangularIndices<D, P>;
	using array_t = std::vector<std::vector<T>>;
	using inner_array_t = std::vector<T>;
	using size_array_t = std::vector<size_t>;
	using state_t = ConservedState<T, D>;
	using array_state_t = ConservedState<array_t, D>;
	static constexpr T zero = T(0), half = T(0.5), one = T(1), two = T(2);
	static constexpr int NF = state_t::NFields;
	static constexpr int P3 = triindices_t::size();
	static constexpr auto slopeLimiter = vanLeer<T>;
	HydroGrid() :
			gridOptions(), N(gridOptions.gridLength), BW(1), dx(one / T(N)), gridDims(N + 2 * BW), gridDomain(Interval<T, D>(one).expand(T(2 * BW) * dx)) {
		N3 = vectorProduct(gridDims);
		gridStrides = size_array_t(D);
		gridStrides[0] = 1;
		for (int k = 0; k < D - 1; k++) {
			gridStrides[k + 1] = gridDims[k] * gridStrides[k];
		}
		for (int f = 0; f < NF; f++) {
			U[f] = std::vector<std::vector<T>>(P3, std::vector<T>(N3));
		}
	}
	void reconstruct() {
		for (int f = 0; f < NF; f++) {
			triindices_t pp;
			for (++pp; pp != pp.end(); pp++) {
				bool firstPass = true;
				auto &Vp = U[f][pp];
				for (int k = 0; k < D; k++) {
					size_t const dN = gridStrides[k];
					T const norm = one / Real(2 * pp.indexAt(k) - 1);
					auto indices = pp.getIndices;
					if (indices[k] == 0) {
						continue;
					} else {
						--indices[k];
					}
					triindices_t const p0(indices);
					auto const &V0 = U[f][p0];
					for (int n = dN; n < N3 - dN; n++) {
						T slope;
						T const v = V0[n];
						slope = norm * slopeLimiter(V0[n + dN] - v, v - V0[n - dN]);
						if (!firstPass) {
							slope = minmod(slope, Vp[n]);
						}
						Vp[n] = slope;
					}
					firstPass = false;
				}
			}
		}
	}
private:
	std::gslice createGslice(Interval<int, D> const &thisBox) {
		size_t start = 0;
		size_array_t sizes(D);
		for (int k = 0; k < D; k++) {
			sizes[k] = thisBox.span(k);
			start += thisBox.begin(k) * gridStrides[k];
		}
		return std::gslice(start, sizes, gridStrides);
	}
	GridOptions<T> gridOptions;
	int N;
	int N3;
	int BW;
	T dx;
	Vector<int, D> gridDims;
	Interval<T, D> gridDomain;
	size_array_t gridStrides;
	array_state_t U;
};

#endif /* INCLUDE_HYDROGRID_HPP_ */
