/*
 * GRGrid.hpp
 *
 *  Created on: Feb 21, 2025
 *      Author: dmarce1
 */

#ifndef INCLUDE_GRGRID_HPP_
#define INCLUDE_GRGRID_HPP_

#include "GridAttributes.hpp"

#include <valarray>

#include <silo.h>

#include "Spacetime.hpp"

template<typename T, typename R>
struct Grid {
	static constexpr R zero = R(0), one = R(1);
	using grid_type = std::valarray<T>;
	Grid(Math::Vector<int, NDIM> const &N, int boundWidth) :
			gAttr(N, boundWidth) {
		U = grid_type(gAttr.extSize);
		for (int dim = 0; dim < NDIM; dim++) {
			std::valarray<size_t> sizes(NDIM);
			int const bwStride = boundWidth * gAttr.gridStrides[dim];
			sizes[dim] = boundWidth;
			for (int k = 0; k < dim; k++) {
				sizes[k] = gAttr.intSizes[k];
			}
			for (int k = dim + 1; k < NDIM; k++) {
				sizes[k] = gAttr.extSizes[k];
			}
			int srcStart1 = zero;
			for (int k = 0; k <= dim; k++) {
				srcStart1 += boundWidth * gAttr.gridStrides[k];
			}
			size_t const dstStart1 = srcStart1 + gAttr.intSizes[dim] * gAttr.gridStrides[dim];
			size_t const dstStart2 = srcStart1 - bwStride;
			size_t const srcStart2 = dstStart1 - bwStride;
			size_t const i0 = 2 * dim;
			size_t const i1 = i0 + 1;
			srcBoundarySlices[i0] = std::gslice(srcStart1, sizes, gAttr.gridStrides);
			dstBoundarySlices[i0] = std::gslice(dstStart1, sizes, gAttr.gridStrides);
			srcBoundarySlices[i1] = std::gslice(srcStart2, sizes, gAttr.gridStrides);
			dstBoundarySlices[i1] = std::gslice(dstStart2, sizes, gAttr.gridStrides);
		}
	}
	T interpolate(Math::Vector<R, NDIM> const &X) const {
		auto const TSC = [](R x) {
			x = abs(x);
			if (x < R(0.5)) {
				return R(0.75) - x * x;
			} else if (x < R(1.5)) {
				return R(0.5) * nSquared(R(1.5) - x);
			} else {
				return R(0.0);
			}
		};
		int const W = 3;
		Math::Vector<int, NDIM> I;
		T sum = zero;
		for (int k = 0; k < NDIM; k++) {
			I[k] = X[k] * gAttr.intSizes[k] + gAttr.boundWidth;
			X[k] -= R(I[k]) + R(0.5);
		}
		int const di = gAttr.gridStrides[0];
		int const dj = gAttr.gridStrides[1];
		int const dk = gAttr.gridStrides[2];
		int const ib = I[0] - 1;
		int const ie = I[0] + 1;
		int const jb = I[1] - 1;
		int const je = I[1] + 1;
		int const kb = I[2] - 1;
		int const ke = I[2] + 1;
		for (int i = ib; i <= ie; i++) {
			for (int j = jb; j <= je; j++) {
				for (int k = kb; k <= ke; k++) {
					R const xTSC = TSC(X[0] - R(i));
					R const yTSC = TSC(X[1] - R(j));
					R const zTSC = TSC(X[2] - R(k));
					sum += U[i * di + j * dj + k * dk] * xTSC * yTSC * zTSC;
				}
			}
		}
		return sum;
	}
	void enforceBoundaryConditions() {
		for (int n = 0; n < 2 * NDIM; n++) {
			U[dstBoundarySlices[n]] = U[srcBoundarySlices[n]];
		}
	}
protected:
	GridAttributes<R> gAttr;
	Math::Vector<std::gslice, 2 * NDIM> srcBoundarySlices;
	Math::Vector<std::gslice, 2 * NDIM> dstBoundarySlices;
	grid_type U;
};

template<typename T>
struct GRGrid: public Grid<SpacetimeState<T>, T> {
	static constexpr int boundWidth = 1;
	static constexpr T zero = T(0), one = T(1);
	using base_type = Grid<SpacetimeState<T>, T>;
	using grid_type = base_type::grid_type;
	using state_type = SpacetimeState<T>;
	GRGrid(Math::Vector<int, NDIM> const &N) :
			base_type(N, boundWidth), gAttr(base_type::gAttr), U(base_type::U) {
		time = T(0);
		frameCount = 0;
	}
	grid_type dUdt() const {
		grid_type dudt(gAttr.extSize);
		for (int n = 0; n < int(gAttr.extSize); n++) {
			dudt[n] = U[n].source();
		}
		for (int dim = 0; dim < NDIM; dim++) {
			int const dN = gAttr.gridStrides[dim];
			T const dxInv = one / gAttr.gridSpacing[dim];
			for (int np = dN; np < int(gAttr.extSize); np++) {
				int const nm = np - dN;
				SpacetimeState<T> const F = riemannFlux(U[nm], U[np], dim);
				dudt[np] += F * dxInv;
				dudt[nm] -= F * dxInv;
			}
		}
		return dudt;
	}
	void initialize() {
		for (int n = 0; n < int(gAttr.extSize); n++) {
			U[n] = zero;
			U[n].alpha = one;
			for (int k = 0; k < NDIM; k++) {
				U[n].gamma_ll[k, k] = one;
			}
		}
	}
	void output() const {
		using namespace Math;
		static constexpr char meshName[] = "Cartesian";
		static constexpr int silo_data_type = DB_DOUBLE;
		std::vector<T> var(gAttr.intSize);
		DBShowErrors(DB_ALL, siloErrorHandler);
		DBoptlist *optList = DBMakeOptlist(1);
		double t = time;
		DBAddOption(optList, DBOPT_DTIME, &t);
		std::string const filename = "X." + std::to_string(frameCount) + ".silo";
		DBfile *db = DBCreate(filename.c_str(), DB_CLOBBER, DB_LOCAL, "Astro-Tiger", DB_HDF5);
		char const *const coordnames[NDIM] = { "x", "y", "z" };
		Vector<std::vector<double>, NDIM> xCoordinates;
		for (int dim = 0; dim < NDIM; dim++) {
			xCoordinates[dim].resize(gAttr.intSizes[dim] + 1);
			for (int n = 0; n <= int(gAttr.intSizes[dim]); n++) {
				xCoordinates[dim][n] = double(n) / gAttr.intSizes[dim];
			}
		}
		void const *const coords[NDIM] = { xCoordinates[XDIM].data(), xCoordinates[YDIM].data(), xCoordinates[ZDIM].data() };
		int dims1[NDIM] = { (int) gAttr.intSizes[XDIM] + 1, (int) gAttr.intSizes[YDIM] + 1, (int) gAttr.intSizes[ZDIM] + 1 };
		int const dims2[NDIM] = { (int) gAttr.intSizes[XDIM], (int) gAttr.intSizes[YDIM], (int) gAttr.intSizes[ZDIM] };
		DBPutQuadmesh(db, meshName, coordnames, coords, dims1, NDIM, silo_data_type, DB_COLLINEAR, optList);
		auto const gSlice = std::gslice(0, gAttr.intSizes, gAttr.gridStrides);
		auto const V = U[gSlice];
		for (int n = 0; n < state_type::NF; n++) {
			for (int k = 0; k < int(gAttr.intSize); k++) {
				var[k] = V[k][n];
			}
			DBPutQuadvar1(db, state_type::fieldNames[n], meshName, var.data(), dims2, NDIM, NULL, 0, silo_data_type, DB_ZONECENT, optList);
		}
		DBClose(db);
		DBFreeOptlist(optList);
		frameCount++;
	}
	void step() {
		T dt = T(0.25) * gAttr.minSpacing;
		auto const dudt = dUdt();
		for (int n = 0; n < int(gAttr.extSize); n++) {
			U[n] += dudt[n] * dt;
		}
		base_type::enforceBoundaryConditions();
		time += dt;
	}
	void evolveToTime(T tEnd, int totalFrameCount) {
		initialize();
		while (time < tEnd) {
			printf("t = %e", double(time));
			if (totalFrameCount * time >= frameCount * tEnd) {
				printf(" (output frame = %i)", frameCount);
				output();
			}
			printf("\n");
			step();
		}
		output();
	}
private:
	static void siloErrorHandler(char *errorString_) {
		std::string const errorString(errorString_);
		std::cout << "SILO returned an error." << "\n";
		std::cout << errorString << "\n";
		std::cout << "Aborting..." << "\n";
		abort();
	}
	T time;
	int mutable frameCount;
	GridAttributes<T> const &gAttr;
	grid_type &U;
};

#endif /* INCLUDE_GRGRID_HPP_ */
