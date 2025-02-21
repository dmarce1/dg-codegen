/*
 * test.cpp
 *
 *  Created on: Feb 19, 2025
 *      Author: dmarce1
 */

#include "Real.hpp"
#include "Relativity.hpp"
#include "Utilities.hpp"

#include <silo.h>
#include <valarray>

template<typename T>
struct GRGrid {
	static constexpr int boundWidth = 3;
	static constexpr T zero = T(0), one = T(1);
	static constexpr int XDIM = 0;
	static constexpr int YDIM = 1;
	static constexpr int ZDIM = 2;
	using grid_type = std::valarray<SpacetimeState<T>>;
	using state_type = SpacetimeState<T>;
	GRGrid(Vector<int, NDIM> const &N) :
			intSizes(NDIM), extSizes(NDIM), gridStrides(NDIM), gridSpacing(NDIM) {
		time = T(0);
		frameCount = 0;
		intSize = 1;
		extSize = 1;
		for (int k = 0; k < NDIM; k++) {
			intSizes[k] = N[k];
			extSizes[k] = intSizes[k] + 2 * boundWidth;
			extSize *= extSizes[k];
			intSize *= intSizes[k];
		}
		for (int k = 0; k < NDIM; k++) {
			gridSpacing[k] = one / T(intSizes[k]);
		}
		gridStrides[XDIM] = extSizes[YDIM] * extSizes[ZDIM];
		gridStrides[YDIM] = extSizes[ZDIM];
		gridStrides[ZDIM] = 1;
		minSpacing = gridSpacing.min();
		U = grid_type(extSize);
		for (int dim = 0; dim < NDIM; dim++) {
			std::valarray<size_t> sizes(NDIM);
			int const bwStride = boundWidth * gridStrides[dim];
			sizes[dim] = boundWidth;
			for (int k = 0; k < dim; k++) {
				sizes[k] = intSizes[k];
			}
			for (int k = dim + 1; k < NDIM; k++) {
				sizes[k] = extSizes[k];
			}
			int srcStart1 = zero;
			for (int k = 0; k <= dim; k++) {
				srcStart1 += boundWidth * gridStrides[k];
			}
			size_t const dstStart1 = srcStart1 + intSizes[dim] * gridStrides[dim];
			size_t const dstStart2 = srcStart1 - bwStride;
			size_t const srcStart2 = dstStart1 - bwStride;
			size_t const i0 = 2 * dim;
			size_t const i1 = i0 + 1;
			srcBoundarySlices[i0] = std::gslice(srcStart1, sizes, gridStrides);
			dstBoundarySlices[i0] = std::gslice(dstStart1, sizes, gridStrides);
			srcBoundarySlices[i1] = std::gslice(srcStart2, sizes, gridStrides);
			dstBoundarySlices[i1] = std::gslice(dstStart2, sizes, gridStrides);
		}
	}
	grid_type dUdt() const {
		grid_type dudt(extSize);
		for (int n = 0; n < int(extSize); n++) {
			dudt[n] = U[n].source();
		}
		for (int dim = 0; dim < NDIM; dim++) {
			int const dN = gridStrides[dim];
			T const dxInv = one / gridSpacing[dim];
			for (int np = dN; np < int(extSize); np++) {
				int const nm = np - dN;
				auto const F = riemannFlux(U[nm], U[np], dim);
				dudt[np] += F * dxInv;
				dudt[nm] -= F * dxInv;
			}
		}
		return dudt;
	}
	void initialize() {
		for (int n = 0; n < int(extSize); n++) {
			U[n] = zero;
			U[n].alpha = one;
			for (int k = 0; k < NDIM; k++) {
				U[n].gamma_ll[k, k] = one;
			}
		}
	}
	void output() const {
		static constexpr char meshName[] = "Cartesian";
		static constexpr int silo_data_type = DB_DOUBLE;
		std::vector<T> var(intSize);
		DBShowErrors(DB_ALL, siloErrorHandler);
		DBoptlist *optList = DBMakeOptlist(1);
		double t = time;
		DBAddOption(optList, DBOPT_DTIME, &t);
		std::string const filename = "X." + std::to_string(frameCount) + ".silo";
		DBfile *db = DBCreate(filename.c_str(), DB_CLOBBER, DB_LOCAL, "Astro-Tiger", DB_HDF5);
		char const *const coordnames[NDIM] = { "x", "y", "z" };
		Vector<std::vector<double>, NDIM> xCoordinates;
		for (int dim = 0; dim < NDIM; dim++) {
			xCoordinates[dim].resize(intSizes[dim] + 1);
			for (int n = 0; n <= int(intSizes[dim]); n++) {
				xCoordinates[dim][n] = double(n) / intSizes[dim];
			}
		}
		void const *const coords[NDIM] = { xCoordinates[XDIM].data(), xCoordinates[YDIM].data(), xCoordinates[ZDIM].data() };
		int dims1[NDIM] = { (int) intSizes[XDIM] + 1, (int) intSizes[YDIM] + 1, (int) intSizes[ZDIM] + 1 };
		int const dims2[NDIM] = { (int) intSizes[XDIM], (int) intSizes[YDIM], (int) intSizes[ZDIM] };
		DBPutQuadmesh(db, meshName, coordnames, coords, dims1, NDIM, silo_data_type, DB_COLLINEAR, optList);
		auto const gSlice = std::gslice(0, intSizes, gridStrides);
		auto const V = U[gSlice];
		for (int n = 0; n < state_type::NF; n++) {
			for (int k = 0; k < int(intSize); k++) {
				var[k] = V[k][n];
			}
			DBPutQuadvar1(db, state_type::fieldNames[n], meshName, var.data(), dims2, NDIM, NULL, 0, silo_data_type, DB_ZONECENT, optList);
		}
		DBClose(db);
		DBFreeOptlist(optList);
		frameCount++;
	}
	void step() {
		T dt = T(0.25) * minSpacing;
		U0 = U;
		auto const dudt = dUdt();
		for (int n = 0; n < int(extSize); n++) {
			U[n] += dudt[n] * dt;
		}
		for (int n = 0; n < 2 * NDIM; n++) {
			U[dstBoundarySlices[n]] = U[srcBoundarySlices[n]];
		}
		time += dt;
	}
	void evolveToTime(T tEnd, int totalFrameCount) {
		initialize();
		while (time < tEnd) {
			printf("t = %e", time);
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
	std::valarray<size_t> intSizes;
	std::valarray<size_t> extSizes;
	std::valarray<size_t> gridStrides;
	std::valarray<T> gridSpacing;
	Vector<std::gslice, 2 * NDIM> srcBoundarySlices;
	Vector<std::gslice, 2 * NDIM> dstBoundarySlices;
	grid_type U0;
	grid_type U;
	size_t extSize;
	size_t intSize;
	T time;
	T minSpacing;
	int mutable frameCount;
};

void test() {
	using T = Real;
	GRGrid<T> test(Vector<int, NDIM> { 32, 32, 32 });
	test.evolveToTime(Real(1.0), 10);
}
