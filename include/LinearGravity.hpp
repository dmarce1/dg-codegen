/*
 * LinearGravity.hpp
 *
 *  Created on: Feb 22, 2025
 *      Author: dmarce1
 */

#ifndef INCLUDE_LINEARGRAVITY_HPP_
#define INCLUDE_LINEARGRAVITY_HPP_

#include <silo.h>

#include "GridAttributes.hpp"
#include "Real.hpp"
#include "Vector.hpp"
#include "Kernels.hpp"

using namespace Math;

static void siloErrorHandler(char *errorString_) {
	std::string const errorString(errorString_);
	std::cout << "SILO returned an error." << "\n";
	std::cout << errorString << "\n";
	std::cout << "Aborting..." << "\n";
	abort();
}

struct LinearGravity {
	using real_type = Real;
	using grid_type = std::valarray<real_type>;
	static constexpr real_type zero = real_type(0), half = real_type(0.5), one = real_type(1);
	static constexpr int BW = 1;
	static constexpr real_type CFL = real_type(0.25);
	LinearGravity(int N_) :
			gAttr(Vector<int, NDIM>( { N_, N_, N_ }), BW), Ni(gAttr.intSize), Ne(gAttr.extSize), sizes(gAttr.intSizes), strides(gAttr.gridStrides), N3(
					Ne * Ne * Ne), start(Ne * (Ne * BW + BW) + BW), U(grid_type(zero, N3)) {
	}
	void step(SymmetricMatrix<grid_type, DIM4> const &S) {
		constexpr real_type G = real_type(-8.0 * M_PI);
		real_type const dx = one / real_type(Ni);
		real_type const dt = CFL * dx;
		Vector<grid_type, DIM4> D;
		for (int i = 0; i < DIM4; i++) {
			for (int j = 0; j <= i; j++) {
				D[0] = (dt * G) * S[i, j];
				for (int k = 1; k < DIM4; k++) {
					int const dk = strides[k - 1];
					grid_type const &Uk = U[k][i, j];
					grid_type const &U0 = U[0][i, j];
					grid_type const Fk = half * ((U0 + U0.cshift(dk)) - (Uk - Uk.cshift(dk)));
					grid_type const F0 = half * ((Uk + Uk.cshift(dk)) - (U0 - U0.cshift(dk)));
					D[0] += (F0.cshift(-dk) - F0) * (one / dx);
					D[k] += (Fk.cshift(-dk) - Fk) * (one / dx);
				}
				for (int k = 0; k < DIM4; k++) {
					U[k][i, j] += dt * D[k];
				}
			}
		}
	}
	void output(DBfile *db, DBoptlist *optList) const {
		using namespace Math;
		static constexpr int XDIM = 0;
		static constexpr int YDIM = 1;
		static constexpr int ZDIM = 2;
		static constexpr char meshName[] = "Cartesian";
		static constexpr int silo_data_type = DB_DOUBLE;
		char const *const coordnames[NDIM] = { "x", "y", "z" };
		Vector<std::vector<double>, NDIM> xCoordinates;
		for (int dim = 0; dim < NDIM; dim++) {
			xCoordinates[dim].resize(Ni + 1);
			for (int n = 0; n <= int(Ni); n++) {
				xCoordinates[dim][n] = double(n) / double(Ni);
			}
		}
		void const *const coords[NDIM] = { xCoordinates[XDIM].data(), xCoordinates[YDIM].data(), xCoordinates[ZDIM].data() };
		int const dims2[NDIM] = { (int) Ni, (int) Ni, (int) Ni };
		int dims1[NDIM] = { dims2[0] + 1, dims2[1] + 1, dims2[2] + 1 };
		DBPutQuadmesh(db, meshName, coordnames, coords, dims1, NDIM, silo_data_type, DB_COLLINEAR, optList);
		auto const gSlice = std::gslice(start, sizes, strides);
		for (int i = 0; i < DIM4; i++) {
			for (int j = 0; j <= i; j++) {
				for (int k = 0; k < DIM4; k++) {
					std::valarray<real_type> const u = U[k][i, j][gSlice];
					std::string fieldname = std::string("g") + std::to_string(i) + std::to_string(j) + std::string("_") + std::to_string(k);
					DBPutQuadvar1(db, fieldname.c_str(), meshName, std::begin(u), dims2, NDIM, NULL, 0, silo_data_type, DB_ZONECENT, optList);
				}
			}
		}
	}
	void enforceBoundaryConditions() {
		for (int n = 0; n < 2 * NDIM; n++) {
			for (int i = 0; i < DIM4; i++) {
				for (int j = 0; j <= i; j++) {
					for (int k = 0; k < DIM4; k++) {
						U[k][i, j][gAttr.dstBoundarySlices[n]] = U[k][i, j][gAttr.srcBoundarySlices[n]];
					}
				}
			}
		}
	}
	auto metricFunction() const {
		auto const f = [this](real_type x0, real_type y0, real_type z0) {
			using return_type = Vector<SymmetricMatrix<real_type, DIM4>, DIM4>;
			return_type u = return_type(zero);
			int const i0 = x0 * real_type(Ni) + BW;
			int const j0 = y0 * real_type(Ni) + BW;
			int const k0 = z0 * real_type(Ni) + BW;
			for (int i = i0 - 1; i <= i0 + 1; i++) {
				for (int j = j0 - 1; j <= j0 + 1; j++) {
					for (int k = k0 - 1; k <= k0 + 1; k++) {
						int index = strides[0] * periodicModulus(i, Ni);
						index += strides[1] * periodicModulus(j, Ni);
						index += strides[2] * periodicModulus(k, Ni);
						real_type const x = real_type(sizes[0]) * x - real_type(i) - half;
						real_type const y = real_type(sizes[1]) * y - real_type(j) - half;
						real_type const z = real_type(sizes[2]) * z - real_type(k) - half;
						for (int l = 0; l < DIM4; l++) {
							for (int m = 0; m < DIM4; m++) {
								for (int n = 0; n <= m; n++) {
									u[l][m, n] += kernelTSC(x) * kernelTSC(y) * kernelTSC(z) * U[l][m, n][index];
								}
							}
						}
					}
				}
			}
			return u;
		};
		return f;
	}
private:
	GridAttributes<real_type> gAttr;
	size_t const &Ni;
	size_t const &Ne;
	std::valarray<size_t> const &sizes;
	std::valarray<size_t> const &strides;
	size_t N3;
	size_t start;
	Vector<SymmetricMatrix<grid_type, DIM4>, DIM4> U;
};

#endif /* INCLUDE_LINEARGRAVITY_HPP_ */
