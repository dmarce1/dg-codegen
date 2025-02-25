/*
 * GRGrid.hpp
 *
 *  Created on: Feb 21, 2025
 *      Author: dmarce1
 */

#ifndef INCLUDE_GRGRID_HPP_
#define INCLUDE_GRGRID_HPP_

#include "GridAttributes.hpp"
#include "Vector.hpp"
#include <silo.h>

template<typename T>
SymmetricMatrix<T, NDIM> matrixInverse(SymmetricMatrix<T, NDIM> const &A) {
	static constexpr T one = T(1);
	SymmetricMatrix<T, NDIM> Ainv;
	Ainv[0, 0] = +(A[1, 1] * A[2, 2] - A[1, 2] * A[2, 1]);
	Ainv[0, 1] = -(A[0, 1] * A[2, 2] - A[0, 2] * A[2, 1]);
	Ainv[0, 2] = +(A[0, 1] * A[1, 2] - A[0, 2] * A[1, 1]);
	Ainv[1, 0] = -(A[1, 0] * A[2, 2] - A[1, 2] * A[2, 0]);
	Ainv[1, 1] = +(A[0, 0] * A[2, 2] - A[0, 2] * A[2, 0]);
	Ainv[1, 2] = -(A[0, 0] * A[1, 2] - A[0, 2] * A[1, 0]);
	Ainv[2, 0] = +(A[1, 0] * A[2, 1] - A[1, 1] * A[2, 0]);
	Ainv[2, 1] = -(A[0, 0] * A[2, 1] - A[0, 1] * A[2, 0]);
	Ainv[2, 2] = +(A[0, 0] * A[1, 1] - A[0, 1] * A[1, 0]);
	T const detInv = one / (A[0, 0] * Ainv[0, 0] + A[0, 1] * Ainv[0, 1] + A[0, 2] * Ainv[0, 2]);
	return Ainv * detInv;
}

template<typename RealType>
struct Einstein {
	using real_type = RealType;
	using grid_type = std::valarray<real_type>;
	using spacetime_type = SymmetricMatrix<grid_type, DIM4>;
	static constexpr real_type zero = real_type(0), half = real_type(0.5), one = real_type(1), two = real_type(2);
	static constexpr int BW = 2;
	Einstein(int N_) :
			gAttr(Vector<int, NDIM>( { N_, N_, N_ }), BW), gamma(grid_type(zero, gAttr.extSize)), U(grid_type(zero, gAttr.extSize)) {
	}
	void step(SymmetricMatrix<grid_type, DIM4> const &T, real_type dt) {
		constexpr real_type G = real_type(-16.0 * M_PI);
		auto const &dxinv = gAttr.intSizes;
		auto const dx = one / dxinv;
		auto const &strides = gAttr.extStrides;
		SymmetricMatrix < grid_type, NDIM > gamma_inv(grid_type(zero, gAttr.extSize));
		Vector<SymmetricMatrix<grid_type, NDIM>, NDIM> D(grid_type(zero, gAttr.extSize));
		Vector<SymmetricMatrix<grid_type, NDIM>, NDIM> Gamma(grid_type(zero, gAttr.extSize));
		SymmetricMatrix < grid_type, NDIM > Ricci(grid_type(zero, gAttr.extSize));
		grid_type R(grid_type(zero, gAttr.extSize));
		for (int n = 0; n < gAttr.extSize; n++) {
			SymmetricMatrix<real_type, NDIM> ginv;
			for (int i = 0; i < NDIM; i++) {
				for (int j = 0; j <= i; j++) {
					ginv[i, j] = gamma[i, j][n];
				}
			}
			ginv = matrixInverse(ginv);
			for (int i = 0; i < NDIM; i++) {
				for (int j = 0; j <= i; j++) {
					gamma_inv[i, j][n] = ginv[i, j];
				}
			}
		}
		for (int k = 0; k < NDIM; k++) {
			size_t const dk = strides[k];
			for (int i = 0; i < NDIM; i++) {
				for (int j = 0; j <= i; j++) {
					for (int n = dk; n < gAttr.extSize - dk; n++) {
						D[k][i, j][n] = half * (gamma[i, j][n + dk] - gamma[i, j][n - dk]) * dxinv;
					}
				}
			}
		}
		for (int k = 0; k < NDIM; k++) {
			for (int i = 0; i < NDIM; i++) {
				for (int j = 0; j <= i; j++) {
					for (int n = 0; n < gAttr.extSize; n++) {
						Gamma[k][i, j] = zero;
						for (int l = 0; l < NDIM; l++) {
							Gamma[k][i, j] += half * gamma_inv[k, l][n] * (D[i][j, l][n] + D[j][i, l][n] - D[l][i, j][n]);
						}
					}
				}
			}
		}
		for (int i = 0; i < NDIM; i++) {
			for (int j = 0; j <= i; j++) {
				int const dj = strides[j];
				for (int n = gAttr.start; n <= gAttr.stop; n++) {
					Ricci[i, j][n] = zero;
					for (int k = 0; k < NDIM; k++) {
						int const dk = strides[k];
						Ricci[i, j][n] += half * (Gamma[k][i, j][n + dk] - Gamma[k][i, j][n - dk]) * dxinv;
						Ricci[i, j][n] += half * (Gamma[k][k, i][n + dj] - Gamma[k][k, i][n - dj]) * dxinv;
					}
					for (int k = 0; k < NDIM; k++) {
						for (int l = 0; l < NDIM; l++) {
							Ricci[i, j][n] += Gamma[k][i, j][n] * Gamma[l][k, l][n] - Gamma[l][i, k][n] * Gamma[k][l, j][n];
						}
					}
				}
			}
		}
		for (int n = gAttr.start; n <= gAttr.stop; n++) {
			R[n] = zero;
			for (int i = 0; i < NDIM; i++) {
				R[n] += gamma_inv[i, i][n] * Ricci[i, i][n];
				for (int j = 0; j < i; j++) {
					R[n] += two * gamma_inv[i, j][n] * Ricci[i, j][n];
				}
			}
		}
	}
	void output(DBfile *db, DBoptlist *optList, SymmetricMatrix<grid_type, DIM4> const *Tptr = nullptr) const {
		using namespace Math;
		static constexpr int silo_data_type = DB_DOUBLE;
		char const *const coordnames[NDIM] = { "x", "y", "z" };
		Vector<std::vector<real_type>, NDIM> xCoordinates;
		for (int dim = 0; dim < NDIM; dim++) {
			xCoordinates[dim].resize(gAttr.intSizes[dim] + 1);
			for (int n = 0; n <= int(gAttr.intSizes[dim]); n++) {
				xCoordinates[dim][n] = real_type(n) / real_type(gAttr.intSizes[dim]);
			}
		}
		void const *const coords[NDIM] = { xCoordinates[XDIM].data(), xCoordinates[YDIM].data(), xCoordinates[ZDIM].data() };
		int const dims2[NDIM] = { (int) gAttr.intSizes[XDIM], (int) gAttr.intSizes[YDIM], (int) gAttr.intSizes[ZDIM] };
		int dims1[NDIM] = { dims2[0] + 1, dims2[1] + 1, dims2[2] + 1 };
		DBPutQuadmesh(db, "quadMesh", coordnames, coords, dims1, NDIM, silo_data_type, DB_COLLINEAR, optList);
		size_t const start = BW * (gAttr.extStrides.sum());
		auto const gSlice = std::gslice(start, gAttr.intSizes, gAttr.extStrides);
		for (int i = 0; i < DIM4; i++) {
			for (int j = 0; j <= i; j++) {
				for (int k = 0; k < DIM4; k++) {
					std::valarray<real_type> const u = U[k][i, j][gSlice];
					std::string fieldname = std::string("g") + std::to_string(i) + std::to_string(j) + std::string("_") + std::to_string(k);
					DBPutQuadvar1(db, fieldname.c_str(), "quadMesh", std::begin(u), dims2, NDIM, NULL, 0, silo_data_type, DB_ZONECENT, optList);
				}
			}
		}
		if (Tptr) {
			auto const &T = *Tptr;
			for (int i = 0; i < DIM4; i++) {
				for (int j = 0; j <= i; j++) {
					std::valarray<real_type> const u = T[i, j][gSlice];
					std::string fieldname = std::string("T") + std::to_string(i) + std::to_string(j);
					DBPutQuadvar1(db, fieldname.c_str(), "quadMesh", std::begin(u), dims2, NDIM, NULL, 0, silo_data_type, DB_ZONECENT, optList);
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
	auto metricDerivatives(real_type x0, real_type y0, real_type z0) const {
		using return_type = Vector<SymmetricMatrix<real_type, DIM4>, DIM4>;
		return_type u = return_type(zero);
		int const i0 = x0 * real_type(gAttr.intSizes[XDIM]) + BW;
		int const j0 = y0 * real_type(gAttr.intSizes[YDIM]) + BW;
		int const k0 = z0 * real_type(gAttr.intSizes[ZDIM]) + BW;
		Vector<int, NDIM> originIndices( { i0, j0, k0 });
		Vector<int, NDIM> indices;
		int &i = indices[XDIM];
		int &j = indices[YDIM];
		int &k = indices[ZDIM];
		Vector<SymmetricMatrix<real_type, DIM4>, DIM4> V;
		Vector<real_type, NDIM> const X0( { x0, y0, z0 });
		size_t const start = (i0 - 1) * gAttr.extStrides[XDIM] + (j0 - 1) * gAttr.extStrides[YDIM] + (k0 - 1) * gAttr.extStrides[ZDIM];
		for (i = 0; i < 3; i++) {
			for (j = 0; j < 3; j++) {
				for (k = 0; k < 3; k++) {
					size_t index = start;
					for (int dim = 0; dim < NDIM; dim++) {
						index += indices[dim] * gAttr.extStrides[dim];
					}
					real_type weight = one;
					weight *= kernelTSC(abs(real_type(i0 + i - 1) - real_type(gAttr.intSizes[XDIM]) * X0[XDIM]));
					weight *= kernelTSC(abs(real_type(j0 + j - 1) - real_type(gAttr.intSizes[YDIM]) * X0[YDIM]));
					weight *= kernelTSC(abs(real_type(k0 + k - 1) - real_type(gAttr.intSizes[ZDIM]) * X0[ZDIM]));
					for (int l = 0; l < DIM4; l++) {
						for (int m = 0; m < DIM4; m++) {
							for (int n = 0; n <= m; n++) {
								u[l][m, n] += weight * U[l][m, n][index];
							}
						}
					}
				}
			}
		}
		for (int k = 0; k < DIM4; k++) {
			u[k] = reverseTrace(u[k]);
		}
		return u;
	}
	auto getStateVars(int xi, int yi, int zi) const {
		Vector<SymmetricMatrix<real_type, DIM4>, DIM4> v;
		int const index = xi * gAttr.extSizes[XDIM] + yi * gAttr.extSizes[YDIM] + zi * gAttr.extSizes[ZDIM];
		for (int k = 0; k < DIM4; k++) {
			for (int i = 0; i < DIM4; i++) {
				for (int j = 0; j <= i; j++) {
					v[k][i, j] = U[k][i, j][index];
				}
			}
		}
		return v;
	}
private:
	GridAttributes<real_type> gAttr;
	SymmetricMatrix<grid_type, NDIM> gamma;
	Vector<spacetime_type, DIM4> U;
}
;
#endif /* INCLUDE_GRGRID_HPP_ */
