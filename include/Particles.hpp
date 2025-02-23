/*
 * Particles.hpp
 *
 *  Created on: Feb 21, 2025
 *      Author: dmarce1
 */

#ifndef INCLUDE_PARTICLES_HPP_
#define INCLUDE_PARTICLES_HPP_

#include "GridAttributes.hpp"
#include "Vector.hpp"
#include "Interval.hpp"
#include "Kernels.hpp"

#include <numeric>
#include <valarray>

template<typename T>
struct Particle {
	Math::Vector<T, NDIM> u;
	Math::Vector<T, NDIM> x;
};

template<typename real_type>
Vector<real_type, DIM4> velocityTo4(Vector<real_type, NDIM> const &u3) {
	static constexpr real_type one = real_type(1);
	Vector<real_type, DIM4> u4;
	real_type const W = sqrt(one + vectorDotProduct(u3, u3));
	for (int n = 0; n < NDIM; n++) {
		u4[n + 1] = u3[n];
	}
	u4[0] = W;
	return u4;
}

double rand1() {
	return (rand() + 0.5) / double(RAND_MAX);
}

template<typename RealType>
struct Particles {
	using real_type = RealType;
	using grid_type = std::valarray<real_type>;
	static constexpr real_type zero = real_type(0), two = real_type(2), half = real_type(0.5), one = real_type(1);
	Particles(size_t count, int Ngrid, int BW) :
			gAttr(Vector<int, NDIM> { Ngrid, Ngrid, Ngrid }, BW), X(grid_type(zero, count)), U(grid_type(zero, count)), particleGrid(gAttr.intSize), M(
					getOptions().totalMass / real_type(count)), time(0.0), N(count) {
		static constexpr real_type r0(0.1);
		for (int n = 0; n < int(count); n++) {
			bool outside;
			do {
				for (int dim = 0; dim < NDIM; dim++) {
					real_type r(rand1());
					r = (two * r - one) * real_type(r0) + half;
					X[dim][n] = r;
				}
				real_type const radius = sqrt(nSquared(X[XDIM][n] - half) + nSquared(X[YDIM][n] - half) + nSquared(X[ZDIM][n] - half));
				outside = radius > real_type(r0);
			} while (outside);
		}
	}
	void swap_particles(size_t i, size_t j) {
		for (int dim = 0; dim < NDIM; dim++) {
			std::swap(X[dim][i], X[dim][j]);
			std::swap(U[dim][i], U[dim][j]);
		}
	}
	void sort_particles(std::pair<size_t, size_t> particleRange, Interval<size_t, NDIM> indexBox) {
		if (indexBox.volume() <= 1) {
			int index = 0;
			for (int dim = 0; dim < NDIM; dim++) {
				index += indexBox.begin(dim) * gAttr.intSizes[dim];
			}
			particleGrid[index] = particleRange;
		} else {
			int const sortDim = indexBox.longest();
			size_t const splitIndex = (indexBox.begin(sortDim) + indexBox.end(sortDim)) >> 1;
			auto indexBoxL = indexBox;
			auto indexBoxR = indexBox;
			auto particleRangeL = particleRange;
			auto particleRangeR = particleRange;
			size_t const midIndex = std::midpoint(indexBox.begin(sortDim), indexBox.end(sortDim));
			indexBoxL.end(sortDim) = indexBoxR.begin(sortDim) = midIndex;
			if (particleRange.second > particleRange.first) {
				real_type const splitPosition = real_type(splitIndex) * gAttr.gridSpacing[sortDim];
				size_t hi = particleRange.second;
				size_t lo = particleRange.first;
				auto const &sortPos = X[sortDim];
				while (lo < hi) {
					if (sortPos[lo] >= splitPosition) {
						while (lo != hi) {
							hi--;
							if (sortPos[hi] < splitPosition) {
								swap_particles(hi, lo);
								break;
							}
						}
					}
					lo++;
				}
				particleRangeL.second = particleRangeR.first = hi;
			}
			sort_particles(particleRangeL, indexBoxL);
			sort_particles(particleRangeR, indexBoxR);
		}

	}
	void sort_particles() {
		Interval<size_t, NDIM> indexBox( { 0, 0, 0 }, { gAttr.intSizes[0], gAttr.intSizes[1], gAttr.intSizes[2] });
		sort_particles(std::pair<size_t, size_t>(0, X[0].size()), indexBox);
	}
	Vector<real_type, DIM4> getX(size_t n) const {
		Vector<real_type, DIM4> x;
		for (int dim = 0; dim < NDIM; dim++) {
			x[dim + 1] = X[dim][n];
		}
		x[0] = time;
		return x;
	}
	Vector<real_type, DIM4> getU(size_t n) const {
		Vector<real_type, DIM4> u;
		real_type W2 = one;
		for (int dim = 0; dim < NDIM; dim++) {
			real_type const v = U[dim][n];
			u[dim + 1] = v;
			W2 += v * v;
		}
		u[0] = sqrt(W2);
		return u;
	}
	void setU(size_t n, Vector<real_type, DIM4> const &u) {
		for (int dim = 0; dim < NDIM; dim++) {
			U[dim][n] = u[dim + 1];
		}
	}
	void setX(size_t n, Vector<real_type, DIM4> const &x) {
		for (int dim = 0; dim < NDIM; dim++) {
			X[dim][n] = x[dim + 1];
		}
	}
	auto stressEnergyTensor() const {
		using namespace Math;
		auto const opts = getOptions();
		int const bw = gAttr.boundWidth;
		Real const dx3inv = gAttr.d3Rinv;
		Real const dx3 = gAttr.d3R;
		SymmetricMatrix<grid_type, DIM4> T(grid_type(zero, gAttr.extSize));
		for (size_t n = 0; n != N; n++) {
			Vector<int, NDIM> I;
			Vector<real_type, DIM4> u = getU(n);
			Vector<real_type, DIM4> x = getX(n);
			SymmetricMatrix<real_type, DIM4> dT;
			for (int j = 0; j < DIM4; j++) {
				for (int k = 0; k <= j; k++) {
					dT[j, k] = M * u[j] * u[k] * dx3inv;
				}
			}
			for (int k = 0; k < NDIM; k++) {
				I[k] = int(x[k + 1] * gAttr.intSizes[k]) + bw;
			}
			for (int ix = I[0] - 1; ix <= I[0] + 1; ix++) {
				for (int iy = I[1] - 1; iy <= I[1] + 1; iy++) {
					for (int iz = I[2] - 1; iz <= I[2] + 1; iz++) {
						auto const K = Vector<int, NDIM>( { ix, iy, iz });
						real_type weight = one;
						int index = 0;
						for (int k = 0; k < NDIM; k++) {
							index += gAttr.extStrides[k] * periodicModulus(K[k], gAttr.intSizes[k]);
							real_type const X = x[k + 1] * real_type(gAttr.intSizes[k]) - real_type(K[k] - bw) - half;
							weight *= kernelTSC(X);
						}
						for (int j = 0; j < DIM4; j++) {
							for (int k = 0; k <= j; k++) {
								T[j, k][index] += weight * dT[j, k];
							}
						}
					}
				}
			}
		}
		for (int j = 0; j < DIM4; j++) {
			for (int k = 0; k <= j; k++) {
				auto const tmp = T[j, k].sum() * dx3;
				T[j, k] -= tmp;
			}
		}
		return T;
	}
	template<typename MetricFunctor>
	void kick(MetricFunctor const &grVars, real_type dt) {
		for (size_t n = 0; n != N; n++) {
			Vector<real_type, DIM4> du = zero;
			auto u = getU(n);
			auto const x = getX(n);
			auto const D = grVars(x[1], x[2], x[3]);
			for (int k = 0; k < DIM4; k++) {
				for (int i = 0; i < DIM4; i++) {
					for (int j = 0; j <= i; j++) {
						du[k] -= half * D[k][i, j] * u[i] * u[j];
					}
				}
			}
			u += du * dt;
			setU(n, u);
		}

	}
	void drift(real_type dt) {
		for (size_t n = 0; n != N; n++) {
			Vector<real_type, DIM4> dx = zero;
			auto const u = getU(n);
			auto x = getX(n);
			real_type const Winv = one / u[0];
			for (int k = 1; k < DIM4; k++) {
				dx[k] += u[k] * Winv;
			}
			x += dx * dt;
			for (int dim = 1; dim < DIM4; dim++) {
				while (x[dim] >= one) {
					x[dim] -= one;
				}
				while (x[dim] < zero) {
					x[dim] += one;
				}
			}
			setX(n, x);
		}
	}
	void output(DBfile *db, DBoptlist *optList) const {
		static constexpr int silo_data_type = DB_DOUBLE;
		static constexpr char meshName[] = "particleMesh";
		real_type const *coords[NDIM] = { &X[XDIM][0], &X[YDIM][0], &X[ZDIM][0] };
		DBPutPointmesh(db, meshName, NDIM, coords, N, silo_data_type, optList);
		for (int dim = 0; dim < NDIM; dim++) {
			auto const name = std::string("u") + std::to_string(dim);
			DBPutPointvar1(db, name.c_str(), meshName, std::begin(U[dim]), N, silo_data_type, optList);
		}
	}
private:
	GridAttributes<real_type> const gAttr;
	Vector<grid_type, NDIM> X;
	Vector<grid_type, NDIM> U;
	std::valarray<std::pair<size_t, size_t>> particleGrid;
	real_type const M;
	real_type time;
	size_t N;
}
;

#endif /* INCLUDE_PARTICLES_HPP_ */
