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
#include "Kernels.hpp"

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

struct Particles {
	using real_type = Real;
	using grid_type = std::valarray<real_type>;
	static constexpr real_type zero = real_type(0), half = real_type(0.5), one = real_type(1);
	Particles(size_t count) :
			X(grid_type(zero, count)), U(grid_type(zero, count)), M(getOptions().totalMass / real_type(count)), time(0.0), N(count) {
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
	auto stressEnergyTensor(GridAttributes<real_type> const &gAttr) const {
		using namespace Math;
		int const bw = gAttr.boundWidth;
		SymmetricMatrix<grid_type, NDIM + 1> T(grid_type(zero, gAttr.extSize));
		for (size_t n = 0; n != N; n++) {
			Vector<int, NDIM> I;
			Vector<real_type, DIM4> u = getU(n);
			Vector<real_type, DIM4> x = getX(n);
			Vector<real_type, NDIM + 1> dT;
			for (int j = 0; j < DIM4; j++) {
				for (int k = 0; k <= j; k++) {
					dT[j, k] = M * u[j] * u[k];
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
							index += gAttr.gridStrides[k] * periodicModulus(K[k], gAttr.intSizes[k]);
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
		return T;
	}
	template<typename MetricFunctor>
	void kick(MetricFunctor const &dgdx, real_type dt) {
		for (size_t n = 0; n != N; n++) {
			Vector<real_type, DIM4> du = zero;
			auto &u = getU(n);
			auto const x = getX(n);
			auto const D = dgdx(x[0], x[1], x[2]);
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
			for (int k = 0; k < NDIM; k++) {
				for (int i = 0; i < DIM4; i++) {
					for (int j = 0; j <= i; j++) {
						dx[k] += u[k] * Winv;
					}
				}
			}
			x += dx * dt;
			setX(n, x);
		}
	}
private:
	Vector<grid_type, NDIM> X;
	Vector<grid_type, NDIM> U;
	real_type const M;
	real_type time;
	size_t N;
};

#endif /* INCLUDE_PARTICLES_HPP_ */
