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
#include "Quadrature.hpp"
#include "RungeKutta.hpp"

#include <valarray>

template<typename, int, int>
struct HydroGrid;

using namespace Hydrodynamics;
using namespace Math;

template<int D>
struct MultiIndex {
	MultiIndex(Vector<size_t, D> const &b, Vector<size_t, D> const &e) {
		start = b;
		stop = e;
		indices = start;
		endIndices[0] = stop[0];
		for (int k = 1; k < D; k++) {
			endIndices[k] = 0;
		}
	}
	MultiIndex(Vector<size_t, D> const &e) {
		MultiIndex(zeroVector<size_t, D>(), e);
	}
	MultiIndex(size_t const &e) {
		MultiIndex(zeroVector<size_t, D>(), Vector<size_t, D>(e));
	}
	MultiIndex& operator++() {
		int k = D - 1;
		if (indices[0] + 1 < stop[0]) {
			indices[k]++;
			while ((k > 0) && (indices[k] == stop[k])) {
				indices[k] = start[k];
				k--;
				indices[k]++;
			}
		}
	}
	MultiIndex operator++(int) {
		auto const rc = *this;
		operator++();
		return rc;
	}
	Vector<size_t, D> begin() const {
		return start;
	}
	Vector<size_t, D> end() const {
		return endIndices;
	}
	Vector<size_t, D> getIndices() const {
		return indices;
	}
	int indexAt(int k) const {
		return indices[k];
	}
	void setIndices(Vector<size_t, D> const &indices_) {
		indices = indices_;
	}
	operator int() const {
		size_t stride = 1;
		size_t index = 0;
		for (int k = D - 1; k >= 0; k--) {
			index += stride * indices[k];
			stride *= (stop[k] - start[k]);
		}
		return index;
	}
private:
	Vector<size_t, D> start;
	Vector<size_t, D> endIndices;
	Vector<size_t, D> stop;
	Vector<size_t, D> indices;
};

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
	using triindices_m1_t = TriangularIndices<D, P - 1>;
	using array_t = std::vector<std::vector<T>>;
	using inner_array_t = std::vector<T>;
	using size_array_t = Vector<size_t, D>;
	using state_t = ConservedState<T, D>;
	using primitive_state_t = PrimitiveState<T, D>;
	using array_state_t = ConservedState<array_t, D>;
	using inner_array_state_t = ConservedState<inner_array_t, D>;
	static T constexpr zero = T(0), half = T(0.5), one = T(1), two = T(2);
	static int constexpr NF = state_t::NFields;
	static int constexpr P3 = triindices_t::size();
	static int constexpr Q = P + 1;
	static int constexpr Q3 = Math::integerPower(Q, D);
	static auto constexpr slopeLimiter = vanLeer<T>;
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
			S[f] = std::vector<std::vector<T>>(P3, std::vector<T>(N3));
			V[f] = std::vector<std::vector<T>>(Q3, std::vector<T>(N3));
			F[f] = std::vector<std::vector<T>>(D, std::vector<T>(N3));
		}
		bottomCorner = vectorSum(gridStrides);
		topCorner = N3 - vectorSum(gridStrides);
	}
	void reconstruct() {
		for (int f = 0; f < NF; f++) {
			auto &u = U[f];
			auto &v = V[f];
			for (int n = 0; n < N; n++) {
				for (MultiIndex<D> q3(Q); q3 != q3.end(); q3++) {
					Vector<T, D> x;
					for (int dim = 0; dim < D; dim++) {
						x[dim] = quadratureRules.x[q3.indexAt(dim)];
					}
					auto const basis = legendreBasis<T, D, P>(x);
					T value = zero;
					for (TriangularIndices<D, P> p3; !p3.end(); p3++) {
						value += u[p3][n] * basis[p3];
					}
					v[q3][n] = value;
				}
			}
		}
	}
	void toCharacteristics(int dim) {
		for (int n = 0; n < N; n++) {
			auto u = getStateAt(n);
			auto const R = matrixInverse(u[0].eigenVectors(dim));
			for (int p3 = 1; p3 < P3; p3++) {
				u[p3] = R * u[0];
			}
			for (int f = 0; f < NF; f++) {
				for (int p3 = 0; p3 < P3; p3++) {
					U[f][p3][n] = u[p3][f];
				}
			}
			setStateAt(n, u);
		}
	}
	void fromCharacteristics(int dim) {
		for (int n = 0; n < N; n++) {
			auto u = getStateAt(n);
			auto const R = u[0].eigenVectors(dim);
			for (int p3 = 1; p3 < P3; p3++) {
				u[p3] = R * u[0];
			}
			for (int f = 0; f < NF; f++) {
				for (int p3 = 0; p3 < P3; p3++) {
					U[f][p3][n] = u[p3][f];
				}
			}
			setStateAt(n, u);
		}
	}
	Vector<state_t, P3> getFieldAt(int field, int n) const {
		Vector<T, P3> u;
		for (int p3 = 0; p3 < P3; p3++) {
			u[p3] = U[field][p3][n];
		}
		return u;
	}
	Vector<state_t, P3> getStateAt(int n) const {
		Vector<state_t, P3> u;
		for (int f = 0; f < NF; f++) {
			for (int p3 = 0; p3 < P3; p3++) {
				u[p3][f] = U[f][p3][n];
			}
		}
		return u;
	}
	void setFieldAt(int field, int n, Vector<T, P3> const &u) const {
		for (int p3 = 0; p3 < P3; p3++) {
			U[field][p3][n] = u[p3];
		}
	}
	void setStateAt(int n, Vector<state_t, P3> const &u) const {
		for (int f = 0; f < NF; f++) {
			for (int p3 = 0; p3 < P3; p3++) {
				U[f][p3][n] = u[p3][f];
			}
		}
	}
	void project(int dim) {
		int const dN = gridStrides[dim];
		for (int n = bottomCorner; n < topCorner; n++) {
			int const O3 = triindices_m1_t::size();
			auto const u0 = getStateAt(n);
			auto const pDelta = getStateAt(n + dN) - u0;
			auto const mDelta = u0 - getStateAt(n - dN);
			Vector<T, P3> delta;
			for (int p3 = 0; p3 < O3; p3++) {
				delta[p3] = slopeLimiter(pDelta[p3], mDelta[p3]);
			}
			std::bitset < P3 > notLimited(false);
			for (triindices_t pp = --triindices_t::end(); pp != triindices_t::begin(); pp--) {
				if (pp[dim] + 1 < P) {
					triindices_t const p0 = pp;
					p0[dim]--;
					T const norm = one / (two * p0[dim] + one);
					T const delta = slopeLimiter(pDelta[p0], mDelta[p0]);
					u0[pp] = minmod(u0[pp], delta[p0]);
				}
			}
			setStateAt(n, u0);
		}

	}
	void computeFlux(int dim) {
		for (int n = bottomCorner; n < topCorner; n++) {
			state_t flux(T(0));
			for (MultiIndex<D - 1> q2(Q); q2 != q2.end(); q2++) {
				size_t rightIndex, leftIndex;
				Vector<size_t, D> rightIndices;
				Vector<size_t, D> leftIndices;
				state_t leftState, rightState;
				const auto q2Indices = q2.getIndices();
				leftIndices[dim] = 0;
				rightIndices[dim] = P - 1;
				for (int l = 0; l < dim; l++) {
					rightIndices[l] = leftIndices[l] = q2Indices[l];
				}
				for (int l = dim; l < D - 1; l++) {
					rightIndices[l + 1] = leftIndices[l + 1] = q2Indices[l];
				}
				leftIndex = MultiIndex<D>().setIndices(leftIndices);
				rightIndex = MultiIndex<D>().setIndices(rightIndices);
				for (int f = 0; f < NF; f++) {
					rightState[f] = V[f][rightIndex][n];
					leftState[f] = V[f][leftIndex][n];
				}
				T const weight = half * quadratureRules.w[leftIndex];
				flux += weight * riemannSolver(leftState, rightState, dim);
			}
			for (int f = 0; f < NF; f++) {
				F[f][dim][n] = flux;
			}
		}
	}
	void computeSource() {
		for (int n = bottomCorner; n < topCorner; n++) {
			Vector<state_t, P3> thisSource = zero;
			for (MultiIndex<D> q3(Q); q3 != q3.end(); q3++) {
				Vector<Vector<T, P3>, D> dBasisDx;
				Vector<T, D> x;
				state_t v;
				T w = one;
				for (int f = 0; f < NF; f++) {
					v[f] = U[f][q3][n];
				}
				for (int dim = 0; dim < D; dim++) {
					dBasisDx[dim] = legendreBasis<T, D, P>(x, dim);
				}
				for (int dim = 0; dim < D; dim++) {
					int const q = q3.indexAt(dim);
					w *= half * quadratureRules.w[q];
					x[dim] = quadratureRules.x[q];
				}
				for (int dim = 0; dim < D; dim++) {
					auto const flux = primState(v[q3]).toFlux(dim);
					for (int p3 = 0; p3 < P3; p3++) {
						thisSource[p3] += w * flux * dBasisDx[dim];
					}
				}
			}
			for (int f = 0; f < NF; f++) {
				for (int p3 = 0; p3 < P3; p3++) {
					S[f][p3][n] = thisSource[p3][f];
				}
			}
		}
	}
private:
	static RungeKutta<T, P> constexpr butcherTable = RungeKutta<T, P>();
	static Quadrature::QuadratureRules<T, Q> const quadratureRules;
	GridOptions<T> gridOptions;
	int N;
	int N3;
	int BW;
	int topCorner;
	int bottomCorner;
	T dx;
	Vector<int, D> gridDims;
	Interval<T, D> gridDomain;
	size_array_t gridStrides;
	array_state_t U;
	array_state_t V;
	array_state_t F;
	array_state_t S;
};

template<typename T, int D, int P>
Quadrature::QuadratureRules<T, HydroGrid<T, D, P>::Q> const HydroGrid<T, D, P>::quadratureRules = Quadrature::gaussLegendreRules<HydroGrid<T, D, P>::Q>();

#endif /* INCLUDE_HYDROGRID_HPP_ */
