#pragma once

#include "Basis.hpp"
#include "Hdf5.hpp"
#include "RungeKutta.hpp"
#include "Vector.hpp"

template<typename T>
inline constexpr T minmod(T const &a, T const &b) {
	static constexpr T half = T(0.5);
	T const sgn = std::copysign(half, a) + std::copysign(half, b);
	T const mag = std::min(std::abs(a), std::abs(b));
	return sgn * mag;
}

template<typename T, int D, int N, int P, typename HyperState>
struct HyperGrid {
	static constexpr Basis<T, D, P> basis { };
	static constexpr RungeKutta<T, P> butcherTable { };
	static constexpr T zero = T(0);
	static constexpr T half = T(0.5);
	static constexpr T one = T(1);
	static constexpr T two = T(2);
	static constexpr int BW = 2;
	static constexpr int NN = N + 2 * BW;
	using hindex_type = IndexTuple<NN, D>;
	using pindex_type = TriangularIndex<P, D>;
	using basis_type = Basis<T, D, P>::basis_type;
	using face_type = Basis<T, D, P>::face_type;
	static constexpr int N3 = hindex_type::elementCount();
	static constexpr int P3 = pindex_type::elementCount();
	static constexpr int Q3 = IndexTuple<P, D>::elementCount();
	static constexpr int NF = HyperState::fieldCount();
	static constexpr int NS = butcherTable.s;
	using initialize_type = HyperState(Math::Vector<T, D>);
	static constexpr auto stride = []() {
		std::array<int, D> stride;
		stride[D - 1] = 1;
		for (int d = D - 2; d >= 0; d--) {
			stride[d] = stride[d + 1] * (N + 2 * BW);
		}
		return stride;
	}();
	static constexpr auto interiorCells(std::array<int, D> const &lb, std::array<int, D> const &ub) {
		static constexpr int M3 = Math::integerPower(N, D - 1) * (N + 1);
		std::array<hindex_type, M3> cells;
		for (int dim = 0; dim < D; dim++) {
			int nextIndex = 0;
			for (auto I = hindex_type::begin(); I != hindex_type::end(); I++) {
				bool flag = true;
				for (int d = 0; d < D; d++) {
					if ((I[d] < lb[d]) || (I[d] >= ub[d])) {
						flag = false;
						break;
					}
				}
				if (flag) {
					cells[nextIndex++] = I;
				}
			}
		}
		return cells;
	}
	static constexpr auto faceCells(int dim) {
		std::array<int, D> lb, ub;
		std::fill(lb.begin(), lb.end(), BW);
		std::fill(ub.begin(), ub.end(), N + BW);
		ub[dim]++;
		return interiorCells(lb, ub);
	}
	static constexpr auto interiorCells(int expand = 0) {
		std::array<int, D> lb, ub;
		std::fill(lb.begin(), lb.end(), BW - expand);
		std::fill(ub.begin(), ub.end(), N + BW + expand);
		return interiorCells(lb, ub);
	}
	HyperGrid(T span = one) :
			currentTime(zero), cellWidth(span / T(N)), stepNumber(0), U(NF) {
	}
	void initialize(Math::Vector<T, D> const &O, initialize_type const &F) {
		static auto const cellIndices = interiorCells();
		origin = O;
		T const h = cellWidth;
		for (auto const &I : cellIndices) {
			auto const hi = I.flatIndex();
			auto const f = [this, h, F, I](Math::Vector<T, D> x) {
				for (int d = 0; d < D; d++) {
					x[d] = h * (origin[d] + (I[d] - BW + half) + half * x[d]);
				}
				return F(x);
			};
			auto const u = basis.template volumeTransform<HyperState>(f);
			for (int fi = 0; fi < NF; fi++) {
				for (int pi = 0; pi < P3; pi++) {
					U[fi][pi][hi] = u[pi][fi];
				}
			}
		}
	}
	auto computeSource() {
		static constexpr auto cellIndices = interiorCells();
		static constexpr auto qpts = basis.getQuadraturePoints();
		std::array<std::array<HyperState, N3>, P3> S;
		for (auto const &hIdx : cellIndices) {
			auto const hi = hIdx.flatIndex();
			std::array<HyperState, Q3> V;
			std::array<std::array<T, Q3>, NF> F;
			for (int fi = 0; fi < NF; fi++) {
				basis_type W;
				for (int pi = 0; pi < P3; pi++) {
					W[pi] = U[fi][pi][hi];
				}
				auto const tmp = basis.inverseVolumeTransform(W);
				for (int qi = 0; qi < Q3; qi++) {
					V[qi][fi] = tmp[qi];
				}
			}
			for (int qi = 0; qi < Q3; qi++) {
				for (int d = 0; d < D; d++) {
					auto const f = V[qi].flux(d);
					for (int fi = 0; fi < NF; fi++) {
						F[fi][qi] = f[fi];
					}
				}
			}
			for (int fi = 0; fi < NF; fi++) {
				auto thisS = basis.weakDivergence(F[fi]);
				for (int pi = 0; pi < P3; pi++) {
					S[hi][pi][fi] = thisS[pi][fi];
				}
			}
		}
		return S;
	}
	auto computeFlux(int dim) {
		static constexpr auto cellIndices = faceCells(dim);
		static constexpr int nSurface = basis.Nf;
		static constexpr int nVolume = basis.Np;
		std::array<HyperState, N3> F;
		auto const str = stride[dim];
		for (auto const &indicesR : cellIndices) {
			auto indicesL = indicesR;
			indicesL[dim]--;
			auto const iL = indicesL.flatIndex();
			auto const iR = indicesR.flatIndex();
			std::array<HyperState, nSurface> uL;
			std::array<HyperState, nSurface> uR;
			std::array<HyperState, nSurface> fRiemann;
			for (int fi = 0; fi < NF; fi++) {
				basis_type phiR, phiL;
				for (int pi = 0; pi < nVolume; pi++) {
					phiL[pi] = U[fi][pi][iL];
					phiR[pi] = U[fi][pi][iR];
				}
				face_type const vL = basis.inverseFaceTransform(phiL);
				face_type const vR = basis.inverseFaceTransform(phiR);
				for (int i = 0; i < nSurface; i++) {
					uR[i][fi] = vR[i];
					uL[i][fi] = vL[i];
				}
			}
			for (int i = 0; i < nSurface; i++) {
				fRiemann[i] = HyperState::riemann(uL[i], uR[i], dim);
			}
			F[iR] = half * basis.faceIntegrate(fRiemann);
		}
		return F;
	}
	void applyLimiter() {
		static constexpr auto cellIndices = interiorCells(1);
		for (int pOrder = P - 1; pOrder > 0; pOrder--) {
			for (auto P1 = pindex_type::begin(); P1 != pindex_type::end(); P1++) {
				auto const Q1 = triangularToLegendre(P1);
				if (Q1[0] == pOrder) {
					for (auto const &hIndices : cellIndices) {
						auto const hFlat = hIndices.flatIndex();
						for (int dim = 0; dim < D; dim++) {
							if (Q1[dim]) {
								HyperState U1, Up, Um;
								auto const y = (two * Q1[dim] - one);
								auto const alpha = (two * Q1[dim] + one) / (y * y);
								auto const str = stride[dim];
								auto Q0 = Q1;
								Q0[dim]--;
								auto const P0 = legendreToTriangular(Q0);
								auto const p0 = P0.flatIndex();
								auto const p1 = P1.flatIndex();
								for (int field = 0; field < NF; field++) {
									auto const u0 = U[field][p0][hFlat];
									auto const &u1 = U[field][p1][hFlat];
									auto const &up = U[field][p0][hFlat + str];
									auto const &um = U[field][p0][hFlat - str];
									U1[field] = u1;
									Up[field] = up - u0;
									Um[field] = u0 - um;
								}
								U1 = U1.toCharacteristic(dim);
								Up = Up.toCharacteristic(dim);
								Um = Um.toCharacteristic(dim);
								for (int field = 0; field < NF; field++) {
									auto const U0 = U1[field];
									auto const Ud = alpha * minmod(Up[field], Um[field]);
									U1[field] = minmod(U0, Ud);
								}
								U1 = U1.fromCharacteristic(dim);
								for (int field = 0; field < NF; field++) {
									U[field][p1][hFlat] = U1[field];
								}
							}
						}
					}
				}
			}
		}
	}
	auto computeTimeDerivatives() {
		static constexpr auto cellIndices = interiorCells(1);
		T const dx = cellWidth;
		T const dxinv = one / dx;
		std::array<std::array<std::array<T, N3>, P3>, NF> dUdt;
		{
			auto const S = computeSource();
			for (auto const &idx : cellIndices) {
				auto const i = idx.flatIndex();
				for (auto qIdx = pindex_type::begin(); qIdx != pindex_type::end(); qIdx++) {
					auto const qi = qIdx.flatIndex();
					for (int fi = 0; fi < NF; fi++) {
						dUdt[fi][qi][i] = S[qi][i][fi];
					}
				}
			}
		}
		for (int dim = 0; dim < D; dim++) {
			auto const F = computeFlux(dim);
			auto const str = stride[dim];
			for (auto const &idx : cellIndices) {
				auto const i0 = idx.flatIndex();
				auto const ip = i0 + str;
				for (auto tQ = pindex_type::begin(); tQ != pindex_type::end(); tQ++) {
					auto const qi = tQ.flatIndex();
					auto const lQ = triangularToLegendre(tQ);
					int const sgn = T(2 * (1 & lQ[dim]) - 1);
					for (int fi = 0; fi < NF; fi++) {
						dUdt[fi][qi][i0] -= dxinv * (F[ip][fi] - sgn * F[i0][fi]);
					}
				}
			}
		}
		return dUdt;
	}
	void substep(int i, T dt) {
		if (i == 0) {
			dUdt.resize(NF);
			return;
		}
		U0 = U;
		static constexpr auto cellIndices = interiorCells(1);
		applyLimiter();
		dUdt[i - 1] = computeTimeDerivatives();
		auto const h = dt / cellWidth;
		for (auto const &hIdx : cellIndices) {
			auto const hi = hIdx.flatIndex();
			for (auto pIdx = pindex_type::begin(); pIdx != pindex_type::end(); pIdx++) {
				auto const pi = pIdx.flatIndex();
				for (int fi = 0; fi < NF; fi++) {
					auto &Yi = U[fi][pi][hi];
					auto const &yn = U0[fi][pi][hi];
					Yi = yn;
					for (int j = 0; j < i; j++) {
						auto const &a_ij = butcherTable.a[i - 1][j];
						auto const &dfdx = dUdt[j][fi][pi][hi];
						Yi += h * a_ij * dfdx;
					}
				}
			}
		}
		if (i == NS - 1) {
			for (auto const &hIdx : cellIndices) {
				auto const hi = hIdx.flatIndex();
				for (auto pIdx = pindex_type::begin(); pIdx != pindex_type::end(); pIdx++) {
					auto const pi = pIdx.flatIndex();
					for (int fi = 0; fi < NF; fi++) {
						auto &ynp1 = U[fi][pi][hi];
						auto const &yn = U0[fi][pi][hi];
						ynp1 = yn;
						for (int j = 0; j < NS; j++) {
							auto const &b_i = butcherTable.b[j];
							auto const &dfdx = dUdt[j][fi][pi][hi];
							ynp1 += h * b_i * dfdx;
						}
					}
				}
			}
			stepNumber++;
			currentTime += dt;
			decltype(U0)().swap(U0);
			decltype(dUdt)().swap(dUdt);
		}
	}
	void output(std::string filename) const {
		//typename T, int D, int N, int P3>
		filename += std::string(".") + std::to_string(stepNumber) + std::string(".hdf5");
		static const auto fnames = HyperState::getFieldNames();
		writeHdf5<T, D, N, P, BW>(filename, cellWidth, U, fnames);
	}
private:
	T currentTime;
	T cellWidth;
	int stepNumber;
	Math::Vector<T, D> origin;
	std::vector<std::array<std::array<std::array<T, N3>, P3>, NF>> dUdt;
	std::vector<std::array<std::array<T, N3>, P3>> U;
	std::vector<std::array<std::array<T, N3>, P3>> U0;
};

