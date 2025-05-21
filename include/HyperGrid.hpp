#pragma once

#include "Basis.hpp"
#include "ContainerArithmetic.hpp"
#include "Quadrature.hpp"

#include <numeric>

template<typename T>
inline constexpr T minmod(T const &a, T const &b) {
	static constexpr T half = T(0.5);
	T const sgn = std::copysign(half, a) + std::copysign(half, b);
	T const mag = std::min(std::abs(a), std::abs(b));
	return sgn * mag;
}

template<typename State, int Span, int Order, typename RungeKutta>
class HyperGrid {
	using T = typename State::value_type;
	using triindex_t = TriIndex<Order, State::dimCount()>;
	static constexpr int D = State::dimCount();
	static constexpr int BW = 2;
	static constexpr int NF = State::fieldCount();
	static constexpr int NS = RungeKutta::stageCount();
	static constexpr int P = Order;
	static constexpr int H = P;
	static constexpr int N = Span + 2 * BW;
	static constexpr int P3 = triindex_t::size();
	static constexpr int H3 = ipow(H, D);
	static constexpr int N3 = ipow(N, D);
	static constexpr Range<int, D> oBox { repeat<D>(-BW), repeat<D>(Span + BW) };
	std::vector<std::array<std::array<std::array<T, N3>, P3>, NF>> k_;
	std::vector<std::array<std::array<T, N3>, P3>> U_;
	std::vector<std::array<std::array<T, N3>, P3>> Un_;
	const T dx;
	const T dxinv;
	static constexpr int stride(int d) {
		return ipow(N, D - 1 - d);
	}
public:
	HyperGrid(T const &xSpan = T(1)) :
			U_(NF), dx(xSpan / T(Span)), dxinv(T(Span) / xSpan) {
	}
	void initialize(std::function<State(std::array<T, D> const&)> const &f) {
		using quad_type = Quadrature<GaussLegendreQuadrature<T, P>, D>;
		constexpr Range<int, D> iBox { repeat<D>(0), repeat<D>(N) };
		constexpr int NH = quad_type::size();
		quad_type const quadrature;
		constexpr Basis<T, P, D> basis;
		using index_type = MultiIndex<oBox, iBox>;
		T const hdx = T(0.5) * dx;
		for (auto I = index_type::begin(); I != index_type::end(); I++) {
			int const ii = I;
			for (int pi = 0; pi < P3; pi++) {
				for (int hi = 0; hi < NH; hi++) {
					auto const phi = basis(quadrature.point(hi));
					T const w = quadrature.weight(hi);
					auto const x0 = quadrature.point(hi);
					std::array<T, D> x;
					for (int d = 0; d < D; d++) {
						x[d] = (T(2 * I[d] + 1) + x0[d]) * hdx;
					}
					auto const dU = w * phi[pi] * f(x);
					for (int fi = 0; fi < NF; fi++) {
						U_[fi][pi][ii] += dU[fi];
					}
				}
			}
		}
	}
	T beginStep() {
		using quad_type = Quadrature<GaussLegendreQuadrature<T, P>, D>;
		constexpr Range<int, D> iBox { repeat<D>(0), repeat<D>(N) };
		constexpr int NH = quad_type::size();
		const quad_type quadrature;
		const RungeKutta rk;
		using index_type = MultiIndex<oBox, iBox>;
		constexpr Basis<T, P, D> basis;
		auto const norm = basis.norm();
		Un_ = U_;
		k_.resize(NS);
		std::array<T, D> maxLambda { T(0) };
		for (auto I = index_type::begin(); I != index_type::end(); I++) {
			int const ii = I;
			for (int hi = 0; hi < NH; hi++) {
				State u;
				auto const phi = basis(quadrature.point(hi));
				for (int fi = 0; fi < NF; fi++) {
					u[fi] = T(0);
					for (int pi = 0; pi < P3; pi++) {
						u[fi] += U_[fi][pi][ii] * norm[pi] * phi[pi];
					}
				}
				auto const max = [](T a, T b) {
					return std::max(a, b);
				};
				for (int d = 0; d < D; d++) {
					auto const lambda = u.eigenValues(d);
					maxLambda[d] = std::accumulate(lambda.begin(), lambda.end(), maxLambda[d], max);
				}
			}
		}
		T const num = dx * rk.cfl();
		T const den = T(2 * P - 1) * std::accumulate(maxLambda.begin(), maxLambda.end(), T(0));
		return num / den;
	}
	void subStep(T const &dt, int i) {
		constexpr Range<int, D> iBox { repeat<D>(0), repeat<D>(N) };
		RungeKutta const rk;
		using index_type = MultiIndex<oBox, iBox>;
		U_ = Un_;
		for (int j = 0; j < i; j++) {
			for (int fi = 0; fi < NF; fi++) {
				for (int pi = 0; pi < P3; pi++) {
					for (auto I = index_type::begin(); I != index_type::end(); I++) {
						int const ii = I;
						U_[fi][pi][ii] += rk.a(i, j) * k_[j][fi][pi][ii];
					}
				}
			}
		}
		auto &dUdt = k_[i];
		T const lambda = dt * dxinv;
		applyLimiter();
		for (int fi = 0; fi < NF; ++fi) {
			applyFlux(fi, lambda, dUdt, std::make_integer_sequence<int, D> { });
		}
		auto const S = computeSource();
		for (int fi = 0; fi < NF; fi++) {
			for (auto I = index_type::begin(); I != index_type::end(); I++) {
				int const ii = I;
				for (auto P = triindex_t::begin(); P != triindex_t::end(); P++) {
					int const pi = P;
					dUdt[fi][pi][ii] -= lambda * S[fi][pi][ii];
				}
			}
		}
	}
	void endStep() {
		constexpr Range<int, D> iBox { repeat<D>(0), repeat<D>(N) };
		using index_type = MultiIndex<oBox, iBox>;
		RungeKutta const rk;
		for (int fi = 0; fi < NF; fi++) {
			for (auto P = triindex_t::begin(); P != triindex_t::end(); P++) {
				int const pi = P;
				for (auto I = index_type::begin(); I != index_type::end(); I++) {
					int const ii = I;
					U_[fi][pi][ii] = Un_[fi][pi][ii];
					for (int i = 0; i < NS; i++) {
						U_[fi][pi][ii] += rk.b(i) * k_[i][fi][pi][ii];
					}
				}
			}
		}
		k_ = { };
		Un_ = { };
	}
	void applyLimiter() {
		static constexpr auto alpha = []() {
			std::array<T, P> a;
			for (int n = 0; n < P; n++) {
				a[n] = T(1) / T(2 * n + 1);
			}
			return a;
		}();
		constexpr Range<int, D> iBox { repeat<D>(-1), repeat<D>(N + 1) };
		using index_type = MultiIndex<oBox, iBox>;
		for (int degree = Order - 1; degree > 0; degree--) {
			for (auto P1 = triindex_t::begin(degree); P1 != triindex_t::end(degree); P1++) {
				for (int d = 0; d < D; d++) {
					if (P1[d] == 0) {
						continue;
					}
					triindex_t P0 = P1.dec(d);
					int const p0 = P0;
					int const p1 = P1;
					int const di = stride(d);
					for (auto I = index_type::begin(); I != index_type::end(); I++) {
						int const ii = I;
						State u1, u0, up, um;
						for (int fi = 0; fi < NF; fi++) {
							u1[fi] = U_[fi][p1][ii];
							u0[fi] = U_[fi][p0][ii];
							up[fi] = U_[fi][p0][ii + di];
							um[fi] = U_[fi][p0][ii - di];
						}
						u1 = u1.toCharacteristic(d);
						u0 = u0.toCharacteristic(d);
						up = up.toCharacteristic(d);
						um = um.toCharacteristic(d);
						T const a = alpha[P0[d]];
						for (int fi = 0; fi < NF; fi++) {
							u1[fi] = minmod(u1[fi], a * minmod(up[fi] - u0[fi], u0[fi] - um[fi]));
						}
						u1 = u1.fromCharacteristic(d);
						for (int fi = 0; fi < NF; fi++) {
							U_[fi][p1][ii] = u1[fi];
						}
					}
				}
			}
		}
	}
	template<int DIM>
	std::vector<std::array<T, N3>> computeFlux() const {
		static_assert(DIM < D);
		using quad_type = Quadrature<GaussLegendreQuadrature<T, P>, D - 1>;
		using index_type = quad_type::index_type;
		quad_type const quadrature;
		constexpr int NH = quad_type::size();
		constexpr Basis<T, P, D> basis;
		auto const norm = basis.norm();
		std::vector<std::array<T, N3>> F(NF);
		for (auto I = index_type::begin(); I != index_type::end(); I++) {
			int const i0 = I;
			int const &iL = i0;
			int const iR = i0 + stride(DIM);
			State uL, uR;
			for (int hi = 0; hi < NH; hi++) {
				auto const x = quadrature.point(hi);
				auto const phiL = basis(insert<D>(T(+1), DIM, x));
				auto const phiR = basis(insert<D>(T(-1), DIM, x));
				for (int fi = 0; fi < NF; fi++) {
					uL[fi] = uR[fi] = T(0);
					for (int pi = 0; pi < P3; pi++) {
						uL[fi] += norm[pi] * phiL[pi] * U_[fi][pi][iL];
						uR[fi] += norm[pi] * phiR[pi] * U_[fi][pi][iR];
					}
				}
				State const flux = solveRiemannProblem(uL, uR, DIM);
				for (int fi = 0; fi < NF; fi++) {
					F[fi][i0] = flux[fi];
				}
			}
		}
		return F;
	}
	template<int DIM>
	void applyFluxByDim(int fi, T const &lambda, std::array<std::array<std::array<T, N3>, P3>, NF> &dUdt) {
		static_assert(DIM < D);
		constexpr Range<int, D> iBox { repeat<D>(0), repeat<D>(N) };
		using index_type = MultiIndex<oBox, iBox>;
		constexpr int di = stride(DIM);
		auto const Fd = computeFlux<DIM>();
		for (auto I = index_type::begin(); I != index_type::end(); ++I) {
			int const ii = I;
			auto const fp = Fd[fi][ii + di];
			auto const fm = Fd[fi][ii];
			for (auto P = triindex_t::begin(); P != triindex_t::end(); ++P) {
				T const sgn = nonepow(P[DIM]);
				dUdt[fi][P][ii] -= lambda * (fp - sgn * fm);
			}
		}
	}
	template<int ... DIM>
	void applyFlux(int fi, T const &lambda, std::array<std::array<std::array<T, N3>, P3>, NF> &dUdt, std::integer_sequence<int, DIM...>) {
		(applyFluxByDim<DIM>(fi, lambda, dUdt), ...);
	}
	std::vector<std::array<std::array<T, N3>, P3>> computeSource() const {
		using quad_type = Quadrature<GaussLegendreQuadrature<T, P>, D>;
		using index_type = quad_type::index_type;
		quad_type const quadrature;
		constexpr int NH = quad_type::size();
		constexpr Basis<T, P, D> basis;
		auto const norm = basis.norm();
		std::vector<std::array<std::array<T, N3>, P3>> S(NF);
		for (auto I = index_type::begin(); I != index_type::end(); I++) {
			int const i0 = I;
			std::array<std::array<State, NH>, D> flux;
			for (int hi = 0; hi < NH; hi++) {
				auto const x = quadrature.point(hi);
				auto const phi = basis(x);
				State u0;
				for (int fi = 0; fi < NF; fi++) {
					u0[fi] = T(0);
					for (int pi = 0; pi < P3; pi++) {
						u0[fi] += norm[pi] * phi[pi] * U_[fi][pi][i0];
					}
				}
				for (int d = 0; d < D; d++) {
					flux[d][hi] = u0.flux(d);
				}
			}
			for (int fi = 0; fi < NF; fi++) {
				for (int pi = 0; pi < P3; pi++) {
					S[fi][pi][i0] = T(0);
				}
			}
			for (int hi = 0; hi < NH; hi++) {
				for (int fi = 0; fi < NF; fi++) {
					auto const x = quadrature.point(hi);
					auto const w = quadrature.weight(hi);
					for (int pi = 0; pi < P3; pi++) {
						for (int d = 0; d < D; d++) {
							auto const dphidx = basis.gradient(d, x);
							S[fi][pi][i0] += w * dphidx[pi] * flux[d][hi][fi];
						}
					}
				}
			}
		}
		return S;
	}
};
