#include "MultiPrecision.hpp"
#include "DoubleReal.hpp"
#include "Integrate.hpp"
#include "FermiDirac.hpp"
#include "Quadrature.hpp"
#include "AutoDifferentiation.hpp"
#include "Polynomial.hpp"
#include "Constants.hpp"
#include <algorithm>
#include <memory>
#include <numeric>
#include <type_traits>
#include <vector>
#include <hpx/hpx_init.hpp>
#include "Octogrid.hpp"
#include "HyperSubgrid.hpp"
#include "dgTransforms.hpp"
#include "Options.hpp"
#include "MultiIndex.hpp"
#include "EulerState.hpp"
#include "Real.hpp"
#include "RungeKutta.hpp"
#include "BiconjugateGradient.hpp"
#include "Polynomial.hpp"
#include <cassert>
void testRadiation();

#include <tuple>
#include <utility>
#include <cstddef>
#include <unordered_map>
#include <ranges>
#include <array>
#include <cstddef>
#include <functional>
#include <map>
#include <array>
#include <cstddef>
#include <functional>

#include <array>
#include <cmath>
#include <algorithm>
#include <span>

template<typename Type, int D>
struct RadiationSource {
	RadiationSource() {
		Constants constants = getCgsConstants();
		mH = constants.u;
		kB = constants.kB;
		c = constants.c;
		aR = constants.aR;
	}
	using AutoType = Auto<Type, D + 1>;
	void solve(Type &Er, std::array<Type, D> &F, Type &Eg, std::array<Type, D> &S, Type const &κ, Type const &χ, Type const &ρ, Type const &γ, Type const &μ,
			Type const &dt) {
		std::array<AutoType, D + 1> R, M, dR;
		for (int d = 0; d < D; d++) {
			R[d] = F[d];
			M[d] = S[d];
		}
		R[D] = Er;
		M[D] = Eg;
		std::function<std::array<AutoType, D + 1>(std::array<AutoType, D + 1> const&)> f = [this, R, M, κ, χ, ρ, γ, μ, dt](
				std::array<AutoType, D + 1> const &dR) {
			return residual(dR, R, M, κ, χ, ρ, γ, μ, dt);
		};
		dR.fill(0);
		return newtonRhapson(f, dR);
	}
	auto residual(std::array<AutoType, D + 1> dR, std::array<AutoType, D + 1> R0, std::array<AutoType, D + 1> M0, Type const &κ, Type const &χ, Type const &ρ,
			Type const &γ, Type const &μ, Type const &dt) {
		std::array<AutoType, D + 1> R, M, G, residual;
		for (int d = 0; d <= D; d++) {
			R[d] = R0[d] + dR[d];
			M[d] = M0[d] - dR[d];
		}
		G = source(R, M, κ, χ, ρ, γ, μ);
		for (int d = 0; d <= D; d++) {
			residual[d] = dR[d] + dt * G[d];
		}
		return residual;
	}
	auto source(std::array<AutoType, D + 1> Rlab, std::array<AutoType, D + 1> M, Type const &κ, Type const &χ, Type const &ρ, Type const &γ, Type const &μ) {
		constexpr Type ε = 1e-20;
		constexpr Type π = M_PI;
		auto const Er = Rlab[D];
		auto const Eg = M[D];
		auto const F = std::span(Rlab.cbegin(), Rlab.cend() - 1);
		auto const S = std::span(M.cbegin(), M.cend() - 1);
		std::array<std::array<AutoType, D + 1>, D + 1> Λ, P, Rcom;
		std::array<AutoType, D + 1> Gcom, Glab;
		std::array<AutoType, D> β, n;
		AutoType W, β2, F2, magF, f, f2, ξ, Dd, Ds, c2, Ek, Ei, p, T, T4;
		c2 = sqr(c);
		F2 = β2 = T(0);
		for (int i = 0; i < D; i++) {
			β[i] = S[i] / (ρ * c);
			β2 += sqr(β[i]);
			F2 += sqr(F[i]);
		}
		Ek = 0.5 * ρ * c2 * β2;
		Ei = Eg - Ek;
		p = (γ - 1) * Ei;
		T = μ * mH * p / (kB * ρ);
		W = 1 / sqrt(1 - β2);
		Λ[D][D] = W;
		for (int j = 0; j < D; j++) {
			Λ[D][j] = Λ[j][D] = -W * β[j];
			for (int k = 0; k < j; k++) {
				Λ[j][k] = Λ[k][j] = (W - 1) * β[j] * β[k] / (β2 + ε);
			}
			Λ[j][j] = 1 + (W - 1) * sqr(β[j]) / (β2 + ε);
		}
		magF = sqrt(F2);
		for (int d = 0; d < D; d++) {
			n[d] = F[d] / (magF + ε);
		}
		f = magF / (Er + ε);
		f2 = f * f;
		ξ = (3 + 4 * f2) / (5 + 2 * sqrt(4 - 3 * f2));
		Dd = 0.5 * (1 - ξ);
		Ds = 0.5 * (3 * ξ - 1);
		P[D][D] = Er;
		for (int j = 0; j < D; j++) {
			P[D][j] = P[j][D] = F[j];
			for (int k = 0; k < j; k++) {
				P[j][k] = P[k][j] = Ds * Er * n[j] * n[k];
			}
			P[j][j] = (Dd + Ds * sqr(n[j])) * Er;
		}
		for (int i = 0; i <= D; i++) {
			Rcom[i] = 0;
			for (int j = 0; j <= D; j++) {
				for (int k = 0; k <= D; k++) {
					Rcom[i] += Λ[i][j] * P[j][k] * Λ[k][D];
				}
			}
		}
		Gcom[D] = c * κ * (Rcom[D] - aR * sqr(sqr(T)));
		for (int j = 0; j < D; j++) {
			Gcom[j] = c * χ * Rcom[j];
		}
		for (int i = 0; i <= D; i++) {
			Glab[i] = 0;
			for (int j = 0; j <= D; j++) {
				Glab[i] -= Λ[i][j] * Gcom[j];
			}
		}
		return Glab;
	}
private:
	Type mH;
	Type kB;
	Type c;
	Type aR;
};

constexpr int P = 3;
constexpr int D = 2;
constexpr int N = 64;
using T = double;
using RK = RungeKutta<T, P>::type;

using SubgridType = HyperSubgrid<T, D, N, P, RK, EulerStateHLLC>;
using OctogridServerType = OctogridServer<SubgridType>;
using OctogridType = hpx::components::component<OctogridServerType>;

HPX_REGISTER_COMPONENT_MODULE();
HPX_REGISTER_COMPONENT (OctogridType);

using GetFaceChildrenAction = typename OctogridServerType::refineAction;
using RefineAction = typename OctogridServerType::getFaceChildrenAction;
using SetAction = typename OctogridServerType::setAction;

HPX_REGISTER_ACTION(GetFaceChildrenAction, getFaceChildrenAction);
HPX_REGISTER_ACTION(RefineAction, refineAction);
HPX_REGISTER_ACTION(SetAction, setAction);

struct FDTest {
	double k;
	double η;
	double β;
	double FD;
};



void polytest();

int hpx_main(int argc, char *argv[]) {
	polytest();
//	for (Real ρ = 1e10; ρ < Real(1e20); ρ *= Real(10)) {
//		printf( " ρ          T          n          ε          p\n");
//		for (Real T = 10; T < Real(1e13); T *= Real(10)) {
//			printf("%e %e %e %e %e\n", (double) eos.ρ, (double) eos.T, (double) eos.n, (double) eos.e, (double) eos.p);
//		}
//		printf("\n");
//	}

//	RadiationSource<double, 3> test;
//	auto tmp = bellPolynomial<6>(6, 3);
//	std::cout << tmp << "\n";

//	std::cout << "/**************************************/\n";
//	fluxExpression<3>(0, { 0, 0, 0, 0, 0 });
//	std::cout << "/**************************************/\n";
//	fluxExpression<3>(0, { 2, 0, 0, 0, 0 });
//	std::cout << "/**************************************/\n";
//	fluxExpression<3>(0, { 0, 2, 0, 0, 0 });
//	std::cout << "/**************************************/\n";
//	fluxExpression<3>(0, { 0, 0, 2, 0, 0 });
//	std::cout << "/**************************************/\n";
//	fluxExpression<3>(0, { 0, 0, 0, 2, 0 });
//	std::cout << "/**************************************/\n";
//	fluxExpression<3>(0, { 0, 0, 0, 0, 2 });
//	std::cout << "/**************************************/\n";
//	constexpr int N = 3;
//	Auto<double> x(10.0, 1.0);
//	Auto<double> dfxdx(10.0, 1.0);
//	auto f = [](auto x) {
//		return x * x * x;
//	};
//	auto g = derivative2(f, 10.0);
//	printf( "%e \n", derivative3(f, 10.0).derivative().derivative().value().value());

//	SymbolicVariable<0> rho("rho");
//	SymbolicVariable<1> S("S");
//	SymbolicVariable<2> E("E");
//	SymbolicConstant gamma(5./3.);
//	SymbolicOne one;
//	SymbolicConstant half(0.5);Constant
//	auto const p = (gamma - one) * (E - half * S * S / rho);
//	auto const U = std::make_tuple(S, S * S / rho + p, (E + p) * S / rho);
//	auto J = jacobian(U);
//	std::cout << jacobianToString(J) << "\n";
//	std::cout << "S = " << static_cast<std::string>(S) << "\n";
//	processOptions(argc, argv);
//	printf("\nPrologue complete\n");
//	RK const rk;
//	SubgridType grid;
//	grid.initialize(initSodShockTube<T, D>);
//	grid.applyLimiter();
//	grid.output("X", 0, Real(0.0));
//	T t = T(0);
//	T tmax = T(.125);
//	T dt;
//	int iter = 0;
//	while (t < tmax) {
//		grid.output("X", iter, t);
//		std::cout << "i = " << std::to_string(iter);
//		std::cout << "  t = " << std::to_string(t);
//		dt = grid.beginStep();
//		std::cout << "  dt = " << dt << std::endl;
//		for (int s = 0; s < rk.stageCount(); s++) {
//			grid.subStep(dt, s);
//		}
//		grid.endStep();
//		iter++;
//		t += dt;
//	}
	printf("\nStopping\n");
	return hpx::local::finalize();
}

thread_local FpeThreadInit _fpeThreadInitGuard;

int main(int argc, char *argv[]) {
	installFpeHandler();
	_fpeThreadInitGuard.touch();
	std::vector<std::string> cfg = { "hpx.commandline.allow_unknown=1" };
	cfg.push_back("hpx.stacks.small_size=1048576");
	hpx::init_params init_params;
	init_params.cfg = std::move(cfg);
	auto rc = hpx::init(argc, argv, init_params);
	return rc;
}
