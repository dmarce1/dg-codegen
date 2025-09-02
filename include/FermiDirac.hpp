#pragma once

#include "Integrate.hpp"
#include "Hypergeometric.hpp"
#include "FermiDirac.hpp"
#include "HalfInteger.hpp"
#include "Util.hpp"

#include <limits>
#include <numeric>

// Α,α,Β,β,Γ,γ,Δ,δ,Ε,ε,Ζ,ζ,Η,η,Θ,θ,ϑ,Ι,ι,Κ,κ,ϰ,Λ,λ,Μ,μ,Ν,ν,Ξ,ξ,Ο,ο,Π,π,ϖ,Ρ,ρ,ϱ,Σ,σ,ς,Τ,τ,Υ,υ,Φ,φ,ϕ,Χ,χ,Ψ,ψ,Ω,ω

#define List( ... ) {__VA_ARGS__}

template<typename Real = double>
struct ElectronEoS {
	ElectronEoS(Real _ρ, Real _T, Real _μe = Real(2)) :
			ρ(_ρ), T(_T), μe(_μe) {
		using std::atan;
		using std::exp;
		using std::log;
		using std::sqrt;
		using std::numeric_limits;
		static Real const zero = Real(0);
		static Real const one = Real(1);
		static Real const two = Real(2);
		static Real const three = Real(3);
		static Real const four = Real(4);
		static Real const five = Real(5);
		static Real const eight = Real(8);
		static Real const half = one / two;
		static Real const threeHalves = three / two;
		static Real const fiveHalves = five / two;
		static Real const infinity = std::numeric_limits<Real>::infinity();
//		static Real const Ξ = std::numeric_limits<Real>::min();
//		static Real const ε = std::numeric_limits<Real>::epsilon();
		static Real const c = 2.997924580000000e+10;
		static Real const h = 6.626070150000000e-27;
		static Real const kB = 3.806490000000000e-16;
		static Real const me = 9.109383701500000e-28;
		static Real const mH = 1.67353271600000e-24;
		static Real const π = four * atan(one);
		static Real const N0 = eight * sqrt(two) * π * ipow(me * c / h, 3);
		static Real const P0 = N0 * sqrt(two) / three * me * sqr(c);
		static Real const E0 = threeHalves * P0;
		static auto const FD = [](Real k, Real η, Real β) {
			std::function<Real(Real const&)> const f = [k, η, β](Real const &x) {
				static Real const xMax = log(numeric_limits<Real>::max());
				if (x < xMax + η) {
					Real const num = pow(x, k) * sqrt(one + half * β * x);
					Real const den = one + exp(x - η);
					return num / den;
				} else {
					return zero;
				}
			};
			return integrate(f, zero, infinity);
		};
//		static auto const dFD_dη = [](Real k, Real η, Real β) {
//			std::function<Real(Real const&)> const f = [k, η, β](Real const &x) {
//				Real constexpr xMax = 0.5 * log(sqrt(std::numeric_limits<Real>::max()) - 1);
//				if (x < xMax + η) {
//					Real const num = pow(x, k) * sqrt(one + half * β * x) * exp(x - η);
//					Real const den = sqr(one + exp(x - η));
//					return num / den;
//				} else {
//					return zero;
//				}
//			};
//			return integrate(f, zero, infinity);
//		};
		Real const β = (kB * T) / (me * sqr(c));
		Real const β32 = pow(β, threeHalves);
		Real const β52 = pow(β, fiveHalves);
		auto fne = [β, β32](Real η) {
			return N0 * β32 * (FD(half, η, β) + FD(threeHalves, η, β));
		};
		auto fnp = [β, β32](Real η) {
			return N0 * β32 * (FD(half, -η - two / β, β) + β * FD(threeHalves, -η - two / β, β));
		};
		auto F = [this, fne, fnp, β, β32](Real η) {
			return ρ / (μe * mH) - fne(η) + fnp(η);
		};
//		auto dFdη = [β, β32](Real η) {
//			Real const dne_dη = +N0 * β32 * (dFD_dη(half, η, β) + dFD_dη(threeHalves, η, β));
//			Real const dnp_dη = -N0 * β32 * (dFD_dη(half, -η - two / β, β) - β * dFD_dη(threeHalves, -η - two / β, β));
//			return dnp_dη - dne_dη;
//		};
		Real ηmax = one;
		while (F(+ηmax) * F(-ηmax) >= zero) {
			printf("1. %e\n", (double) ηmax);
			ηmax *= two;
		}
		Real ηmin = -ηmax;
		Real η;
		do {
			η = half * (ηmin + ηmax);
			if (F(η) * F(ηmin) >= zero) {
				ηmin = η;
			} else {
				ηmax = η;
			}
			printf("4. %e\n", (double) (ηmax - ηmin));
		} while (ηmax > nexttoward(ηmin, ηmax));
//		do {
//			Real δη = -F(η) / dFdη(η);
//			η += δη;
//			δ = abs(δη / (η + Ξ));
//		} while (δ > ε);
		Real const ne = fne(η);
		Real const np = fnp(η);
		Real const F15e = FD(threeHalves, η, β);
		Real const F25e = FD(fiveHalves, η, β);
		Real const F15p = FD(threeHalves, -η - two / β, β);
		Real const F25p = FD(fiveHalves, -η - two / β, β);
		Real const ee = E0 * β52 * (F15e + β * F25e);
		Real const ep = E0 * β52 * (F15p + β * F25p);
		Real const pe = P0 * β52 * (F15e + half * β * F25e);
		Real const pp = P0 * β52 * (F15p + half * β * F25p);
		n = ne + np;
		p = pe + pp;
		e = ee + ep;
	}
	Real ρ;
	Real T;
	Real μe;
	Real η;
	Real n;
	Real e;
	Real p;
};

