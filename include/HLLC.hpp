/*
 * HLLC.hpp
 *
 *  Created on: Dec 17, 2024
 *      Author: dmarce1
 */

#ifndef INCLUDE_HLLC_HPP_
#define INCLUDE_HLLC_HPP_

#include <cmath>

#include "Constants.hpp"
#include "Options.hpp"
#include "Matrix.hpp"
#include "Vector.hpp"

namespace Hydrodynamics {

using namespace Math;

template<typename, int>
struct ConservedState;

template<typename, int>
struct PrimitiveState;

template<typename Type>
struct HydrodynamicsOptions {
	Type fluidGamma;
	Type dualEnergyPressureSwitch;
	Type dualEnergyUpdateSwitch;
	HydrodynamicsOptions() {
		auto const &opts = getOptions();
		fluidGamma = Type(opts.fluidGamma);
		dualEnergyPressureSwitch = Type(opts.dualEnergyPressureSwitch);
		dualEnergyUpdateSwitch = Type(opts.dualEnergyUpdateSwitch);
	}
};

template<typename Type>
Type energy2Entropy(Type rho, Type e) {
	constexpr Type one(1);
	static HydrodynamicsOptions<Type> const opts;
	static Type const gam = opts.fluidGamma;
	static Type const igamm1 = one / (gam - one);
	return rho * igamm1 * (log(e) - gam * log(rho));
}

template<typename Type>
Type entropy2Energy(Type rho, Type S) {
	constexpr Type one(1);
	static HydrodynamicsOptions<Type> const opts;
	static Type const gam = opts.fluidGamma;
	static Type const gamm1 = gam - one;
	return exp(gamm1 * S / rho) * pow(rho, gam);
}

template<typename Type, int Ndim>
struct ConservedState: public Vector<Type, Ndim + 3> {
	static constexpr int NFields = Ndim + 3;
	using base_type = Vector<Type, NFields>;
	Type &D;
	Type &E;
	Type &tau;
	Vector<Type, Ndim> &S;
	ConservedState() :
			D(base_type::operator[](0)), E(base_type::operator[](1)), tau(base_type::operator[](2)), S(
					(Vector<Type, Ndim>&) base_type::operator[](3)) {
	}
	ConservedState(ConservedState<Type, Ndim> const &U) :
			D(base_type::operator[](0)), E(base_type::operator[](1)), tau(base_type::operator[](2)), S(
					(Vector<Type, Ndim>&) base_type::operator[](3)) {
		*((base_type*) this) = (base_type const&) U;
	}
	ConservedState& operator=(PrimitiveState<Type, Ndim> const &V) {
		constexpr Type const half(0.5);
		D = V.rho;
		E = V.rho * (V.eps + half * vectorNorm(V.v));
		S = V.rho * V.v;
		tau = V.rho * V.s;
		return *this;
	}
	ConservedState& operator=(ConservedState<Type, Ndim> const &U) {
		*((base_type*) this) = (base_type const&) U;
		return *this;
	}
	ConservedState& operator=(base_type const &U) {
		(*(base_type*) this) = U;
		return *this;
	}
	void dualEnergyUpdate(Real maxEnergy) {
		constexpr Type const half(0.5);
		static HydrodynamicsOptions<Type> const opts;
		static Type const energySwitch = opts.dualEnergyUpdateSwitch;
		Type const eThermal = E - half * vectorNorm(S) / D;
		if (eThermal > energySwitch * maxEnergy) {
			tau = energy2Entropy(D, eThermal);
		}
	}
}
;

template<typename Type, int Ndim>
struct PrimitiveState: public Vector<Type, Ndim + 5> {
	static constexpr int NFields = Ndim + 5;
	using base_type = Vector<Type, NFields>;
	Type &rho;
	Type &eps;
	Type &s;
	Type &p;
	Type &c;
	Vector<Type, Ndim> &v;
	PrimitiveState() :
			rho(base_type::operator[](0)), eps(base_type::operator[](1)), s(base_type::operator[](2)), p(
					base_type::operator[](3)), c(base_type::operator[](4)), v(
					(Vector<Type, Ndim>&) base_type::operator[](5)) {
	}
	PrimitiveState& operator=(PrimitiveState<Type, Ndim> const &U) {
		*((base_type*) this) = (base_type const&) U;
		return *this;
	}
	PrimitiveState(ConservedState<Type, Ndim> const &U) :
			rho(base_type::operator[](0)), eps(base_type::operator[](1)), s(base_type::operator[](2)), p(
					base_type::operator[](3)), c(base_type::operator[](4)), v(
					(Vector<Type, Ndim>&) base_type::operator[](5)) {
		constexpr Type zero(0);
		constexpr Type half(0.5);
		constexpr Type one(1);
		static HydrodynamicsOptions<Type> const opts;
		static Type const gam = opts.fluidGamma;
		static Type const gamm1 = gam - one;
		static Type const energySwitch = opts.dualEnergyPressureSwitch;
		Type const iD = one / U.D;
		v = U.S * iD;
		rho = U.D;
		eps = U.E - half * vectorNorm(U.S) * iD;
		s = U.tau * iD;
		if (eps < energySwitch * U.E) {
			eps = entropy2Energy(U.D, U.tau);
		}
		p = gamm1 * eps;
		eps *= iD;
		c = sqrt(gam * gamm1 * eps);
	}
	ConservedState<Type, Ndim> toFlux(int k) const {
		static Type const half(0.5);
		ConservedState<Type, Ndim> F;
		F.D = rho * v[k];
		F.S = F.D * v;
		F.E = (p + rho * eps + half * rho * vectorNorm(v)) * v[k];
		F.tau = F.D * s;
		F.S[k] += p;
		return F;
	}
	Vector<Type, NFields> eigenValues(int k) const {
		static Type const half(0.5), one(1);
		Vector<Type, NFields> lambda;
		auto u = v;
		std::swap(u[0], u[k]);
		Type const h = eps + p / rho + half * vectorNorm(u);
		lambda[0] = u[0] - c;
		lambda[1] = u[0] + c;
		lambda[2] = u[0];
		lambda[3] = u[0];
		for (int l = 1; l < Ndim; l++) {
			lambda[3 + l] = u[l];
		}
		std::swap(lambda[3], lambda[3 + k]);
		return lambda;
	}
	SquareMatrix<Type, NFields> eigenVectors(int k) const {
		static Type const half(0.5), zero(0), one(1);
		Vector<Type, NFields> R;
		auto u = v;
		std::swap(u[0], u[k]);
		Type const h = gamma * eps + half * vectorNorm(u);
		R[0, 0] = one;
		R[1, 0] = h - u[0] * c;
		R[2, 0] = u[0] - c;
		R[3, 0] = zero;
		for (int l = 1; l < Ndim; l++) {
			R[3 + l, 0] = zero;
		}
		R[0, 1] = one;
		R[1, 1] = h + u[0] * c;
		R[2, 1] = u[0] + c;
		R[3, 1] = zero;
		for (int l = 1; l < Ndim; l++) {
			R[3 + l, 1] = u[l];
		}
		R[0, 2] = one;
		R[1, 2] = half * vectorNorm(u);
		R[2, 2] = u[0];
		R[3, 2] = zero;
		for (int l = 1; l < Ndim; l++) {
			R[3 + l, 2] = zero;
		}
		R[0, 3] = zero;
		R[1, 3] = zero;
		R[2, 3] = zero;
		R[3, 3] = one;
		for (int l = 1; l < Ndim; l++) {
			R[3 + l, 3] = zero;
		}
		for (int m = 1; m < Ndim; m++) {
			R[3 + m, 0] = zero;
			R[3 + m, 1] = u[m];
			R[3 + m, 2] = zero;
			R[3 + m, 3] = zero;
			for (int l = 1; l < Ndim; l++) {
				R[3 + m, 3 + l] = ((m == l) ? one : zero);
			}
		}
		for (int l = 0; l < NFields; l++) {
			std::swap(R[3, l], R[3 + k, l]);
		}
		return R;
	}
};

template<typename Type, int Ndim>
struct RiemannReturn {
	ConservedState<Type, Ndim> flux;
	Type signalSpeed;
};

template<typename Type, int Ndim>
RiemannReturn<Type, Ndim> riemannSolver(ConservedState<Type, Ndim> UL, ConservedState<Type, Ndim> UR, int k = 0) {
	static Type const zero(0), half(0.5), one(1), two(2);
	RiemannReturn<Type, Ndim> rc;
	ConservedState<Type, Ndim> U0;
	std::swap(UL.S[0], UL.S[k]);
	std::swap(UR.S[0], UR.S[k]);
	PrimitiveState<Type, Ndim> const VL(UL);
	PrimitiveState<Type, Ndim> const VR(UR);
	Type const sL = min(VL.v[0] - VL.c, VR.v[0] - VR.c);
	Type const sR = max(VL.v[0] + VL.c, VR.v[0] + VR.c);
	Type const wL = VL.rho * (sL - VL.v[0]);
	Type const wR = VR.rho * (sR - VR.v[0]);
	Type const s0 = (VR.p - VL.p + wL * VL.v[0] - wR * VR.v[0]) / (wL - wR);
	assert(sL <= s0);
	assert(s0 <= sR);
	if (s0 >= zero) {
		rc.flux = VL.toFlux(0);
		if (zero > sL) {
			U0.D = wL / (sL - s0);
			U0.E = U0.D * (UL.E / VL.rho + (s0 - VL.v[0]) * (s0 + VL.p / wL));
			U0.tau = UL.tau * (sL - VL.v[0]) / (sL - s0);
			U0.S[0] = U0.D * s0;
			for (int k = 1; k < Ndim; k++) {
				U0.S[k] = U0.D * VL.v[k];
			}
			rc.flux += sL * (U0 - UL);
		}
	} else {
		rc.flux = VR.toFlux(0);
		if (sR > zero) {
			U0.D = wR / (sR - s0);
			U0.E = U0.D * (UR.E / VR.rho + (s0 - VR.v[0]) * (s0 + VR.p / wR));
			U0.tau = UR.tau * (sR - VR.v[0]) / (sR - s0);
			U0.S[0] = U0.D * s0;
			for (int k = 1; k < Ndim; k++) {
				U0.S[k] = U0.D * VR.v[k];
			}
			rc.flux += sR * (U0 - UR);
		}
	}
	rc.signalSpeed = max(abs(sR), abs(sL));
	std::swap(rc.flux.S[0], rc.flux.S[k]);
	return rc;
}

}

#endif /* INCLUDE_HLLC_HPP_ */
