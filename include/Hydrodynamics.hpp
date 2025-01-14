/*
 * HLLC.hpp
 *
 *  Created on: Dec 17, 2024
 *      Author: dmarce1
 */

#ifndef INCLUDE_HYDRODYNAMICS_HPP_
#define INCLUDE_HYDRODYNAMICS_HPP_

#include <cmath>

#include "Constants.hpp"
#include "Options.hpp"
#include "Matrix.hpp"
#include "Vector.hpp"

namespace Hydrodynamics {

using namespace Math;

template<typename Type>
struct HydrodynamicsOptions;

template<typename Type>
struct EquationOfState;

template<typename, int>
struct ConservedState;

template<typename, int>
struct PrimitiveState;

template<typename Type, int Ndim>
struct RiemannReturn;

template<typename Type, int Ndim>
RiemannReturn<Type, Ndim> riemannSolver(ConservedState<Type, Ndim> UL, ConservedState<Type, Ndim> UR, int k = 0);

template<typename Type>
Type energy2Entropy(Type rho, Type e);

template<typename Type>
Type entropy2Energy(Type rho, Type S);

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
		static constexpr Type half(0.5);
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
		static constexpr Type half(0.5);
		static HydrodynamicsOptions<Type> const opts;
		static Type const energySwitch = opts.dualEnergyUpdateSwitch;
		Type const eThermal = E - half * vectorNorm(S) / D;
		if (eThermal > energySwitch * maxEnergy) {
			tau = energy2Entropy(D, eThermal);
		}
	}
};


template<typename Type>
struct EquationOfState {
	Type P;
	Type dPdRho;
	Type dPdEps;
	EquationOfState(Type rho, Type eps) {
		static HydrodynamicsOptions<Type> const opts;
		static constexpr Type one(1);
		static Type const gam = opts.fluidGamma;
		static Type const gamm1 = gam - one;
		P = gamm1 * rho * eps;
		dPdRho = gamm1 * eps;
		dPdEps = gamm1 * rho;
	}
	Type soundSpeed(Type rho, Type eps) const {
		Type const c2 = dPdRho + P * dPdEps / sqr(rho);
		return sqrt(c2);
	}
};

template<typename Type, int Ndim>
struct PrimitiveState: public Vector<Type, Ndim + 3> {
	static constexpr int NFields = Ndim + 3;
	using base_type = Vector<Type, NFields>;
	Type &rho;
	Type &eps;
	Type &s;
	Vector<Type, Ndim> &v;
	PrimitiveState() :
			rho(base_type::operator[](0)), eps(base_type::operator[](1)), s(base_type::operator[](2)), v(
					(Vector<Type, Ndim>&) base_type::operator[](3)) {
	}
	PrimitiveState& operator=(PrimitiveState<Type, Ndim> const &U) {
		*((base_type*) this) = (base_type const&) U;
		return *this;
	}
	PrimitiveState(ConservedState<Type, Ndim> const &U) :
			rho(base_type::operator[](0)), eps(base_type::operator[](1)), s(base_type::operator[](2)), v(
					(Vector<Type, Ndim>&) base_type::operator[](3)) {
		static constexpr Type one(1);
		static constexpr Type half(0.5);
		static HydrodynamicsOptions<Type> const opts;
		static Type const energySwitch = opts.dualEnergyPressureSwitch;
		Type const iD = one / U.D;
		v = U.S * iD;
		rho = U.D;
		eps = U.E - half * vectorNorm(U.S) * iD;
		s = U.tau * iD;
		if (eps < energySwitch * U.E) {
			eps = entropy2Energy(U.D, U.tau);
		}
		eps *= iD;
	}
	ConservedState<Type, Ndim> toFlux(int k) const {
		static Type const half(0.5);
		EquationOfState eos(rho, eps);
		ConservedState<Type, Ndim> F;
		F.D = v[k] * rho;
		F.S = v[k] * rho * v;
		F.S[k] += eos.P;
		F.E = v[k] * (eos.P + rho * (eps + half * vectorNorm(v)));
		F.tau = v[k] * rho * s;
		return F;
	}
	Vector<Type, NFields> eigenValues(int k) const {
		Vector<Type, NFields> lambda;
		EquationOfState eos(rho, eps);
		Type const a = eos.soundSpeed(rho, eps);
		Type const u = v[k];
		lambda[0] = u - a;
		lambda[1] = u + a;
		lambda[2] = u;
		for (int n = 0; n < Ndim; n++) {
			lambda[3 + n] = u;
		}
		return lambda;
	}
	SquareMatrix<Type, NFields> eigenVectors(int k) const {
		static constexpr Type zero(0);
		static constexpr Type half(0.5);
		static constexpr Type one(1);
		EquationOfState<Type> eos(rho, eps);
		Vector<Type, NFields> R;
		Type const h = gamma * eps + half * vectorNorm(v);
		Type const a = eos.soundSpeed(rho, eps);
		Type const u = v[k];
		Type const psi = one - (rho * eos.dPdRho) / (eps * eos.dPdEps);
		R[0, 0] = one;
		R[1, 0] = h - u * a;
		R[2, 0] = s;
		for (int n = 0; n < Ndim; n++) {
			R[3 + n, 0] = u[n];
		}
		R[3 + k, 0] -= a;
		R[0, 1] = one;
		R[1, 1] = h + u * a;
		R[2, 1] = s;
		for (int n = 0; n < Ndim; n++) {
			R[3 + n, 1] = u[n];
		}
		R[3 + k, 1] += a;
		R[0, 2] = zero;
		R[1, 2] = zero;
		R[2, 2] = one;
		for (int n = 0; n < Ndim; n++) {
			R[3 + n, 2] = zero;
		}
		R[0, 3 + k] = one;
		R[1, 3 + k] = half * vectorNorm(u) + psi * eps;
		R[2, 3 + k] = zero;
		for (int n = 0; n < Ndim; n++) {
			R[3 + n, 3 + k] = u[n];
		}
		for (int m = 0; m < Ndim; m++) {
			if (m != k) {
				R[0, 3 + m] = zero;
				R[1, 3 + m] = u[m];
				R[2, 3 + m] = zero;
				for (int n = 0; n < Ndim; n++) {
					R[3 + n, 3 + m] = ((m == n) ? one : zero);
				}
			}
		}
		return R;
	}
};

template<typename Type, int Ndim>
struct RiemannReturn {
	ConservedState<Type, Ndim> flux;
	Type signalSpeed;
};

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

template<typename Type, int Ndim>
RiemannReturn<Type, Ndim> riemannSolver(ConservedState<Type, Ndim> UL, ConservedState<Type, Ndim> UR, int k) {
	static constexpr Type zero(0);
	std::swap(UL.S[0], UL.S[k]);
	std::swap(UR.S[0], UR.S[k]);
	RiemannReturn<Type, Ndim> rc;
	ConservedState<Type, Ndim> U0;
	PrimitiveState<Type, Ndim> const VL(UL);
	PrimitiveState<Type, Ndim> const VR(UR);
	EquationOfState<Type> eosL(VL.rho, VL.eps);
	EquationOfState<Type> eosR(VR.rho, VR.eps);
	Type const pL = eosL.P;
	Type const pR = eosR.P;
	Type const aL = eosL.soundSpeed(VL.rho, VL.eps);
	Type const aR = eosR.soundSpeed(VL.rho, VL.eps);
	Type const sL = min(VL.v[0] - aL, VR.v[0] - aR);
	Type const sR = max(VL.v[0] + aL, VR.v[0] + aR);
	Type const wL = VL.rho * (sL - VL.v[0]);
	Type const wR = VR.rho * (sR - VR.v[0]);
	Type const s0 = (pR - pL + wL * VL.v[0] - wR * VR.v[0]) / (wL - wR);
	assert(sL <= s0);
	assert(s0 <= sR);
	if (s0 >= zero) {
		rc.flux = VL.toFlux(0);
		if (zero > sL) {
			U0.D = wL / (sL - s0);
			U0.E = U0.D * (UL.E / VL.rho + (s0 - VL.v[0]) * (s0 + pL / wL));
			U0.tau = U0.D * VL.s;
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
			U0.E = U0.D * (UR.E / VR.rho + (s0 - VR.v[0]) * (s0 + pR / wR));
			U0.tau = U0.D * VR.s;
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

template<typename Type>
Type energy2Entropy(Type rho, Type e) {
	static HydrodynamicsOptions<Type> const opts;
	static constexpr Type one(1);
	static Type const gam = opts.fluidGamma;
	static Type const igamm1 = one / (gam - one);
	return rho * igamm1 * (log(e) - gam * log(rho));
}

template<typename Type>
Type entropy2Energy(Type rho, Type S) {
	static HydrodynamicsOptions<Type> const opts;
	static constexpr Type one(1);
	static Type const gam = opts.fluidGamma;
	static Type const gamm1 = gam - one;
	return exp(gamm1 * S / rho) * pow(rho, gam);
}

}

#endif /* INCLUDE_HYDRODYNAMICS_HPP_ */
