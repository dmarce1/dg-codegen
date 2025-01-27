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
#include "Utilities.hpp"
#include "Vector.hpp"

namespace Hydrodynamics {

using namespace Math;

static constexpr int scalarFieldCount = 3;

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

template<typename >
struct EquationOfState;

template<typename T, int N, typename = Vector<T, N + scalarFieldCount>>
struct ConservedState;

template<typename T, int N, typename = Vector<T, N + scalarFieldCount>>
struct PrimitiveState;

template<typename, int, typename >
struct RiemannReturn;

template<typename Type, int Ndim, typename Container>
RiemannReturn<Type, Ndim, Container> riemannSolver(ConservedState<Type, Ndim, Container> UL, ConservedState<Type, Ndim, Container> UR, int k = 0);

template<typename Type>
Type energy2Entropy(Type rho, Type e);

template<typename Type>
Type entropy2Energy(Type rho, Type S);

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
		Type const c2 = dPdRho + P * dPdEps / nSquared(rho);
		return sqrt(c2);
	}
};

template<typename Type, int Ndim, typename Container>
struct ConservedState: public Container {
	static constexpr int NFields = Ndim + scalarFieldCount;
	ConservedState() :
			createContainer(*this), D(Container::operator[](0)), E(Container::operator[](1)), tau(Container::operator[](2)), S(
					(Vector<Type, Ndim>&) Container::operator[](3)) {
	}
	ConservedState(ConservedState const &U) :
			createContainer(*this), D(Container::operator[](0)), E(Container::operator[](1)), tau(Container::operator[](2)), S(
					(Vector<Type, Ndim>&) Container::operator[](3)) {
		*((Container*) this) = (Container const&) U;
	}
	ConservedState& operator=(PrimitiveState<Type, Ndim, Container> const &V) {
		static constexpr Type half = Type(0.5);
		D = V.rho;
		E = V.rho * (V.eps + half * vectorNorm(V.v));
		S = V.rho * V.v;
		tau = V.rho * V.s;
		return *this;
	}
	ConservedState& operator=(ConservedState const &U) {
		*((Container*) this) = (Container const&) U;
		return *this;
	}
	ConservedState& operator=(Container const &U) {
		(*(Container*) this) = U;
		return *this;
	}
	SquareMatrix<Type, NFields> eigenVectors(int k) const {
		return PrimitiveState<Type, Ndim, Container>().eigenVectors(k);
	}
	void dualEnergyUpdate(Real maxEnergy) {
		static constexpr Type half = Type(0.5);
		Type const eThermal = E - half * vectorNorm(S) / D;
		if (eThermal > hydroOptions.dualEnergyUpdateSwitch * maxEnergy) {
			tau = energy2Entropy(D, eThermal);
		}
	}
	static HydrodynamicsOptions<Type> const hydroOptions;
	ContainerResizer<Container, NFields> const createContainer;
	Type &D;
	Type &E;
	Type &tau;
	Vector<Type, Ndim> &S;
};

template<typename Type, int Ndim, typename Container>
struct PrimitiveState: public Container {
	static constexpr int NFields = Ndim + scalarFieldCount;
	PrimitiveState() :
			createContainer(*this), rho(Container::operator[](0)), eps(Container::operator[](1)), s(Container::operator[](2)), v(
					(Vector<Type, Ndim>&) Container::operator[](3)) {
	}
	PrimitiveState(ConservedState<Type, Ndim, Container> const &U) :
			createContainer(*this), rho(Container::operator[](0)), eps(Container::operator[](1)), s(Container::operator[](2)), v(
					(Vector<Type, Ndim>&) Container::operator[](3)) {
		static constexpr Type half = Type(0.5), one = Type(1);
		Type const iD = one / U.D;
		v = U.S * iD;
		rho = U.D;
		eps = U.E - half * vectorNorm(U.S) * iD;
		s = U.tau * iD;
		if (eps < hydroOptions.dualEnergyPressureSwitch * U.E) {
			eps = entropy2Energy(U.D, U.tau);
		}
		eps *= iD;
	}
	PrimitiveState& operator=(PrimitiveState const &U) {
		*((Container*) this) = (Container const&) U;
		return *this;
	}
	ConservedState<Type, Ndim, Container> toFlux(int k) const {
		static constexpr Type half = Type(0.5);
		ConservedState<Type, Ndim, Container> F;
		EquationOfState eos(rho, eps);
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
		static constexpr Type zero = Type(0), half = Type(0.5), one = Type(1);
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
	static HydrodynamicsOptions<Type> const hydroOptions;
	ContainerResizer<Container, NFields> const createContainer;
	Type &rho;
	Type &eps;
	Type &s;
	Vector<Type, Ndim> &v;
};

template<typename Type, int Ndim, typename Container>
struct RiemannReturn {
	ConservedState<Type, Ndim, Container> flux;
	Type signalSpeed;
};

template<typename Type, int Ndim, typename Container>
RiemannReturn<Type, Ndim, Container> riemannSolver(ConservedState<Type, Ndim, Container> UL, ConservedState<Type, Ndim, Container> UR, int k) {
	static constexpr Type zero(0);
	std::swap(UL.S[0], UL.S[k]);
	std::swap(UR.S[0], UR.S[k]);
	RiemannReturn<Type, Ndim, Container> rc;
	ConservedState<Type, Ndim, Container> U0;
	PrimitiveState<Type, Ndim, Container> const VL(UL);
	PrimitiveState<Type, Ndim, Container> const VR(UR);
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

template<typename Type, int Ndim, typename Container>
HydrodynamicsOptions<Type> const ConservedState<Type, Ndim, Container>::hydroOptions { };

template<typename Type, int Ndim, typename Container>
HydrodynamicsOptions<Type> const PrimitiveState<Type, Ndim, Container>::hydroOptions { };

}

#endif /* INCLUDE_HYDRODYNAMICS_HPP_ */
