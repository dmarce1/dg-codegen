/*
 * HLLC.hpp
 *
 *  Created on: Dec 17, 2024
 *      Author: dmarce1
 */

#ifndef INCLUDE_HLLC_HPP_
#define INCLUDE_HLLC_HPP_

#include <cmath>

#include "Matrix.hpp"
#include "Vector.hpp"

namespace Hydrodynamics {

using namespace Math;

template<typename, double>
struct EquationOfState;

template<typename, int, typename >
struct ConservedState;

template<typename, int, typename >
struct PrimitiveState;

template<typename Type, double Gamma>
struct EquationOfState {
	Type pressure(Type rho, Type eps) const {
		return Type(Gamma - 1.0) * rho * eps;
	}
	Type soundSpeed(Type rho, Type eps) const {
		return sqrt(Type(Gamma * (Gamma - 1.0)) * eps);
	}
};

template<typename Type, int Ndim, typename EoS>
struct ConservedState: public Vector<Type, Ndim + 2> {
	static constexpr int NFields = Ndim + 2;
	using base_type = Vector<Type, NFields>;
	Type &D;
	Type &Tau;
	Vector<Type, Ndim> &S;
	ConservedState() :
			D(base_type::operator[](0)), Tau(base_type::operator[](1)), S(
					(Vector<Type, Ndim>&) base_type::operator[](2)) {
	}
	ConservedState(ConservedState<Type, Ndim, EoS> const &U) :
			D(base_type::operator[](0)), Tau(base_type::operator[](1)), S(
					(Vector<Type, Ndim>&) base_type::operator[](2)) {
		*((base_type*) this) = (base_type const&) U;
	}
	ConservedState& operator=(PrimitiveState<Type, Ndim, EoS> const &V) {
		static Type const half(0.5);
		D = V.rho;
		Tau = V.rho * (V.eps + half * vectorNorm(V.v));
		S = V.rho * V.v;
		return *this;
	}
	ConservedState& operator=(ConservedState<Type, Ndim, EoS> const &U) {
		*((base_type*) this) = (base_type const&) U;
		return *this;
	}
	ConservedState& operator=(base_type const &U) {
		(*(base_type*) this) = U;
		return *this;
	}
};

template<typename Type, int Ndim, typename Eos>
struct PrimitiveState: public Vector<Type, Ndim + 4> {
	static constexpr int NFields = Ndim + 4;
	using base_type = Vector<Type, NFields>;
	Type &rho;
	Type &eps;
	Type &p;
	Type &c;
	Vector<Type, Ndim> &v;
	PrimitiveState() :
			rho(base_type::operator[](0)), eps(base_type::operator[](1)), p(base_type::operator[](2)), c(
					base_type::operator[](3)), v((Vector<Type, Ndim>&) base_type::operator[](4)) {
	}
	PrimitiveState& operator=(PrimitiveState<Type, Ndim, Eos> const &U) {
		*((base_type*) this) = (base_type const&) U;
		return *this;
	}
	PrimitiveState(ConservedState<Type, Ndim, Eos> const &U) :
			rho(base_type::operator[](0)), eps(base_type::operator[](1)), p(base_type::operator[](2)), c(
					base_type::operator[](3)), v((Vector<Type, Ndim>&) base_type::operator[](4)) {
		Eos eos;
		static Type const zero(0);
		static Type const half(0.5);
		static Type const one(1);
		Type const iDensity = one / U.D;
		v = U.S * iDensity;
		rho = U.D;
		eps = max((U.Tau - half * vectorNorm(U.S) * iDensity), zero) * iDensity;
		p = eos.pressure(rho, eps);
		c = eos.soundSpeed(rho, eps);
	}
	ConservedState<Type, Ndim, Eos> toFlux(int k) const {
		static Type const half(0.5);
		ConservedState<Type, Ndim, Eos> F;
		F.D = rho * v[k];
		F.S = rho * v * v[k];
		F.Tau = (p + rho * eps + half * rho * vectorNorm(v)) * v[k];
		F.S[k] += p;
		return F;
	}
/*	Vector<Type, NFields> eigenValues(int k) const {
		static Type const half(0.5), one(1);
		Vector<Type, NFields> lambda;
		auto u = v;
		std::swap(u[0], u[k]);
		Type const h = eps + p / rho + half * vectorNorm(u);
		lambda[0] = u[0] - c;
		lambda[1] = u[0] + c;
		lambda[2] = u[0];
		for (int l = 1; l < Ndim; l++) {
			lambda[2 + l] = u[l];
		}
		std::swap(lambda[2], lambda[2 + k]);
		return lambda;
	}
	SquareMatrix<Type, NFields> eigenVectors(int k) const {
		static Type const half(0.5), zero(0), one(1);
		Vector<Type, NFields> R;
		auto u = v;
		std::swap(u[0], u[k]);
		Type const h = gamma * eps + half * vectorNorm(u);
		R[0, 0] = one;
		R[0, 1] = one;
		R[0, 2] = one;
		for (int l = 1; l < Ndim; l++) {
			R[0, 2 + l] = zero;
		}
		R[1, 0] = h + u[0] * c;
		R[1, 1] = h - u[0] * c;
		R[1, 2] = half * vectorNorm(u);
		for (int l = 1; l < Ndim; l++) {
			R[1, 2 + l] = u[l];
		}
		R[2, 0] = u[0] - c;
		R[2, 1] = u[0] + c;
		R[2, 2] = u[0];
		for (int l = 1; l < Ndim; l++) {
			R[2, 2 + l] = zero;
		}
		for (int l = 1; l < Ndim; l++) {
			R[2 + l, 0] = R[2 + l, 1] = R[21 + l, 2] = u[l];
			for (int m = 1; m < Ndim; l++) {
				R[2 + l, 2 + m] = ((m == l) ? one : zero);
			}
		}
		for (int l = 0; l < NFields; l++) {
			std::swap(R[l, 2], R[l, 2 + k]);
		}
		for (int l = 0; l < NFields; l++) {
			std::swap(R[2, l], R[2 + k, l]);
		}
		return R;
	}*/
};

template<typename Type, int Ndim, typename Eos>
struct RiemannReturn {
	ConservedState<Type, Ndim, Eos> flux;
	Type signalSpeed;
};

template<typename Type, int Ndim, typename Eos>
RiemannReturn<Type, Ndim, Eos> riemannSolver(ConservedState<Type, Ndim, Eos> UL, ConservedState<Type, Ndim, Eos> UR,
		int k = 0) {
	static Type const zero(0), half(0.5), one(1), two(2);
	RiemannReturn<Type, Ndim, Eos> rc;
	ConservedState<Type, Ndim, Eos> U0;
	std::swap(UL.S[0], UL.S[k]);
	std::swap(UR.S[0], UR.S[k]);
	PrimitiveState<Type, Ndim, Eos> const VL(UL);
	PrimitiveState<Type, Ndim, Eos> const VR(UR);
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
			U0.S[0] = U0.D * s0;
			U0.Tau = U0.D * (UL.Tau / VL.rho + (s0 - VL.v[0]) * (s0 + VL.p / wL));
			for (int k = 1; k < Ndim; k++) {
				U0.S[k] = U0.D * VL.v[k];
			}
			rc.flux += sL * (U0 - UL);
		}
	} else {
		rc.flux = VR.toFlux(0);
		if (sR > zero) {
			U0.D = wR / (sR - s0);
			U0.S[0] = U0.D * s0;
			U0.Tau = U0.D * (UR.Tau / VR.rho + (s0 - VR.v[0]) * (s0 + VR.p / wR));
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
