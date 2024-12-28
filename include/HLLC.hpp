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

using namespace Math;

template<typename Type, int Ndim>
struct PrimitiveState;

template<typename Type, int Ndim>
struct ConservedState {
	static constexpr int Nfield = Ndim + 2;
	ConservedState() :
			stateVector(), massDensity(stateVector[0]), energyDensity(stateVector[2]), momentumDensity(
					*((Vector<Type, Ndim>*) (&stateVector[1]))) {

	}
	ConservedState(Vector<Type, Nfield> const &vec) :
			stateVector(vec), massDensity(stateVector[0]), energyDensity(stateVector[2]), momentumDensity(
					*((Vector<Type, Ndim>*) (&stateVector[1]))) {
	}
	ConservedState& operator=(PrimitiveState<Type, Ndim> const &other);
	ConservedState& operator=(ConservedState<Type, Ndim> const &other) {
		stateVector = other.stateVector;
		return *this;
	}
	operator Vector<Type, Nfield> const&() const {
		return stateVector;
	}
	operator Vector<Type, Nfield>&() {
		return stateVector;
	}
	ConservedState operator+(ConservedState const &B) const {
		ConservedState const &A = *this;
		ConservedState C;
		C.stateVector = A.stateVector + B.stateVector;
		return C;
	}
	ConservedState operator-(ConservedState const &B) const {
		ConservedState const &A = *this;
		ConservedState C;
		C.stateVector = A.stateVector - B.stateVector;
		return C;
	}
	ConservedState& operator+=(ConservedState const &other) {
		*this = *this + other;
		return *this;
	}
	ConservedState& operator-=(ConservedState const &other) {
		*this = *this - other;
		return *this;
	}
	friend ConservedState operator*(Type A, ConservedState const &B) {
		ConservedState C;
		C.stateVector = A * B.stateVector;
		return C;
	}
	Vector<Type, Nfield> stateVector;
	Type &massDensity;
	Type &energyDensity;
	Vector<Type, Ndim> &momentumDensity;
};

template<typename Type, int Ndim>
struct EquationOfState {
	static const Type fluidGamma;
	Type pressure(PrimitiveState<Type, Ndim> const &state) const {
		return (fluidGamma - Type(1)) * state.massDensity * state.specificEnergy;
	}
	Type soundSpeed(PrimitiveState<Type, Ndim> const &state) const {
		return sqrt(fluidGamma * pressure(state) / state.massDensity);
	}
};

template<typename Type, int Ndim>
const Type EquationOfState<Type, Ndim>::fluidGamma = Type(5.0 / 3.0);

template<typename Type, int Ndim>
struct PrimitiveState {
	EquationOfState<Type, Ndim> eOS;
	static constexpr int Nfield = Ndim + 2;
	PrimitiveState(ConservedState<Type, Ndim> const &con) {
		massDensity = con.massDensity;
		auto const volumeDensity = Type(1) / massDensity;
		velocity = con.momentumDensity * volumeDensity;
		specificEnergy = volumeDensity * con.energyDensity - Type(0.5) * vectorDotProduct(velocity, velocity);
	}
	ConservedState<Type, Ndim> getFlux() const {
		ConservedState<Type, Ndim> flux;
		flux.massDensity = massDensity * velocity[0];
		for (int dim = 0; dim < Ndim; dim++) {
			flux.momentumDensity[dim] = flux.massDensity * velocity[dim];
		}
		auto const pressure = eOS.pressure(*this);
		flux.momentumDensity[0] += pressure;
		flux.energyDensity = flux.massDensity * specificEnergy;
		flux.energyDensity += pressure * velocity[0];
		flux.energyDensity += Type(0.5) * flux.massDensity * vectorDotProduct(velocity, velocity);
		return flux;
	}
	union {
		struct {
			Type massDensity;
			Vector<Type, Ndim> velocity;
			Type specificEnergy;
		};
		Vector<Type, Nfield> stateVector;
	};
};

template<typename Type, int Ndim>
ConservedState<Type, Ndim>& ConservedState<Type, Ndim>::operator=(PrimitiveState<Type, Ndim> const &prim) {
	massDensity = prim.massDensity;
	momentumDensity = prim.velocity * massDensity;
	energyDensity = massDensity * prim.specificEnergy + Type(0.5) * vectorDotProduct(prim.velocity, momentumDensity);
	return *this;
}

template<typename Type, int Ndim>
struct RiemannReturn {
	RiemannReturn(ConservedState<Type, Ndim> const &f) :
			flux(f) {
	}
	ConservedState<Type, Ndim> flux;
	Type soundSpeed;
};

template<typename Type, int Ndim>
struct HLLCSolver {
	RiemannReturn<Type, Ndim> operator()(ConservedState<Type, Ndim> const &lState,
			ConservedState<Type, Ndim> const &rState) const {
		PrimitiveState<Type, Ndim> const lPrimitives(lState);
		PrimitiveState<Type, Ndim> const rPrimitives(rState);
		ConservedState<Type, Ndim> flux;
		Type const uR = rPrimitives.velocity[0];
		Type const uL = lPrimitives.velocity[0];
		Type const aR = eOS.soundSpeed(rPrimitives);
		Type const aL = eOS.soundSpeed(lPrimitives);
		Type const pR = eOS.pressure(rPrimitives);
		Type const pL = eOS.pressure(lPrimitives);
		Type const sR = std::max(uR + aR, uL + aL);
		Type const sL = std::min(uR - aR, uL - aL);
		Type const rhoR = rState.massDensity;
		Type const rhoL = lState.massDensity;
		Type const eGasR = rState.energyDensity;
		Type const eGasL = lState.energyDensity;
		Type const weightL = rhoL * (sL - uL);
		Type const weightR = rhoR * (sR - uR);
		Type const weightSum = weightL - weightR;
		Type const weightInv = Type(1) / (weightSum + Type(1e-100));
		Type const sStar = ((pR - pL) + weightL * uL - weightR * uR) * weightInv;
		if ((Type(0) <= sL) || (Type(0) >= sR)) {
			if (Type(0) <= sL) {
				flux = lPrimitives.getFlux();
			} else {
				flux = rPrimitives.getFlux();
			}
		} else {
			if (Type(0) <= sStar) {
				Type const rhoStarL = rhoL * (sL - uL) / (sL - sStar);
				Type const momStarL = rhoStarL * sStar;
				Type const energyStarL = rhoStarL * (eGasL / rhoL + (sStar - uL) * (sStar + pL / (rhoL * (sL - uL))));
				Type const momL = lState.momentumDensity[0];
				flux = lPrimitives.getFlux();
				flux.massDensity += sL * (rhoStarL - rhoL);
				flux.momentumDensity += sL * (momStarL - momL);
				flux.energyDensity += sL * (energyStarL - eGasL);
			} else {
				Type const rhoStarR = rhoR * (sR - uR) / (sR - sStar);
				Type const momStarR = rhoStarR * sStar;
				Type const energyStarR = rhoStarR * (eGasR / rhoR + (sStar - uR) * (sStar + pR / (rhoR * (sR - uR))));
				Type const momR = rState.momentumDensity[0];
				flux = rPrimitives.getFlux();
				flux.massDensity += sR * (rhoStarR - rhoR);
				flux.momentumDensity += sR * (momStarR - momR);
				flux.energyDensity += sR * (energyStarR - eGasR);
			}
		}
		RiemannReturn<Type, Ndim> rc(flux);
		rc.soundSpeed = max(abs(uR) + aR, abs(uL) + aL);
		return rc;
	}
	EquationOfState<Type, Ndim> eOS;
};

template<typename Type, int Ndim>
struct RoeSolver {
	RiemannReturn<Type, Ndim> operator()(ConservedState<Type, Ndim> const &lState,
			ConservedState<Type, Ndim> const &rState) const {
		static constexpr int Nfield = Ndim + 2;
		PrimitiveState<Type, Ndim> lPrimitives(lState);
		PrimitiveState<Type, Ndim> rPrimitives(rState);
		ConservedState<Type, Ndim> flux;
		assert(rPrimitives.specificEnergy >= Type(0));
		assert(lPrimitives.specificEnergy >= Type(0));
		Type const one = Type(1.0);
		Type const half = Type(0.5);
		Type const fGamma = eOS.fluidGamma;
		Type const rhoL = lState.massDensity;
		Type const rhoR = rState.massDensity;
		Type const egasL = lState.energyDensity;
		Type const egasR = rState.energyDensity;
		Type const pL = eOS.pressure(lPrimitives);
		Type const pR = eOS.pressure(rPrimitives);
		Type const vL = lPrimitives.velocity[0];
		Type const vR = rPrimitives.velocity[0];
		Type const hL = (pL + egasL) / rhoL;
		Type const hR = (pR + egasR) / rhoR;
		Type const wR = sqrt(rhoR);
		Type const wL = sqrt(rhoL);
		Type const v_ = (wR * vR + wL * vL) / (wR + wL);
		Type const h_ = (wR * hR + wL * hL) / (wR + wL);
		Type const kin_ = half * sqr(v_);
		assert((h_ - kin_) > Type(0));
		Type const c2_ = (fGamma - one) * (h_ - kin_);
		Type const c_ = sqrt(c2_);
		Math::SquareMatrix<Type, Nfield> R;
		Math::SquareMatrix<Type, Nfield> Rinv;
		Vector<Type, Nfield> lambda;
		Matrix<Type, Nfield, 1> dU;
		dU[0, 0] = rhoR - rhoL;
		dU[1, 0] = rState.momentumDensity[0] - lState.momentumDensity[0];
		dU[2, 0] = egasR - egasL;
		auto Lambda = identityMatrix<Type, Nfield>();
		Lambda[0, 0] = abs(v_ - c_);
		Lambda[1, 1] = abs(v_);
		Lambda[2, 2] = abs(v_ + c_);
		R[0, 0] = one;
		R[1, 0] = v_ - c_;
		R[2, 0] = h_ - v_ * c_;
		R[0, 1] = one;
		R[1, 1] = v_;
		R[2, 1] = kin_;
		R[0, 2] = one;
		R[1, 2] = v_ + c_;
		R[2, 2] = h_ + v_ * c_;
		Rinv = matrixInverse(R);
		flux = half * (rPrimitives.getFlux() + lPrimitives.getFlux());
		flux.stateVector -= toVector(half * R * Lambda * (Rinv * dU));
		RiemannReturn<Type, Ndim> rc(flux);
		rc.soundSpeed = abs(v_) + c_;
		return rc;
	}
	EquationOfState<Type, Ndim> eOS;
};

#endif /* INCLUDE_HLLC_HPP_ */
