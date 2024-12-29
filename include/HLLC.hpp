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
using HydroState = Vector<Type, Ndim + 2>;

constexpr int D_i = 0;
constexpr int S_i = 1;
constexpr int E_i = 2;

constexpr int rho_i = 0;
constexpr int vel_i = 1;
constexpr int eps_i = 2;

#define fGamma (7.0 / 5.0)
#define Nfield 3

template<typename Type, int Ndim>
struct RiemannReturn {
	HydroState<Type, Ndim> flux;
	Type signalSpeed;
};

template<typename Type, int Ndim>
HydroState<Type, Ndim> conToPrim(HydroState<Type, Ndim> const &con) {
	HydroState<Type, Ndim> prim;
	auto const half = Type(0.5);
	auto const rho = con[D_i];
	auto const rhoInv = Type(1) / rho;
	auto const v = con[S_i] * rhoInv;
	auto const eKin = half * rho * sqr(v);
	auto const eps = (con[E_i] - eKin) * rhoInv;
	prim[rho_i] = rho;
	prim[vel_i] = v;
	prim[eps_i] = eps;
	return prim;
}

template<typename Type, int Ndim>
HydroState<Type, Ndim> primToCon(HydroState<Type, Ndim> const &prim) {
	HydroState<Type, Ndim> con;
	auto const half = Type(0.5);
	auto const D = prim[rho_i];
	auto const S = D * prim[vel_i];
	auto const E = D * prim[eps_i] + half * sqr(S) / D;
	con[D_i] = D;
	con[S_i] = S;
	con[E_i] = E;
	return con;
}

template<typename Type, int Ndim>
HydroState<Type, Ndim> primToFlux(const HydroState<Type, Ndim> &prim) {
	HydroState<Type, Ndim> F;
	auto const one = Type(1);
	auto const half = Type(0.5);
	auto const gam = Type(fGamma);
	auto const rho = prim[rho_i];
	auto const v = prim[vel_i];
	auto const eps = prim[eps_i];
	auto const v2 = sqr(v);
	auto const e = rho * eps + half * rho * v2;
	auto const p = (gam - one) * rho * eps;
	F[D_i] = rho * v;
	F[S_i] = rho * v2 + p;
	F[E_i] = v * (e + p);
	return F;
}

template<typename Type, int Ndim>
SquareMatrix<Type, Nfield> primToEigenvectors(const HydroState<Type, Ndim> &prim) {
	SquareMatrix<Type, Nfield> R;
	auto const one = Type(1);
	auto const half = Type(0.5);
	auto const gam = Type(fGamma);
	auto const v = prim[vel_i];
	auto const eps = prim[eps_i];
	auto const v2 = sqr(v);
	auto const h = gam * eps + half * v2;
	auto const c2 = gam * (gam - one) * eps;
	auto const c = sqrt(c2);
	R[0, 0] = one;
	R[1, 0] = v + c;
	R[2, 0] = h + v * c;
	R[0, 1] = one;
	R[1, 1] = v;
	R[2, 1] = half * v2;
	R[0, 2] = one;
	R[1, 2] = v - c;
	R[2, 2] = h - v * c;
	return R;
}

template<typename Type, int Ndim>
SquareMatrix<Type, Nfield> roeToEigenvectors(const HydroState<Type, Ndim> &roe) {
	SquareMatrix<Type, Nfield> R;
	auto const one = Type(1);
	auto const half = Type(0.5);
	auto const gam = Type(fGamma);
	auto const v = roe[vel_i];
	auto const h = roe[eps_i];
	auto const c = sqrt((gam - one) * (h - half * sqr(v)));
	R[0, 0] = one;
	R[1, 0] = v + c;
	R[2, 0] = h + v * c;
	R[0, 1] = one;
	R[1, 1] = v;
	R[2, 1] = half * sqr(v);
	R[0, 2] = one;
	R[1, 2] = v - c;
	R[2, 2] = h - v * c;
	return R;
}

template<typename Type, int Ndim>
HydroState<Type, Ndim> stateToChar(const HydroState<Type, Ndim> &dU, const HydroState<Type, Ndim> &prim) {
	return toVector(matrixInverse(primToEigenvectors<Type, Ndim>(prim)) * toColumnVector(dU));
}
template<typename Type, int Ndim>
HydroState<Type, Ndim> primToRoe(const HydroState<Type, Ndim> &primL, const HydroState<Type, Ndim> &primR) {
	auto const half = Type(0.5);
	auto const gam = Type(fGamma);
	auto const rhoR = primR[rho_i];
	auto const vR = primR[vel_i];
	auto const epsR = primR[eps_i];
	auto const hR = gam * epsR + half * sqr(vR);
	auto const rhoL = primL[rho_i];
	auto const vL = primL[vel_i];
	auto const epsL = primL[eps_i];
	auto const hL = gam * epsL + half * sqr(vL);
	auto const wR = sqrt(rhoR);
	auto const wL = sqrt(rhoL);
	auto const rho_ = wR * wL;
	auto const v_ = (wR * vR + wL * vL) / (wR + wL);
	auto const h_ = (wR * hR + wL * hL) / (wR + wL);
	HydroState<Type, Ndim> dW = { rho_, v_, h_ };
	return dW;
}

template<typename Type, int Ndim>
HydroState<Type, Ndim> charToState(const HydroState<Type, Ndim> &dW, const HydroState<Type, Ndim> &prim) {
	return toVector(primToEigenvectors<Type, Ndim>(prim) * toColumnVector(dW));
}

template<typename Type, int Ndim>
struct HLLCSolver {
	RiemannReturn<Type, Ndim> operator()(HydroState<Type, Ndim> const &UL, HydroState<Type, Ndim> const &UR) const {
		RiemannReturn<Type, Ndim> rc;
		HydroState<Type, Ndim> U0;
		auto const gam = Type(fGamma);
		auto const two = Type(2);
		auto const one = Type(1);
		auto const half = Type(0.5);
		auto const zero = Type(0);
		auto const VL = conToPrim<Type, Ndim>(UL);
		auto const VR = conToPrim<Type, Ndim>(UR);
		Type const rhoL = VL[rho_i];
		Type const rhoR = VR[rho_i];
		Type const uL = VL[vel_i];
		Type const uR = VR[vel_i];
		Type const epsL = VL[eps_i];
		Type const epsR = VR[eps_i];
		Type const eL = UL[E_i];
		Type const eR = UR[E_i];
		Type const rhoInvL = one / rhoL;
		Type const rhoInvR = one / rhoR;
		Type const pL = (gam - one) * rhoL * epsL;
		Type const pR = (gam - one) * rhoR * epsR;
		Type const aL = sqrt(gam * pL * rhoInvL);
		Type const aR = sqrt(gam * pR * rhoInvR);
		Type const eta0 = (gam + one) / (two * gam);
		Type const z = (gam - one) / (two * gam);
		Type const pRz = pow(pR, z);
		Type const pLz = pow(pL, z);
		Type const pz = (aL + aR - half * (gam - one) * (uR - uL)) * pLz * pRz / (aL * pRz + aR * pLz);
		Type const p0 = pow(max(zero, pz), one / z);
		Type const qL = (p0 < pL) ? one : sqrt(one + eta0 * (p0 / pL - one));
		Type const qR = (p0 < pR) ? one : sqrt(one + eta0 * (p0 / pR - one));
		Type const sL = uL - qL * aL;
		Type const sR = uR + qR * aR;
		Type const wL = rhoL * (sL - uL);
		Type const wR = rhoR * (sR - uR);
		Type const s0 = (pR - pL + wL * uL - wR * uR) / (wL - wR);
		assert(sL <= s0);
		assert(s0 <= sR);
		if (s0 >= zero) {
			rc.flux = primToFlux<Type, Ndim>(VL);
			if (zero > sL) {
				U0[D_i] = rhoL * (sL - uL) / (sL - s0);
				U0[S_i] = U0[D_i] * s0;
				U0[E_i] = U0[D_i] * (eL / rhoL + (s0 - uL) * (s0 + pL / wL));
				rc.flux += sL * (U0 - UL);
			}
		} else {
			rc.flux = primToFlux<Type, Ndim>(VR);
			if (sR > zero) {
				U0[D_i] = rhoR * (sR - uR) / (sR - s0);
				U0[S_i] = U0[D_i] * s0;
				U0[E_i] = U0[D_i] * (eR / rhoR + (s0 - uR) * (s0 + pR / wR));
				rc.flux += sR * (U0 - UR);
			}
		}
		rc.signalSpeed = max(abs(sR), abs(sL));
		return rc;
	}
};
/*
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
 */
#endif /* INCLUDE_HLLC_HPP_ */
