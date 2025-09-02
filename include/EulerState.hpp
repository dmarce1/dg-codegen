#pragma once

#include "ContainerArithmetic.hpp"
#include "Matrix.hpp"
#include "Util.hpp"
#include <array>
#include <cmath>
#include <vector>
#include <string>

#define EULERS_CONSTRUCTION : \
	S  (*std::launder(reinterpret_cast<std::array<T, D>*>(base_type::data()))), \
	rho(*std::launder(reinterpret_cast<T*>(base_type::data() + D))), \
	eg (*std::launder(reinterpret_cast<T*>(base_type::data() + 1 + D)))

enum class RiemannSolver : int {
	LLF, HLL, HLLC
};

template<typename T, int D, RiemannSolver riemannSolver = RiemannSolver::LLF>
struct EulerState: public std::array<T, 2 + D> {
	static_assert(D <= 3);
	using EleType = typename ElementType<T>::type;
	static constexpr int NF = 2 + D;
	static constexpr EleType zero = EleType(0);
	static constexpr EleType one = EleType(1);
	static constexpr EleType two = EleType(2);
	static constexpr EleType three = EleType(3);
	static constexpr EleType four = EleType(4);
	static constexpr EleType five = EleType(5);
	static constexpr EleType quarter = one / four;
	static constexpr EleType third = one / three;
	static constexpr EleType half = one / two;
	static constexpr EleType gamma = five * third;
	static constexpr EleType gamm1 = gamma - one;
	static constexpr EleType gamp1o2gam = (gamma + one) / (two * gamma);
	static constexpr EleType igamm1 = one / gamm1;
	static constexpr int sx_i = 0;
	static constexpr int rho_i = D;
	static constexpr int eg_i = D + 1;
	using base_type = std::array<T, NF>;
	using value_type = T;
	using eigensys_type = std::pair<std::array<T, NF>, SquareMatrix<T, NF>>;
	static int dimCount() noexcept {
		return D;
	}
	static constexpr int fieldCount() noexcept {
		return NF;
	}
	EulerState() EULERS_CONSTRUCTION {
	}
	EulerState(base_type const &other) EULERS_CONSTRUCTION {
		((base_type&)(*this)).operator=(other);
	}
	EulerState(EulerState const &other) EULERS_CONSTRUCTION {
		*this = other;
	}
	EulerState(T const &other) EULERS_CONSTRUCTION {
		for( auto& u : *this) {
			u = other;
		}
	}
	EulerState(EulerState &&other) EULERS_CONSTRUCTION {
		*this = std::move(other);
	}
	EulerState& operator=(EulerState const &other) {
		static_cast<base_type&>(*this) = static_cast<base_type const&>(other);
		return *this;
	}
	EulerState& operator=(EulerState &&other) {
		static_cast<base_type&>(*this) = std::move(static_cast<base_type&&>(other));
		return *this;
	}
	std::array<T, NF> eigenvalues(int dim) const {
		std::array<T, NF> lambda;
		T const irho = one / rho;
		T const v = S[dim] * irho;
		T ek = half * irho * S[0] * S[0];
		for(int d2 = 1; d2 < D; d2++) {
			ek += half * irho* S[d2] * S[d2];
		}
		T const ei = eg - ek;
		T const p = gamm1 * ei;
		T const a = sqrt(gamma * p * irho);
		std::fill(lambda.begin(), lambda.end(), v);
		lambda[dim] -= a;
		lambda[D + 1] += a;
		return lambda;

	}
	eigensys_type eigenSystem(int dim) const {
		eigensys_type rc;
		auto& eigenvalues = rc.first;
		auto& eigenvectors = rc.second;
		for(int r = 0; r < NF; r++) {
			eigenvalues[r] = allocateLike((*this)[0]);
			for(int c = 0; c < NF; c++) {
				eigenvectors(r, c) = allocateLike((*this)[0]);
			}
		}
		T const irho = one / rho;
		std::array<T, D> u;
		for(int d2 = 0; d2 < D; d2++) {
			u[d2] = S[d2] * irho;
		}
		if(dim != D - 1) {
			std::swap(u[dim], u[D - 1]);
		}
		T ek = half * u[0] * u[0];
		for(int d2 = 1; d2 < D; d2++) {
			ek += half * u[d2] * u[d2];
		}
		T const ei = eg * irho - ek;
		T const p = gamm1 * rho * ei;
		T const a = sqrt(gamma * p * irho);
		T const h = (p + eg) * irho;
		T const v = u[D - 1];
		std::fill(eigenvalues.begin(), eigenvalues.end(), v);
		eigenvalues[D - 1] -= a;
		eigenvalues[D + 1] += a;
		for(int r = 0; r < D; r++) {
			for(int c = 0; c < D - 1; c++) {
				eigenvectors(r, c) = (r == c) ? one : zero;
			}
			eigenvectors(r, D - 1) = u[r];
			eigenvectors(r, D) = u[r];
			eigenvectors(r, D + 1) = u[r];
		}
		for(int c = 0; c < D - 1; c++) {
			eigenvectors(D, c) = zero;
			eigenvectors(D + 1, c) = u[c];
		}
		for(int c = D - 1; c <= D + 1; c++) {
			eigenvectors(D, c) = one;
		}
		eigenvectors(D - 1, D - 1) -= a;
		eigenvectors(D - 1, D + 1) += a;
		eigenvectors(D + 1, D - 1) = h - a * v;
		eigenvectors(D + 1, D) = ek;
		eigenvectors(D + 1, D + 1) = h + a * v;
		if( dim != D - 1 ) {
			for(int c = 0; c < NF; c++) {
				std::swap(eigenvectors(D - 1, c), eigenvectors(dim, c));
			}
		}
		return rc;
	}
	EulerState flux(int dim) const {
		EulerState F;
		T const irho = one / rho;
		std::array<T, D> v;
		for(int d2 = 0; d2 < D; d2++) {
			v[d2] = S[d2] * irho;
		}
		T ek = half * v[0] * S[0];
		for(int d2 = 1; d2 < D; d2++) {
			ek += half * v[d2] * S[d2];
		}
		T const ei = eg - ek;
		T const p = gamm1 * ei;
		T const u = v[dim];
		F.rho = S[dim];
		F.eg = u * (eg + p);
		for(int d2 = 0; d2 < D; d2++) {
			F.S[d2] = u * S[d2];
		}
		F.S[dim] += p;
		return F;
	}
	friend EulerState solveRiemannProblem(const EulerState &uL, const EulerState &uR, int dim) {
		using namespace Math;
		EulerState F;
		T const irhoR = one / uR.rho;
		T const irhoL = one / uL.rho;
		T const vR = irhoR * uR.S[dim];
		T const vL = irhoL * uL.S[dim];
		T ekR, ekL;
		ekR = half * irhoR * uR.S[0] * uR.S[0];
		ekL = half * irhoL * uL.S[0] * uL.S[0];
		for(int d2 = 1; d2 < D; d2++) {
			ekR += half * irhoR * uR.S[d2] * uR.S[d2];
			ekL += half * irhoL * uL.S[d2] * uL.S[d2];
		}
		T const eiR = uR.eg - ekR;
		T const eiL = uL.eg - ekL;
		T const pR = gamm1 * eiR;
		T const pL = gamm1 * eiL;
		T const aR = sqrt(gamma * pR * irhoR);
		T const aL = sqrt(gamma * pL * irhoL);
		EulerState const fL = uL.flux(dim);
		EulerState const fR = uR.flux(dim);
		if constexpr(riemannSolver == RiemannSolver::LLF) {
			T const sStar = max(aL + abs(vL), aR + abs(vR));
			for(int fi = 0; fi < NF; fi++) {
				F[fi] = (fL[fi] + fR[fi] - sStar * (uR[fi] - uL[fi])) * half;
			}
		} else if constexpr(riemannSolver == RiemannSolver::HLL) {
			T const sL = min(zero, min(T(vL - aL), T(vR - aR)));
			T const sR = max(zero, max(T(vL + aL), T(vR + aR)));
			T const w = EleType(1) / (sR - sL);
			for(int fi = 0; fi < NF; fi++) {
				F[fi] = w * (sR * fL[fi] - sL * fR[fi] + sL * sR * (uR[fi] - uL[fi]));
			}
		} else/*if constexpr(riemannSolver == RiemannSolver::HLLC)*/{
			EulerState fStar, uK, fK;
			T const pStar = max(zero, T(half * (pR + pL - quarter * (vR - vL) * (aR + aL) * (uR.rho + uL.rho))));
			T const qL = sqrt(one + gamp1o2gam * max(zero, T(pStar / pL - one)));
			T const qR = sqrt(one + gamp1o2gam * max(zero, T(pStar / pR - one)));
			T const sL = vL - qL * aL;
			T const sR = vR + qR * aR;
			T const num = pR - pL + vL * uL.rho * (sL - vL) - vR * uR.rho * (sR - vR);
			T const den = uL.rho * (sL - vL) - uR.rho * (sR - vR);
			T const sStar = num / den;
			T const pLR = half * (pR + pL + uL.rho * (sL - vL) * (sStar - vL) + uR.rho * (sR - vR) * (sStar - vR));
			T const wL = half + copysign(half, sStar);
			T const wR = half - copysign(half, sStar);
			T const sK = wL * sL + wR * sR;
			for(int fi = 0; fi < NF; fi++) {
				uK[fi] = wL * uL[fi] + wR * uR[fi];
				fK[fi] = wL * fL[fi] + wR * fR[fi];
			}
			T const iden = one / (sK - sStar);
			for(int fi = 0; fi < NF; fi++) {
				fStar[fi] = sStar * (sK * uK[fi] - fK[fi]) * iden;
			}
			T const tmp = sK * pLR * iden;
			fStar[dim] += tmp;
			fStar.back() += sStar * tmp;
			T const sStarK = sStar * sK;
			T const wD = (half - copysign(half, T(half * sStar * sK)));
			for(int fi = 0; fi < NF; fi++) {
				F[fi] = wD * fStar[fi] + (one - wD) * fK[fi];
			}
		}
		return F;
	}
	bool sanityCheck() const {
		return true;
	}
	friend T findPositivityPreservingTheta(EulerState const &u0, EulerState const &uh) {
		using EleType = typename ElementType<T>::type;
		constexpr EleType epsRho = EleType(1e-13);
		constexpr EleType tiny = EleType(sqrt(std::numeric_limits<EleType>::min()));
		auto const s0 = u0.S;
		auto const s1 = uh.S;
		auto const dS = s1 - s0;
		T const rho0 = u0.rho;
		T const rho1 = uh.rho;
		T const dRho = rho1 - rho0;
		T const E0 = u0.eg;
		T const E1 = uh.eg;
		auto const dE = E1 - E0;
		auto const maskRho = dRho < zero;
		T thetaRho = safeDiv(T(epsRho - rho0), dRho);
		thetaRho[!maskRho] = one;
		T const a = dE * dRho - half * dot(dS, dS);
		T const b = E0 * dRho + dE * rho0 - dot(s0, dS);
		T const c = E0 * rho0 - half * dot(s0, s0);
		T const g1 = a + b + c;
		auto const needP = (g1 < zero);
		auto const skipP = !needP;
		T disc = b * b - four * a * c;
		disc[skipP || (disc < zero)] = zero;
		T const sqrtDisc = sqrt(disc);
		T const thetaLin = -safeDiv(c, b);
		T const denom = -b - copysign(sqrtDisc, b);
		T const thetaQuad = two * safeDiv(c, denom);
		T thetaP = thetaQuad;
		auto const linFlag = (abs(a) < tiny) && needP;
		thetaP[linFlag] = thetaLin[linFlag];
		thetaP[skipP] = one;
		T const theta = clamp(zero, T(min(thetaRho, thetaP)), one);
		return theta;
	}
	T const& getDensity() const {
		return rho;
	}
	T const& getEnergy() const {
		return eg;
	}
	std::array<T, D> const& getMomentum() const {
		return S;
	}
	T const& getMomentum(int d) const {
		return S[d];
	}
	void setDensity(T const& value) {
		rho = value;
	}
	void setEnergy(T const& value) {
		eg = value;
	}
	void setMomentum(std::array<T, D> const& value) {
		S = value;
	}
	void setMomentum(int d, T const& value) {
		S[d] = value;
	}
	static std::vector<std::string> getFieldNames() {
		static std::vector<std::string> const fieldNames = []() {
			std::vector<std::string> names;
			for(int d = 0; d < D; d++) {
				std::string const name = std::string("s_") + std::string(1, 'x' + d);
				names.push_back(name);
			}
			names.push_back("rho");
			names.push_back("eg");
			return names;
		}();
		return fieldNames;
	}
private:
	std::array<T, D>& S;
	T& rho;
	T& eg;
};

template<typename T, int D>
struct CanDoArithmetic<EulerState<T, D>> {
	static constexpr bool value = true;
};

template<typename T, int D>
EulerState<T, D> initSodShockTube(std::array<T, D> x) {
	using EleType = typename ElementType<T>::type;
	static constexpr EleType half = EleType(0.5);
	/*********************************************************/
	static constexpr T rhoL = EleType(1.000);
	static constexpr T pL = EleType(1.000);
	static constexpr T vL = EleType(0.000);
	static constexpr T rhoR = EleType(0.125);
	static constexpr T pR = EleType(0.100);
	static constexpr T vR = EleType(0.000);
	/*********************************************************/
	static constexpr T c0 = EulerState<T, D>::igamm1;
	static constexpr T eL = c0 * pL + half * rhoL * sqr(vL);
	static constexpr T eR = c0 * pR + half * rhoR * sqr(vR);
	EulerState<T, D> u;
	T const r = x[0];
	u.setMomentum(zeroArray<T, D>());
	if (r < 0.25 || 0.75 < r) {
		u.setDensity(rhoL);
		u.setEnergy(eL);
		u.setMomentum(0, vL * rhoL);
	} else {
		u.setDensity(rhoR);
		u.setEnergy(eR);
		u.setMomentum(0, vR * rhoR);
	}
	return u;
}

template<typename T, int D>
using EulerStateLLF = EulerState<T, D, RiemannSolver::LLF>;

template<typename T, int D>
using EulerStateHLL = EulerState<T, D, RiemannSolver::HLL>;

template<typename T, int D>
using EulerStateHLLC = EulerState<T, D, RiemannSolver::HLLC>;

