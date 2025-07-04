#pragma once

#include "ContainerArithmetic.hpp"
#include "Matrix.hpp"
#include "Util.hpp"
#include "ValArray.hpp"

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
		using std::sqrt;
		T const irho = one / rho;
		T const v = S[dim] * irho;
		T const ek = half * irho * dot(S, S);
		T const ei = max(eg - ek, zero);
		T const p = gamm1 * ei;
		T const a = sqrt(gamma * p * irho);
		std::fill(lambda.begin(), lambda.end(), v);
		lambda[D - 1] -= a;
		lambda[D + 1] += a;
		return lambda;

	}
	eigensys_type eigenSystem(int dim) const {
		eigensys_type rc;
		auto& eigenvalues = rc.first;
		auto& eigenvectors = rc.second;
		T const irho = one / rho;
		std::array<T, D> u = S * irho;
		if(dim != D - 1) {
			std::swap(u[dim], u[D - 1]);
		}
		T const ek = half * dot(u, u);
		T const ei = max(zero, eg * irho - ek);
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
	SquareMatrix<T, NF> rightEigenvectors() const {
		std::array<std::array<T, D>, D> n;
		T const magS = sqrt(dot(S, S));
		n[0] = S / (magS + EleType(std::numeric_limits<double>::min()));
		T theta = dot(n[0], n[0]);
		n[0][0] += one - theta;
		if constexpr (D == 2) {
			n[1] = n[0];
			std::swap(n[1][0], n[1][1]);
			n[1][0] = -n[1][0];
		} else if constexpr (D == 3) {
			auto const maskT = (sqr(n[0][0]) + sqr(n[0][1]) > two * third);
			auto const maskF = !maskT;
			n[1][0][maskT] = zero;
			n[1][0][maskF] = one;
			n[1][1] = zero;
			n[1][2][maskT] = one;
			n[1][2][maskF] = zero;
			n[1] -= n[0] * dot(n[1], n[0]);
			n[1] = n[1] / norm(n[1]);
			n[2][0] = +n[0][1] * n[1][2] - n[0][2] * n[1][1];
			n[2][1] = -n[0][0] * n[1][2] + n[0][2] * n[1][0];
			n[2][2] = +n[0][0] * n[1][1] - n[0][1] * n[1][0];
			n[2] = n[2] / norm(n[2]);
		}
		T const irho = one / rho;
		T const v = magS * irho;
		T const ek = half * v * magS;
		T const ei = max(zero, eg - ek);
		T const p = gamm1 * ei;
		T const a = sqrt(gamma * p * irho);
		T const h = (eg + p) * irho;
		SquareMatrix<T, NF> R;
		for(int d = 0; d < D; d++) {
			R(d, 0) = (v - a) * n[0][d];
			R(d, 1) = v * n[0][d];
			R(d, 2) = (v + a) * n[0][d];
		}
		R(D + 0, 0) = R(D + 0, 1) = R(D + 0, 2) = one;
		R(D + 1, 0) = h - a * v;
		R(D + 1, 1) = sqr(v) * half;
		R(D + 1, 2) = h + a * v;
		if constexpr(D == 2) {
			for(int d = 0; d < D; d++) {
				R(d, 3) = n[1][d];
			}
			R(D + 0, 3) = R(D + 1, 3) = zero;
		} else if constexpr(D == 3) {
			for(int d = 0; d < D; d++) {
				R(d, 4) = n[2][d];
			}
			R(D + 0, 4) = R(D + 1, 4) = zero;
		}
		return R;
	}
	EulerState flux(int d1) const {
		EulerState F;
		T const irho = one / rho;
		auto const v = S * irho;
		T const ek = half * dot(v, S);
		T const ei = max(zero, eg - ek);
		T const p = gamm1 * ei;
		T const u = v[d1];
		F.rho = S[d1];
		F.eg = u * (eg + p);
		for(int d2 = 0; d2 < D; d2++) {
			F.S[d2] = u * S[d2];
		}
		F.S[d1] += p;
		return F;
	}
	friend EulerState solveRiemannProblem(const EulerState &uL, const EulerState &uR, int dim) {
		EulerState F;
		T const irhoR = one / uR.rho;
		T const irhoL = one / uL.rho;
		T const vR = irhoR * uR.S[dim];
		T const vL = irhoL * uL.S[dim];
		T const ekR = half * irhoR * dot(uR.S, uR.S);
		T const ekL = half * irhoL * dot(uL.S, uL.S);
		T const eiR = max(zero, uR.eg - ekR);
		T const eiL = max(zero, uL.eg - ekL);
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
			T const sL = min(zero, min(vL - aL, vR - aR));
			T const sR = max(zero, max(vL + aL, vR + aR));
			T const w = EleType(1) / (sR - sL);
			for(int fi = 0; fi < NF; fi++) {
				F[fi] = w * (sR * fL[fi] - sL * fR[fi] + sL * sR * (uR[fi] - uL[fi]));
			}
		} else/*if constexpr(riemannSolver == RiemannSolver::HLLC)*/{
			EulerState fStar, uK, fK;
			T const pStar = max(zero, half * (pR + pL - quarter * (vR - vL) * (aR + aL) * (uR.rho + uL.rho)));
			T const qL = sqrt(one + gamp1o2gam * max(zero, pStar / pL - one));
			T const qR = sqrt(one + gamp1o2gam * max(zero, pStar / pR - one));
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
			T const wD = (one - copysign(one, sStar * sK)) * half;
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
		return one;
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
		T const dE = E1 - E0;
		T const maskRho = (dRho < zero).template cast<EleType>();
		T const thetaRho = maskRho * (safeDiv(T(epsRho - rho0), dRho) - one) + one;
		T const a = dE * dRho - half * dot(dS, dS);
		T const b = E0 * dRho + dE * rho0 - dot(s0, dS);
		T const c = E0 * rho0 - half * dot(s0, s0);
		T const g1 = a + b + c;
		auto const needP = (g1 < zero).template cast<EleType>();
		auto const skipP = one - needP;
		T const disc = max(zero, needP * (b * b - four * a * c));
		T const sqrtDisc = sqrt(disc);
		auto const maskLin = needP * (abs(a) < tiny).template cast<EleType>();;
		auto const maskQuad = needP - maskLin;
		T const thetaLin = -safeDiv(c, b);
		T const denom = -b - copysign(sqrtDisc, b);
		T const thetaQuad = two * safeDiv(c, denom);
		T const thetaP = skipP * one + maskLin * thetaLin + maskQuad * thetaQuad;
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
	static constexpr EleType one = EleType(1);
	/*********************************************************/
	static constexpr T rhoL = one;
	static constexpr T pL = one;
	static constexpr T rhoR = T(0.125);
	static constexpr T pR = T(0.1);
	/*********************************************************/
	static constexpr T c0 = one / (EulerState<T, D>::gamma - one);
	static constexpr T eL = c0 * pL;
	static constexpr T eR = c0 * pR;
	EulerState<T, D> u;
	T const r = x[0];
	//T const r = T(1.0/2.0) * (x[0] - x[1] + T(1));
	u.setMomentum(zeroArray<T, D>());
//	printf( "%e %e %e \n", r, x[0],  x[1]);
	if (r < 0.25 || 0.75 < r) {
		u.setDensity(rhoL);
		u.setEnergy(eL);
	} else {
		u.setDensity(rhoR);
		u.setEnergy(eR);
	}
	return u;
}

template<typename T, int D>
using EulerStateLLF = EulerState<T, D, RiemannSolver::LLF>;

template<typename T, int D>
using EulerStateHLL = EulerState<T, D, RiemannSolver::HLL>;

template<typename T, int D>
using EulerStateHLLC = EulerState<T, D, RiemannSolver::HLLC>;

