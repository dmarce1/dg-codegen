#pragma once

#include "Vector.hpp"
#include <array>
#include <cmath>

#define EULERS_CONSTRUCTION       :                    \
	rho(*((T*) (base_type::data() + 0))),              \
	eg(*((T*) (base_type::data() + 1))),               \
	tau(*((T*) (base_type::data() + 2))),              \
	S(*((Math::Vector<T, D>*) (base_type::data() + 3)))

template<typename T, int D, T G = 5.0 / 3.0>
struct EulerState: public Math::Vector<T, 3 + D> {
	static constexpr T zero = T(0);
	static constexpr T half = T(0.5);
	static constexpr T one = T(1);
	static constexpr T gamma = G;
	static constexpr T igamma = one / gamma;
	static constexpr T gamm1 = gamma - one;
	static constexpr T des1 = T(1e-1);
	static constexpr T des2 = T(1e-3);
	using base_type = Math::Vector<T, 3 + D>;
	using eigensys_type = std::pair<std::array<T, 3 + D>, Math::SquareMatrix<T, 3 + D>>;
	static constexpr int fieldCount() noexcept {
		return 3 + D;
	}
	constexpr EulerState() EULERS_CONSTRUCTION {
	}
	constexpr EulerState(EulerState const &other) EULERS_CONSTRUCTION {
		*this = other;
	}
	constexpr EulerState(T const &other) EULERS_CONSTRUCTION {
		for( auto& u : *this) {
			u = other;
		}
	}
	constexpr EulerState(EulerState &&other) EULERS_CONSTRUCTION {
		*this = std::move(other);
	}
	constexpr EulerState& operator=(EulerState const &other) {
		static_cast<base_type&>(*this) = static_cast<base_type const&>(other);
		return *this;
	}
	constexpr EulerState& operator=(EulerState &&other) {
		static_cast<base_type&>(*this) = std::move(static_cast<base_type&&>(other));
		return *this;
	}
	constexpr eigensys_type eigenSystem(int dim) const {
		eigensys_type rc { {zero}, {zero}};
		auto& lambda = rc.first;
		auto& R = rc.second;
		T const irho = one / rho;
		auto const v = S * irho;
		T const ek = half * Math::vectorDotProduct(v, S);
		T const ei = ((eg - ek) < des2 * eg) ? (eg - ek) : std::pow(tau, gamma);
		T const p = gamm1 * ei;
		T const a = std::sqrt(gamma * p * irho);
		T const h = (ei + p) * irho;
		T const u = v[dim];
		for( int d = 0; d < D + 3; d++) {
			lambda[d] = u;
		}
		lambda[D + 0] -= a;
		lambda[D + 1] += a;
		for( int d = 0; d < D; d++) {
			R[d][D + 0] = v[d];
			R[d][D + 1] = v[d];
		}
		for( int d = 0; d < D; d++) {
			R[d][d] = one;
		}
		R[dim][D + 0] -= a;
		R[dim][D + 1] += a;
		R[D + 0][D + 0] = one;
		R[D + 0][D + 1] = one;
		R[D + 1][D + 0] = h - u * a;
		R[D + 1][D + 1] = h + u * a;
		for( int i = 0; i < D; i++) {
			R[D + 1][i] = v[i];
		}
		R[D + 2][D + 2] = one;
		return rc;
	}
	constexpr EulerState flux(int d) const noexcept {
		EulerState F;
		T const irho = one / rho;
		auto const v = S * irho;
		T const ek = half * Math::vectorDotProduct(v, S);
		T const ei = ((eg - ek) < des2 * eg) ? (eg - ek) : std::pow(tau, gamma);
		T const p = gamm1 * ei;
		T const u = v[d];
		F.rho = S[d];
		F.eg = u * (eg + p);
		F.tau = u * tau;
		F.S = u * S;
		F.S[d] += p;
		return F;
	}
	constexpr void syncEntropy() noexcept {
		T const irho = one / rho;
		auto const v = S * irho;
		T const ekin = half * Math::vectorDotProduct(v, S);
		T const eint = eg - ekin;
		if (eint > des1 * eg) {
			tau = std::pow(eint, igamma);
		}
	}
	constexpr EulerState toCharacteristic(int dim) const {
		return EulerState {matrixInverse(eigenSystem(dim).second) * (base_type const&)(*this)};
	}
	constexpr EulerState fromCharacteristic(int dim) const {
		return EulerState {eigenSystem(dim).second * (base_type const&)(*this)};
	}
	friend constexpr EulerState riemann(const EulerState &uL, const EulerState &uR, int d) noexcept {
		using state_type = EulerState;
		T const irhoL = one / uL.rho;
		T const irhoR = one / uR.rho;
		T const vL = uL.S[d] * irhoL;
		T const vR = uR.S[d] * irhoR;
		T const ekL = half * irhoL * Math::vectorDotProduct(uL.S, uL.S);
		T const ekR = half * irhoR * Math::vectorDotProduct(uR.S, uR.S);
		T const eiL = ((uL.eg - ekL) < des2 * uL.eg) ? (uL.eg - ekL) : std::pow(uL.tau, gamma);
		T const eiR = ((uR.eg - ekR) < des2 * uR.eg) ? (uR.eg - ekR) : std::pow(uR.tau, gamma);
		T const pL = gamm1 * eiL * uL.rho;
		T const pR = gamm1 * eiR * uR.rho;
		T const aL = std::sqrt(gamma * pL * irhoL);
		T const aR = std::sqrt(gamma * pR * irhoR);
		T const sL = std::min(vL - aL, vR - aR);
		T const sR = std::max(vL + aL, vR + aR);
		T const num = pR - pL + uL.S[d] * (sL - vL) - uR.S[d] * (sR - vR);
		T const den = uL.rho * (sL - vL) - uR.rho * (sR - vR);
		T const sM = num / den;
		if (sM > zero) {
			auto const f = uL.flux(d);
			if (sL > zero) {
				return f;
			} else {
				state_type u;
				T const rho = uL.rho * (sL - vL) / (sL - sM);
				u.rho = rho;
				u.S = uL.S;
				u.S[d] += rho * (sM - vL);
				u.eg = uL.eg + (sM - vL) * (rho * sM + pL / (sL - vL));
				u.tau = uL.tau * (sL - vL) / (sL - sM);
				return f + sL * (u - uL) / (sL - sM);
			}
		} else {
			auto const f = uR.flux(d);
			if (sR > zero) {
				state_type u;
				T const rho = uR.rho * (sR - vR) / (sR - sM);
				u.rho = rho;
				u.S = uR.S;
				u.S[d] += rho * (sM - vR);
				u.eg = uR.eg + (sM - vR) * (rho * sM + pR / (sR - vR));
				u.tau = uR.tau * (sR - vR) / (sR - sM);
				return f + sR * (u - uR) / (sR - sM);
			} else {
				return f;
			}
		}
	}
	constexpr T const& getDensity() const {
		return rho;
	}
	constexpr T const& getEnergy() const {
		return eg;
	}
	T const& getEntropy() const {
		return tau;
	}
	constexpr Math::Vector<T, D> const& getMomentum() const {
		return S;
	}
	constexpr T const& getMomentum(int d) const {
		return S[d];
	}
	constexpr void setDensity(T const& value) {
		rho = value;
	}
	constexpr void setEnergy(T const& value) {
		eg = value;
	}
	constexpr void setEntropy(T const& value) {
		tau = value;
	}
	constexpr void setMomentum(Math::Vector<T, D> const& value) {
		S = value;
	}
	constexpr void setMomentum(int d, T const& value) {
		S[d] = value;
	}
	static std::vector<std::string> getFieldNames() {
		static std::vector<std::string> const fieldNames = []() {
			std::vector<std::string> names;
			names.push_back("rho");
			names.push_back("eg");
			names.push_back("tau");
			for(int d = 0; d < D; d++) {
				std::string const name = std::string("s_") + std::string(1, 'x' + d);
				names.push_back(name);
			}
			return names;
		}();
		return fieldNames;
	}
private:
	T& rho;
	T& eg;
	T& tau;
	Math::Vector<T, D>& S;
};

template<typename T, int D>
EulerState<T, D> initSodShockTube(Math::Vector<T, D> x) {
	/*********************************************************/
	static constexpr T rhoL = T(1.0);
	static constexpr T pL = T(1.0);
	static constexpr T rhoR = T(0.125);
	static constexpr T pR = T(0.1);
	/*********************************************************/
	static constexpr T c0 = T(1) / (EulerState<T, D>::gamma - T(1));
	static constexpr T c1 = T(1) / EulerState<T, D>::gamma;
	static constexpr T eL = c0 * pL;
	static constexpr T eR = c0 * pR;
	static constexpr T tauL = std::pow(eL, c1);
	static constexpr T tauR = std::pow(eR, c1);
	EulerState<T, D> u;
	u.setMomentum(T(0));
	if (x[0] < T(0.5)) {
		u.setDensity(rhoL);
		u.setEnergy(eL);
		u.setEntropy(tauL);
	} else if (x[0] > T(0.5)) {
		u.setDensity(rhoR);
		u.setEnergy(eR);
		u.setEntropy(tauR);
	} else {
		u.setDensity(0.5 * (rhoL + rhoR));
		u.setEnergy(0.5 * (eL + eR));
		u.setEntropy(0.5 * (tauL + tauR));
	}
	return u;
}

