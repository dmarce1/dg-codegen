#pragma once

#include "ContainerArithmetic.hpp"
#include "Matrix.hpp"
#include "Util.hpp"

#include <array>
#include <cmath>
#include <vector>
#include <string>

#define EULERS_CONSTRUCTION : \
	rho(*std::launder(reinterpret_cast<T*>(base_type::data() + 0))), \
	eg (*std::launder(reinterpret_cast<T*>(base_type::data() + 1))), \
	S  (*std::launder(reinterpret_cast<std::array<T, D>*>(base_type::data() + 2)))

enum class RiemannSolver : int {
	LLF, HLL, HLLC
};

template<typename T, int D, RiemannSolver riemannSolver = RiemannSolver::LLF>
struct EulerState: public std::array<T, 2 + D> {
	using EleType = typename ElementType<T>::type;
	static constexpr int NF = 2 + D;
	static constexpr EleType zero = EleType(0);
	static constexpr EleType half = EleType(1) / EleType(2);
	static constexpr EleType one = EleType(1);
	static constexpr EleType gamma = EleType(5) / EleType(3);
	static constexpr EleType gamm1 = gamma - one;
	static constexpr EleType igamm1 = one / gamm1;
	using base_type = std::array<T, NF>;
	using value_type = EleType;
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
		lambda[0] -= a;
		lambda[2] += a;
		return lambda;

	}
	eigensys_type eigenSystem(int dim) const {
		eigensys_type rc;
		auto& eigenvalues = rc.first;
		auto& rightEigenvectors = rc.second;
		T const invρ = one / rho;
		std::array<T, D> u = S * invρ;
		if(dim) {
			std::swap(u[dim], u[0]);
		}
		T const ek = half * dot(u, u);
		T const ei = max(zero, eg - rho * ek);
		T const p = gamm1 * ei;
		T const a = sqrt(gamma * p * invρ);
		T const h0 = (p + eg) * invρ;
		std::fill(eigenvalues.begin(), eigenvalues.end(), u[0]);
		eigenvalues[0] -= a;
		eigenvalues[2] += a;
		rightEigenvectors(0, 0) = one;
		rightEigenvectors(0, 1) = one;
		rightEigenvectors(0, 2) = one;
		rightEigenvectors(1, 0) = u[0] - a;
		rightEigenvectors(1, 1) = u[0];
		rightEigenvectors(1, 2) = u[0] + a;
		for(int r = 2; r < NF - 1; r++) {
			rightEigenvectors(r, 0) = u[r - 1];
			rightEigenvectors(r, 1) = u[r - 1];
			rightEigenvectors(r, 2) = u[r - 1];
		}
		rightEigenvectors(NF - 1, 0) = h0 - a * u[0];
		rightEigenvectors(NF - 1, 1) = ek;
		rightEigenvectors(NF - 1, 2) = h0 + a * u[0];
		for(int c = 3; c < NF; c++) {
			rightEigenvectors(0, c) = zero;
			for(int r = 1; r < NF - 1; r++) {
				rightEigenvectors(r, c) = (r + 1 == c) ? one : zero;
			}
			rightEigenvectors(NF - 1, c) = u[c - 2];
		}
		if( dim ) {
			for(int c = 0; c < NF; c++) {
				std::swap(rightEigenvectors(1, c), rightEigenvectors(dim + 1, c));
			}
		}
		return rc;
	}
	EulerState flux(int d1) const noexcept {
		EulerState F;
		T const irho = one / rho;
		auto const v = S * irho;
		T const ek = half * dot(v, S);
		T const ei = max(zero, eg - ek);
		T const p = (gamma - one) * ei;
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
		using std::sqrt;
		using std::abs;
		if constexpr(riemannSolver == RiemannSolver::LLF) {
			constexpr int N2 = 2;
			constexpr int L = 0;
			constexpr int R = 1;
			std::array<EulerState, N2> u = {uL, uR};
			std::array<EulerState, N2> f;
			std::array<T, N2> irho, v, ek, ei, a;
			std::array<T, N2> p;
			for(int i = 0; i < N2; i++) {
				irho[i] = one / u[i].rho;
				v[i] = irho[i] * u[i].S[dim];
				ek[i] = zero;
				for(int d = 0; d < D; d++) {
					ek[i] += half * irho[i] * sqr(u[i].S[d]);
				}
				ei[i] = max(zero, u[i].eg - ek[i]);
				p[i] = gamm1 * ei[i];
				a[i] = sqrt(gamma * p[i] * irho[i]);
			}
			T const s = max(T(a[L] + abs(v[L])), T(a[R] + abs(v[R])));
			for(int i = 0; i < N2; i++) {
				f[i] = u[i].flux(dim);
			}
			EulerState flux;
			for(int fi = 0; fi < NF; fi++) {
				flux[fi] = (f[L][fi] + f[R][fi] - s * (uR[fi] - uL[fi])) * half;
			}
			return flux;
		} else if constexpr(riemannSolver == RiemannSolver::HLL) {
			constexpr int N2 = 2;
			constexpr int L = 0;
			constexpr int R = 1;
			std::array<EulerState, N2> u = {uL, uR};
			std::array<EulerState, N2> f;
			std::array<T, N2> irho, v, ek, ei, a;
			std::array<T, N2> s, p;
			for(int i = 0; i < N2; i++) {
				irho[i] = one / u[i].rho;
				v[i] = irho[i] * u[i].S[dim];
				ek[i] = zero;
				for(int d = 0; d < D; d++) {
					ek[i] += half * irho[i] * sqr(u[i].S[d]);
				}
				ei[i] = std::max(zero, u[i].eg - ek[i]);
				p[i] = (gamma - one) * ei[i];
				a[i] = sqrt(gamma * p[i] * irho[i]);
			}
			s[L] = min(v[L] - a[L], v[R] - a[R]);
			s[R] = max(v[R] + a[R], v[L] + a[L]);
			s[L] = min(s[L], zero);
			s[R] = max(s[R], zero);
			for(int i = 0; i < N2; i++) {
				f[i] = u[i].flux(dim);
			}
			return (s[R] * f[L] - s[L] * f[R] + s[L] * s[R] *(uR - uL)) / (s[R] - s[L]);
		} else/*if constexpr(riemannSolver == RiemannSolver::HLLC)*/{
			constexpr int N2 = 2;
			constexpr int N3 = 3;
			constexpr int L = 0;
			constexpr int R = 1;
			constexpr int STAR = 2;
			std::array<EulerState, N3> u = {uL, uR};
			std::array<EulerState, N2> f;
			std::array<T, N2> irho, v, ek, ei, a;
			std::array<T, N3> s, p;
			for(int i = 0; i < N2; i++) {
				irho[i] = one / u[i].rho;
				v[i] = irho[i] * u[i].S[dim];
				ek[i] = zero;
				for(int d = 0; d < D; d++) {
					ek[i] += half * irho[i] * sqr(u[i].S[d]);
				}
				ei[i] = std::max(zero, u[i].eg - ek[i]);
				p[i] = (gamma - one) * ei[i];
				a[i] = sqrt(gamma * p[i] * irho[i]);
			}
			s[L] = min(v[L] - a[L], v[R] - a[R]);
			s[R] = max(v[L] + a[L], v[R] + a[R]);
			T const num = p[R] - p[L] + u[L].rho * v[L] * (s[L] - v[L]) - u[R].rho * v[R] * (s[R] - v[R]);
			T const den = u[L].rho * (s[L] - v[L]) - u[R].rho * (s[R] - v[R]);
			s[STAR] = num / den;
			for(int i = 0; i < N2; i++) {
				f[i] = u[i].flux(dim);
			}
			if (zero < s[L]) {
				return f[L];
			} else if (zero > s[R]) {
				return f[R];
			} else {
				int const i = ((s[STAR] > zero) ? L : R);
				u[STAR].rho = u[i].rho * (s[i] - v[i]) / (s[i] - s[STAR]);
				for(int dir = 0; dir < D; dir++) {
					if(dim != dir ) {
						u[STAR].S[dir] = u[STAR].rho * irho[i] * u[i].S[dir];
					} else {
						u[STAR].S[dir] = u[STAR].rho * s[STAR];
					}
				}
				p[STAR] = p[i] + u[i].rho * (s[i] - v[i]) * (s[STAR] - v[i]);
				u[STAR].eg = ((s[i] - v[i]) * u[i].eg - p[i] * v[i] + p[STAR] * s[STAR]) / (s[i] - s[STAR]);
				return f[i] + s[i] * (u[STAR] - u[i]);
			}
		}
	}
	bool sanityCheck() const {
		return true;
	}
	friend T findPositivityPreservingTheta(const EulerState &uBar, const EulerState &uNode) noexcept {
		return one;
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
			names.push_back("rho");
			names.push_back("eg");
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
	std::array<T, D>& S;
};

template<typename T, int D>
struct CanDoArithmetic<EulerState<T, D>> {
	static constexpr bool value = true;
};

template<typename T, int D>
EulerState<T, D> initSodShockTube(std::array<T, D> x) {
	using EleType = typename ElementType<T>::type;
	static constexpr EleType half = EleType(1) / EleType(2);
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
	u.setMomentum(zeroArray<T, D>());
	if (x[0] < half) {
		u.setDensity(rhoL);
		u.setEnergy(eL);
	} else if (x[0] > half) {
		u.setDensity(rhoR);
		u.setEnergy(eR);
	} else {
		u.setDensity(half * (rhoL + rhoR));
		u.setEnergy(half * (eL + eR));
	}
	return u;
}

