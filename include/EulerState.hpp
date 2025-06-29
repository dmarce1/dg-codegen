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
	static constexpr int NF = 2 + D;
	static constexpr double gamma = 5.0 / 3.0;
	static constexpr double gamm1 = 2.0 / 3.0;
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
		T const irho = T(1) / rho;
		T const v = S[dim] * irho;
		T const ek = T(0.5) * irho * dot(S, S);
		T const ei = max(eg - ek, 0.0);
		T const p = gamm1 * ei;
		T const a = sqrt(gamma * p * irho);
		std::fill(lambda.begin(), lambda.end(), v);
		lambda.front() -= a;
		lambda.back() += a;
		return lambda;

	}
	eigensys_type eigenSystem(int dim) const {
		eigensys_type rc;
		using std::sqrt;
		auto& eigenvalues = rc.first;
		auto& rightEigenvectors = rc.second;
		int column = 0;
		std::fill(rightEigenvectors.begin(), rightEigenvectors.end(), T(0));
		T const invρ = T(1) / rho;
		auto const v = S * invρ;
		T ke = T(0);
		for(int d = 0; d < D; d++) {
			ke += T(0.5) * v[d] * S[d];
		}
		T const ei = max(T(0), eg - ke);
		T const p = gamm1 * ei;
		T const c = sqrt(gamma * p * invρ);
		T const h = (eg + p)*invρ;
		T const u = v[dim];
		for(int i = 0; i < NF; i++) {
			eigenvalues[i] = u;
		}
		eigenvalues.front() -= c;
		eigenvalues.back() += c;
		rightEigenvectors(0, column) = T(1);
		for(int thisDimension = 0; thisDimension < D; thisDimension++) {
			rightEigenvectors(1 + thisDimension, column) = ((thisDimension == dim) ? (u - c) : v[thisDimension]);
		}
		rightEigenvectors(D + 1, column) = h - u * c;
		column++;
		for(int thisDimension = 0; thisDimension < D; thisDimension++) {
			if (thisDimension == dim) {
				continue;
			}
			rightEigenvectors(0, column) = T(0);
			rightEigenvectors(1 + thisDimension, column) = T(1);
			rightEigenvectors(D + 1, column) = T(0);
			column++;
		}
		rightEigenvectors(0, column) = T(1);
		for(int thisDimension = 0; thisDimension < D; thisDimension++) {
			rightEigenvectors(1 + thisDimension, column) = v[thisDimension];
		}
		rightEigenvectors(D + 1, column) = T(0);
		column++;
		rightEigenvectors(0, column) = T(1);
		for(int thisDimension = 0; thisDimension < D; thisDimension++) {
			rightEigenvectors(1 + thisDimension, column) = ((thisDimension == dim) ? (u + c) : v[thisDimension]);
		}
		rightEigenvectors(D + 1, column) = h + u * c;
		return rc;
	}
	EulerState flux(int d1) const noexcept {
		EulerState F;
		T const irho = T(1) / rho;
		auto const v = S * irho;
		T const ek = T(0.5) * dot(v, S);
		T const ei = max(T(0), eg - ek);
		T const p = (gamma - T(1)) * ei;
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
				irho[i] = T(1) / u[i].rho;
				v[i] = irho[i] * u[i].S[dim];
				ek[i] = T(0);
				for(int d = 0; d < D; d++) {
					ek[i] += T(0.5) * irho[i] * sqr(u[i].S[d]);
				}
				ei[i] = max(T(0), u[i].eg - ek[i]);
				p[i] = (gamma - T(1)) * ei[i];
				a[i] = sqrt(gamma * p[i] * irho[i]);
			}
			T const s =  max(T(a[L] + abs(v[L])), T(a[R] + abs(v[R])));
			for(int i = 0; i < N2; i++) {
				f[i] = u[i].flux(dim);
			}
			EulerState flux;
			for(int fi = 0; fi < NF; fi++) {
				flux[fi] = (f[L][fi] + f[R][fi] - s * (uR[fi] - uL[fi])) * 0.5;
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
				irho[i] = T(1) / u[i].rho;
				v[i] = irho[i] * u[i].S[dim];
				ek[i] = T(0);
				for(int d = 0; d < D; d++) {
					ek[i] += T(0.5) * irho[i] * sqr(u[i].S[d]);
				}
				ei[i] = std::max(T(0), u[i].eg - ek[i]);
				p[i] = (gamma - T(1)) * ei[i];
				a[i] = sqrt(gamma * p[i] * irho[i]);
			}
			s[L] = min(v[L] - a[L], v[R] - a[R]);
			s[R] = max(v[R] + a[R], v[L] + a[L]);
			s[L] = min(s[L], T(0));
			s[R] = max(s[R], T(0));
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
				irho[i] = T(1) / u[i].rho;
				v[i] = irho[i] * u[i].S[dim];
				ek[i] = T(0);
				for(int d = 0; d < D; d++) {
					ek[i] += T(0.5) * irho[i] * sqr(u[i].S[d]);
				}
				ei[i] = std::max(T(0), u[i].eg - ek[i]);
				p[i] = (gamma - T(1)) * ei[i];
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
			if (T(0) < s[L]) {
				return f[L];
			} else if (T(0) > s[R]) {
				return f[R];
			} else {
				int const i = ((s[STAR] > T(0)) ? L : R);
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
		return T(1.0);
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
	/*********************************************************/
	static constexpr T rhoL = T(1.0);
	static constexpr T pL = T(1.0);
	static constexpr T rhoR = T(0.125);
	static constexpr T pR = T(0.1);
	/*********************************************************/
	static constexpr T c0 = T(1) / (EulerState<T, D>::gamma - T(1));
	static constexpr T eL = c0 * pL;
	static constexpr T eR = c0 * pR;
	EulerState<T, D> u;
	u.setMomentum(zero<T, D>());
	if (x[0] < T(0.5)) {
		u.setDensity(rhoL);
		u.setEnergy(eL);
	} else if (x[0] > T(0.5)) {
		u.setDensity(rhoR);
		u.setEnergy(eR);
	} else {
		u.setDensity(T(0.5) * (rhoL + rhoR));
		u.setEnergy(T(0.5) * (eL + eR));
	}
	return u;
}

