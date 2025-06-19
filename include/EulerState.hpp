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

template<typename T, int D, T G = T(5.0 / 3.0)>
struct EulerState: public std::array<T, 2 + D> {
	static constexpr int NF = 2 + D;
	using base_type = std::array<T, NF>;
	static constexpr T gamma = G;
	using value_type = T;
	using eigensys_type = std::pair<std::array<T, NF>, SquareMatrix<T, NF>>;
	static constexpr int dimCount() noexcept {
		return D;
	}
	static constexpr int fieldCount() noexcept {
		return NF;
	}
	constexpr EulerState() EULERS_CONSTRUCTION {
	}
	constexpr EulerState(base_type const &other) EULERS_CONSTRUCTION {
		((base_type&)(*this)).operator=(other);
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
	constexpr std::array<T, NF> eigenvalues(int dim) const {
		std::array<T, NF> lambda;
		using std::sqrt;
		T const irho = T(1) / rho;
		T const v = S[dim] * irho;
		T const ek = T(0.5) * irho * dot(S, S);
		T const ei = std::max(eg - ek, T(0));
		T const p = (gamma - T(1)) * ei;
		T const a = sqrt(gamma * p * irho);
		std::fill(lambda.begin(), lambda.end(), v);
		lambda.front() -= a;
		lambda.back() += a;
		return lambda;

	}
	constexpr eigensys_type eigenSystem(int dim) const {
		eigensys_type rc;
		using std::sqrt;

		auto& λ = rc.first;     // eigenvalues
		auto& R = rc.second;// right‐eigenvector matrix (NF×NF)
		std::fill(R.begin(), R.end(), T(0));

		// primitive variables
		T const invρ = T(1)/rho;
		auto const v = S * invρ;// velocity vector
		T ke = T(0);
		for(int d = 0; d < D; ++d)
		ke += T(0.5)*v[d]*S[d];
		T const ei = std::max(T(0), eg - ke);
		T const p = (gamma - T(1))*ei;
		T const c = sqrt(gamma*p*invρ);
		T const h = (eg + p)*invρ;
		T const u_n = v[dim];// normal velocity

		// fill eigenvalues: [u-c, u, …, u, u+c]
		for(int i = 0; i < NF; ++i) λ[i] = u_n;
		λ[0] = u_n - c;
		λ[NF-1] = u_n + c;

		int col = 0;

		// 1) slow acoustic wave (u-c)
		R(0, col) = T(1);
		for(int i = 0; i < D; ++i) {
			R(1 + i, col) = (i == dim ? u_n - c : v[i]);
		}
		R(D+1, col) = h - u_n*c;
		++col;

		// 2) D−1 shear waves (speed = u_n), each with jump only in one transverse momentum
		for(int i = 0; i < D; ++i) {
			if (i == dim) continue;
			R(0, col) = T(0);         // no density jump
			R(1 + i, col) = T(1);// δS_i = 1
			R(D+1, col) = T(0);// no energy jump
			++col;
		}

		// 3) entropy/contact wave (speed = u_n)
		R(0, col) = T(1);
		for(int i=0; i<D; ++i)
		R(1+i, col) = v[i];
		R(D+1, col) = T(0);// no pressure/energy jump
		++col;

		// 4) fast acoustic wave (u+c)
		R(0, col) = T(1);
		for(int i=0; i<D; ++i) {
			R(1 + i, col) = (i == dim ? u_n + c : v[i]);
		}
		R(D+1, col) = h + u_n*c;

		return rc;
	}
	constexpr EulerState flux(int d) const noexcept {
		EulerState F;
		T const irho = T(1) / rho;
		auto const v = S * irho;
		T const ek = T(0.5) * dot(v, S);
		T const ei = std::max(T(0), eg - ek);
		T const p = (gamma - T(1)) * ei;
		T const u = v[d];
		F.rho = S[d];
		F.eg = u * (eg + p);
		F.S = u * S;
		F.S[d] += p;
		return F;
	}
	friend constexpr EulerState solveRiemannProblem(const EulerState &uL, const EulerState &uR, int dim) noexcept {
		using std::sqrt;
		using state_t = EulerState;
		constexpr int N2 = 2;
		constexpr int N3 = 3;
		constexpr int L = 0;
		constexpr int R = 1;
		constexpr int STAR = 2;
		std::array<state_t, N3> u = {uL, uR};
		std::array<state_t, N2> f;
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
		s[L] = std::min(v[L] - a[L], v[R] - a[R]);
		s[R] = std::max(v[L] + a[L], v[R] + a[R]);
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
	bool sanityCheck() const {
		return true;
	}
	friend T findPositivityPreservingTheta(const EulerState &uBar, const EulerState &uNode) noexcept {
		return T(1.0);
	}
	constexpr T const& getDensity() const {
		return rho;
	}
	constexpr T const& getEnergy() const {
		return eg;
	}
	constexpr std::array<T, D> const& getMomentum() const {
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
	constexpr void setMomentum(std::array<T, D> const& value) {
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

