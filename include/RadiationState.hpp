/******************************************************************************
 Copyright (C) 2024  Dominic C. Marcello
*******************************************************************************/

#ifndef INCLUDE_RADIATIONSTATE_HPP_
#define INCLUDE_RADIATIONSTATE_HPP_

#include "Constants.hpp"
#include "Matrix.hpp"

#include <vector>


template<typename T, int D = 3>
struct RadiationState: public std::array<T, 1 + D> {
	static constexpr int fx_i = 0;
	static constexpr int fy_i = 1;
	static constexpr int fz_i = 2;
	static constexpr int er_i = 3;
	static constexpr int NF = 1 + D;
	using base_type = std::array<T, NF>;
	using value_type = T;
	using eigensys_type = std::pair<std::array<T, NF>, SquareMatrix<T, NF>>;
	static constexpr int dimCount() noexcept {
		return D;
	}
	static constexpr int fieldCount() noexcept {
		return NF;
	}
	RadiationState() {
	}
	RadiationState(base_type const &other) {
		((base_type&) (*this)).operator=(other);
	}
	RadiationState(RadiationState const &other) {
		*this = other;
	}
	RadiationState(T const &other) {
		for (auto &u : *this) {
			u = other;
		}
	}
	RadiationState(RadiationState &&other) {
		*this = std::move(other);
	}
	RadiationState& operator=(RadiationState const &other) {
		static_cast<base_type&>(*this) = static_cast<base_type const&>(other);
		return *this;
	}
	RadiationState& operator=(RadiationState &&other) {
		static_cast<base_type&>(*this) = std::move(static_cast<base_type&&>(other));
		return *this;
	}
	std::array<T, NF> eigenvalues(int dim) const {
		std::array<T, NF> lambda { };
		auto const cons = getCodeConstants();
		using std::sqrt;
		auto const &U = *this;
		T const F = sqrt(sqr(U[fx_i]) + sqr(U[fy_i]) + sqr(U[fz_i]));
		T const f = F / U[er_i];
		if ((1 > f) && (f > 0)) {
			T const Finv = T(1) / F;
			T const mu = U[dim] * Finv;
			T const sqrt4m3f2 = sqrt(T(4.0) - T(3.0) * sqr(f));
			T const lambda1 = T(cons.c) * f * mu / sqrt4m3f2;
			T const lambda2 =
					T(cons.c) * sqrt((T(2.0) / T(3.0) * (T(4.0) - T(3.0) * sqr(f) - sqrt4m3f2)) + T(2.0) * (T(2.0) - sqr(f) - sqrt4m3f2) * sqr(mu)) / sqrt4m3f2;
			lambda.front() = lambda1 + lambda2;
			lambda.back() = lambda1 - lambda2;
			for (int i = 1; i < NF - 1; i++) {
				lambda[i] = T(cons.c) * (T(2.0) - sqrt4m3f2) * mu / f;
			}
		} else if (f == 1) {
			T const Finv = T(1) / F;
			T const mu = U[dim] * Finv;
			lambda.front() = T(cons.c) * mu;
			lambda.back() = T(cons.c) * mu;
			for (int i = 1; i < NF - 1; i++) {
				lambda[i] = T(cons.c) * mu;
			}
		} else if (f == 0) {
			lambda.front() = T(cons.c) * sqrt(T(1.0 / 3.0));
			lambda.back() = -T(cons.c) * sqrt(T(1.0 / 3.0));
			for (int i = 1; i < NF - 1; i++) {
				lambda[i] = T(0.0);
			}
		} else {
			printf("\n%e\n", f - 1.0);
			assert(false);
			abort();
		}
		return lambda;

	}
	eigensys_type eigenSystem(int d) const {
		eigensys_type rc;
		auto &lambda = rc.first;
		auto &R = rc.second;
		std::fill(R.begin(), R.end(), T(0));
		lambda = eigenvalues(d);
		R(0, 0) = lambda.front();
		R(D, 0) = T(1);
		R(0, D) = lambda.back();
		R(D, D) = T(1);
		R(1, 1) = T(1);
		R(2, 2) = T(1);
		for (int n = 0; n <= D; n++) {
			std::swap(R(d, n), R(0, n));
		}
		return rc;
	}
	RadiationState flux(int d1) const noexcept {
		auto const &U = *this;
		RadiationState flux;
		T const F2 = sqr(U[fx_i]) + sqr(U[fy_i]) + sqr(U[fz_i]);
		T const F = sqrt(F2);
		if (F) {
			T const F2inv = T(1) / F2;
			T const f = F / U[er_i];
			T const Xi = (T(5.0) - T(2.0) * sqrt(T(4.0) - T(3.0) * sqr(f))) / T(3.0);
			T const cd = (T(1.0) - Xi) / T(2.0);
			T const cs = (T(3.0) * Xi - T(1.0)) / T(2.0);
			flux[er_i] = U[fx_i + d1];
			for (int d2 = 0; d2 < D; d2++) {
				flux[fx_i + d2] = cs * U[fx_i + d1] * U[fx_i + d2] * F2inv * U[er_i];
			}
			flux[fx_i + d1] += cd * U[er_i];
		} else {
			flux.fill(T(0));
			flux[fx_i + d1] = T(1.0 / 3.0) * U[er_i];
		}
		return flux;
	}
	bool sanityCheck() const {
		auto const &U = *this;
		T const F = sqrt(sqr(U[fx_i]) + sqr(U[fy_i]) + sqr(U[fz_i]));
		T const f = F / U[er_i];
		if (f > T(1)) {
			printf("f = %e is gt 1 (%e, %e)\n", f, F, U[er_i]);
			return false;
		}
		return true;
	}
	friend T findPositivityPreservingTheta(const RadiationState &uBar, const RadiationState &uNode) {
		using std::max;
		T constexpr almostOne = T(1.0 - cbrt(std::numeric_limits<double>::epsilon()));
		T const fx0 = uBar[fx_i];
		T const fy0 = uBar[fy_i];
		T const fz0 = uBar[fz_i];
		T const er0 = uBar[er_i];
		T const fx = uNode[fx_i];
		T const fy = uNode[fy_i];
		T const fz = uNode[fz_i];
		T const er = uNode[er_i];
		T theta = T(1);
		if (sqrt(sqr(fx) + sqr(fy) + sqr(fz)) > er) {
			T const num1 = fx0 * (-fx + fx0) + fy0 * (-fy + fy0) - fz * fz0 + sqr(fz0) + (er - er0) * er0;
			T const num2 =
					sqrt(
							-sqr(fy0) * sqr(fz) + T(2.0) * fy * fy0 * fz * fz0 - sqr(fy) * sqr(fz0) + sqr(fy0) * sqr(er) + sqr(fz0) * sqr(er) - sqr(fx0) * (sqr(fy) + sqr(
									fz) - sqr(er)) - T(2.0) * (fy * fy0 + fz * fz0) * er * er0 + (sqr(fy) + sqr(fz)) * sqr(er0) + T(2.0) * fx * fx0 * (fy * fy0 + fz * fz0 - er * er0) - sqr(
									fx) * (sqr(fy0) + sqr(fz0) - sqr(er0)));
			T const den_inv = T(1) / (sqr(fx - fx0) + sqr(fy - fy0) + (fz - fz0 + er - er0) * (fz - fz0 - er + er0));
			T const thetaP = (num1 + num2) * den_inv;
			T const thetaM = (num1 - num2) * den_inv;
			bool const validP = (T(1) > thetaP) && (thetaP >= T(0));
			bool const validM = (T(1) > thetaM) && (thetaM >= T(0));
			if (validP) {
				if (validM) {
					theta = (thetaP > thetaM) ? thetaP : thetaM;
				} else {
					theta = thetaP;
				}
			} else {
				theta = thetaM;
			}
			theta *= almostOne;
		}
		return theta;
	}
	friend RadiationState solveRiemannProblem(const RadiationState &uL, const RadiationState &uR, int dim) noexcept {
		RadiationState fR, fL;
		fL = uL.flux(dim);
		fR = uR.flux(dim);
		auto lambdaL = uL.eigenvalues(dim);
		auto lambdaR = uR.eigenvalues(dim);
		T sR = T(0);
		T sL = T(0);
		for (int l = 0; l < NF; l++) {
			sR = std::max(sR, lambdaR[l]);
			sR = std::max(sR, lambdaL[l]);
			sL = std::min(sL, lambdaR[l]);
			sL = std::min(sL, lambdaL[l]);
		}
		return (sR * fL - sL * fR + sR * sL * (uR - uL)) / (sR - sL);
	}
	constexpr T const& getEnergy() const {
		return (*this)[D];
	}
	constexpr T const& getFlux(int d) const {
		return (*this)[d];
	}
	constexpr void setEnergy(T const &value) {
		(*this)[D] = value;
	}
	constexpr void setFlux(int d, T const &value) {
		(*this)[d] = value;
	}
	static std::vector<std::string> getFieldNames() {
		static std::vector<std::string> const fieldNames = []() {
			std::vector<std::string> names;
			for (int d = 0; d < D; d++) {
				std::string const name = std::string("F") + std::string(1, 'x' + d);
				names.push_back(name);
			}
			names.push_back("E");
			return names;
		}();
		return fieldNames;
	}
};


#endif /* INCLUDE_RADIATIONSTATE_HPP_ */
