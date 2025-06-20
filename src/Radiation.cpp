#include "Constants.hpp"
#include "EulerState.hpp"
#include "HyperGrid.hpp"
#include "Real.hpp"
#include "RungeKutta.hpp"
#include <array>
#include <iostream>
using T = double;

static const PhysicalUnits<T> physicalUnits(T(1), T(1), T(1));
static const auto cons = physicalUnits.getPhysicalConstants();

constexpr int NDIM = 3;

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
		for (int n = 0; n <= NDIM; n++) {
			std::swap(R(n, d), R(n, 0));
		}
		for (int n = 0; n <= NDIM; n++) {
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
		return (*this)[NDIM];
	}
	constexpr T const& getFlux(int d) const {
		return (*this)[d];
	}
	constexpr void setEnergy(T const &value) {
		(*this)[NDIM] = value;
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

template<typename T, int D>
struct CanDoArithmetic<RadiationState<T, D>> {
	static constexpr bool value = true;
};

template<typename T, int D>
RadiationState<T, D> initRadiationFront(std::array<T, D> x) {
	/*********************************************************/
	RadiationState<T, D> u;
	T Tr = T(5000.0);
	T Er = T(cons.a) * sqr(sqr(Tr));
	u.setFlux(0, T(0));
	u.setFlux(1, T(0));
	u.setFlux(2, T(0));
	u.setEnergy(Er);
	if (T(0.75) > x[0] && x[0] > T(0.25)) {
		u.setFlux(0, Er * T(0.9));
	} else {
		u.setEnergy(Er / T(1000));
	}
	return u;
}

//auto computeG(T dEr, std::array<T, NDIM> dF, T Er0, std::array<T, NDIM> F0, T Eg0, std::array<T, NDIM> Beta0, T rho, T mu, T kappa, T chi, T gamma, T dt) {
//	constexpr T tiny = 1e-50;
//	std::array<T, NDIM> Beta;
//	for (int n = 0; n < NDIM; n++) {
//		Beta[n] = Beta0[n] - dF[n] / (rho * cons.c * cons.c);
//	}
//	T Ek = 0;
//	std::array<T, NDIM> dEk_dF;
//	for (int n = 0; n < NDIM; n++) {
//		Ek += 0.5 * rho * sqr(cons.c * Beta[n]);
//		dEk_dF[n] = -Beta[n];
//	}
//	T const Eg = Eg0 - dEr;
//	T const eps = Eg - Ek;
//	std::array<T, NDIM> const deps_dF = -dEk_dF;
//	T F2 = 0;
//	for (int n = 0; n < NDIM; n++) {
//		F2 += sqr(F0[n] + dF[n]);
//	}
//	T const absF = sqrt(F2);
//	T const iCv = (mu * cons.amu) * (gamma - 1.0) / (cons.kB * rho);
//	T const temp = eps * iCv;
//	T const dtemp_dEr = -iCv;
//	std::array<T, NDIM> const dtemp_dF = deps_dF * iCv;
//	T const temp2 = sqr(temp);
//	T const temp4 = sqr(temp2);
//	std::array<std::array<T, NDIM>, NDIM> dBeta_dF;
//	for (int n = 0; n < NDIM; n++) {
//		for (int m = 0; m < NDIM; m++) {
//			dBeta_dF[n][m] = -T(n == m) / (rho * cons.c * cons.c);
//		}
//	}
//	std::array<T, NDIM> N;
//	for (int n = 0; n < NDIM; n++) {
//		N[n] = (F0[n] + dF[n]) / (absF + tiny);
//	}
//	T const temp3 = temp * temp2;
//	std::array<std::array<T, NDIM>, NDIM> dN_dF;
//	for (int n = 0; n < NDIM; n++) {
//		for (int m = 0; m < NDIM; m++) {
//			dN_dF[n][m] = (T(n == m) - N[n] * N[m]) / (absF + tiny);
//		}
//	}
//	std::array<T, NDIM> df_dF;
//	T const f = absF / (Er0 + dEr);
//	T const df_dEr = -f / (Er0 + dEr);
//	for (int n = 0; n < NDIM; n++) {
//		df_dF[n] = N[n] / (Er0 + dEr);
//	}
//	T const Xi = (5 - 2 * sqrt(4 - 3 * sqr(f))) / 3;
//	T const dXi_df = 2 * f / sqrt(4 - 3 * sqr(f));
//	T const dXi_dEr = dXi_df * df_dEr;
//	std::array<T, NDIM> dXi_dF;
//	for (int n = 0; n < NDIM; n++) {
//		dXi_dF[n] = dXi_df * df_dF[n];
//	}
//	std::array<std::array<T, NDIM>, NDIM> D;
//	for (int n = 0; n < NDIM; n++) {
//		for (int m = 0; m < NDIM; m++) {
//			D[n][m] = ((1 - Xi) / 2 * T(n == m) + (3 * Xi - 1) / 2 * N[n] * N[m]);
//		}
//	}
//	std::array<std::array<T, NDIM>, NDIM> dD_dEr;
//	for (int n = 0; n < NDIM; n++) {
//		for (int m = 0; m < NDIM; m++) {
//			dD_dEr[n][m] = -dXi_dEr / 2 * T(n == m);
//			dD_dEr[n][m] += 3 * dXi_dEr / 2 * N[n] * N[m];
//		}
//	}
//	std::array<std::array<std::array<T, NDIM>, NDIM>, NDIM> dD_dF;
//	for (int n = 0; n < NDIM; n++) {
//		for (int m = 0; m < NDIM; m++) {
//			for (int l = 0; l < NDIM; l++) {
//				dD_dF[n][m][l] = -dXi_dF[l] / 2 * T(n == m);
//				dD_dF[n][m][l] += 3 * dXi_dF[l] / 2 * N[n] * N[m];
//				dD_dF[n][m][l] += (3 * Xi - 1) / 2 * (dN_dF[n][l] * N[m] + dN_dF[m][l] * N[n]);
//			}
//		}
//	}
//	std::array<std::array<T, NDIM>, NDIM> P;
//	for (int n = 0; n < NDIM; n++) {
//		for (int m = 0; m < NDIM; m++) {
//			P[n][m] = (Er0 + dEr) * D[n][m];
//		}
//	}
//	std::array<std::array<T, NDIM>, NDIM> dP_dEr;
//	for (int n = 0; n < NDIM; n++) {
//		for (int m = 0; m < NDIM; m++) {
//			dP_dEr[n][m] = (Er0 + dEr) * dD_dEr[n][m] + D[n][m];
//		}
//	}
//	std::array<std::array<std::array<T, NDIM>, NDIM>, NDIM> dP_dF;
//	for (int n = 0; n < NDIM; n++) {
//		for (int m = 0; m < NDIM; m++) {
//			for (int l = 0; l < NDIM; l++) {
//				dP_dF[l][n][m] = (Er0 + dEr) * dD_dF[l][n][m];
//			}
//		}
//	}
//	T gk = kappa * (Er0 + dEr - cons.a * temp4);
//	for (int k = 0; k < NDIM; k++) {
//		gk -= 2 * kappa * Beta[k] * (F0[k] + dF[k]);
//	}
//	T const dgk_dEr = kappa * (1 - 4 * cons.a * temp3 * dtemp_dEr);
//	std::array<T, NDIM> dgk_dF;
//	for (int n = 0; n < NDIM; n++) {
//		dgk_dF[n] = -kappa * (2 * Beta[n] + 4 * cons.a * temp3 * dtemp_dF[n]);
//		for (int k = 0; k < NDIM; k++) {
//			dgk_dF[n] -= 2 * kappa * (F0[k] + dF[k]) * dBeta_dF[k][n];
//		}
//	}
//	std::array<T, NDIM> Gx { };
//	std::array<T, NDIM> dGx_dEr { };
//	std::array<std::array<T, NDIM>, NDIM> dGx_dF { };
//	for (int n = 0; n < NDIM; n++) {
//		Gx[n] = chi * (F0[n] + dF[n] - Beta[n] * (Er0 + dEr));
//		for (int k = 0; k < NDIM; k++) {
//			Gx[n] -= chi * P[k][n] * Beta[k];
//		}
//	}
//	for (int n = 0; n < NDIM; n++) {
//		dGx_dEr[n] = -chi * Beta[n];
//		for (int k = 0; k < NDIM; k++) {
//			dGx_dEr[n] -= chi * dP_dEr[k][n] * Beta[k];
//		}
//	}
//	for (int n = 0; n < NDIM; n++) {
//		for (int m = 0; m < NDIM; m++) {
//			dGx_dF[n][m] = chi * (T(n == m) - dBeta_dF[n][m] * (Er0 + dEr));
//			for (int k = 0; k < NDIM; k++) {
//				dGx_dF[n][m] -= chi * dP_dF[k][n][m] * Beta[k];
//				dGx_dF[n][m] -= chi * P[k][n] * dBeta_dF[k][m];
//			}
//		}
//	}
//	T gx = 0;
//	T dgx_dEr = 0;
//	std::array<T, NDIM> dgx_dF { };
//	for (int k = 0; k < NDIM; k++) {
//		gx += Gx[k] * Beta[k];
//	}
//	for (int k = 0; k < NDIM; k++) {
//		dgx_dEr += dGx_dEr[k] * Beta[k];
//	}
//	for (int n = 0; n < NDIM; n++) {
//		dgx_dF[n] = 0;
//		for (int k = 0; k < NDIM; k++) {
//			dgx_dF[n] += dGx_dF[k][n] * Beta[k] + dBeta_dF[k][n] * Gx[k];
//		}
//	}
//	std::array<T, NDIM> Gk { };
//	for (int n = 0; n < NDIM; n++) {
//		Gk[n] = gk * Beta[n];
//	}
//	std::array<T, NDIM> dGk_dEr;
//	for (int n = 0; n < NDIM; n++) {
//		dGk_dEr[n] = dgk_dEr * Beta[n];
//	}
//	std::array<std::array<T, NDIM>, NDIM> dGk_dF;
//	for (int n = 0; n < NDIM; n++) {
//		for (int m = 0; m < NDIM; m++) {
//			dGk_dF[n][m] = Beta[n] * dgk_dF[m] + gk * dBeta_dF[n][m];
//		}
//	}
//	T const hr = dEr + dt * cons.c * (gk + gx);
//	T const dhr_dEr = 1 + dt * cons.c * (dgk_dEr + dgx_dEr);
//	std::array<T, NDIM> dhr_dF;
//	for (int n = 0; n < NDIM; n++) {
//		dhr_dF[n] = dt * cons.c * (dgk_dF[n] + dgx_dF[n]);
//	}
//	std::array<T, NDIM> Hr;
//	for (int n = 0; n < NDIM; n++) {
//		Hr[n] = dF[n] + dt * cons.c * (Gk[n] + Gx[n]);
//	}
//	std::array<T, NDIM> dHr_dEr;
//	for (int n = 0; n < NDIM; n++) {
//		dHr_dEr[n] = dt * cons.c * (dGk_dEr[n] + dGx_dEr[n]);
//	}
//	std::array<std::array<T, NDIM>, NDIM> dHr_dF;
//	for (int n = 0; n < NDIM; n++) {
//		for (int m = 0; m < NDIM; m++) {
//			dHr_dF[n][m] = T(n == m) + dt * cons.c * (dGk_dF[n][m] + dGx_dF[n][m]);
//		}
//	}
//
//	std::pair<std::array<T, NDIM + 1>, std::array<std::array<T, NDIM + 1>, NDIM + 1>> rc;
//	std::array<T, NDIM + 1> &F4 = rc.first;
//	std::array<std::array<T, NDIM + 1>, NDIM + 1> &dF4 = rc.second;
//	for (int k = 0; k < NDIM; k++) {
//		F4[k] = Hr[k];
//		for (int n = 0; n < NDIM; n++) {
//			dF4[k][n] = dHr_dF[k][n];
//		}
//		dF4[NDIM][k] = dhr_dF[k];
//		dF4[k][NDIM] = dHr_dEr[k];
//	}
//	F4[NDIM] = hr;
//	dF4[NDIM][NDIM] = dhr_dEr;
//	return rc;
//}

#include <random>

double randomLogNormal(double mean = 1.0, double stddev = 1.0) {
	static thread_local std::mt19937_64 rng(42);
	std::lognormal_distribution<double> dist(mean, stddev);
	return dist(rng);
}

double random(double a, double b) {
	static thread_local std::mt19937_64 rng(42);
	double const lower = std::nextafter(a, b);
	double const upper = std::nextafter(b, a);
	std::uniform_real_distribution<double> dist(lower, upper);
	return dist(rng);
}

double rand1() {
	return random(0.0, 1.0);
}

double randomLog(double a, double b) {
	static thread_local std::mt19937_64 rng(42);
	double const logA = std::log(std::nextafter(a, b));
	double const logB = std::log(std::nextafter(b, a));
	std::uniform_real_distribution<double> dist(logA, logB);
	return std::exp(dist(rng));
}

double randomNormal(double mean = 1.0, double stddev = 1.0) {
	static thread_local std::mt19937_64 rng(42);
	std::normal_distribution<double> dist(mean, stddev);
	return dist(rng);
}

//Marsaglia (1972)
std::array<double, 3> randomUnitVector3D() {
	static thread_local std::mt19937_64 rng(42);
	static thread_local std::uniform_real_distribution<double> dist(-1.0, 1.0);
	double x1, x2, s;
	do {
		x1 = dist(rng);
		x2 = dist(rng);
		s = x1 * x1 + x2 * x2;
	} while (s >= 1.0 || s == 0.0);
	double factor = std::sqrt(1 - s);
	return {2 * x1 * factor, 2 * x2 * factor, 1 - 2 * s};
}

template<typename T>
std::array<T, 3> crossProduct(std::array<T, 3> const &a, std::array<T, 3> const &b) {
	std::array<T, 3> c;
	c[0] = a[1] * b[2] - a[2] * b[1];
	c[1] = a[2] * b[0] - a[0] * b[2];
	c[2] = a[0] * b[1] - a[1] * b[0];
	return c;
}

template<typename T, auto D>
T dot(std::array<T, D> const &a, std::array<T, D> const &b) {
	T sum = T(0);
	for (int i = 0; i < int(D); i++) {
		sum += a[i] * b[i];
	}
	return sum;
}

template<typename T, auto D>
T norm(std::array<T, D> const &a) {
	return sqrt(dot(a, a));
}

template<typename T, auto D>
std::array<T, D> normalize(std::array<T, D> const &a) {
	std::array<T, D> n = a;
	n *= T(1) / norm(a);
	return n;
}

template<typename T>
SquareMatrix<T, 3> crossProductMatrix(std::array<T, 3> const &a) {
	SquareMatrix<T, 3> R;
	R(0, 0) = R(1, 1) = R(2, 2) = T(0);
	R(0, 1) = -a[2];
	R(0, 2) = +a[1];
	R(1, 2) = -a[0];
	R(1, 0) = -R(0, 1);
	R(2, 0) = -R(0, 2);
	R(2, 1) = -R(1, 2);
	return R;
}

template<int D>
SquareMatrix<T, D> rotateToPlane(std::array<T, D> const &a, std::array<T, D> const &b) {
	auto const u = normalize(crossProduct(a, b));
	std::array n = { T(0), T(0), T(1) };
	auto const v = crossProduct(u, n);
	T const s = norm(v);
	T const c = dot(u, n);
	auto const vx = crossProductMatrix(v);
	auto const R = SquareMatrix<T, D>::identity() + vx + (vx * vx) * ((1 - c) / (s * s));
	return R;
}

struct ImplicitEqs {
	ImplicitEqs(T Er0_, T F0_x_, T F0_y_, T Eg0_, T Beta0_x_, T Beta0_y_, T rho_, T mu_, T kappa_, T chi_, T gamma_, T dt_) {
		Er0 = Er0_;
		F0_x = F0_x_;
		F0_y = F0_y_;
		Eg0 = Eg0_;
		Beta0_x = Beta0_x_;
		Beta0_y = Beta0_y_;
		rho = rho_;
		mu = mu_;
		kappa = kappa_;
		chi = chi_;
		gamma = gamma_;
		dt = dt_;

	}
	auto jacobian(std::array<T, 3> const &X) const {
		T const dEr = X[2];
		T const dF_x = X[0];
		T const dF_y = X[1];
		T const Beta_x = Beta0_x - dF_x / (rho * cons.c * cons.c);
		T const Beta_y = Beta0_y - dF_y / (rho * cons.c * cons.c);
		T const dBeta_dF = -T(1.0) / (rho * cons.c * cons.c);
		T const Ek = rho * (sqr(cons.c * Beta_x) + sqr(cons.c * Beta_y)) / T(2);
		T const Eg = Eg0 - dEr;
		T const eps = Eg - Ek;
		T const iCv = (eps > T(0)) ? (mu * cons.amu * (gamma - T(1)) / (cons.kB * rho)) : T(0);
		T const temp = eps * iCv;
		if (temp < T(0)) {
			THROW("Invalid temperature");
		}
		T const dtemp_dF_x = Beta_x * iCv;
		T const dtemp_dF_y = Beta_y * iCv;
		T const temp2 = sqr(temp);
		T const temp3 = temp * temp2;
		T const temp4 = sqr(temp2);
		T const F_x = F0_x + dF_x;
		T const F_y = F0_y + dF_y;
		T const F2 = sqr(F_x) + sqr(F_y);
		T const F = sqrt(F2);
		T const Finv = (F ? T(1) : T(0)) / (F ? F : T(1));
		T const N_x = F_x * Finv;
		T const N_y = F_y * Finv;
		T const dNx_dx = +N_y * N_y * Finv;
		T const dNy_dy = +N_x * N_x * Finv;
		T const dNx_dy = -N_x * N_y * Finv;
		T const Er = Er0 + dEr;
		T const Er_inv = T(1) / Er;
		T const f = F * Er_inv;
		T const df_dEr = -f * Er_inv;
		T const df_dF_x = N_x * Er_inv;
		T const df_dF_y = N_y * Er_inv;
		if ((f > T(1)) || (f < T(0))) {
			THROW(std::string("Unphysical flux in implicit solver: f = " + std::to_string(f)));
		}
		T const sqrt4m3f2 = sqrt(T(4) - T(3) * sqr(f));
		T const Xi = (T(5) / T(3) - T(2) / T(3) * sqrt4m3f2);
		T const dXi_df = T(2) * f / sqrt4m3f2;
		T const dXi_dEr = dXi_df * df_dEr;
		T const dXi_dF_x = dXi_df * df_dF_x;
		T const dXi_dF_y = dXi_df * df_dF_y;
		T const strCo = (T(3) / T(2) * Xi - T(1) / T(2));
		T const difCo = (T(1) / T(2) - T(1) / T(2) * Xi);
		T const D_xx = strCo * N_x * N_x + difCo;
		T const D_yy = strCo * N_y * N_y + difCo;
		T const D_xy = strCo * N_x * N_y;
		T const dD_xx_dEr = T(3) / T(2) * dXi_dEr * N_x * N_x - T(1) / T(2) * dXi_dEr;
		T const dD_yy_dEr = T(3) / T(2) * dXi_dEr * N_y * N_y - T(1) / T(2) * dXi_dEr;
		T const dD_xy_dEr = T(3) / T(2) * dXi_dEr * N_x * N_y;
		T const dD_xx_dF_x = T(3) / T(2) * dXi_dF_x * N_x * N_x + strCo * T(2) * (dNx_dx * N_x) - T(1) / T(2) * dXi_dF_x;
		T const dD_xx_dF_y = T(3) / T(2) * dXi_dF_y * N_x * N_x + strCo * T(2) * (dNx_dy * N_x) - T(1) / T(2) * dXi_dF_y;
		T const dD_yy_dF_x = T(3) / T(2) * dXi_dF_x * N_y * N_y + strCo * T(2) * (dNx_dy * N_y) - T(1) / T(2) * dXi_dF_x;
		T const dD_yy_dF_y = T(3) / T(2) * dXi_dF_y * N_y * N_y + strCo * T(2) * (dNy_dy * N_y) - T(1) / T(2) * dXi_dF_y;
		T const dD_xy_dF_x = T(3) / T(2) * dXi_dF_x * N_x * N_y + strCo * (dNx_dx * N_y + dNx_dy * N_x);
		T const dD_xy_dF_y = T(3) / T(2) * dXi_dF_y * N_x * N_y + strCo * (dNx_dy * N_y + dNy_dy * N_x);
		T const P_xx = Er * D_xx;
		T const P_xy = Er * D_xy;
		T const P_yy = Er * D_yy;
		T const dP_xx_dEr = Er * dD_xx_dEr + D_xx;
		T const dP_xy_dEr = Er * dD_xy_dEr + D_xy;
		T const dP_yy_dEr = Er * dD_yy_dEr + D_yy;
		T const dP_xx_dF_x = Er * dD_xx_dF_x;
		T const dP_xx_dF_y = Er * dD_xx_dF_y;
		T const dP_xy_dF_x = Er * dD_xy_dF_x;
		T const dP_xy_dF_y = Er * dD_xy_dF_y;
		T const dP_yy_dF_x = Er * dD_yy_dF_x;
		T const dP_yy_dF_y = Er * dD_yy_dF_y;
		T const aT3 = cons.a * temp3;
		T const gk = kappa * (Er0 + dEr - cons.a * temp4 - T(2) * (Beta_x * F_x + Beta_y * F_y));
		T const dgk_dEr = kappa * (T(1) + T(4) * aT3 * iCv);
		T const dgk_dF_x = -kappa * T(2) * (T(2) * aT3 * dtemp_dF_x + F_x * dBeta_dF + Beta_x);
		T const dgk_dF_y = -kappa * T(2) * (T(2) * aT3 * dtemp_dF_y + F_y * dBeta_dF + Beta_y);
		T const dGx_dEr_x = -chi * (Beta_x + dP_xx_dEr * Beta_x + dP_xy_dEr * Beta_y);
		T const dGx_dEr_y = -chi * (Beta_y + dP_xy_dEr * Beta_x + dP_yy_dEr * Beta_y);
		T const Gx_x = chi * (F_x - Beta_x * Er - P_xx * Beta_x - P_xy * Beta_y);
		T const Gx_y = chi * (F_y - Beta_y * Er - P_xy * Beta_x - P_yy * Beta_y);
		T const dG_dF = dBeta_dF * Er - T(1);
		T const dGx_x_dF_x = -chi * (dP_xx_dF_x * Beta_x + P_xx * dBeta_dF + dP_xy_dF_x * Beta_y + dG_dF);
		T const dGx_y_dF_y = -chi * (dP_xy_dF_y * Beta_x + P_yy * dBeta_dF + dP_yy_dF_y * Beta_y + dG_dF);
		T const dGx_x_dF_y = -chi * (dP_xx_dF_y * Beta_x + P_xy * dBeta_dF + dP_xy_dF_y * Beta_y);
		T const dGx_y_dF_x = -chi * (dP_xy_dF_x * Beta_x + P_xy * dBeta_dF + dP_yy_dF_x * Beta_y);
		T const dgx_dEr = dGx_dEr_x * Beta_x + dGx_dEr_y * Beta_y;
		T const dgx_dF_x = dGx_x_dF_x * Beta_x + dGx_y_dF_x * Beta_y + Gx_x * dBeta_dF;
		T const dgx_dF_y = dGx_x_dF_y * Beta_x + dGx_y_dF_y * Beta_y + Gx_y * dBeta_dF;
		T const dGk_x_dEr = dgk_dEr * Beta_x;
		T const dGk_y_dEr = dgk_dEr * Beta_y;
		T const dGk_x_dF_y = dgk_dF_y * Beta_x;
		T const dGk_y_dF_x = dgk_dF_x * Beta_y;
		T const dGk_x_dF_x = dgk_dF_x * Beta_x + gk * dBeta_dF;
		T const dGk_y_dF_y = dgk_dF_y * Beta_y + gk * dBeta_dF;
		T const cdt = cons.c * dt;
		T const dhr_dEr = T(1) + cdt * (dgk_dEr + dgx_dEr);
		T const dhr_dF_x = cdt * (dgk_dF_x + dgx_dF_x);
		T const dhr_dF_y = cdt * (dgk_dF_y + dgx_dF_y);
		T const dHr_x_dEr = cdt * (dGk_x_dEr + dGx_dEr_x);
		T const dHr_y_dEr = cdt * (dGk_y_dEr + dGx_dEr_y);
		T const dHr_x_dF_x = cdt * (dGk_x_dF_x + dGx_x_dF_x) + T(1);
		T const dHr_y_dF_y = cdt * (dGk_y_dF_y + dGx_y_dF_y) + T(1);
		T const dHr_x_dF_y = cdt * (dGk_x_dF_y + dGx_x_dF_y);
		T const dHr_y_dF_x = cdt * (dGk_y_dF_x + dGx_y_dF_x);
		SquareMatrix<T, NDIM> J;
		J(0, 0) = dHr_x_dF_x;
		J(0, 1) = dHr_x_dF_y;
		J(1, 0) = dHr_y_dF_x;
		J(1, 1) = dHr_y_dF_y;
		J(0, 2) = dHr_x_dEr;
		J(1, 2) = dHr_y_dEr;
		J(2, 0) = dhr_dF_x;
		J(2, 1) = dhr_dF_y;
		J(2, 2) = dhr_dEr;
		return J;
	}
	auto residual(std::array<T, 3> const &X) const {
		T const dEr = X[2];
		T const dF_x = X[0];
		T const dF_y = X[1];
		T const Beta_x = Beta0_x - dF_x / (rho * cons.c * cons.c);
		T const Beta_y = Beta0_y - dF_y / (rho * cons.c * cons.c);
		T const Ek = T(1) / T(2) * rho * (sqr(cons.c * Beta_x) + sqr(cons.c * Beta_y));
		T const Eg = Eg0 - dEr;
		T const eps = Eg - Ek;
		T const iCv = (eps > T(0)) ? (mu * cons.amu * (gamma - T(1)) / (cons.kB * rho)) : T(0);
		T const temp = eps * iCv;
		if (temp < 0) {
			THROW("Invalid temperature");
		}
		T const temp2 = sqr(temp);
		T const temp4 = sqr(temp2);
		T const F_x = F0_x + dF_x;
		T const F_y = F0_y + dF_y;
		T const F2 = sqr(F_x) + sqr(F_y);
		T const F = sqrt(F2);
		T const Finv = (F ? T(1) : T(0)) / (F ? F : T(1));
		T const N_x = F_x * Finv;
		T const N_y = F_y * Finv;
		T const Er = Er0 + dEr;
		T const Er_inv = T(1) / Er;
		T const f = F * Er_inv;
		if ((f > T(1)) || (f < T(0))) {
			THROW(std::string("Unphysical flux in implicit solver: f = " + std::to_string(f)));
		}
		T const sqrt4m3f2 = sqrt(T(4) - T(3) * sqr(f));
		T const Xi = (T(5) / T(3) - T(2) / T(3) * sqrt4m3f2);
		T const strCo = (T(3) / T(2) * Xi - T(1) / T(2));
		T const difCo = (T(1) / T(2) - T(1) / T(2) * Xi);
		T const D_xx = strCo * N_x * N_x + difCo;
		T const D_yy = strCo * N_y * N_y + difCo;
		T const D_xy = strCo * N_x * N_y;
		T const P_xx = Er * D_xx;
		T const P_xy = Er * D_xy;
		T const P_yy = Er * D_yy;
		T const gk = kappa * (Er0 + dEr - cons.a * temp4 - T(2) * (Beta_x * F_x + Beta_y * F_y));
		T const Gx_x = chi * (F_x - Beta_x * Er - P_xx * Beta_x - P_xy * Beta_y);
		T const Gx_y = chi * (F_y - Beta_y * Er - P_xy * Beta_x - P_yy * Beta_y);
		T const gx = Gx_x * Beta_x + Gx_y * Beta_y;
		T const Gk_x = gk * Beta_x;
		T const Gk_y = gk * Beta_y;
		T const cdt = cons.c * dt;
		T const hr = dEr + cdt * (gk + gx);
		T const Hr_x = dF_x + cdt * (Gk_x + Gx_x);
		T const Hr_y = dF_y + cdt * (Gk_y + Gx_y);
		std::array<T, NDIM> dH;
		dH[0] = Hr_x;
		dH[1] = Hr_y;
		dH[2] = hr;
		return dH;
	}
	T Er0;
	T F0_x;
	T F0_y;
	T Eg0;
	T Beta0_x;
	T Beta0_y;
	T rho;
	T mu;
	T kappa;
	T chi;
	T gamma;
	T dt;
};

int solveImplicitRadiation(std::array<T, NDIM + 1> &uR, std::array<T, NDIM + 1> &uG, T rho, T mu, T kappa, T chi, T gamma, T dt) {
	std::array<T, NDIM> X;
	std::array<T, NDIM> Fr;
	std::array<T, NDIM> Beta;
	T betaNorm = T(1) / (rho * cons.c);
	T Ek1 = T(0);
	for (int d = 0; d < NDIM; d++) {
		X[d] = T(0);
		Fr[d] = uR[d];
		Beta[d] = uG[d] * betaNorm;
		Ek1 += sqr(uG[d]) / rho / T(2);
	}
	X[NDIM] = T(0);
	T Er = uR[NDIM];
	T Eg = uG[NDIM];
	T Ei = Eg - Ek1;
	auto const Rot = rotateToPlane<NDIM>(Fr, Beta);
	Fr = Rot * Fr;
	Beta = Rot * Beta;
	T Ek = T(0);
	for (int d = 0; d < NDIM; d++) {
		X[d] = T(0);
		Fr[d] = uR[d];
		Beta[d] = uG[d] * betaNorm;
		Ek += sqr(Beta[d] / betaNorm) * T(0.5) / rho;
	}
	Eg = Ei + Ek;
	T err = 1;
	int n = 0;
	ImplicitEqs const eqs(Er, Fr[0], Fr[1], Eg, Beta[0], Beta[1], rho, mu, kappa, chi, gamma, dt);
	auto dH = eqs.residual(X);
	auto J = eqs.jacobian(X);
	auto Jinv = matrixInverse(J);
	auto RR = dot(dH, dH);
	auto R0 = dot(dH, dH);
	T R, f;
	do {
		T alpha = T(1);
		auto const X0 = X;
		while (1) {
			bool flag = false;
			X = X0 - alpha * Jinv * dH;
			if ((Er + X[2]) > T(0)) {
				f = sqrt(sqr(Fr[0] + X[0]) + sqr(Fr[1] + X[1])) / (Er + X[2]);
				if (f < T(1)) {
					auto const dH1 = eqs.residual(X);
					R = dot(dH1, dH1);
					err = sqrt(R / RR);
					if (R < (T(1) - (1e-3) * alpha) * R0) {
						printf("%i    %e %e %e %e\n", n, err, alpha, R, R0);
						R0 = R;
						flag = true;
					}
				}
			}
			alpha *= T(1) / T(2);
			if (flag) {
				J = eqs.jacobian(X);
				Jinv = matrixInverse(J);
				dH = eqs.residual(X);
				break;
			}
		}
		++n;
		if (n > 100) {
			T Ek = (Beta[0] * Beta[0] + Beta[1] * Beta[1]) * (0.5 * rho * cons.c * cons.c);
			T Ek1 = (0.5 * rho * cons.c * cons.c) * (sqr(Beta[0] - X[0] / (rho * cons.c * cons.c)) + sqr(Beta[1] - X[1] / (rho * cons.c * cons.c)));
			Jinv = matrixInverse(eqs.jacobian(X));
			dH = eqs.residual(X);
			auto dX = Jinv * dH;
			printf("X  = %e %e %e\n", X[0], X[1], X[2]);
			printf("dX = %e %e %e\n", dX[0], dX[1], dX[2]);
			printf("Er    = %e -> %e\n", Er, Er + X[2]);
			printf("Fx    = %e -> %e\n", Fr[0], Fr[0] + X[0]);
			printf("Fy    = %e -> %e\n", Fr[1], Fr[1] + X[0]);
			printf("Eg    = %e -> %e\n", Eg, Eg - X[2]);
			printf("Ek    = %e -> %e\n", Ek, Ek1);
			printf("Betax = %e -> %e\n", Beta[0], Beta[0] - X[0] / (rho * cons.c * cons.c));
			printf("Betay = %e -> %e\n", Beta[1], Beta[1] - X[1] / (rho * cons.c * cons.c));
			THROW("Solver failed to converge.");
		}
	} while (err > 1e-2);
	auto const RotInv = matrixInverse(Rot);
	Fr[0] += X[0];
	Fr[1] += X[1];
	Fr[2] = T(0);
	Beta[0] -= X[0] / (rho * cons.c * cons.c);
	Beta[1] -= X[1] / (rho * cons.c * cons.c);
	Beta[2] = T(0);
	Fr = RotInv * Fr;
	Beta = RotInv * Beta;
	for (int d = 0; d < NDIM; d++) {
		uR[d] = Fr[d];
		uG[d] = Beta[d] * rho * cons.c;
	}
	uR[NDIM] = Er + X[2];
	uG[NDIM] = Eg - X[2];
	return n;
}

void testRadiation() {
//	constexpr int D = 3;
//	constexpr int N = 32;
//	using T = Real;
//	using RK = typename RungeKutta<T, P>::type;
//	using S = RadiationState<T>;
//	HyperGrid<S, N, P, RK> grid;
//	std::cout << cons;
//	std::cout << physicalUnits;
//	grid.initialize(initRadiationFront<T, D>);
//	grid.enforceBoundaryConditions();
//	T t = T(0);
//	T tmax = T(.15);
//	T dt;
//	RK const rk;
//	int iter = 0;
//	while (t < tmax) {
//		grid.output("X", iter, t);
//		std::cout << "i = " << std::to_string(iter);
//		std::cout << "  t = " << std::to_string(t);
//		dt = grid.beginStep();
//		std::cout << "  dt = " << dt << std::endl;
//		for (int s = 0; s < rk.stageCount(); s++) {
//			grid.subStep(dt, s);
//			grid.enforceBoundaryConditions();
//		}
//		grid.endStep();
//		grid.enforceBoundaryConditions();
//		iter++;
//		t += dt;
//	}
	std::cout << physicalUnits << "\n";
	int ntrial = 100000;
	double err_max = 0.0;
	for (int trialNum = 0; trialNum < ntrial; trialNum++) {
		printf("%i\n", trialNum);
		std::array<T, NDIM + 1> uR;
		std::array<T, NDIM + 1> uG;
		T Tr = randomLog(10, 1000000);
		T Tg = Tr * randomLogNormal(1.0, .1);
		auto const v1 = randomUnitVector3D() * rand1();
		auto const v2 = randomUnitVector3D() * rand1();
		uR[NDIM] = cons.a * sqr(sqr(Tr));
		T rho = 1.0e-3;
		T mu = 4.0 / 3.0;
		T chi = 1.0;
		T kappa = 1.0;
		T dt = 1.0;
		for (int d = 0; d < NDIM; d++) {
			uR[d] = v1[d] * uR[NDIM];
			uG[d] = v2[d] * cons.c * rho;
		}
		constexpr T gamma = 5.0 / 3.0;
		T const eint = cons.kB * rho * Tg / ((gamma - 1.0) * mu * cons.amu);
		uG[NDIM] = eint;
		uG[NDIM] += 0.5 * sqr(uG[0]) / rho;
		uG[NDIM] += 0.5 * sqr(uG[1]) / rho;
		uG[NDIM] += 0.5 * sqr(uG[2]) / rho;
		T f = sqrt(sqr(uR[0]) + sqr(uR[1]) + sqr(uR[2])) / uR[NDIM];
		printf("Er = %e ", uR[NDIM]);
		printf("Eg = %e ", uG[NDIM]);
		printf("Ei = %e ", eint);
		printf("Tr = %e ", Tr);
		printf("Tg = %e ", Tg);
		printf("f = %e\n", f);
		printf("Fr = %e %e %e ", uR[0], uR[1], uR[2]);
		printf("Beta = %e %e %e \n", uG[0] / (rho * cons.c), uG[1] / (rho * cons.c), uG[2] / (rho * cons.c));
		auto const eps = sqrt(std::numeric_limits<double>::epsilon());
		auto sol = solveImplicitRadiation(uR, uG, rho, mu, kappa, chi, gamma, dt);

		//		std::array<std::array<T, NDIM + 1>, NDIM + 1> dJ2;
//		auto dJ1 = computeG2D(0, 0, 0, Er0, F0x, F0y, Eg0, Beta0x, Beta0y, rho, mu, kappa, chi, gamma, dt).second;
//		for (int l = 0; l <= NDIM; l++) {
//			if (l == 2) {
//				continue;
//			}
//			T dErp = T(0), dErm = T(0);
//			T dFmx = T(0);
//			T dFmy = T(0);
//			T dFpx = T(0);
//			T dFpy = T(0);
//			if (l == NDIM) {
//				dErp = 0.5 * dEr;
//				dErm = -0.5 * dEr;
//			} else {
//				if (l == 0) {
//					dFpx += 0.5 * dEr;
//					dFmx -= 0.5 * dEr;
//				} else {
//					dFpy += 0.5 * dEr;
//					dFmy -= 0.5 * dEr;
//				}
//			}
//			//	(T dEr, T dF_x, T dF_y, T Er0, T F0_x, T F0_y, T Eg0, T Beta0_x, T Beta0_y, T rho, T mu, T kappa, T chi, T gamma, T dt)
//			auto dJp = computeG2D(dErp, dFpx, dFpy, Er0, F0x, F0y, Eg0, Beta0x, Beta0y, rho, mu, kappa, chi, gamma, dt).first;
//			auto dJm = computeG2D(dErm, dFmx, dFmy, Er0, F0x, F0y, Eg0, Beta0x, Beta0y, rho, mu, kappa, chi, gamma, dt).first;
//			for (int m = 0; m <= NDIM; m++) {
//				if (m == 2) {
//					continue;
//				}
//				dJ2[m][l] = (dJp[m] - dJm[m]) / (eps * Er);
//			}
//		}
//		double err = 0.0;
//		double norm = 0.0;
//		for (int l = 0; l < NDIM + 1; l++) {
//			if (l == 2) {
//				continue;
//			}
//			for (int m = 0; m < NDIM + 1; m++) {
//				if (m == 2) {
//					continue;
//				}
//				err += sqr(dJ1[l][m] - dJ2[l][m]);
//				norm += 0.5 * (std::abs(dJ1[l][m]) + std::abs(dJ2[l][m]));
////				printf("n = %i m = %i analytic %16.8e numerical %16.8e abs error %16.8e rel error %16.8e\n", l, m, dF1[l][m], dF2[l][m], dF1[l][m] - dF2[l][m],
////						(dF1[l][m] - dF2[l][m]) / (dF1[l][m] + 1e-20));
//			}
//		}
//		err = sqrt(err);
//		err /= norm;
//		printf("%e\n", err);
//		if (err > 0.1) {
//			printf("Tr = %e\n", Tr);
//			printf("Tg = %e\n", Tg);
//			printf("f = %e\n", f);
//			printf("F0x = %e\n", F0x);
//			printf("F0y = %e\n", F0y);
//			printf("Beta0x = %e\n", Beta0x);
//			printf("Beta0y = %e\n", Beta0y);
//			for (int l = 0; l < NDIM + 1; l++) {
//				if (l == 2) {
//					continue;
//				}
//				for (int m = 0; m < NDIM + 1; m++) {
//					if (m == 2) {
//						continue;
//					}
//					printf("n = %i m = %i analytic %16.8e numerical %16.8e abs error %16.8e rel error %16.8e\n", l, m, dJ1[l][m], dJ2[l][m], dJ1[l][m] - dJ2[l][m],
//							(dJ1[l][m] - dJ2[l][m]) / (dJ1[l][m] + 1e-20));
//				}
//			}
//			abort();
//		}
//
//		err_max = std::max(err, err_max);
	}
}
