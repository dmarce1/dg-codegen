#include "Constants.hpp"
#include "EulerState.hpp"
#include <array>

using T = double;
constexpr int NDIM = 3;

auto computeG(T Er, std::array<T, NDIM> F, T Er0, std::array<T, NDIM> F0, T Eg0, std::array<T, NDIM> Beta0, T rho, T mu, T kappa, T chi, T gamma, T dt) {
	constexpr T tiny = 1e-50;
	std::array<T, NDIM> Beta;
	for (int n = 0; n < NDIM; n++) {
		Beta[n] = Beta0[n] + (F0[n] - F[n]) / (rho * cgs::c * cgs::c);
	}
	T Ek = 0;
	std::array<T, NDIM> dEk_dF;
	for (int n = 0; n < NDIM; n++) {
		Ek += 0.5 * rho * sqr(cgs::c * Beta[n]);
		dEk_dF[n] = -Beta[n];
	}
	T const Eg = Eg0 + Er0 - Er;
	T const eps = Eg - Ek;
	std::array<T, NDIM> const deps_dF = -dEk_dF;
	T F2 = 0;
	for (int n = 0; n < NDIM; n++) {
		F2 += sqr(F[n]);
	}
	T const absF = sqrt(F2);
	T const iCv = (mu * cgs::amu) * (gamma - 1.0) / (cgs::kB * rho);
	T const temp = eps * iCv;
	T const dtemp_dEr = -iCv;
	std::array<T, NDIM> const dtemp_dF = deps_dF * iCv;
	T const temp2 = sqr(temp);
	T const temp4 = sqr(temp2);
	std::array<std::array<T, NDIM>, NDIM> dBeta_dF;
	for (int n = 0; n < NDIM; n++) {
		for (int m = 0; m < NDIM; m++) {
			dBeta_dF[n][m] = -T(n == m) / (rho * cgs::c * cgs::c);
		}
	}
	std::array<T, NDIM> N;
	for (int n = 0; n < NDIM; n++) {
		N[n] = F[n] / (absF + tiny);
	}
	T const temp3 = temp * temp2;
	std::array<std::array<T, NDIM>, NDIM> dN_dF;
	for (int n = 0; n < NDIM; n++) {
		for (int m = 0; m < NDIM; m++) {
			dN_dF[n][m] = (T(n == m) - N[n] * N[m]) / (absF + tiny);
		}
	}
	std::array<T, NDIM> df_dF;
	T const f = absF / Er;
	T const df_dEr = -f / Er;
	for (int n = 0; n < NDIM; n++) {
		df_dF[n] = N[n] / Er;
	}
	T const Xi = (5 - 2 * sqrt(4 - 3 * sqr(f))) / 3;
	T const dXi_df = 2 * f / sqrt(4 - 3 * sqr(f));
	T const dXi_dEr = dXi_df * df_dEr;
	std::array<T, NDIM> dXi_dF;
	for (int n = 0; n < NDIM; n++) {
		dXi_dF[n] = dXi_df * df_dF[n];
	}
	std::array<std::array<T, NDIM>, NDIM> D;
	for (int n = 0; n < NDIM; n++) {
		for (int m = 0; m < NDIM; m++) {
			D[n][m] = ((1 - Xi) / 2 * T(n == m) + (3 * Xi - 1) / 2 * N[n] * N[m]);
		}
	}
	std::array<std::array<T, NDIM>, NDIM> dD_dEr;
	for (int n = 0; n < NDIM; n++) {
		for (int m = 0; m < NDIM; m++) {
			dD_dEr[n][m] = -dXi_dEr / 2 * T(n == m);
			dD_dEr[n][m] += 3 * dXi_dEr / 2 * N[n] * N[m];
		}
	}
	std::array<std::array<std::array<T, NDIM>, NDIM>, NDIM> dD_dF;
	for (int n = 0; n < NDIM; n++) {
		for (int m = 0; m < NDIM; m++) {
			for (int l = 0; l < NDIM; l++) {
				dD_dF[n][m][l] = -dXi_dF[l] / 2 * T(n == m);
				dD_dF[n][m][l] += 3 * dXi_dF[l] / 2 * N[n] * N[m];
				dD_dF[n][m][l] += (3 * Xi - 1) / 2 * (dN_dF[n][l] * N[m] + dN_dF[m][l] * N[n]);
			}
		}
	}
	std::array<std::array<T, NDIM>, NDIM> P;
	for (int n = 0; n < NDIM; n++) {
		for (int m = 0; m < NDIM; m++) {
			P[n][m] = Er * D[n][m];
		}
	}
	std::array<std::array<T, NDIM>, NDIM> dP_dEr;
	for (int n = 0; n < NDIM; n++) {
		for (int m = 0; m < NDIM; m++) {
			dP_dEr[n][m] = Er * dD_dEr[n][m] + D[n][m];
		}
	}
	std::array<std::array<std::array<T, NDIM>, NDIM>, NDIM> dP_dF;
	for (int n = 0; n < NDIM; n++) {
		for (int m = 0; m < NDIM; m++) {
			for (int l = 0; l < NDIM; l++) {
				dP_dF[l][n][m] = Er * dD_dF[l][n][m];
			}
		}
	}
	T gk = kappa * (Er - cgs::a * temp4);
	for (int k = 0; k < NDIM; k++) {
		gk -= 2 * kappa * Beta[k] * F[k];
	}
	T const dgk_dEr = kappa * (1 - 4 * cgs::a * temp3 * dtemp_dEr);
	std::array<T, NDIM> dgk_dF;
	for (int n = 0; n < NDIM; n++) {
		dgk_dF[n] = -kappa * (2 * Beta[n] + 4 * cgs::a * temp3 * dtemp_dF[n]);
		for (int k = 0; k < NDIM; k++) {
			dgk_dF[n] -= 2 * kappa * F[k] * dBeta_dF[k][n];
		}
	}
	std::array<T, NDIM> Gx { };
	std::array<T, NDIM> dGx_dEr { };
	std::array<std::array<T, NDIM>, NDIM> dGx_dF { };
	for (int n = 0; n < NDIM; n++) {
		Gx[n] = chi * (F[n] - Beta[n] * Er);
		for (int k = 0; k < NDIM; k++) {
			Gx[n] -= chi * P[k][n] * Beta[k];
		}
	}
	for (int n = 0; n < NDIM; n++) {
		dGx_dEr[n] = -chi * Beta[n];
		for (int k = 0; k < NDIM; k++) {
			dGx_dEr[n] -= chi * dP_dEr[k][n] * Beta[k];
		}
	}
	for (int n = 0; n < NDIM; n++) {
		for (int m = 0; m < NDIM; m++) {
			dGx_dF[n][m] = chi * (T(n == m) - dBeta_dF[n][m] * Er);
			for (int k = 0; k < NDIM; k++) {
				dGx_dF[n][m] -= chi * dP_dF[k][n][m] * Beta[k];
				dGx_dF[n][m] -= chi * P[k][n] * dBeta_dF[k][m];
			}
		}
	}
	T gx = 0;
	T dgx_dEr = 0;
	std::array<T, NDIM> dgx_dF { };
	for (int k = 0; k < NDIM; k++) {
		gx += Gx[k] * Beta[k];
	}
	for (int k = 0; k < NDIM; k++) {
		dgx_dEr += dGx_dEr[k] * Beta[k];
	}
	for (int n = 0; n < NDIM; n++) {
		dgx_dF[n] = 0;
		for (int k = 0; k < NDIM; k++) {
			dgx_dF[n] += dGx_dF[k][n] * Beta[k] + dBeta_dF[k][n] * Gx[k];
		}
	}
	std::array<T, NDIM> Gk { };
	for (int n = 0; n < NDIM; n++) {
		Gk[n] = gk * Beta[n];
	}
	std::array<T, NDIM> dGk_dEr;
	for (int n = 0; n < NDIM; n++) {
		dGk_dEr[n] = dgk_dEr * Beta[n];
	}
	std::array<std::array<T, NDIM>, NDIM> dGk_dF;
	for (int n = 0; n < NDIM; n++) {
		for (int m = 0; m < NDIM; m++) {
			dGk_dF[n][m] = Beta[n] * dgk_dF[m] + gk * dBeta_dF[n][m];
		}
	}
	T const hr = Er - Er0 + dt * cgs::c * (gk + gx);
	T const dhr_dEr = 1 + dt * cgs::c * (dgk_dEr + dgx_dEr);
	std::array<T, NDIM> dhr_dF;
	for (int n = 0; n < NDIM; n++) {
		dhr_dF[n] = dt * cgs::c * (dgk_dF[n] + dgx_dF[n]);
	}
	std::array<T, NDIM> Hr;
	for (int n = 0; n < NDIM; n++) {
		Hr[n] = F[n] - F0[n] + dt * cgs::c * (Gk[n] + Gx[n]);
	}
	std::array<T, NDIM> dHr_dEr;
	for (int n = 0; n < NDIM; n++) {
		dHr_dEr[n] = dt * cgs::c * (dGk_dEr[n] + dGx_dEr[n]);
	}
	std::array<std::array<T, NDIM>, NDIM> dHr_dF;
	for (int n = 0; n < NDIM; n++) {
		for (int m = 0; m < NDIM; m++) {
			dHr_dF[n][m] = T(n == m) + dt * cgs::c * (dGk_dF[n][m] + dGx_dF[n][m]);
		}
	}

	std::pair<std::array<T, NDIM + 1>, std::array<std::array<T, NDIM + 1>, NDIM + 1>> rc;
	std::array<T, NDIM + 1> &F4 = rc.first;
	std::array<std::array<T, NDIM + 1>, NDIM + 1> &dF4 = rc.second;
	for (int k = 0; k < NDIM; k++) {
		F4[k] = Hr[k];
		for (int n = 0; n < NDIM; n++) {
			dF4[k][n] = dHr_dF[k][n];
		}
		dF4[NDIM][k] = dhr_dF[k];
		dF4[k][NDIM] = dHr_dEr[k];
	}
	F4[NDIM] = hr;
	dF4[NDIM][NDIM] = dhr_dEr;
	return rc;
}

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
void testRadiation() {
	int ntrial = 100000;
	double err_max = 0.0;
	for (int trialNum = 0; trialNum < ntrial; trialNum++) {

		T Tr = randomLog(10, 1000000);
		T Tg = Tr * randomLogNormal(1.0, 1);
		T Er = cgs::a * sqr(sqr(Tr));
		T const f = rand1();
		std::array<T, NDIM> const Nr = randomUnitVector3D();
		std::array<T, NDIM> const Ng = randomUnitVector3D();
		std::array<T, NDIM> F0 { };
		std::array<T, NDIM> Beta0 = Ng;
		for (int i = 0; i < NDIM; i++) {
			F0[i] = Nr[i] * f * Er;
		}
		T rho = 1.0;
		T mu = 4.0 / 3.0;
		T chi = 1.0;
		T kappa = 1.0;
		T dt = 1.0;
		T Er0 = Er;
		std::array<T, NDIM> F = F0;
		constexpr T gamma = 5.0 / 3.0;
		T Eg0 = cgs::kB * rho * Tg / ((gamma - 1.0) * mu * cgs::amu);
		for (int l = 0; l < NDIM; l++) {
			Eg0 += 0.5 * rho * sqr(cgs::c * Beta0[l]);
		}
		auto const eps = sqrt(std::numeric_limits<double>::epsilon());
	//	auto const eps = 1e-4;
		T dEr = eps * Er;
		T dF = eps * Er;

		std::array<std::array<T, NDIM + 1>, NDIM + 1> dF2;
		auto dF1 = computeG(Er, F, Er0, F0, Eg0, Beta0, rho, mu, kappa, chi, gamma, dt).second;
		for (int l = 0; l <= NDIM; l++) {
			T Erp = Er, Erm = Er;
			std::array<T, NDIM> Fp = F, Fm = F;
			if (l == NDIM) {
				Erp = Er + 0.5 * dEr;
				Erm = Er - 0.5 * dEr;
			} else {
				Fp[l] += 0.5 * dEr;
				Fm[l] -= 0.5 * dEr;
			}
			auto dFp = computeG(Erp, Fp, Er0, F0, Eg0, Beta0, rho, mu, kappa, chi, gamma, dt).first;
			auto dFm = computeG(Erm, Fm, Er0, F0, Eg0, Beta0, rho, mu, kappa, chi, gamma, dt).first;
			for (int m = 0; m <= NDIM; m++) {
				dF2[m][l] = (dFp[m] - dFm[m]) / (eps * Er);
			}
		}
		double err = 0.0;
		double norm = 0.0;
		for (int l = 0; l < NDIM + 1; l++) {
			for (int m = 0; m < NDIM + 1; m++) {
				err += sqr(dF1[l][m] - dF2[l][m]);
				norm += 0.5 * (std::abs(dF1[l][m]) + std::abs(dF2[l][m]));
//				printf("n = %i m = %i analytic %16.8e numerical %16.8e abs error %16.8e rel error %16.8e\n", l, m, dF1[l][m], dF2[l][m], dF1[l][m] - dF2[l][m],
//						(dF1[l][m] - dF2[l][m]) / (dF1[l][m] + 1e-20));
			}
		}
		err = sqrt(err);
		err /= norm;
		printf("%e\n", err);
		if (err > 0.5) {
			printf("Tr = %e\n", Tr);
			printf("Tg = %e\n", Tg);
			printf("f = (%e) ", f);
			for (int i = 0; i < NDIM; i++) {
				F0[i] = Nr[i] * f * Er;
				printf("%e ", Nr[i] * f);
			}
			printf("\nBeta = ");
			for (int i = 0; i < NDIM; i++) {
				printf("%e ", Ng[i]);
			}
			printf("\n");
			for (int l = 0; l < NDIM + 1; l++) {
				for (int m = 0; m < NDIM + 1; m++) {
					printf("n = %i m = %i analytic %16.8e numerical %16.8e abs error %16.8e rel error %16.8e\n", l, m, dF1[l][m], dF2[l][m], dF1[l][m] - dF2[l][m],
							(dF1[l][m] - dF2[l][m]) / (dF1[l][m] + 1e-20));
				}
			}
			abort();
		}

		err_max = std::max(err, err_max);
	}
	printf("-->> %e\n", err_max);
}
