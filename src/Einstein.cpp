/*
 * Einstein.cpp
 *
 *  Created on: Feb 25, 2025
 *      Author: dmarce1
 */

#include <algorithm>

static constexpr int NDIM2 = NDIM * (NDIM + 1) / 2;

template<typename T>
struct SpacetimeFlux {
private:
	T Theta_;
	T Z_i[NDIM];
	T K_ij[NDIM2];
	T A_k[NDIM];
	T B_kI[NDIM][NDIM];
	T D_kij[NDIM][NDIM2];
};

template<typename T>
struct Spacetime {
	static constexpr T _zeta = T(0);
	static constexpr T _zeta1 = T(0);
	static constexpr T _f = T(1);
	static constexpr T _m = T(2);
	static constexpr T _xi = T(1);
private:
	T alpha_;
	T beta_I[NDIM];
	T gamma_ij[NDIM2];
	T Theta_;
	T Z_i[NDIM];
	T K_ij[NDIM2];
	T A_k[NDIM];
	T B_kI[NDIM][NDIM];
	T D_kij[NDIM][NDIM2];

	static constexpr int si1(int i, int j) {
		return ((i * i + i) >> 1) + j;
	}

	static constexpr int si2(int i, int j) {
		return si1(std::max(i, j), std::min(i, j));
	}

	static constexpr T delta(int i, int j) {
		return T(i == j);
	}

	Spacetime flux(int k) {

		SpacetimeFlux<T> F;
		T gamma_, gamma_inv, K_, B_, Q_;
		T LambdaK_i[NDIM], Q_i[NDIM], V_i[NDIM];
		T gamma_IJ[NDIM2], lambdaK_ij[NDIM2], Q_ij[NDIM2];

		gamma_IJ[si1(0, 0)] = gamma_ij[2] * gamma_ij[5] - gamma_ij[4] * gamma_ij[4];
		gamma_IJ[si1(1, 0)] = gamma_ij[3] * gamma_ij[4] - gamma_ij[1] * gamma_ij[5];
		gamma_IJ[si1(2, 0)] = gamma_ij[1] * gamma_ij[4] - gamma_ij[2] * gamma_ij[3];
		gamma_IJ[si1(1, 1)] = gamma_ij[0] * gamma_ij[5] - gamma_ij[3] * gamma_ij[3];
		gamma_IJ[si1(2, 1)] = gamma_ij[1] * gamma_ij[3] - gamma_ij[0] * gamma_ij[4];
		gamma_IJ[si1(2, 2)] = gamma_ij[0] * gamma_ij[2] - gamma_ij[1] * gamma_ij[1];
		gamma_ = gamma_ij[0] * gamma_IJ[si1(0, 0)] + gamma_ij[1] * gamma_IJ[si1(1, 0)] + gamma_ij[3] * gamma_IJ[si1(2, 0)];
		gamma_inv = T(1) / gamma_;
		for (int i = 0; i < NDIM; i++) {
			for (int j = i; j < NDIM; j++) {
				gamma_IJ[si1(j, i)] *= gamma_inv;
			}
		}
		B_ = T(0);
		for (int m = 0; m < NDIM; m++) {
			B_ + B_kI[m][m];
		}
		K_ = T(0);
		for (int m = 0; m < NDIM; m++) {
			for (int n = 0; n <= m; n++) {
				K_ += (T(2) - delta(n, m)) * gamma_IJ[si1(m, n)] * K_ij[si1(m, n)];
			}
		}
		for (int i = 0; i < NDIM; i++) {
			V_i[i] = -Z_i[i];
			for (int m = 0; m < NDIM; m++) {
				for (int n = 0; n < NDIM; n++) {
					if (i != n) {
						V_i[i] += gamma_IJ[si2(n, m)] * (D_kij[i][si2(m, n)] - D_kij[n][si2(m, i)]);
					}
				}
			}
		}
		Q_ = _f * (K_ - _m * Theta_);
		for (int i = 0; i < NDIM; i++) {
			Q_i[i] = LambdaK_i[i] = A_k[i] - T(2) * Z_i[i];
			for (int m = 0; m < NDIM; m++) {
				for (int n = m; n < NDIM; n++) {
					auto const termD = D_kij[i][si2(n, m)] * gamma_IJ[si1(m, n)];
					auto const termE = D_kij[n][si2(m, i)] * gamma_IJ[si1(m, n)];
					LambdaK_i[i] += (T(2) - delta(n, m)) * (termD - (T(1) - _zeta) * termE);
					Q_i[i] += (T(2) - delta(n, m)) * (termD - T(2) * termE);
				}
			}
			Q_i[i] *= alpha_;
		}

		for (int i = 0; i < NDIM; i++) {
			F.A_k[i] = -beta_I[k] * A_k[i] + delta(i, k) * alpha_ * Q_;
		}

		for (int i = 0; i < NDIM; i++) {
			for (int j = 0; j < NDIM; j++) {
				F.B_kI[i][j] = -beta_I[k] * B_kI[i][j];
				for (int m = 0; m < NDIM; m++) {
					F.B_kI[i][j] += delta(i, k) * alpha_ * Q_i[m] * gamma_IJ[si2(m, j)];
				}
			}
		}

		for (int l = 0; l < NDIM; l++) {
			for (int i = 0; i < NDIM; i++) {
				for (int j = 0; j <= i; j++) {
					int const ij = si1(i, j);
					F.D_kij[l][ij] = -beta_I[k] * D_kij[l][ij];
				}
			}
		}
		for (int i = 0; i < NDIM; i++) {
			for (int j = 0; j <= i; j++) {
				int const ij = si1(i, j);
				F.D_kij[k][ij] += alpha_ * K_ij[si1(i, j)];
				for (int m = 0; m < NDIM; m++) {
					F.D_kij[k][ij] += T(0.5) * (B_kI[i][m] * gamma_ij[si2(m, j)] + B_kI[j][m] * gamma_ij[si2(m, i)]);
				}
			}
		}

		for (int i = 0; i < NDIM; i++) {
			for (int j = 0; j <= i; j++) {
				int const ij = si1(i, j);
				lambdaK_ij[ij] = T(0.5) * (delta(k, i) * LambdaK_i[j] + delta(k, j) * LambdaK_i[i]);
				for (int m = 0; m < NDIM; m++) {
					lambdaK_ij[ij] += gamma_IJ[si2(k, m)] * (D_kij[m][ij] - T(0.5) * (T(1) + _zeta) * (D_kij[j][si2(i, m)] + D_kij[i][si2(j, m)]));
				}
				F.K_ij[si1(i, j)] = -beta_I[k] * K_ij[ij] + alpha_ * lambdaK_ij[ij];
			}
		}

		F.Theta_ = -beta_I[k] * Theta_ + alpha_ * V_i[k];

		for (int i = 0; i < NDIM; i++) {
			F.Z_i[i] = delta(i, k) * (K_ - Theta_);
			for (int m = 0; m < NDIM; m++) {
				F.Z_i[i] -= K_ij[si2(i, m)] * gamma_IJ[si2(m, k)];
			}
			F.Z_i[i] = -beta_I[k] * Z_i[i] + alpha_ * F.Z_i[i] + _zeta1 * (B_kI[i][k] - delta(k, i) * B_);
		}

		return F;
	}

	Spacetime source() {

		SpacetimeFlux<T> S;
		T gamma_, gamma_inv;
		T gamma_IJ[NDIM2];
		T D_Kij[NDIM][NDIM2];
		T D_kiJ[NDIM][NDIM][NDIM];
		T D_k[NDIM];
		T E_k[NDIM];
		T Gamma_Kij[NDIM][NDIM2];

		gamma_IJ[si1(0, 0)] = gamma_ij[2] * gamma_ij[5] - gamma_ij[4] * gamma_ij[4];
		gamma_IJ[si1(1, 0)] = gamma_ij[3] * gamma_ij[4] - gamma_ij[1] * gamma_ij[5];
		gamma_IJ[si1(2, 0)] = gamma_ij[1] * gamma_ij[4] - gamma_ij[2] * gamma_ij[3];
		gamma_IJ[si1(1, 1)] = gamma_ij[0] * gamma_ij[5] - gamma_ij[3] * gamma_ij[3];
		gamma_IJ[si1(2, 1)] = gamma_ij[1] * gamma_ij[3] - gamma_ij[0] * gamma_ij[4];
		gamma_IJ[si1(2, 2)] = gamma_ij[0] * gamma_ij[2] - gamma_ij[1] * gamma_ij[1];
		gamma_ = gamma_ij[0] * gamma_IJ[si1(0, 0)] + gamma_ij[1] * gamma_IJ[si1(1, 0)] + gamma_ij[3] * gamma_IJ[si1(2, 0)];
		gamma_inv = T(1) / gamma_;
		for (int i = 0; i < NDIM; i++) {
			for (int j = i; j < NDIM; j++) {
				gamma_IJ[si1(j, i)] *= gamma_inv;
			}
		}
		for (int k = 0; k < NDIM; k++) {
			for (int i = 0; i < NDIM; i++) {
				for (int j = 0; j <= i; j++) {
					D_Kij[k][si1(i, j)] = T(0);
					D_kiJ[k][i][j] = D_kiJ[k][j][i] = T(0);
					for (int m = 0; m < NDIM; m++) {
						D_Kij[k][si2(i, j)] += gamma_IJ[si2(k, m)] * D_kij[m][si2(i, j)];
						D_kiJ[k][i][j] += D_kij[k][si2(i, m)] * gamma_IJ[si2(m, j)];
						if (i != j) {
							D_kiJ[k][j][i] += D_kij[k][si2(j, m)] * gamma_IJ[si2(m, i)];
						}
					}
				}
			}
		}
		for (int k = 0; k < NDIM; k++) {
			for (int i = 0; i < NDIM; i++) {
				for (int j = 0; j <= i; j++) {
					Gamma_Kij[k][si2(i, j)] = T(0.5) * (D_Kij[i][si2(k, j)] + D_Kij[j][si2(i, k)] - D_Kij[k][si2(i, j)]);
				}
			}
		}
		for (int i = 0; i < NDIM; i++) {
			D_k[i] = T(0);
			E_k[i] = T(0);
			for (int k = 0; k < NDIM; k++) {
				D_k[i] += D_kiJ[i][k][k];
				E_k[i] += D_Kij[k][si2(i, k)];
			}
		}

		for (int i = 0; i < NDIM; i++) {
			for (int j = 0; j <= i; j++) {
				S.K_ij[si1(i, j)] = T(0.5) * (A_k[i] * D_k[j] + A_k[j] * D_k[i]);
				S.K_ij[si1(i, j)] -= T(0.5) * (T(1) - _xi) * (A_k[i] * E_k[j] + A_k[j] * E_k[i]);
				for (int k = 0; k < NDIM; k++) {
					S.K_ij[si1(i, j)] -= T(0.5) * (T(1) + _xi) * A_k[k] * Gamma_Kij[k][si1(i, j)];
					S.K_ij[si1(i, j)] += T(0.5) * (T(1) - _xi) * A_k[k] * D_Kij[k][si1(i, j)];
					S.K_ij[si1(i, j)] -= (T(1) - _xi) * E_k[k] * (D_kiJ[i][j][k] + D_kiJ[j][i][k]);
					for (int m = 0; m < NDIM; m++) {
						S.K_ij[si1(i, j)] -= (T(1) + _xi) * (D_kiJ[i][k][m] * D_kiJ[k][si2(m, j)] + D_kiJ[j][k][m] * D_kiJ[k][si2(m, i)]);
					}
				}
			}
		}

		return S;
	}
};
;

void testEinstein() {
}

