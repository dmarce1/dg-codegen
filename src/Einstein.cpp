/*
 * Einstein.cpp
 *
 *  Created on: Feb 25, 2025
 *      Author: dmarce1
 */

#include <numbers>

#include "Tensor.hpp"

template<typename T>
struct Z4SpaceTime {
	static constexpr int NDIM2 = NDIM * (NDIM + 1) / 2;
	static constexpr T zero = T(0);
	static constexpr T half = T(0.5);
	static constexpr T one = T(1);
	static constexpr T two = T(2);
	using tensor0 = T;
	using tensor1 = std::array<T, NDIM>;
	using tensor2 = std::array<std::array<T, NDIM>, NDIM>;
	using tensorS = std::array<T, NDIM2>;
	using tensor3 = std::array<std::array<T, NDIM2>, NDIM>;
private:
	tensor0 alpha0;
	tensor1 beta1;
	tensorS gamma2;
	tensor1 A1;
	tensor2 B2;
	tensorS K2;
	tensor3 D3;
	tensor0 Theta0;
	tensor1 Z1;
	static inline T f() {
		return one;
	}
	tensorS matrix_inverse(tensorS const &g) {
		tensorS ginv;
		T g00 = (g[2] * g[5] - g[4] * g[4]);
		T g01 = (g[1] * g[5] - g[4] * g[3]);
		T g02 = (g[1] * g[4] - g[2] * g[3]);
		T const det = g00 * g[0] - g01 * g[1] + g02 * g[2];
		T const det_inv = T(1) / det;
		g00 *= det_inv;
		g01 *= det_inv;
		g02 *= det_inv;
		T const g11 = (g[0] * g[5] - g[3] * g[3]) * det_inv;
		T const g12 = (g[0] * g[4] - g[1] * g[3]) * det_inv;
		T const g22 = (g[0] * g[2] - g[1] * g[1]) * det_inv;
		ginv[0] = g00;
		ginv[1] = g01;
		ginv[2] = g11;
		ginv[3] = g02;
		ginv[4] = g12;
		ginv[5] = g22;
		return ginv;
	}
public:
	Z4SpaceTime flux(int l) {
		constexpr tensor2 sym = { { 0, 3, 5 }, { 3, 1, 4 }, { 5, 4, 2 } };
		Z4SpaceTime F = { zero, { zero, zero, zero }, { zero, zero, zero, zero, zero, zero }, { zero, zero, zero }, { { zero, zero, zero },
				{ zero, zero, zero }, { zero, zero, zero } }, { zero, zero, zero, zero, zero, zero }, { { zero, zero, zero, zero, zero, zero }, { zero, zero,
				zero, zero, zero, zero }, { zero, zero, zero, zero, zero, zero } }, zero, { zero, zero, zero } };
		tensor0 const ialpha0 = one / alpha0;
		tensorS const igamma2 = matrix_inverse(gamma2);
		tensor0 B0 = zero;
		tensor0 K0 = zero;
		tensor0 Q0 = zero;
		tensor1 D1 = { zero, zero, zero };
		tensor1 E1 = { zero, zero, zero };
		tensor1 V1 = { zero, zero, zero };
		tensor1 Q1 = { zero, zero, zero };
		tensorS Q2 = { zero, zero, zero, zero, zero, zero };
		tensorS lambda2 = { zero, zero, zero, zero, zero, zero };
		for (int i = 0; i < NDIM; i++) {
			B0 += B2[i][i];
		}
		for (int i = 0; i < NDIM; i++) {
			int const ii = sym[i][i];
			K0 += K2[ii] * igamma2[ii];
			for (int j = 0; j < i; j++) {
				int const ij = sym[i][j];
				K0 += K2[ij] * igamma2[ij] * two;
			}
		}
		Q0 = f() * (K0 - two * Theta0);
		for (int i = 0; i < NDIM; i++) {
			for (int j = 0; j < NDIM; j++) {
				for (int k = 0; k < NDIM; k++) {
					D1[i] += D3[i][sym[j][k]] * igamma2[sym[j][k]];
				}
			}
		}
		for (int i = 0; i < NDIM; i++) {
			for (int j = 0; j < NDIM; j++) {
				for (int k = 0; k < NDIM; k++) {
					E1[i] += D3[j][sym[i][k]] * igamma2[sym[k][j]];
				}
			}
		}
		for (int i = 0; i < NDIM; i++) {
			V1[i] = D1[i] - E1[i] - Z1[i];
		}
		for (int i = 0; i < NDIM; i++) {
			Q1[i] = A1[i] - D1[i] + two * V1[i];
		}
		for (int i = 0; i < NDIM; i++) {
			for (int j = 0; j <= i; j++) {
				int const ij = sym[i][j];
				for (int k = 0; k < NDIM; k++) {
					Q2[ij] = K2[ij] - half * ialpha0 * (B2[i][k] * gamma2[sym[k][j]] + B2[j][k] * gamma2[sym[k][i]]);
				}
			}
		}
		for (int i = 0; i < l; i++) {
			int const il = sym[i][l];
			lambda2[l][il] += half * ialpha0 * Q1[l];
			for (int j = 0; j <= i; j++) {
				int const ij = sym[i][j];
				for (int k = 0; k < NDIM; k++) {
					lambda2[l][ij] += igamma2[sym[l][k]] * D3[k][ij];
				}
			}
		}
		int const ll = sym[l][l];
		lambda2[l][ll] += ialpha0 * Q1[l];
		for (int j = 0; j <= l; j++) {
			int const lj = sym[l][j];
			for (int k = 0; k < NDIM; k++) {
				lambda2[l][lj] += igamma2[sym[l][k]] * D3[k][lj];
			}
		}
		for (int i = l + 1; i < NDIM; i++) {
			int const il = sym[i][l];
			lambda2[l][il] += half * ialpha0 * Q1[l];
			for (int j = 0; j <= i; j++) {
				int const ij = sym[i][j];
				for (int k = 0; k < NDIM; k++) {
					lambda2[l][ij] += igamma2[sym[l][k]] * D3[k][ij];
				}
			}
		}
		F.A1[l] += alpha0 * Q0;
		for (int i = 0; i < NDIM; i++) {
			for (int j = 0; j < NDIM; j++) {
				int const ij = sym[i][j];
				F.B2[l][i] += alpha0 * Q1[j] * igamma2[ij];
			}
		}
		for (int i = 0; i < NDIM; i++) {
			for (int j = 0; j <= i; j++) {
				int const ij = sym[i][j];
				F.D3[l][ij] += alpha0 * Q2[ij];
			}
		}
		for (int i = 0; i < NDIM; i++) {
			for (int j = 0; j <= i; j++) {
				int const ij = sym[i][j];
				F.K2[ij] = alpha0 * lambda2[l][ij];
			}
		}
		for (int i = 0; i < NDIM; i++) {
			F.Theta0 += alpha0 * V1[i] * igamma2[sym[i][l]];
		}
		for (int i = 0; i < NDIM; i++) {
			F.Z1[i] += alpha0 * (K0 - Theta0) + B0;
			for (int j = 0; j < NDIM; j++) {
				F.Z1[i] += alpha0 * K2[sym[i][j]] * igamma2[sym[j][l]] - B2[i][l];
			}
		}
		for (int i = 0; i < NDIM; i++) {
			F.A1[i] -= beta1[l] * A1[i];
		}
		for (int i = 0; i < NDIM; i++) {
			for (int j = 0; j < NDIM; j++) {
				F.B2[i][j] -= beta1[l] * B2[i][j];
			}
		}
		for (int i = 0; i < NDIM; i++) {
			for (int j = 0; j < NDIM; j++) {
				for (int k = 0; k <= j; k++) {
					int const jk = sym[j][k];
					F.D3[i][jk] -= beta1[l] * D3[i][jk];
				}
			}
		}
		for (int i = 0; i < NDIM; i++) {
			for (int j = 0; j <= i; j++) {
				int const ij = sym[i][j];
				F.K2[ij] -= beta1[l] * K2[ij];
			}
		}
		F.Theta0 -= beta1[l] * Theta0;
		for (int i = 0; i < NDIM; i++) {
			F.Z1[i] -= beta1[l] * Z1[i];
		}
		return F;
	}
	Z4SpaceTime vacuum_source() {
		constexpr tensor2 sym = { { 0, 3, 5 }, { 3, 1, 4 }, { 5, 4, 2 } };
		Z4SpaceTime S = { zero, { zero, zero, zero }, { zero, zero, zero, zero, zero, zero }, { zero, zero, zero }, { { zero, zero, zero },
				{ zero, zero, zero }, { zero, zero, zero } }, { zero, zero, zero, zero, zero, zero }, { { zero, zero, zero, zero, zero, zero }, { zero, zero,
				zero, zero, zero, zero }, { zero, zero, zero, zero, zero, zero } }, zero, { zero, zero, zero } };
		tensorS const igamma2 = matrix_inverse(gamma2);
		tensor0 B0 = zero;
		tensor0 K0 = zero;
		tensor0 Q0 = zero;
		tensor1 D1 = { zero, zero, zero };
		tensor1 E1 = { zero, zero, zero };
		tensor1 Q1 = { zero, zero, zero };
		tensor3 Gamma3 = { { zero, zero, zero, zero, zero, zero }, { zero, zero, zero, zero, zero, zero }, { zero, zero, zero, zero, zero, zero } };
		for (int i = 0; i < NDIM; i++) {
			B0 += B2[i][i];
		}
		for (int i = 0; i < NDIM; i++) {
			int const ii = sym[i][i];
			K0 += K2[ii] * igamma2[ii];
			for (int j = 0; j < i; j++) {
				int const ij = sym[i][j];
				K0 += K2[ij] * igamma2[ij] * two;
			}
		}
		Q0 = f() * (K0 - two * Theta0);
		for (int i = 0; i < NDIM; i++) {
			Q1[i] = A1[i] + D1[i] - two * (E1[i] + Z1[i]);
		}
		for (int i = 0; i < NDIM; i++) {
			for (int j = 0; j < NDIM; j++) {
				for (int k = 0; k < NDIM; k++) {
					D1[i] += D3[i][sym[j][k]] * igamma2[sym[j][k]];
				}
			}
		}
		for (int i = 0; i < NDIM; i++) {
			for (int j = 0; j < NDIM; j++) {
				for (int k = 0; k < NDIM; k++) {
					E1[i] += D3[j][sym[i][k]] * igamma2[sym[k][j]];
				}
			}
		}
		for (int i = 0; i < NDIM; i++) {
			for (int j = 0; j < NDIM; j++) {
				for (int k = 0; k < NDIM; k++) {
					int const jk = sym[j][k];
					for (int l = 0; l < NDIM; l++) {
						Gamma3[i][jk] = half * igamma2[sym[i][l]] * (D3[j][sym[l][k]] + D3[k][sym[l][j]] - D3[l][sym[j][k]]);
					}
				}
			}
		}
		for (int i = 0; i < NDIM; i++) {
			S.alpha0 += alpha0 * (beta1[i] * A1[i] - alpha0 * Q0);
		}
		for (int i = 0; i < NDIM; i++) {
			for (int j = 0; j < NDIM; j++) {
				int const ij = sym[i][j];
				S.beta1[i] += beta1[j] * B2[j][i] - alpha0 * igamma2[ij] * Q1[j];
			}
		}
		for (int i = 0; i < NDIM; i++) {
			for (int j = 0; j < NDIM; j++) {
				int const ij = sym[i][j];
				S.gamma2[ij] -= two * alpha0 * K2[ij];
				for (int k = 0; k < NDIM; k++) {
					int const ik = sym[i][k];
					int const jk = sym[j][k];
					S.gamma2[ij] += two * beta1[k] * D3[k][ij];
					S.gamma2[ij] += B2[i][k] * igamma2[jk] + B2[j][k] * igamma2[ik];
				}
			}
		}
		for (int i = 0; i < NDIM; i++) {
			S.A1[i] -= B0 * A1[i];
			for (int j = 0; j < NDIM; j++) {
				S.A1[i] += B2[i][j] * A1[j];
			}
		}
		for (int i = 0; i < NDIM; i++) {
			for (int j = 0; j < NDIM; j++) {
				S.B2[i][j] -= B0 * B2[i][j];
				for (int k = 0; k < NDIM; k++) {
					S.B2[i][j] += B2[i][k] * B2[k][j];
				}
			}
		}
		for (int i = 0; i < NDIM; i++) {
			for (int j = 0; j < NDIM; j++) {
				for (int k = 0; k < NDIM; k++) {
					int const jk = sym[i][j];
					S.D3[i][jk] -= B0 * D3[i][jk];
					for (int l = 0; l < NDIM; l++) {
						S.D3[i][jk] += B2[i][l] * D3[l][jk];
					}
				}
			}
		}
		for (int i = 0; i < NDIM; i++) {
			for (int j = 0; j <= i; j++) {
				int const ij = sym[i][j];
				S.K2[ij] -= half * alpha0 * (A1[i] * (two * E1[j] - D1[j]) + A1[j] * (two * E1[i] - D1[i]));
				S.K2[ij] -= K2[ij] * B0;
				S.K2[ij] -= alpha0 * (A1[i] * Z1[j] + A1[j] * Z1[i]);
				S.K2[ij] += alpha0 * (K0 - two * Theta0) * K2[ij];
				for (int k = 0; k < NDIM; k++) {
					int const ik = sym[i][k];
					int const jk = sym[j][k];
					S.K2[ij] += K2[ik] * B2[j][k] + K2[jk] * B2[i][k];
					for (int l = 0; l < NDIM; l++) {
						int const kl = sym[k][l];
						int const jl = sym[j][l];
						int const il = sym[i][l];
						S.K2[ij] += alpha0 * A1[k] * igamma2[kl] * D3[l][ij];
						S.K2[ij] -= two * alpha0 * igamma2[kl] * K2[il] * K2[jk];
						for (int m = 0; m < NDIM; m++) {
							int const jm = sym[j][m];
							int const km = sym[k][m];
							S.K2[ij] -= two * alpha0 * E1[k] * (D3[i][jl] + D3[j][il]) * igamma2[kl];
							S.K2[ij] += alpha0 * (D1[k] + A1[k] - two * Z1[k]) * igamma2[kl] * Gamma3[l][ij];
							for (int n = 0; n < NDIM; n++) {
								int const kn = sym[k][n];
								S.K2[ij] += two * alpha0 * igamma2[sym[n][l]] * igamma2[km] * (D3[i][kn] * D3[m][jl] + D3[j][kn] * D3[m][il]);
								S.K2[ij] -= alpha0 * igamma2[kl] * Gamma3[l][jm] * igamma2[sym[m][n]] * Gamma3[n][ik];
							}
						}
					}
				}
			}
		}
		for (int i = 0; i < NDIM; i++) {
			S.Z1[i] -= Z1[i] * B0;
			S.Z1[i] += alpha0 * A1[i] * (K0 - two * Theta0);
			for (int j = 0; j < NDIM; j++) {
				int const ij = sym[i][j];
				S.Z1[i] -= Z1[j] * B2[i][j];
				for (int k = 0; k < NDIM; k++) {
					int const jk = sym[j][k];
					int const ik = sym[i][k];
					S.Z1[i] -= alpha0 * A1[j] * igamma2[jk] * K2[ik];
					S.Z1[i] += alpha0 * igamma2[jk] * K2[ik] * (D1[j] - two * Z1[j]);
					for (int l = 0; l < NDIM; l++) {
						S.Z1[i] -= alpha0 * igamma2[jk] * K2[sym[k][l]] * Gamma3[l][ij];
					}
				}
			}
		}
		S.Theta -= Theta0 * B0;
		S.Theta += K0 * (K0 - two * Theta0);
		for (int i = 0; i < NDIM; i++) {
			for (int j = 0; j < NDIM; j++) {
				int const ij = sym[i][j];
				S.Theta += alpha0 * A1[j] * igamma2[ij] * (D1[i] - E1[i] - two * Z1[i]);
				S.Theta -= half * alpha0 * D1[j] * igamma2[ij] * (D1[i] - two * Z1[i]);
				for (int k = 0; k < NDIM; k++) {
					int const ik = sym[i][k];
					int const jk = sym[j][k];
					for (int l = 0; l < NDIM; l++) {
						int const il = sym[i][l];
						int const jl = sym[j][l];
						S.Theta -= half * alpha0 * K2[il] * igamma2[jl] * igamma2[jk] * K2[ik];
						for (int m = 0; m < NDIM; m++) {
							int const km = sym[k][m];
							int const lm = sym[l][m];
							S.Theta += half * alpha0 * igamma2[jl] * D3[i][lm] * igamma2[km] + Gamma3[i][jk];
						}
					}
				}
			}
		}
		return S;
	}
	Z4SpaceTime matter_source(tensor0 const &tau, tensor1 const &s, tensorS const &S) {
		Z4SpaceTime src;
		return src;
	}

};

void testEinstein() {
//	Z4SpaceTime<double> st;
//	st.flux(0);
//	st.vacuum_source();
//	st.matter_source(Z4SpaceTime<double>::tensor0(), Z4SpaceTime<double>::tensor1(), Z4SpaceTime<double>::tensorS());
}

