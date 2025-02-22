/*
 * Spacetime.hpp
 *
 *  Created on: Feb 21, 2025
 *      Author: dmarce1
 */

#ifndef INCLUDE_SPACETIME_HPP_
#define INCLUDE_SPACETIME_HPP_

#include "Vector.hpp"
/*
 template<typename T>
 struct GRState: public Math::Vector<T, 50> {
 static constexpr int NF = 50;
 static constexpr int alpha_index = 0;
 static constexpr int beta_index = 1;
 static constexpr int gamma_index = 4;
 static constexpr int A_index = 10;
 static constexpr int B_index = 13;
 static constexpr int D_index = 22;
 static constexpr int K_index = 40;
 static constexpr int Theta_index = 46;
 static constexpr int Z_index = 47;
 using base_type = Math::Vector<T, NF>;
 using Math::Vector<T, NDIM> = Math::Vector<T, NDIM>;
 using Tensor2 = Math::SquareMatrix<T, NDIM>;
 using SymmetricTensor2 = Math::SymmetricMatrix<T, NDIM>;
 GRState() :
 U(*reinterpret_cast<base_type*>(this)) {
 }
 GRState(GRState const &other) :
 U(other.U) {
 }
 GRState(GRState &&other) :
 U(std::move(other.U)) {
 }
 GRState& operator=(GRState const &other) {
 U = other.U;
 return *this;
 }
 GRState& operator=(GRState &&other) {
 U = std::move(other.U);
 return *this;
 }
 T& alpha() {
 return U[alpha_index];
 }
 T& beta(int i) {
 return U[beta_index + i];
 }
 T& gamma(int i, int j) {
 int const k = std::max(i, j);
 int const l = std::min(i, j);
 return U[gamma_index + ((k * (k + 1)) >> 1) + l];
 }
 T& A(int i) {
 return U[A_index + i];
 }
 T& B(int i, int j) {
 return U[B_index + NDIM * i + j];
 }
 T& D(int i, int j, int k) {
 int const l = std::max(j, k);
 int const m = std::min(j, k);
 return U[D_index + 6 * i + ((l * (l + 1)) >> 1) + m];
 }
 T& K(int i, int j) {
 int const k = std::max(i, j);
 int const l = std::min(i, j);
 return U[K_index + ((k * (k + 1)) >> 1) + l];
 }
 T& Theta() {
 return U[Theta_index];
 }
 T& Z(int i) {
 return U[beta_index + i];
 }
 T const& alpha() const {
 return U[alpha_index];
 }
 T const& beta(int i) const {
 return U[beta_index + i];
 }
 T const& gamma(int i, int j) const {
 int const k = std::max(i, j);
 int const l = std::min(i, j);
 return U[gamma_index + ((k * (k + 1)) >> 1) + l];
 }
 T const& A(int i) const {
 return U[A_index + i];
 }
 T const& B(int i, int j) const {
 return U[B_index + NDIM * i + j];
 }
 T const& D(int i, int j, int k) const {
 int const l = std::max(j, k);
 int const m = std::min(j, k);
 return U[D_index + 6 * i + ((l * (l + 1)) >> 1) + m];
 }
 T const& K(int i, int j) const {
 int const k = std::max(i, j);
 int const l = std::min(i, j);
 return U[K_index + ((k * (k + 1)) >> 1) + l];
 }
 T const& Theta() const {
 return U[Theta_index];
 }
 T const& Z(int i) const {
 return U[beta_index + i];
 }
 private:
 base_type &U;
 };
 */

template<typename T>
struct SpacetimeState: public Math::Vector<T, 50> {
	static constexpr int NF = 50;
	using base_type = Math::Vector<T, NF>;
	static constexpr T zero = T(0), half = T(0.5), one = T(1), two = T(2);
	struct params_type {
		T f;
		T zeta;
		T zeta1;
		T xi;
		T m;
	};
	static constexpr params_type params = { one, zero, zero, zero, two };
	base_type &U;
	T &alpha;
	T &Theta;
	Math::Vector<T, NDIM> &beta_u;
	Math::Vector<T, NDIM> &A_l;
	Math::Vector<T, NDIM> &Z_l;
	Math::Vector<Math::Vector<T, NDIM>, NDIM> &B_lu;
	Math::SymmetricMatrix<T, NDIM> &gamma_ll;
	Math::SymmetricMatrix<T, NDIM> &K_ll;
	Math::Vector<Math::SymmetricMatrix<T, NDIM>, NDIM> &D_lll;
public:
	static constexpr int eigenvalueCount = 5;
	static constexpr char const *fieldNames[NF] = { "alpha", "Theta", "beta_x", "beta_y", "beta_z", "A_x", "A_y", "A_z", "Z_x", "Z_y", "Z_z", "B_xx", "B_xy",
			"B_xz", "B_yx", "B_yy", "B_yz", "B_zx", "B_zy", "B_zz", "gamma_xx", "gamma_xy", "gamma_xz", "gamma_yy", "gamma_yz", "gamma_zz", "K_xx", "K_xy",
			"K_xz", "K_yy", "K_yz", "K_zz", "D_xxx", "D_xxy", "D_xxz", "D_xyy", "D_xyz", "D_xzz", "D_yxx", "D_yxy", "D_yxz", "D_yyy", "D_yyz", "D_yzz", "D_zxx",
			"D_zxy", "D_zxz", "D_zyy", "D_zyz", "D_zzz" };
	SpacetimeState() :
			U(*this), alpha(U[0]), Theta(U[1]), beta_u((Math::Vector<T, NDIM>&) U[2]), A_l((Math::Vector<T, NDIM>&) U[5]), Z_l((Math::Vector<T, NDIM>&) U[8]), B_lu(
					(Math::Vector<Math::Vector<T, NDIM>, NDIM>&) U[11]), gamma_ll((Math::SymmetricMatrix<T, NDIM>&) U[20]), K_ll(
					(Math::SymmetricMatrix<T, NDIM>&) U[26]), D_lll((Math::Vector<Math::SymmetricMatrix<T, NDIM>, NDIM>&) U[32]) {
	}
	SpacetimeState& operator=(SpacetimeState const &other) {
		base_type::operator=((base_type const&) other);
		return *this;
	}
	SpacetimeState& operator=(T const &other) {
		for (int n = 0; n < NF; n++) {
			U[n] = other;
		}
		return *this;
	}
	SpacetimeState flux(int dim) const {
		using namespace Math;
		SpacetimeState F;
		Vector<Math::SymmetricMatrix<T, NDIM>, NDIM> D_llu, D_ull;
		Vector<Vector<T, NDIM>, NDIM> B_ll, K_lu;
		Math::SymmetricMatrix<T, NDIM> gamma_uu, Q_ll, lambdak_ll;
		Vector<T, NDIM> V_l, Q_l, D_l, E_l, Q_u;
		T Q, TrB, TrK, Vk;
		gamma_uu = matrixInverse(gamma_ll);
		for (int i = 0; i < NDIM; i++) {
			D_l[i] = zero;
			for (int k = 0; k < NDIM; k++) {
				B_ll[i][k] = K_lu[i][k] = zero;
				for (int n = 0; n < NDIM; n++) {
					D_llu[i][n, k] = zero;
					D_ull[i][n, k] = zero;
					for (int m = 0; m < NDIM; m++) {
						D_llu[i][n, k] += D_lll[i][n, m] * gamma_uu[m, k];
						D_ull[i][n, k] += gamma_uu[i, m] * D_lll[m][n, k];
					}
				}
				D_l[i] += D_llu[i][k, k];
			}
		}
		for (int i = 0; i < NDIM; i++) {
			for (int k = 0; k < NDIM; k++) {
				for (int n = 0; n < NDIM; n++) {
					B_ll[k][n] += B_ll[k][i] * gamma_ll[i, n];
					K_lu[k][n] += K_ll[k, i] * gamma_uu[i, n];
				}
			}
		}
		TrB = TrK = Vk = zero;
		for (int i = 0; i < NDIM; i++) {
			E_l[i] = zero;
			TrB += B_lu[i][i];
			TrK += K_lu[i][i];
			for (int k = 0; k <= i; k++) {
				E_l[i] += D_llu[k][k, i];
				Q_ll[i, k] = K_ll[i, k] - (half / alpha) * (B_ll[i][k] + B_ll[k][i]);
			}
			V_l[i] = D_l[i] - E_l[i] - Z_l[i];
			Q_l[i] = alpha * (A_l[i] - D_l[i] + two * V_l[i]);
			Vk += V_l[i] * gamma_uu[i, dim];
		}
		Q = params.f * (TrK - params.m * Theta);
		for (int i = 0; i < NDIM; i++) {
			Q_u[i] = zero;
			for (int k = 0; k <= i; k++) {
				Q_u[i] += gamma_uu[i, k] * Q_l[k];
				lambdak_ll[i, k] = D_ull[dim][i, k] - half * (one + params.zeta) * (D_llu[i][k, dim] + D_llu[k][i, dim]);
			}
		}
		for (int i = 0; i < NDIM; i++) {
			T const tmp = half * (A_l[i] + D_l[i] + (params.zeta - one) * E_l[i]) - Z_l[i];
			lambdak_ll[dim, i] += tmp;
			lambdak_ll[i, dim] += tmp;
		}
		F.Theta = -beta_u[dim] * Theta + alpha * Vk;
		F.alpha = zero;
		for (int i = 0; i < NDIM; i++) {
			F.beta_u[i] = zero;
		}
		for (int i = 0; i < NDIM; i++) {
			for (int k = 0; k <= i; k++) {
				F.gamma_ll[i, k] = zero;
			}
		}
		for (int i = 0; i < NDIM; i++) {
			F.A_l[i] = -beta_u[dim] * A_l[i];
		}
		F.A_l[dim] += alpha * Q;
		for (int i = 0; i < NDIM; i++) {
			F.Z_l[i] = -beta_u[dim] * Z_l[i] - alpha * K_lu[i][dim] + params.zeta1 * B_lu[i][dim];
		}
		F.Z_l[dim] += alpha * (TrK - Theta) - params.zeta1 * TrB;
		for (int i = 0; i < NDIM; i++) {
			for (int k = 0; k < NDIM; k++) {
				F.B_lu[k][i] = -beta_u[dim] * B_lu[k][i];
			}
			F.B_lu[dim][i] += alpha * Q_u[i];
		}
		for (int i = 0; i < NDIM; i++) {
			for (int k = 0; k <= i; k++) {
				F.K_ll[i, k] = -beta_u[dim] * K_ll[i, k] + alpha * lambdak_ll[i, k];
			}
		}
		for (int i = 0; i < NDIM; i++) {
			for (int k = 0; k <= i; k++) {
				for (int n = 0; n < NDIM; n++) {
					F.D_lll[n][i, k] = -beta_u[dim] * D_lll[n][i, k];
				}
				F.D_lll[dim][i, k] += alpha * Q_ll[i, k];
			}
		}
		return F;
	}
	SpacetimeState source(Math::SymmetricMatrix<T, NDIM> const &S_ll, Math::Vector<T, NDIM> const &S_l, T tau) const {
		using namespace Math;
		static constexpr T eightPi = T(8.0 * M_PI);
		auto const gamma_uu = matrixInverse(gamma_ll);
		auto S = source(gamma_uu, true);
		T TrS = zero;
		for (int i = 0; i < NDIM; i++) {
			for (int k = 0; k < NDIM; k++) {
				TrS += gamma_uu[k, i] * S_ll[i, k];
			}
		}
		S.Theta -= eightPi * alpha * tau;
		for (int i = 0; i < NDIM; i++) {
			S.Z_l[i] -= eightPi * alpha * S_l[i];
			for (int k = 0; k <= i; k++) {
				S.K_ll[i, k] -= eightPi * alpha * (S_ll[i, k] - half * (TrS - tau) * gamma_ll[i, k]);
			}
		}
		return S;
	}
	SpacetimeState source(Math::SymmetricMatrix<T, NDIM> gamma_uu = zero, bool has_gamma_uu = false) const {
		using namespace Math;
		SpacetimeState S;
		Vector<Math::SymmetricMatrix<T, NDIM>, NDIM> D_luu, D_ull, Gamma_ull;
		Vector<Vector<Vector<T, NDIM>, NDIM>, NDIM> D_llu;
		Vector<Vector<T, NDIM>, NDIM> B_ll, K_lu, K_ul;
		Math::SymmetricMatrix<T, NDIM> Q_ll;
		Vector<T, NDIM> V_l, Q_l, D_l, E_l, D_u, E_u, Z_u, Q_u;
		T Q, TrB, TrK;
		if (!has_gamma_uu) {
			gamma_uu = matrixInverse(gamma_ll);
		}
		for (int i = 0; i < NDIM; i++) {
			for (int k = 0; k < NDIM; k++) {
				D_l[i] = B_ll[i][k] = K_lu[i][k] = K_ul[i][k] = zero;
				for (int n = 0; n < NDIM; n++) {
					D_llu[i][n][k] = zero;
					B_ll[i][k] += B_lu[i][n] * gamma_ll[n, k];
					K_lu[i][k] += K_ll[i, n] * gamma_uu[n, k];
					K_ul[i][k] += gamma_uu[i, n] * K_ll[n, k];
					for (int m = 0; m < NDIM; m++) {
						D_llu[i][n][k] += D_lll[i][n, m] * gamma_uu[m, k];
					}
				}
				for (int n = 0; n <= k; n++) {
					Gamma_ull[i][n, k] = D_ull[i][n, k] = D_luu[i][n, k] = zero;
					for (int m = 0; m < NDIM; m++) {
						Gamma_ull[i][n, k] += half * gamma_uu[i, n] * (D_lll[n][k, m] + D_lll[k][m, n] - D_lll[m][k, n]);
						D_ull[i][n, k] += gamma_uu[i, m] * D_lll[m][n, k];
						D_luu[i][n, k] += D_llu[i][m][k] * gamma_uu[m, n];
					}
				}
				D_l[i] += D_llu[i][k][k];
			}
		}
		TrB = TrK = zero;
		for (int i = 0; i < NDIM; i++) {
			TrB += B_lu[i][i];
			TrK += K_lu[i][i];
			E_l[i] = zero;
			for (int k = 0; k <= i; k++) {
				E_l[i] += D_llu[k][k, i];
				Q_ll[i, k] = K_ll[i, k] - (half / alpha) * (B_ll[i][k] + B_ll[k][i]);
			}
			V_l[i] = D_l[i] - E_l[i] - Z_l[i];
			Q_l[i] = alpha * (A_l[i] - D_l[i] + two * V_l[i]);
		}
		Q = params.f * (TrK - params.m * Theta);
		for (int i = 0; i < NDIM; i++) {
			Q_u[i] = D_u[i] = E_u[i] = Z_u[i] = zero;
			for (int k = 0; k < NDIM; k++) {
				Q_u[i] += gamma_uu[i, k] * Q_l[k];
				D_u[i] += gamma_uu[i, k] * D_l[k];
				E_u[i] += gamma_uu[i, k] * E_l[k];
				Z_u[i] += gamma_uu[i, k] * Z_l[k];
			}
		}

		S.alpha = -nSquared(alpha) * Q;
		for (int i = 0; i < NDIM; i++) {
			S.alpha += alpha * beta_u[i] * A_l[i];
		}
		for (int i = 0; i < NDIM; i++) {
			S.beta_u[i] = -alpha * Q_u[i];
			for (int k = 0; k < NDIM; k++) {
				S.beta_u[i] += beta_u[k] * B_lu[k][i];
			}
		}
		for (int i = 0; i < NDIM; i++) {
			for (int k = 0; k <= i; k++) {
				S.gamma_ll[i, k] = B_lu[i][k] + B_lu[k][i] - two * alpha * K_ll[k, i];
				for (int n = 0; n < NDIM; n++) {
					S.gamma_ll[i, k] += two * beta_u[n] * D_lll[n][i, k];
				}
			}
		}
		for (int i = 0; i < NDIM; i++) {
			S.A_l[i] = -TrB * A_l[i];
			for (int k = 0; k <= i; k++) {
				S.A_l[i] += B_lu[i][k] * A_l[k];
			}
		}
		for (int i = 0; i < NDIM; i++) {
			for (int k = 0; k < NDIM; k++) {
				S.B_lu[i][k] = -TrB * B_lu[i][k];
				for (int n = 0; n < NDIM; n++) {
					S.B_lu[i][k] += B_lu[i][n] * B_lu[n][k];
				}
			}
		}
		for (int i = 0; i < NDIM; i++) {
			for (int k = 0; k < NDIM; k++) {
				for (int n = 0; n <= k; n++) {
					S.D_lll[i][k, n] = -TrB * D_lll[i][k, n];
					for (int m = 0; m < NDIM; m++) {
						S.D_lll[i][k, n] += B_lu[i][m] * D_lll[m][k, n];
					}
				}
			}
		}
		for (int i = 0; i < NDIM; i++) {
			for (int k = 0; k <= i; k++) {
				T K0, K1, K2, K3;
				K0 = -TrB * K_ll[i, k];
				K1 = half * (A_l[i] * D_l[k] + A_l[k] * D_l[i]);
				K2 = -A_l[i] * (E_l[k] - half * D_l[k]) + A_l[k] * (E_l[i] - half * D_l[i]);
				K3 = -A_l[i] * Z_l[k] - A_l[k] * Z_l[i] + (TrK - Theta) * K_ll[i, k];
				for (int n = 0; n < NDIM; n++) {
					K0 += B_lu[i][n] * K_ll[n, k];
					K0 += B_lu[k][n] * K_ll[n, i];
					K1 -= A_l[n] * Gamma_ull[n][i, k];
					K2 += A_l[n] * D_ull[n][i, k];
					K2 -= two * (D_llu[i][k, n] + D_llu[k][i, n]) * E_l[n];
					K3 += (D_l[n] + A_l[n] - two * Z_l[n]) * Gamma_ull[n][i, k];
					K3 -= K_ul[n][i] * K_ll[n, k];
					for (int m = 0; m < NDIM; m++) {
						K2 += two * (D_llu[i][n, m] * D_ull[n][m, k] + D_llu[k][n, m] * D_ull[n][m, i]);
						K3 -= Gamma_ull[n][m, k] * Gamma_ull[m][n, i];
					}
				}
				S.K_ll[i, k] = K0 + alpha * (half * ((one + params.xi) * K1 + (one - params.xi) * K2) + K3);
			}
		}
		for (int i = 0; i < NDIM; i++) {
			T Z0;
			S.Z_l[i] = -TrB * Z_l[i];
			Z0 = A_l[i] * (TrK - two * Theta);
			for (int k = 0; k < NDIM; k++) {
				S.Z_l[i] += Z_l[k] * B_lu[i][k];
				Z0 += K_ul[k][i] * (D_l[k] - A_l[k] - two * Z_l[k]);
				for (int n = 0; n < NDIM; n++) {
					Z0 -= K_ul[k][n] * Gamma_ull[n][k, i];
				}
			}
			S.Z_l[i] += alpha * Z0;
		}
		T Theta0 = TrK * (TrK - two * Theta);
		for (int i = 0; i < NDIM; i++) {
			Theta0 += two * A_l[i] * (D_u[i] - E_u[i] - two * Z_u[i]);
			Theta0 -= D_u[i] * (D_l[i] - two * Z_l[i]);
			for (int k = 0; k < NDIM; k++) {
				Theta0 -= K_ul[i][k] * K_ul[k][i];
				for (int n = 0; n < NDIM; n++) {
					Theta0 += D_luu[i][k, n] * Gamma_ull[i][k, n];
				}
			}
		}
		S.Theta = half * alpha * Theta0 - TrB * Theta;
		return S;
	}
	Math::Vector<T, eigenvalueCount> eigenvalues(int dim) const {
		Math::Vector<T, eigenvalueCount> lambda;
		lambda[0] = -params.f * alpha - beta_u[dim];
		lambda[1] = -alpha - beta_u[dim];
		lambda[2] = -beta_u[dim];
		lambda[3] = +alpha - beta_u[dim];
		lambda[4] = +params.f * alpha - beta_u[dim];
		return lambda;
	}
	friend SpacetimeState riemannFlux(SpacetimeState const &UL, SpacetimeState const &UR, int dim) {
		SpacetimeState F;
		SpacetimeState const FR = UR.flux(dim);
		SpacetimeState const FL = UL.flux(dim);
		auto const lambdaL = UL.eigenvalues(dim);
		auto const lambdaR = UR.eigenvalues(dim);
		T sL = zero, sR = zero;
		for (int i = 0; i < eigenvalueCount; i++) {
			sL = min(sL, min(lambdaL[i], lambdaR[i]));
			sR = max(sR, max(lambdaL[i], lambdaR[i]));
		}
		for (int n = 0; n < NF; n++) {
			F[n] = (sR * FL[n] - sL * FR[n] + sR * sL * (UR[n] - UL[n])) / (sR - sL);
		}
		return F;
	}
};

#endif /* INCLUDE_SPACETIME_HPP_ */
