/*
 * Einstein.cpp
 *
 *  Created on: Feb 25, 2025
 *      Author: dmarce1
 */

#include <algorithm>
#include <functional>

#include "SparseMatrix.hpp"
#include "Vector.hpp"

static constexpr int NDIM2 = NDIM * (NDIM + 1) / 2;

using T = double;

constexpr int symmtetricMap(int k, int j) {
	static constexpr int map[NDIM][NDIM] = { { 0, 1, 3 }, { 1, 2, 4 }, { 3, 4, 5 } };
	return map[k][j];
}

template<typename >
class Z4Metric;

template<typename >
struct Z4DynamicVariables;

template<typename T>
struct Z4Metric: public Math::Vector<T, 10> {
	static constexpr int NF = 10;
	static constexpr T zero = T(0);
	static constexpr T one = T(1);
	T alpha() const {
		return (*this)[alphaIndex];
	}
	T beta_I(int i) const {
		return (*this)[betaIndex + i];
	}
	T gamma_ij(int i, int j) const {
		const int k = symmtetricMap(i, j);
		return (*this)[gammaIndex + k];
	}
	T& alpha() {
		return (*this)[alphaIndex];
	}
	T& beta_I(int i) {
		return (*this)[betaIndex + i];
	}
	T& gamma_ij(int i, int j) {
		const int k = symmtetricMap(i, j);
		return (*this)[gammaIndex + k];
	}
	void reset() {
		std::fill(this->begin(), this->end(), 0);
	}
	static Z4Metric genMinkowskiMetric() {
		Z4Metric M;
		M.alpha() = one;
		for (int i = 0; i < NDIM; i++) {
			M.beta_I(i) = zero;
			M.gamma_ij(i, i) = one;
			for (int j = i + 1; j < NDIM; j++) {
				M.gamma_ij(i, j) = zero;
			}
		}
		return M;
	}
	auto genIndexLower() const {
		Math::Vector<T, NF> const gam { *this };
		return [gam](int j, int k) {
			return gam[symmtetricMap(j, k)];
		};
	}
	auto genIndexRaise() const {
		T const *gam = &(*this)[gammaIndex];
		Math::Vector<T, NF> igam;
		igam[0] = gam[2] * gam[5] - gam[4] * gam[4];
		igam[1] = gam[4] * gam[3] - gam[1] * gam[5];
		igam[3] = gam[1] * gam[4] - gam[2] * gam[3];
		T const idet = T(1) / (gam[0] * igam[0] + gam[1] * igam[1] + gam[3] * igam[3]);
		igam[0] *= idet;
		igam[1] *= idet;
		igam[2] = (gam[0] * gam[5] - gam[3] * gam[3]) * idet;
		igam[3] *= idet;
		igam[4] = (gam[1] * gam[3] - gam[0] * gam[4]) * idet;
		igam[5] = (gam[0] * gam[2] - gam[1] * gam[1]) * idet;
		return [igam](int j, int k) {
			return igam[symmtetricMap(j, k)];
		};
	}
private:
	static constexpr int alphaIndex = 0;
	static constexpr int betaIndex = 1;
	static constexpr int gammaIndex = 4;
};

template<typename T>
struct Z4DynamicVariables: public Math::Vector<T, 40> {
	static constexpr int NF = 40;
	static constexpr int xi = 0;
	static constexpr int zeta = 0;
	static constexpr int zeta1 = 0;
	static constexpr T half = T(0.5);
	static constexpr T one = T(1);
	static constexpr T two = T(2);
	T K_ij(int i, int j) const {
		const int k = symmtetricMap(i, j);
		return (*this)[indexK + k];
	}
	T A_i(int i) const {
		return (*this)[indexA + i];
	}
	T B_iJ(int i, int j) const {
		return (*this)[indexB + NDIM * i + j];
	}
	T D_kij(int k, int i, int j) const {
		const int m = symmtetricMap(i, j);
		return (*this)[indexD + NDIM2 * k + m];
	}
	T Theta() const {
		return (*this)[indexTheta];
	}
	T Z_i(int i) const {
		return (*this)[indexZ + i];
	}
	T& K_ij(int i, int j) {
		const int k = symmtetricMap(i, j);
		return (*this)[indexK + k];
	}
	T& A_i(int i) {
		return (*this)[indexA + i];
	}
	T& B_iJ(int i, int j) {
		return (*this)[indexB + NDIM * i + j];
	}
	T& D_kij(int k, int i, int j) {
		const int m = symmtetricMap(i, j);
		return (*this)[indexD + NDIM2 * k + m];
	}
	T& Theta() {
		return (*this)[indexTheta];
	}
	T& Z_i(int i) {
		return (*this)[indexZ + i];
	}
	void reset() {
		std::fill(this->begin(), this->end(), 0);
	}

	Z4DynamicVariables<T> flux(Z4Metric<T> const &V, int k) const;

	static constexpr bool isTransverseField(int fieldNumber, int dim) {
		static constexpr bool flags[NDIM][NF] = { {
		/*Kxx*/false,
		/*Kxy*/false,
		/*Kyy*/false,
		/*Kxz*/false,
		/*Kyz*/false,
		/*Kzz*/false,
		/*Ax*/false,
		/*Ay*/true,
		/*Az*/true,
		/*Bxx*/false,
		/*Bxy*/false,
		/*Bxz*/false,
		/*Byx*/true,
		/*Byy*/true,
		/*Byz*/true,
		/*Bzx*/true,
		/*Bzy*/true,
		/*Bzz*/true,
		/*Dxxx*/false,
		/*Dxxy*/false,
		/*Dxyy*/false,
		/*Dxxz*/false,
		/*Dxyz*/false,
		/*Dxzz*/false,
		/*Dyxx*/true,
		/*Dyxy*/true,
		/*Dyyy*/true,
		/*Dyxz*/true,
		/*Dyyz*/true,
		/*Dyzz*/true,
		/*Dzxx*/true,
		/*Dzxy*/true,
		/*Dzyy*/true,
		/*Dzxz*/true,
		/*Dzyz*/true,
		/*Dzzz*/true,
		/*Theta*/false,
		/*Zx*/false,
		/*Zy*/false,
		/*Zz*/false }, {
		/*Kxx*/false,
		/*Kxy*/false,
		/*Kyy*/false,
		/*Kxz*/false,
		/*Kyz*/false,
		/*Kzz*/false,
		/*Ax*/true,
		/*Ay*/false,
		/*Az*/true,
		/*Bxx*/true,
		/*Bxy*/true,
		/*Bxz*/true,
		/*Byx*/false,
		/*Byy*/false,
		/*Byz*/false,
		/*Bzx*/true,
		/*Bzy*/true,
		/*Bzz*/true,
		/*Dxxx*/true,
		/*Dxxy*/true,
		/*Dxyy*/true,
		/*Dxxz*/true,
		/*Dxyz*/true,
		/*Dxzz*/true,
		/*Dyxx*/false,
		/*Dyxy*/false,
		/*Dyyy*/false,
		/*Dyxz*/false,
		/*Dyyz*/false,
		/*Dyzz*/false,
		/*Dzxx*/true,
		/*Dzxy*/true,
		/*Dzyy*/true,
		/*Dzxz*/true,
		/*Dzyz*/true,
		/*Dzzz*/true,
		/*Theta*/false,
		/*Zx*/false,
		/*Zy*/false,
		/*Zz*/false }, {
		/*Kxx*/false,
		/*Kxy*/false,
		/*Kyy*/false,
		/*Kxz*/false,
		/*Kyz*/false,
		/*Kzz*/false,
		/*Ax*/true,
		/*Ay*/true,
		/*Az*/false,
		/*Bxx*/true,
		/*Bxy*/true,
		/*Bxz*/true,
		/*Byx*/true,
		/*Byy*/true,
		/*Byz*/true,
		/*Bzx*/false,
		/*Bzy*/false,
		/*Bzz*/false,
		/*Dxxx*/true,
		/*Dxxy*/true,
		/*Dxyy*/true,
		/*Dxxz*/true,
		/*Dxyz*/true,
		/*Dxzz*/true,
		/*Dyxx*/true,
		/*Dyxy*/true,
		/*Dyyy*/true,
		/*Dyxz*/true,
		/*Dyyz*/true,
		/*Dyzz*/true,
		/*Dzxx*/false,
		/*Dzxy*/false,
		/*Dzyy*/false,
		/*Dzxz*/false,
		/*Dzyz*/false,
		/*Dzzz*/false,
		/*Theta*/false,
		/*Zx*/false,
		/*Zy*/false,
		/*Zz*/false } };
		return flags[fieldNumber][dim];
	}
private:
	static constexpr int indexK = 0;
	static constexpr int indexA = 6;
	static constexpr int indexB = 9;
	static constexpr int indexD = 18;
	static constexpr int indexTheta = 36;
	static constexpr int indexZ = 37;
};

template<typename T>
struct Z4AuxiliaryVariables: public Math::Vector<T, 11> {
	static constexpr int NF = 11;
	Z4AuxiliaryVariables(std::function<T(int, int)> const &raiseIndex, Z4DynamicVariables<T> const &st);
	T D_i(int i) const {
		return (*this)[indexD + i];
	}
	T E_i(int i) const {
		return (*this)[indexE + i];
	}
	T V_i(int i) const {
		return (*this)[indexV + i];
	}
	T B() const {
		return (*this)[indexB];
	}
	T K() const {
		return (*this)[indexK];
	}
private:
	void reset() {
		std::fill(this->begin(), this->end(), 0);
	}
	static constexpr int indexK = 0;
	static constexpr int indexB = 1;
	static constexpr int indexD = 2;
	static constexpr int indexE = 5;
	static constexpr int indexV = 8;
};

template<typename T>
Z4AuxiliaryVariables<T>::Z4AuxiliaryVariables(std::function<T(int, int)> const &gRaise, Z4DynamicVariables<T> const &U) {
	T *D1 = &(*this)[indexD];
	T *E1 = &(*this)[indexE];
	T *V1 = &(*this)[indexV];
	T *K0 = &(*this)[indexK];
	T *B0 = &(*this)[indexB];
	reset();
	for (int i = 0; i < NDIM; i++) {
		B0[0] += U.B_iJ(i, i);
		for (int j = 0; j <= i; j++) {
			T w = 1.0 + T(i != j);
			K0[0] += w * U.K_ij(i, j) * gRaise(j, i);
			for (int k = 0; k < NDIM; k++) {
				D1[k] += w * U.D_kij(k, i, j) * gRaise(j, i);
				E1[k] += U.D_kij(i, k, j) * gRaise(j, i);
			}
		}
		for (int j = i + 1; j < NDIM; j++) {
			for (int k = 0; k < NDIM; k++) {
				E1[k] += U.D_kij(i, k, j) * gRaise(j, i);
			}
		}
	}
	for (int i = 0; i < NDIM; i++) {
		V1[i] = D1[i] - E1[i] - U.Z_i(i);
	}
}

template<typename T>
Z4DynamicVariables<T> Z4DynamicVariables<T>::flux(Z4Metric<T> const &V, int k) const {
	static constexpr T zero = T(0);
	auto const lower = V.genIndexLower();
	auto const raise = V.genIndexRaise();
	Z4DynamicVariables<T> F;
	F.reset();
	T const ahpla = one / V.alpha();
	Math::Vector<T, NDIM> D_i(zero);
	Math::Vector<T, NDIM> E_i(zero);
	Math::Vector<T, NDIM> V_i(zero);
	T K0 = zero;
	T B0 = zero;
	for (int i = 0; i < NDIM; i++) {
		B0 += B_iJ(i, i);
		for (int j = 0; j <= i; j++) {
			T w = 1.0 + T(i != j);
			K0 += w * K_ij(i, j) * raise(j, i);
			for (int k = 0; k < NDIM; k++) {
				D_i[k] += w * D_kij(k, i, j) * raise(j, i);
				E_i[k] += D_kij(i, k, j) * raise(j, i);
			}
		}
		for (int j = i + 1; j < NDIM; j++) {
			for (int k = 0; k < NDIM; k++) {
				E_i[k] += D_kij(i, k, j) * raise(j, i);
			}
		}
	}
	for (int i = 0; i < NDIM; i++) {
		V_i[i] = D_i[i] - E_i[i] - Z_i(i);
	}

	for (int i = 0; i < NDIM; i++) {
		T const wt = T(i == k) + T(1);
		F.Theta() += raise(k, i) * V_i[i];
	}

	F.A_i(k) = K0 - two * Theta();

	F.Z_i(k) += K0 - Theta();
	F.Z_i(k) -= zeta1 * B0;
	for (int i = 0; i < NDIM; i++) {
		T const wt = T(i == k) + T(1);
		F.Z_i(i) += zeta1 * B_iJ(i, k);
		for (int j = 0; j < NDIM; j++) {
			F.Z_i(i) -= raise(k, j) * K_ij(j, i);
		}
	}

	for (int i = 0; i < NDIM; i++) {
		T const wt = T(i == k) + T(1);
		F.K_ij(k, i) += half * wt * (A_i(i) + D_i[i] - two * Z_i(i) - (one - zeta) * E_i[i]);
		for (int j = 0; j <= i; j++) {
			for (int m = 0; m < NDIM; m++) {
				F.K_ij(i, j) += raise(k, m) * (D_kij(m, i, j) - half * (one + zeta) * (D_kij(j, i, m) + D_kij(i, j, m)));
			}
		}
	}

	for (int i = 0; i < NDIM; i++) {
		for (int j = 0; j < NDIM; j++) {
			F.B_iJ(k, i) += raise(i, j) * (A_i(j) - D_i[j] + two * V_i[j]);
		}
	}

	for (int i = 0; i < NDIM; i++) {
		for (int j = 0; j <= i; j++) {
			F.D_kij(k, i, j) += K_ij(i, j);
			for (int m = 0; m < NDIM; m++) {
				F.D_kij(k, i, j) -= half * ahpla * (lower(j, m) * B_iJ(i, m) + lower(i, m) * B_iJ(j, m));
			}
		}
	}

	auto itF = F.begin();
	auto itU = this->begin();
	for (; itU != this->end(); itU++, itF++) {
		auto &f = *itF;
		auto &u = *itU;
		*itF = f * V.alpha() - V.beta_I(k) * u;
	}
	return F;

}

template<typename T>
Z4DynamicVariables<T> Z4RiemannSolver(Z4DynamicVariables<T> &UR, Z4DynamicVariables<T> const &UL, Z4Metric<T> const &g, int dim) {
	using namespace std;
	auto const transverseField = Z4DynamicVariables<T>::isTransverseField;
	static constexpr int NF = Z4DynamicVariables<T>::NF;
	static constexpr T zero = T(0);
	static constexpr T half = T(0.5);
	static constexpr T one = T(1);
	Z4DynamicVariables<T> F;
	T const alpha = g.alpha();
	T const beta = g.beta(dim);
	T const sR = max(-beta + alpha, zero);
	T const sL = min(-beta - alpha, zero);
	T const wInv = one / (sR - sL);
	T const wR = -sL * wInv;
	T const wL = +sR * wInv;
	T const c0 = -sR * sL * wInv;
	if ((sL < zero) && (zero < sR)) {
		auto const FL = UL.flux(dim);
		auto const FR = UR.flux(dim);
		auto const &F0 = (sL < zero) ? FR : FL;
		for (int fi = 0; fi < NF; fi++) {
			if (transverseField(fi, dim)) {
				F[fi] = F0[fi];
			} else {
				F[fi] = (wL * FL[fi] + wR * FR[fi]) - c0 * (UR[fi] - UL[fi]);
			}
		}
	} else if (zero <= sL) {
		F = UL.flux(dim);
	} else /*if (sR <= zero)*/{
		F = UR.flux(dim);
	}
	return F;
}

struct Z4source_type {
	Z4Metric<T> g;
	Z4DynamicVariables<T> u;
};

template<typename Type>
auto Z4StressEnergySource(Z4Metric<Type> const &g, Math::SymmetricMatrix<Type, NDIM + 1> T) {
	static constexpr Type zero = Type(0);
	static constexpr Type half = Type(0.5);
	Type constexpr g0 = 8.0 * M_PI;
	auto const lower = g.genIndexLower();
	auto const raise = g.genIndexRaise();
	Z4source_type S;
	S.g.reset();
	S.u.reset();
	Type TrS = zero;
	for (int i = 0; i < NDIM; i++) {
		for (int j = 0; j < NDIM; j++) {
			TrS += raise(j, i) * T[i + 1, j + 1];
		}
	}
	S.u.Theta() -= g0 * g.alpha() * T[0, 0];
	for (int i = 0; i < NDIM; i++) {
		S.u.Z_i(i) -= g0 * g.alpha() * T[0, i + 1];
		for (int j = 0; j <= i; j++) {
			S.u.K_ij(i, j) -= g0 * g.alpha() * (T[i + 1, j + 1] - half * (TrS - T[0, 0]) * lower(i, j));
		}
	}
	return S;
}

template<typename T>
auto Z4VacuumSource(Z4Metric<T> const &g, Z4DynamicVariables<T> const &u) {
	static constexpr T zero = T(0);
	static constexpr T half = T(0.5);
	static constexpr T one = T(1);
	static constexpr T two = T(2);
	static constexpr T xi = Z4DynamicVariables<T>::xi;
	auto const lower = g.genIndexLower();
	auto const raise = g.genIndexRaise();
	Z4source_type S;
	S.g.reset();
	S.u.reset();
	Math::Vector<T, NDIM> D_i(zero);
	Math::Vector<T, NDIM> E_i(zero);
	Math::Vector<T, NDIM> V_i(zero);
	Math::Vector<Math::SymmetricMatrix<T, NDIM>, NDIM> Gamma_Kij(zero);
	Math::Vector<Math::SymmetricMatrix<T, NDIM>, NDIM> D_Kij(zero);
	Math::Vector<Math::SquareMatrix<T, NDIM>, NDIM> D_kiJ(zero);
	Math::SquareMatrix<T, NDIM> K_Ij(zero);

	T K0 = zero;
	T B0 = zero;
	for (int i = 0; i < NDIM; i++) {
		B0 += u.B_iJ(i, i);
		for (int j = 0; j <= i; j++) {
			T w = 1.0 + T(i != j);
			K0 += w * u.K_ij(i, j) * raise(j, i);
			for (int k = 0; k < NDIM; k++) {
				K_Ij[i, j] += raise(i, k) * u.K_ij(k, j);
			}
		}
		for (int j = 0; j < NDIM; j++) {
			for (int k = 0; k < NDIM; k++) {
				K_Ij[i, j] += raise(i, k) * u.K_ij(k, j);
				E_i[k] += u.D_kij(i, k, j) * raise(j, i);
			}
		}
	}
	for (int k = 0; k < NDIM; k++) {
		V_i[k] = D_i[k] - E_i[k] - u.Z_i(k);
		for (int j = 0; j < NDIM; j++) {
			for (int i = 0; i <= j; i++) {
				for (int m = 0; m < NDIM; m++) {
					Gamma_Kij[k][i, j] += half * raise(k, m) * (u.D_kij(i, m, j) + u.D_kij(j, i, m) - u.D_kij(m, i, j));
					D_Kij[k][i, j] += raise(k, m) * u.D_kij(m, i, j);
				}
			}
			for (int i = 0; i < NDIM; i++) {
				for (int m = 0; m < NDIM; m++) {
					D_kiJ[k][i, j] += raise(j, m) * u.D_kij(k, i, m);
				}
			}
		}
	}
	T const ahpla = one / g.alpha();
	T const alpha2 = g.alpha() * g.alpha();
	auto &Salpha = S.g.alpha();
	Salpha -= alpha2 * (K0 - two * u.Theta());
	for (int k = 0; k < NDIM; k++) {
		auto &Sbeta = S.g.beta_I(k);
		Salpha += g.alpha() * g.beta_I(k) * u.A_i(k);
		for (int j = 0; j <= k; j++) {
			Sbeta += g.beta_I(j) * u.B_iJ(j, k);
			for (int m = 0; m < NDIM; m++) {
				Sbeta -= alpha2 * raise(j, m) * (u.A_i(m) + D_i[m] - two * (E_i[m] + u.Z_i(m)));
			}
		}
		for (int j = 0; j < NDIM; j++) {
			auto &Sgamma = S.g.gamma_ij(j, k);
			Sgamma -= two * g.alpha() * u.K_ij(j, k);
			for (int m = 0; m < NDIM; m++) {
				Sgamma += two * g.beta_I(m) * u.D_kij(m, j, k);
				Sgamma += u.B_iJ(j, m) * lower(m, k) + u.B_iJ(k, m) * lower(m, j);
			}
		}
	}

	Salpha -= alpha2 * (K0 - two * u.Theta());
	for (int k = 0; k < NDIM; k++) {
		Salpha += g.alpha() * g.beta_I(k) * u.A_i(k);
	}

	for (int k = 0; k < NDIM; k++) {
		auto &Sbeta = S.g.beta_I(k);
		Salpha += g.alpha() * g.beta_I(k) * u.A_i(k);
		for (int j = 0; j < NDIM; j++) {
			Sbeta -= alpha2 * raise(j, j) * (u.A_i(j) + D_i[j] - two * (E_i[j] + u.Z_i(j)));
		}
	}

	for (int k = 0; k < NDIM; k++) {
		for (int j = 0; j < NDIM; j++) {
			auto &Sgamma = S.g.gamma_ij(j, k);
			Sgamma -= two * g.alpha() * u.K_ij(j, k);
			for (int m = 0; m < NDIM; m++) {
				Sgamma += two * g.beta_I(m) * u.D_kij(m, j, k);
				Sgamma += u.B_iJ(j, m) * lower(m, k) + u.B_iJ(k, m) * lower(m, j);
			}
		}
	}

	for (int k = 0; k < NDIM; k++) {
		for (int j = 0; j < NDIM; j++) {
			S.u.A_i(k) += u.B_iJ(k, j) * u.A_i(j) - B0 * u.A_i(k);
		}
	}

	for (int k = 0; k < NDIM; k++) {
		for (int i = 0; i < NDIM; i++) {
			for (int j = 0; j < NDIM; j++) {
				S.u.B_iJ(k, i) += u.B_iJ(k, j) * u.B_iJ(j, i) - B0 * u.B_iJ(k, i);
			}
		}
	}

	for (int k = 0; k < NDIM; k++) {
		for (int i = 0; i < NDIM; i++) {
			for (int j = 0; j <= i; j++) {
				for (int m = 0; m < NDIM; m++) {
					S.u.D_kij(k, i, j) += u.B_iJ(k, m) * u.D_kij(m, i, j) - B0 * u.D_kij(k, i, j);
				}
			}
		}
	}

	constexpr T c0 = half * (one + xi);
	constexpr T c1 = half * (one - xi);
	for (int i = 0; i < NDIM; i++) {
		for (int j = 0; j <= i; j++) {
			auto &SK = S.u.K_ij(i, j);
			SK += u.K_ij(i, j) * (K0 - two * u.Theta());
			SK -= u.A_i(i) * u.Z_i(j) + u.A_i(j) * u.Z_i(i);
			SK += u.A_i(i) * D_i[j] + u.A_i(j) * D_i[i];
			SK -= c1 * u.A_i(i) * E_i[j] + u.A_i(j) * E_i[i];
			for (int k = 0; k < NDIM; k++) {
				SK += (D_i[k] + u.A_i(k) - two * u.Z_i(k)) * Gamma_Kij[k][i, j];
				SK -= c0 * u.A_i(k) * Gamma_Kij[k][i, j];
				SK += c1 * u.A_i(k) * u.D_kij(k, i, j);
				SK -= c1 * two * E_i[k] * (D_kiJ[i][j, k] + D_kiJ[j][i, k]);
				for (int m = 0; m < NDIM; m++) {
					SK += c1 * two * (D_kiJ[i][k, m] * D_Kij[k][m, j] + D_kiJ[j][k, m] * D_Kij[k][m, i]);
					SK -= two * raise(k, m) * u.K_ij(m, i) * u.K_ij(k, j);
					SK -= Gamma_Kij[k][m, j] * Gamma_Kij[m][k, i];
				}
			}
			SK = g.alpha() * SK - u.K_ij(i, j) * B0;
			for (int k = 0; k < NDIM; k++) {
				SK += u.K_ij(i, k) * u.B_iJ(j, k) + u.K_ij(j, k) * u.B_iJ(i, k);
			}
		}
	}

	for (int i = 0; i < NDIM; i++) {
		auto &SZ = S.u.Z_i(i);
		SZ += u.A_i(i) * (K0 - two * u.Theta());
		for (int k = 0; k < NDIM; k++) {
			SZ += (D_i[k] - u.A_i(k) - two * u.Z_i(k)) * K_Ij[k, i];
			for (int m = 0; m < NDIM; m++) {
				SZ -= K_Ij[k, m] * Gamma_Kij[m][k, i];
			}
			SZ = g.alpha() * SZ - B0 * u.Z_i(i);
			for (int k = 0; k < NDIM; k++) {
				SZ += u.Z_i(k) * u.B_iJ(i, k);
			}
		}

	}

	auto &STheta = S.u.Theta();
	STheta += K0 * (K0 - two * u.Theta());
	for (int k = 0; k < NDIM; k++) {
		for (int r = 0; r < NDIM; r++) {
			STheta += raise(k, r) * (u.A_i(k) * (D_i[r] - E_i[r] - two * u.Z_i(r)) - half * D_i[r] * (D_i[r] - two * u.Z_i(r)));
			STheta -= K_Ij[k, r] * K_Ij[r, k];
			for (int s = 0; s < NDIM; s++) {
				for (int m = 0; m < NDIM; m++) {
					STheta += half * u.D_kij(k, m, s) * raise(r, m) * Gamma_Kij[k][r, s];
				}
			}
		}
	}
	STheta = g.alpha() * STheta - B0 * u.Theta();

	return S;
}

#include <iostream>
static constexpr int indexK = 0;
static constexpr int indexA = 6;
static constexpr int indexB = 9;
static constexpr int indexD = 18;
static constexpr int indexTheta = 36;
static constexpr int indexZ = 37;

template<int D>
auto TensorMapping(std::vector<int> outputIndices, std::vector<int> inputIndices) {
	int const R = outputIndices.size();
	int const size = std::pow(D, R);
	std::vector<int> oMap(size);
	std::vector<int> iMap(size);
	std::vector<int> oStrides(R);
	std::vector<int> iStrides(R);
	oStrides[0] = 1;
	for (int r = 1; r < R; r++) {
		oStrides[r] = D * oStrides[r - 1];
	}
	iStrides = oStrides;
	auto const iTmp = iStrides;
	auto const oTmp = oStrides;
	for (int r = 0; r < R; r++) {
		iStrides[r] = iTmp[inputIndices[r]];
		oStrides[r] = oTmp[outputIndices[r]];
		printf("%i %i %i\n", r, oStrides[r], iStrides[r]);
	}
	printf( "\n");
	std::vector<int> I(R, 0);
	int iIndex = 0;
	int oIndex = 0;
	std::vector<int> ioMapping(size);
	for (int count = 0; count < size; count++) {
		printf("%i %i %i\n", count, iIndex, oIndex);
		ioMapping[iIndex] = oIndex;
		int dim = 0;
		while (++I[dim] == D) {
			iIndex -= (D - 1) * iStrides[dim];
			oIndex -= (D - 1) * oStrides[dim];
			I[dim++] = 0;
		}
		iIndex += iStrides[dim];
		oIndex += oStrides[dim];
	}
}
;

void testEinstein() {
	TensorMapping<NDIM>(std::vector<int> { 2, 1, 0 }, std::vector<int> { 1, 0, 2 });
//	Z4DynamicVariables<double> st;
//	const char *name[40] = { "Kxx", "Kxy", "Kyy", "Kxz", "Kxy", "Kzz", "Ax", "Ay", "Az", "Bxx", "Bxy", "Bxz", "Byx", "Byy", "Byz", "Bzx", "Bzy", "Bzz", "Dxxx",
//			"Dxxy", "Dxyy", "Dxxz", "Dxyz", "Dxzz", "Dyxx", "Dyxy", "Dyyy", "Dyxz", "Dyyz", "Dyzz", "Dzxx", "Dzxy", "Dzyy", "Dzxz", "Dzyz", "Dzzz", "The0",
//			"Zx", "Zy", "Zz" };
}

