/*
 * Einstein.cpp
 *
 *  Created on: Feb 25, 2025
 *      Author: dmarce1
 */

#include <numbers>

#include "Tensor.hpp"

using namespace Tensors;

template<typename T>
struct Z4SpaceTime {
	using tensor0 = T;
	using tensor1 = Tensor<T, NDIM, 1>;
	using tensor2 = Tensor<T, NDIM, 2>;
	using tensorS = Tensor<T, NDIM, 2, Sym<0, 1>>;
	using tensor3 = Tensor<T, NDIM, 3, Sym<0, 1>>;
private:
	static constexpr int zeta = 0;
	static constexpr int zeta1 = 0;
	static constexpr int chi = 0;
	static constexpr Index<'i'> i { };
	static constexpr Index<'j'> j { };
	static constexpr Index<'k'> k { };
	static constexpr Index<'m'> m { };
	static constexpr Index<'n'> n { };
	static constexpr Index<'o'> o { };
	static constexpr Index<'p'> p { };
	static constexpr Index<'q'> q { };
	static constexpr Delta<T, NDIM> delta { };
	static constexpr T Pi = std::numbers::pi_v<T>;
	tensor0 alpha;
	tensor0 Theta;
	tensor1 beta;
	tensor1 A;
	tensor1 Z;
	tensor2 B;
	tensorS gamma;
	tensorS K;
	tensor3 D;
public:
	Z4SpaceTime flux(int l) {
		tensor0 trK, Q0;
		tensor1 trD, E, V, Q1;
		tensorS gamma_inv, Q2, lambda;
		Z4SpaceTime F;
		gamma_inv = gamma.inverse();
		trD(i) = D(i, j, m) * gamma_inv(j, m);
		E(i) = D(j, i, m) * gamma_inv(j, m);
		V(i) = trD(i) - E(i) - Z(i);
		trK = K(i, j) * gamma_inv(i, j);
		Q0 = trK - 2 * Theta;
		Q1(i) = A(i) + trD(i) - 2 * E(i) - 2 * Z(i);
		Q2(i, j) = K(i, j) - (0.5 / alpha) * (B(i, m) * gamma(m, j) + B(j, m) * gamma(m, i));
		lambda(i, j) = gamma_inv(l, m) * D(m, i, j) + 0.5 * (delta(l, i) * Q1(j) + delta(l, j) * Q1(i));
		if constexpr (zeta != -1) {
			lambda(i, j) += -0.5 * (zeta + 1) * ((D(i, j, m) + D(j, i, m)) * gamma_inv(m, l) + delta(l, i) * E(j) + delta(l, j) * E(i));
		}
		F.alpha = 0;
		F.beta(i) = 0;
		F.gamma(i, j) = 0;
		F.A(i) = -beta(l) * A(i) + alpha * delta(l, i) * Q0;
		F.B(i, j) = -beta(l) * B(i, j) + alpha * delta(l, i) * Q1(m) * gamma_inv(m, j);
		F.D(k, i, j) = -beta(l) * D(k, i, j) + alpha * delta(l, k) * Q2(i, j);
		F.K(i, j) = -beta(l) * K(i, j) + alpha * lambda(i, j);
		F.Theta = -beta(l) * Theta + alpha * V(m) * gamma_inv(m, l);
		F.Z(i) = -beta(l) * Z(i) + alpha * (delta(l, i) * (trK - Theta) - K(i, m) * gamma_inv(m, l));
		if constexpr (zeta1 != 0) {
			F.Z(i) = F.Z(i) + zeta1 * (B(i, l) - delta(i, l) * delta(m, k) * B(m, k));
		}
		return F;
	}
	Z4SpaceTime vacuum_source() {
		tensor0 trK, Q0, trB;
		tensor1 Q1, trD, E;
		tensorS gamma_inv, Q2;
		tensor3 Gamma;
		Z4SpaceTime src;
		gamma_inv = gamma.inverse();
		trB = delta(m, k) * B(m, k);
		trK = K(m, k) * gamma_inv(m, k);
		trD(i) = D(i, k, m) * gamma_inv(k, m);
		E(i) = D(k, i, m) * gamma_inv(k, m);
		Q0 = trK - 2 * Theta;
		Q1(i) = A(i) + trD(i) - 2 * E(i) - 2 * Z(i);
		Q2(i, j) = K(i, j) - (0.5 / alpha) * (B(i, m) * gamma(m, j) + B(j, m) * gamma(m, i));
		Gamma(k, i, j) = 0.5 * gamma_inv(k, m) * (D(i, j, m) + D(j, i, m) - D(m, j, i));
		src.alpha = alpha * (beta(m) * A(m) - alpha * Q0);
		src.beta(i) = beta(m) * B(m, i) - alpha * Q1(m) * gamma_inv(m, i);
		src.gamma(i, j) = 2 * (beta(m) * D(m, i, j) - alpha * Q2(i, j));
		src.A(i) = B(i, m) * A(m) - trB * A(i);
		src.B(i, j) = B(i, m) * B(m, j) - trB * B(i, j);
		src.D(k, i, j) = B(k, m) * D(m, i, j) - trB * D(k, i, j);
		src.K(i, j) = alpha * ((trD(k) + A(k) - 2 * Z(k)) * Gamma(k, i, j) - Gamma(k, m, j) * Gamma(m, k, i) - A(i) * Z(j) + A(j) * Z(i));
		src.K(i, j) -= trB * K(i, j) - K(i, k) * B(j, k) - K(j, k) * B(i, k) + alpha * (2 * gamma_inv(k, m) * K(m, i) * K(k, j) - (trK - 2 * Theta) * K(i, j));
		if constexpr (chi != -1) {
			src.K(i, j) += alpha * (0.5 * (1 + chi) * (-A(k) * Gamma(k, i, j) + 0.5 * (A(i) * trD(j) + A(j) * trD(i))));
		}
		if constexpr (chi != +1) {
			T const co = alpha * 0.5 * (1 - chi);
			src.K(i, j) += co * (A(k) * gamma_inv(k, m) * D(m, i, j) - 0.5 * (A(i) * (2 * E(j) - trD(j)) + A(j) * (2 * E(i) - trD(i))));
			src.K(i, j) += 2 * co * (D(i, p, k) * gamma_inv(k, m) * gamma_inv(p, q) * D(q, m, j) + D(j, p, k) * gamma_inv(k, m) * gamma_inv(p, q) * D(q, m, i));
			src.K(i, j) -= 2 * co * E(k) * (D(i, j, m) * gamma_inv(m, k) + D(j, i, m) * gamma_inv(m, k));
		}
		src.Z(i) = -Z(i) * trB + Z(k) * B(i, k) + alpha * (A(i) * (trK - 2 * Theta) - A(k) * gamma_inv(k, m) * K(m, i));
		src.Z(i) += alpha * (gamma_inv(k, m) * K(m, p) * Gamma(p, k, i) + gamma_inv(k, m) * K(m, i) * (trD(k) - 2 * Z(k)));
		src.Theta = -Theta * trB + 0.5 * alpha * (2 * A(k) * gamma_inv(k, m) * (trD(m) - E(m) - 2 * Z(m)));
		src.Theta += 0.5 * alpha * (D(i, p, q) * gamma_inv(p, j) * gamma_inv(q, k) * Gamma(i, j, k) - gamma_inv(k, m) * trD(m) * (trD(k) - 2 * Z(k)));
		src.Theta -= 0.5 * alpha * (gamma_inv(p, m) * K(m, q) * gamma_inv(q, k) * K(k, p) - trK * (trK - 2 * Theta));
		return src;
	}
	Z4SpaceTime matter_source(tensor0 const &tau, tensor1 const &s, tensorS const &S) {
		tensorS gamma_inv;
		Z4SpaceTime src;
		gamma_inv = gamma.inverse();
		src.alpha = 0;
		src.beta(i) = 0;
		src.gamma(i, j) = 0;
		src.A(i) = 0;
		src.B(i, j) = 0;
		src.D(k, i, j) = 0;
		src.Theta = -(8.0 * Pi) * tau;
		src.Z(i) = -(8.0 * Pi) * s(i);
		src.K(i, j) = -(8.0 * Pi) * (S(i, j) - 0.5 * (S(k, m) * gamma_inv(k, m) - tau) * gamma(i, j));
		return src;
	}

};

void testEinstein() {
	Z4SpaceTime<double> st;
	st.flux(0);
	st.vacuum_source();
	st.matter_source(Z4SpaceTime<double>::tensor0(), Z4SpaceTime<double>::tensor1(), Z4SpaceTime<double>::tensorS());
}

