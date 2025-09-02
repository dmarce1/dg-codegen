/******************************************************************************
 Copyright (C) 2024  Dominic C. Marcello
 *******************************************************************************/

#ifndef INCLUDE_BICONJUGATEGRADIENT_HPP_
#define INCLUDE_BICONJUGATEGRADIENT_HPP_

#include "Matrix.hpp"
//r0 = b − Ax0
//Choose an arbitrary vector r̂0 such that (r̂0, r0) ≠ 0, e.g., r̂0 = r0
//ρ0 = (r̂0, r0)
//p0 = r0
//For i = 1, 2, 3, …
//v = Api−1
//α = ρi−1/(r̂0, v)
//h = xi−1 + αpi−1
//s = ri−1 − αv
//If h is accurate enough, i.e., if s is small enough, then set xi = h and quit
//t = As
//ω = (t, s)/(t, t)
//xi = h + ωs
//ri = s − ωt
//If xi is accurate enough, i.e., if ri is small enough, then quit
//ρi = (r̂0, ri)
//β = (ρi/ρi−1)(α/ω)
//pi = ri + β(pi−1 − ωv)

template<typename T, int N>
std::array<T, N> biconjugateGradientSolve(SquareMatrix<T, N> const &A, std::array<T, N> const &b) {
	constexpr T toler = 1e-5;
	std::array<T, N> pn, pnp1, rn, rnp1, xn, xnp1;
	T rhon, rhonp1;
	xnp1.fill(0);
	rnp1 = b;
	rhonp1 = dot(rnp1, rnp1);
	pnp1 = rnp1;
	auto const r0 = rnp1;
	do {
		pn = pnp1;
		rn = rnp1;
		xn = xnp1;
		rhon = rhonp1;
		auto const v = A * pn;
		auto const alpha = rhon / dot(r0, v);
		auto const h = xn + alpha * pn;
		auto const s = rn - alpha * v;
		auto const t = A * s;
		auto const omega = dot(t, s) / dot(t, t);
		xnp1 = h + omega * s;
		rnp1 = s - omega * t;
		rhonp1 = dot(r0, rnp1);
		auto const beta = (rhonp1 / rhon) * (alpha / omega);
		pnp1 = rnp1 + beta * (pn - omega * v);
		for (int d = 0; d < N; d++) {
			printf("%e ", xnp1[d]);
		}
		printf(" - %e\n", norm(xnp1 - xn));
	} while (norm(xnp1 - xn) > T(0));
	return xn;
}

#endif /* INCLUDE_BICONJUGATEGRADIENT_HPP_ */
