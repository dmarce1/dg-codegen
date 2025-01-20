/*
 * LegendreP.hpp
 *
 *  Created on: Dec 15, 2024
 *      Author: dmarce1
 */

#ifndef INCLUDE_LEGENDREP_HPP_
#define INCLUDE_LEGENDREP_HPP_

#include "TriangularArray.hpp"

namespace Math {

template<typename T, int P>
Vector<T, P> legendreBasis(T x) {
	static constexpr T one(1);
	Vector<T, P> Pn;
	Pn[0] = one;
	if (P > 1) {
		Pn[1] = x;
	}
	for (int p = 0; p < P - 1; p++) {
		Pn[p + 1] = (T(2 * p + 1) * x * Pn[p] - T(p) * Pn[p - 1]) * invInteger<T>(p + 1);
	}
	return Pn;
}

template<typename T, int D, int P>
TriangularArray<T, D, P> legendreBasis(Vector<T, D> const &x) {
	static constexpr T one(1);
	TriangularArray<T, D, P> Qn;
	Vector<Vector<T, P>, D> Pn;
	for (int k = 0; k < D; k++) {
		Pn[k] = legendreBasis(x[k]);
	}
	for (TriangularIndices<D, P> I; I != I.end(); I++) {
		Qn[I] = one;
		for (int k = 0; k < D; k++) {
			Qn[I] *= Pn[k];
		}
	}
	return Qn;
}

template<typename Type>
Type legendreP(int n, Type x) {
	static constexpr Type one(1);
	if (n == 0) {
		return one;
	} else if (n == 1) {
		return x;
	} else {
		Type Pnm1, Pn, Pnp1;
		Pnm1 = one;
		Pn = x;
		for (int l = 1; l < n; l++) {
			Pnp1 = (Type(2 * l + 1) * Pn * x - Type(l) * Pnm1) * invInteger<Type>(l + 1);
			Pnm1 = Pn;
			Pn = Pnp1;
		}
		return Pn;
	}
}

template<typename Type>
Type dLegendrePdX(int n, Type x) {
	static constexpr Type zero(0), one(1);
	if (n == 0) {
		return zero;
	} else if (n == 1) {
		return one;
	} else {
		Type Pnm1, Pn, Pnp1;
		Type dPnm1dX, dPndX, dPnp1dX;
		Pnm1 = one;
		Pn = x;
		dPndX = one;
		for (int l = 1; l < n - 1; l++) {
			Pnp1 = (Type(2 * l + 1) * Pn * x - Type(l) * Pnm1) * invInteger<Type>(l + 1);
			dPnp1dX = Type(l + 1) * Pn + x * dPndX;
			Pnm1 = Pn;
			Pn = Pnp1;
			dPndX = dPnp1dX;
		}
		return Type(n) * Pn + x * dPndX;
	}
}

template<typename Type>
Type d2LegendrePdX2(int n, Type x) {
	static constexpr Type zero(0), one(1), three(3);
	if (n <= 1) {
		return zero;
	} else if (n == 2) {
		return three;
	} else {
		Type Pnm1, Pn, Pnp1;
		Type dPnm1dX, dPndX, dPnp1dX;
		Type d2Pnm1dX2, d2PndX2, d2Pnp1dX2;
		Pnm1 = one;
		Pn = x;
		dPndX = one;
		d2PndX2 = three;
		for (int l = 1; l < n - 2; l++) {
			Pnp1 = (Type(2 * l + 1) * Pn * x - Type(l) * Pnm1) * invInteger<Type>(l + 1);
			dPnp1dX = Type(l + 1) * Pn + x * dPndX;
			d2Pnp1dX2 = Type(l + 2) * dPndX + x * d2PndX2;
			Pnm1 = Pn;
			Pn = Pnp1;
			dPndX = dPnp1dX;
			d2PndX2 = d2Pnp1dX2;
		}
		dPnp1dX = Type(n - 1) * Pn + x * dPndX;
		d2PndX2 = Type(n) * dPndX + x * d2PndX2;
		dPndX = dPnp1dX;
		return Type(n + 1) * dPndX + x * d2PndX2;
	}
}

}

#endif /* INCLUDE_LEGENDREP_HPP_ */
