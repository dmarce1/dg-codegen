/******************************************************************************
 Copyright (C) 2024  Dominic C. Marcello
 *******************************************************************************/

#include "Indent.hpp"
#include "Symbolic.hpp"
#include "TriIndex.hpp"

static Indent indent { };

constexpr int eulerFieldCount(int dimCount) {
	return dimCount + 2;
}

template<int dimCount, int maxOrder>
auto eulerFlux(int dim) {
	using namespace SymEngine;
	constexpr int fieldCount = dimCount + 2;
	SymbolicExpression const half = rational(1, 2);
	const int sx_i = 0;
	const int rho_i = dimCount;
	const int eg_i = dimCount + 1;
	using IndexType = TriIndex<maxOrder, fieldCount>;
	vec_basic Fin(fieldCount * IndexType::count());
	std::vector<SymbolicSymbol> U(fieldCount);
	for (int d = 0; d < dimCount; d++) {
		auto const name = std::string("S") + std::string(1, 'x' + d);
		U[sx_i + d] = symbol(name.c_str());
	}
	U[rho_i] = symbol("rho");
	U[eg_i] = symbol("eg");
	indent++;
	for(IndexType Dn = IndexType::begin(); Dn != IndexType::end(); Dn++) {
		vec_basic S(dimCount);
		vec_basic v(dimCount);
		auto* F = Fin.data() + fieldCount * int(Dn);
		for (int d = 0; d < dimCount; d++) {
			v[d] = U[sx_i + d] / U[rho_i];
		}
		for (int f = 0; f < fieldCount; f++) {
			F[f] = U[f] * v[dim];
		}
		SymbolicExpression ek = U[sx_i] * v[0];
		for (int d = 1; d < dimCount; d++) {
			ek += U[sx_i + d] * v[d];
		}
		ek *= rational(1, 2);
		SymbolicExpression const p = rational(2, 3) * (U[eg_i] - ek);
		F[sx_i + dim] += p;
		F[eg_i] += p * v[dim];

		for (int f = 0; f < fieldCount; f++) {
			int const n = Dn[f];
			for (int i = 0; i < n; i++) {
				for (int f2 = 0; f2 < fieldCount; f2++) {
					F[f2] = derivative(F[f2], U[f]);
				}
			}
		}
	}
	vec_pair subexprs;
	vec_basic Fout;
	for (int f = 0; f < fieldCount * IndexType::count(); f++) {
		Fin[f] = expand(Fin[f]);
	}
	cse(subexprs, Fout, Fin);
	factorPowersFromSubstitutions(subexprs);
	for (int f = 0; f < int(subexprs.size()); f++) {
		std::cout << indent << "T const " << *subexprs[f].first << " = " << toCxxCode(subexprs[f].second) << ";\n";
	}
	for(IndexType Dn = IndexType::begin(); Dn != IndexType::end(); Dn++) {
		int const dFlat = int(Dn);
		auto* F = Fout.data() + fieldCount * dFlat;
		for (int f = 0; f < fieldCount; f++) {
			if(toCxxCode(F[f]) != "0") {
				std::cout << indent << "F[" << std::to_string(f) << "][" << std::to_string(Dn) << "] " << toCxxCode(F[f]) << ";\n";
			}
		}
	}
	indent--;
}

constexpr int maxOrder = 4;

void stateCodegen() {
	eulerFlux<3, maxOrder>(0);
//	constexpr int dimCount = 1;
//	constexpr int fieldCount = eulerFieldCount(dimCount);
//	using IndexType = TriIndex<maxOrder, fieldCount>;
//	for(IndexType index = IndexType::begin(); index < IndexType::end(); index++) {
//		std::vector<int> indices(fieldCount);
//		for(int i = 0; i < fieldCount; i ++) {
//			indices[i] = index[i];
//		}
//		printf( "\n" );
//		eulerFlux(dimCount, 0, indices);
//	}
}
