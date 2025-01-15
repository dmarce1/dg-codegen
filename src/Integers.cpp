/*
 * Integers.cpp
 *
 *  Created on: Jan 14, 2025
 *      Author: dmarce1
 */

#include <cassert>
#include <cmath>
#include <cstdio>
#include <set>
#include <vector>

using Integer = long long;

bool isPrime(Integer testNumber) {
	static constexpr Integer zero = Integer(0);
	static thread_local Integer maxTested;
	static thread_local std::set<Integer> primes;
	if (primes.size() == 0) {
		primes.insert(2);
		maxTested = 2;
	}
	while (testNumber > maxTested) {
		Integer const candidatePrime = maxTested + 1;
		Integer const maxFactor = sqrt(candidatePrime);
		bool candidatePasses = true;
		for (auto i = primes.begin(); *i <= maxFactor; ++i) {
			if (candidatePrime % *i == zero) {
				candidatePasses = false;
				break;
			}
		}
		if (candidatePasses) {
			primes.insert(candidatePrime);
		}
		maxTested++;
	}
	return bool(primes.find(testNumber) != primes.end());
}

Integer factorial(Integer n) {
	assert(n <= 20);
	static constexpr Integer one = Integer(1);
	static thread_local std::vector<Integer> fact;
	while (n >= Integer(fact.size())) {
		if (fact.size() > 0) {
			fact.push_back(fact.back() * fact.size());
		} else {
			fact.push_back(one);
		}
	}
	return fact[n];
}

Integer doubleFactorial(Integer n) {
	assert(n <= 33);
	static constexpr Integer one = Integer(1);
	static thread_local std::vector<Integer> dFact;
	while (n >= Integer(dFact.size())) {
		if (dFact.size() > 1) {
			Integer const m = dFact.size();
			dFact.push_back(m * dFact[m - 2]);
		} else {
			dFact.push_back(one);
		}
	}
	return dFact[n];
}

Integer binomialCoefficient(Integer j, Integer k) {
	assert(j <= 66);
	static constexpr Integer one = Integer(1);
	static thread_local std::vector<std::vector<Integer>> binCo;
	while (j >= Integer(binCo.size())) {
		if (binCo.size() > 0) {
			Integer const n = binCo.size();
			binCo.resize(n + 1);
			binCo[n].resize(n + 1);
			binCo[n][0] = one;
			for (int k = 1; k < n; k++) {
				binCo[n][k] = binCo[n - 1][k - 1] + binCo[n - 1][k];
			}
			binCo[n][n] = one;
		} else {
			binCo.resize(1);
			binCo[0].resize(1);
			binCo[0][0] = 1;
		}
	}
	return binCo[j][k];
}

void testIntegers() {
	constexpr Integer N = 100;
	int cnt = 0;
	isPrime(1LL << 26LL);
}
