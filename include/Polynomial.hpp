/******************************************************************************
 Copyright (C) 2024  Dominic C. Marcello
 *******************************************************************************/

#ifndef INCLUDE_POLYNOMIAL_HPP_
#define INCLUDE_POLYNOMIAL_HPP_

#include <algorithm>
#include <array>
#include <complex>
#include <iostream>
#include <ostream>
#include <unordered_map>
#include <vector>
#include "Util.hpp"

template<typename T, size_t dimensionCount>
struct ArrayHash {
	size_t operator()(std::array<T, dimensionCount> const &arr) const {
		size_t seed = dimensionCount;
		for (T const &val : arr) {
			seed ^= std::hash<T> { }(val) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
		}
		return seed;
	}
};

template<typename T>
struct Monomial {
	T coefficient;
	int degree;
	Monomial() = default;
	Monomial(T c, int d) :
			coefficient(c), degree(d) {
	}
	Monomial operator/(Monomial<T> const &divisor) const {
		Monomial quotient;
		quotient.degree = degree - divisor.degree;
		quotient.coefficient = coefficient / divisor.coefficient;
		return quotient;
	}
};

template<typename T, int dimensionCount>
struct MultivariatePolynomial;

template<typename T, int dimensionCount>
class CoefficientReference {
	MultivariatePolynomial<T, dimensionCount> *ptr;
	std::array<int, dimensionCount> n;
	void checkZero() {
		if (std::abs(ptr->coefficients_[n]) < std::abs(std::numeric_limits<T>::epsilon())) {
			ptr->coefficients_.erase(n);
		}
	}
	void initZero() {
		if (ptr->coefficients_.find(n) == ptr->coefficients_.end()) {
			ptr->coefficients_[n] = 0;
		}
	}
public:
	CoefficientReference(std::array<int, dimensionCount> _n, MultivariatePolynomial<T, dimensionCount> &ref) {
		n = _n;
		ptr = &ref;
	}
	operator T() const {
		if (ptr->coefficients_.find(n) == ptr->coefficients_.end()) {
			return 0;
		}
		return ptr->coefficients_[n];
	}
	CoefficientReference& operator=(T const &value) {
		ptr->coefficients_[n] = value;
		checkZero();
		return *this;
	}
	CoefficientReference& operator-=(T const &value) {
		initZero();
		ptr->coefficients_[n] -= value;
		checkZero();
		return *this;
	}
	CoefficientReference& operator+=(T const &value) {
		initZero();
		ptr->coefficients_[n] += value;
		checkZero();
		return *this;
	}
	CoefficientReference& operator*=(T const &value) {
		initZero();
		ptr->coefficients_[n] *= value;
		checkZero();
		return *this;
	}
	CoefficientReference& operator/=(T const &value) {
		initZero();
		ptr->coefficients_[n] /= value;
		checkZero();
		return *this;
	}
};

template<typename T, int dimensionCount>
struct MultivariatePolynomial {
	MultivariatePolynomial() = default;
	MultivariatePolynomial(Monomial<T> const &term) {
		(*this)[term.degree] = term.coefficient;
	}
	CoefficientReference<T, dimensionCount> operator[](std::array<int, dimensionCount> indices) {
		auto it = coefficients_.find(indices);
		if (it == coefficients_.end()) {
			it = coefficients_.insert(std::pair(indices, T(0))).first;
		}
		return CoefficientReference<T, dimensionCount>(it->first, *this);
	}
	T const operator[](std::array<int, dimensionCount> indices) const {
		auto it = coefficients_.find(indices);
		if (it == coefficients_.end()) {
			return T(0);
		}
		if (std::abs(it->second) < std::abs(std::numeric_limits<T>::epsilon())) {
			return T(0);
		}
		return it->second;
	}
	template<typename T2>
	T2 const operator()(std::array<T2, dimensionCount> const &x) const {
		T2 sum = T2(0);
		for (auto it = coefficients_.begin(); it != coefficients_.end(); it++) {
			T2 product = T2(it->second);
			for (unsigned r = 0; r < it->first.size(); r++) {
				int const exponent = it->first[r];
				product *= ipow(x[r], exponent);
			}
			sum += product;
		}
		return sum;
	}
	template<int one = dimensionCount, typename = std::enable_if_t<one == 1>>
	auto operator[](int i) {
		std::array<int, dimensionCount> indices = { i };
		return (*this)[indices];
	}
	template<int one = dimensionCount, typename = std::enable_if_t<one == 1>>
	T const operator[](int i) const {
		std::array<int, dimensionCount> indices = { i };
		return (*this)[indices];
	}
	template<typename T2, int one = dimensionCount, typename = std::enable_if_t<one == 1>>
	T2 const operator()(T2 x) const {
		std::array<T2, dimensionCount> xs = { x };
		return (*this)(xs);
	}
	friend auto operator+(MultivariatePolynomial const &a) {
		return a;
	}
	friend auto operator-(MultivariatePolynomial const &b) {
		MultivariatePolynomial a;
		for (auto &ele : b.coefficients_) {
			a[ele.first] = -ele.second;
		}
		return a;
	}
	friend auto operator+(MultivariatePolynomial const &b, MultivariatePolynomial const &c) {
		MultivariatePolynomial a = b;
		for (auto &ele : c.coefficients_) {
			a[ele.first] += ele.second;
		}
		return a;
	}
	friend auto operator-(MultivariatePolynomial const &b, MultivariatePolynomial const &c) {
		MultivariatePolynomial a = b;
		for (auto &ele : c.coefficients_) {
			a[ele.first] -= ele.second;
		}
		return a;
	}
	friend auto operator*(MultivariatePolynomial const &b, MultivariatePolynomial const &c) {
		MultivariatePolynomial a;
		for (auto &eleB : b.coefficients_) {
			for (auto &eleC : c.coefficients_) {
				std::array<int, dimensionCount> indices = eleB.first;
				for (int r = 0; r < dimensionCount; ++r) {
					indices[r] += eleC.first[r];
				}
				a[indices] += eleB.second * eleC.second;
			}
		}
		return a;
	}
	friend auto operator*(MultivariatePolynomial const &b, T const &c) {
		MultivariatePolynomial a = b;
		for (auto &eleA : a.coefficients_) {
			eleA.second *= c;
		}
		return a;
	}
	friend auto operator*(T const &b, MultivariatePolynomial const &c) {
		return c * b;
	}
	MultivariatePolynomial& operator+=(MultivariatePolynomial const &other) {
		*this = *this + other;
		return *this;
	}
	MultivariatePolynomial& operator-=(MultivariatePolynomial const &other) {
		*this = *this - other;
		return *this;
	}
	MultivariatePolynomial& operator*=(MultivariatePolynomial const &other) {
		*this = *this * other;
		return *this;
	}
	MultivariatePolynomial& operator*=(T const &other) {
		*this = *this * other;
		return *this;
	}
	friend std::ostream& operator<<(std::ostream &os, MultivariatePolynomial const &p) {
		bool firstTerm = true;
		for (auto const& [indices, coeff] : p.coefficients_) {
			if (std::abs(coeff) < 1e-12) {
				continue;
			}
			if (!firstTerm) {
				os << " + ";
			}
			bool needCoeff = std::abs(coeff) != T(1) || std::all_of(indices.begin(), indices.end(), [](int x) {
				return x == 0;
			});
			bool printedSomething = false;
			if (needCoeff) {
				os << coeff;
				printedSomething = true;
			}
			for (int j = 0; j < dimensionCount; ++j) {
				if (indices[j] != 0) {
					if (printedSomething) {
						os << "*";
					}
					if (dimensionCount > 1) {
						os << "x" << j;
					} else {
						os << "x";
					}
					if (indices[j] != 1) {
						os << "^" << indices[j];
					}
					printedSomething = true;
				}
			}
			firstTerm = false;
		}
		if (firstTerm) {
			os << "0";
		}
		return os;
	}
	template<typename V>
	MultivariatePolynomial(MultivariatePolynomial<V, dimensionCount> const &other) {
		for (auto it = other.coefficients_.begin(); it != other.coefficients_.end(); it++) {
			coefficients_[it->first] = it->second;
		}
	}
	template<int one = dimensionCount, typename = std::enable_if_t<one == 1>>
	Monomial<T> leadingTerm() const {
		Monomial<T> term;
		auto &maxDegree = term.degree;
		auto &coefficient = term.coefficient;
		maxDegree = -1;
		coefficient = 0;
		for (auto it = coefficients_.begin(); it != coefficients_.end(); it++) {
			if (std::abs(it->second) && maxDegree < it->first[0]) {
				maxDegree = it->first[0];
				coefficient = it->second;
			}
		}
		return term;
	}
	MultivariatePolynomial truncate(int maxDegree) const {
		MultivariatePolynomial P;
		for (auto it = coefficients_.begin(); it != coefficients_.end(); it++) {
			auto const index = it->first;
			int const deg = std::accumulate(index.begin(), index.end(), 0);
			if (deg <= maxDegree) {
				P[index] = it->second;
			}
		}
		return P;
	}
	template<int one = dimensionCount, typename = std::enable_if_t<one == 1>>
	MultivariatePolynomial const operator()(MultivariatePolynomial const &q) const {
		MultivariatePolynomial const &p = *this;
		MultivariatePolynomial qn, r;
		qn[0] = T(1);
		int const N = degree();
		for (int n = 0; n <= N; n++) {
			r += p[n] * qn;
			qn *= q;
		}
		return r;
	}
	template<int one = dimensionCount, typename = std::enable_if_t<one == 1>>
	T leadingCoefficient() const {
		int maxDegree = -1;
		T coefficient = 0;
		for (auto it = coefficients_.begin(); it != coefficients_.end(); it++) {
			if (std::abs(it->second) && maxDegree < it->first[0]) {
				maxDegree = it->first[0];
				coefficient = it->second;
			}
		}
		return coefficient;
	}
	template<int one = dimensionCount, typename = std::enable_if_t<one == 1>>
	int degree() const {
		int maxDegree = 0;
		for (auto it = coefficients_.begin(); it != coefficients_.end(); it++) {
			if (std::abs(it->second)) {
				maxDegree = std::max(maxDegree, it->first[0]);
			}
		}
		return maxDegree;
	}
	template<int one = dimensionCount, typename = std::enable_if_t<one == 1>>
	MultivariatePolynomial& operator+=(Monomial<T> const &b) {
		(*this)[b.degree] += b.coefficient;
		return *this;
	}
	template<int one = dimensionCount, typename = std::enable_if_t<one == 1>>
	MultivariatePolynomial operator+(Monomial<T> const &b) {
		auto a = *this;
		a += b;
		return a;
	}
	template<int one = dimensionCount, typename = std::enable_if_t<one == 1>>
	MultivariatePolynomial& operator-=(Monomial<T> const &b) {
		(*this)[b.degree] -= b.coefficient;
		return *this;
	}
	template<int one = dimensionCount, typename = std::enable_if_t<one == 1>>
	MultivariatePolynomial operator-(Monomial<T> const &b) {
		auto a = *this;
		a -= b;
		return a;
	}
	template<int one = dimensionCount, typename = std::enable_if_t<one == 1>>
	MultivariatePolynomial operator*(Monomial<T> const &b) const {
		MultivariatePolynomial c;
		for (auto it = coefficients_.begin(); it != coefficients_.end(); it++) {
			c[it->first[0] + b.degree] = it->second * b.coefficient;
		}
		return c;
	}
	template<int one = dimensionCount, typename = std::enable_if_t<one == 1>>
	friend MultivariatePolynomial operator*(Monomial<T> const &a, MultivariatePolynomial const &b) {
		return b * a;
	}
	template<int one = dimensionCount, typename = std::enable_if_t<one == 1>>
	MultivariatePolynomial& operator*=(Monomial<T> const &a) {
		*this = *this * a;
		return a;
	}
	friend MultivariatePolynomial polynomialAntiderivative(MultivariatePolynomial const &p, int d = 0) {
		MultivariatePolynomial P;
		for (auto it = p.coefficients_.begin(); it != p.coefficients_.end(); it++) {
			int const power = it->first[d];
			P[power + 1] = it->second / T(1 + power);
		}
		return P;
	}
	friend MultivariatePolynomial polynomialDerivative(MultivariatePolynomial const &p, int d = 0) {
		MultivariatePolynomial dpdx;
		for (auto it = p.coefficients_.begin(); it != p.coefficients_.end(); it++) {
			int const power = it->first[d];
			if (power > 0) {
				dpdx[power - 1] = T(power) * it->second;
			}
		}
		return dpdx;
	}
	friend class CoefficientReference<T, dimensionCount> ;
	template<typename, int>
	friend class MultivariatePolynomial;
private:
	std::unordered_map<std::array<int, dimensionCount>, T, ArrayHash<int, dimensionCount>> coefficients_;
};

template<typename T>
using Polynomial = MultivariatePolynomial<T, 1>;


template<typename T>
Polynomial<T> operator+(Monomial<T> const &a, Monomial<T> const &b) {
	Polynomial<T> c;
	c += a;
	c += b;
	return c;
}

template<typename T>
Polynomial<T> operator-(Monomial<T> const &a, Monomial<T> const &b) {
	Polynomial<T> c;
	c += a;
	c -= b;
	return c;
}

template<typename T>
struct PolynomialDivision {
	Polynomial<T> quotient;
	Polynomial<T> remainder;
};

template<typename T>
PolynomialDivision<T> divide(Polynomial<T> const &dividend, Polynomial<T> const &divisor) {
	PolynomialDivision<T> rc;
	auto &quotient = rc.quotient;
	auto &remainder = rc.remainder;
	remainder = dividend;
	auto remainderDegree = remainder.degree();
	auto divisorDegree = divisor.degree();
	while (remainderDegree >= divisorDegree) {
		auto const quotientLeadingTerm = remainder.leadingTerm() / divisor.leadingTerm();
		quotient += quotientLeadingTerm;
		remainder -= quotientLeadingTerm * divisor;
		remainderDegree = remainder.degree();
		divisorDegree = divisor.degree();
	}
	return rc;
}

template<typename T>
Polynomial<T> operator/(Polynomial<T> const &dividend, Polynomial<T> const &divisor) {
	return divide(dividend, divisor).quotient;
}

template<typename T>
Polynomial<T>& operator/=(Polynomial<T> &a, Polynomial<T> const &b) {
	a = a / b;
	return a;
}

template<typename T>
Polynomial<T>& operator/=(Polynomial<T> &a, Monomial<T> const &b) {
	a /= Polynomial<T>(b);
	return a;
}

template<typename T>
Polynomial<T> operator%(Polynomial<T> const &dividend, Polynomial<T> const &divisor) {
	return divide(dividend, divisor).remainder;
}

template<typename T>
std::vector<std::complex<T>> polynomialRoots(Polynomial<T> pReal) {
	using MonomialType = Monomial<std::complex<T>>;
	Polynomial<std::complex<T>> p(pReal);
	std::vector<std::complex<T>> roots;
	std::complex<T> z(1, 1);
	while (p.degree()) {
		T err;
		for (int iter = 0; iter < 20; iter++) {
			auto const pz = p(z);
			auto const dpdz = derivative(p)(z);
			auto const dz = -pz / dpdz;
			z += dz;
			err = std::abs(dz / z);
			printf("%i %e\n", iter, err);
		}
		std::cout << z << "\n";
		roots.push_back(z);
		auto const pd = MonomialType(1, 1) - MonomialType(z, 0);
		p /= pd;
	}
	return roots;
}

#endif /* INCLUDE_POLYNOMIAL_HPP_ */
