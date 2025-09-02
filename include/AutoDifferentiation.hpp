#pragma once
#pragma once

#include <array>
#include <cmath>
#include <iostream>

#include "ContainerArithmetic.hpp"
#include "Polynomial.hpp"

template<typename T, size_t N>
class Auto {
	T value_;
	std::array<T, N> gradient_;

public:

	constexpr Auto() :
			value_(T(0)) {
		gradient_.fill(T(0));
	}

	constexpr Auto(T val, size_t varIndex) :
			value_(val) {
		gradient_.fill(T(0));
		gradient_[varIndex] = T { 1 };
	}

	constexpr operator T() const {
		return value_;
	}

	constexpr auto gradient() const {
		return gradient_;
	}

	constexpr T gradient(int i) const {
		return gradient_[i];
	}

	constexpr Auto operator+() const {
		return *this;
	}

	constexpr Auto operator-() const {
		Auto res;
		res.value_ = -value_;
		for (size_t i = 0; i < N; ++i) {
			res.gradient_[i] = -gradient_[i];
		}
		return res;
	}

	constexpr Auto operator+(Auto const &other) const {
		Auto res;
		res.value_ = value_ + other.value_;
		for (size_t i = 0; i < N; ++i) {
			res.gradient_[i] = gradient_[i] + other.gradient_[i];
		}
		return res;
	}

	constexpr Auto operator-(Auto const &other) const {
		Auto res;
		res.value_ = value_ - other.value_;
		for (size_t i = 0; i < N; ++i) {
			res.gradient_[i] = gradient_[i] - other.gradient_[i];
		}
		return res;
	}

	constexpr Auto operator*(Auto const &other) const {
		Auto res;
		res.value_ = value_ * other.value_;
		for (size_t i = 0; i < N; ++i) {
			res.gradient_[i] = value_ * other.gradient_[i] + gradient_[i] * other.value_;
		}
		return res;
	}

	constexpr Auto operator/(Auto const &other) const {
		Auto res;
		res.value_ = value_ / other.value_;
		for (size_t i = 0; i < N; ++i) {
			res.gradient_[i] = (gradient_[i] * other.value_ - value_ * other.gradient_[i]) / (other.value_ * other.value_);
		}
		return res;
	}

	constexpr Auto& operator+=(Auto const &other) {
		value_ += other.value_;
		for (size_t i = 0; i < N; ++i) {
			gradient_[i] += other.gradient_[i];
		}
		return *this;
	}

	constexpr Auto& operator-=(Auto const &other) {
		value_ -= other.value_;
		for (size_t i = 0; i < N; ++i) {
			gradient_[i] -= other.gradient_[i];
		}
		return *this;
	}

	constexpr Auto& operator*=(Auto const &other) {
		T newVal = value_ * other.value_;
		for (size_t i = 0; i < N; ++i) {
			gradient_[i] = value_ * other.gradient_[i] + gradient_[i] * other.value_;
		}
		value_ = newVal;
		return *this;
	}

	constexpr Auto& operator/=(Auto const &other) {
		T denom = other.value_ * other.value_;
		T newVal = value_ / other.value_;
		for (size_t i = 0; i < N; ++i) {
			gradient_[i] = (gradient_[i] * other.value_ - value_ * other.gradient_[i]) / denom;
		}
		value_ = newVal;
		return *this;
	}

	constexpr Auto operator*(T scalar) const {
		Auto res;
		res.value_ = value_ * scalar;
		for (size_t i = 0; i < N; ++i) {
			res.gradient_[i] = gradient_[i] * scalar;
		}
		return res;
	}

	constexpr Auto operator/(T scalar) const {
		Auto res;
		res.value_ = value_ / scalar;
		for (size_t i = 0; i < N; ++i) {
			res.gradient_[i] = gradient_[i] / scalar;
		}
		return res;
	}

	constexpr Auto& operator*=(T scalar) {
		value_ *= scalar;
		for (size_t i = 0; i < N; ++i) {
			gradient_[i] *= scalar;
		}
		return *this;
	}

	constexpr Auto& operator/=(T scalar) {
		value_ /= scalar;
		for (size_t i = 0; i < N; ++i) {
			gradient_[i] /= scalar;
		}
		return *this;
	}

	friend constexpr Auto sin(Auto const &x) {
		Auto res;
		res.value_ = std::sin(x.value_);
		T deriv = std::cos(x.value_);
		for (size_t i = 0; i < N; ++i) {
			res.gradient_[i] = deriv * x.gradient_[i];
		}
		return res;
	}

	friend constexpr Auto cos(Auto const &x) {
		Auto res;
		res.value_ = std::cos(x.value_);
		T deriv = -std::sin(x.value_);
		for (size_t i = 0; i < N; ++i) {
			res.gradient_[i] = deriv * x.gradient_[i];
		}
		return res;
	}

	friend constexpr Auto sqrt(Auto const &x) {
		Auto res;
		res.value_ = std::sqrt(x.value_);
		T deriv = T { 0.5 } / res.value_;
		for (size_t i = 0; i < N; ++i) {
			res.gradient_[i] = deriv * x.gradient_[i];
		}
		return res;
	}

	friend constexpr Auto exp(Auto const &x) {
		Auto res;
		res.value_ = std::exp(x.value_);
		for (size_t i = 0; i < N; ++i) {
			res.gradient_[i] = res.value_ * x.gradient_[i];
		}
		return res;
	}

	friend constexpr Auto pow(Auto const &base, Auto const &exponent) {
		Auto res;
		res.value_ = std::pow(base.value_, exponent.value_);
		T logBase = std::log(base.value_);
		for (size_t i = 0; i < N; ++i) {
			res.gradient_[i] = res.value_ * (exponent.gradient_[i] * logBase + exponent.value_ * base.gradient_[i] / base.value_);
		}
		return res;
	}

	friend constexpr Auto pow(Auto const &base, T exponent) {
		Auto res;
		res.value_ = std::pow(base.value_, exponent);
		T coeff = exponent * std::pow(base.value_, exponent - T { 1 });
		for (size_t i = 0; i < N; ++i) {
			res.gradient_[i] = coeff * base.gradient_[i];
		}
		return res;
	}

	friend constexpr Auto pow(T base, Auto const &exponent) {
		Auto res;
		res.value_ = std::pow(base, exponent.value_);
		T logBase = std::log(base);
		for (size_t i = 0; i < N; ++i) {
			res.gradient_[i] = res.value_ * logBase * exponent.gradient_[i];
		}
		return res;
	}

	friend constexpr Auto operator*(T scalar, Auto const &x) {
		return x * scalar;
	}

	friend constexpr Auto operator/(T scalar, Auto const &x) {
		Auto res;
		T const inv = T(1) / x.value_;
		res.value_ = scalar * inv;
		T coeff = T(-1) / (inv * inv);
		for (size_t i = 0; i < N; ++i) {
			res.gradient_[i] = coeff * x.gradient_[i];
		}
		return res;
	}
};

#include "Matrix.hpp"
#include <functional>

template<typename Real, typename Auto, int rank>
std::array<Real, rank> newtonRhapson(std::function<std::array<Auto, rank>(std::array<Auto, rank> const&)> const &function,
		std::array<Real, rank> const &initialGuess, Real tolerance = 1e-12, int maxIter = 50) {
	std::array<Real, rank> x = initialGuess;
	std::array<Auto, rank> fx;
	SquareMatrix<Real, rank> J;
	for (int iter = 0; iter < maxIter; ++iter) {
		for (int i = 0; i < rank; ++i) {
			fx[i] = Auto(x[i], i);
		}
		auto const fval = function(fx);
		for (int n = 0; n < rank; ++n) {
			for (int m = 0; m < rank; ++m) {
				J[n][m] = fval[n].gradient(m);
			}
		}
		auto const delta = matrixInverse(J) * fval;

		Real normDeltaSq = 0;
		Real normXSq = 0;
		for (int i = 0; i < rank; ++i) {
			x[i] -= delta[i];
			normDeltaSq += delta[i] * delta[i];
			normXSq += x[i] * x[i];
		}

		Real const relativeChange = std::sqrt(normDeltaSq) / (std::sqrt(normXSq) + Real(1e-14));
		if (relativeChange < tolerance) {
			return x;
		}
	}
	THROW("max iters reached\n");
	return x;
}

class Partitions {
	void genPartitions(int n, int i, std::vector<std::vector<int>> &currentPartition) {
		if (i > n) {
			partitions_[n].push_back(currentPartition);
			return;
		}
		for (unsigned index = 0; index < currentPartition.size(); index++) {
			currentPartition[index].push_back(i);
			genPartitions(n, i + 1, currentPartition);
			currentPartition[index].pop_back();
		}
		currentPartition.push_back(std::vector<int> { i });
		genPartitions(n, i + 1, currentPartition);
		currentPartition.pop_back();
	}
	std::vector<std::vector<std::vector<std::vector<int>>>> partitions_;
public:
	size_t size(int n) const {
		return partitions_[n].size();
	}
	auto begin(int n) {
		return partitions_[n].begin();
	}
	auto end(int n) {
		return partitions_[n].end();
	}
	auto begin(int n) const {
		return partitions_[n].begin();
	}
	auto end(int n) const {
		return partitions_[n].end();
	}
	std::vector<std::vector<std::vector<int>>> const& get(int n) const {
		return partitions_[n];
	}
	Partitions(int N) {
		partitions_.resize(N + 1);
		std::vector<std::vector<int>> workspace;
		for (int n = 0; n <= N; n++) {
			genPartitions(n, 1, workspace);
		}
	}
};

template<int maxOrder>
struct BellPolynomials {
	MultivariatePolynomial<double, maxOrder + 1>& operator()(int n, int k) {
		return B_[n * (n + 1) / 2 + k];
	}
	std::array<MultivariatePolynomial<double, maxOrder + 1>, (maxOrder + 2) * (maxOrder + 1) / 2> B_;
public:
	BellPolynomials() {
		auto &B = *this;
		std::array<int, maxOrder + 1> index;
		MultivariatePolynomial<double, maxOrder + 1> zero { };
		MultivariatePolynomial<double, maxOrder + 1> one;
		std::array<MultivariatePolynomial<double, maxOrder + 1>, maxOrder + 1> Xi;
		index.fill(0);
		one[index] = 1.0;
		for (int n = 0; n <= maxOrder; n++) {
			index[n] = 1;
			Xi[n][index] = 1.0;
			index[n] = 0;
		}
		B(0, 0) = one;
		for (int n = 1; n <= maxOrder; n++) {
			B(n, 0) = zero;
			B(0, n) = zero;
		}
		for (int n = 0; n < maxOrder; n++) {
			for (int k = 0; k <= n; k++) {
				for (int i = 0; i <= n - k; i++) {
					B(n + 1, k + 1) += binco(n, i) * Xi[i + 1] * B(n - i, k);
				}
			}
		}

	}
	MultivariatePolynomial<double, maxOrder + 1> const& operator()(int n, int k) const {
		return B_[n * (n + 1) / 2 + k];
	}
};

// Α,α,Β,β,Γ,γ,Δ,δ,Ε,ε,Ζ,ζ,Η,η,Θ,θ,ϑ,Ι,ι,Κ,κ,ϰ,Λ,λ,Μ,μ,Ν,ν,Ξ,ξ,Ο,ο,Π,π,ϖ,Ρ,ρ,ϱ,Σ,σ,ς,Τ,τ,Υ,υ,Φ,φ,ϕ,Χ,χ,Ψ,ψ,Ω,ω
template<typename Type, int maxOrder, int varCount = 1>
class Jet {
	using ElementType = typename std::conditional<varCount == 1, Type, Jet<Type, maxOrder, varCount - 1>>::type;
	std::array<ElementType, maxOrder + 1> d_;
public:
	Jet() {
		d_.fill(Type(0));
	}
	Jet& operator=(Type cons) {
		d_.fill(Type(0));
		d_[0] = cons;
		return *this;
	}
	Type const& operator[](int i) const {
		return d_[i];
	}
	Type& operator[](int i) {
		return d_[i];
	}
	friend Jet faàDiBruno(std::function<Type(Type const&, int)> const &F, Jet const &g) {
		BellPolynomials<maxOrder> Bnk;
		Jet<Type, maxOrder + 1> dFdg;
		Jet dFdx;
		std::array<double, maxOrder + 1> x = g.d_;
		x[0] = 0.0;
		dFdx.d_.fill(0);
		for (int n = 0; n <= maxOrder; n++) {
			dFdg[n] = F(g[0], n);
		}
		dFdx[0] = dFdg[0];
		for (int n = 1; n <= maxOrder; n++) {
			dFdx[n] = 0.0;
			for (int k = 1; k <= n; k++) {
				dFdx[n] += dFdg[k] * Bnk(n, k)(x);
				fflush(stdout);
			}
		}
		return dFdx;
	}
	static Jet inverse(Jet const &G) {
		return faàDiBruno([](Type const &x, int n) {
			if (n == 0) {
				return 1.0 / x;
			} else {
				double const co = factorial(n + 1) * nonepow(n);
				double const xn = pow(x, -(n + 1));
				return co * xn;
			}
		}, G);
	}
	friend Jet operator+(Jet const &F) {
		return F;
	}
	friend Jet operator-(Jet const &F) {
		Jet const P;
		for (int n = 0; n <= maxOrder; n++) {
			P[n] = -F[n];
		}
		return P;
	}
	friend Jet operator+(Jet const &F, Jet const &G) {
		Jet const P;
		for (int n = 0; n <= maxOrder; n++) {
			P[n] = F[n] + G[n];
		}
		return P;
	}
	friend Jet operator-(Jet const &F, Jet const &G) {
		Jet const P;
		for (int n = 0; n <= maxOrder; n++) {
			P[n] = F[n] - G[n];
		}
		return P;
	}
	friend Jet operator*(Jet const &F, Jet const &G) {
		Jet P;
		for (int n = 0; n <= maxOrder; n++) {
			P[n] = 0.0;
			for (int k = 0; k <= n; k++) {
				P[n] += binco(n, k) * F[k] * G[n - k];
			}
		}
		return P;
	}
	friend Jet operator/(Jet const &B, Jet const &C) {
		auto iC = inverse(C);
		return B * iC;
	}
	Jet& operator+=(Jet const &F) {
		*this = *this + F;
		return *this;
	}
	Jet& operator-=(Jet const &F) {
		*this = *this - F;
		return *this;
	}
	Jet& operator*=(Jet const &F) {
		*this = *this * F;
		return *this;
	}
	Jet& operator/=(Jet const &F) {
		*this = *this / F;
		return *this;
	}
	friend Jet exp(Jet const &G) {
		return faàDiBruno([](Type const &x, int n) {
			return exp(x);
		}, G);
	}
	static Jet genVar(Type value, int dim = 0) {
		Jet x;
		if constexpr (varCount > 1) {
			using NextJet = Jet<Type, maxOrder, varCount - 1>;
			if (dim > 0) {
				x[0] = NextJet::genVar(value, dim - 1);
				return x;
			}
		}
		x[0] = value;
		x[1] = 1.0;
		return x;
	}
	friend std::ostream& operator<<(std::ostream &os, Jet const &j) {
		os << "{";
		for (int n = 0; n < maxOrder; n++) {
			os << j[n] << ", ";
		}
		os << j[maxOrder] << "}\n";
		return os;
	}
};
