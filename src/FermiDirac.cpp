#include "AutoDifferentiation.hpp"
#include "PadeApproximant.hpp"
#include <cmath>
#include <limits>

#include "Integrate.hpp"
#include "HalfInteger.hpp"
#include "FermiDirac.hpp"
#include "Polynomial.hpp"
#include "DifferentialEquation.hpp"

static constexpr int nInterp = 3;
// ρ = g(η(ρ, β, β)
// δρ/δη = δg(η(ρ, β, β)/δη δη/
// Α,α,Β,β,Γ,γ,Δ,δ,Ε,ε,Ζ,ζ,Η,η,Θ,θ,ϑ,Ι,ι,Κ,κ,ϰ,Λ,λ,Μ,μ,Ν,ν,Ξ,ξ,Ο,ο,Π,π,ϖ,Ρ,ρ,ϱ,Σ,σ,ς,Τ,τ,Υ,υ,Φ,φ,ϕ,Χ,χ,Ψ,ψ,Ω,ω
static double constexpr infinity = std::numeric_limits<double>::infinity();
static HalfInteger<int> constexpr half = HalfInteger<int>::half();
static double constexpr c = 2.997924580000000e+10;
static double constexpr h = 6.626070150000000e-27;
static double constexpr kB = 3.806490000000000e-16;
static double constexpr me = 9.109383701500000e-28;
static double constexpr mH = 1.67353271600000e-24;
static double constexpr π = 4.0 * atan(1.0);
static int constexpr dCount = 3;


constexpr int stirling1(int n, int k) {
	if (k == 0) {
		return n == 0;
	}
	int const sum = k * stirling1(n - 1, k) + stirling1(n - 1, k - 1);
	return sum;
}

constexpr int stirling2(int n, int k) {
	if ((n == 0) || (k == 0) || (k > n)) {
		return 0;
	}
	if (n == k) {
		return 1;
	}
	int const sum = k * stirling2(n - 1, k) + stirling2(n - 1, k - 1);
	return sum;
}

double polylog(int n, double z) {
	constexpr double infinity = std::numeric_limits<double>::infinity();
	constexpr double δz = pow(std::numeric_limits<double>::epsilon(), 1.0 / 3.0);
	if (n < 0) {
		double sum = 0.0;
		double t = z / (1 - z);
		double tk = 1;
		for (int k = 0; k <= -n; k++) {
			tk *= t;
			sum += factorial(k) * stirling2(-n + 1, k + 1) * tk;
		}
		return nonepow(n + 1) * sum;
	} else if (n == 0) {
		if (z == -infinity) {
			return -1.0;
		}
		if (-z < δz) {
			double y = 1.0;
			y = fma(z, y, 1.0);
			y = fma(z, y, 1.0);
			y = fma(z, y, 1.0);
			y *= z;
			return y;
		}
		double const t = -1.0 / z;
		if (t < δz) {
			double y = 1.0;
			y = fma(t, y, -1.0);
			y = fma(t, y, +1.0);
			y = fma(t, y, -1.0);
			return y;
		}
		return z / (1 - z);
	} else {
		assert(false);
	}
	return 1.0 / 0.0;
}

double logisticCDF(double z) {
	return -polylog(0, -z);
}

//std::function<double(double)> compositionDerivative(std::function<double(double)> const &g, std::function<double(double)> const &h, double x, int n, int i) {
//	if (i == 0) {
//		return [n](double) {
//			return (n == 0);
//		};
//	}
//	auto const childA = compositionDerivative(g, h, x, n - 1, i);
//	auto const childB = compositionDerivative(g, h, x, n - 1, i - 1);
//	std::function<double(double)> f;
//	return f;
//}

template<int N>
std::array<double, N + 1> faáDiBruno(std::function<double(double, int)> const &f, Polynomial<double> g, double x) {
	static BellPolynomials<N> Bnk;
	std::array<double, N + 1> dfdx;
	std::array<double, N> dgdx;
	for (int n = 0; n < N; n++) {
		dgdx[n] = g(x);
		g = polynomialDerivative(g);
	}
	dfdx[0] = f(g(x), 0);
	for (int m = 1; m <= N; m++) {
		dfdx[m] = 0.0;
		for (int n = 0; n <= m; n++) {
			for (int k = 0; k <= n; k++) {
				dfdx[m] += binco(n, k) * f(g(x), k) * Bnk(n, k)(dgdx);
			}
		}
	}
	return dfdx;
}

Polynomial<double> polynomialPower(Polynomial<double> const &p, int m) {
	Polynomial<double> pn;
	Polynomial<double> pm = p;
	pn[0] = 1;
	while (m) {
		if (m & 1) {
			pn *= pm;
		}
		m >>= 1;
		if (m) {
			pm *= pm;
		}
	}
	return pn;
}

double FD(HalfInteger<int> const &k, double η, double β, double α = 1.0);

double Nele(double η, double β, int dηN = 0, int dβN = 0) {
	static double const N0 = 8.0 * sqrt(2.0) * π * ipow(me * c / h, 3);
	return N0 * FD(half, η, β);
}

double Npos(double η, double β, int dηN = 0, int dβN = 0) {
	static double const N0 = 8.0 * sqrt(2.0) * π * ipow(me * c / h, 3);
	η += 2 / β;
	return N0 * FD(half, -η, β);
}

double Pele(double η, double β, int dηN = 0, int dβN = 0) {
	static double const N0 = 8.0 * sqrt(2.0) * π * ipow(me * c / h, 3);
	return N0 * FD(3 * half, η, β);
}

double Ppos(double η, double β, int dηN = 0, int dβN = 0) {
	static double const N0 = 8.0 * sqrt(2.0) * π * ipow(me * c / h, 3);
	η += 2 / β;
	return N0 * FD(3 * half, -η, β);

}

auto Derivative(int n, int k) {
	return [n, k](std::function<double(double, double, int, int)> const &F) {
		return [n, k, F](double x, double y) {
			return F(x, y, n, k);
		};
	};
}

#define Power(a, b) ipow(a, b)
#define List(...) {__VA_ARGS__}

//double F12(double η, double β, int n = 0, int k = 0) {
//	return FD(half, η, β, n, k);
//}
//;
//
//double F32(double η, double β, int n = 0, int k = 0) {
//	return FD(3 * half, η, β, n, k);
//}
//;
//
//double F52(double η, double β, int n = 0, int k = 0) {
//	return FD(5 * half, η, β, n, k);
//}
//;

double solve4η(double ρ, double μe, double T) {

	constexpr double mH = 1.67353271600000e-24;
	constexpr double To = me * sqr(c) / kB;
	double const β = T / To;
	double const N = ρ / (μe * mH);
	auto F = [N, β](double η) {
		return N + Npos(η, β) - Nele(η, β);
	};
	double η = 0.0;
	double ηa = -1.0;
	double ηb = 1.0;
	while (F(ηa) * F(ηb) > 0.0) {
		ηa *= 10.0;
		ηb *= 10.0;
	}
	while (ηb > std::nextafter(ηa, ηb)) {
		η = 0.5 * (ηa + ηb);
		if (F(ηa) * F(η) > 0.0) {
			ηa = η;
		} else {
			ηb = η;
		}
	}
	return η;
}

double helmholtz(double ρ, double μe, double T) {
	constexpr double To = me * sqr(c) / kB;
	double const n = ρ / (μe * mH);
	double const β = T / To;
	const double η = solve4η(ρ, μe, T);
	const double μ = β * η;
	const double P = Pele(η, β) - Ppos(η, β);

	return μ - P / n;
}
//struct Adaptivηble {
//	Adaptivηble() {
//		ready = false;
//		leaf = true;
//		zmin.fill(-1);
//		zmax.fill(+1);
//		bmin.fill(true);
//		bmax.fill(true);
//	}
//	bool isBound() {
//		for (int d = 0; d < D; d++) {
//			if (bmin[d]) {
//				return true;
//			}
//			if (bmax[d]) {
//				return true;
//			}
//		}
//		return false;
//	}
//	void refine(int d) {
//		chilren.resize(2);
//		Adaptivηble l = *this;
//		Adaptivηble r = *this;
//		double const zmid = 0.5 * (zmin[d] + zmax[d]);
//		l.zmax[d] = r.zmin[d] = zmid;
//		l.bmax[d] = r.bmax[d] = false;
//		leaf = false;
//	}
//	double lookup(std::array<double, D> const &z) {
//
//	}
//	double operator()(std::array<double, D> const &x) {
//		std::array<double, D> z;
//		for (int d = 0; d < D; d++) {
//			z[d] = tanh(sinh(x / x0));
//		}
//		return lookup(z);
//	}
//private:
//	bool leaf;
//	bool ready;
//	std::vector<Adaptivηble> children;
//	std::array<bool, D> bmin;
//	std::array<bool, D> bmax;
//	std::array<double, D> zmin;
//	std::array<double, D> zmax;
//	std::array<double, ipow(2 * P, D)> C;
//};

#include <bits/stdc++.h>

void printPartition(std::vector<std::vector<int>> const &p) {
	std::cout << "{ ";
	for (size_t i = 0; i < p.size(); ++i) {
		std::cout << "{";
		for (size_t j = 0; j < p[i].size(); ++j) {
			std::cout << p[i][j];
			if (j + 1 < p[i].size()) {
				std::cout << ",";
			}
		}
		std::cout << "}";
		if (i + 1 < p.size()) {
			std::cout << ", ";
		}
	}
	std::cout << " }" << '\n';
}


//template<int N, int D, int M>
//double faàDiBruno(std::function<double(std::array<double, M>, std::array<int, M> n)>const &f,
//		std::array<std::function<double(std::array<double, D>, std::array<int, D>)>, M> const &g, std::array<double, D> x) {
//	double sum = 0.0;
//	printf("F\n");
//	auto const incrementIndices = [](std::array<int, D> &I) {
//		int d;
//		d = 0;
//		while (I[d]++ == N) {
//			if (d == D) {
//				break;
//			}
//			I[d] = 0;
//			d++;
//		}
//		if (d == D) {
//			return -1;
//		}
//		return std::accumulate(I.begin(), I.end(), 0);
//	};
//	for (int m = 0; m < M; m++) {
//		std::array<int, D> δD;
//		std::array<int, 2 * N + 1> δN;
//		δD.fill(0);
//		δN.fill(0);
//		static Partitions<2 * N + 1> Π { };
//		int isum = incrementIndices(δD);
//		while (isum >= 0) {
//			printf("\n");
//			for (int n = 0, d = 0; d < D; d++) {
//				for (int n0 = n; n < n0 + δD[d]; n++) {
//					δN[n] = d;
//				}
//			}
//			std::array<int, D> dn;
//			for (auto const &ϖ : Π.get(isum)) {
//				int const Cϖ = ϖ.size();
//				double dngdxn = 1.0;
//				bool first = true;
//				for (auto const &B : ϖ) {
//					dn.fill(0);
//					for (auto j : B) {
//						dn[δN[j - 1]]++;
//					}
//					if (!first) {
//						printf(" * ");
//					}
//					printf("η^(");
//					for (unsigned i = 0; i < dn.size(); i++) {
//						printf("%i", dn[i]);
//						if (i + 1 < dn.size()) {
//							printf(",");
//						}
//					}
//					printf(") ");
//					//		dngdxn *= g(x, dn);
//					first = false;
//				}
//				printf(" * F^%i\n", Cϖ);
//				dn.fill(0);
//				sum += f(g(x, dn), Cϖ) * dngdxn;
//			}
//			isum = incrementIndices(δD);
//		}
//	}
//	return sum;
//}

template<typename T, int N>
struct PowerSeries {
	PowerSeries(std::function<T(int)> const &f) {
		for (int n = 0; n <= N; n++) {
			an[N - n] = f(n);
		}
	}
	T operator()(double x) const {
		double y = an[0];
		for (int n = 1; n <= N; n++) {
			y = x * y + an[n];
		}
		return y;
	}
private:
	std::array<double, N + 1> an;
};

//std::complex<double> iBeta(double x, HalfInteger<int> a, HalfInteger<int> b) {
//	constexpr std::complex<double> i(0, 1);
//	int const n = (2 * a - 1).toInteger();
//	int const m = (2 * b - 1).toInteger();
//	std::complex<double> const z = asin(sqrt(std::complex<double>(x)));
//	std::complex<double> sum = 0.0;
//	for (int k = 0; k <= n; k++) {
//		double const nk = binco(n, k);
//		for (int l = 0; l <= m; l++) {
//			double const ml = binco(m, l);
//			int const q = m + n - 2 * (k + l);
//			std::complex<double> num, den;
//			if (q) {
//				num = std::exp(i * double(q) * z) - 1.0;
//				den = i * double(q);
//			} else {
//				num = z;
//				den = 1.0;
//			}
//			num *= nk * ml;
//			den *= double(1 << (m + n)) * ipow(i, m) * ipow(-1.0, k);
//			sum += num / den;
//		}
//	}
//	return 2.0 * sum;
//}

std::complex<double> iBeta(double x, HalfInteger<int> a, HalfInteger<int> b) {
	constexpr std::complex<double> i(0, 1);
	auto const denCo = [a, b, x](int m) {
		if (m == 0) {
			return 1.0;
		}
		if (m & 1) {
			m = (m - 1) >> 1;
			auto const num = -double(a + m) * double(a + b + m) * x;
			auto const den = double(a + 2 * m) * double(a + 2 * m + 1);
			return num / den;
		}
		m >>= 1;
		auto const num = m * (b - m) * x;
		auto const den = double(a + 2 * m - 1) * double(a + 2 * m);
		return num / den;
	};
	int constexpr N = 32;
	double den = 1;
	for (int n = N; n > 0; n--) {
		den = denCo(n) / den;
		den += 1;
	}
	den *= a;
	std::complex<double> z(x, 0);
	std::complex<double> const loga = log(z);
	std::complex<double> const logb = log(std::complex<double>(1) - z);
	return exp(double(a) * loga) * exp(double(b) * logb) / den;
}

double FD(HalfInteger<int> const &k, double η, double β, double α) {
	using std::pow;
	using std::sqrt;
	constexpr int qOrder = 24;
	static auto const gaussLaguerre = gaussLaguerreQuadraturePoints<double>(qOrder);
	static auto const gaussLegendre = gaussLegendreQuadraturePoints<double>(qOrder);
	constexpr double D = 3.3609;
	constexpr double σ = 9.1186e-2;
	constexpr double a1 = 6.7774;
	constexpr double b1 = 1.1418;
	constexpr double c1 = 2.9826;
	constexpr double a2 = 3.7601;
	constexpr double b2 = 9.3719e-2;
	constexpr double c2 = 2.1063e-2;
	constexpr double d2 = 3.1084e1;
	constexpr double e2 = 1.0056;
	constexpr double a3 = 7.5669;
	constexpr double b3 = 1.1695;
	constexpr double c3 = 7.5416e-1;
	constexpr double d3 = 6.6558;
	constexpr double e3 = -1.2819e-1;
	auto const Exp = [](double x) {
		constexpr double zero = 0.0;
		constexpr double χmax = std::log(std::numeric_limits<double>::max());
		if (x > χmax) {
			return zero;
		} else if (x < -χmax) {
			return zero;
		} else {
			return std::exp(x);
		}
	};
	double const ξ = log1p(Exp(σ * (η - D))) / σ;
	double const Xa = ((a1 + b1 * ξ + c1 * ξ * ξ) / (1.0 + c1 * ξ));
	double const Xb = ((a2 + b2 * ξ + c2 * d2 * ξ * ξ) / (1.0 + e2 * ξ + c2 * ξ * ξ));
	double const Xc = ((a3 + b3 * ξ + c3 * d3 * ξ * ξ) / (1.0 + e3 * ξ + c3 * ξ * ξ));
	double const S1 = Xa - Xb;
	double const S2 = Xa;
	double const S3 = Xa + Xc;
	double const ha = sqrt(S1);
	double const hb = (S2 - S1);
	double const hc = (S3 - S2);
	KahnSum<double> Ia, Ib, Ic, Id;
	Ia = 0.0;
	Ib = 0.0;
	Ic = 0.0;
	Id = 0.0;
	for (int q = 0; q < qOrder; q++) {
		double const x = gaussLegendre[q].position;
		double const wp = 0.5 * gaussLegendre[q].weight;
		double const xd = (S3 + gaussLaguerre[q].position);
		double const wl = gaussLaguerre[q].weight;
		double const τ = 0.5 * (x + 1.0);
		double const xa = τ * ha;
		double const xb = (S1 + τ * hb);
		double const xc = (S2 + τ * hc);
		double const x2 = xa * xa;
		Ia += wp * pow(x2, k) * (1 + α * β * x2) * sqrt(1 + 0.5 * β * x2) / (1 + Exp(x2 - η)) * ha * 2.0 * xa;
		Ib += wp * pow(xb, k) * (1 + α * β * xb) * sqrt(1 + 0.5 * β * xb) / (1 + Exp(xb - η)) * hb;
		Ic += wp * pow(xc, k) * (1 + α * β * xc) * sqrt(1 + 0.5 * β * xc) / (1 + Exp(xc - η)) * hc;
		Id += wl * pow(xd, k) * (1 + α * β * xd) * sqrt(1 + 0.5 * β * xd) / (1 + Exp(xd - η)) * Exp(xd - S3);
	}
//	printf("%e %e %e %e\n", double(Ia), double(Ib), double(Ic), double(Id));
	return (Ia + Ib + Ic + Id);
}

#include <vector>
#include <cmath>
#include <stdexcept>
#include <limits>

template<class T>
auto computeY(int N, T const &x, T const &k) -> std::vector<T> {
	if (N < 0) {
		return {};
	}
	if (x == T(0) || x == T(-2)) {
		throw std::domain_error("Denominator (2+n)*x*(2+x) vanishes for x=0 or x=-2.");
	}
	if (x + T(2) <= T(0)) {
		throw std::domain_error("sqrt(2+x) requires x > -2.");
	}
	// Note: if x<0 and k is non-integer, std::pow(x,k) is not real.
	std::vector<T> y(static_cast<std::size_t>(N + 1));
	auto const sqrt2 = std::sqrt(T(2));
	auto const sqrt2px = std::sqrt(T(2) + x);

	y[0] = std::pow(x, k) * (sqrt2px / sqrt2);
	if (N == 0) {
		return y;
	}
	y[1] = std::pow(x, k - T(1)) * (x + T(2) * k * (T(2) + x)) / (T(2) * sqrt2 * sqrt2px);

	for (int n = 0; n <= N - 2; ++n) {
		auto const a = (-T(1) + T(2) * n - T(2) * k);
		auto const b = (T(4) + T(4) * n - T(4) * k + T(3) * x + T(4) * n * x - T(2) * k * x);
		auto const c = T(2) * (T(2) + n) * x * (T(2) + x);
		y[static_cast<std::size_t>(n + 2)] = -(b * y[static_cast<std::size_t>(n + 1)] + a * y[static_cast<std::size_t>(n)]) / c;
	}
	return y;
}

template<int O>
double fermiDirac(HalfInteger<int> s, double η, double β) {
	constexpr int N = 3;
	constexpr double α = 1.0;
	constexpr double π = M_PI;
	constexpr std::array<double, O + 1> Im = []() {
		std::array<double, O + 1> Im;
		for (int m = 0; m <= O; m++) {
			if (m % 2 == 0) {
				double B2 = 0.0;
				for (int k = 0; k <= m; k++) {
					for (int j = 0; j <= k; j++) {
						B2 += double(binco(k, j) * nonepow(j) * ipow(j + 1, m)) / double(k + 1);
					}
				}
				Im[m] = std::abs(((1 << m) - 2) * ipow(π, m) * B2);
			} else {
				Im[m] = 0.0;
			}
		}
		return Im;
	}();
	auto const fψ = [α, s](double ε) {
		constexpr std::complex<double> i(0, 1);
		constexpr std::complex<double> one(1, 0);
		constexpr std::complex<double> two(2, 0);
		auto const ψ = i * pow(two, s + 1) * ipow(-one, s) * (iBeta(-0.5 * ε, s + 1, 3 * half) - 2.0 * α * iBeta(-0.5 * ε, s + 2, 3 * half));
		return ψ.real();
	};
	std::array<double, O + 1> Yn;
	std::array<double, O + 1> ψn;
	constexpr double sqrt2 = sqrt(2);
	constexpr double isqrt2 = 1 / sqrt2;
	double const μ = η * β;
	auto yn = computeY(O, μ, 1.5);
	ψn[0] = fψ(μ);
	ψn[1] = ipow(μ, s) * (1 + α * μ) * sqrt(1 + 0.5 * μ);
	for (int n = 1; n < O; n++) {
		ψn[n + 1] = factorial(n) * (α * Yn[n] + (1 + α * μ) * Yn[n + 2]);
	}
	double sum = 0.0;
	sum = ψn[0];
	for (int m = 2; m <= O; m += 2) {
		auto const c = pow(β, s + 1 - m);
		printf("%e %e %e %e\n", sum, Im[m], ψn[m], c);
		sum += Im[m] * ψn[m] / factorial(m) * c;
	}
	return sum;
}


#include <cstdint>
auto test() {

//	std::complex<double> I(0, 1);
//	for (double x = 10.0; x < 101.0; x *= 1.1) {
//		printf("%e %e %e %e\n", x, (I * iBeta(-0.5 * x, 3 * half, 3 * half)).real(), (I * iBeta(-0.5 * x, 5 * half, 3 * half)).real(),
//				(I * iBeta(-0.5 * x, 7 * half, 3 * half)).real());
//	}
//	return;
	double β = 1e-6;
	HalfInteger<int> k = half * 3;
	PadeApproximant<16> pade([k](long double x, int n) {
		double const arr[33] = {66.03346919316417,151.83272672423314,372.5892539872859,928.6137153408156,
				   2319.424125055658,5795.195148371597,14502.08872272234,36254.118005093675,90315.18570453653,
				   226427.78660756798,575607.5146336989,1.4086630584417323e6,3.1082090920605077e6,
				   8.896740796184724e6,4.756492715955744e7,1.1399453026039888e8,-1.6890141681888695e9,
				   -1.2090409978887466e10,1.407866098280436e11,2.0760130019224863e12,-8.537240490744178e12,
				   -3.250372665992431e14,-3.264586214853177e14,4.93526789471258e16,3.2041012812560666e17,
				   -7.158018902280004e18,-1.0281900584085494e20,9.406948652791451e20,2.8185218104198913e22,
				   -9.150120925819329e22,-7.784494928026343e24,-6.148371657113851e24,2.392619084675169e27};
		return arr[n];
	}, std::exp(1));
	for (double η = 1; η < 100; η += 0.01) {
		auto fd1 = FD(k, η, 0.0);
		auto fd2 = pade(log(η)- 2);
		double e1 = fd2 - fd1;
		printf("%e %e %e %e\n", η, fd1, fd2, e1);
	}
	return;
	//fermiDirac(half, 1.0, 1.0);
	return;
	constexpr double s = 1.5;
	PowerSeries<double, 1024> polylog([s](int n) {
		if (n > 0) {
			return std::pow(n, -s);
		}
		return 0.0;
	});
	for (double z = -0.5; z < 0.5; z += 0.01) {
		printf("%e %e\n", z, polylog(z));
	}
	return;
	constexpr int N = 32;
	double ρa = 1e-12;
	double ρb = 1e+12;
	double dρ = (log(ρb) - log(ρa)) / N;
	double Ta = 1e3;
	double Tb = 1e12;
	double dT = (log(Tb) - log(Ta)) / N;
	for (double logρ = log(ρa); logρ <= log(ρb); logρ += dρ) {
		double const ρ = exp(logρ);
		for (double logT = log(Ta); logT <= log(Tb); logT += dT) {
			double const T = exp(logT);
			auto const η = solve4η(ρ, 2.0, T);
			printf("%i %e %e %e\n", N, exp(logρ), exp(logT), η);
			fflush(stdout);
		}
	}
	return;
}
;

void polytest() {
	constexpr int N = 8;
	auto f = Jet<double, N - 1>::genVar(1);
	std::cout << exp(f * f);
	return;
//	test();
}
