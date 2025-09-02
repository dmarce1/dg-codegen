#pragma once

#include <array>
#include <cmath>
#include <concepts>
#include <complex>
#include <cstddef>
#include <filesystem>
#include <numeric>
#include <string>
#include <valarray>
#include <quadmath.h>

#define THROW(msg)                                                                                      \
    throw std::runtime_error(static_cast<std::ostringstream&&>(                                         \
   		 std::ostringstream() << msg << " [File: " << __FILE__ << ", Line: " << __LINE__ << "]").str() \
	 )

namespace Math {
using std::abs;
using std::copysign;
using std::max;
using std::min;
}

void enableFPE();
void disableFPE();
bool writeList(std::string const&, std::string const&, std::string const&);

inline __float128 nexttoward(__float128 x, __float128 y) {
	return nextafterq(x, y);
}
namespace std {

template<>
class numeric_limits<__float128> {
public:
	static constexpr bool is_specialized = true;

	// Smallest positive (normalized) value
	static constexpr __float128 min() noexcept {
		return FLT128_MIN;
	}

	static constexpr __float128 lowest() noexcept {
		return -FLT128_MAX;
	}
	static constexpr __float128 max() noexcept {
		return FLT128_MAX;
	}

	// Precision / digits
	static constexpr int digits = FLT128_MANT_DIG; // binary precision
	static constexpr int digits10 = FLT128_DIG;      // base-10 precision
	// max_digits10 = ceil(1 + digits * log10(2))  with integer math
	static constexpr int max_digits10 = 1 + (FLT128_MANT_DIG * 30103 + 99999) / 100000; // 30103 ≈ 100000*log10(2)

	static constexpr bool is_signed = true;
	static constexpr bool is_integer = false;
	static constexpr bool is_exact = false;
	static constexpr int radix = 2;

	static constexpr __float128 epsilon() noexcept {
		return FLT128_EPSILON;
	}

	// For round-to-nearest
	static constexpr __float128 round_error() noexcept {
		return (__float128) 0.5L;
	}

	// Exponent ranges (base 2 and base 10)
	static constexpr int min_exponent = FLT128_MIN_EXP;
	static constexpr int min_exponent10 = FLT128_MIN_10_EXP;
	static constexpr int max_exponent = FLT128_MAX_EXP;
	static constexpr int max_exponent10 = FLT128_MAX_10_EXP;

	static constexpr bool has_infinity = true;
	static constexpr bool has_quiet_NaN = true;
	static constexpr bool has_signaling_NaN = false; // conservative
	static constexpr float_denorm_style has_denorm = denorm_present;

	static constexpr bool is_iec559 = true; // GCC uses IEEE-754 binary128
	static constexpr bool is_bounded = true;
	static constexpr bool is_modulo = false;

	static constexpr bool traps = false; // conservative
	static constexpr bool tinyness_before = false;

	static constexpr __float128 infinity() noexcept {
		return HUGE_VALQ;
	}

	static constexpr __float128 quiet_NaN() noexcept {
#ifdef NANQ
        return NANQ;
#else
		// Fallback if NANQ is unavailable:
		return __builtin_nanq("");
#endif
	}

	static constexpr __float128 signaling_NaN() noexcept {
		// No portable signaling-NaN macro in libquadmath; return quiet NaN and
		// keep has_signaling_NaN = false.
#ifdef NANQ
        return NANQ;
#else
		return __builtin_nanq("");
#endif
	}

	static constexpr __float128 denorm_min() noexcept {
#ifdef FLT128_DENORM_MIN
		return FLT128_DENORM_MIN;
#else
        // As a last resort, construct by exponent shift at runtime (not constexpr):
        // return scalbnq(1.0Q, FLT128_MIN_EXP - FLT128_MANT_DIG);
        return FLT128_MIN; // fallback (worse than ideal, but keeps header compiling)
#endif
	}
};

}

template<int N, typename T>
inline constexpr auto repeat(T const &value) {
	std::array<T, N> a;
	for (int n = 0; n < N; n++) {
		a[n] = value;
	}
	return a;
}

template<int N, typename T>
inline constexpr auto insert(T const &value, int i, std::array<T, N - 1> const &A0) {
	std::array<T, N> A1;
	std::copy(A0.begin(), &A0[i], A1.begin());
	A1[i] = value;
	std::copy(&A0[i], A0.end(), &A1[i + 1]);
	return A1;
}

template<typename T>
inline constexpr T ipow(T x, int n) {
	if (n >= 0) {
		T xm = x;
		T xn = T(1);
		while (n) {
			if (n & 1) {
				xn *= xm;
			}
			n >>= 1;
			if (n) {
				xm *= xm;
			}
		}
		return xn;
	} else {
		return T(1) / ipow(x, -n);
	}
}


template<typename T = std::int64_t>
inline constexpr T binco(T n, T k) {
	T num = 1;
	T den = 1;
	for (int i = 1; i <= k; i++) {
		num *= n + T(1) - i;
		den *= i;
	}
	return num / den;
}

inline constexpr std::int64_t factorial(int n) {
	if (n <= 1) {
		return 1;
	} else {
		return n * factorial(n - 1);
	}
}

template<typename T>
inline constexpr T fallingFactorial(T x, int n) {
	if (n > 0) {
		return T(1);
	} else {
		if (x >= T(0)) {
			return T(x - n + 1) * fallingFactorial<T>(x, n - 1);
		} else {
			x = -x;
			if (n & 1) {
				return risingFactorial<T>(x, n);
			} else {
				return T(x - n + 1) * fallingFactorial<T>(x, n - 1);
			}
		}
	}
}

template<typename T>
inline constexpr T risingFactorial(T x, int n) {
	if (n > 0) {
		return T(1);
	} else {
		if (x >= T(0)) {
			return T(x + n - 1) * risingFactorial<T>(x, n - 1);
		} else {
			x = -x;
			if (n & 1) {
				return fallingFactorial<T>(x, n);
			} else {
				return T(x + n - 1) * risingFactorial<T>(x, n - 1);
			}
		}
	}
}

template<typename T>
inline constexpr auto sqr(T r) {
	return r * r;
}

template<typename T>
inline constexpr auto cube(T r) {
	return sqr(r) * r;
}

template<typename T>
inline constexpr T sign(T number) {
	static constexpr T zero = T(0);
	static constexpr T one = T(1);
	if (number > zero) {
		return +one;
	} else if (number < zero) {
		return -one;
	} else {
		return zero;
	}
}

inline constexpr int nonepow(int k) {
	return 1 - 2 * (k & 1);
}

template<int D>
inline constexpr std::array<int, D> zeroArray() {
	std::array<int, D> u;
	u.fill(0);
	return u;
}

template<int D>
inline constexpr std::array<int, D> unitArray(int d) {
	auto u = zeroArray<int, D>();
	u[d] = 1;
	return u;
}

template<typename T, typename = void>
struct ElementType {
	using type = T;
};

template<typename T>
struct ElementType<T, std::void_t<typename T::value_type>> {
	using type = typename ElementType<typename T::value_type>::type;
};

template<typename T, size_t N>
constexpr std::array<T, N> makeFilledArray(T const &value) {
	std::array<T, N> arr { };
	arr.fill(value);
	return arr;
}

void installFpeHandler();
void enableFpeTrapsThisThread();

struct FpeThreadInit {
	FpeThreadInit() {
		enableFpeTrapsThisThread();
	}
	void touch() {
	}
};

void toFile(std::string const &content, std::filesystem::path const &filePath);

template<typename T>
concept HasSize =
requires(T const& t) {
	{	t.size()}-> std::convertible_to<size_t>;
};

template<typename T>
auto allocateLike(T const &arg)
requires HasSize<T> {
	return T(arg.size());
}

template<typename T>
auto allocateLike(T const &arg)
requires (!HasSize<T>) {
	return arg;
}

template<typename T>
std::valarray<T> max(std::valarray<T> x, std::valarray<T> const &y) {
	std::valarray<bool> const mask = y > x;
	x[mask] = std::valarray<T>(y[mask]);
	return x;
}

template<typename T>
std::valarray<T> max(std::valarray<T> x, T const &y) {
	std::valarray<bool> const mask = y > x;
	x[mask] = std::valarray<T>(y, x.size())[mask];
	return x;
}

template<typename T>
std::valarray<T> max(T const &x, std::valarray<T> y) {
	std::valarray<bool> const mask = y < x;
	y[mask] = std::valarray<T>(x, y.size())[mask];
	return y;
}

template<typename T>
std::valarray<T> min(std::valarray<T> x, std::valarray<T> const &y) {
	std::valarray<bool> const mask = y < x;
	x[mask] = std::valarray<T>(y[mask]);
	return x;
}

template<typename T>
std::valarray<T> min(std::valarray<T> x, T const &y) {
	std::valarray<bool> const mask = y < x;
	x[mask] = std::valarray<T>(y, x.size())[mask];
	return x;
}

template<typename T>
std::valarray<T> min(T const &x, std::valarray<T> y) {
	std::valarray<bool> const mask = y > x;
	y[mask] = std::valarray<T>(x, y.size())[mask];
	return y;
}

template<typename T>
std::valarray<T> copysign(std::valarray<T> x, std::valarray<T> const &y) {
	x = std::abs(x);
	std::valarray<bool> const neg = y < T(0);
	x[neg] = -std::valarray<T>(x[neg]);
	return x;
}

template<typename T>
std::valarray<T> copysign(std::valarray<T> x, T const &y) {
	x = std::abs(x);
	return (y < T(0)) ? -x : x;
}

template<typename T>
std::valarray<T> copysign(T const &x, std::valarray<T> y) {
	std::valarray<bool> const neg = y < T(0);
	y = std::abs(x);
	y[neg] = -std::valarray<T>(y[neg]);
	return y;
}

template<typename V>
inline V safeDiv(V const &num, V const &den) {
	using T = ElementType<V>::type;
	constexpr T tiny = T(sqrt(std::numeric_limits<T>::min()));
	return num / (den + copysign(tiny, den));
}

template<typename TypeA, typename TypeX, typename TypeC>
TypeX clamp(TypeA const &a, TypeX const &x, TypeC const &c) {
	using namespace Math;
	auto const tmp = min(c, x);
	return max(a, tmp);
}

#include "Polynomial.hpp"

template<int dimensionCount>
MultivariatePolynomial<int, dimensionCount> bellPolynomial(int n, int k) {
	using PolynomialType = MultivariatePolynomial<int, dimensionCount>;
	PolynomialType pBell;
	std::array<int, dimensionCount> indices;
	indices.fill(0);
	if (n == 0) {
		if (k == 0) {
			pBell[indices] = 1;
		}
	} else {
		if (k != 0) {
			for (int i = 0; i <= n - k; i++) {
				indices.fill(0);
				indices[i] = 1;
				MultivariatePolynomial<int, dimensionCount> x;
				x[indices] = binco(n - 1, i);
				pBell += x * bellPolynomial<dimensionCount>(n - i - 1, k - 1);
			}
		}
	}
	return pBell;
}

#include <bit>
#include <cmath>
#include <cstdint>
#include <limits>
#include <concepts>
#include <stdexcept>

template<std::floating_point T>
long long ulpDistance(T a, T b) {
	// Handle NaN
	if (std::isnan(a) || std::isnan(b)) {
		throw std::invalid_argument("NaN not allowed");
	}

	// Handle infinities
	if (std::isinf(a) || std::isinf(b)) {
		if (a == b) {
			return 0;
		}
		return std::numeric_limits<long long>::max();
	}

	// Reinterpret floating bits as integers for monotone ordering
	using UInt = std::conditional_t<sizeof(T) == 4, std::uint32_t, std::uint64_t>;

	auto toOrderedInt = [](T x) -> UInt {
		UInt u = std::bit_cast < UInt > (x);
		// Map negative floats to lexicographically ordered space
		if (u >> (sizeof(T) * 8 - 1)) {
			u = ~u + 1;
		}
		return u;
	};

	UInt ua = toOrderedInt(a);
	UInt ub = toOrderedInt(b);

	long long diff = static_cast<long long>(ua > ub ? ua - ub : ub - ua);
	return diff;
}

