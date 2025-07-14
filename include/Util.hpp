#pragma once

#include <array>
#include <cmath>
#include <concepts>
#include <cstddef>
#include <filesystem>
#include <numeric>
#include <string>
#include <valarray>

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
	static constexpr T one = T(1);
	if (n >= 0) {
		T xm = x;
		T xn = one;
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
		return one / ipow(x, -n);
	}
}

template<typename T>
inline constexpr T binco(T n, T k) {
	static constexpr T one = T(1);
	T num = one;
	T den = one;
	for (int i = 1; i <= k; i++) {
		num *= T(n + 1 - i);
		den *= T(i);
	}
	return num / den;
}

template<typename T>
inline constexpr T factorial(int n) {
	static constexpr T one = T(1);
	if (n <= 1) {
		return one;
	} else {
		return T(n) * nFactorial < T > (n - 1);
	}
}

template<typename T>
inline constexpr T sqr(T r) {
	return r * r;
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

template<typename T, int D>
inline constexpr std::array<T, D> zeroArray() {
	std::array<T, D> u;
	u.fill(T(0));
	return u;
}

template<int D>
inline constexpr std::array<int, D> unit(int d) {
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

