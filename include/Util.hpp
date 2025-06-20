#pragma once

#include <array>
#include <cmath>
#include <filesystem>
#include <string>

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
inline constexpr std::array<T, D> zero() {
	std::array<T, D> u;
	u.fill(T(0));
	return u;
}

template<int D>
inline constexpr std::array<int, D> unit(int d) {
	auto u = zero<int, D>();
	u[d] = 1;
	return u;
}

void installFpeHandler();

void toFile(std::string const &content, std::filesystem::path const &filePath);


