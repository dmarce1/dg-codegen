#pragma once

#include <array>
#include <cmath>
#include <filesystem>
#include <numeric>
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

template<typename V>
inline V safeDiv(V const &num, V const &den) {
	using T = ElementType<V>::type;
	constexpr T tiny = T(sqrt(std::numeric_limits<T>::min()));
    return num / (den + copysign(tiny, den));
}

template<typename TypeA, typename TypeX, typename TypeC>
TypeX clamp(TypeA const& a, TypeX const& x, TypeC const& c) {
	return max(a, min(c, x));
}

template<typename T, size_t N>
constexpr std::array<T, N> makeFilledArray(T const& value) {
    std::array<T, N> arr{};
    arr.fill(value);
    return arr;
}

constexpr auto allSevens = makeFilledArray<int, 4>(7);    // {7,7,7,7}

void installFpeHandler();

void toFile(std::string const &content, std::filesystem::path const &filePath);

//template<typename T1, typename T2, typename = void>
//struct Multiplies {
//	static constexpr bool value = false;
//};
//
//template<typename T1, typename T2>
//struct Multiplies<T1, T2, std::void_t<decltype(std::declval<T1>() * std::declval<T2>())>>  {
//	static constexpr bool value = true;
//};
//
//template<typename T, auto D>
//std::array<std::array<T, D>, D> gramSchmidt(std::array<T, D> v) {
//	using U = typename ElementType<T>::type;
//	constexpr U half = U(1) / U(2);
//	std::array<std::array<T, D>, D> u;
//	std::array<T, D> k;
//	std::iota(k.begin(), k.end(), T(0));
//	auto const theta = dot(v, v);
//	for(int d  = 0; d < int(D - 1); d++) {
//		v[d] *= theta;
//	}
//	v.back() = (U(1) - theta) + theta * v.back();
//	u[0] = v;
//	for (int n = 0; n < D; n++) {
//		auto &un = u[n];
//		auto const u0 = un;
//		for (int m = 0; m < n; m++) {
//			un -= u[m][k[n]] * u[m];
//		}
//		un = normalize(un);
//	}
//	return u;
//}
//
