#pragma once

#include <cmath>
#include <filesystem>

namespace Math {
using std::abs;
using std::copysign;
using std::max;
using std::min;
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
inline constexpr auto sqr(T r) {
	return r * r;
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


template<typename T, typename = void>
struct ElementType {
	using type = T;
};

template<typename T>
struct ElementType<T, std::void_t<typename T::value_type>> {
	using type = typename ElementType<typename T::value_type>::type;
};


void toFile(std::string const &content, std::filesystem::path const &filePath);

