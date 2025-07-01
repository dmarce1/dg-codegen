#pragma once

#include <type_traits>
#include <array>
#include <utility>

template<typename T>
struct CanDoArithmetic {
	static constexpr bool value = false;
};

template<typename T, auto N>
struct CanDoArithmetic<std::array<T, N>> {
	static constexpr bool value = true;
};


template<typename T, typename U, std::enable_if_t<CanDoArithmetic<T>::value && CanDoArithmetic<U>::value, int> = 0>
inline constexpr T& operator+=(T &A, U const &B) {
	int const cnt = A.size();
	for (int i = 0; i < cnt; i++) {
		A[i] += B[i];
	}
	return A;
}

template<typename T, std::enable_if_t<CanDoArithmetic<T>::value, int> = 0>
inline constexpr T& operator-=(T &A, T const &B) {
	int const cnt = A.size();
	for (int i = 0; i < cnt; i++) {
		A[i] -= B[i];
	}
	return A;
}

template<typename T, std::enable_if_t<CanDoArithmetic<T>::value, int> = 0>
inline constexpr T& operator*=(T &A, T const &B) {
	int const cnt = A.size();
	for (int i = 0; i < cnt; i++) {
		A[i] *= B[i];
	}
	return A;
}

template<typename T, std::enable_if_t<CanDoArithmetic<T>::value, int> = 0>
inline constexpr T& operator/=(T &A, T const &B) {
	int const cnt = A.size();
	for (int i = 0; i < cnt; i++) {
		A[i] /= B[i];
	}
	return A;
}

template<typename T, std::enable_if_t<CanDoArithmetic<T>::value, int> = 0>
inline constexpr T& operator%=(T &A, T const &B) {
	int const cnt = A.size();
	for (int i = 0; i < cnt; i++) {
		A[i] %= B[i];
	}
	return A;
}

template<typename T, std::enable_if_t<CanDoArithmetic<T>::value, int> = 0>
inline constexpr T operator+(T const &A) {
	return A;
}

template<typename T, std::enable_if_t<CanDoArithmetic<T>::value, int> = 0>
inline constexpr T operator-(T const &A) {
	auto B = A;
	for (auto &b : B) {
		b = -b;
	}
	return B;
}

template<typename T, typename U, std::enable_if_t<CanDoArithmetic<T>::value && CanDoArithmetic<U>::value, int> = 0>
inline constexpr T operator+(T const &A, U const &B) {
	auto C = A;
	C += B;
	return C;
}

template<typename T, std::enable_if_t<CanDoArithmetic<T>::value, int> = 0>
inline constexpr T operator-(T const &A, T const &B) {
	auto C = A;
	C -= B;
	return C;
}

template<typename T, typename U, std::enable_if_t<CanDoArithmetic<T>::value && CanDoArithmetic<U>::value, int> = 0>
inline constexpr T operator*(T const &A, U const &B) {
	auto C = A;
	C *= B;
	return C;
}

template<typename T, std::enable_if_t<CanDoArithmetic<T>::value, int> = 0>
inline constexpr T operator/(T const &A, T const &B) {
	auto C = A;
	C /= B;
	return C;
}

template<typename T, std::enable_if_t<CanDoArithmetic<T>::value, int> = 0>
inline constexpr T operator%(T const &A, T const &B) {
	auto C = A;
	C %= B;
	return C;
}

template<typename T, std::enable_if_t<CanDoArithmetic<T>::value, int> = 0>
inline constexpr T& operator*=(T &A, typename T::value_type const &b) {
	int const cnt = A.size();
	for (int i = 0; i < cnt; i++) {
		A[i] *= b;
	}
	return A;
}

template<typename T, std::enable_if_t<CanDoArithmetic<T>::value, int> = 0>
inline constexpr T& operator/=(T &A, typename T::value_type const &b) {
	int const cnt = A.size();
	for (int i = 0; i < cnt; i++) {
		A[i] /= b;
	}
	return A;
}

template<typename T, std::enable_if_t<CanDoArithmetic<T>::value, int> = 0>
inline constexpr T operator*(T const &A, typename T::value_type const &b) {
	auto C = A;
	C *= b;
	return C;
}

template<typename T, std::enable_if_t<CanDoArithmetic<T>::value, int> = 0>
inline constexpr T operator*(typename T::value_type const &a, T const &B) {
	return B * a;
}

template<typename T, std::enable_if_t<CanDoArithmetic<T>::value, int> = 0>
inline constexpr T operator/(T const &A, typename T::value_type const &b) {
	auto C = A;
	C /= b;
	return C;
}

template<typename T, std::enable_if_t<CanDoArithmetic<T>::value, int> = 0>
inline constexpr auto dot(T const& a, T const& b) {
	typename T::value_type sum = a[0] * b[0];
	for (int d = 1; d < (int) a.size(); d++) {
		sum += a[d] * b[d];
	}
	return sum;
}


template<typename T, std::enable_if_t<CanDoArithmetic<T>::value, int> = 0>
inline constexpr auto norm(T const& a) {
	return sqrt(dot(a, a));
}

template<typename T, std::enable_if_t<CanDoArithmetic<T>::value, int> = 0>
inline constexpr auto normalize(T const& a) {
	using U = typename ElementType<T>::type;
	return a / (norm(a) + (U(2 * std::numeric_limits<U>::min())));
}
