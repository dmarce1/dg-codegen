#pragma once

#include <type_traits>
#include <utility>

template<typename T, typename = void>
struct isArithmeticable {
	static constexpr bool value = false;
};

template<typename T>
struct isArithmeticable<T, std::void_t<decltype(std::declval<T&>()[std::declval<std::size_t>()]), decltype(std::declval<T&>().size())>> {
	static constexpr bool value = std::is_integral<decltype(std::declval<T&>().size())>::value;
};

template<typename T, std::enable_if_t<isArithmeticable<T>::value, int> = 0>
inline constexpr T& operator+=(T &A, T const &B) {
	int const cnt = A.size();
	for (int i = 0; i < cnt; i++) {
		A[i] += B[i];
	}
	return A;
}

template<typename T, std::enable_if_t<isArithmeticable<T>::value, int> = 0>
inline constexpr T& operator-=(T &A, T const &B) {
	int const cnt = A.size();
	for (int i = 0; i < cnt; i++) {
		A[i] -= B[i];
	}
	return A;
}

template<typename T, std::enable_if_t<isArithmeticable<T>::value, int> = 0>
inline constexpr T& operator*=(T &A, T const &B) {
	int const cnt = A.size();
	for (int i = 0; i < cnt; i++) {
		A[i] *= B[i];
	}
	return A;
}

template<typename T, std::enable_if_t<isArithmeticable<T>::value, int> = 0>
inline constexpr T& operator/=(T &A, T const &B) {
	int const cnt = A.size();
	for (int i = 0; i < cnt; i++) {
		A[i] /= B[i];
	}
	return A;
}

template<typename T, std::enable_if_t<isArithmeticable<T>::value, int> = 0>
inline constexpr T& operator%=(T &A, T const &B) {
	int const cnt = A.size();
	for (int i = 0; i < cnt; i++) {
		A[i] %= B[i];
	}
	return A;
}

template<typename T, std::enable_if_t<isArithmeticable<T>::value, int> = 0>
inline constexpr T operator+(T const &A) {
	return A;
}

template<typename T, std::enable_if_t<isArithmeticable<T>::value, int> = 0>
inline constexpr T operator-(T const &A) {
	auto B = A;
	for (auto &b : B) {
		b = -b;
	}
	return B;
}

template<typename T, std::enable_if_t<isArithmeticable<T>::value, int> = 0>
inline constexpr T operator+(T const &A, T const &B) {
	auto C = A;
	C += B;
	return C;
}

template<typename T, std::enable_if_t<isArithmeticable<T>::value, int> = 0>
inline constexpr T operator-(T const &A, T const &B) {
	auto C = A;
	C -= B;
	return C;
}

template<typename T, std::enable_if_t<isArithmeticable<T>::value, int> = 0>
inline constexpr T operator*(T const &A, T const &B) {
	auto C = A;
	C *= B;
	return C;
}

template<typename T, std::enable_if_t<isArithmeticable<T>::value, int> = 0>
inline constexpr T operator/(T const &A, T const &B) {
	auto C = A;
	C /= B;
	return C;
}

template<typename T, std::enable_if_t<isArithmeticable<T>::value, int> = 0>
inline constexpr T operator%(T const &A, T const &B) {
	auto C = A;
	C %= B;
	return C;
}

template<typename T, std::enable_if_t<isArithmeticable<T>::value, int> = 0>
inline constexpr T& operator*=(T &A, typename T::value_type const &b) {
	int const cnt = A.size();
	for (int i = 0; i < cnt; i++) {
		A[i] *= b;
	}
	return A;
}

template<typename T, std::enable_if_t<isArithmeticable<T>::value, int> = 0>
inline constexpr T& operator/=(T &A, typename T::value_type const &b) {
	int const cnt = A.size();
	for (int i = 0; i < cnt; i++) {
		A[i] /= b;
	}
	return A;
}

template<typename T, std::enable_if_t<isArithmeticable<T>::value, int> = 0>
inline constexpr T operator*(T const &A, typename T::value_type const &b) {
	auto C = A;
	C *= b;
	return C;
}

template<typename T, std::enable_if_t<isArithmeticable<T>::value, int> = 0>
inline constexpr T operator*(typename T::value_type const &a, T const &B) {
	return B * a;
}

template<typename T, std::enable_if_t<isArithmeticable<T>::value, int> = 0>
inline constexpr T operator/(T const &A, typename T::value_type const &b) {
	auto C = A;
	C /= b;
	return C;
}

