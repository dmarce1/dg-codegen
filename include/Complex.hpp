/*
 * Complex.hpp
 *
 *  Created on: Dec 14, 2024
 *      Author: dmarce1
 */

#ifndef INCLUDE_COMPLEX_HPP_
#define INCLUDE_COMPLEX_HPP_

#include "Real.hpp"

namespace Math {

template<typename >
struct Complex;

template<typename T>
Complex<T> conj(Complex<T>);

template<typename T>
T norm(Complex<T> const&);

template<typename T>
T abs(Complex<T> const&);

template<typename T>
Complex<T> exp(Complex<T> const&);

template<typename Type>
struct Complex {
	constexpr Complex() {
		if constexpr (CHECK_REALS) {
			y = x = std::numeric_limits<Type>::signaling_NaN();

		}
	}
	constexpr Complex(Type const &_value) :
			x(_value), y(0) {
	}
	constexpr Complex(Type const &_x, Type const &_y) :
			x(_x), y(_y) {
	}
	std::string toString() const {
		return std::to_string(x) + " + " + std::to_string(y) + "I";

	}
	Type real() const {
		return x;
	}
	Type imaginary() const {
		return y;
	}
	Type& real() {
		return x;
	}
	Type& imaginary() {
		return y;
	}
	Complex operator+() const {
		return *this;
	}
	Complex operator+(Complex const &B) const {
		auto const &A = *this;
		Complex C;
		C.x = A.x + B.x;
		C.y = A.y + B.y;
		return C;
	}
	Complex operator-() const {
		auto const &A = *this;
		Complex B;
		B.x = -A.x;
		B.y = -A.y;
		return B;
	}
	Complex operator-(Complex const &B) const {
		auto const &A = *this;
		Complex C;
		C.x = A.x - B.x;
		C.y = A.y - B.y;
		return C;
	}
	Complex operator*(Complex const &B) const {
		auto const &A = *this;
		Complex C;
		C.x = A.x * B.x - A.y * B.y;
		C.y = A.x * B.y + A.y * B.x;
		return C;
	}
	Complex operator*(Type const &B) const {
		auto C = *this;
		C.x *= B;
		C.y *= B;
		return C;
	}
	Complex operator/(Complex const &B) const {
		Type const denInv = Type(1) / (B * conj(B)).real();
		auto const num = (*this * conj(B));
		return num * denInv;
	}
	Complex& operator+=(Complex const &B) {
		auto &A = *this;
		A.x += B.x;
		A.y += B.y;
		return *this;
	}
	Complex& operator-=(Complex const &B) {
		auto &A = *this;
		A.x -= B.x;
		A.y -= B.y;
		return *this;
	}
	Complex& operator*=(Complex const &other) {
		*this = *this * other;
		return *this;
	}
	Complex& operator/=(Complex const &other) {
		*this = *this / other;
		return *this;
	}
	friend Complex<Type> exp<Type>(Complex<Type> const &z);
	friend Complex<Type> conj<Type>(Complex<Type> z);
	friend Type norm<Type>(Complex<Type> const &z);
	friend Type abs<Type>(Complex<Type> const &z);
private:
	Type x;
	Type y;
};

template<typename T>
Complex<T> exp(Complex<T> const &z) {
	T const theta = T(2.0 * M_PI) * z.imaginary();
	T const radius = z.real();
	T const expR = std::exp(radius);
	T const expCos = std::cos(theta);
	T const expSin = std::sin(theta);
	return Complex<T>(expR * expCos, expR * expSin);
}

template<typename T>
Complex<T> conj(Complex<T> z) {
	z.y = -z.y;
	return z;
}

template<typename T>
T norm(Complex<T> const &z) {
	return z.x * z.x + z.y * z.y;
}

template<typename T>
T abs(Complex<T> const &z) {
	auto const result = sqrt(z.x * z.x + z.y * z.y);
	return result;
}

template<typename Type>
Complex<Type> operator*(Type const &A, Complex<Type> const &B) {
	return B * A;
}

template<typename Type>
Complex<Type> operator+(Type const &R, Complex<Type> const &C) {
	Complex<Type> result = C;
	result += R;
	return result;
}

template<typename Type>
Complex<Type> operator+(Complex<Type> const &C, Type const &R) {
	return R + C;
}

template<typename Type>
Complex<Type> operator-(Type const &R, Complex<Type> const &C) {
	return R + (-C);
}

template<typename Type>
Complex<Type> operator-(Complex<Type> const &C, Type const &R) {
	return -(R - C);
}

template<typename Type>
Complex<Type> operator+=(Complex<Type> &C, Type const &R) {
	C += Complex<Type>(R, 0);
	return C;
}

template<typename Type>
Complex<Type> operator-=(Complex<Type> &C, Type const &R) {
	C.x -= R;
	return C;
}

}

constexpr Math::Complex<double> operator""_Imag(long double a) {
	return Math::Complex<double>(0, a);
}

#endif /* INCLUDE_COMPLEX_HPP_ */
