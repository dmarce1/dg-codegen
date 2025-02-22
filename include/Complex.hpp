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

struct Complex;

Complex conj(Complex);

Real norm(Complex const&);

Real abs(Complex const&);

Complex exp(Complex const&);

struct Complex {
	using Type = Real;
	Complex() :
			x(), y() {
	}
	explicit Complex(Type const &_value) :
			x(_value), y(0) {
	}
	constexpr Complex(Type const &_x, Type const &_y) :
			x(_x), y(_y) {
	}
	std::string toString() const {
		return to_string(x) + " + " + to_string(y) + "I";

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
	friend Complex exp(Complex const &z);
	friend Complex conj(Complex z);
	friend Type norm(Complex const &z);
	friend Type abs(Complex const &z);
private:
	Type x;
	Type y;
};


inline Complex exp(Complex const &z) {
	Real const theta = z.imaginary();
	Real const radius = z.real();
	Real const expR = exp(radius);
	Real const expCos = cos(theta);
	Real const expSin = sin(theta);
	return Complex(expR * expCos, expR * expSin);
}



inline Complex conj(Complex z) {
	z.y = -z.y;
	return z;
}


inline Real norm(Complex const &z) {
	return z.x * z.x + z.y * z.y;
}


inline Real abs(Complex const &z) {
	auto const result = sqrt(z.x * z.x + z.y * z.y);
	return result;
}

template<typename Type>
inline Complex operator*(Type const &A, Complex const &B) {
	return B * A;
}

template<typename Type>
inline Complex operator+(Type const &R, Complex const &C) {
	Complex result = C;
	result += R;
	return result;
}

template<typename Type>
inline Complex operator+(Complex const &C, Type const &R) {
	return R + C;
}

template<typename Type>
inline Complex operator-(Type const &R, Complex const &C) {
	return R + (-C);
}

template<typename Type>
inline Complex operator-(Complex const &C, Type const &R) {
	return -(R - C);
}

template<typename Type>
inline Complex operator+=(Complex &C, Type const &R) {
	C += Complex(R, Real(0));
	return C;
}

template<typename Type>
inline Complex operator-=(Complex &C, Type const &R) {
	C.x -= R;
	return C;
}

}

constexpr Math::Complex operator""_Imag(long double a) {
	return Math::Complex(0.0_Real, Math::Real(a));
}


#endif /* INCLUDE_COMPLEX_HPP_ */
