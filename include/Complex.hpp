/*
 * Complex.hpp
 *
 *  Created on: Dec 14, 2024
 *      Author: dmarce1
 */

#ifndef INCLUDE_COMPLEX_HPP_
#define INCLUDE_COMPLEX_HPP_

namespace Math {

template<typename T>
struct Complex;

template<typename T>
Complex<T> conj(Complex<T>);

template<typename T>
T norm(Complex<T> const&);

template<typename T>
T abs(Complex<T> const&);

template<typename T>
Complex<T> exp(Complex<T> const&);

template<typename T>
struct Complex {
	Complex() :
			x(), y() {
	}
	explicit Complex(T const &_value) :
			x(_value), y(0) {
	}
	constexpr Complex(T const &_x, T const &_y) :
			x(_x), y(_y) {
	}
	std::string toString() const {
		return to_string(x) + " + " + to_string(y) + "I";

	}
	T real() const {
		return x;
	}
	T imaginary() const {
		return y;
	}
	T& real() {
		return x;
	}
	T& imaginary() {
		return y;
	}
	Complex<T> operator+() const {
		return *this;
	}
	Complex<T> operator+(Complex<T> const &B) const {
		auto const &A = *this;
		Complex<T> C;
		C.x = A.x + B.x;
		C.y = A.y + B.y;
		return C;
	}
	Complex<T> operator-() const {
		auto const &A = *this;
		Complex<T> B;
		B.x = -A.x;
		B.y = -A.y;
		return B;
	}
	Complex<T> operator-(Complex const &B) const {
		auto const &A = *this;
		Complex<T> C;
		C.x = A.x - B.x;
		C.y = A.y - B.y;
		return C;
	}
	Complex<T> operator*(Complex<T> const &B) const {
		auto const &A = *this;
		Complex<T> C;
		C.x = A.x * B.x - A.y * B.y;
		C.y = A.x * B.y + A.y * B.x;
		return C;
	}
	Complex<T> operator*(T const &B) const {
		auto C = *this;
		C.x *= B;
		C.y *= B;
		return C;
	}
	Complex<T> operator/(Complex<T> const &B) const {
		T const denInv = T(1) / (B * conj(B)).real();
		auto const num = (*this * conj(B));
		return num * denInv;
	}
	Complex<T>& operator+=(Complex<T> const &B) {
		auto &A = *this;
		A.x += B.x;
		A.y += B.y;
		return *this;
	}
	Complex<T>& operator-=(Complex<T> const &B) {
		auto &A = *this;
		A.x -= B.x;
		A.y -= B.y;
		return *this;
	}
	Complex<T>& operator*=(Complex<T> const &other) {
		*this = *this * other;
		return *this;
	}
	Complex<T>& operator/=(Complex<T> const &other) {
		*this = *this / other;
		return *this;
	}
	friend Complex exp<>(Complex const &z);
	friend Complex conj<>(Complex z);
	friend T norm<>(Complex const &z);
	friend T abs<>(Complex const &z);
private:
	T x;
	T y;
};


template<typename T>
inline Complex<T> exp(Complex<T> const &z) {
	T const theta = z.imaginary();
	T const radius = z.real();
	T const expR = exp(radius);
	T const expCos = cos(theta);
	T const expSin = sin(theta);
	return Complex<T>(expR * expCos, expR * expSin);
}



template<typename T>
inline Complex<T> conj(Complex<T> z) {
	z.y = -z.y;
	return z;
}


template<typename T>
inline T norm(Complex<T> const &z) {
	return z.x * z.x + z.y * z.y;
}


template<typename T>
inline T abs(Complex<T> const &z) {
	auto const result = sqrt(z.x * z.x + z.y * z.y);
	return result;
}

template<typename T>
inline Complex<T> operator*(T const &A, Complex<T> const &B) {
	return B * A;
}

template<typename T>
inline Complex<T> operator+(T const &R, Complex<T> const &C) {
	Complex<T> result = C;
	result += R;
	return result;
}

template<typename T>
inline Complex<T> operator+(Complex<T> const &C, T const &R) {
	return R + C;
}

template<typename T>
inline Complex<T> operator-(T const &R, Complex<T> const &C) {
	return R + (-C);
}

template<typename T>
inline Complex<T> operator-(Complex<T> const &C, T const &R) {
	return -(R - C);
}

template<typename T>
inline Complex<T> operator+=(Complex<T> &C, T const &R) {
	C += Complex<T>(R, T(0));
	return C;
}

template<typename T>
inline Complex<T> operator-=(Complex<T> &C, T const &R) {
	C.x -= R;
	return C;
}

}



#endif /* INCLUDE_COMPLEX_HPP_ */
