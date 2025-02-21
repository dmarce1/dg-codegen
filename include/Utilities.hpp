/*
 * Utilities.hpp
 *
 *  Created on: Dec 12, 2024
 *      Author: dmarce1
 */

#ifndef INCLUDE_UTILITIES_HPP_
#define INCLUDE_UTILITIES_HPP_

#include <cctype>

#define BINARY_OPERATOR_MEMBER(CLASS, OPERATOR)      \
	CLASS operator OPERATOR (CLASS const &B) const { \
        CLASS C;                                     \
		CLASS const &A = *this;                      \
		for (int n = 0; n != int(size()); n++) {     \
			C[n] = A[n] OPERATOR B[n];               \
		}                                            \
		return C;                                    \
	}                                                \
	CLASS& operator OPERATOR##=(CLASS const &B) {    \
    	CLASS &A = *this;                            \
		for (int n = 0; n != int(size()); n++) {     \
			A[n] OPERATOR##= B[n];                   \
		}                                            \
		return *this;                                \
	}

#define BINARY_OPERATOR_MEMBER_SCALAR(CLASS, SCALAR, OPERATOR)               \
	CLASS operator OPERATOR (SCALAR const &B) const {                        \
        CLASS C;                                                             \
		CLASS const &A = *this;                                              \
		for (int n = 0; n != int(size()); n++) {                             \
			C[n] = A[n] OPERATOR B;                                          \
		}                                                                    \
		return C;                                                            \
	}                                                                        \
	friend CLASS operator OPERATOR (SCALAR const &A, CLASS const &B) {       \
        CLASS C;                                                             \
		for (int n = 0; n != int(size()); n++) {                             \
			C[n] = A OPERATOR B[n];                                          \
		}                                                                    \
		return C;                                                            \
	}                                                                        \
	CLASS operator OPERATOR##=(SCALAR const &B) {                            \
    	CLASS &A = *this;                                                    \
		for (int n = 0; n < int(size()); n++) {                              \
			A[n] OPERATOR##= B;                                              \
		}                                                                    \
		return *this;                                                        \
	}

#define UNARY_OPERATOR_MEMBER(CLASS, OPERATOR)            \
	CLASS operator OPERATOR () const {                    \
        CLASS C;                                          \
		CLASS const &A = *this;                           \
		for (int n = 0; n < int(size()); n++) {           \
			C[n] = OPERATOR A[n];                         \
		}                                                 \
		return C;                                         \
	}
#define USE_STANDARD_ARITHMETIC(CLASS, TYPE)      \
	BINARY_OPERATOR_MEMBER(CLASS, +)              \
	BINARY_OPERATOR_MEMBER(CLASS, -)              \
	UNARY_OPERATOR_MEMBER(CLASS, +)               \
	UNARY_OPERATOR_MEMBER(CLASS, -)               \
    BINARY_OPERATOR_MEMBER_SCALAR(CLASS, TYPE, *) \
    BINARY_OPERATOR_MEMBER_SCALAR(CLASS, TYPE, /)

#define USE_STANDARD_DEFAULTS(CLASS)          \
	constexpr CLASS() = default;              \
	CLASS(CLASS const&) = default;            \
	CLASS(CLASS&&) = default;                 \
	CLASS& operator=(CLASS const&) = default; \
	CLASS& operator=(CLASS&&) = default;

#define STANDARD_COMPARISON_OPERATORS(CLASS) \
	auto operator !=(CLASS const& other) const { \
		return (other < *this) || (*this < other); \
	} \
	auto operator ==(CLASS const& other) const { \
		return !(other != *this); \
	} \
	auto operator >(CLASS const& other) const { \
		return other < *this; \
	} \
	auto operator >=(CLASS const& other) const { \
		return !(*this < other); \
	} \
	auto operator <=(CLASS const& other) const { \
		return !(*this > other); \
	} \

template<auto Init>
struct AutoInit {
	using value_type = decltype(Init);
	constexpr AutoInit() :
			value(Init) {
	}
	AutoInit(value_type _value) :
			value(_value) {
	}
	AutoInit(AutoInit const &init) :
			value(init) {
	}
	AutoInit(AutoInit &&init) :
			value(std::move(init)) {
	}
	AutoInit& operator=(AutoInit const &_value) {
		value = _value;
	}
	AutoInit& operator=(AutoInit &&_value) {
		value = std::move(_value);
	}
	operator const value_type() const {
		return value;
	}
	operator value_type&() {
		return value;
	}
private:
	value_type value;
};

#include <deque>
#include <valarray>
#include <vector>

template<typename, int>
struct containerHasResize: public std::false_type {
	static constexpr bool value = false;
};

template<typename T, int N>
struct containerHasResize<std::deque<T>, N> : public std::true_type {
	static constexpr bool value = true;
};

template<typename T, int N>
struct containerHasResize<std::valarray<T>, N> : public std::true_type {
	static constexpr bool value = true;
};

template<typename T, int N>
struct containerHasResize<std::vector<T>, N> : public std::true_type {
	static constexpr bool value = true;
};

template<typename Container, int N, bool = containerHasResize<Container, N>::value>
struct ContainerResizer;

template<typename Container, int N>
struct ContainerResizer<Container, N, false> {
	ContainerResizer(Container&) {
	}
};

template<typename Container, int N>
struct ContainerResizer<Container, N, true> {
	ContainerResizer(Container &container) {
		container = Container(N);
	}
};

template <typename T, typename = void>
struct Underlying {
    using type = T;
};

template <typename T>
struct Underlying<T, std::void_t<typename T::value_type>> {
    using type = typename T::value_type;
};

template <typename T>
using underlying_type = typename Underlying<T>::type;

#endif /* INCLUDE_UTILITIES_HPP_ */
