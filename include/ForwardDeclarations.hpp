/*
 * ForwardDeclarations.hpp
 *
 *  Created on: Jan 14, 2025
 *      Author: dmarce1
 */

#ifndef INCLUDE_FORWARDDECLARATIONS_HPP_
#define INCLUDE_FORWARDDECLARATIONS_HPP_

namespace Math {


template<typename T, typename U>
constexpr T integerPower(T, U);

template<typename T>
constexpr T kroneckerDelta(T, T);

template<typename T>
constexpr T negativeOne2Power(T);

template<typename T>
constexpr T nChooseK(T, T);

template<typename T>
constexpr T nFactorial(T);

template<typename T>
constexpr T nSign(T);

template<typename T>
constexpr T nSquared(T);

}

#endif /* INCLUDE_FORWARDDECLARATIONS_HPP_ */
