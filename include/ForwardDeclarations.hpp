/*
 * ForwardDeclarations.hpp
 *
 *  Created on: Jan 14, 2025
 *      Author: dmarce1
 */

#ifndef INCLUDE_FORWARDDECLARATIONS_HPP_
#define INCLUDE_FORWARDDECLARATIONS_HPP_

namespace Math {


template<typename T>
constexpr T integerPower(T, int);

template<typename T>
constexpr T kroneckerDelta(int, int);

template<typename T>
constexpr T negativeOne2Power(int);

template<typename T>
constexpr T nChooseK(int, int);

template<typename T>
constexpr T nFactorial(int);

template<typename T>
constexpr T nSign(T);

template<typename T>
constexpr T nSquared(T);

}

#endif /* INCLUDE_FORWARDDECLARATIONS_HPP_ */
