/*
 * Metric.hpp
 *
 *  Created on: Mar 6, 2025
 *      Author: dmarce1
 */

#ifndef INCLUDE_METRIC_HPP_
#define INCLUDE_METRIC_HPP_

#include "Matrix.hpp"

enum class tensor_t : int {
	covariant, contravariant, mixed
};

template<typename T, int N>
struct MetricTensor {
	MetricTensor(Math::SymmetricMatrix<T, N> const &metric_) :
			metric(metric_), inverseMetric(matrixInverse(metric)), raiseTransform(inverseMetric), lowerTransform(metric) {
		using namespace Math;
		matrixLUDecompose(raiseTransform);
		matrixLUDecompose(lowerTransform);
	}
	template<int M>
	Math::Matrix<T, N, M> raiseIndex(Math::Matrix<T, N, M> &coVector) const {
		using namespace Math;
		matrixLUMultiply(raiseTransform, coVector);
	}
	template<int M>
	Math::Matrix<T, N, M> lowerIndex(Math::Matrix<T, N, M> &contraVector) const {
		using namespace Math;
		matrixLUMultiply(lowerTransform, contraVector);
	}
	template<tensor_t I = tensor_t::covariant>
	T tensorTrace(Math::SquareMatrix<T, N> &A) const {
		T trace = T(0);
		if constexpr (I == tensor_t::mixed) {
			for (int n = 0; n < N; n++) {
				trace += A[n, n];
			}
		} else {
			auto const G = (I == tensor_t::contravariant) ? metric : inverseMetric;
			for (int n = 0; n < N; n++) {
				for (int k = 0; k < n; k++) {
					trace += G[n, k] * (A[n, k] + A[k, n]);
				}
				trace += G[n, n] * A[n, n];
			}

		}
		return trace;
	}
	template<tensor_t I = tensor_t::covariant>
	T tensorTrace(Math::SymmetricMatrix<T, N> &A) const {
		T trace = T(0);
		if constexpr (I == tensor_t::mixed) {
			for (int n = 0; n < N; n++) {
				trace += A[n, n];
			}
		} else {
			auto const G = (I == tensor_t::contravariant) ? metric : inverseMetric;
			for (int n = 0; n < N; n++) {
				for (int k = 0; k < n; k++) {
					trace += T(2) * G[n, k] * A[n, k];
				}
				trace += G[n, n] * A[n, n];
			}
		}
		return trace;
	}
private:
	Math::SquareMatrix<T, N> metric;
	Math::SquareMatrix<T, N> inverseMetric;
	Math::SquareMatrix<T, N> raiseTransform;
	Math::SquareMatrix<T, N> lowerTransform;
};

#endif /* INCLUDE_METRIC_HPP_ */
