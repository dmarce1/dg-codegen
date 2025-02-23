/*
 * GridAttributes.hpp
 *
 *  Created on: Feb 21, 2025
 *      Author: dmarce1
 */

#ifndef INCLUDE_GRIDATTRIBUTES_HPP_
#define INCLUDE_GRIDATTRIBUTES_HPP_

#include "Vector.hpp"

#include <valarray>

template<typename T>
struct GridAttributes {
	static constexpr T zero = T(0), one = T(1);
	std::valarray<size_t> intSizes;
	std::valarray<size_t> extSizes;
	std::valarray<size_t> intStrides;
	std::valarray<size_t> extStrides;
	std::valarray<T> gridSpacing;
	T d3Rinv;
	T d3R;
	size_t intSize;
	size_t extSize;
	size_t boundWidth;
	T minSpacing;
	Math::Vector<std::gslice, 2 * NDIM> srcBoundarySlices;
	Math::Vector<std::gslice, 2 * NDIM> dstBoundarySlices;
	GridAttributes(Math::Vector<int, NDIM> const &N, int bw) :
			intSizes(NDIM), extSizes(NDIM), intStrides(NDIM), extStrides(NDIM), gridSpacing(NDIM), boundWidth(bw) {
		intSize = 1;
		extSize = 1;
		for (int k = 0; k < NDIM; k++) {
			intSizes[k] = N[k];
			extSizes[k] = intSizes[k] + 2 * boundWidth;
			extSize *= extSizes[k];
			intSize *= intSizes[k];
		}
		d3Rinv = one;
		d3R = one;
		for (int k = 0; k < NDIM; k++) {
			d3Rinv *= T(intSizes[k]);
			gridSpacing[k] = one / T(intSizes[k]);
			d3R *= gridSpacing[k];
		}
		extStrides[ZDIM] = intStrides[ZDIM] = 1;
		for (int k = NDIM - 1; k > 0; k--) {
			extStrides[k - 1] = extStrides[k] * extSizes[k];
			intStrides[k - 1] = intStrides[k] * intSizes[k];
		}
		minSpacing = gridSpacing.min();
		for (int dim = 0; dim < NDIM; dim++) {
			std::valarray<size_t> sizes(NDIM);
			int const bwStride = boundWidth * extStrides[dim];
			sizes[dim] = boundWidth;
			for (int k = 0; k < dim; k++) {
				sizes[k] = intSizes[k];
			}
			for (int k = dim + 1; k < NDIM; k++) {
				sizes[k] = extSizes[k];
			}
			int srcStart1 = zero;
			for (int k = 0; k <= dim; k++) {
				srcStart1 += boundWidth * extStrides[k];
			}
			size_t const dstStart1 = srcStart1 + intSizes[dim] * extStrides[dim];
			size_t const dstStart2 = srcStart1 - bwStride;
			size_t const srcStart2 = dstStart1 - bwStride;
			size_t const i0 = 2 * dim;
			size_t const i1 = i0 + 1;
			srcBoundarySlices[i0] = std::gslice(srcStart1, sizes, extStrides);
			dstBoundarySlices[i0] = std::gslice(dstStart1, sizes, extStrides);
			srcBoundarySlices[i1] = std::gslice(srcStart2, sizes, extStrides);
			dstBoundarySlices[i1] = std::gslice(dstStart2, sizes, extStrides);
		}
	}
};

#endif /* INCLUDE_GRIDATTRIBUTES_HPP_ */
