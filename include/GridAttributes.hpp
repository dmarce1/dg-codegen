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
	std::valarray<size_t> gridStrides;
	std::valarray<T> gridSpacing;
	size_t intSize;
	size_t extSize;
	size_t boundWidth;
	T minSpacing;
	Math::Vector<std::gslice, 2 * NDIM> srcBoundarySlices;
	Math::Vector<std::gslice, 2 * NDIM> dstBoundarySlices;
	GridAttributes(Math::Vector<int, NDIM> const &N, int bw) :
			intSizes(NDIM), extSizes(NDIM), gridStrides(NDIM), gridSpacing(NDIM), boundWidth(bw) {
		intSize = 1;
		extSize = 1;
		for (int k = 0; k < NDIM; k++) {
			intSizes[k] = N[k];
			extSizes[k] = intSizes[k] + 2 * boundWidth;
			extSize *= extSizes[k];
			intSize *= intSizes[k];
		}
		for (int k = 0; k < NDIM; k++) {
			gridSpacing[k] = one / T(intSizes[k]);
		}
		gridStrides[XDIM] = extSizes[YDIM] * extSizes[ZDIM];
		gridStrides[YDIM] = extSizes[ZDIM];
		gridStrides[ZDIM] = 1;
		minSpacing = gridSpacing.min();
		for (int dim = 0; dim < NDIM; dim++) {
			std::valarray<size_t> sizes(NDIM);
			int const bwStride = boundWidth * gridStrides[dim];
			sizes[dim] = boundWidth;
			for (int k = 0; k < dim; k++) {
				sizes[k] = intSizes[k];
			}
			for (int k = dim + 1; k < NDIM; k++) {
				sizes[k] = extSizes[k];
			}
			int srcStart1 = zero;
			for (int k = 0; k <= dim; k++) {
				srcStart1 += boundWidth * gridStrides[k];
			}
			size_t const dstStart1 = srcStart1 + intSizes[dim] * gridStrides[dim];
			size_t const dstStart2 = srcStart1 - bwStride;
			size_t const srcStart2 = dstStart1 - bwStride;
			size_t const i0 = 2 * dim;
			size_t const i1 = i0 + 1;
			srcBoundarySlices[i0] = std::gslice(srcStart1, sizes, gridStrides);
			dstBoundarySlices[i0] = std::gslice(dstStart1, sizes, gridStrides);
			srcBoundarySlices[i1] = std::gslice(srcStart2, sizes, gridStrides);
			dstBoundarySlices[i1] = std::gslice(dstStart2, sizes, gridStrides);
		}
	}
};




#endif /* INCLUDE_GRIDATTRIBUTES_HPP_ */
