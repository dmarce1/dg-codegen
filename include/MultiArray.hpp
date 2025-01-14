/*
 * MultiArray.hpp
 *
 *  Created on: Jan 13, 2025
 *      Author: dmarce1
 */

#ifndef INCLUDE_MULTIARRAY_HPP_
#define INCLUDE_MULTIARRAY_HPP_

#include <vector>

#include "Interval.hpp"

template<typename Type, int Ndim, int Nfield = 1>
struct MultiArray {
	MultiArray(Interval<size_t, Ndim> const &box_) {
		box = box_;
		size_t const vol = box.volume();
		for( int f = 0; f < Nfield; f++) {
			data[f].resize(vol);
		}
	}
private:
	Interval<size_t, Ndim> box;
	std::array<std::vector<Type>, Nfield> data;
};

#endif /* INCLUDE_MULTIARRAY_HPP_ */
