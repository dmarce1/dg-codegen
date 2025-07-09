#include "Valarray.hpp"


Slice::Slice(size_t start, size_t size, size_t stride) :
		start_(start), size_(size), stride_(stride) {
}


size_t Slice::start() const {
	return start_;
}

size_t Slice::size() const {
	return size_;
}

size_t Slice::stride() const {
	return stride_;
}

GSlice::GSlice(size_t start, Valarray<size_t> const &sizes, Valarray<size_t> const &strides) :
		start_(start), sizes_(sizes), strides_(strides) {
}

size_t GSlice::start() const {
	return start_;
}

Valarray<size_t> const& GSlice::size() const {
	return sizes_;
}

Valarray<size_t> const& GSlice::stride() const {
	return strides_;
}
