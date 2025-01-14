/*
 * MultiArray.hpp
 *
 *  Created on: Jan 13, 2025
 *      Author: dmarce1
 */

#ifndef INCLUDE_MULTIARRAY_HPP_
#define INCLUDE_MULTIARRAY_HPP_

#include <iterator>
#include <vector>

#include "Interval.hpp"

template<typename, int>
struct MultiArray;

template<typename, int>
struct MultiArrayView;

template<typename Type, int Ndim>
struct MultiArray {
	using iterator = std::vector<Type>::iterator;
	using const_iterator = std::vector<Type>::const_iterator;
	MultiArray(Math::Vector<int, Ndim> const&);
	MultiArray& operator=(MultiArray const&);
	MultiArray& operator=(MultiArrayView<Type, Ndim> const&);
	MultiArrayView<Type, Ndim> getView(Interval<int, Ndim> const&);
	MultiArrayView<Type, Ndim> getDefaultView();
	const_iterator begin() const;
	const_iterator end() const;
	const_iterator cbegin() const;
	const_iterator cend() const;
	iterator begin();
	iterator end();
	size_t size() const;
private:
	std::vector<Type> data;
	Math::Vector<int, Ndim> dimension;
	Math::Vector<int, Ndim> stride;
};

template<typename Type, int Ndim>
struct MultiArrayView {
	template<typename >
	struct base_iterator;
	struct iterator;
	struct const_iterator;
	MultiArrayView(Type*, Interval<int, Ndim> const&, Math::Vector<int, Ndim> const&);
	MultiArrayView();
	MultiArrayView(MultiArrayView const&);
	MultiArrayView(MultiArrayView&&) = delete;
	MultiArrayView& operator=(MultiArrayView const&);
	MultiArrayView& operator=(MultiArray<Type, Ndim> const&);
	MultiArrayView& operator=(MultiArrayView&&) = delete;
	iterator begin();
	iterator end();
	const_iterator begin() const;
	const_iterator end() const;
	const_iterator cbegin() const;
	const_iterator cend() const;
	MultiArrayView<Type, Ndim> getView(Interval<int, Ndim> const&);
	int size() const;
	void setGhostWidth(int);
	template<typename Derived>
	struct base_iterator {
		friend class MultiArrayView;
		using derived_type = Derived;
		bool operator==(Derived const&) const;
		bool operator!=(Derived const&) const;
		derived_type& operator++();
		derived_type& operator--();
		derived_type operator++(int);
		derived_type operator--(int);
		derived_type increment(int, int = 1) const;
		derived_type decrement(int, int = 1) const;
	protected:
		int scalarIndex;
	private:
		Math::Vector<int, Ndim> vectorIndices;
		MultiArrayView const *viewPointer;
	};
	struct iterator: public base_iterator<iterator> {
		friend class MultiArrayView;
		static constexpr std::bidirectional_iterator_tag iterator_category { };
		using base_type = base_iterator<iterator>;
		using difference_type = int;
		using pointer = Type*;
		using reference = Type&;
		using value_type = Type;
		reference operator*();
		pointer operator->();
	private:
		MultiArrayView *viewPointer;
	};
	struct const_iterator: public base_iterator<const_iterator> {
		friend class MultiArrayView;
		static constexpr std::bidirectional_iterator_tag iterator_category { };
		using base_type = base_iterator<const_iterator>;
		using difference_type = int;
		using pointer = Type const*;
		using reference = Type const&;
		using value_type = Type;
		reference operator*();
		pointer operator->();
	private:
		MultiArrayView const *viewPointer;
	};
private:
	Type *dataPtr;
	Interval<int, Ndim> subBox;
	Math::Vector<int, Ndim> stride;
	int start;
	int ghostWidth;
};

template<typename Type, int Ndim>
MultiArray<Type, Ndim>::MultiArray(Math::Vector<int, Ndim> const &dims) :
		data(std::accumulate(dims.begin(), dims.end(), 1, std::multiplies<int>())) {
	stride[Ndim - 1] = 1;
	for (int k = 0; k < Ndim; k++) {
		dimension[k] = dims[k];
	}
	for (int k = Ndim - 1; k > 0; k--) {
		stride[k - 1] = stride[k] * dims[k];
	}
}

template<typename Type, int Ndim>
MultiArray<Type, Ndim>& MultiArray<Type, Ndim>::operator=(MultiArrayView<Type, Ndim> const &other) {
	assert(size() == other.size());
	std::copy(other.begin(), other.end(), begin());
	return *this;
}

template<typename Type, int Ndim>
MultiArray<Type, Ndim>& MultiArray<Type, Ndim>::operator=(MultiArray const &other) {
	assert(size() == other.size());
	std::copy(other.begin(), other.end(), begin());
	return *this;
}

template<typename Type, int Ndim>
MultiArrayView<Type, Ndim> MultiArray<Type, Ndim>::getDefaultView() {
	return MultiArrayView<Type, Ndim>(data.data(), dimension, stride);
}

template<typename Type, int Ndim>
MultiArrayView<Type, Ndim> MultiArray<Type, Ndim>::getView(Interval<int, Ndim> const& subBox) {
	return MultiArrayView<Type, Ndim>(data.data(), subBox, stride);
}

template<typename Type, int Ndim>
MultiArray<Type, Ndim>::iterator MultiArray<Type, Ndim>::begin() {
	return data.begin();
}

template<typename Type, int Ndim>
MultiArray<Type, Ndim>::iterator MultiArray<Type, Ndim>::end() {
	return data.end();
}

template<typename Type, int Ndim>
MultiArray<Type, Ndim>::const_iterator MultiArray<Type, Ndim>::begin() const {
	return data.cbegin();
}

template<typename Type, int Ndim>
MultiArray<Type, Ndim>::const_iterator MultiArray<Type, Ndim>::end() const {
	return data.cend();
}

template<typename Type, int Ndim>
MultiArray<Type, Ndim>::const_iterator MultiArray<Type, Ndim>::cbegin() const {
	return data.cbegin();
}

template<typename Type, int Ndim>
MultiArray<Type, Ndim>::const_iterator MultiArray<Type, Ndim>::cend() const {
	return data.cend();
}

template<typename Type, int Ndim>
size_t MultiArray<Type, Ndim>::size() const {
	return data.size();
}

template<typename Type, int Ndim>
MultiArrayView<Type, Ndim>::MultiArrayView() {
}

template<typename Type, int Ndim>
MultiArrayView<Type, Ndim>::MultiArrayView(Type *dPtr, Interval<int, Ndim> const &subBox_,
		Math::Vector<int, Ndim> const &parentStride) {
	subBox = subBox_;
	ghostWidth = 0;
	start = 0;
	for (int k = 0; k < Ndim; k++) {
		start += parentStride[k] * subBox.begin(k);
		stride[k] = parentStride[k];
	}
	dataPtr = dPtr;
}

template<typename Type, int Ndim>
MultiArrayView<Type, Ndim>::MultiArrayView(MultiArrayView const &other) :
		dataPtr(other.dataPtr), subBox(other.subBox), stride(other.stride), ghostWidth(other.ghostWidth) {
}

template<typename Type, int Ndim>
MultiArrayView<Type, Ndim>& MultiArrayView<Type, Ndim>::operator=(MultiArrayView const &other) {
	assert(size() == other.size());
	std::copy(other.begin(), other.end(), begin());
	return *this;
}

template<typename Type, int Ndim>
MultiArrayView<Type, Ndim>& MultiArrayView<Type, Ndim>::operator=(MultiArray<Type, Ndim> const &other) {
	assert(size() == other.size());
	std::copy(other.begin(), other.end(), begin());
	return *this;
}

template<typename Type, int Ndim>
MultiArrayView<Type, Ndim>::iterator MultiArrayView<Type, Ndim>::begin() {
	iterator I;
	I.viewPointer = this;
	I.scalarIndex = start;
	for (int k = 0; k < Ndim; k++) {
		I.vectorIndices[k] = ghostWidth;
		I.scalarIndex += ghostWidth * stride[k];
	}
	return I;
}

template<typename Type, int Ndim>
MultiArrayView<Type, Ndim>::iterator MultiArrayView<Type, Ndim>::end() {
	iterator I;
	I.viewPointer = this;
	I.vectorIndices[0] = subBox.end(0) - ghostWidth;
	for (int k = 1; k < Ndim; k++) {
		I.vectorIndices[k] = subBox.begin(k) + ghostWidth;
	}
	I.scalarIndex = start;
	for (int k = 0; k < Ndim; k++) {
		I.index += I.vectorIndices[k] * stride[k];
	}
	return I;
}

template<typename Type, int Ndim>
MultiArrayView<Type, Ndim>::const_iterator MultiArrayView<Type, Ndim>::begin() const {
	return cbegin();
}

template<typename Type, int Ndim>
MultiArrayView<Type, Ndim>::const_iterator MultiArrayView<Type, Ndim>::end() const {
	return cend();
}

template<typename Type, int Ndim>
MultiArrayView<Type, Ndim>::const_iterator MultiArrayView<Type, Ndim>::cbegin() const {
	const_iterator I;
	I.viewPointer = this;
	I.scalarIndex = start;
	for (int k = 0; k < Ndim; k++) {
		I.vectorIndices[k] = ghostWidth;
		I.scalarIndex += ghostWidth * stride[k];
	}
	return I;
}

template<typename Type, int Ndim>
MultiArrayView<Type, Ndim>::const_iterator MultiArrayView<Type, Ndim>::cend() const {
	const_iterator I;
	I.viewPointer = this;
	I.vectorIndices[0] = subBox.end(0) - ghostWidth;
	for (int k = 1; k < Ndim; k++) {
		I.vectorIndices[k] = subBox.begin(k) + ghostWidth;
	}
	I.scalarIndex = start;
	for (int k = 0; k < Ndim; k++) {
		I.scalarIndex += I.vectorIndices[k] * stride[k];
	}
	return I;
}

template<typename Type, int Ndim>
MultiArrayView<Type, Ndim> MultiArrayView<Type, Ndim>::getView(Interval<int, Ndim> const& subBox) {
	return MultiArrayView<Type, Ndim>(dataPtr, subBox, stride);
}

template<typename Type, int Ndim>
int MultiArrayView<Type, Ndim>::size() const {
	return subBox.volume();
}

template<typename Type, int Ndim>
void MultiArrayView<Type, Ndim>::setGhostWidth(int gw) {
	ghostWidth = gw;
}

template<typename Type, int Ndim>
template<typename Derived>
bool MultiArrayView<Type, Ndim>::base_iterator<Derived>::operator==(Derived const &other) const {
	return bool((other.scalarIndex == scalarIndex) && (other.viewPointer == viewPointer));
}

template<typename Type, int Ndim>
template<typename Derived>
bool MultiArrayView<Type, Ndim>::base_iterator<Derived>::operator!=(Derived const &other) const {
	return !operator==(other);
}

template<typename Type, int Ndim>
template<typename Derived>
Derived& MultiArrayView<Type, Ndim>::base_iterator<Derived>::operator++() {
	auto const &box = viewPointer->subBox;
	auto const &str = viewPointer->stride;
	auto const &gw = viewPointer->ghostWidth;
	int dim = Ndim - 1;
	while (vectorIndices[dim] == box.end(dim) - gw - 1) {
		vectorIndices[dim] = gw;
		scalarIndex -= str[dim] * (box.span(dim) - 2 * gw - 1);
		dim--;
		if (dim == 0) {
			break;
		}
	}
	vectorIndices[dim]++;
	scalarIndex += str[dim];
	return *(Derived*) this;
}

template<typename Type, int Ndim>
template<typename Derived>
Derived& MultiArrayView<Type, Ndim>::base_iterator<Derived>::operator--() {
	auto const &box = viewPointer->subBox;
	auto const &str = viewPointer->stride;
	auto const &gw = viewPointer->ghostWidth;
	int dim = 0;
	while (vectorIndices[dim] == box.begin(dim) + gw) {
		vectorIndices[dim] = box.end(dim) - 1 - gw;
		scalarIndex += str[dim] * (box.span(dim) - 2 * gw - 1);
		dim++;
		if (dim == Ndim - 1) {
			break;
		}
	}
	vectorIndices[dim]++;
	scalarIndex += str[dim];
	return *(Derived*) this;
}

template<typename Type, int Ndim>
template<typename Derived>
Derived MultiArrayView<Type, Ndim>::base_iterator<Derived>::operator++(int) {
	Derived rc = *this;
	operator++();
	return rc;
}

template<typename Type, int Ndim>
template<typename Derived>
Derived MultiArrayView<Type, Ndim>::base_iterator<Derived>::operator--(int) {
	Derived rc = *this;
	operator--();
	return rc;
}

template<typename Type, int Ndim>
template<typename Derived>
Derived MultiArrayView<Type, Ndim>::base_iterator<Derived>::increment(int k, int cnt) const {
	Derived I = *this;
	I.vectorIndices[k] += cnt;
	I.scalarIndex += cnt * viewPointer->stride[k];
	return I;
}

template<typename Type, int Ndim>
template<typename Derived>
Derived MultiArrayView<Type, Ndim>::base_iterator<Derived>::decrement(int k, int cnt) const {
	Derived I = *this;
	I.vectorIndices[k] -= cnt;
	I.scalarIndex -= cnt * viewPointer->stride[k];
	return I;
}

template<typename Type, int Ndim>
MultiArrayView<Type, Ndim>::iterator::reference MultiArrayView<Type, Ndim>::iterator::operator*() {
	return viewPointer->dataPtr[base_type::scalarIndex];
}

template<typename Type, int Ndim>
MultiArrayView<Type, Ndim>::iterator::pointer MultiArrayView<Type, Ndim>::iterator::operator->() {
	return &viewPointer->dataPtr[base_type::scalarIndex];
}

template<typename Type, int Ndim>
MultiArrayView<Type, Ndim>::const_iterator::reference MultiArrayView<Type, Ndim>::const_iterator::operator*() {
	return viewPointer->dataPtr[base_type::scalarIndex];
}

template<typename Type, int Ndim>
MultiArrayView<Type, Ndim>::const_iterator::pointer MultiArrayView<Type, Ndim>::const_iterator::operator->() {
	return &viewPointer->dataPtr[base_type::scalarIndex];
}

#endif /* INCLUDE_MULTIARRAY_HPP_ */
