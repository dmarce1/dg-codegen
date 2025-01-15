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
#include "Numbers.hpp"

template<typename, int>
struct MultiArray;

template<typename, int>
struct MultiArrayView;

template<typename Type, int Ndim>
struct MultiArray {
	using iterator = std::vector<Type>::iterator;
	using const_iterator = std::vector<Type>::const_iterator;
	MultiArray();
	MultiArray(int);
	MultiArray(Math::Vector<int, Ndim> const&);
	MultiArray(Interval<int, Ndim> const&);
	MultiArray(MultiArray const&);
	MultiArray(MultiArray&&);
	MultiArray(MultiArrayView<Type, Ndim> const&);
	MultiArray& operator=(MultiArray&&);
	MultiArray& operator=(MultiArray const&);
	MultiArray& operator=(MultiArrayView<Type, Ndim> const&);
	Type& operator[](Math::Vector<int, Ndim> const&);
	Type operator[](Math::Vector<int, Ndim> const&) const;
	MultiArrayView<Type, Ndim> getView(Interval<int, Ndim> const&);
	MultiArrayView<Type, Ndim> getDefaultView();
	const_iterator begin() const;
	const_iterator end() const;
	const_iterator cbegin() const;
	const_iterator cend() const;
	iterator begin();
	iterator end();
	size_t size() const;
	void swapDimensions(int, int);
private:
	int computeIndex(Math::Vector<int, Ndim> const&);
	std::vector<Type> dataVector;
	Interval<int, Ndim> domainBox;
	Math::Vector<int, Ndim> dataStride;
};

template<typename Type, int Ndim>
struct MultiArrayView {
	struct const_iterator;
	struct iterator;
	MultiArrayView(Type*, Interval<int, Ndim> const&, Math::Vector<int, Ndim> const&);
	MultiArrayView();
	MultiArrayView(MultiArrayView const&);
	MultiArrayView(MultiArrayView&&) = delete;
	MultiArrayView& operator=(MultiArrayView const&);
	MultiArrayView& operator=(MultiArrayView&&) = delete;
	MultiArrayView& operator=(MultiArray<Type, Ndim> const&);
	Type& operator[](Math::Vector<int, Ndim> const&);
	Type operator[](Math::Vector<int, Ndim> const&) const;
	iterator begin();
	iterator end();
	const_iterator begin() const;
	const_iterator end() const;
	const_iterator cbegin() const;
	const_iterator cend() const;
	MultiArrayView<Type, Ndim> getView(Interval<int, Ndim> const&);
	int size() const;
	void setGhostWidth(int);
	void swapDimensions(int, int);
	template<typename Derived>
	struct base_iterator {
		using derived_type = Derived;
		bool operator==(Derived const&) const;
		bool operator!=(Derived const&) const;
		derived_type& operator++();
		derived_type& operator--();
		derived_type operator++(int);
		derived_type operator--(int);
		derived_type increment(int, int = 1) const;
		derived_type decrement(int, int = 1) const;
		friend class MultiArrayView;
	protected:
		int scalarIndex;
	private:
		Math::Vector<int, Ndim> vectorIndices;
	};
	struct const_iterator: public base_iterator<const_iterator> {
		static constexpr std::bidirectional_iterator_tag iterator_category { };
		using base_type = base_iterator<const_iterator>;
		using difference_type = int;
		using pointer = Type const*;
		using reference = Type const&;
		using value_type = Type;
		reference operator*();
		pointer operator->();
		friend class MultiArrayView;
		friend class base_iterator<const_iterator> ;
	private:
		MultiArrayView const* getViewPointer() const;
		MultiArrayView const *viewPointer;
	};
	struct iterator: public base_iterator<iterator> {
		static constexpr std::bidirectional_iterator_tag iterator_category { };
		using base_type = base_iterator<iterator>;
		using difference_type = int;
		using pointer = Type*;
		using reference = Type&;
		using value_type = Type;
		reference operator*();
		pointer operator->();
		friend class MultiArrayView;
		friend class base_iterator<const_iterator> ;
	private:
		MultiArrayView const* getViewPointer() const;
		MultiArrayView *viewPointer;
	};
private:
	int computeIndex(Math::Vector<int, Ndim> const&);
	Type *dataPtr;
	Interval<int, Ndim> domainBox;
	Math::Vector<int, Ndim> dataStride;
	int dataStart;
	int ghostWidth;
};

template<typename Type, int Ndim>
MultiArray<Type, Ndim>::MultiArray() :
		dataVector(), domainBox(nullInterval<int, Ndim>()), dataStride(Math::zeroVector<int, Ndim>()) {
}

template<typename Type, int Ndim>
MultiArray<Type, Ndim>::MultiArray(int length) :
		dataVector(Math::integerPower(length, Ndim)) {
	for (int k = 0; k < Ndim; k++) {
		domainBox.begin(k) = 0;
		domainBox.end(k) = length;
	}
	dataStride[Ndim - 1] = 1;
	for (int k = Ndim - 2; k >= 0; k--) {
		dataStride[k] = length * dataStride[k + 1];
	}
}

template<typename Type, int Ndim>
MultiArray<Type, Ndim>::MultiArray(Math::Vector<int, Ndim> const &dims) :
		dataVector(std::accumulate(dims.begin(), dims.end(), 1, std::multiplies<int>())) {
	for (int k = 0; k < Ndim; k++) {
		domainBox.begin(k) = 0;
		domainBox.end(k) = dims[k];
	}
	dataStride[Ndim - 1] = 1;
	for (int k = Ndim - 2; k >= 0; k--) {
		dataStride[k] = dataStride[k + 1] * domainBox.span(k + 1);
	}
}

template<typename Type, int Ndim>
MultiArray<Type, Ndim>::MultiArray(Interval<int, Ndim> const &interval) :
		dataVector(interval.volume()) {
	domainBox = interval;
	dataStride[Ndim - 1] = 1;
	for (int k = Ndim - 2; k >= 0; k--) {
		dataStride[k] = dataStride[k + 1] * domainBox.span(k + 1);
	}
}

template<typename Type, int Ndim>
MultiArray<Type, Ndim>::MultiArray(MultiArray const &other) {
	*this = other;
}

template<typename Type, int Ndim>
MultiArray<Type, Ndim>::MultiArray(MultiArray &&other) {
	*this = std::move(other);
}

template<typename Type, int Ndim>
MultiArray<Type, Ndim>& MultiArray<Type, Ndim>::operator=(MultiArray &&other) {
	assert(size() == other.size());
	dataVector = std::move(other.dataVector);
	domainBox = std::move(other.domainBox);
	dataStride = std::move(other.dataStride);
	return *this;
}

template<typename Type, int Ndim>
MultiArray<Type, Ndim>& MultiArray<Type, Ndim>::operator=(MultiArray const &other) {
	assert(size() == other.size());
	dataVector = other.dataVector;
	domainBox = other.domainBox;
	dataStride = other.dataStride;
	return *this;
}

template<typename Type, int Ndim>
MultiArray<Type, Ndim>::MultiArray(MultiArrayView<Type, Ndim> const &other) {
	*this = other;
}

template<typename Type, int Ndim>
MultiArray<Type, Ndim>& MultiArray<Type, Ndim>::operator=(MultiArrayView<Type, Ndim> const &other) {
	domainBox = other.domainBox;
	dataVector.resize(domainBox.volume());
	dataStride[Ndim - 1] = 1;
	for (int k = Ndim - 2; k >= 0; k--) {
		dataStride[k] = dataStride[k + 1] * domainBox.span(k + 1);
	}
	assert(size() == other.size());
	std::copy(other.begin(), other.end(), begin());
	return *this;
}

template<typename Type, int Ndim>
Type& MultiArray<Type, Ndim>::operator[](Math::Vector<int, Ndim> const &indices) {
	return dataVector[computeIndex(indices)];
}

template<typename Type, int Ndim>
Type MultiArray<Type, Ndim>::operator[](Math::Vector<int, Ndim> const &indices) const {
	return dataVector[computeIndex(indices)];
}

template<typename Type, int Ndim>
MultiArrayView<Type, Ndim> MultiArray<Type, Ndim>::getDefaultView() {
	return MultiArrayView<Type, Ndim>(dataVector.data(), domainBox, dataStride);
}

template<typename Type, int Ndim>
MultiArrayView<Type, Ndim> MultiArray<Type, Ndim>::getView(Interval<int, Ndim> const &subBox) {
	return MultiArrayView<Type, Ndim>(dataVector.data(), subBox, dataStride);
}

template<typename Type, int Ndim>
MultiArray<Type, Ndim>::iterator MultiArray<Type, Ndim>::begin() {
	return dataVector.begin();
}

template<typename Type, int Ndim>
MultiArray<Type, Ndim>::iterator MultiArray<Type, Ndim>::end() {
	return dataVector.end();
}

template<typename Type, int Ndim>
MultiArray<Type, Ndim>::const_iterator MultiArray<Type, Ndim>::begin() const {
	return dataVector.cbegin();
}

template<typename Type, int Ndim>
MultiArray<Type, Ndim>::const_iterator MultiArray<Type, Ndim>::end() const {
	return dataVector.cend();
}

template<typename Type, int Ndim>
MultiArray<Type, Ndim>::const_iterator MultiArray<Type, Ndim>::cbegin() const {
	return dataVector.cbegin();
}

template<typename Type, int Ndim>
MultiArray<Type, Ndim>::const_iterator MultiArray<Type, Ndim>::cend() const {
	return dataVector.cend();
}

template<typename Type, int Ndim>
size_t MultiArray<Type, Ndim>::size() const {
	return dataVector.size();
}

template<typename Type, int Ndim>
void MultiArray<Type, Ndim>::swapDimensions(int k, int n) {
	std::swap(dataStride[k], dataStride[n]);
	std::swap(domainBox[k], domainBox[n]);
}

template<typename Type, int Ndim>
int MultiArray<Type, Ndim>::computeIndex(Math::Vector<int, Ndim> const &indices) {
	int i = 0;
	for (int k = 0; k < Ndim; k++) {
		i += (indices[k] - domainBox.begin(k)) * dataStride[k];
	}
	return i;
}

template<typename Type, int Ndim>
MultiArrayView<Type, Ndim>::MultiArrayView() {
}

template<typename Type, int Ndim>
MultiArrayView<Type, Ndim>::MultiArrayView(Type *dPtr, Interval<int, Ndim> const &subBox_,
		Math::Vector<int, Ndim> const &parentStride) {
	domainBox = subBox_;
	ghostWidth = 0;
	dataStart = 0;
	for (int k = 0; k < Ndim; k++) {
		dataStart += parentStride[k] * domainBox.begin(k);
		dataStride[k] = parentStride[k];
	}
	dataPtr = dPtr;
}

template<typename Type, int Ndim>
MultiArrayView<Type, Ndim>::MultiArrayView(MultiArrayView const &other) :
		dataPtr(other.dataPtr), domainBox(other.domainBox), dataStride(other.dataStride), ghostWidth(other.ghostWidth) {
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
int MultiArrayView<Type, Ndim>::computeIndex(Math::Vector<int, Ndim> const &indices) {
	int i = 0;
	for (int k = 0; k < Ndim; k++) {
		i += (indices[k] - domainBox.begin(k)) * dataStride[k];
	}
	return i;
}

template<typename Type, int Ndim>
MultiArrayView<Type, Ndim>::iterator MultiArrayView<Type, Ndim>::begin() {
	iterator I;
	I.viewPointer = this;
	I.scalarIndex = dataStart;
	for (int k = 0; k < Ndim; k++) {
		I.vectorIndices[k] = ghostWidth;
		I.scalarIndex += ghostWidth * dataStride[k];
	}
	return I;
}

template<typename Type, int Ndim>
MultiArrayView<Type, Ndim>::iterator MultiArrayView<Type, Ndim>::end() {
	iterator I;
	I.viewPointer = this;
	I.vectorIndices[0] = domainBox.end(0) - ghostWidth;
	for (int k = 1; k < Ndim; k++) {
		I.vectorIndices[k] = domainBox.begin(k) + ghostWidth;
	}
	I.scalarIndex = dataStart;
	for (int k = 0; k < Ndim; k++) {
		I.index += I.vectorIndices[k] * dataStride[k];
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
	I.scalarIndex = dataStart;
	for (int k = 0; k < Ndim; k++) {
		I.vectorIndices[k] = ghostWidth;
		I.scalarIndex += ghostWidth * dataStride[k];
	}
	return I;
}

template<typename Type, int Ndim>
MultiArrayView<Type, Ndim>::const_iterator MultiArrayView<Type, Ndim>::cend() const {
	const_iterator I;
	I.viewPointer = this;
	I.vectorIndices[0] = domainBox.end(0) - ghostWidth;
	for (int k = 1; k < Ndim; k++) {
		I.vectorIndices[k] = domainBox.begin(k) + ghostWidth;
	}
	I.scalarIndex = dataStart;
	for (int k = 0; k < Ndim; k++) {
		I.scalarIndex += I.vectorIndices[k] * dataStride[k];
	}
	return I;
}

template<typename Type, int Ndim>
MultiArrayView<Type, Ndim> MultiArrayView<Type, Ndim>::getView(Interval<int, Ndim> const &subBox) {
	return MultiArrayView<Type, Ndim>(dataPtr, subBox, dataStride);
}

template<typename Type, int Ndim>
int MultiArrayView<Type, Ndim>::size() const {
	return domainBox.volume();
}

template<typename Type, int Ndim>
void MultiArrayView<Type, Ndim>::setGhostWidth(int gw) {
	ghostWidth = gw;
}

template<typename Type, int Ndim>
void MultiArrayView<Type, Ndim>::swapDimensions(int k, int n) {
	std::swap(dataStride[k], dataStride[n]);
	domainBox.swapDimensions(k, n);
}

template<typename Type, int Ndim>
template<typename Derived>
bool MultiArrayView<Type, Ndim>::base_iterator<Derived>::operator==(Derived const &other) const {
	return bool((other.scalarIndex == scalarIndex) && (other.getViewPointer() == ((Derived*) this)->getViewPointer()));
}

template<typename Type, int Ndim>
template<typename Derived>
bool MultiArrayView<Type, Ndim>::base_iterator<Derived>::operator!=(Derived const &other) const {
	return !operator==(other);
}

template<typename Type, int Ndim>
template<typename Derived>
Derived& MultiArrayView<Type, Ndim>::base_iterator<Derived>::operator++() {
	auto const *viewPointer = ((Derived*) this)->getViewPointer();
	auto const &box = viewPointer->domainBox;
	auto const &str = viewPointer->dataStride;
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
	auto const *viewPointer = ((Derived*) this)->getViewPointer();
	auto const &box = viewPointer->domainBox;
	auto const &str = viewPointer->dataStride;
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
	I.scalarIndex += cnt * ((Derived*) this)->getViewPointer()->dataStride[k];
	return I;
}

template<typename Type, int Ndim>
template<typename Derived>
Derived MultiArrayView<Type, Ndim>::base_iterator<Derived>::decrement(int k, int cnt) const {
	Derived I = *this;
	I.vectorIndices[k] -= cnt;
	I.scalarIndex -= cnt * ((Derived*) this)->getViewPointer()->dataStride[k];
	return I;
}

template<typename Type, int Ndim>
MultiArrayView<Type, Ndim> const* MultiArrayView<Type, Ndim>::iterator::getViewPointer() const {
	return viewPointer;
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
MultiArrayView<Type, Ndim> const* MultiArrayView<Type, Ndim>::const_iterator::getViewPointer() const {
	return viewPointer;
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
