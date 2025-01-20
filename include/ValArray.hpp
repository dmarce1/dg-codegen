/*
 * ValArray.hpp
 *
 *  Created on: Jan 19, 2025
 *      Author: dmarce1
 */

#ifndef INCLUDE_VALARRAY_HPP_
#define INCLUDE_VALARRAY_HPP_

#include <algorithm>
#include <vector>

#include <Interval.hpp>

template<typename >
struct ValArrayBase;

template<typename T, typename Allocator = std::allocator<T>>
struct ValArray;

template<typename T, typename O>
struct UnaryOp: public ValArrayBase<T> {
	using base_type = ValArrayBase<T>;
	using op_type = O;
	using value_type = base_type::value_type;
	using pointer = base_type::pointer;
	using reference = base_type::reference;
	using const_pointer = base_type::const_pointer;
	using const_reference = base_type::const_reference;
	using size_type = base_type::size_type;
	using difference_type = base_type::difference_type;
	UnaryOp(base_type const &a, op_type op) :
			Op(op), A(a) {
	}
	virtual ~UnaryOp() = default;
	virtual value_type operator[](size_type pos) const {
		return Op(A[pos]);
	}
	virtual size_type size() const {
		return A.size();
	}
private:
	op_type const &Op;
	base_type const &A;
};

template<typename T, typename O>
struct BinaryOp: public ValArrayBase<T> {
	using base_type = ValArrayBase<T>;
	using op_type = O;
	using value_type = base_type::value_type;
	using pointer = base_type::pointer;
	using reference = base_type::reference;
	using const_pointer = base_type::const_pointer;
	using const_reference = base_type::const_reference;
	using size_type = base_type::size_type;
	using difference_type = base_type::difference_type;
	BinaryOp(base_type const &a, base_type const &b, op_type op) :
			Op(op), A(a), B(b) {
		assert(A.size() == B.size());
	}
	virtual ~BinaryOp() = default;
	virtual value_type operator[](size_type pos) const {
		return Op(A[pos], B[pos]);
	}
	virtual size_type size() const {
		return A.size();
	}
private:
	op_type const &Op;
	base_type const &A;
	base_type const &B;
};

template<typename T>
struct Shift: public ValArrayBase<T> {
	using base_type = ValArrayBase<T>;
	using value_type = base_type::value_type;
	using pointer = base_type::pointer;
	using reference = base_type::reference;
	using const_pointer = base_type::const_pointer;
	using const_reference = base_type::const_reference;
	using size_type = base_type::size_type;
	using difference_type = base_type::difference_type;
	Shift(base_type const &a, size_type s) :
			A(a), shift(s) {
	}
	virtual ~Shift() = default;
	virtual value_type operator[](size_type pos) const {
		static constexpr T zero = T(0);
		difference_type const index = difference_type(pos) - difference_type(shift);
		if (index >= 0 && index < size()) {
			return A[index];
		} else {
			return zero;
		}
	}
	virtual size_type size() const {
		return A.size();
	}
private:
	base_type const &A;
	size_type shift;
};

template<typename T>
struct Cshift: public ValArrayBase<T> {
	using base_type = ValArrayBase<T>;
	using value_type = base_type::value_type;
	using pointer = base_type::pointer;
	using reference = base_type::reference;
	using const_pointer = base_type::const_pointer;
	using const_reference = base_type::const_reference;
	using size_type = base_type::size_type;
	using difference_type = base_type::difference_type;
	Cshift(base_type const &a, size_type s) :
			A(a), shift(s) {
	}
	virtual ~Cshift() = default;
	virtual value_type operator[](size_type pos) const {
		static constexpr T zero = T(0);
		difference_type index = difference_type(pos) - difference_type(shift);
		if (index < 0) {
			index += size();
		} else if (index >= size()) {
			index -= size();
		}
		return A[index];
	}
	virtual size_type size() const {
		return A.size();
	}
private:
	base_type const &A;
	size_type shift;
};

template<typename T>
struct ValArrayBase {
	using value_type = T;
	using pointer = T*;
	using reference = T&;
	using const_pointer = T const*;
	using const_reference = T const&;
	using size_type = std::size_t;
	using difference_type = std::ptrdiff_t;
	ValArrayBase() = default;
	ValArrayBase<T> operator+() const {
		return *this;
	}
	UnaryOp<T, std::plus<T>> operator-() const {
		return UnaryOp<T, std::negate<T>>(*this);
	}
	BinaryOp<T, std::plus<T>> operator+(ValArrayBase<T> const &other) const {
		return BinaryOp<T, std::plus<T>>(*this, other);
	}
	BinaryOp<T, std::minus<T>> operator-(ValArrayBase<T> const &other) const {
		return BinaryOp<T, std::minus<T>>(*this, other);
	}
	BinaryOp<T, std::multiplies<T>> operator*(ValArrayBase<T> const &other) const {
		return BinaryOp<T, std::multiplies<T>>(*this, other);
	}
	BinaryOp<T, std::divides<T>> operator/(ValArrayBase<T> const &other) const {
		return BinaryOp<T, std::divides<T>>(*this, other);
	}
	BinaryOp<T, std::modulus<T>> operator%(ValArrayBase<T> const &other) const {
		return BinaryOp<T, std::modulus<T>>(*this, other);
	}
	Shift<T> shift(size_type k) const {
		return Shift<T>(*this, k);
	}
	Cshift<T> cshift(size_type k) const {
		return Cshift<T>(*this, k);
	}
	virtual ~ValArrayBase() = default;
	virtual T operator[](size_type pos) const = 0;
	virtual size_type size() const = 0;
};

template<typename T, typename Allocator>
struct ValArray: public ValArrayBase<T> {
	using base_type = ValArrayBase<T>;
	using value_type = base_type::value_type;
	using pointer = base_type::pointer;
	using reference = base_type::reference;
	using const_pointer = base_type::const_pointer;
	using const_reference = base_type::const_reference;
	using size_type = base_type::size_type;
	using difference_type = base_type::difference_type;
	using iterator = std::vector<T>::iterator;
	using const_iterator = std::vector<T>::const_iterator;
	ValArray(size_type count = 0, const T &value = T(), const Allocator &alloc = Allocator()) :
			V(count, value, alloc) {
	}
	ValArray(const Allocator &alloc) :
			V(alloc) {
	}
	ValArray(const ValArray &other, const Allocator &alloc) :
			V(other.V, alloc) {
	}
	ValArray(ValArray &&other, const Allocator &alloc) :
			V(std::move(other.V), alloc) {
	}
	ValArray(std::initializer_list<T> init, const Allocator &alloc = Allocator()) :
			V(init, alloc) {
	}
	template<class Iterator>
	ValArray(Iterator first, Iterator last, const Allocator &alloc = Allocator()) :
			V(first, last, alloc) {
	}
	ValArray(ValArrayBase<T> const &op) {
		*this = op;
	}
	ValArray(ValArrayBase<T> &&op) {
		*this = std::move(op);
	}
	ValArray& operator=(ValArrayBase<T> const &op) {
		V.resize(op.size());
		for (size_type i = 0; i != op.size(); i++) {
			V[i] = op[i];
		}
		return *this;
	}
	ValArray& operator=(ValArrayBase<T> &&op) {
		V.resize(op.size());
		for (size_type i = 0; i != op.size(); i++) {
			V[i] = std::move(op[i]);
		}
		return *this;
	}
	value_type operator[](size_type pos) const {
		return V[pos];
	}
	iterator begin() {
		return V.begin();
	}
	iterator end() {
		return V.end();
	}
	iterator rbegin() {
		return V.rbegin();
	}
	iterator rend() {
		return V.rend();
	}
	const_iterator begin() const {
		return V.begin();
	}
	const_iterator end() const {
		return V.end();
	}
	const_iterator rbegin() const {
		return V.rbegin();
	}
	const_iterator rend() const {
		return V.rend();
	}
	virtual ~ValArray() = default;
	virtual reference operator[](size_type pos) {
		return V[pos];
	}
	virtual size_type size() const {
		return V.size();
	}
private:
	std::vector<T> V;
};

#endif /* INCLUDE_VALARRAY_HPP_ */
