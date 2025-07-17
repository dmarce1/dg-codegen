/******************************************************************************
 Copyright (C) 2024  Dominic C. Marcello
 *******************************************************************************/

#ifndef INCLUDE_FACE_HPP_
#define INCLUDE_FACE_HPP_

#include <cstdint>

template<int8_t D>
struct Face {
	static constexpr int8_t one = int8_t(1);
	Face() = default;
	Face(int8_t f) : face_(f) {}
	Face(int8_t dim, int8_t dir) {
		face_ = (dim << one) ^ ((dir + one) >> one);
	}
	Face& operator=(int8_t f) {
		face_ = f;
		return *this;
	}
	Face<D>& operator++() {
		face_++;
		return *this;
	}
	Face<D> operator++(int) {
		auto rc = *this;
		face_++;
		return rc;
	}
	Face<D>& operator--() {
		face_--;
		return *this;
	}
	Face<D> operator--(int) {
		auto rc = *this;
		face_--;
		return rc;
	}
	Face flip() const {
		Face f;
		f.face_ = face_ ^ one;
		return f;
	}
	int8_t dimension() const {
		return face_ >> one;
	}
	int8_t direction() const {
		return ((face_ & one) << one) - one;
	}
	operator int8_t() const {
		return face_;
	}
	static constexpr int8_t count() {
		return 2 * D;
	}
	static constexpr Face begin() {
		Face f;
		f.face_ = 0;
		return f;
	}
	static constexpr Face end() {
		Face f;
		f.face_ = count();
		return f;
	}
private:
	int8_t face_;
};

template<int8_t D>
struct Child {
	static constexpr int8_t zero = int8_t(0);
	static constexpr int8_t one = int8_t(1);
	Child() = default;
	Child(int8_t bits) :
			child_(bits) {
	}
	Child(Child<D - 1> child, Face<D> face) {
		child_ = zero;
		for (int8_t dim = 0; dim < D; dim++) {
			if (dim != face.dimension()) {
				child_ &= (child.child_ & (one << dim)) << int8_t(dim > face.dimension());
			} else {
				child_ &= int8_t(face.direction() > zero) << dim;
			}
		}
	}
	Child flip(int8_t dim = -1) const {
		Child c;
		if (dim >= 0) {
			c.child_ = child_ ^ (one << dim);
		} else {
			c.child_ = child_ ^ (count() - one);
		}
		return c;
	}
	int8_t direction(int8_t dim) const {
		return (((child_ >> dim) & one) << one) - one;
	}
	bool onFace(Face<D> face) const {
		return direction(face.dimension()) == face.direction();
	}
	operator int8_t() const {
		return child_;
	}
	Child<D>& operator++() {
		child_++;
		return *this;
	}
	Child<D> operator++(int) {
		auto rc = *this;
		child_++;
		return rc;
	}
	Child<D>& operator--() {
		child_--;
		return *this;
	}
	Child<D> operator--(int) {
		auto rc = *this;
		child_--;
		return rc;
	}
	static constexpr int8_t count() {
		return one << D;
	}
	static constexpr Child begin() {
		Child c;
		c.child_ = 0;
		return c;
	}
	static constexpr Child end() {
		Child c;
		c.child_ = count();
		return c;
	}
	template<int8_t D1>
	friend class Child;
private:
	int8_t child_;
};

#endif /* INCLUDE_FACE_HPP_ */
