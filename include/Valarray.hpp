#pragma once

#include <algorithm>
#include <cassert>
#include <cmath>
#include <iostream>
#include <memory>
#include <numeric>
#include <sstream>
#include <string>
#include <type_traits>
#include <vector>

#include "Util.hpp"

#define INDIRECTARRAY_COMPOUND_ASSIGNMENT(op)           \
	IndirectArray& operator op##= (auto const& _src) {     \
		auto const src = toExpression(_src);            \
		size_t k = 0;                                   \
		for(size_t j = 0; j != indices_->size(); j++) { \
			(*ptr_)[(*indices_)[j]] op##= src[j];                   \
		}                                               \
		return *this;                                   \
	}

#define GSLICEARRAY_COMPOUND_ASSIGNMENT(op)           \
	GSliceArray& operator op##= (auto const& _src) {     \
		auto const src = toExpression(_src);\
		auto const &strides = slc_.stride();\
		auto const &sizes = slc_.size();\
		size_t totalSize = sizes[0];\
		int const ndim = sizes.size();\
		for (int i = 1; i < ndim; i++) {\
			totalSize *= sizes[i];\
		}\
		Valarray<size_t> indices(size_t(0), ndim);\
		size_t k = slc_.start();\
		(*ptr_)[k] op##= src[0];\
		for (size_t j = 1; j < totalSize; j++) {\
			int dim = ndim - 1;\
			while (++indices[dim] == sizes[dim]) {\
				indices[dim] = 0;\
				k -= (sizes[dim] - 1) * strides[dim];\
				dim--;\
			}\
			k += strides[dim];\
			(*ptr_)[k] op##= src[j];\
		}\
		return *this; \
	}\

#define MASKARRAY_COMPOUND_ASSIGNMENT(op)\
	MaskArray& operator op##= (auto const& _src) { \
		auto const src = toExpression(_src);    \
		size_t const cnt = ptr_.size();         \
		size_t k = 0;                           \
		for(size_t j = 0; j != cnt; j++) {      \
			if((*mask_)[j]) {                   \
				(*ptr_)[j] op##= src[k];           \
				k++;                            \
			}                                   \
		}                                       \
		return *this;                           \
	}

#define SLICEARRAY_COMPOUND_ASSIGNMENT(op)                    \
	SliceArray& operator op##= (auto const& _src) {              \
		auto const src = toExpression(_src);                  \
		size_t const beg = slc_.start();                      \
		size_t const str = slc_.stride();                     \
		size_t const cnt = slc_.size();                       \
		for(size_t j = 0, k = beg; j != cnt; j++, k += str) { \
			(*ptr_)[k] op##= src[j];                             \
		}                                                     \
		return *this;                                         \
	}

#define VALARRAY_BINARY(_class, name, op)                                        \
	template<typename Ea, typename Eb>                                           \
	struct name##Expression: public Expression<T, name##Expression<Ea, Eb>> {       \
		using base_type = Expression<T, name##Expression<Ea, Eb>>;\
		name##Expression(Ea const& e1, Eb const& e2) : expr1_(e1), expr2_(e2) {  \
		}                                                                        \
		auto access(size_t index) const {                                        \
			return expr1_[index] op expr2_[index];                               \
		}                                                                        \
		size_t count() const {\
			return std::max(expr1_.size(), expr2_.size()); \
		}\
	private:                                                                     \
		typename std::remove_reference<Ea>::type expr1_;                                                               \
		typename std::remove_reference<Eb>::type expr2_;                                                               \
	};                                                                           \
	template<typename Derived1>                                                  \
	auto operator op (Expression<T, Derived1> const& expr) const {                       \
		auto const expr1 = toExpression(*this); 								 \
		auto const expr2 = toExpression(expr); 								     \
		return name##Expression(expr1, ((Derived1 const&) expr2));                                   \
	}                                                                            \
	auto operator op (Valarray<T> const& varr) const {                                \
		auto const expr1 = toExpression(*this); 								 \
		auto const expr2 = toExpression(varr); 								     \
		return name##Expression(expr1, expr2);                                   \
	}                                                                            \
	auto operator op (T const& val) const {                                           \
		auto const expr1 = toExpression(*this); 								 \
		auto const expr2 = toExpression(val); 								     \
		return name##Expression(expr1, expr2);                                   \
	}                                                                            \
	friend auto operator op (T const& val, _class const& obj) {                                           \
		auto const expr1 = toExpression(val); 								     \
		auto const expr2 = toExpression(obj); 								 \
		return name##Expression(expr1, expr2);                                   \
	}                                                                            \


#define VALARRAY_BINARY_FUNCTION2(name, func)                                    \
	template<typename Ea, typename Eb>                                           \
	struct name##Expression: public Expression<typename ElementType<Ea>::type, name##Expression<Ea, Eb>> {       \
		using base_type = Expression<typename ElementType<Ea>::type, name##Expression<Ea, Eb>>;\
		name##Expression(Ea const& e1, Eb const& e2) : expr1_(e1), expr2_(e2) {  \
		}                                                                        \
		auto access(size_t index) const {                                  \
			using std::min; \
			using std::max; \
			return func(expr1_[index], expr2_[index]);                           \
		}                                                                        \
		size_t count() const {\
			size_t const count = std::max(expr1_.size(), expr2_.size());\
			return count; \
		}\
	private:                                                                     \
		typename std::remove_reference<Ea>::type expr1_;                                                               \
		typename std::remove_reference<Eb>::type expr2_;                                                               \
	};                                                                           \
	template<typename T, typename Derived1, typename Derived2>\
	auto name (Expression<T, Derived1> const& arg1, Expression<T, Derived2> const& arg2) {                                  \
		auto const expr1 = toExpression(arg1); 									 \
		auto const expr2 = toExpression(arg2);                                   \
		return name##Expression(expr1, expr2);                                   \
	}                                                                            \
	template<typename T, typename Derived>\
	auto name (Expression<T, Derived> const& arg1, Valarray<T> const& arg2) {                                  \
		auto const expr1 = toExpression(arg1); 									 \
		auto const expr2 = toExpression(arg2);                                   \
		return name##Expression((Derived const&) expr1, expr2);                                   \
	}                                                                            \
	template<typename T, typename Derived>\
	auto name (Valarray<T> const& arg1, Expression<T, Derived> const& arg2) {                                  \
		auto const expr1 = toExpression(arg1); 									 \
		auto const expr2 = toExpression(arg2);                                   \
		return name##Expression(expr1, expr2);                                   \
	}                                                                            \
	template<typename T>\
	auto name (Valarray<T> const& arg1, Valarray<T> const& arg2) {                                  \
		auto const expr1 = toExpression(arg1); 									 \
		auto const expr2 = toExpression(arg2);                                   \
		return name##Expression(expr1, expr2);                                   \
	}                                                                            \
	template<typename T, typename Derived>\
	auto name (Expression<T, Derived> const& arg1, T const& arg2) {                                  \
		auto const expr1 = toExpression(arg1); 									 \
		auto const expr2 = toExpression(arg2);                                   \
		return name##Expression(expr1, expr2);                                   \
	}                                                                            \
	template<typename T, typename Derived>\
	auto name (T const& arg1, Expression<T, Derived>  const& arg2) {                                  \
		auto const expr1 = toExpression(arg1); 									 \
		auto const expr2 = toExpression(arg2);                                   \
		return name##Expression(expr1, expr2);                                   \
	}                                                                            \
	template<typename T>\
	auto name (T const& arg1, Valarray<T>  const& arg2) {                                  \
		auto const expr1 = toExpression(arg1); 									 \
		auto const expr2 = toExpression(arg2);                                   \
		return name##Expression(expr1, expr2);                                   \
	}                                                                            \
	template<typename T>\
	auto name (Valarray<T> const& arg1, T const& arg2) {                                  \
		auto const expr1 = toExpression(arg1); 									 \
		auto const expr2 = toExpression(arg2);                                   \
		return name##Expression(expr1, expr2);                                   \
	}                                                                            \

#define VALARRAY_BINARY_FUNCTION1(func) VALARRAY_BINARY_FUNCTION2(func, func)

#define VALARRAY_COMPARE(_class, name, op)                                              \
	template<typename Ea, typename Eb>                                          \
	struct name##Expression: public Expression<T, name##Expression<Ea, Eb>> {      \
		using base_type = Expression<T, name##Expression<Ea, Eb>>;\
		name##Expression(Ea const& e1, Eb const& e2) : expr1_(e1), expr2_(e2) { \
		}                                                                       \
		auto access(size_t index) const {                                       \
			return expr1_[index] op expr2_[index];                              \
		}                                                                       \
		size_t count() const {\
			return std::max(expr1_.size(), expr2_.size()); \
		}\
	private:                                                                    \
		typename std::remove_reference<Ea>::type expr1_;                                                               \
		typename std::remove_reference<Eb>::type expr2_;                                                               \
	};                                                                          \
	template<typename Derived1>                                                  \
	auto operator op (Expression<T, Derived1> const& expr2)const  {                       \
		auto const expr1 = toExpression(*this); 								 \
		return name##Expression(expr1, (Derived1 const&) expr2);                                   \
	}                                                                            \
	auto operator op (Valarray<T> const& varr) const  {                                \
		auto const expr1 = toExpression(*this); 								 \
		auto const expr2 = toExpression(varr); 								     \
		return name##Expression(expr1, expr2);                                   \
	}                                                                            \
	auto operator op (T const& val) const {                                           \
		auto const expr1 = toExpression(*this); 								 \
		auto const expr2 = toExpression(val); 								     \
		return name##Expression(expr1, expr2);                                   \
	}                                                                            \
	friend auto operator op (T const& val, _class const& obj) {                                           \
		auto const expr1 = toExpression(val); 								     \
		auto const expr2 = toExpression(obj); 								 \
		return name##Expression(expr1, expr2);                                   \
	}                                                                            \

#define VALARRAY_COMPOUND_ASSIGNMENT(op)                                       	\
	Valarray<T>& operator op##= (Valarray<T> const& _src) {                                 	\
		auto const src = ValarrayExpression<T>(_src);                                   	\
		size_t const count = src.size();                                       	\
		for(size_t i = 0; i != count; i++) {                                   	\
			(*ptr_)[i] op##= src[i];                                              	\
		}                                                                      	\
		return *this;  \
	}\
	template<typename Derived> \
	Valarray<T>& operator op##= (Expression<T, Derived> const& src) {                                 	\
		size_t const count = src.size();                                           	\
		for(size_t i = 0; i != count; i++) {                                   	\
			(*ptr_)[i] op##= src[i];                                              	\
		}                                                                      	\
		return *this;                                                          	\
	}                                                                          	\
	Valarray<T>& operator op##= (T const& src) {                                 	\
		size_t const count = size();                                           	\
		for(size_t i = 0; i != count; i++) {                                   	\
			(*ptr_)[i] op##= src;                                              	\
		}                                                                      	\
		return *this;                                                          	\
	}                                                                          	\

#define VALARRAY_UNARY(_class, name, op)                                                  \
	template<typename E>                                                          \
	struct name##Expression: public Expression<T, name##Expression<E>> {             \
		using base_type = Expression<T, name##Expression<E>>;\
		name##Expression(E const& e) : expr_(e) {                                 \
		}                                                                         \
		auto access(size_t index) const {                                   \
			return op expr_[index];                                               \
		}                                                                         \
		size_t count() const {      \
			return expr_.size(); \
		}\
	private:                                                                      \
		typename std::remove_reference<E>::type expr_;                                                               \
	};                                                                            \
	auto operator op () const {                                           \
		auto const expr = toExpression(*this); 								 \
		return name##Expression(expr);                                   \
	}                                                                            \

#define VALARRAY_UNARY_FUNCTION2(name, func)                                      \
	template<typename E>                                                          \
	struct name##Expression: public Expression<typename ElementType<E>::type, name##Expression<E>> {             \
		using base_type = Expression<typename ElementType<E>::type, name##Expression<E>>; \
		name##Expression(E const& e) : expr_(e) {                                 \
		}                                                                         \
		auto access(size_t index) const {                                   \
			return func (expr_[index]);                                           \
		}                                                                         \
		size_t count() const {      \
			return expr_.size(); \
		}\
	private:                                                                      \
		typename std::remove_reference<E>::type expr_;                                                               \
	};                                                                            \
	template<typename T1, typename Derived> \
	auto name (Expression<T1, Derived> const& arg) {                                                   \
		auto const expr = toExpression(arg); 									  \
		return name##Expression(expr);                                            \
	}                                                                             \
	template<typename T1> \
	auto name (Valarray<T1> const& arg) {                                                   \
		auto const expr = toExpression(arg); 									  \
		return name##Expression(expr);                                            \
	}                                                                             \

#define VALARRAY_UNARY_FUNCTION1(func) VALARRAY_UNARY_FUNCTION2(func, func)

struct GSlice;

struct Slice;

template<typename, typename >
struct Expression;

template<typename T>
auto toExpression(T const&);

template<typename >
struct IsExpression;

template<typename >
struct IsValarray;

template<class T>
struct IsValarrayExpr;

template<typename >
struct Valarray;

template<typename >
struct IndirectArray;

template<typename >
struct MaskArray;

template<typename >
struct GSliceArray;

template<typename >
struct SliceArray;

template<typename >
struct ValarrayExpression;

template<typename T>
struct Valarray {
	using value_type = std::vector<T>::value_type;
	using const_pointer = std::vector<T>::const_pointer;
	using const_reference = std::vector<T>::const_reference;
	using pointer = std::vector<T>::pointer;
	using reference = std::vector<T>::reference;
	using iterator = std::vector<T>::iterator;
	using const_iterator = std::vector<T>::const_iterator;
	Valarray();
	Valarray(size_t);
	Valarray(const_reference, size_t);
	Valarray(const_pointer, size_t);
	Valarray(const Valarray&);
	Valarray(Valarray&&);
	Valarray(std::initializer_list<T>);
	template<size_t count>
	Valarray(std::array<T, count> const&);
	Valarray(GSliceArray<T> const&);
	Valarray(SliceArray<T> const&);
	Valarray(IndirectArray<T> const&);
	Valarray(MaskArray<T> const&);
	template<typename U>
	Valarray(Expression<T, U> const&);
	virtual ~Valarray();
	Valarray& operator=(const Valarray&);
	Valarray& operator=(Valarray&&);
	Valarray& operator=(std::initializer_list<T>);
	Valarray& operator=(GSliceArray<T> const&);
	Valarray& operator=(SliceArray<T> const&);
	Valarray& operator=(IndirectArray<T> const&);
	Valarray& operator=(MaskArray<T> const&);
	reference operator[](size_t);
	T operator[](size_t) const;
	SliceArray<T> operator[](Slice const&) const;
	GSliceArray<T> operator[](GSlice const&) const;
	MaskArray<T> operator[](Valarray<bool> const&) const;
	IndirectArray<T> operator[](Valarray<size_t> const&) const;
	void swap(Valarray&&);
	value_type sum() const;
	value_type min() const;
	value_type max() const;
	Valarray shift(int shft) const;
	Valarray cshift(int shft) const;
	Valarray apply(T(T)) const;
	Valarray apply(T(T const&)) const;
	size_t size() const;
	void fill(T const&);
	iterator begin();
	iterator end();
	const_iterator begin() const;
	const_iterator end() const;
	template<typename T1>
	Valarray<T1> cast() const;
	template<typename Derived>
	Valarray<T>& operator=(Expression<T, Derived> const&);
	Valarray<T>& operator=(T const&);
	/**/
	VALARRAY_COMPOUND_ASSIGNMENT(+)
	VALARRAY_COMPOUND_ASSIGNMENT(-)
	VALARRAY_COMPOUND_ASSIGNMENT(*)
	VALARRAY_COMPOUND_ASSIGNMENT(/)
	VALARRAY_COMPOUND_ASSIGNMENT(%)
	VALARRAY_COMPOUND_ASSIGNMENT(&)
	VALARRAY_COMPOUND_ASSIGNMENT(|)
	VALARRAY_COMPOUND_ASSIGNMENT(^)
	VALARRAY_COMPOUND_ASSIGNMENT(>>)
	VALARRAY_COMPOUND_ASSIGNMENT(<<)
	VALARRAY_BINARY(Valarray, Add, +)
	VALARRAY_BINARY(Valarray, Subtract, -)
	VALARRAY_BINARY(Valarray, Multiply, *)
	VALARRAY_BINARY(Valarray, Divide, /)
	VALARRAY_BINARY(Valarray, Modulus, %)
	VALARRAY_BINARY(Valarray, LogicalAnd, &&)
	VALARRAY_BINARY(Valarray, LogicalOr, ||)
	VALARRAY_BINARY(Valarray, BitwiseAnd, &)
	VALARRAY_BINARY(Valarray, BitwiseOr, |)
	VALARRAY_BINARY(Valarray, BitwiseXor, ^)
	VALARRAY_BINARY(Valarray, ShiftRight, >>)
	VALARRAY_BINARY(Valarray, ShiftLeft, <<)
	VALARRAY_COMPARE(Valarray, Equal, ==)
	VALARRAY_COMPARE(Valarray, NotEqual, !=)
	VALARRAY_COMPARE(Valarray, Greater, >)
	VALARRAY_COMPARE(Valarray, Less, <)
	VALARRAY_COMPARE(Valarray, GreaterOrEqual, >=)
	VALARRAY_COMPARE(Valarray, LessOrEqual, <=)
	VALARRAY_UNARY(Valarray, Positive, +)
	VALARRAY_UNARY(Valarray, Negative, -)
	VALARRAY_UNARY(Valarray, LogicalNot, !)
	VALARRAY_UNARY(Valarray, BitwiseNot, ~)
	friend class ValarrayExpression<T> ;
private:
	std::shared_ptr<std::vector<T>> ptr_;
};

template<typename T>
Valarray<T>& Valarray<T>::operator=(Valarray<T> const &src) {
	size_t const count = src.size();
	ptr_ = std::make_shared<std::vector<T>>(std::vector<T>(count));
	for (size_t i = 0; i != count; i++) {
		(*ptr_)[i] = src[i];
	}
	return *this;
}

template<typename T>
template<typename Derived>
Valarray<T>& Valarray<T>::operator=(Expression<T, Derived> const &src) {
	size_t const count = src.size();
	ptr_ = std::make_shared<std::vector<T>>(count);
	for (size_t i = 0; i != count; i++) {
		(*ptr_)[i] = src[i];
	}
	return *this;
}

template<typename T>
Valarray<T>& Valarray<T>::operator=(T const &val) {
	size_t const count = size();
	for (size_t i = 0; i != count; i++) {
		(*ptr_)[i] = val;
	}
	return *this;
}

struct Slice {
	Slice() = default;
	Slice(Slice const&) = default;
	Slice(size_t, size_t, size_t);
	size_t start() const;
	size_t size() const;
	size_t stride() const;
private:
	size_t start_;
	size_t size_;
	size_t stride_;
};

class GSlice {
public:
	GSlice() = default;
	GSlice(GSlice const&) = default;
	GSlice(size_t, Valarray<size_t> const&, Valarray<size_t> const&);
	size_t start() const;
	Valarray<size_t> const& size() const;
	Valarray<size_t> const& stride() const;
private:
	size_t start_;
	Valarray<size_t> sizes_;
	Valarray<size_t> strides_;
};

template<typename T>
class SliceArray {
	SliceArray(std::shared_ptr<std::vector<T>> const&, Slice const&);
	size_t size() const;
public:
	SliceArray() = default;
	SliceArray(SliceArray const&) = default;
	friend class Valarray<T> ;
	SliceArray& operator=(auto const&);
	SLICEARRAY_COMPOUND_ASSIGNMENT(+)
	SLICEARRAY_COMPOUND_ASSIGNMENT(-)
	SLICEARRAY_COMPOUND_ASSIGNMENT(*)
	SLICEARRAY_COMPOUND_ASSIGNMENT(/)
	SLICEARRAY_COMPOUND_ASSIGNMENT(%)
	SLICEARRAY_COMPOUND_ASSIGNMENT(&)
	SLICEARRAY_COMPOUND_ASSIGNMENT(|)
	SLICEARRAY_COMPOUND_ASSIGNMENT(^)
	SLICEARRAY_COMPOUND_ASSIGNMENT(>>)
	SLICEARRAY_COMPOUND_ASSIGNMENT(<<)
private:
	std::shared_ptr<std::vector<T>> ptr_;
	Slice slc_;
};

template<typename T>
class GSliceArray {
	GSliceArray(std::shared_ptr<std::vector<T>> const&, GSlice const&);
	size_t size() const;
public:
	GSliceArray() = default;
	GSliceArray(GSliceArray const&) = default;
	friend class Valarray<T> ;
	GSliceArray& operator=(auto const&);
	GSLICEARRAY_COMPOUND_ASSIGNMENT(+)
	GSLICEARRAY_COMPOUND_ASSIGNMENT(-)
	GSLICEARRAY_COMPOUND_ASSIGNMENT(*)
	GSLICEARRAY_COMPOUND_ASSIGNMENT(/)
	GSLICEARRAY_COMPOUND_ASSIGNMENT(%)
	GSLICEARRAY_COMPOUND_ASSIGNMENT(&)
	GSLICEARRAY_COMPOUND_ASSIGNMENT(|)
	GSLICEARRAY_COMPOUND_ASSIGNMENT(^)
	GSLICEARRAY_COMPOUND_ASSIGNMENT(>>)
	GSLICEARRAY_COMPOUND_ASSIGNMENT(<<)
private:
	std::shared_ptr<std::vector<T>> ptr_;
	GSlice slc_;
};

template<typename T>
class MaskArray {
	MaskArray(Valarray<T> const&, Valarray<bool> const&);
	size_t size() const;
public:
	MaskArray() = default;
	MaskArray(MaskArray<T> const&) = default;
	friend class Valarray<T> ;
	MaskArray& operator =(auto const&);
	MASKARRAY_COMPOUND_ASSIGNMENT(=)
	MASKARRAY_COMPOUND_ASSIGNMENT(+)
	MASKARRAY_COMPOUND_ASSIGNMENT(-)
	MASKARRAY_COMPOUND_ASSIGNMENT(*)
	MASKARRAY_COMPOUND_ASSIGNMENT(/)
	MASKARRAY_COMPOUND_ASSIGNMENT(%)
	MASKARRAY_COMPOUND_ASSIGNMENT(&)
	MASKARRAY_COMPOUND_ASSIGNMENT(|)
	MASKARRAY_COMPOUND_ASSIGNMENT(^)
	MASKARRAY_COMPOUND_ASSIGNMENT(>>)
	MASKARRAY_COMPOUND_ASSIGNMENT(<<)
private:
	std::shared_ptr<std::vector<T>> ptr_;
	std::shared_ptr<std::vector<bool>> mask_;
};
template<typename T>
class IndirectArray {
	IndirectArray(Valarray<T> const&, Valarray<size_t> const&);
	size_t size() const;
public:
	IndirectArray() = default;
	IndirectArray(IndirectArray<T> const&) = default;
	friend class Valarray<T> ;
	IndirectArray& operator=(auto const&);
	INDIRECTARRAY_COMPOUND_ASSIGNMENT(+)
	INDIRECTARRAY_COMPOUND_ASSIGNMENT(-)
	INDIRECTARRAY_COMPOUND_ASSIGNMENT(*)
	INDIRECTARRAY_COMPOUND_ASSIGNMENT(/)
	INDIRECTARRAY_COMPOUND_ASSIGNMENT(%)
	INDIRECTARRAY_COMPOUND_ASSIGNMENT(&)
	INDIRECTARRAY_COMPOUND_ASSIGNMENT(|)
	INDIRECTARRAY_COMPOUND_ASSIGNMENT(^)
	INDIRECTARRAY_COMPOUND_ASSIGNMENT(>>)
	INDIRECTARRAY_COMPOUND_ASSIGNMENT(<<)
private:
	std::shared_ptr<std::vector<T>> ptr_;
	std::shared_ptr<std::vector<size_t>> indices_;
};

template<typename T>
size_t IndirectArray<T>::size() const {
	return indices_->size();
}

template<typename T>
struct ValarrayExpression: public Expression<T, ValarrayExpression<T>> {
	using base_type = Expression<T, ValarrayExpression<T>>;
	ValarrayExpression() = delete;
	ValarrayExpression(ValarrayExpression<T> const&) = default;
	ValarrayExpression(ValarrayExpression<T>&&) = delete;
	ValarrayExpression& operator=(ValarrayExpression<T> const&) = delete;
	ValarrayExpression& operator=(ValarrayExpression<T>&&) = delete;
	ValarrayExpression(Valarray<T> const&);
	virtual ~ValarrayExpression();
	T access(size_t) const;
	size_t count() const;
private:
	std::shared_ptr<std::vector<T>> ptr_;
};

template<typename T>
struct ScalarExpression: public Expression<T, ScalarExpression<T>> {
	using base_type = Expression<T, ScalarExpression<T>>;
	T access(size_t) const;
	ScalarExpression(T const&);
	virtual ~ScalarExpression();
	size_t count() const;
	operator T() const;
private:
	std::remove_reference_t<T> value_;
};

template<typename T, typename Derived>
struct Expression {
	using value_type = std::vector<T>::value_type;
	using const_pointer = std::vector<T>::const_pointer;
	using const_reference = std::vector<T>::const_reference;
	using pointer = std::vector<T>::pointer;
	using reference = std::vector<T>::reference;
	using derived_type = Derived;
	Expression();
	virtual ~Expression();
	auto operator[](size_t) const;
	size_t size() const;
	template<typename T2>
	Valarray<T2> cast() const;VALARRAY_BINARY(Derived, Add, +)
	VALARRAY_BINARY(Derived, Subtract, -)
	VALARRAY_BINARY(Derived, Multiply, *)
	VALARRAY_BINARY(Derived, Divide, /)
	VALARRAY_BINARY(Derived, Modulus, %)
	VALARRAY_BINARY(Derived, LogicalAnd, &&)
	VALARRAY_BINARY(Derived, LogicalOr, ||)
	VALARRAY_BINARY(Derived, BitwiseAnd, &)
	VALARRAY_BINARY(Derived, BitwiseOr, |)
	VALARRAY_BINARY(Derived, BitwiseXor, ^)
	VALARRAY_BINARY(Derived, ShiftRight, >>)
	VALARRAY_BINARY(Derived, ShiftLeft, <<)
	VALARRAY_COMPARE(Derived, Equal, ==)
	VALARRAY_COMPARE(Derived, NotEqual, !=)
	VALARRAY_COMPARE(Derived, Greater, >)
	VALARRAY_COMPARE(Derived, Less, <)
	VALARRAY_COMPARE(Derived, GreaterOrEqual, >=)
	VALARRAY_COMPARE(Derived, LessOrEqual, <=)
	VALARRAY_UNARY(Derived, Positive, +)
	VALARRAY_UNARY(Derived, Negative, -)
	VALARRAY_UNARY(Derived, LogicalNot, !)
	VALARRAY_UNARY(Derived, BitwiseNot, ~)
};

template<typename T>
struct IsExpression {
private:
	template<typename U, typename V, typename D>
	static std::true_type probe(Expression<V, D> const*);
	template<typename >
	static std::false_type probe(...);

public:
	static constexpr bool value = decltype(probe<T>(static_cast<T const*>(nullptr)))::value;
};

template<typename T>
struct IsValarray {
	static constexpr bool value = false;
};

template<typename T>
struct IsValarray<Valarray<T>> {
	static constexpr bool value = true;
};

template<class T>
struct IsValarrayExpr {
	static constexpr bool value = IsExpression<T>::value || IsValarray<T>::value;
};
;

template<typename T>
IndirectArray<T>::IndirectArray(Valarray<T> const &ptr, Valarray<size_t> const &indices) :
		ptr_(ptr), indices_(indices) {
}

template<typename T>
MaskArray<T>::MaskArray(Valarray<T> const &ptr, Valarray<bool> const &mask) :
		ptr_(ptr), mask_(mask) {
}

template<typename T>
SliceArray<T>::SliceArray(std::shared_ptr<std::vector<T>> const &varray, Slice const &slc) :
		ptr_(varray), slc_(slc) {
}

template<typename T>
GSliceArray<T>::GSliceArray(std::shared_ptr<std::vector<T>> const &ptr, GSlice const &slc) :
		ptr_(ptr), slc_(slc) {
}

template<typename T>
Valarray<T>::Valarray() {
}

template<typename T>
Valarray<T>::Valarray(size_t count) :
		ptr_(new std::vector<T>(count)) {
}

template<typename T>
Valarray<T>::Valarray(const_reference val, size_t count) :
		ptr_(std::make_shared<std::vector<T>>(count, val)) {
}

template<typename T>
Valarray<T>::Valarray(const_pointer ptr, size_t count) :
		ptr_(std::make_shared<std::vector<T>>(ptr, ptr + count)) {
}

template<typename T>
Valarray<T>::Valarray(const Valarray &other) :
		ptr_(std::make_shared<std::vector<T>>(other.ptr_->begin(), other.ptr_->end())) {
}

template<typename T>
Valarray<T>::Valarray(Valarray &&other) :
		ptr_(nullptr) {
	*this = other;
}

template<typename T>
Valarray<T>::Valarray(std::initializer_list<T> list) :
		ptr_(new std::vector<T>(list.begin(), list.end())) {
}

template<typename T>
template<size_t count>
Valarray<T>::Valarray(std::array<T, count> const &arr) {
	ptr_->resize(count);
	std::copy(arr.begin(), arr.end(), ptr_->begin());
}

template<typename T>
template<typename U>
Valarray<T>::Valarray(Expression<T, U> const &expr) {
	*this = expr;
}

template<typename T>
Valarray<T>::Valarray(GSliceArray<T> const &other) {
	*this = other;
}

template<typename T>
Valarray<T>::Valarray(SliceArray<T> const &other) {
	*this = other;
}

template<typename T>
Valarray<T>::Valarray(IndirectArray<T> const &other) {
	*this = other;
}

template<typename T>
Valarray<T>::Valarray(MaskArray<T> const &other) {
	*this = other;
}

template<typename T>
Valarray<T>::~Valarray() {
}

template<typename T>
Valarray<T>& Valarray<T>::operator=(IndirectArray<T> const &src) {
	ptr_ = std::make_shared<std::vector<T>>(src.size());
	size_t k = 0;
	for (size_t j = 0; j != src.indices_->size(); j++) {
		(*this)[j] = (*src.ptr_)[k];
		k++;
	}
	return *this;
}

template<typename T>
Valarray<T>& Valarray<T>::operator=(GSliceArray<T> const &src) {
	size_t count = src.size();
	ptr_ = std::make_shared<std::vector<T>>(count);
	auto const &strides = src.slc_.stride();
	auto const &sizes = src.slc_.size();
	size_t const ndim = sizes.size();
	Valarray<size_t> indices(size_t(0), ndim);
	size_t k = src.slc_.start();
	(*this)[0] = (*src.ptr_)[k];
	for( size_t  j = 1; j < count; j++) {
		int dim = ndim - 1;
		while (++indices[dim] == sizes[dim]) {
			indices[dim] = 0;
			k -= (sizes[dim] - 1) * strides[dim];
			dim--;
		}
		k += strides[dim];
		(*this)[j] = (*src.ptr_)[k];
	}
	return *this;
}

template<typename T>
Valarray<T>& Valarray<T>::operator=(MaskArray<T> const &src) {
	size_t const cnt = src.size();
	ptr_ = std::make_shared<std::vector<T>>(src.size());
	size_t k = 0;
	for (size_t j = 0; j != cnt; j++) {
		if ((*this)[j]) {
			(*ptr_)[j] = (*src.ptr_)[k];
			k++;
		}
	}
	return *this;
}

template<typename T>
Valarray<T>& Valarray<T>::operator=(SliceArray<T> const &src) {
	ptr_ = std::make_shared<std::vector<T>>(src.size());
	size_t const beg = src.slc_.start();
	size_t const str = src.slc_.stride();
	size_t const cnt = src.slc_.size();
	for (size_t j = 0, k = beg; j != cnt; j++, k += str) {
		(*this)[j] = (*src.ptr_)[k];
	}
	return *this;
}

template<typename T>
SliceArray<T>& SliceArray<T>::operator=(auto const &_src) {
	auto const src = toExpression(_src);
	size_t const count = src.size();
//	ptr_ = std::make_shared<std::vector<T>>(std::vector<T>(count));
	size_t const beg = slc_.start();
	size_t const str = slc_.stride();
	size_t const cnt = slc_.size();
	for (size_t j = 0, k = beg; j != cnt; j++, k += str) {
		(*ptr_)[k] = src[j];
	}
	return *this;
}

template<typename T>
GSliceArray<T>& GSliceArray<T>::operator=(auto const &_src) {
	auto const src = toExpression(_src);
	auto const &strides = slc_.stride();
	auto const &sizes = slc_.size();
	size_t count = sizes[0];
	int const ndim = sizes.size();
	for (int i = 1; i < ndim; i++) {
		count *= sizes[i];
	}
//	ptr_ = std::make_shared<std::vector<T>>(std::vector<T>(count));
	Valarray<size_t> indices(size_t(0), ndim);
	size_t k = slc_.start();
	(*ptr_)[k] = src[0];
	for (size_t j = 1; j < count; j++) {
		int dim = ndim - 1;
		while (++indices[dim] == sizes[dim]) {
			indices[dim] = 0;
			k -= (sizes[dim] - 1) * strides[dim];
			dim--;
		}
		k += strides[dim];
		(*ptr_)[k] = src[j];
	}
	return *this;
}

template<typename T>
MaskArray<T>& MaskArray<T>::operator=(auto const &_src) {
	auto const src = toExpression(_src);
	size_t const count = ptr_.size();
//	ptr_ = std::make_shared<std::vector<T>>(std::vector<T>(count));
	size_t k = 0;
	for (size_t j = 0; j != count; j++) {
		if ((*mask_)[j]) {
			(*ptr_)[j] = src[k];
			k++;
		}
	}
	return *this;
}

template<typename T>
IndirectArray<T>& IndirectArray<T>::operator=(auto const &_src) {
	auto const src = toExpression(_src);
	size_t const count = src.size();
//	ptr_ = std::make_shared<std::vector<T>>(std::vector<T>(count));
	for (size_t j = 0; j != count; j++) {
		(*ptr_)[(*indices_)[j]] = src[j];
	}
	return *this;
}

//template<typename T>
//Valarray<T>& Valarray<T>::operator=(const_reference val) {
//	std::fill(ptr_->begin(), ptr_->end(), val);
//	return *this;
//}

//template<typename T>
//Valarray<T>& Valarray<T>::operator=(const Valarray &other) {
//	*ptr_ = *other.ptr_;
//	return *this;
//}

template<typename T>
Valarray<T>& Valarray<T>::operator=(Valarray &&other) {
	ptr_ = other.ptr_;
	other.ptr_ = nullptr;
	return *this;
}

template<typename T>
Valarray<T>& Valarray<T>::operator=(std::initializer_list<T> list) {
	*ptr_ = list;
	return *this;
}

//template<typename T>
//template<typename U>
//Valarray<T>& Valarray<T>::operator=(Expression<T, U> const &expr) {
//	size_t const count = expr.size();
//	ptr_->resize(count);
//	for (size_t i = 0; i != count; i++) {
//		(*ptr_)[i] = expr[i];
//	}
//	return *this;
//}

template<typename T>
typename Valarray<T>::reference Valarray<T>::operator[](size_t index) {
	return (*ptr_)[index];
}

template<typename T>
T Valarray<T>::operator[](size_t index) const {
	return (*ptr_)[index];
}

template<typename T>
SliceArray<T> Valarray<T>::operator[](Slice const &slc) const {
	return SliceArray<T>(ptr_, slc);
}

template<typename T>
GSliceArray<T> Valarray<T>::operator[](GSlice const &slc) const {
	return GSliceArray<T>(ptr_, slc);
}

template<typename T>
MaskArray<T> Valarray<T>::operator[](Valarray<bool> const &mask) const {
	return MaskArray<T>(ptr_, mask);
}

template<typename T>
void Valarray<T>::swap(Valarray<T> &&other) {
	std::swap(other.ptr_, ptr_);
}

template<typename T>
Valarray<T>::value_type Valarray<T>::sum() const {
	return std::accumulate(ptr_->begin(), ptr_->end(), T(0));
}

template<typename T>
Valarray<T>::value_type Valarray<T>::min() const {
	return *std::min_element(ptr_->begin(), ptr_->end());
}

template<typename T>
Valarray<T>::value_type Valarray<T>::max() const {
	return *std::max_element(ptr_->begin(), ptr_->end());
}

template<typename T>
Valarray<T> Valarray<T>::shift(int shft) const {
	Valarray arry(*this);
	std::shift_left(arry.begin(), arry.end(), shft);
}

template<typename T>
Valarray<T> Valarray<T>::cshift(int shft) const {
	Valarray arry(*this);
	if (shft > 0) {
		std::copy(ptr_->begin() + shft, ptr_->end(), arry.ptr_->begin());
		std::copy(ptr_->begin(), ptr_->begin() + shft, arry.ptr_->end() - shft);
	} else {
		std::copy(ptr_->begin(), ptr_->end() + shft, arry.ptr_->begin() - shft);
		std::copy(ptr_->end() + shft, ptr_->end(), arry.ptr_->begin());
	}
	return arry;
}

template<typename T>
Valarray<T> Valarray<T>::apply(T f(T)) const {
	Valarray arry(*this);
	std::transform(arry.begin(), arry.end(), f);
	return *this;
}

template<typename T>
Valarray<T> Valarray<T>::apply(T f(T const&)) const {
	Valarray arry(*this);
	std::transform(arry.begin(), arry.end(), f);
	return *this;
}

template<typename T>
Valarray<T>::iterator Valarray<T>::begin() {
	return ptr_->begin();
}

template<typename T>
Valarray<T>::iterator Valarray<T>::end() {
	return ptr_->end();
}

template<typename T>
Valarray<T>::const_iterator Valarray<T>::begin() const {
	return ptr_->begin();
}

template<typename T>
Valarray<T>::const_iterator Valarray<T>::end() const {
	return ptr_->end();
}

template<typename T1>
template<typename T2>
Valarray<T2> Valarray<T1>::cast() const {
	size_t const count = ptr_->size();
	Valarray<T2> results(count);
	for (size_t i = 0; i < count; i++) {
		results[i] = T2((*this)[i]);
	}
	return results;
}

template<typename T>
void Valarray<T>::fill(T const &val) {
	std::fill(begin(), end(), val);
}

template<typename T>
size_t Valarray<T>::size() const {
	return ptr_->size();
}

template<typename T>
ValarrayExpression<T>::ValarrayExpression(Valarray<T> const &array) :
		ptr_(array.ptr_) {
}

template<typename T>
ValarrayExpression<T>::~ValarrayExpression() {
}

template<typename T>
T ValarrayExpression<T>::access(size_t index) const {
	return (*ptr_)[index];
}

template<typename T>
size_t ValarrayExpression<T>::count() const {
	//assert((long) ptr_->size() > 0);
	return ptr_->size();
}

template<typename T>
ScalarExpression<T>::ScalarExpression(T const &value) :
		value_(value) {
}

template<typename T>
ScalarExpression<T>::~ScalarExpression() {
}

template<typename T>
T ScalarExpression<T>::access(size_t index) const {
	return value_;
}

template<typename T>
size_t ScalarExpression<T>::count() const {
	return 1;
}

template<typename T>
auto toExpression(T const &obj) {
	if constexpr (IsExpression<T>::value) {
		return (typename T::derived_type const&) obj;
	} else if constexpr (IsValarray<T>::value) {
		return ValarrayExpression<typename T::value_type>(obj);
	} else {
		using Type = typename std::remove_reference<T>::type;
		return ScalarExpression<Type>(obj);
	}
}

template<typename T, typename Derived>
Expression<T, Derived>::Expression() {
}

template<typename T1, typename Derived>
template<typename T2>
Valarray<T2> Expression<T1, Derived>::cast() const {
	size_t const count = size();
	Valarray<T2> results(count);
	for (size_t i = 0; i < count; i++) {
		results[i] = T2((*this)[i]);
	}
	return results;
}

template<typename T, typename Derived>
Expression<T, Derived>::~Expression() {
}

template<typename T, typename Derived>
auto Expression<T, Derived>::operator[](size_t index) const {
	return static_cast<derived_type const*>(this)->access(index);
}

template<typename T, typename Derived>
size_t Expression<T, Derived>::size() const {
	return static_cast<derived_type const*>(this)->count();
}

template<typename T>
size_t SliceArray<T>::size() const {
	return slc_.size();
}

template<typename T>
size_t MaskArray<T>::size() const {
	size_t sum = 0;
	for (auto const &bit : *mask_) {
		if (bit) {
			sum++;
		}
	}
	return sum;
}

template<typename T>
size_t GSliceArray<T>::size() const {
	auto const &sizes = slc_.size();
	if (sizes.size()) {
		constexpr std::multiplies<size_t> multOp;
		return std::accumulate(sizes.begin(), sizes.end(), size_t(1), multOp);
	} else {
		return 0;
	}
}

VALARRAY_BINARY_FUNCTION1(copysign)
VALARRAY_BINARY_FUNCTION1(max)
VALARRAY_BINARY_FUNCTION1(min)
VALARRAY_UNARY_FUNCTION1(abs)
VALARRAY_UNARY_FUNCTION1(sqrt)

template<typename T>
auto toString(Valarray<T> const &values) {
	std::ostringstream oss { };
	oss << '[';
	for (size_t i = 0; i != values.size(); i++) {
		oss << values[i];
		if (i + 1 != values.size()) {
			oss << ", ";
		}
	}
	oss << ']';
	return oss.str();
}

