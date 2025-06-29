#ifndef INCLUDE_VAL_ARRAY_HPP_
#define INCLUDE_VAL_ARRAY_HPP_

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <initializer_list>
#include <functional>
#include <memory>
#include <numeric>
#include <stack>
#include <type_traits>
#include <utility>

#define VAL_ARRAY_UNARY_OPERATOR(op)                                             \
	auto operator op() const {                                                    \
		using Handle = decltype(this->getHandle());                                \
		struct ThisExpression : public Expression<Type, elementCount, ThisExpression> { \
			ThisExpression(Handle const &handle) : handle_(handle) {                                      \
			}                                                                       \
			Type access(int64_t index) const {                                      \
				return op handle_[index];                                            \
			}                                                                       \
		private:                                                                   \
			Handle handle_;                                                         \
		};                                                                         \
		return ThisExpression(this->getHandle());                          \
	}

#define VAL_ARRAY_BINARY_OPERATOR1(op)                                                  \
	auto operator op(auto const& other) const {                                          \
		using HandleA = decltype(this->getHandle());                                      \
		using HandleB = decltype(other.getHandle());                                      \
		struct ThisExpression : public Expression<Type, elementCount, ThisExpression> {        \
			ThisExpression(HandleA const &handleA, HandleB const &handleB) : \
				handleA_(handleA), handleB_(handleB){                         \
			}                                                                              \
			Type access(int64_t index) const {                                             \
				return handleA_[index] op handleB_[index];                                  \
			}                                                                              \
		private:                                                                          \
			HandleA handleA_;                                                              \
			HandleB handleB_;                                                              \
		};                                                                                \
		return ThisExpression(this->getHandle(), other.getHandle());        \
	}                                                                                    \
	auto operator op(Type const& value) const {                                          \
		using Handle = decltype(this->getHandle());                                       \
		struct ThisExpression : public Expression<Type, elementCount, ThisExpression> {        \
			ThisExpression(Handle const &handle, Type const &value) :        \
				handle_(handle), value_(value) {                               \
			}                                                                              \
			Type access(int64_t index) const {                                             \
				return handle_[index] op value_;                                            \
			}                                                                              \
		private:                                                                          \
			Handle handle_;                                                                \
			Type value_;                                                                   \
		};                                                                                \
		return ThisExpression(this->getHandle(), value);                          \
	}                                                                                    \

#define VAL_ARRAY_BINARY_OPERATOR2(op)                                                  \
	template<typename Type, int64_t elementCount, typename Derived> \
	auto operator op(Type const& value, Expression<Type, elementCount, Derived> const& other) {                                          \
		using Handle = decltype(other.getHandle());                                      \
		struct ThisExpression : public Expression<Type, elementCount, ThisExpression> {        \
			ThisExpression(Type const& value, Handle const &handle) : \
				value_(value), handle_(handle){                         \
			}                                                                              \
			Type access(int64_t index) const {                                             \
				return value_ op handle_[index];                                  \
			}                                                                              \
		private:                                                                          \
			Type const& value_;\
			Handle handle_;                                                              \
		};                                                                                \
		return ThisExpression(value, other.getHandle());        \
	}                                                                                    \
	template<typename Type, int64_t elementCount> \
	auto operator op(Type const& value, ValArray<Type, elementCount> const& other) {                                          \
		using Handle = decltype(other.getHandle());                                      \
		struct ThisExpression : public Expression<Type, elementCount, ThisExpression> {        \
			ThisExpression(Type const& value, Handle const &handle) : \
				value_(value), handle_(handle){                         \
			}                                                                              \
			Type access(int64_t index) const {                                             \
				return value_ op handle_[index];                                  \
			}                                                                              \
		private:                                                                          \
			Type const& value_;\
			Handle handle_;                                                              \
		};                                                                                \
		return ThisExpression(value, other.getHandle());        \
	}                                                                                    \

#define VAL_ARRAY_COMPARE_OPERATOR(op)                                                  \
	auto operator op(auto const& other) const {                                          \
		using HandleA = decltype(this->getHandle());                                      \
		using HandleB = decltype(other.getHandle());                                      \
		struct ThisExpression : public Expression<bool, elementCount, ThisExpression> {        \
			ThisExpression(HandleA const &handleA, HandleB const &handleB) : \
				handleA_(handleA), handleB_(handleB) {                         \
			}                                                                              \
			bool operator[](int64_t index) {                                               \
				return handleA_[index] op handleB_[index];                                  \
			}                                                                              \
		private:                                                                          \
			HandleA handleA_;                                                              \
			HandleB handleB_;                                                              \
		};                                                                                \
		return ThisExpression(this->getHandle(), other.getHandle());        \
	}

#define VAL_ARRAY_UNARY_FUNCTION(class_, name)                                   \
	friend auto name(class_ const& object) {                                      \
		using Handle = decltype(object.getHandle());                               \
		struct ThisExpression : public Expression<Type, elementCount, ThisExpression> { \
			ThisExpression(Handle const &handle) :                    \
				handle_(handle) {                                       \
			}                                                                       \
			Type access(int64_t index) const {                                      \
				return name(handle_[index]);                                         \
			} 	                                                                     \
		private:                                                                   \
			Handle handle_;                                                         \
		};                                                                         \
		return ThisExpression(object.getHandle());                  \
	}

#define VAL_ARRAY_BINARY_FUNCTION(class_, name, func)                                   \
	friend auto name(class_ const& objectA, auto const& objectB) {                       \
		using HandleA = decltype(objectA.getHandle());                                    \
		using HandleB = decltype(objectB.getHandle());                                    \
		struct ThisExpression : public Expression<Type, elementCount, ThisExpression> {        \
			ThisExpression(HandleA const &handleA, HandleB const &handleB) : \
				handleA_(handleA), handleB_(handleB) {                         \
			}                                                                              \
			Type access(int64_t index) const {                                             \
				return func(handleA_[index], handleB_[index]);                              \
			}                                                                              \
		private:                                                                          \
			HandleA handleA_;                                                              \
			HandleB handleB_;                                                              \
	};                                                                                   \
		return ThisExpression(objectA.getHandle(), objectB.getHandle());  \
	}                                                                                    \
	friend auto name(class_ const& object, Type const& value) {                          \
		using Handle = decltype(object.getHandle());                                      \
		struct ThisExpression : public Expression<Type, elementCount, ThisExpression> {        \
			ThisExpression(Handle const &handle,  Type const &value) :       \
				handle_(handle), value_(value) {                               \
			}                                                                              \
			Type access(int64_t index) const {                                             \
				return func(handle_[index], value_);                                        \
			}                                                                              \
		private:                                                                          \
			Handle handle_;                                                                \
			Type value_;                                                                   \
		};                                                                                \
		return ThisExpression(object.getHandle(), value);                  \
	}                                                                                    \
	friend auto name(Type const& value, class_ const& object) {                          \
		using Handle = decltype(object.getHandle());                                      \
		struct ThisExpression : public Expression<Type, elementCount, ThisExpression> {        \
			ThisExpression(Handle const &handle, Type const &value) :        \
				handle_(handle), value_(value) {                               \
			}                                                                              \
			Type access(int64_t index) const {                                             \
				return func(value_, handle_[index]);                                        \
			}                                                                              \
		private:                                                                          \
			Handle handle_;                                                                \
			Type value_;                                                                   \
		};                                                                                \
		return ThisExpression(object.getHandle(), value);                  \
	}                                                                                    \


#define VAL_ARRAY_COMPOUND_ASSIGNMENT_OPERATOR(OP)                                   \
	ValArray& operator OP##= (ValArray<Type, elementCount> const &array) {                 \
		for (int64_t i  = 0; i < elementCount; i++) {                          \
			data_[i] OP##= array[i];                                                    \
		}                                                                              \
		return *this;                                                                  \
	}                                                                                 \
	template<typename Derived>                                                        \
	ValArray& operator OP##= (Expression<Type, elementCount, Derived> const &expression) { \
		for (int64_t i  = 0; i < elementCount; i++) {                          \
			data_[i] OP##= expression[i];                                               \
		}                                                                              \
		return *this;                                                                  \
	}                                                                                 \
	ValArray& operator OP##= (Type const &value) {                                    \
		for (int64_t i  = 0; i < elementCount; i++) {                          \
			data_[i] OP##= value;                                                       \
		}                                                                              \
		return *this;                                                                  \
	}                                                                                 \

#define SLICE_ARRAY_COMPOUND_ASSIGNMENT_OPERATOR(OP)                                           \
	SliceArray& operator OP##= (auto const &other) {                                            \
		for (int64_t i  = slice_.start(), j = 0; j != slice_.size(); i += slice_.stride(), j++) { \
			(*pointer_)[i] OP##= other[j];                                                        \
		}                                                                                        \
		return *this;                                                                            \
	}                                                                                           \
	SliceArray& operator OP##= (Type const &value) {                                            \
		for (int64_t i  = slice_.start(), j = 0; j != slice_.size(); i += slice_.stride(), j++) { \
			(*pointer_)[i] OP##= value;                                                           \
		}                                                                                        \
		return *this;                                                                            \
	}                                                                                           \

#define GSLICE_ARRAY_COMPOUND_ASSIGNMENT_OPERATOR(OP)                                                   \
	GSliceArray& operator OP##= (auto const &other) {                                                    \
		int64_t const start = gSlice_.start();                                                            \
		auto const &strides = gSlice_.strides();                                                          \
		auto const &sizes = gSlice_.sizes();                                                              \
		ValArray<int64_t, elementCount> indices;                                              \
		indices.fill(0);\
		int64_t i  = start;                                                                                \
		(*pointer_)[start] OP##= other[0];                                                                \
		for (int64_t j = 1; j < elementCount; j++) {                                                         \
			int dimension = dimensionCount - 1;                                                          \
			while (++indices[dimension] == sizes[dimension]) {                                             \
				indices[dimension] = 0;                                                                     \
				i -= (sizes[dimension] - 1) * strides[dimension];                                           \
				dimension--;                                                                                \
			}                                                                                              \
			i += strides[dimension];                                                                       \
			(*pointer_)[i] OP##= other[j];                                                                 \
		}                                                                                                 \
		return *this;                                                                                     \
	}                                                                                                    \
	GSliceArray& operator OP##= (Type const &value) {                                                    \
		int64_t const start = gSlice_.start();                                                            \
		auto const &strides = gSlice_.strides();                                                          \
		auto const &sizes = gSlice_.sizes();                                                              \
		ValArray<int64_t, elementCount> indices(0, sizes.size());                                              \
		indices.fill(0);\
		int64_t i  = start;                                                                                \
		(*pointer_)[start] OP##= value;                                                                   \
		for (int64_t j = 1; j < elementCount; j++) {                                                         \
			int64_t dimension = dimensionCount - 1;                                                          \
			while (++indices[dimension] == sizes[dimension]) {                                             \
				indices[dimension] = 0;                                                                     \
				i -= (sizes[dimension] - 1) * strides[dimension];                                           \
					dimension--;                                                                                \
			}                                                                                              \
			i += strides[dimension];                                                                       \
			(*pointer_)[i] OP##= value;                                                                    \
		}                                                                                                 \
		return *this;                                                                                     \
	}

#define MASK_ARRAY_COMPOUND_ASSIGNMENT_OPERATOR(OP) \
	MaskArray& operator OP##= (auto const &other) {  \
		int64_t i  = 0;                                \
		for (int64_t j = 0; j < elementCount; j++) { \
			if(mask_[j]) {                             \
				(*pointer_)[i++] OP##= other[j];        \
			}                                          \
		}                                             \
		return *this;                                 \
}                                                   \
	MaskArray& operator OP##= (Type const &value) {  \
		int64_t i  = 0;                                \
		for (int64_t j = 0; j < elementCount; j++) { \
			if(mask_[j]) {                             \
				(*pointer_)[i++] OP##= value;           \
			}                                          \
		}                                             \
		return *this;                                 \
}                                                   \


#define INDIRECT_ARRAY_COMPOUND_ASSIGNMENT_OPERATOR(OP) \
	IndirectArray& operator OP##= (auto const &other){   \
		for (int64_t i  = 0; i < elementCount; i++) {  \
			(*pointer_)[indices_[i]] OP##= other[i];       \
		}                                                 \
		return *this;                                     \
	}                                                    \
	IndirectArray& operator OP##= (Type const &value){   \
		for (int64_t i  = 0; i < elementCount; i++) {  \
			(*pointer_)[indices_[i]] OP##= value;          \
		}                                                 \
		return *this;                                     \
	}                                                    \

template<typename, int64_t, typename >
struct Expression;

template<typename, int64_t>
struct ValArray;

template<typename >
struct UnaryPlus;

template<typename >
struct ShiftLeft;

template<typename >
struct ShiftRight;

struct Slice;

template<int>
class GSlice;

template<typename, int64_t>
struct SliceArray;

template<typename, int64_t, int>
struct GSliceArray;

#include <cstddef>

class Slice {
	int64_t start_;
	int64_t size_;
	int64_t stride_;
public:
	Slice();
	Slice(int64_t start, int64_t size, int64_t stride) :
			start_(start), size_(size), stride_(stride) {
	}
	Slice(Slice const &other) :
			start_(other.start_), size_(other.size_), stride_(other.stride_) {
	}
	int64_t start() const {
		return start_;
	}
	int64_t size() const {
		return size_;
	}
	int64_t stride() const {
		return stride_;
	}
};

template<int dimensionCount>
class GSlice {
	int64_t start_;
	ValArray<int64_t, dimensionCount> sizes_;
	ValArray<int64_t, dimensionCount> strides_;
public:
	GSlice() {
	}
	GSlice(int64_t start, ValArray<int64_t, dimensionCount> sizes, ValArray<int64_t, dimensionCount> strides) :
			start_(start), sizes_(sizes), strides_(strides) {
	}
	GSlice(GSlice const &other) :
			start_(other.start_), sizes_(other.sizes_), strides_(other.strides_) {
	}
	int64_t start() const {
		return start_;
	}
	ValArray<int64_t, dimensionCount> const& sizes() const {
		return sizes_;
	}
	ValArray<int64_t, dimensionCount> const& strides() const {
		return strides_;
	}
};

template<typename Type, int64_t elementCount>
class SliceArray {
	friend class ValArray<Type, elementCount> ;
	ValArray<Type, elementCount> *pointer_;
	Slice slice_;
	SliceArray(ValArray<Type, elementCount> *pointer, Slice slice) :
			pointer_(pointer), slice_(slice) {
	}
public:
	SliceArray(SliceArray const &other) :
			pointer_(other.pointer_), slice_(other.slice_) {
	}
	SliceArray& operator=(auto const &other) {
		for (int64_t i = slice_.start(), j = 0; j != slice_.size(); i += slice_.stride(), j++) {
			(*pointer_)[i] = other[j];
		}
		return *this;
	}
	SliceArray& operator=(Type const &value) {
		for (int64_t i = slice_.start(), j = 0; j != slice_.size(); i += slice_.stride(), j++) {
			(*pointer_)[i] = value;
		}
		return *this;
	}
	SLICE_ARRAY_COMPOUND_ASSIGNMENT_OPERATOR(+)
	SLICE_ARRAY_COMPOUND_ASSIGNMENT_OPERATOR(-)
	SLICE_ARRAY_COMPOUND_ASSIGNMENT_OPERATOR(*)
	SLICE_ARRAY_COMPOUND_ASSIGNMENT_OPERATOR(/)
	SLICE_ARRAY_COMPOUND_ASSIGNMENT_OPERATOR(%)
	SLICE_ARRAY_COMPOUND_ASSIGNMENT_OPERATOR(&)
	SLICE_ARRAY_COMPOUND_ASSIGNMENT_OPERATOR(|)
	SLICE_ARRAY_COMPOUND_ASSIGNMENT_OPERATOR(^)
	SLICE_ARRAY_COMPOUND_ASSIGNMENT_OPERATOR(<<)
	SLICE_ARRAY_COMPOUND_ASSIGNMENT_OPERATOR(>>)
};

template<typename Type, int64_t elementCount, int dimensionCount>
class GSliceArray {
	friend class ValArray<Type, elementCount> ;
	ValArray<Type, elementCount> *pointer_;
	GSlice<dimensionCount> gSlice_;
	GSliceArray(ValArray<Type, elementCount> *pointer, GSlice<dimensionCount> const &slice) :
			pointer_(pointer), gSlice_(slice) {
	}
public:
	GSliceArray(GSliceArray const &other) :
			pointer_(other.pointer_), gSlice_(other.gSlice_) {
	}
	GSliceArray& operator=(auto const &other) {
		constexpr std::multiplies<int64_t> integerMultiply { };
		int64_t const start = gSlice_.start();
		auto const &strides = gSlice_.strides();
		auto const &sizes = gSlice_.sizes();
		int64_t const totalSize = std::accumulate(sizes.begin(), sizes.end(), int64_t(1), integerMultiply);
		ValArray<int64_t, elementCount> indices(0, sizes.size());
		int64_t i = start;
		(*pointer_)[start] = other[0];
		for (int64_t j = 1; j < totalSize; j++) {
			int dimension = dimensionCount - 1;
			while (++indices[dimension] == sizes[dimension]) {
				indices[dimension] = 0;
				i -= (sizes[dimension] - 1) * strides[dimension];
				dimension--;
			}
			i += strides[dimension];
			(*pointer_)[i] = other[j];
		}
		return *this;
	}
	GSliceArray& operator=(Type const &value) {
		constexpr std::multiplies<int64_t> integerMultiply { };
		int64_t const start = gSlice_.start();
		auto const &strides = gSlice_.strides();
		auto const &sizes = gSlice_.sizes();
		int64_t const totalSize = std::accumulate(sizes.begin(), sizes.end(), int64_t(1), integerMultiply);
		ValArray<int64_t, dimensionCount> indices;
		indices.fill(0);
		int64_t i = start;
		(*pointer_)[start] = value;
		for (int64_t j = 1; j < totalSize; j++) {
			int dimension = dimensionCount - 1;
			while (++indices[dimension] == sizes[dimension]) {
				indices[dimension] = 0;
				i -= (sizes[dimension] - 1) * strides[dimension];
				dimension--;
			}
			i += strides[dimension];
			(*pointer_)[i] = value;
		}
		return *this;
	}
	GSLICE_ARRAY_COMPOUND_ASSIGNMENT_OPERATOR(+)
	GSLICE_ARRAY_COMPOUND_ASSIGNMENT_OPERATOR(-)
	GSLICE_ARRAY_COMPOUND_ASSIGNMENT_OPERATOR(*)
	GSLICE_ARRAY_COMPOUND_ASSIGNMENT_OPERATOR(/)
	GSLICE_ARRAY_COMPOUND_ASSIGNMENT_OPERATOR(%)
	GSLICE_ARRAY_COMPOUND_ASSIGNMENT_OPERATOR(&)
	GSLICE_ARRAY_COMPOUND_ASSIGNMENT_OPERATOR(|)
	GSLICE_ARRAY_COMPOUND_ASSIGNMENT_OPERATOR(^)
	GSLICE_ARRAY_COMPOUND_ASSIGNMENT_OPERATOR(<<)
	GSLICE_ARRAY_COMPOUND_ASSIGNMENT_OPERATOR(>>)
};

template<typename Type, int64_t elementCount>
class MaskArray {
	friend class ValArray<Type, elementCount> ;
	ValArray<Type, elementCount> *pointer_;
	ValArray<bool, elementCount> mask_;
	MaskArray(ValArray<Type, elementCount> *pointer, ValArray<bool, elementCount> const &mask) :
			pointer_(pointer), mask_(mask) {
	}
public:
	int64_t sum() const {
		int64_t result = 0;
		for (int64_t j = 0; j < mask_.size(); j++) {
			result = mask_[j] ? (result + 1) : result;
		}
		return result;
	}
	MaskArray(MaskArray const &other) :
			pointer_(other.pointer_), mask_(other.mask_) {
	}
	MaskArray& operator=(auto const &other) {
		int64_t i = 0;
		for (int64_t j = 0; j < mask_.size(); j++) {
			if (mask_[j]) {
				(*pointer_)[i++] = other[j];
			}
		}
		return *this;
	}
	MaskArray& operator=(Type const &value) {
		int64_t i = 0;
		for (int64_t j = 0; j < elementCount; j++) {
			if (mask_[j]) {
				(*pointer_)[i++] = value;
			}
		}
		return *this;
	}
	MASK_ARRAY_COMPOUND_ASSIGNMENT_OPERATOR(+)
	MASK_ARRAY_COMPOUND_ASSIGNMENT_OPERATOR(-)
	MASK_ARRAY_COMPOUND_ASSIGNMENT_OPERATOR(*)
	MASK_ARRAY_COMPOUND_ASSIGNMENT_OPERATOR(/)
	MASK_ARRAY_COMPOUND_ASSIGNMENT_OPERATOR(%)
	MASK_ARRAY_COMPOUND_ASSIGNMENT_OPERATOR(&)
	MASK_ARRAY_COMPOUND_ASSIGNMENT_OPERATOR(|)
	MASK_ARRAY_COMPOUND_ASSIGNMENT_OPERATOR(^)
	MASK_ARRAY_COMPOUND_ASSIGNMENT_OPERATOR(<<)
	MASK_ARRAY_COMPOUND_ASSIGNMENT_OPERATOR(>>)
};

template<typename Type, int64_t elementCount>
class IndirectArray {
	friend class ValArray<Type, elementCount> ;
	ValArray<Type, elementCount> *pointer_;
	ValArray<int64_t, elementCount> indices_;
	IndirectArray(ValArray<Type, elementCount> *pointer, ValArray<int64_t, elementCount> const &indices) :
			pointer_(pointer), indices_(indices) {
	}
public:
	IndirectArray(IndirectArray const &other) :
			pointer_(other.pointer_), indices_(other.indices_) {
	}
	IndirectArray& operator=(auto const &other) {
		for (int64_t i = 0; i < elementCount; i++) {
			(*pointer_)[indices_[i]] = other[i];
		}
		return *this;
	}
	IndirectArray& operator=(Type const &value) {
		for (int64_t i = 0; i < elementCount; i++) {
			(*pointer_)[indices_[i]] = value;
		}
		return *this;
	}
	INDIRECT_ARRAY_COMPOUND_ASSIGNMENT_OPERATOR(+)
	INDIRECT_ARRAY_COMPOUND_ASSIGNMENT_OPERATOR(-)
	INDIRECT_ARRAY_COMPOUND_ASSIGNMENT_OPERATOR(*)
	INDIRECT_ARRAY_COMPOUND_ASSIGNMENT_OPERATOR(/)
	INDIRECT_ARRAY_COMPOUND_ASSIGNMENT_OPERATOR(%)
	INDIRECT_ARRAY_COMPOUND_ASSIGNMENT_OPERATOR(&)
	INDIRECT_ARRAY_COMPOUND_ASSIGNMENT_OPERATOR(|)
	INDIRECT_ARRAY_COMPOUND_ASSIGNMENT_OPERATOR(^)
	INDIRECT_ARRAY_COMPOUND_ASSIGNMENT_OPERATOR(<<)
	INDIRECT_ARRAY_COMPOUND_ASSIGNMENT_OPERATOR(>>)
};

template<typename Type, int64_t elementCount>
struct ValArray {
	using value_type = Type;
	ValArray() {
		data_.resize(elementCount);
	}
	ValArray(const Type &value) {
		data_.resize(elementCount, value);
	}
	ValArray(const Type *valuesPointer) {
		data_.resize(elementCount);
		std::copy(valuesPointer, valuesPointer + elementCount, data_.begin());
	}
	ValArray(const ValArray &other) {
		*this = other;
	}
	ValArray(ValArray &&other) {
		*this = std::move(other);
	}
	ValArray(std::initializer_list<Type> valuesList) {
		data_.resize(elementCount);
		std::copy(valuesList.begin(), valuesList.end(), data_.begin());
	}
	ValArray(SliceArray<Type, elementCount> const &sliceArray) {
		data_.resize(elementCount);
		*this = sliceArray;
	}
	template<int dimensionCount>
	ValArray(GSliceArray<Type, elementCount, dimensionCount> const &gSliceArray) {
		data_.resize(elementCount);
		*this = gSliceArray;
	}
	ValArray(MaskArray<Type, elementCount> const &maskArray) {
		data_.resize(elementCount);
		*this = maskArray;
	}
	ValArray(IndirectArray<Type, elementCount> const &indirectArray) {
		data_.resize(elementCount);
		*this = indirectArray;
	}
	template<typename Derived>
	ValArray(Expression<Type, elementCount, Derived> const &expression) {
		data_.resize(elementCount);
		*this = expression;
	}
	virtual ~ValArray() {
	}
	ValArray& operator=(const Type &value) {
		std::fill(data_.begin(), data_.end(), value);
		return *this;
	}
	ValArray& operator=(const ValArray &other) {
		data_ = other.data_;
		return *this;
	}
	ValArray& operator=(ValArray &&other) {
		swap(std::move(other));
		return *this;
	}
	ValArray& operator=(std::initializer_list<Type> valuesList) {
		*this = ValArray { std::move(valuesList) };
		return *this;
	}
	ValArray& operator=(SliceArray<Type, elementCount> const &sliceArray) {
		auto const &slice = sliceArray.slice_;
		auto const &other = *(sliceArray.pointer_);
		for (int64_t i = 0, j = slice.start(); i != slice.size(); i++, j += slice.stride()) {
			data_[i] = other[j];
		}
		return *this;
	}
	template<int dimensionCount>
	ValArray& operator=(GSliceArray<Type, elementCount, dimensionCount> const &gSliceArray) {
		int64_t const start = gSliceArray.gSlice_.start();
		auto const &strides = gSliceArray.gSlice_.strides();
		auto const &sizes = gSliceArray.gSlice_.sizes();
		auto const &other = *(gSliceArray.pointer_);
		ValArray<int64_t, dimensionCount> indices;
		indices.fill(0);
		int64_t j = start;
		data_[0] = other[start];
		for (int64_t i = 1; i < elementCount; i++) {
			int dimension = dimensionCount - 1;
			while (++indices[dimension] == sizes[dimension]) {
				indices[dimension] = 0;
				j -= (sizes[dimension] - 1) * strides[dimension];
				dimension--;
			}
			j += strides[dimension];
			data_[i] = other[j];
		}
		return *this;
	}
	ValArray& operator=(MaskArray<Type, elementCount> const &maskArray) {
		int64_t k = 0;
		auto const &other = *(maskArray.pointer_);
		for (int64_t i = 0; i < maskArray.mask_.size(); i++) {
			if (maskArray.mask_[i]) {
				data_[i] = other[k++];
			}
		}
		return *this;
	}
	ValArray& operator=(IndirectArray<Type, elementCount> const &indirectArray) {
		auto const &other = *(indirectArray.pointer_);
		auto const &indices = indirectArray.indices_;
		for (int64_t i = 0; i < indices.size(); i++) {
			data_[indices[i]] = other[i];
		}
		return *this;
	}
	template<typename Derived>
	ValArray& operator=(Expression<Type, elementCount, Derived> const &expression) {
		for (int64_t i = 0; i < elementCount; i++) {
			data_[i] = expression[i];
		}
		return *this;
	}
	Type& operator[](int64_t index) {
		return data_[index];
	}
	Type const& operator[](int64_t index) const {
		return data_[index];
	}
	SliceArray<Type, elementCount> operator[](Slice slice) {
		return SliceArray<Type, elementCount> { this, slice };
	}
	template<int dimensionCount>
	GSliceArray<Type, elementCount, dimensionCount> operator[](GSlice<dimensionCount> const &gSlice) {
		return GSliceArray<Type, elementCount, dimensionCount> { this, std::move(gSlice) };
	}
	MaskArray<Type, elementCount> operator[](ValArray<bool, elementCount> const &mask) {
		return MaskArray<Type, elementCount> { this, std::move(mask) };
	}
	IndirectArray<Type, elementCount> operator[](ValArray<int64_t, elementCount> const &indirect) {
		return IndirectArray<Type, elementCount> { this, std::move(indirect) };
	}
	auto operator[](Slice slice) const {
		return const_cast<ValArray*>(this)->operator[](slice);
	}
	template<int dimensionCount>
	auto operator[](GSlice<dimensionCount> const &gSlice) const {
		return const_cast<ValArray*>(this)->operator[](gSlice);
	}
	auto operator[](ValArray<bool, elementCount> const &mask) const {
		return const_cast<ValArray*>(this)->operator[](mask);
	}
	auto operator[](ValArray<int64_t, elementCount> const &indirect) const {
		return const_cast<ValArray*>(this)->operator[](indirect);
	}
	void swap(ValArray<Type, elementCount> &&other) {
		std::swap(data_, other.data_);
	}
	Type sum() const {
		Type result = data_[0];
		for (int64_t i = 1; i < elementCount; i++) {
			result += data_[i];
		}
		return result;
	}
	Type min() const {
		Type result = data_[0];
		for (int64_t i = 1; i < elementCount; i++) {
			Type const &thisElement = data_[i];
			if (thisElement < result) {
				result = thisElement;
			}
		}
		return result;
	}
	Type max() const {
		Type result = data_[0];
		for (int64_t i = 1; i < elementCount; i++) {
			Type const &thisElement = data_[i];
			if (result < thisElement) {
				result = thisElement;
			}
		}
		return result;
	}
	ValArray shift(int thisShift) const {
		ValArray result(elementCount);
		if (thisShift > 0) {
			std::copy_n(data_.begin() + thisShift, elementCount - thisShift, result.begin());
		} else if (thisShift < 0) {
			std::copy_n(data_.begin(), elementCount + thisShift, result.begin() - thisShift);
		}
		return result;
	}
	ValArray cshift(int thisShift) const {
		ValArray result(elementCount);
		if (thisShift > 0) {
			std::copy_n(data_.begin() + thisShift, elementCount - thisShift, result.begin());
			std::copy_n(data_.begin(), thisShift, result.begin() + elementCount - thisShift);
		} else if (thisShift < 0) {
			std::copy_n(data_.begin(), elementCount + thisShift, result.begin() - thisShift);
			std::copy_n(data_.begin() + elementCount + thisShift, -thisShift, result.begin());
		}
		return result;
	}
	ValArray apply(Type thisFunction(Type)) const {
		ValArray result(elementCount);
		for (int64_t i = 0; i < elementCount; i++) {
			result.data_[i] = thisFunction(data_);
		}
		return result;
	}
	ValArray apply(Type thisFunction(const Type&)) const {
		ValArray result(elementCount);
		for (int64_t i = 0; i < elementCount; i++) {
			result.data_[i] = thisFunction(data_);
		}
		return result;
	}
	VAL_ARRAY_UNARY_OPERATOR(+)
	VAL_ARRAY_UNARY_OPERATOR(-)
	VAL_ARRAY_UNARY_OPERATOR(!)
	VAL_ARRAY_UNARY_OPERATOR(~)
	VAL_ARRAY_BINARY_OPERATOR1(+)
	VAL_ARRAY_BINARY_OPERATOR1(-)
	VAL_ARRAY_BINARY_OPERATOR1(*)
	VAL_ARRAY_BINARY_OPERATOR1(/)
	VAL_ARRAY_BINARY_OPERATOR1(%)
	VAL_ARRAY_BINARY_OPERATOR1(&&)
	VAL_ARRAY_BINARY_OPERATOR1(||)
	VAL_ARRAY_BINARY_OPERATOR1(&)
	VAL_ARRAY_BINARY_OPERATOR1(|)
	VAL_ARRAY_BINARY_OPERATOR1(^)
	VAL_ARRAY_BINARY_OPERATOR1(>>)
	VAL_ARRAY_BINARY_OPERATOR1(<<)
	VAL_ARRAY_COMPARE_OPERATOR(<)
	VAL_ARRAY_COMPARE_OPERATOR(>)
	VAL_ARRAY_COMPARE_OPERATOR(<=)
	VAL_ARRAY_COMPARE_OPERATOR(>=)
	VAL_ARRAY_COMPARE_OPERATOR(==)
	VAL_ARRAY_COMPARE_OPERATOR(!=)
	VAL_ARRAY_UNARY_FUNCTION(ValArray, abs)
	VAL_ARRAY_UNARY_FUNCTION(ValArray, exp)
	VAL_ARRAY_UNARY_FUNCTION(ValArray, log)
	VAL_ARRAY_UNARY_FUNCTION(ValArray, log10)
	VAL_ARRAY_UNARY_FUNCTION(ValArray, sqrt)
	VAL_ARRAY_UNARY_FUNCTION(ValArray, asin)
	VAL_ARRAY_UNARY_FUNCTION(ValArray, acos)
	VAL_ARRAY_UNARY_FUNCTION(ValArray, atan)
	VAL_ARRAY_UNARY_FUNCTION(ValArray, sin)
	VAL_ARRAY_UNARY_FUNCTION(ValArray, cos)
	VAL_ARRAY_UNARY_FUNCTION(ValArray, tan)
	VAL_ARRAY_UNARY_FUNCTION(ValArray, sinh)
	VAL_ARRAY_UNARY_FUNCTION(ValArray, cosh)
	VAL_ARRAY_UNARY_FUNCTION(ValArray, tanh)
	VAL_ARRAY_BINARY_FUNCTION(ValArray, atan2, atan2)
	VAL_ARRAY_BINARY_FUNCTION(ValArray, pow, pow)
	VAL_ARRAY_BINARY_FUNCTION(ValArray, copysign, copysign)
	VAL_ARRAY_BINARY_FUNCTION(ValArray, min, std::min)
	VAL_ARRAY_BINARY_FUNCTION(ValArray, max, std::max)
	VAL_ARRAY_COMPOUND_ASSIGNMENT_OPERATOR(+)
	VAL_ARRAY_COMPOUND_ASSIGNMENT_OPERATOR(-)
	VAL_ARRAY_COMPOUND_ASSIGNMENT_OPERATOR(*)
	VAL_ARRAY_COMPOUND_ASSIGNMENT_OPERATOR(/)
	VAL_ARRAY_COMPOUND_ASSIGNMENT_OPERATOR(%)
	VAL_ARRAY_COMPOUND_ASSIGNMENT_OPERATOR(&)
	VAL_ARRAY_COMPOUND_ASSIGNMENT_OPERATOR(|)
	VAL_ARRAY_COMPOUND_ASSIGNMENT_OPERATOR(^)
	VAL_ARRAY_COMPOUND_ASSIGNMENT_OPERATOR(<<)
	VAL_ARRAY_COMPOUND_ASSIGNMENT_OPERATOR(>>)
	ValArray& fill(Type const &value) {
		std::fill_n(data_.begin(), elementCount, value);
		return *this;
	}
	auto begin() {
		return data_.begin();
	}
	auto end() {
		return data_.end();
	}
	auto begin() const {
		return data_.begin();
	}
	auto end() const {
		return data_.end();
	}
	std::vector<Type>::iterator getHandle() {
		return data_.begin();
	}
	std::vector<Type>::const_iterator getHandle() const {
		return data_.begin();
	}
private:
	std::vector<Type> data_;
};

template<typename Type, int64_t elementCount, typename Derived>
struct Handle {
	Handle(Expression<Type, elementCount, Derived> const &reference) :
			reference_(reference) {
	}
	Type operator[](int64_t index) const {
		return reference_[index];
	}
private:
	Expression<Type, elementCount, Derived> const &reference_;
};

template<typename Type, int64_t elementCount, typename Derived>
struct Expression {
	VAL_ARRAY_UNARY_OPERATOR(+)
	VAL_ARRAY_UNARY_OPERATOR(-)
	VAL_ARRAY_UNARY_OPERATOR(!)
	VAL_ARRAY_UNARY_OPERATOR(~)
	VAL_ARRAY_BINARY_OPERATOR1(+)
	VAL_ARRAY_BINARY_OPERATOR1(-)
	VAL_ARRAY_BINARY_OPERATOR1(*)
	VAL_ARRAY_BINARY_OPERATOR1(/)
	VAL_ARRAY_BINARY_OPERATOR1(%)
	VAL_ARRAY_BINARY_OPERATOR1(&&)
	VAL_ARRAY_BINARY_OPERATOR1(||)
	VAL_ARRAY_BINARY_OPERATOR1(&)
	VAL_ARRAY_BINARY_OPERATOR1(|)
	VAL_ARRAY_BINARY_OPERATOR1(^)
	VAL_ARRAY_BINARY_OPERATOR1(>>)
	VAL_ARRAY_BINARY_OPERATOR1(<<)
	VAL_ARRAY_COMPARE_OPERATOR(<)
	VAL_ARRAY_COMPARE_OPERATOR(>)
	VAL_ARRAY_COMPARE_OPERATOR(<=)
	VAL_ARRAY_COMPARE_OPERATOR(>=)
	VAL_ARRAY_COMPARE_OPERATOR(==)
	VAL_ARRAY_COMPARE_OPERATOR(!=)
	VAL_ARRAY_UNARY_FUNCTION(Derived, abs)
	VAL_ARRAY_UNARY_FUNCTION(Derived, exp)
	VAL_ARRAY_UNARY_FUNCTION(Derived, log)
	VAL_ARRAY_UNARY_FUNCTION(Derived, log10)
	VAL_ARRAY_UNARY_FUNCTION(Derived, sqrt)
	VAL_ARRAY_UNARY_FUNCTION(Derived, asin)
	VAL_ARRAY_UNARY_FUNCTION(Derived, acos)
	VAL_ARRAY_UNARY_FUNCTION(Derived, atan)
	VAL_ARRAY_UNARY_FUNCTION(Derived, sin)
	VAL_ARRAY_UNARY_FUNCTION(Derived, cos)
	VAL_ARRAY_UNARY_FUNCTION(Derived, tan)
	VAL_ARRAY_UNARY_FUNCTION(Derived, sinh)
	VAL_ARRAY_UNARY_FUNCTION(Derived, cosh)
	VAL_ARRAY_UNARY_FUNCTION(Derived, tanh)
	VAL_ARRAY_BINARY_FUNCTION(Derived, pow, pow)
	VAL_ARRAY_BINARY_FUNCTION(Derived, copysign, copysign)
	VAL_ARRAY_BINARY_FUNCTION(Derived, atan2, atan2)
	VAL_ARRAY_BINARY_FUNCTION(Derived, min, std::min)
	VAL_ARRAY_BINARY_FUNCTION(Derived, max, std::max)
	Type operator[](int64_t index) const {
		return static_cast<Derived const*>(this)->access(index);
	}
	Handle<Type, elementCount, Derived> getHandle() const {
		return Handle<Type, elementCount, Derived> { *this };
	}
	int64_t size() const {
		return static_cast<Derived const*>(this)->count();
	}
};

template<typename Type>
struct UnaryPlus {
	Type operator()(Type const &value) const {
		return value;
	}
};

template<typename Type>
struct ShiftLeft {
	Type operator()(Type const &value, Type const &shift) const {
		return value << shift;
	}
};

template<typename Type>
struct ShiftRight {
	Type operator()(Type const &value, Type const &shift) const {
		return value >> shift;
	}
};

VAL_ARRAY_BINARY_OPERATOR2(+)
VAL_ARRAY_BINARY_OPERATOR2(-)
VAL_ARRAY_BINARY_OPERATOR2(*)
VAL_ARRAY_BINARY_OPERATOR2(/)
VAL_ARRAY_BINARY_OPERATOR2(%)
VAL_ARRAY_BINARY_OPERATOR2(&&)
VAL_ARRAY_BINARY_OPERATOR2(||)
VAL_ARRAY_BINARY_OPERATOR2(&)
VAL_ARRAY_BINARY_OPERATOR2(|)
VAL_ARRAY_BINARY_OPERATOR2(^)
VAL_ARRAY_BINARY_OPERATOR2(>>)
VAL_ARRAY_BINARY_OPERATOR2(<<)

#endif /* INCLUDE_VAL_ARRAY_HPP_ */
