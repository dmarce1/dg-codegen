#ifndef INCLUDE_VAL_ARRAY_HPP_
#define INCLUDE_VAL_ARRAY_HPP_

#include <cstdlib>
#include <initializer_list>
#include <functional>
#include <memory>
#include <numeric>
#include <stack>
#include <type_traits>
#include <utility>

#define CONCAT_IMPL(x, y) x##y
#define CONCAT(x, y) CONCAT_IMPL(x, y)
#define UNIQUE_NAME(base) CONCAT(base, __LINE__)

#define VAL_ARRAY_UNARY_OPERATOR(op)                                                               \
	auto operator op() const {                                                                      \
		using Handle = decltype(this->getHandle());                                                  \
		struct UNIQUE_NAME(Expression) : public Expression<Type, Integer, UNIQUE_NAME(Expression)> { \
		UNIQUE_NAME(Expression)(Handle const &handle) :                                              \
			handle_(handle) {                                                                         \
		}                                                                                            \
		Type operator[](Integer index) const {                                                       \
			return op handle_[index];                                                                 \
		}                                                                                            \
		private:                                                                                     \
			Handle handle_;                                                                           \
		};                                                                                           \
		return UNIQUE_NAME(Expression)(this->getHandle());                                           \
	}

#define VAL_ARRAY_BINARY_OPERATOR(op)                                                              \
	auto operator op(auto const& other) const {                                                     \
		using HandleA = decltype(this->getHandle());                                                 \
		using HandleB = decltype(other.getHandle());                                                 \
		struct UNIQUE_NAME(Expression) : public Expression<Type, Integer, UNIQUE_NAME(Expression)> { \
		UNIQUE_NAME(Expression)(HandleA const &handleA, HandleB const &handleB) :                    \
			handleA_(handleA), handleB_(handleB) {                                                    \
		}                                                                                            \
		Type operator[](Integer index) const {                                                       \
			return handleA_[index] op handleB_[index];                                                \
		}                                                                                            \
		private:                                                                                     \
			HandleA handleA_;                                                                         \
			HandleB handleB_;                                                                         \
		};                                                                                           \
		return UNIQUE_NAME(Expression)(this->getHandle(), other.getHandle());                       \
	}

#define VAL_ARRAY_COMPARE_OPERATOR(op)                                                             \
	auto operator op(auto const& other) const {                                                     \
		using HandleA = decltype(this->getHandle());                                                 \
		using HandleB = decltype(other.getHandle());                                                 \
		struct UNIQUE_NAME(Expression) : public Expression<bool, Integer, UNIQUE_NAME(Expression)> { \
		UNIQUE_NAME(Expression)(HandleA const &handleA, HandleB const &handleB) :                    \
			handleA_(handleA), handleB_(handleB) {                                                    \
		}                                                                                            \
		bool operator[](Integer index) const {                                                       \
			return handleA_[index] op handleB_[index];                                                \
		}                                                                                            \
		private:                                                                                     \
			HandleA handleA_;                                                                         \
			HandleB handleB_;                                                                         \
		};                                                                                           \
		return UNIQUE_NAME(Expression)(this->getHandle(), other->getHandle());                       \
	}

#define VAL_ARRAY_UNARY_FUNCTION(class_, name)                                       \
	friend auto name(class_ const& object) {                                          \
		using Handle = decltype(object.getHandle());                                   \
		struct Expression##name : public Expression<Type, Integer, Expression##name> { \
		Expression##name(Handle const &handle) :                                       \
			handle_(handle) {                                                           \
		}                                                                              \
		Type operator[](Integer index) const {                                         \
			using namespace std;                                                        \
			return name(handle_[index]);                                                \
		}                                                                              \
		private:                                                                       \
			Handle handle_;                                                             \
		};                                                                             \
		return Expression##name(object.getHandle());                                   \
	}

#define VAL_ARRAY_BINARY_FUNCTION(class_, name)                                        \
	friend auto name(class_ const& objectA, auto const& objectB) {                      \
		using HandleA = decltype(objectA.getHandle());                                   \
		using HandleB = decltype(objectB.getHandle());                                   \
		struct Expression##name : public Expression<Type, Integer, Expression##name> {   \
		Expression##name(HandleA const &handleA, HandleB const &handleB) :               \
			handleA_(handleA), handleB_(handleB) {                                        \
		}                                                                                \
		Type operator[](Integer index) const {                                           \
			using namespace std;                                                          \
			return name(handleA_[index], handleB_[index]);                                \
		}                                                                                \
		private:                                                                         \
			HandleA handleA_;                                                             \
			HandleB handleB_;                                                             \
		};                                                                               \
		return Expression##name(objectA.getHandle(), objectB.getHandle());               \
	}                                                                                   \
	friend auto name(class_ const& object, Type const& value) {                         \
		using Handle = decltype(object.getHandle());                                     \
		struct Expression##name : public Expression<Type, Integer, Expression##name> {   \
		Expression##name(Handle const &handle,  Type const &value) :                     \
			handle_(handle), value_(value) {                                              \
		}                                                                                \
		Type operator[](Integer index) const {                                           \
			using namespace std;                                                          \
			return name(handle_[index], value_);                                          \
		}                                                                                \
		private:                                                                         \
			Handle handle_;                                                               \
			Type value_;                                                                  \
		};                                                                               \
		return Expression##name(object.getHandle(), value);                              \
	}                                                                                   \
	friend auto name(Type const& value, class_ const& object) {                         \
		using Handle = decltype(object.getHandle());                                     \
		struct Expression##name : public Expression<Type, Integer, Expression##name> {   \
		Expression##name(Handle const &handle,  Type const &value) :                     \
			handle_(handle), value_(value) {                                              \
		}                                                                                \
		Type operator[](Integer index) const {                                           \
			using namespace std;                                                          \
			return name(value_, handle_[index]);                                          \
		}                                                                                \
		private:                                                                         \
			Handle handle_;                                                               \
			Type value_;                                                                  \
		};                                                                               \
		return Expression##name(object.getHandle(), value);                              \
	}                                                                                   \

#define VAL_ARRAY_COMPOUND_ASSIGNMENT_OPERATOR(OP)                                   \
	template<typename Derived>                                                        \
	ValArray& operator OP##= (Expression<Type, Integer, Derived> const &expression) { \
		for (Integer i = 0; i < Integer(data_.size()); i++) {                          \
			data_[i] OP##= expression[i];                                               \
		}                                                                              \
		return *this;                                                                  \
	}                                                                                 \
	ValArray& operator OP##= (Type const &value) {                                    \
		for (Integer i = 0; i < Integer(data_.size()); i++) {                          \
			data_[i] OP##= value;                                                       \
		}                                                                              \
		return *this;                                                                  \
	}                                                                                 \

#define SLICE_ARRAY_COMPOUND_ASSIGNMENT_OPERATOR(OP)                                           \
	SliceArray& operator OP##= (auto const &other) {                                            \
		for (Integer i = slice_.start(), j = 0; j != slice_.size(); i += slice_.stride(), j++) { \
			(*pointer_)[i] OP##= other[j];                                                        \
		}                                                                                        \
		return *this;                                                                            \
	}                                                                                           \
	SliceArray& operator OP##= (Type const &value) {                                            \
		for (Integer i = slice_.start(), j = 0; j != slice_.size(); i += slice_.stride(), j++) { \
			(*pointer_)[i] OP##= value;                                                           \
		}                                                                                        \
		return *this;                                                                            \
	}                                                                                           \

#define GSLICE_ARRAY_COMPOUND_ASSIGNMENT_OPERATOR(OP)                                                  \
	GSliceArray& operator OP##= (auto const &other) {                                                   \
		constexpr std::multiplies<Integer> integerMultiply{};                                            \
		Integer const start = gSlice_.start();                                                           \
		auto const &strides = gSlice_.strides();                                                         \
		auto const &sizes = gSlice_.sizes();                                                             \
		Integer const totalSize = std::accumulate(begin(sizes), end(sizes), Integer(1), integerMultiply);\
		ValArray<Integer, Integer> indices(0, sizes.size());                                             \
		Integer i = start;                                                                               \
		(*pointer_)[start] OP##= other[0];                                                               \
		for (Integer j = 1; j < totalSize; j++) {                                                        \
			Integer dimension = sizes.size();                                                             \
			while (++indices[--dimension] == sizes[dimension]) {                                          \
				i += strides[dimension];                                                                   \
				indices[dimension] = 0;                                                                    \
				i -= sizes[dimension] * strides[dimension];                                                \
			}                                                                                             \
			i += strides[dimension];                                                                      \
			(*pointer_)[i] OP##= other[j];                                                                \
		}                                                                                                \
		return *this;                                                                                    \
	}                                                                                                   \
	GSliceArray& operator OP##= (Type const &value) {                                                   \
		constexpr std::multiplies<Integer> integerMultiply{};                                            \
		Integer const start = gSlice_.start();                                                           \
		auto const &strides = gSlice_.strides();                                                         \
		auto const &sizes = gSlice_.sizes();                                                             \
		Integer const totalSize = std::accumulate(begin(sizes), end(sizes), Integer(1), integerMultiply);\
		ValArray<Integer, Integer> indices(0, sizes.size());                                             \
		Integer i = start;                                                                               \
		(*pointer_)[start] OP##= value;                                                                  \
		for (Integer j = 1; j < totalSize; j++) {                                                        \
			Integer dimension = sizes.size();                                                             \
			while (++indices[--dimension] == sizes[dimension]) {                                          \
				i += strides[dimension];                                                                   \
				indices[dimension] = 0;                                                                    \
				i -= sizes[dimension] * strides[dimension];                                                \
			}                                                                                             \
			i += strides[dimension];                                                                      \
			(*pointer_)[i] OP##= value;                                                                   \
		}                                                                                                \
		return *this;                                                                                    \
	}

#define MASK_ARRAY_COMPOUND_ASSIGNMENT_OPERATOR(OP)\
	MaskArray& operator OP##= (auto const &other) { \
		Integer i = 0;                               \
		for (Integer j = 0; j < mask_->size(); j++) {\
			if(mask_[j]) {                            \
				(*pointer_)[i++] OP##= other[j];       \
			}                                         \
		}                                            \
		return *this;                                \
}                                                  \
	MaskArray& operator OP##= (Type const &value) { \
		Integer i = 0;                               \
		for (Integer j = 0; j < mask_->size(); j++) {\
			if(mask_[j]) {                            \
				(*pointer_)[i++] OP##= value;          \
			}                                         \
		}                                            \
		return *this;                                \
}                                                  \


#define INDIRECT_ARRAY_COMPOUND_ASSIGNMENT_OPERATOR(OP)\
	IndirectArray& operator OP##= (auto const &other){  \
		for (Integer i = 0; i < indices_->size(); i++) { \
			(*pointer_)[indices_[i]] OP##= other[i];      \
		}                                                \
		return *this;                                    \
	}                                                   \
	IndirectArray& operator OP##= (Type const &value){  \
		for (Integer i = 0; i < indices_->size(); i++) { \
			(*pointer_)[indices_[i]] OP##= value;         \
		}                                                \
		return *this;                                    \
	}                                                   \

template<typename, typename Integer, typename Derived>
struct Expression;

template<typename, typename, typename >
struct UnaryExpression;

template<typename, typename, typename, typename >
struct ComparisonExpression;

template<typename, typename, typename, typename >
struct BinaryExpression;

template<typename, typename = int>
struct ValArray;

template<typename >
struct UnaryPlus;

template<typename >
struct ShiftLeft;

template<typename >
struct ShiftRight;

template<typename = int>
struct Slice;

template<typename >
class GSlice;

template<typename, typename = int>
struct SliceArray;

template<typename, typename = int>
struct GSliceArray;

#include <cstddef>

template<typename Integer>
class Slice {
	Integer start_;
	Integer size_;
	Integer stride_;
public:
	Slice();
	Slice(Integer start, Integer size, Integer stride) :
			start_(start), size_(size), stride_(stride) {
	}
	Slice(Slice const &other) :
			start_(other.start_), size_(other.size_), stride_(other.stride_) {
	}
	Integer start() const {
		return start_;
	}
	Integer size() const {
		return size_;
	}
	Integer stride() const {
		return stride_;
	}
};

template<typename Integer>
class GSlice {
	Integer start_;
	ValArray<Integer, Integer> sizes_;
	ValArray<Integer, Integer> strides_;
public:
	GSlice();
	GSlice(Integer start, ValArray<Integer, Integer> sizes, ValArray<Integer, Integer> strides) :
			start_(start), sizes_(sizes), strides_(strides) {
	}
	GSlice(GSlice const &other) :
			start_(other.start_), sizes_(other.sizes_), strides_(other.strides_) {
	}
	Integer start() const {
		return start_;
	}
	ValArray<Integer, Integer> const& sizes() const {
		return sizes_;
	}
	ValArray<Integer, Integer> const& strides() const {
		return strides_;
	}
};

template<typename Type, typename Integer>
class SliceArray {
	friend class ValArray<Type, Integer> ;
	ValArray<Type, Integer> *pointer_;
	Slice<Integer> slice_;
	SliceArray(ValArray<Type, Integer> *pointer, Slice<Integer> slice) :
			pointer_(pointer), slice_(slice) {
	}
public:
	SliceArray(SliceArray const &other) :
			pointer_(other.pointer_), slice_(other.slice_) {
	}
	SliceArray& operator=(auto const &other) {
		for (Integer i = slice_.start(), j = 0; j != slice_.size(); i += slice_.stride(), j++) {
			(*pointer_)[i] = other[j];
		}
		return *this;
	}
	SliceArray& operator=(Type const &value) {
		for (Integer i = slice_.start(), j = 0; j != slice_.size(); i += slice_.stride(), j++) {
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

template<typename Type, typename Integer>
class GSliceArray {
	friend class ValArray<Type, Integer> ;
	ValArray<Type, Integer> *pointer_;
	GSlice<Integer> gSlice_;
	GSliceArray(ValArray<Type, Integer> *pointer, GSlice<Integer> const &slice) :
			pointer_(pointer), gSlice_(slice) {
	}
public:
	GSliceArray(GSliceArray const &other) :
			pointer_(other.pointer_), gSlice_(other.gSlice_) {
	}
	GSliceArray& operator=(auto const &other) {
		constexpr std::multiplies<Integer> integerMultiply { };
		Integer const start = gSlice_.start();
		auto const &strides = gSlice_.strides();
		auto const &sizes = gSlice_.sizes();
		Integer const totalSize = std::accumulate(sizes.begin(), sizes.end(), Integer(1), integerMultiply);
		ValArray<Integer, Integer> indices(0, sizes.size());
		Integer i = start;
		(*pointer_)[start] = other[0];
		for (Integer j = 1; j < totalSize; j++) {
			Integer dimension = sizes.size() - 1;
			while (++indices[dimension] == sizes[dimension]) {
				dimension--;
				indices[dimension] = 0;
				i -= (sizes[dimension] - 1) * strides[dimension];
			}
			i += strides[dimension];
			(*pointer_)[i] = other[j];
		}
		return *this;
	}
	GSliceArray& operator=(Type const &value) {
		constexpr std::multiplies<Integer> integerMultiply { };
		Integer const start = gSlice_.start();
		auto const &strides = gSlice_.strides();
		auto const &sizes = gSlice_.sizes();
		Integer const totalSize = std::accumulate(sizes.begin(), sizes.end(), Integer(1), integerMultiply);
		ValArray<Integer, Integer> indices(0, sizes.size());
		Integer i = start;
		(*pointer_)[start] = value;
		for (Integer j = 1; j < totalSize; j++) {
			Integer dimension = sizes.size() - 1;
			while (++indices[dimension] == sizes[dimension]) {
				dimension--;
				indices[dimension] = 0;
				i -= (sizes[dimension] - 1) * strides[dimension];
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

template<typename Type, typename Integer>
class MaskArray {
	friend class ValArray<Type, Integer> ;
	ValArray<Type, Integer> *pointer_;
	ValArray<bool, Integer> mask_;
	MaskArray(ValArray<Type, Integer> *pointer, ValArray<bool, Integer> const &mask) :
			pointer_(pointer), mask_(mask) {
	}
public:
	MaskArray(MaskArray const &other) :
			pointer_(other.pointer_), mask_(other.mask_) {
	}
	MaskArray& operator=(auto const &other) {
		Integer i = 0;
		for (Integer j = 0; j < mask_->size(); j++) {
			if (mask_[j]) {
				(*pointer_)[i++] = other[j];
			}
		}
		return *this;
	}
	MaskArray& operator=(Type const &value) {
		Integer i = 0;
		for (Integer j = 0; j < mask_->size(); j++) {
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

template<typename Type, typename Integer>
class IndirectArray {
	friend class ValArray<Type, Integer> ;
	ValArray<Type, Integer> *pointer_;
	ValArray<Integer, Integer> indices_;
	IndirectArray(ValArray<Type, Integer> *pointer, ValArray<Integer, Integer> const &indices) :
			pointer_(pointer), indices_(indices) {
	}
public:
	IndirectArray(IndirectArray const &other) :
			pointer_(other.pointer_), indices_(other.indices_) {
	}
	IndirectArray& operator=(auto const &other) {
		for (Integer i = 0; i < indices_->size(); i++) {
			(*pointer_)[indices_[i]] = other[i];
		}
		return *this;
	}
	IndirectArray& operator=(Type const &value) {
		for (Integer i = 0; i < indices_->size(); i++) {
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

template<typename Type, typename Integer>
struct ValArray {
	ValArray() :
			data_() {
	}
	explicit ValArray(Integer count) :
			data_(count) {
	}
	ValArray(const Type &value, Integer count) :
			data_(count, value) {
	}
	ValArray(const Type *valuesPointer, Integer count) :
			data_(valuesPointer, valuesPointer + count) {
	}
	ValArray(const ValArray &other) {
		*this = other;
	}
	ValArray(ValArray &&other) :
			data_() {
		*this = std::move(other);
	}
	ValArray(std::initializer_list<Type> valuesList) :
			data_(valuesList) {
	}
	ValArray(SliceArray<Type, Integer> const &sliceArray) {
		*this = sliceArray;
	}
	ValArray(GSliceArray<Type, Integer> const &gSliceArray) {
		*this = gSliceArray;
	}
	ValArray(MaskArray<Type, Integer> const &maskArray) {
		*this = maskArray;
	}
	ValArray(IndirectArray<Type, Integer> const &indirectArray) {
		*this = indirectArray;
	}
	template<typename Derived>
	ValArray(Expression<Type, Integer, Derived> const &expression) {
		*this = expression;
	}
	virtual ~ValArray() {
	}
	ValArray& operator=(const Type &value) {
		std::fill_n(data_, data_.size(), value);
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
		deallocate(data_);
		*this = ValArray { std::move(valuesList) };
		return *this;
	}
	ValArray& operator=(SliceArray<Type, Integer> const &sliceArray) {
		auto const &slice = sliceArray.slice_;
		auto const &other = *(sliceArray.pointer_);
		data_.resize(slice.size());
		for (Integer i = 0, j = slice.start(); i != slice.size(); i++, j += slice.stride()) {
			data_[i] = other[j];
		}
		return *this;
	}
	ValArray& operator=(GSliceArray<Type, Integer> const &gSliceArray) {
		constexpr std::multiplies<Integer> integerMultiply { };
		Integer const start = gSliceArray.gSlice_.start();
		auto const &strides = gSliceArray.gSlice_.strides();
		auto const &sizes = gSliceArray.gSlice_.sizes();
		auto const &other = *(gSliceArray.pointer_);
		Integer const totalSize = std::accumulate(sizes.begin(), sizes.end(), Integer(1), integerMultiply);
		ValArray<Integer, Integer> indices(0, sizes.size());
		Integer j = start;
		data_.resize(totalSize);
		data_[0] = other[start];
		for (Integer i = 1; i < totalSize; i++) {
			Integer dimension = sizes.size() - 1;
			while (++indices[dimension] == sizes[dimension]) {
				dimension--;
				indices[dimension] = 0;
				j -= (sizes[dimension] - 1) * strides[dimension];
			}
			j += strides[dimension];
			data_[i] = other[j];
		}
		return *this;
	}
	ValArray& operator=(MaskArray<Type, Integer> const &maskArray) {
		Integer k = 0;
		auto const &other = *(maskArray.pointer_);
		data_.resize(maskArray.sum());
		for (Integer i = 0; i < maskArray.mask_->size(); i++) {
			if (maskArray.mask_[i]) {
				data_[i] = other[k++];
			}
		}
		return *this;
	}
	ValArray& operator=(IndirectArray<Type, Integer> const &indirectArray) {
		auto const &other = *(indirectArray.pointer_);
		auto const& indices = indirectArray.indices_;
		data_.resize(indices.size());
		for (Integer i = 0; i < indices.size(); i++) {
			data_[indices[i]] = other[i];
		}
		return *this;
	}
	template<typename Derived>
	ValArray& operator=(Expression<Type, Integer, Derived> const &expression) {
		data_.resize(expression.size());
		for (Integer i = 0; i < Integer(data_.size()); i++) {
			data_[i] = expression[i];
		}
		return *this;
	}
	Type& operator[](Integer index) {
		return data_[index];
	}
	const Type& operator[](Integer index) const {
		return data_[index];
	}
	SliceArray<Type, Integer> operator[](Slice<Integer> slice) {
		return SliceArray<Type, Integer> { this, slice };
	}
	GSliceArray<Type, Integer> operator[](GSlice<Integer> const &gSlice) {
		return GSliceArray<Type, Integer> { this, std::move(gSlice) };
	}
	MaskArray<Type, Integer> operator[](ValArray<bool, Integer> const &mask) {
		return MaskArray<Type, Integer> { this, std::move(mask) };
	}
	IndirectArray<Type, Integer> operator[](ValArray<Integer, Integer> const &indirect) {
		return IndirectArray<Type, Integer> { this, std::move(indirect) };
	}
	auto operator[](Slice<Integer> slice) const {
		return const_cast<ValArray*>(this)->operator[](slice);
	}
	auto operator[](GSlice<Integer> const &gSlice) const {
		return const_cast<ValArray*>(this)->operator[](gSlice);
	}
	auto operator[](ValArray<bool, Integer> const &mask) const {
		return const_cast<ValArray*>(this)->operator[](mask);
	}
	auto operator[](ValArray<Integer, Integer> const &indirect) const {
		return const_cast<ValArray*>(this)->operator[](indirect);
	}
	void swap(ValArray<Type, Integer> &&other) {
		std::swap(data_, other.data_);
	}
	Integer size() const {
		return data_.size();
	}
	void resize(Integer count, Type value = Type { }) {
		*this = ValArray<Type, Integer>(value, count);
	}
	Type sum() const {
		Type result = data_[0];
		for (Integer i = 1; i < data_.size(); i++) {
			result += data_[i];
		}
		return result;
	}
	Type min() const {
		Type result = data_[0];
		for (Integer i = 1; i < data_.size(); i++) {
			Type const &thisElement = data_[i];
			if (thisElement < result) {
				result = thisElement;
			}
		}
		return result;
	}
	Type max() const {
		Type result = data_[0];
		for (Integer i = 1; i < data_.size(); i++) {
			Type const &thisElement = data_[i];
			if (result < thisElement) {
				result = thisElement;
			}
		}
		return result;
	}
	ValArray shift(int thisShift) const {
		ValArray result(data_.size());
		if (thisShift > 0) {
			std::copy_n(data_.begin() + thisShift, data_.size() - thisShift, result.begin());
		} else if (thisShift < 0) {
			std::copy_n(data_.begin(), data_.size() + thisShift, result.begin() - thisShift);
		}
		return result;
	}
	ValArray cshift(int thisShift) const {
		ValArray result(data_.size());
		if (thisShift > 0) {
			std::copy_n(data_.begin() + thisShift, data_.size() - thisShift, result.begin());
			std::copy_n(data_.begin(), thisShift, result.begin() + data_.size() - thisShift);
		} else if (thisShift < 0) {
			std::copy_n(data_.begin(), data_.size() + thisShift, result.begin() - thisShift);
			std::copy_n(data_.begin() + data_.size() + thisShift, -thisShift, result.begin());
		}
		return result;
	}
	ValArray apply(Type thisFunction(Type)) const {
		ValArray result(data_.size());
		for (Integer i = 0; i < data_.size(); i++) {
			result.data_[i] = thisFunction(data_);
		}
		return result;
	}
	ValArray apply(Type thisFunction(const Type&)) const {
		ValArray result(data_.size());
		for (Integer i = 0; i < data_.size(); i++) {
			result.data_[i] = thisFunction(data_);
		}
		return result;
	}
	VAL_ARRAY_UNARY_OPERATOR(+)
	VAL_ARRAY_UNARY_OPERATOR(-)
	VAL_ARRAY_UNARY_OPERATOR(!)
	VAL_ARRAY_UNARY_OPERATOR(~)
	VAL_ARRAY_BINARY_OPERATOR(+)
	VAL_ARRAY_BINARY_OPERATOR(-)
	VAL_ARRAY_BINARY_OPERATOR(*)
	VAL_ARRAY_BINARY_OPERATOR(/)
	VAL_ARRAY_BINARY_OPERATOR(%)
	VAL_ARRAY_BINARY_OPERATOR(&&)
	VAL_ARRAY_BINARY_OPERATOR(||)
	VAL_ARRAY_BINARY_OPERATOR(&)
	VAL_ARRAY_BINARY_OPERATOR(|)
	VAL_ARRAY_BINARY_OPERATOR(^)
	VAL_ARRAY_BINARY_OPERATOR(>>)
	VAL_ARRAY_BINARY_OPERATOR(<<)
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
	VAL_ARRAY_BINARY_FUNCTION(ValArray, pow)
	VAL_ARRAY_BINARY_FUNCTION(ValArray, min)
	VAL_ARRAY_BINARY_FUNCTION(ValArray, max)
	VAL_ARRAY_BINARY_FUNCTION(ValArray, copysign)
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
		std::fill_n(data_.begin(), data_.size(), value);
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
	auto const& getHandle() const {
		return data_;
	}
private:
	template<typename Derived>
	operator Expression<Type, Integer, Derived>() const {
		return operator+();
	}
	static constexpr bool isPod = std::is_standard_layout<Type>::value && std::is_trivial<Type>::value;
	std::vector<Type> data_;
};

template<typename Type>
void swap(ValArray<Type> &&valarrayA, ValArray<Type> &&valarrayB) {
	valarrayA.swap(std::move(valarrayB));
}

template<typename Type>
auto begin(ValArray<Type> &valarray) {
	return valarray.begin();
}

template<typename Type>
auto end(ValArray<Type> &valarray) {
	return valarray.end();
}

template<typename Type>
auto begin(ValArray<Type> const &valarray) {
	return valarray.begin();
}

template<typename Type>
auto end(ValArray<Type> const &valarray) {
	return valarray.end();
}


template<typename Type, typename Integer, typename Derived>
struct Handle {
	Handle(Expression<Type, Integer, Derived> const &reference) :
			reference_(reference) {
	}
	Type operator[](Integer index) const {
		return reference_[index];
	}
private:
	Expression<Type, Integer, Derived> const &reference_;
};

template<typename Type, typename Integer, typename Derived>
struct Expression {
	VAL_ARRAY_UNARY_OPERATOR(+)
	VAL_ARRAY_UNARY_OPERATOR(-)
	VAL_ARRAY_UNARY_OPERATOR(!)
	VAL_ARRAY_UNARY_OPERATOR(~)
	VAL_ARRAY_BINARY_OPERATOR(+)
	VAL_ARRAY_BINARY_OPERATOR(-)
	VAL_ARRAY_BINARY_OPERATOR(*)
	VAL_ARRAY_BINARY_OPERATOR(/)
	VAL_ARRAY_BINARY_OPERATOR(%)
	VAL_ARRAY_BINARY_OPERATOR(&&)
	VAL_ARRAY_BINARY_OPERATOR(||)
	VAL_ARRAY_BINARY_OPERATOR(&)
	VAL_ARRAY_BINARY_OPERATOR(|)
	VAL_ARRAY_BINARY_OPERATOR(^)
	VAL_ARRAY_BINARY_OPERATOR(>>)
	VAL_ARRAY_BINARY_OPERATOR(<<)
	VAL_ARRAY_COMPARE_OPERATOR(<)
	VAL_ARRAY_COMPARE_OPERATOR(>)
	VAL_ARRAY_COMPARE_OPERATOR(<=)
	VAL_ARRAY_COMPARE_OPERATOR(>=)
	VAL_ARRAY_COMPARE_OPERATOR(==)
	VAL_ARRAY_COMPARE_OPERATOR(!=)
	VAL_ARRAY_UNARY_FUNCTION(Expression, abs)
	VAL_ARRAY_UNARY_FUNCTION(Expression, exp)
	VAL_ARRAY_UNARY_FUNCTION(Expression, log)
	VAL_ARRAY_UNARY_FUNCTION(Expression, log10)
	VAL_ARRAY_UNARY_FUNCTION(Expression, sqrt)
	VAL_ARRAY_UNARY_FUNCTION(Expression, asin)
	VAL_ARRAY_UNARY_FUNCTION(Expression, acos)
	VAL_ARRAY_UNARY_FUNCTION(Expression, atan)
	VAL_ARRAY_UNARY_FUNCTION(Expression, sin)
	VAL_ARRAY_UNARY_FUNCTION(Expression, cos)
	VAL_ARRAY_UNARY_FUNCTION(Expression, tan)
	VAL_ARRAY_UNARY_FUNCTION(Expression, sinh)
	VAL_ARRAY_UNARY_FUNCTION(Expression, cosh)
	VAL_ARRAY_UNARY_FUNCTION(Expression, tanh)
	VAL_ARRAY_BINARY_FUNCTION(Expression, pow)
	VAL_ARRAY_BINARY_FUNCTION(Expression, atan2)
	VAL_ARRAY_BINARY_FUNCTION(Expression, min)
	VAL_ARRAY_BINARY_FUNCTION(Expression, max)
	VAL_ARRAY_BINARY_FUNCTION(Expression, copysign)
	Type operator[](Integer index) const {
		return static_cast<Derived const*>(this)->operator[](index);
	}
	auto getHandle() const {
		return Handle<Type, Integer, Derived> { *this };
	}
	Integer size() const {
		return static_cast<Derived const*>(this)->size();
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

#endif /* INCLUDE_VAL_ARRAY_HPP_ */
