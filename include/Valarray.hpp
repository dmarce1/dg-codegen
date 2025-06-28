#ifndef INCLUDE_VALARRAY_HPP_
#define INCLUDE_VALARRAY_HPP_

#include <cstdlib>
#include <initializer_list>
#include <functional>
#include <type_traits>
#include <utility>

#define CONCAT_IMPL(x, y) x##y
#define CONCAT(x, y) CONCAT_IMPL(x, y)
#define UNIQUE_NAME(base) CONCAT(base, __LINE__)

#define VALARRY_UNARY_OPERATOR(op)                                                                 \
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

#define VALARRY_BINARY_OPERATOR(op)                                                                \
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
		return UNIQUE_NAME(Expression)(this->getHandle(), other->getHandle());                       \
	}


#define VALARRY_COMPARE_OPERATOR(op)                                                               \
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

#define VALARRY_COMPOUND_ASSIGNMENT_OPERATOR(OP)                                     \
	template<typename Derived>                                                        \
	Valarray& operator OP##= (Expression<Type, Integer, Derived> const &expression) { \
		for (Integer i = 0; i < size_; i++) {                                          \
			data_[i] OP##= expression[i];                                               \
		}                                                                              \
		return *this;                                                                  \
	}                                                                                 \

#define VALARRY_UNARY_FUNCTION(class_, name)                                         \
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

#define VALARRY_BINARY_FUNCTION(class_, name)                                          \
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

template<typename Type, typename Integer, typename Derived>
struct Expression;

template<typename Type, typename Integer, typename Handle>
struct UnaryExpression;

template<typename Type, typename Integer, typename HandleA, typename HandleB>
struct ComparisonExpression;

template<typename Type, typename Integer, typename HandleA, typename HandleB>
struct BinaryExpression;

template<typename Type, typename Integer = int>
struct Valarray;

template<typename Type>
struct UnaryPlus;

template<typename Type>
struct ShiftLeft;

template<typename Type>
struct ShiftRight;

template<typename Type, typename Integer>
struct Valarray {
	Valarray() :
			data_(allocate(0)), size_(0) {
	}
	explicit Valarray(Integer count) :
			data_(allocate(count)), size_(count) {
		for (Integer i = 0; i < size_; i++) {
			new (data_ + i) Type { };
		}
	}
	Valarray(const Type &value, Integer count) :
			data_(allocate(count)), size_(count) {
		if constexpr (isPod) {
			*this = value;
		} else {
			for (Integer i = 0; i < size_; i++) {
				new (data_ + i) Type { value };
			}
		}
	}
	Valarray(const Type *valuesPointer, Integer count) :
			data_(allocate(count)), size_(count) {
		if constexpr (isPod) {
			std::copy_n(valuesPointer, size_, data_);
		} else {
			for (Integer i = 0; i < size_; i++) {
				new (data_ + i) Type { valuesPointer[i] };
			}
		}
	}
	Valarray(const Valarray &other) {
		*this = other;
	}
	Valarray(Valarray &&other) :
			data_(nullptr), size_(0) {
		*this = std::move(other);
	}
	Valarray(std::initializer_list<Type> valuesList) :
			data_(valuesList.size()), size_(valuesList.size()) {
		if constexpr (isPod) {
			std::copy_n(valuesList.data_, size_, data_);
		} else {
			auto const *vals = valuesList.data();
			for (Integer i = 0; i < size_; i++) {
				new (data_ + i) Type { vals[i] };
			}
		}
	}
	virtual ~Valarray() {
		deallocate(data_);
	}
	Valarray& operator=(const Type &value) {
		std::fill_n(data_, size_, value);
		return *this;
	}
	Valarray& operator=(const Valarray &other) {
		deallocate(data_);
		data_ = allocate(other.size_);
		size_ = other.size_;
		if constexpr (isPod) {
			std::copy_n(other.data_, size_, data_);
		} else {
			for (Integer i = 0; i < size_; i++) {
				new (data_ + i) Type { other.data_[i] };
			}
		}
		return *this;
	}
	Valarray& operator=(Valarray &&other) {
		swap(std::move(other));
	}
	Valarray& operator=(std::initializer_list<Type> valuesList) {
		deallocate(data_);
		*this = Valarray { std::move(valuesList) };
	}
	const Type& operator[](Integer index) const {
		return data_[index];
	}
	Type& operator[](Integer index) {
		return data_[index];
	}
	void swap(Valarray<Type, Integer> &&other) {
		std::swap(data_, other.data_);
		std::swap(size_, other.size_);
	}
	Integer size() const {
		return size_;
	}
	void resize(Integer count, Type value = Type { }) {
		*this = Valarray<Type, Integer>(value);
	}
	Type sum() const {
		Type result = data_[0];
		for (Integer i = 1; i < size_; i++) {
			result += data_[i];
		}
		return result;
	}
	Type min() const {
		Type result = data_[0];
		for (Integer i = 1; i < size_; i++) {
			Type const &thisElement = data_[i];
			if (thisElement < result) {
				result = thisElement;
			}
		}
		return result;
	}
	Type max() const {
		Type result = data_[0];
		for (Integer i = 1; i < size_; i++) {
			Type const &thisElement = data_[i];
			if (result < thisElement) {
				result = thisElement;
			}
		}
		return result;
	}
	Valarray shift(int thisShift) const {
		Valarray result(size_);
		if (thisShift > 0) {
			std::copy_n(data_ + thisShift, size_ - thisShift, result.data_);
		} else if (thisShift < 0) {
			std::copy_n(data_, size_ + thisShift, result.data_ - thisShift);
		}
		return result;
	}
	Valarray cshift(int thisShift) const {
		Valarray result(size_);
		if (thisShift > 0) {
			std::copy_n(data_ + thisShift, size_ - thisShift, result.data_);
			std::copy_n(data_, thisShift, result.data_ + size_ - thisShift);
		} else if (thisShift < 0) {
			std::copy_n(data_, size_ + thisShift, result.data_ - thisShift);
			std::copy_n(data_ + size_ + thisShift, -thisShift, result.data_);
		}
		return result;
	}
	Valarray apply(Type thisFunction(Type)) const {
		Valarray result(size_);
		for (Integer i = 0; i < size_; i++) {
			result.data_[i] = thisFunction(data_);
		}
		return result;
	}
	Valarray apply(Type thisFunction(const Type&)) const {
		Valarray result(size_);
		for (Integer i = 0; i < size_; i++) {
			result.data_[i] = thisFunction(data_);
		}
		return result;
	}
	template<typename Derived>
	Valarray& operator=(Expression<Type, Integer, Derived> const &expression) {
		for (Integer i = 0; i < size_; i++) {
			data_[i] = expression[i];
		}
		return *this;
	}
	VALARRY_UNARY_OPERATOR(+)
	VALARRY_UNARY_OPERATOR(-)
	VALARRY_UNARY_OPERATOR(!)
	VALARRY_UNARY_OPERATOR(~)
	VALARRY_BINARY_OPERATOR(+)
	VALARRY_BINARY_OPERATOR(-)
	VALARRY_BINARY_OPERATOR(*)
	VALARRY_BINARY_OPERATOR(/)
	VALARRY_BINARY_OPERATOR(%)
	VALARRY_BINARY_OPERATOR(&&)
	VALARRY_BINARY_OPERATOR(||)
	VALARRY_BINARY_OPERATOR(&)
	VALARRY_BINARY_OPERATOR(|)
	VALARRY_BINARY_OPERATOR(^)
	VALARRY_BINARY_OPERATOR(>>)
	VALARRY_BINARY_OPERATOR(<<)
	VALARRY_COMPARE_OPERATOR(<)
	VALARRY_COMPARE_OPERATOR(>)
	VALARRY_COMPARE_OPERATOR(<=)
	VALARRY_COMPARE_OPERATOR(>=)
	VALARRY_COMPARE_OPERATOR(==)
	VALARRY_COMPARE_OPERATOR(!=)
   VALARRY_UNARY_FUNCTION(Valarray, abs)
   VALARRY_UNARY_FUNCTION(Valarray, exp)
   VALARRY_UNARY_FUNCTION(Valarray, log)
	VALARRY_UNARY_FUNCTION(Valarray, log10)
   VALARRY_UNARY_FUNCTION(Valarray, sqrt)
	VALARRY_UNARY_FUNCTION(Valarray, asin)
   VALARRY_UNARY_FUNCTION(Valarray, acos)
   VALARRY_UNARY_FUNCTION(Valarray, atan)
	VALARRY_UNARY_FUNCTION(Valarray, sin)
   VALARRY_UNARY_FUNCTION(Valarray, cos)
   VALARRY_UNARY_FUNCTION(Valarray, tan)
	VALARRY_UNARY_FUNCTION(Valarray, sinh)
   VALARRY_UNARY_FUNCTION(Valarray, cosh)
   VALARRY_UNARY_FUNCTION(Valarray, tanh)
   VALARRY_BINARY_FUNCTION(Valarray, pow)
   VALARRY_BINARY_FUNCTION(Valarray, atan2)
	VALARRY_COMPOUND_ASSIGNMENT_OPERATOR(+)
	VALARRY_COMPOUND_ASSIGNMENT_OPERATOR(-)
	VALARRY_COMPOUND_ASSIGNMENT_OPERATOR(*)
	VALARRY_COMPOUND_ASSIGNMENT_OPERATOR(/)
	VALARRY_COMPOUND_ASSIGNMENT_OPERATOR(%)
	VALARRY_COMPOUND_ASSIGNMENT_OPERATOR(&)
	VALARRY_COMPOUND_ASSIGNMENT_OPERATOR(|)
	VALARRY_COMPOUND_ASSIGNMENT_OPERATOR(^)
	VALARRY_COMPOUND_ASSIGNMENT_OPERATOR(<<)
	VALARRY_COMPOUND_ASSIGNMENT_OPERATOR(>>)
private:
	template<typename Derived>
	operator Expression<Type, Integer, Derived>() const {
		return operator+();
	}
	auto const* getHandle() const {
		return data_;
	}
	static Type* allocate(Integer count) {
		return (Type*) malloc(sizeof(Type) * count);
	}
	static void deallocate(Type *pointer) {
		if (pointer) {
			free(pointer);
		}
	}
	static constexpr bool isPod = std::is_standard_layout<Type>::value && std::is_trivial<Type>::value;
	Type *data_;
	Integer size_;
};

template<typename Type>
void swap(Valarray<Type> &&valarrayA, Valarray<Type> &&valarrayB) {
	valarrayA.swap(std::move(valarrayB));
}

template<typename Type>
Type* begin(Valarray<Type> &valarray) {
	return &(valarray[0]);
}

template<typename Type>
Type* end(Valarray<Type> &valarray) {
	return &(valarray[0]) + valarray.size();
}

template<typename Type>
Type const* begin(Valarray<Type> const &valarray) {
	return &(valarray[0]);
}

template<typename Type>
Type const* end(Valarray<Type> const &valarray) {
	return &(valarray[0]) + valarray.size();
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
	VALARRY_UNARY_OPERATOR(+)
	VALARRY_UNARY_OPERATOR(-)
	VALARRY_UNARY_OPERATOR(!)
	VALARRY_UNARY_OPERATOR(~)
	VALARRY_BINARY_OPERATOR(+)
	VALARRY_BINARY_OPERATOR(-)
	VALARRY_BINARY_OPERATOR(*)
	VALARRY_BINARY_OPERATOR(/)
	VALARRY_BINARY_OPERATOR(%)
	VALARRY_BINARY_OPERATOR(&&)
	VALARRY_BINARY_OPERATOR(||)
	VALARRY_BINARY_OPERATOR(&)
	VALARRY_BINARY_OPERATOR(|)
	VALARRY_BINARY_OPERATOR(^)
	VALARRY_BINARY_OPERATOR(>>)
	VALARRY_BINARY_OPERATOR(<<)
	VALARRY_COMPARE_OPERATOR(<)
	VALARRY_COMPARE_OPERATOR(>)
	VALARRY_COMPARE_OPERATOR(<=)
	VALARRY_COMPARE_OPERATOR(>=)
	VALARRY_COMPARE_OPERATOR(==)
	VALARRY_COMPARE_OPERATOR(!=)
   VALARRY_UNARY_FUNCTION(Expression, abs)
   VALARRY_UNARY_FUNCTION(Expression, exp)
   VALARRY_UNARY_FUNCTION(Expression, log)
	VALARRY_UNARY_FUNCTION(Expression, log10)
   VALARRY_UNARY_FUNCTION(Expression, sqrt)
	VALARRY_UNARY_FUNCTION(Expression, asin)
   VALARRY_UNARY_FUNCTION(Expression, acos)
   VALARRY_UNARY_FUNCTION(Expression, atan)
	VALARRY_UNARY_FUNCTION(Expression, sin)
   VALARRY_UNARY_FUNCTION(Expression, cos)
   VALARRY_UNARY_FUNCTION(Expression, tan)
	VALARRY_UNARY_FUNCTION(Expression, sinh)
   VALARRY_UNARY_FUNCTION(Expression, cosh)
   VALARRY_UNARY_FUNCTION(Expression, tanh)
   VALARRY_BINARY_FUNCTION(Expression, pow)
   VALARRY_BINARY_FUNCTION(Expression, atan2)
 	Type operator[](Integer index) const {
		return Derived::operator[](index);
	}
	auto getHandle() const {
		return Handle { this };
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

#endif /* INCLUDE_VALARRAY_HPP_ */
