#pragma once

#include <algorithm>
#include <memory>
#include <numeric>
#include <type_traits>
#include <vector>

#pragma once
#include <type_traits>

template<typename T, typename = int>
struct ValueType {
	using type = T;
};

template<typename T>
struct ValueType<T, std::void_t<typename T::value_type>> {
	using type = typename T::value_type;
};

#define VALARRAY_BINARY(name, op)                                               \
	template<typename Ea, typename Eb>                                          \
	struct name##Expression: public Expression<name##Expression<Ea, Eb>> {      \
		using value_type = typename ValueType<Ea>::type;                        \
		using const_pointer = value_type const*;                                \
		using const_reference = value_type const&;                              \
		using pointer = value_type*;                                            \
		using reference = value_type&;                                          \
		name##Expression(Ea const& e1, Eb const& e2) : expr1_(e1), expr2_(e2) { \
		}                                                                       \
		value_type access(size_t index) const {                                 \
			return expr1_[index] op expr2_[index];                              \
		}                                                                       \
	private:                                                                    \
		Ea expr1_;                                                              \
		Eb expr2_;                                                              \
	};                                                                          \
	auto operator op (auto const& arg1, auto const& arg2) {                     \
		auto const expr1 = toExpression(arg1); 									\
		auto const expr2 = toExpression(arg2);                                  \
		return name##Expression(expr1, expr2);                                  \
	}                                                                           \

#define VALARRAY_COMPARE(name, op)                                              \
	template<typename Ea, typename Eb>                                          \
	struct name##Expression: public Expression<name##Expression<Ea, Eb>> {      \
		using value_type = bool;                                                \
		using const_pointer = bool const*;                                      \
		using const_reference = bool const&;                                    \
		using pointer = bool;                                                   \
		using reference = bool&;                                                \
		name##Expression(Ea const& e1, Eb const& e2) : expr1_(e1), expr2_(e2) { \
		}                                                                       \
		bool access(size_t index) const {                                       \
			return expr1_[index] op expr2_[index];                              \
		}                                                                       \
	private:                                                                    \
		Ea expr1_;                                                              \
		Eb expr2_;                                                              \
	};                                                                          \
	auto operator op (auto const& arg1, auto const& arg2) {                     \
		auto const expr1 = toExpression(arg1); 									\
		auto const expr2 = toExpression(arg2);                                  \
		return name##Expression(expr1, expr2);                                  \
	}                                                                           \

#define VALARRAY_COMPOUND_ASSIGNMENT(_class, op)                           \
	template<typename T>                                                        \
	_class<T>& operator op##= (_class<T> const& dst, auto const& _src) {        \
		auto const src = toExpression(_src);                                    \
		size_t const count = dst.size();                                        \
		for(size_t i = 0; i != count; i++) {                                    \
			dst[i] op##= src[i];                                                \
		}                                                                       \
		return dst;                                                             \
	}                                                                           \

#define VALARRAY_UNARY(name, op)                                                \
	template<typename E>                                                        \
	struct name##Expression: public Expression<name##Expression<E>> {           \
		using value_type = typename ValueType<E>::type;                         \
		using const_pointer = value_type const*;                                \
		using const_reference = value_type const&;                              \
		using pointer = value_type*;                                            \
		using reference = value_type&;                                          \
		name##Expression(E const& e) : expr_(e) {                               \
		}                                                                       \
		value_type access(size_t index) const {                                 \
			return op expr_[index];                                             \
		}                                                                       \
	private:                                                                    \
		E expr_;                                                                \
	};                                                                          \
	auto operator op (auto const& arg) {                                        \
		auto const expr = toExpression(arg); 									\
		return name##Expression(expr);                                          \
	}                                                                           \


template<typename >
struct Expression;

template<typename >
struct IsExpression;

template<typename >
struct IsValarray;

template<typename >
struct Valarray;

template<typename T>
struct ValarrayExpression {
	using value_type = T;
	using const_pointer = T const*;
	using const_reference = T const&;
	using pointer = T*;
	using reference = T&;
	ValarrayExpression(Valarray<T> const&);
	virtual ~ValarrayExpression();
	const_reference access(size_t) const;
	size_t size() const;
private:
	std::shared_ptr<std::vector<T>> ptr_;
};

template<typename T>
struct ValueExpression {
	using value_type = T;
	using const_pointer = T const*;
	using const_reference = T const&;
	using pointer = T*;
	using reference = T&;
	ValueExpression(T const&);
	virtual ~ValueExpression();
	const_reference access(size_t) const;
	size_t size() const;
private:
	T value_;
};

template<typename T>
struct Valarray {
	using value_type = T;
	using const_pointer = T const*;
	using const_reference = T const&;
	using pointer = T*;
	using reference = T&;
	Valarray();
	Valarray(size_t);
	Valarray(const_reference, size_t);
	Valarray(const_pointer, size_t);
	Valarray(const Valarray&);
	Valarray(Valarray&&);
	Valarray(std::initializer_list<T>);
	virtual ~Valarray();
	Valarray& operator=(const_reference);
	Valarray& operator=(const Valarray&);
	Valarray& operator=(Valarray&&);
	Valarray& operator=(std::initializer_list<T>);
	reference operator[](size_t);
	const_reference operator[](size_t) const;
	void swap(Valarray&&);
	value_type sum() const;
	value_type min() const;
	value_type max() const;
	Valarray shift(int shft) const;
	Valarray cshift(int shft) const;
	Valarray apply(value_type(T)) const;
	Valarray apply(value_type(const_reference)) const;
	size_t size() const;
	friend class ValarrayExpression<T>;
private:
	std::shared_ptr<std::vector<T>> ptr_;
};

template<typename Derived>
struct Expression {
	Expression();
	virtual ~Expression();
	auto operator[](size_t) const;
	size_t size() const;
private:
	Derived const &derived_;
};

template<typename T>
auto toExpression(T const&);

template<typename >
struct IsExpression {
	static constexpr bool value = false;
};

template<typename T>
struct IsExpression<Expression<T>> {
	static constexpr bool value = true;
};

template<typename T>
struct IsValarray {
	static constexpr bool value = false;
};

template<typename T>
struct IsValarray<Valarray<T>> {
	static constexpr bool value = true;
};

template<typename T>
Valarray<T>::Valarray() {
}

template<typename T>
Valarray<T>::Valarray(size_t count) :
		ptr_(new std::vector<T>(count)) {
}

template<typename T>
Valarray<T>::Valarray(const_reference val, size_t count) :
		ptr_(new std::vector<T>(count, val)) {
}

template<typename T>
Valarray<T>::Valarray(const_pointer ptr, size_t count) :
		ptr_(new std::vector<T>(ptr, ptr + count)) {
}

template<typename T>
Valarray<T>::Valarray(const Valarray &other) :
		ptr_(new std::vector<T>(other.ptr_.begin(), other.ptr_.end())) {
}

template<typename T>
Valarray<T>::Valarray(Valarray &&other) :
		ptr_(new std::vector<T>(std::move(other))) {
}

template<typename T>
Valarray<T>::Valarray(std::initializer_list<T> list) :
		ptr_(new std::vector<T>(list.begin(), list.end())) {
}

template<typename T>
Valarray<T>::~Valarray() {
}

template<typename T>
Valarray<T>& Valarray<T>::operator=(const_reference val) {
	std::fill(ptr_->begin(), ptr_->end(), val);
}

template<typename T>
Valarray<T>& Valarray<T>::operator=(const Valarray &other) {
	*ptr_ = *other.ptr_;
}

template<typename T>
Valarray<T>& Valarray<T>::operator=(Valarray &&other) {
	ptr_ = std::move(other.ptr_);
}

template<typename T>
Valarray<T>& Valarray<T>::operator=(std::initializer_list<T> list) {
	*ptr_ = list;
}

template<typename T>
Valarray<T>::reference Valarray<T>::operator[](size_t index) {
	return (*ptr_)[index];
}

template<typename T>
Valarray<T>::const_reference Valarray<T>::operator[](size_t index) const {
	return (*ptr_)[index];
}

template<typename T>
void Valarray<T>::swap(Valarray<T> &&other) {
	other.swap(ptr_);
}

template<typename T>
Valarray<T>::value_type Valarray<T>::sum() const {
	return std::accumulate(ptr_->begin(), ptr_->end(), T(0));
}

template<typename T>
Valarray<T>::value_type Valarray<T>::min() const {
	return std::min_element(ptr_->begin(), ptr_->end());
}

template<typename T>
Valarray<T>::value_type Valarray<T>::max() const {
	return std::max_element(ptr_->begin(), ptr_->end());
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
Valarray<T> Valarray<T>::apply(value_type f(T)) const {
	Valarray arry(*this);
	std::transform(arry.begin(), arry.end(), f);
	return *this;
}

template<typename T>
Valarray<T> Valarray<T>::apply(value_type f(const_reference)) const {
	Valarray arry(*this);
	std::transform(arry.begin(), arry.end(), f);
	return *this;
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
ValarrayExpression<T>::const_reference ValarrayExpression<T>::access(size_t index) const {
	return (*ptr_)[index];
}

template<typename T>
size_t ValarrayExpression<T>::size() const {
	return ptr_->size();
}

template<typename T>
ValueExpression<T>::ValueExpression(T const &value) :
		value_(value) {
}

template<typename T>
ValueExpression<T>::~ValueExpression() {
}

template<typename T>
ValueExpression<T>::const_reference ValueExpression<T>::access(size_t index) const {
	return value_;
}

template<typename T>
size_t ValueExpression<T>::size() const {
	return 1;
}

template<typename T>
auto toExpression(T const &object) {
	if constexpr (IsExpression<T>::value) {
		return object;
	} else {
		if constexpr (IsValarray<T>::value) {
			return ValarrayExpression<typename T::value_type>(object);
		} else {
			return ValueExpression<T>(object);
		}
	}
}

template<typename Derived>
Expression<Derived>::Expression() :
		derived_(*static_cast<Derived*>(this)) {
}

template<typename Derived>
Expression<Derived>::~Expression() {
}

template<typename Derived>
auto Expression<Derived>::operator[](size_t index) const {
	return derived_.access(index);
}

template<typename Derived>
size_t Expression<Derived>::size() const {
	return derived_.size();
}

VALARRAY_BINARY(Add, +)
VALARRAY_BINARY(Subtract, -)
VALARRAY_BINARY(Multiply, *)
VALARRAY_BINARY(Divide, /)
VALARRAY_BINARY(Modulus, %)
VALARRAY_BINARY(LogicalAnd, &&)
VALARRAY_BINARY(LogicalOr, ||)
VALARRAY_BINARY(BitwiseAnd, &)
VALARRAY_BINARY(BitwiseOr, |)
VALARRAY_BINARY(BitwiseXor, ^)
VALARRAY_BINARY(ShiftRight, >>)
VALARRAY_BINARY(ShiftLeft, <<)
VALARRAY_COMPARE(Equal, ==)
VALARRAY_COMPARE(NotEqual, !=)
VALARRAY_COMPARE(Greater, >)
VALARRAY_COMPARE(Less, <)
VALARRAY_COMPARE(GreaterOrEqual, >=)
VALARRAY_COMPARE(LessOrEqual, <=)
VALARRAY_COMPOUND_ASSIGNMENT(Valarray, +)
VALARRAY_COMPOUND_ASSIGNMENT(Valarray, -)
VALARRAY_COMPOUND_ASSIGNMENT(Valarray, *)
VALARRAY_COMPOUND_ASSIGNMENT(Valarray, /)
VALARRAY_COMPOUND_ASSIGNMENT(Valarray, %)
VALARRAY_COMPOUND_ASSIGNMENT(Valarray, &)
VALARRAY_COMPOUND_ASSIGNMENT(Valarray, |)
VALARRAY_COMPOUND_ASSIGNMENT(Valarray, ^)
VALARRAY_COMPOUND_ASSIGNMENT(Valarray, >>)
VALARRAY_COMPOUND_ASSIGNMENT(Valarray, <<)
VALARRAY_UNARY(Positive, +)
VALARRAY_UNARY(Negative, -)
VALARRAY_UNARY(LogicalNot, !)
VALARRAY_UNARY(BitwiseNot, ~)


//#include <hpx/hpx_init.hpp>
//#include "Valarray.hpp"
//#include "dgTransforms.hpp"
//#include "Options.hpp"
//#include "MultiIndex.hpp"
//#include "HyperGrid.hpp"
//#include "EulerState.hpp"
//#include "Real.hpp"
//#include "RungeKutta.hpp"

void testRadiation();

int hpx_main(int argc, char *argv[]) {
	Valarray<double> v1;
	Valarray<double> v2;
	auto const v3 = v1 < v2;
//	printf("\nStarting\n");
//	enableFPE();
//	processOptions(argc, argv);
//	printf("\nPrologue complete\n");
//	constexpr int P = 4;
//	constexpr int D = 2;
//	constexpr int N = 128;
//	using T = Real;
//	using RK = Tname RungeKutta<T, P>::T;
//	HyperGrid<T, D, N, P, RK, EulerStateHLLC> grid;
//	grid.initialize(initSodShockTube<T, D>);
//	grid.output("X", 0, Real(0.0));
//	grid.applyLimiter();
//	T t = T(0);
//	T tmax = T(.125);
//	T dt;
//	RK const rk;
//	int iter = 0;
//	while (t < tmax) {
//		grid.output("X", iter, t);
//		std::cout << "i = " << std::to_string(iter);
//		std::cout << "  t = " << std::to_string(t);
//		dt = grid.beginStep();
//		std::cout << "  dt = " << dt << std::endl;
//		for (int s = 0; s < rk.stageCount(); s++) {
//			grid.subStep(dt, s);
//		}
//		grid.endStep();
//		iter++;
//		t += dt;
//	}
//	printf("\nStopping\n");
//	return hpx::local::finalize();
}

int main(int argc, char *argv[]) {
//#ifndef NDEBUG
//	installFpeHandler();
//#endif
//	std::vector<std::string> cfg = { "hpx.commandline.allow_unknown=1" };
//	cfg.push_back("hpx.stacks.small_size=1048576");
//	hpx::init_params init_params;
//	init_params.cfg = std::move(cfg);
//	auto rc = hpx::init(argc, argv, init_params);
	return 0;
//	return rc;
}
