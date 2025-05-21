#pragma once

#include <cassert>
#include <concepts>
#include <cstddef>
#include <experimental/simd>
#include <vector>

template<typename Derived>
struct Expression {
	inline double operator[](int i) const {
		return static_cast<Derived const&>(*this)[i];
	}
	inline int size() const {
		return static_cast<Derived const&>(*this).size();
	}
};

//template<typename LHS, typename RHS>
//struct AddExpr: Expression<AddExpr<LHS, RHS>> {
//	LHS const &lhs;
//	RHS const &rhs;
//
//	AddExpr(LHS const &lhs_, RHS const &rhs_) :
//			lhs(lhs_), rhs(rhs_) {
//	}
//
//	inline double operator[](size_t i) const {
//		return lhs[i] + rhs[i];
//	}
//
//	inline size_t size() const {
//		return lhs.size(); // assume same size
//	}
//};

template<typename Type>
struct SimdArray: Expression<SimdArray> {
	using simd_type = std::experimental::native_simd<Type>;
	static constexpr size_t simdWidth = simd_type::size();
	static constexpr size_t simdWidthBytes = sizeof(Type) * simdWidth;
	static constexpr bool isPod = std::is_trivial<Type>::value && std::is_standard_layout<Type>::value;
	SimdArray(size_t n) {
		logicalSize = 0;
		allocate(n);
	}
	~SimdArray() {
		deallocate();
	}
	double& operator[](size_t i) {
		return data[i];
	}
	double operator[](size_t i) const {
		return data[i];
	}
	size_t size() const {
		return data.size();
	}

	template<typename Expr>
	SimdArray& operator=(Expr const &expr) {
		for (size_t i = 0; i < size(); ++i) {
			data[i] = expr[i];
		}
		return *this;
	}
private:
	void allocate(size_t n, Type initVal = Type()) {
		size_t newLogicalSize = n;
		size_t const newPaddedSize = simdWidth * ((newLogicalSize + simdWidth - 1) / simdWidth);
		if (newPaddedSize > PaddedSize) {
			raw = malloc(newPaddedSize + simdWidth - 1);
			aligned = (T*) ((uintptr_t(raw) + simdWidthBytes - 1) & ~(simdWidthBytes - 1));
			paddedSize = newPaddedSize;
		}
		if constexpr (!isPod) {
			if (newLogicalSize != logicalSize) {
				if (newLogicalSize < logicalSize) {
					for (size_t i = newLogicalSize; i < logicalSize; i++) {
						(aligned + i)->~Type();
					}
				} else {
					for (size_t i = logicalSize; i < newLogicalSize; i++) {
						new (aligned + i) Type(initVal);
					}
				}
				logicalSize = newLogicalSize;
			}
		}
	}
	void deallocate() {
		if (raw) {
			if constexpr (!isPod) {
				(aligned + i)->~Type();
			}
			free(raw);
		}
	}
	T *aligned = nullptr;
	void *raw = nullptr;
	size_t logicalSize = 0;
	size_t paddedSize = 0;
};

// Operator overload
template<typename LHS, typename RHS>
inline AddExpr<LHS, RHS> operator+(Expression<LHS> const &a, Expression<RHS> const &b) {
	return AddExpr<LHS, RHS>(a, b);
}
#endif /* INCLUDE_ARRAY_HPP_ */
