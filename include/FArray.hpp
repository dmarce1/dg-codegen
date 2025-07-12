//#ifndef INCLUDE_FARRAY_HPP_
//#define INCLUDE_FARRAY_HPP_
//
//#include <array>
//#include <initializer_list>
//#include <limits>
//#include <utility>
//
//using FArrayInt = int64_t;
//
//struct Bounds {
//	FArrayInt lowerBound;
//	FArrayInt upperBound;
//};
//
//template<FArrayInt D>
//struct FRange {
//	FRange() {
//		for (int d = 0; d < D; d++) {
//			bounds_[d].lowerBound = bounds_[d].upperBound = 0;
//		}
//	}
//	FRange(FArrayInt n) {
//		for (int d = 0; d < D; d++) {
//			bounds_[d].lowerBound = 0;
//			bounds_[d].upperBound = n;
//		}
//	}
//	FRange(FArrayInt b, FArrayInt e) {
//		for (int d = 0; d < D; d++) {
//			bounds_[d].lowerBound = b;
//			bounds_[d].upperBound = e;
//		}
//	}
//	FRange(std::array<FArrayInt, D> const &b, std::array<FArrayInt, D> const &e) {
//		for (int d = 0; d < D; d++) {
//			bounds_[d].lowerBound = b[d];
//			bounds_[d].upperBound = e[d];
//		}
//	}
//	FArrayInt span(int d) const {
//		return bounds_[d].upperBound - bounds_[d].lowerBound;
//	}
//	FArrayInt size() const {
//		FArrayInt vol = 1;
//		for (int d = 0; d < D; d++) {
//			vol *= span(d);
//		}
//		return vol;
//	}
//	FArrayInt const& lower(int d) const {
//		return bounds_[d].lowerBound;
//	}
//	FArrayInt const& upper(int d) const {
//		return bounds_[d].upperBound;
//	}
//	FArrayInt& lower(int d) {
//		return bounds_[d].lowerBound;
//	}
//	FArrayInt& upper(int d) {
//		return bounds_[d].upperBound;
//	}
//private:
//	std::array<Bounds, D> bounds_;
//};
//
//template<int D>
//class FArrayView {
//	constexpr FArrayInt nullEntry = std::numeric_limits<FArrayInt>::min();
//public:
//	template<std::same_as<std::initializer_list<FArrayInt>> ... Lists>
//	FArrayView(Lists ... lists) {
//		static_assert(sizeof...(Lists) == D);
//		std::array<std::initializer_list<FArrayInt>, D> listArray { lists... };
//		for (std::size_t d = 0; d < D; ++d) {
//			auto const &list = listArray[d];
//			switch (l.size()) {
//			case 0:
//				innerBox_.lower(d) = nullEntry;
//				innerBox_.upper(d) = nullEntry;
//				skipStrides_[d] = nullEntry;
//				break;
//			case 1:
//				innerBox_.lower(d) = 0;
//				innerBox_.upper(d) = list[0];
//				skipStrides_[d] = nullEntry;
//				break;
//				break;
//			case 2:
//				break;
//			default:
//				break;
//			}
//		}
//	}
//private:
//	FRange<D> innerBox_;
//	std::array<FArrayInt, D> skipStrides_;
//};
//
//template<typename T, FArrayInt D>
//class FArray {
//	using const_pointer = typename std::vector<T>::const_pointer;
//	using const_reference = typename std::vector<T>::const_reference;
//	using pointer = typename std::vector<T>::pointer;
//	using reference = typename std::vector<T>::reference;
//	using value_type = typename std::vector<T>::value_type;
//	void init() {
//		for (int d = 0; d < D; d++) {
//			sizes_[d] = box_.span(d);
//		}
//		strides_[D - 1] = 1;
//		for (int d = D - 2; d >= 0; d--) {
//			strides_[d] = strides_[d + 1] * sizes_[d + 1];
//		}
//	}
//	FArrayInt flatIndex(std::array<FArrayInt, D> const &indices) const {
//		size_t index = indices[D - 1] - box_[D - 1].lowerBound;
//		for (int d = D - 2; d >= 0; d--) {
//			index += strides_[d] * (indices[d] - box_.lower(d));
//		}
//		return index;
//	}
//public:
//	FArray() {
//		ptr_ = std::make_shared<std::vector<T>>(0);
//		init();
//	}
//	FArray(FRange<D> const &box) {
//		box_ = box;
//		ptr_ = std::make_shared<std::vector<T>>(box_.size());
//		init();
//	}
//	FArray(FArray const &other) {
//		*this = other;
//	}
//	FArray(FArray &&other) {
//		*this = std::move(other);
//	}
//	FArray& operator=(FArray const &other) {
//		box_ = other.box;
//		ptr_ = std::make_shared<std::vector<T>>(box_.size());
//		init();
//		*ptr_ = *other.ptr_;
//		return *this;
//	}
//	FArray& operator=(FArray &&other) {
//		box_ = std::move(other.box);
//		ptr_ = std::make_shared<std::vector<T>>(box_.size());
//		init();
//		*ptr_ = std::move(*other.ptr_);
//		return *this;
//	}
//	const_reference operator[](std::array<FArrayInt, D> const &indices) const {
//		return (*ptr_)[flatIndex(indices)];
//	}
//	reference operator[](std::array<FArrayInt, D> const &indices) {
//		return (*ptr_)[flatIndex(indices)];
//	}
//private:
//	FRange<D> box_;
//	std::array<FArrayInt, D> sizes_;
//	std::array<FArrayInt, D> strides_;
//	std::shared_ptr<std::vector<T>> ptr_;
//};
//
//#endif
