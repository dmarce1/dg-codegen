#pragma once

template<int M>
struct Octant {
	constexpr Octant(int j) :
			i(j) {
	}
	constexpr Octant() {
	}
	inline constexpr Octant& operator++() {
		i++;
		return *this;
	}
	inline constexpr Octant& operator--() {
		i--;
		return *this;
	}
	inline constexpr Octant operator++(int) {
		auto const rc = *this;
		operator++();
		return rc;
	}
	inline constexpr Octant operator--(int) {
		auto const rc = *this;
		operator--();
		return rc;
	}
	inline constexpr bool operator[](int j) const {
		return (1 & (i >> j));
	}
	inline static constexpr Octant begin() {
		Octant o;
		o.i = 0;
		return o;
	}
	inline static constexpr Octant end() {
		Octant o;
		o.i = 1 << M;
		return o;
	}
	inline constexpr bool operator==(Octant const &other) const {
		return i == other.i;
	}
	inline constexpr bool operator!=(Octant const &other) const {
		return i != other.i;
	}
	inline constexpr bool operator<(Octant const &other) const {
		return i < other.i;
	}
	inline constexpr bool operator>(Octant const &other) const {
		return i > other.i;
	}
	inline constexpr bool operator<=(Octant const &other) const {
		return i <= other.i;
	}
	inline constexpr bool operator>=(Octant const &other) const {
		return i >= other.i;
	}
	inline constexpr static int size() {
		return 1 << M;
	}
	inline constexpr operator int() const {
		return i;
	}
private:
	int i;
};



