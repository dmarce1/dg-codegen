#ifndef PHYSICAL_CONSTANTS_CGS_HPP
#define PHYSICAL_CONSTANTS_CGS_HPP

#include <ostream>
#include <iomanip>

template<typename Type>
struct PhysicalConstants {
	Type c;
	Type G;
	Type sigma;
	Type kB;
	Type amu;
	Type a;
};

// pretty-print PhysicalConstants in scientific notation
template<typename Type>
inline std::ostream& operator<<(std::ostream &os, PhysicalConstants<Type> const &pc) {
	os << std::scientific << std::setprecision(8) << "PhysicalConstants {\n" << "  c     = " << pc.c << "\n" << "  G     = " << pc.G << "\n" << "  sigma = "
			<< pc.sigma << "\n" << "  kB    = " << pc.kB << "\n" << "  amu   = " << pc.amu << "\n" << "  a     = " << pc.a << "\n" << "}";
	return os;
}

template<typename Type>
class PhysicalUnits {
	static constexpr Type cgs_sigma = 5.670374419e-5;
	static constexpr Type cgs_kB = 1.380649e-16;
	static constexpr Type cgs_c = 2.99792458e10;
	static constexpr Type cgs_amu = 1.66053906660e-24;
	static constexpr Type cgs_G = 6.67430e-8;

	Type code2cm_;
	Type code2s_;
	Type code2g_;

public:
	PhysicalUnits(Type Lo = 1.0) {
		code2cm_ = Lo;
		code2s_ = code2cm_ / cgs_c;
		code2g_ = code2cm_ * cgs_c * cgs_c / cgs_G;
	}

	PhysicalUnits(Type Lo, Type To, Type Mo) {
		code2cm_ = Lo;
		code2s_ = To;
		code2g_ = Mo;
	}

	Type getCode2cm() const {
		return code2cm_;
	}
	Type getCode2s() const {
		return code2s_;
	}
	Type getCode2g() const {
		return code2g_;
	}

	PhysicalConstants<Type> getPhysicalConstants() const {
		PhysicalConstants<Type> constants;
		constants.c = cgs_c * code2s_ / code2cm_;
		constants.G = cgs_G * code2s_ * code2s_ * code2g_ / (code2cm_ * code2cm_ * code2cm_);
		constants.sigma = cgs_sigma * code2s_ * code2s_ * code2s_ / code2g_;
		constants.kB = cgs_kB * code2s_ * code2s_ / (code2cm_ * code2cm_ * code2g_);
		constants.amu = cgs_amu / code2g_;
		constants.a = 4 * constants.sigma / constants.c;
		return constants;
	}
};

// pretty-print PhysicalUnits in scientific notation
template<typename Type>
inline std::ostream& operator<<(std::ostream &os, PhysicalUnits<Type> const &pu) {
	auto pc = pu.getPhysicalConstants();
	os << std::scientific << std::setprecision(8) << "PhysicalUnits {\n" << "  code2cm = " << pu.getCode2cm() << " [cm]\n" << "  code2s  = " << pu.getCode2s()
			<< " [s]\n" << "  code2g  = " << pu.getCode2g() << " [g]\n\n" << "  " << pc << "\n" << "}";
	return os;
}
#endif /* INCLUDE_CONSTANTS_HPP_ */
