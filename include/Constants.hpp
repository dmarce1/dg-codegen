#ifndef PHYSICAL_CONSTANTS_CGS_HPP
#define PHYSICAL_CONSTANTS_CGS_HPP

#include <ostream>
#include <iomanip>

struct PhysicalConstants {
	double c;
	double G;
	double sigma;
	double kB;
	double amu;
	double a;
};

// pretty-print PhysicalConstants in scientific notation
inline std::ostream& operator<<(std::ostream &os, PhysicalConstants const &pc) {
	os << std::scientific << std::setprecision(8) << "PhysicalConstants {\n" << "  c     = " << pc.c << "\n" << "  G     = " << pc.G << "\n" << "  sigma = "
			<< pc.sigma << "\n" << "  kB    = " << pc.kB << "\n" << "  amu   = " << pc.amu << "\n" << "  a     = " << pc.a << "\n" << "}";
	return os;
}

class PhysicalUnits {
	static constexpr double cgs_sigma = 5.670374419e-5;
	static constexpr double cgs_kB = 1.380649e-16;
	static constexpr double cgs_c = 2.99792458e10;
	static constexpr double cgs_amu = 1.66053906660e-24;
	static constexpr double cgs_G = 6.67430e-8;

	double code2cm_;
	double code2s_;
	double code2g_;

public:
	PhysicalUnits(double Lo = 1.0) {
		code2cm_ = Lo;
		code2s_ = code2cm_ / cgs_c;
		code2g_ = code2cm_ * cgs_c * cgs_c / cgs_G;
	}

	PhysicalUnits(double Lo, double To, double Mo) {
		code2cm_ = Lo;
		code2s_ = To;
		code2g_ = Mo;
	}

	double getCode2cm() const {
		return code2cm_;
	}
	double getCode2s() const {
		return code2s_;
	}
	double getCode2g() const {
		return code2g_;
	}

	PhysicalConstants getPhysicalConstants() const {
		PhysicalConstants constants;
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
inline std::ostream& operator<<(std::ostream &os, PhysicalUnits const &pu) {
	auto pc = pu.getPhysicalConstants();
	os << std::scientific << std::setprecision(8) << "PhysicalUnits {\n" << "  code2cm = " << pu.getCode2cm() << " [cm]\n" << "  code2s  = " << pu.getCode2s()
			<< " [s]\n" << "  code2g  = " << pu.getCode2g() << " [g]\n\n" << "  " << pc << "\n" << "}";
	return os;
}
#endif /* INCLUDE_CONSTANTS_HPP_ */
