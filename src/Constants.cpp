/******************************************************************************
 Copyright (C) 2024  Dominic C. Marcello
 *******************************************************************************/

#include "Constants.hpp"

Units::Units(double code2cm) {
	Constants const cgsConstants = getCgsConstants();
	code2cm_ = code2cm;
	code2g_ = code2cm_ * sqr(cgsConstants.c) / cgsConstants.G;
	code2s_ = code2cm_ / cgsConstants.c;
	code2K_ = 1.0;
}

double Units::code2cm() const {
	return code2cm_;
}

double Units::code2s() const {
	return code2s_;
}

double Units::code2g() const {
	return code2g_;
}

double Units::code2K() const {
	return code2K_;
}

double Units::code2dyne() const {
	return code2g_ * code2cm_ / sqr(code2s_);
}

double Units::code2erg() const {
	return code2dyne() * code2cm_;
}

double Units::code2Hz() const {
	return 1.0 / code2s_;
}

double Units::cm2code() const {
	return 1.0 / code2cm_;
}

double Units::s2code() const {
	return 1.0 / code2s_;
}

double Units::g2code() const {
	return 1.0 / code2g_;
}

double Units::K2code() const {
	return 1.0 / code2K();
}

double Units::dyne2code() const {
	return 1.0 / code2dyne();
}

double Units::erg2code() const {
	return 1.0 / code2erg();
}

double Units::Hz2code() const {
	return 1.0 / code2Hz();
}

std::ostream& operator<<(std::ostream &os, Units const &u) {
	os << std::setprecision(10);
	os << "Units (code → CGS):\n";
	os << "  code2cm = " << u.code2cm_ << " cm\n";
	os << "  code2s  = " << u.code2s_ << " s\n";
	os << "  code2g  = " << u.code2g_ << " g\n";
	os << "  code2Hz   = " << 1.0 / u.code2s_ << " Hz\n";
	os << "  code2dyne = " << u.code2g_ * u.code2cm_ / (u.code2s_ * u.code2s_) << " dyn\n";
	os << "  code2erg  = " << u.code2g_ * std::pow(u.code2cm_, 2) / std::pow(u.code2s_, 2) << " erg\n";
	return os;
}

static Units units { };

void setLengthScale(double code2cm) {
	units = Units(code2cm);
}

Units const& getUnits() {
	return units;
}

Constants getCgsConstants() {
	Constants constants;
	constants.c = 2.99792458e10;
	constants.G = 6.67259e-8;
	constants.amu = 1.6605402e-24;
	constants.kB = 1.380658e-16;
	constants.h = 6.6260755e-27;
	constants.hbar = constants.h / (2.0 * constants.pi);
	constants.sigma = (sqr(constants.pi) * ipow(constants.kB, 4)) / (60.0 * ipow(constants.hbar, 3) * sqr(constants.c));
	constants.aR = 4.0 * constants.sigma / constants.c;
	constants.Msol = 1.99e33;
	constants.Rsol = 6.96e10;
	constants.Lsol = 3.9e33;
	constants.AU = 1.496e13;
	constants.pc = constants.AU * (360.0 * 60.0 * 60.0) / (2.0 * constants.pi);
	constants.ly = constants.c * (365.25 * 24.0 * 60.0 * 60.0);
	return constants;
}

Constants getCodeConstants() {
	auto const cgs = getCgsConstants();
	Constants constants;
	constants.c = cgs.c * units.cm2code() / units.s2code();
	constants.G = cgs.G * ipow(units.cm2code(), 3) / (units.g2code() * sqr(units.s2code()));
	constants.amu = cgs.amu * units.g2code();
	constants.kB = cgs.kB * units.erg2code() / units.K2code();
	constants.h = cgs.h * units.erg2code() / units.Hz2code();
	constants.hbar = constants.h / (2.0 * constants.pi);
	constants.sigma = (sqr(constants.pi) * ipow(constants.kB, 4)) / (60.0 * ipow(constants.hbar, 3) * sqr(constants.c));
	constants.aR = 4.0 * constants.sigma / constants.c;
	constants.Msol = cgs.Msol * units.g2code();
	constants.Rsol = cgs.Rsol * units.cm2code();
	constants.Lsol = cgs.Lsol * units.erg2code() / units.s2code();
	constants.AU = cgs.AU * units.cm2code();
	constants.pc = constants.AU * (360.0 * 60.0 * 60.0) / (2.0 * constants.pi);
	constants.ly = constants.c * (365.25 * 24.0 * 60.0 * 60.0 * units.s2code());
	return constants;
}


std::ostream& operator<<(std::ostream &os, Constants const &constants) {
	os << std::setprecision(10);
	os << "  c     = " << constants.c << "\n";
	os << "  G     = " << constants.G << "\n";
	os << "  kB    = " << constants.kB << "\n";
	os << "  h     = " << constants.h << "\n";
	os << "  hbar  = " << constants.hbar << "\n";
	os << "  sigma = " << constants.sigma << "\n";
	os << "  aR    = " << constants.aR << "\n";
	os << "  Msol  = " << constants.Msol << "\n";
	os << "  Rsol  = " << constants.Rsol << "\n";
	os << "  Lsol  = " << constants.Lsol << "\n";
	os << "  AU    = " << constants.AU << "\n";
	os << "  pc    = " << constants.pc << "\n";
	os << "  ly    = " << constants.ly << "\n";
	return os;
}

