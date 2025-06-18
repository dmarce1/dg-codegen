/******************************************************************************
 Copyright (C) 2024  Dominic C. Marcello
 *******************************************************************************/

#ifndef INCLUDE_CONSTANTS_HPP_
#define INCLUDE_CONSTANTS_HPP_

// physical_constants_cgs.hpp
#ifndef PHYSICAL_CONSTANTS_CGS_HPP
#define PHYSICAL_CONSTANTS_CGS_HPP

namespace cgs {
//constexpr double sigma = 1;
//constexpr double kB = 1;
//constexpr double c = 1;
//constexpr double amu = 1;
constexpr double sigma = 5.670374419e-5;
constexpr double kB = 1.380649e-16;
constexpr double c = 2.99792458e10;
constexpr double amu = 1.66053906660e-24;
constexpr double a = 4 * sigma / c;
constexpr double G = 6.67430e-8;
constexpr double h = 6.62607015e-27;
constexpr double hbar = 1.054571817e-27;
constexpr double Rsol = 6.957e10;
constexpr double Lsol = 3.828e33;
constexpr double Msol = 1.98847e33;
constexpr double H = 2.184e-18;
constexpr double Mpc = 3.08567758149137e24;
constexpr double ly = 9.4607304725808e17;
constexpr double minute = 60;
constexpr double hour = 60 * minute;
constexpr double day = 86400.0;
constexpr double year = 365.25 * day;
constexpr double mp = 1.67262192369e-24;
constexpr double me = 9.1093837015e-28;
constexpr double mH = mp + me;
}

#endif

#endif /* INCLUDE_CONSTANTS_HPP_ */
