/*
 * Constants.hpp
 *
 *  Created on: Jan 12, 2025
 *      Author: dmarce1
 */

#ifndef INCLUDE_CONSTANTS123_HPP_
#define INCLUDE_CONSTANTS123_HPP_

#include "Util.hpp"

#include <numbers>
#include <ostream>

class Units;

struct Constants {
	static constexpr double pi = 3.14159265358979323846;
	double c;
	double G;
	double u;
	double me;
	double kB;
	double h;
	double hbar;
	double sigma;
	double aR;
	double Msol;
	double Rsol;
	double Lsol;
	double AU;
	double pc;
	double ly;
};

void setLengthScale(double);
Constants getCodeConstants();
Constants getCgsConstants();
Units const& getUnits();

class Units {
	double code2cm_;
	double code2g_;
	double code2s_;
	double code2K_;
public:
	Units(double = 1.0);
	double code2cm() const;
	double code2s() const;
	double code2g() const;
	double code2K() const;
	double code2dyne() const;
	double code2erg() const;
	double code2Hz() const;
	double cm2code() const;
	double s2code() const;
	double g2code() const;
	double K2code() const;
	double dyne2code() const;
	double erg2code() const;
	double Hz2code() const;
	friend std::ostream& operator<<(std::ostream &os, Units const &u);
};

std::ostream& operator<<(std::ostream &, Constants const &);


#endif /* INCLUDE_CONSTANTS_HPP_ */
