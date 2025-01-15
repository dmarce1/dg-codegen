/*
 * Options.hpp
 *
 *  Created on: Jan 13, 2025
 *      Author: dmarce1
 */

#ifndef INCLUDE_OPTIONS_HPP_
#define INCLUDE_OPTIONS_HPP_

#include <string>

struct Options {
	std::string configFile;
	double dualEnergyPressureSwitch;
	double dualEnergyUpdateSwitch;
	double fluidGamma;
	double gridScale;
	int gridLength;
};

const Options& getOptions();
bool processOptions(int, char*[]);

#endif /* INCLUDE_OPTIONS_HPP_ */
