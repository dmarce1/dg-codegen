/*
 * test.cpp
 *
 *  Created on: Feb 19, 2025
 *      Author: dmarce1
 */

#include "LinearGravity.hpp"
#include "Interpolate.hpp"
#include "Real.hpp"
#include "Relativity.hpp"
#include "Options.hpp"
#include "Particles.hpp"

struct CosmicGR {
	CosmicGR() :
			GR(getOptions().gridLength), DM(getOptions().particleCount), frameCount(0) {
	}
	void output() const {
		using namespace Math;
		DBShowErrors(DB_ALL, siloErrorHandler);
		DBoptlist *optList = DBMakeOptlist(0);
		std::string const filename = "X." + std::to_string(frameCount) + ".silo";
		DBfile *db = DBCreate(filename.c_str(), DB_CLOBBER, DB_LOCAL, "Astro-Tiger", DB_HDF5);
		GR.output(db, optList);
		DBClose(db);
		DBFreeOptlist(optList);
		frameCount++;
	}
private:
	LinearGravity GR;
	Particles DM;
	mutable int frameCount;
};

void test() {
	using namespace Math;
	LinearGravity test(32);
//TricubicSpline<double> test;
}
