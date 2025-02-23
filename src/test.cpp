/*
 * test.cpp
 *
 *  Created on: Feb 19, 2025
 *      Author: dmarce1
 */

#include "Timer.hpp"
#include "LinearGravity.hpp"
#include "Interpolate.hpp"
#include "Real.hpp"
#include "Relativity.hpp"
#include "Options.hpp"
#include "Particles.hpp"

#include <unordered_map>

static std::unordered_map<std::string, Timer> timerDatabase;
static std::string currentTimer = "";

static std::string readTimers() {
	std::string output;
	output += std::string("Timing Results:\n");
	for (auto iter = timerDatabase.begin(); iter != timerDatabase.end(); iter++) {
		output += iter->first + std::string(" ") + std::to_string(iter->second.read());
	}
	return output;
}

static void resetTimers() {
	timerDatabase.clear();
}

static void startTimer(char const *name) {
	timerDatabase[name].start();
}

static void stopTimer() {
	timerDatabase[currentTimer].stop();
}

template<typename RealType>
struct CosmicGR {
	using real_type = RealType;
	static constexpr real_type outputTimeInterval = real_type(0.01);
	static constexpr real_type half = real_type(0.5), one = real_type(1);
	static constexpr real_type c = one;
	static constexpr real_type cfl = real_type(0.25);
	CosmicGR() :
			GR(getOptions().gridLength), DM(getOptions().particleCount, getOptions().gridLength, LinearGravity<Real>::BW), dx(one / getOptions().gridLength), time(
					0), lastOutputTime(-real_type(1.000001) * outputTimeInterval), frameCount(0) {
	}
	void output(SymmetricMatrix<std::valarray<real_type>, DIM4> const *Tptr = nullptr) const {
		using namespace Math;
		DBShowErrors(DB_ALL, siloErrorHandler);
		DBoptlist *optList = DBMakeOptlist(1);
		double thisTime = time;
		DBAddOption(optList, DBOPT_DTIME, &thisTime);
		std::string const filename = "X." + std::to_string(frameCount) + ".silo";
		DBfile *db = DBCreate(filename.c_str(), DB_CLOBBER, DB_LOCAL, "Astro-Tiger", DB_HDF5);
		GR.output(db, optList, Tptr);
		DM.output(db, optList);
		DBClose(db);
		DBFreeOptlist(optList);
		frameCount++;
	}
	void execute(real_type maxTime) {
		real_type const dt = cfl * dx / c;
		GR.enforceBoundaryConditions();
		auto const dgdx = [this](real_type x, real_type y, real_type z) {
			return GR.metricDerivatives(x, y, z);
		};
		int step = 0;
		while (time < maxTime) {
			printf("step = %i | time = %e | dt = %e\n", step, time, dt);
			if (time - lastOutputTime >= outputTimeInterval) {
				printf("outputing frame %i\n", frameCount);
				auto const T = DM.stressEnergyTensor();
				output(&T);
				lastOutputTime = time;
			}
			printf("kick\n");
			DM.kick(dgdx, half * dt);
			printf("drift\n");
			DM.drift(half * dt);
			printf("stress-energy tensor\n");
			auto const T = DM.stressEnergyTensor();
			printf("GR step\n");
			GR.step(T, dt);
			printf("BC\n");
			GR.enforceBoundaryConditions();
			printf("drift\n");
			DM.drift(half * dt);
			printf("kick\n");
			DM.kick(dgdx, half * dt);
			time += dt;
			step++;
		}
		output();
	}
private:
	LinearGravity<Real> GR;
	Particles<Real> DM;
	real_type dx;
	real_type time;
	real_type lastOutputTime;
	mutable int frameCount;
};

void test() {
	using namespace Math;
	CosmicGR<Real> test;
	test.execute(Real(1));
//TricubicSpline<double> test;
}
