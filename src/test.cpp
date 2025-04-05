/*
 * test.cpp
 *
 *  Created on: Feb 19, 2025
 *      Author: dmarce1
 */

#include "Timer.hpp"
#include "LinearGravity.hpp"
#include "Interpolate.hpp"
#include "GRGrid.hpp"
#include "Real.hpp"
#include "Relativity.hpp"
#include "Options.hpp"
#include "Particles.hpp"
#include "Polynomial.hpp"
#include "Tensor.hpp"

#include <unordered_map>

static std::unordered_map<std::string, Timer> timerDatabase;
static std::string currentTimer = "";

static std::string readTimers() {
	std::string output;
	output += std::string("Timing Results:\n");
	for (auto iter = timerDatabase.begin(); iter != timerDatabase.end(); iter++) {
		output += iter->first + std::string(" ") + std::to_string(iter->second.read()) + " s\n";
	}
	return output;
}


static void startTimer(char const *name) {
	timerDatabase[name].start();
	currentTimer = name;
}

static void stopTimer() {
	timerDatabase[currentTimer].stop();
	currentTimer = "";
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
		auto const opts = getOptions();
		real_type const dt = cfl * dx / c;
		GR.enforceBoundaryConditions();
//		auto const grVars = [this](int x, int  y, int  z) {
//			return GR.getStateVars(x, y, z);
//		};
		auto const grVars = [this](real_type x, real_type y, real_type z) {
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

			startTimer("sort");
			printf("sort\n");
			DM.sort_particles();
			stopTimer();

			startTimer("kick");
			printf("kick\n");
			DM.kick(grVars, half * dt);
			stopTimer();

			startTimer("drift");
			printf("drift\n");
			DM.drift(half * dt);
			stopTimer();

			startTimer("stress-energy tensor");
			printf("stress-energy tensor\n");
			auto const T = DM.stressEnergyTensor();
			stopTimer();

			startTimer("GR step");
			printf("GR step\n");
			GR.step(T, dt);
			stopTimer();

			startTimer("BC");
			printf("BC\n");
			GR.enforceBoundaryConditions();
			stopTimer();

			startTimer("drift");
			printf("drift\n");
			DM.drift(half * dt);
			stopTimer();

			startTimer("kick");
			printf("kick\n");
			DM.kick(grVars, half * dt);
			stopTimer();

			time += dt;
			step++;
			std::cout << "\n\n" << readTimers() << "\n\n";
		}
		output();
	}
private:
	LinearGravity<real_type> GR;
	Particles<real_type> DM;
	real_type dx;
	real_type time;
	real_type lastOutputTime;
	mutable int frameCount;
};

void testEinstein();

// a0 x^0 + a1 x^1 + a2 x ^2 + ... + an x^n

// a0 x^0 + ... + a(n/2) x^(2*n)
// b0 x^0 + ... + b(n/2) x^(2*n)

void FPT(std::vector<double> &Y, int N) {
	if (N == 2) {
		Y[1] -= Y[0];
		return;
	}
	int const N2 = N / 2;
	std::vector<double> Ye(N2), Yo(N2);
	for (int n2 = 0; n2 < N2; n2++) {
		Ye[n2] = Y[2 * n2];
		Yo[n2] = Y[2 * n2 + 1];
	}
	FPT(Ye, N2);
	FPT(Yo, N2);

}

using sparse_polynomial = std::vector<std::pair<int, double>>;

using polynomial = std::vector<double>;

polynomial polynomialReduce(polynomial const &D, sparse_polynomial const &I) {
	int const N = D.size() - 1;
	int const M = I.size() - 1;
	polynomial Q(N - M + 1, 0.0);
	polynomial R = D;
	for (int n = N; n >= M; n--) {
		for (int m = 1; m <= M; m++) {
			R[n + I[m].first - M] += R[n] * I[m].second;
		}
		Q[n - M] = R[n];
		R[n] = 0.0;
	}
	return Q;
}


void test() {
	for(size_t d = 0; d < 5; d++) {
		Polynomial<double> P;
		P[d] = 1.0;
		std::cout << toString(P) << " ";
		P = polynomialDiscreteAntiDerivative(P, d);
		std::cout << toString(P) << "\n";
	}
	//	auto DD = A(i) * B(j);
}

/***/
/***/
/***/
/***/
/***/
/***/
/***/
/***/
/***/
/***/
/***/
/***/
/***/
/***/
/***/
/***/
/***/
