/*
 * Timer.hpp
 *
 *  Created on: Feb 23, 2025
 *      Author: dmarce1
 */

#ifndef INCLUDE_TIMER_HPP_
#define INCLUDE_TIMER_HPP_

#include <chrono>
#include <stack>

struct Timer {
	using time_type = std::chrono::time_point<std::chrono::high_resolution_clock>;
	Timer() :
			running(false), elapsedTime(0.0) {
	}
	bool isRunning() const {
		return running;
	}
	double read() {
		double readTime;
		if (running) {
			stop();
			readTime = elapsedTime;
			start();
		} else {
			readTime = elapsedTime;
		}
		return readTime;
	}
	double reset() {
		if (running) {
			stop();
		}
		double const readTime = elapsedTime;
		elapsedTime = 0.0;
		return readTime;
	}
	void start() {
		startTime = readSystemClock();
		running = true;
	}
	void stop() {
		time_type const stopTime = readSystemClock();
		running = false;
		elapsedTime += std::chrono::duration<double>(stopTime - startTime).count();
	}
private:
	time_type readSystemClock() {
		return std::chrono::high_resolution_clock::now();
	}
	bool running;
	double elapsedTime;
	time_type startTime;
};

#endif /* INCLUDE_TIMER_HPP_ */
