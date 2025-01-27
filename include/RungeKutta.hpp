/*
 * RungeKutta.hpp
 *
 *  Created on: Jan 25, 2025
 *      Author: dmarce1
 */

#ifndef INCLUDE_RUNGEKUTTA_HPP_
#define INCLUDE_RUNGEKUTTA_HPP_

template<typename T, int O>
struct RungeKutta;

template<typename T>
struct RungeKutta<T, 1> {
	static constexpr int s = 1;
	static constexpr std::array<T, s> a = { { } };
	static constexpr std::array<T, s> b = { T(1.00000000000000) };
	static constexpr std::array<T, s> c = { T(0.00000000000000) };
};

template<typename T>
struct RungeKutta<T, 2> {
	static constexpr int s = 2;
	static constexpr std::array<T, s> a = { { }, //
			{ T(1.00000000000000) } //
	};
	static constexpr std::array<T, s> b = { T(0.50000000000000), T(0.50000000000000) };
	static constexpr std::array<T, s> c = { T(0.00000000000000), T(1.00000000000000) };
};

template<typename T>
struct RungeKutta<T, 3> {
	static constexpr int s = 3;
	static constexpr std::array<T, s> a = { { }, //
			{ T(1.00000000000000) }, //
			{ T(0.25000000000000), T(0.25000000000000) } }; //
	static constexpr std::array<T, s> b = { T(0.16666666666667), T(0.16666666666667), T(0.66666666666667) }; //
	static constexpr std::array<T, s> c = { T(0.00000000000000), T(1.00000000000000), T(0.50000000000000) };
};

template<typename T>
struct RungeKutta<T, 4> {
	static constexpr int s = 5;
	static constexpr std::array<T, s> a = { { }, //
			{ T(0.39175222700392) }, //
			{ T(0.21766909633821), T(0.36841059262959) }, //
			{ T(0.08269208670950), T(0.13995850206999), T(0.25189177424738) }, //
			{ T(0.06796628370320), T(0.11503469844438), T(0.20703489864929), T(0.50000000000000) } //
	};
	static constexpr std::array<T, s> b = { T(0.14681187618661), T(0.24848290924556), T(0.10425883036650), T(0.27443890091960), T(0.22600748319395) };
	static constexpr std::array<T, s> c = { T(0.00000000000000), T(0.39175222700392), T(0.58607968896779), T(0.47454236302687), T(0.93501063100924) };
};

template<typename T>
struct RungeKutta<T, 5> {
	static constexpr int s = 7;
	static constexpr std::array<T, s> a = { { }, //
			{ T(+0.392382208054010) }, //
			{ T(+0.310348765296963), T(+0.523846724909595) }, //
			{ T(+0.114817342432177), T(+0.248293597111781), T(+0.000000000000000) }, //
			{ T(+0.136041285050893), T(+0.163250087363657), T(+0.000000000000000), T(+0.557898557725281) }, //
			{ T(+0.135252145083336), T(+0.207274083097540), T(-0.180995372278096), T(+0.326486467604174), T(+0.348595427190109) }, //
			{ T(+0.082675687408986), T(+0.146472328858960), T(-0.160507707995237), T(+0.161924299217425), T(+0.028864227879979), T(+0.070259587451358) }  //
	};
	static constexpr std::array<T, s> b = //
			{ T(+0.110184169931401), T(+0.122082833871843), T(-0.117309105328437), T(+0.169714358772186), T(+0.143346980044187), T(+0.348926696469455), T(
					+0.223054066239366) }; //
	static constexpr std::array<T, s> c = //
			{ T(+0.000000000000000), T(+0.392382208054010), T(+0.310348765296963), T(+0.523846724909595), T(+0.114817342432177), T(+0.248293597111781), T(
					+0.557898557725281) }; //
};

#endif /* INCLUDE_RUNGEKUTTA_HPP_ */
