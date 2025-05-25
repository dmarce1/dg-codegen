#include "Util.hpp"

#include <cstring>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <stacktrace>

#include <fenv.h>
#include <signal.h>

void enableFPE() {
	feclearexcept(FE_ALL_EXCEPT);
	feenableexcept(FE_INVALID | FE_DIVBYZERO | FE_OVERFLOW);
}

void disableFPE() {
	feclearexcept(FE_ALL_EXCEPT);
	fedisableexcept(FE_INVALID | FE_DIVBYZERO | FE_OVERFLOW);
}

bool writeList(std::string const &filename, std::string const &header, std::string const &item) {
	namespace fs = std::filesystem;
	bool created;
	std::ofstream file;
	if (!fs::exists(filename)) {
		created = true;
		file.open(filename, std::ios::out);
		if (!file) {
			throw std::runtime_error("Failed to create file: " + filename);
		}
		file << header << '\n';
	} else {
		created = false;
		file.open(filename, std::ios::app);
		if (!file) {
			throw std::runtime_error("Failed to open file for appending: " + filename);
		}
	}
	file << item << std::endl;
	return created;
}

#ifndef NDEBUG

void fpeHandler(int, siginfo_t*, void*) {
	std::cerr << "SIGFPE!\n" << std::stacktrace::current() << std::endl;
	std::_Exit(1);
}

void installFpeHandler() {
	struct sigaction sa;
	memset(&sa, 0, sizeof(sa));
	sa.sa_flags = SA_SIGINFO;
	sa.sa_sigaction = fpeHandler;
	sigemptyset(&sa.sa_mask);
	if (sigaction(SIGFPE, &sa, nullptr) != 0) {
		perror("sigaction");
		std::exit(1);
	}
}

#endif
