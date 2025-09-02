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

void toFile(std::string const &content, std::filesystem::path const &filePath) {
	namespace fs = std::filesystem;
	auto parentPath = filePath.parent_path();
	if (!parentPath.empty() && !fs::exists(parentPath)) {
		fs::create_directories(parentPath);
	}
	std::ofstream ofs(filePath, std::ios::out | std::ios::trunc);
	if (!ofs) {
		throw std::runtime_error("Failed to open file: " + filePath.string());
	}
	ofs << content;
}

#include <mutex>

static void fpeHandler(int, siginfo_t* info, void*) noexcept {
    const char* reason = "Unknown FPE";
    switch (info->si_code) {
        case FPE_INTDIV: reason = "Integer divide by zero"; break;
        case FPE_INTOVF: reason = "Integer overflow"; break;
        case FPE_FLTDIV: reason = "Floating-point divide by zero"; break;
        case FPE_FLTOVF: reason = "Floating-point overflow"; break;
        case FPE_FLTUND: reason = "Floating-point underflow"; break;
        case FPE_FLTRES: reason = "Floating-point inexact result"; break;
        case FPE_FLTINV: reason = "Invalid floating-point operation"; break;
        case FPE_FLTSUB: reason = "Subscript out of range"; break;
    }

    // Write to stderr using async-signal-safe functions
    (void)!write(STDERR_FILENO, reason, strlen(reason));
    (void)!write(STDERR_FILENO, "\n", 1);

    _Exit(EXIT_FAILURE);
}
void installFpeHandler() {
	static std::once_flag once;
	std::call_once(once, [] {
		struct sigaction sa { };
		sa.sa_flags = SA_SIGINFO;
		sa.sa_sigaction = fpeHandler;
		sigemptyset(&sa.sa_mask);
		if (sigaction(SIGFPE, &sa, nullptr) != 0) {
			perror("sigaction");
			std::exit(1);
		}
	});
}

void enableFpeTrapsThisThread() {
	sigset_t set;
	sigemptyset(&set);
	sigaddset(&set, SIGFPE);
	pthread_sigmask(SIG_UNBLOCK, &set, nullptr);
	feenableexcept(FE_DIVBYZERO | FE_INVALID | FE_OVERFLOW);
}


