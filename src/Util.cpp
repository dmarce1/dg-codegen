#include "Util.hpp"

#include <fstream>



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

