///******************************************************************************
// Copyright (C) 2024  Dominic C. Marcello
// *******************************************************************************/
//
//#ifndef SRC_SYMBOLIC_CPP_
//#define SRC_SYMBOLIC_CPP_
//
//#include <array>
//#include <string>
//#include <map>
//#include <set>
//#include <vector>
//
//inline void hashCombine(std::size_t &seed, std::size_t value) {
//	seed ^= value + 0x9e3779b97f4a7c15ULL + (seed << 6) + (seed >> 2);
//}
//
//struct FunctionType {
//	FunctionType() {
//	}
//	FunctionType(std::string n, std::vector<std::string> p) {
//		name = n;
//		parameters = p;
//	}
//	std::string name;
//	std::vector<std::string> parameters;
//	bool operator==(FunctionType other) const {
//		if (other.name != name) {
//			return false;
//		}
//		if (other.parameters != parameters) {
//			return false;
//		}
//		return true;
//	}
//	bool operator<(FunctionType const &other) const {
//		return name < other.name;
//	}
//};
//
//struct Derivatives {
//	Derivatives(std::vector<FunctionType> const &funcs) {
//		for (auto const &f : funcs) {
//			deps_[f.name] = f;
//		}
//		for (auto const &f : funcs) {
//			vars_.merge(findIndependents(f));
//		}
//		indexCount = 0;
//		for (auto const &dep : deps_) {
//			std::pair<int, int> rng;
//			rng.first = indexCount;
//			indexCount += dep.second.parameters.size() + 1;
//			rng.second = indexCount;
//			index_[dep.first] = rng;
//		}
//
//	}
//	std::set<std::string> findIndependents(FunctionType f) const {
//		std::set<std::string> indes;
//		if (f.parameters.size() == 0) {
//			indes.insert(f.name);
//		} else {
//			for (auto const &p : f.parameters) {
//				FunctionType const g = deps_.find(p)->second;
//				indes.merge(findIndependents(g));
//			}
//		}
//		return indes;
//	}
//	Derivatives differentiate(int k) const {
//		if (k < vars_.size()) {
//			Derivatives dydz;
//		} else {
//			return *this;
//		}
//	}
//	int indexCount;
//	std::map<std::string, std::pair<int, int>> index_;
//	std::map<std::string, FunctionType> deps_;
//	std::set<std::string> vars_;
//};
//
//void testD() {
//	std::vector<FunctionType> f;
//	Derivatives D(f);
//}
//#endif /* SRC_SYMBOLIC_CPP_ */
