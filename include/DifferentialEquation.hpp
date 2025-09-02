/******************************************************************************
 Copyright (C) 2024  Dominic C. Marcello
 *******************************************************************************/

#ifndef INCLUDE_DIFFERENTIALEQUATION_HPP_
#define INCLUDE_DIFFERENTIALEQUATION_HPP_

#include <memory>
#include <string>
#include <vector>

//struct Variable {
//	Variable(std::string const &name, int dimCount = 1) {
//		name_ = name;
//		exponent_ = 0;
//		dimCount_ = 1;
//		derivativeOrder_.resize(dimCount_, 0);
//		dependencies_.resize(dimCount_, nullptr);
//	}
//	std::vector<Variable> differentiate(int dimension) const {
//		std::vector<Variable> factors;
//		Variable F = *this;
//		F.derivativeOrder_[dimension]++;
//		factors.push_back(F);
//		if (dependencies_[dimension] != nullptr) {
//			auto const U = *dependencies_[dimension];
//			auto const ufactors = U.differentiate(dimension);
//			for(auto const& factor : ufactors) {
//				factors.push_back(factor);
//			}
//		}
//		return factors;
//	}
//	friend std::vector<std::vector<Variable>> differentiate(std::vector<Variable> const& factors, int dimension) const {
//		return factors;
//	}
//	std::string name_;
//	std::vector<int> derivativeOrder_;
//	std::vector<std::shared_ptr<Variable>> dependencies_;
//	int exponent_;
//	int dimCount_;
//};
//
#endif /* INCLUDE_DIFFERENTIALEQUATION_HPP_ */
