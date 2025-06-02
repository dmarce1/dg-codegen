/******************************************************************************
 Copyright (C) 2024  Dominic C. Marcello
*******************************************************************************/

#ifndef INCLUDE_INDENT_HPP_
#define INCLUDE_INDENT_HPP_


#include <cassert>
#include <string>



struct Indention {
	static constexpr char tab_[] = "   ";
	Indention() :
			count_(0) {
	}
	Indention& operator++() {
		assert(count_ >= 0);
		count_++;
		return *this;
	}
	Indention& operator--() {
		count_--;
		assert(count_ >= 0);
		return *this;
	}
	Indention operator++(int) {
		auto rc = *this;
		operator++();
		return rc;
	}
	Indention operator--(int) {
		auto rc = *this;
		operator--();
		return rc;
	}
	operator std::string() const {
		std::string tabs;
		for (int i = 0; i < count_; i++) {
			tabs += tab_;
		}
		return tabs;
	}
private:
	int count_;
};


#endif /* INCLUDE_INDENT_HPP_ */
