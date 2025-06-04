/******************************************************************************
 Copyright (C) 2024  Dominic C. Marcello
 *******************************************************************************/

#ifndef INCLUDE_INDENT_HPP_
#define INCLUDE_INDENT_HPP_

#include <cassert>
#include <string>

struct Indent {
	static constexpr int spacesPerIndention_ = 3;
	Indent() {
	}
	Indent& operator++() {
		for (int n = 0; n < spacesPerIndention_; n++) {
			tabString_ += " ";
		}
		return *this;
	}
	Indent& operator--() {
		assert(!tabString_.empty());
		for (int n = 0; n < spacesPerIndention_; n++) {
			tabString_.pop_back();
		}
		return *this;
	}
	Indent operator++(int) {
		auto rc = *this;
		operator++();
		return rc;
	}
	Indent operator--(int) {
		auto rc = *this;
		operator--();
		return rc;
	}
	operator std::string() const {
		return tabString_;
	}
	operator char const*() const {
		return tabString_.c_str();
	}
	std::string operator+(std::string const &other) const {
		return std::string(tabString_) + other;
	}
private:
	std::string tabString_;
};

#endif /* INCLUDE_INDENT_HPP_ */
