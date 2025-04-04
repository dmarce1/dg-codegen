#pragma once

#include <concepts>
#include <string>
#include <cstddef>
#include <type_traits>

template<typename T>
concept Indexable = requires(T a, size_t i) {
	{a[i]}-> std::convertible_to<typename T::value_type>;
	{a.size()}-> std::convertible_to<size_t>;
};

template<Indexable Container>
std::string toString(Container const &p) {
	size_t const N = p.size();
	std::string str = "(";
	for (size_t n = 0; n < N; n++) {
		str += std::to_string(p[n]);
		str += (n + 1 < N) ? "," : "";
	}
	str += ")";
	return str;
}
