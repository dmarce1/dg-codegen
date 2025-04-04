/*
 * Index.hpp
 *
 *  Created on: Apr 3, 2025
 *      Author: dmarce1
 */

#ifndef INCLUDE_INDEX_HPP_
#define INCLUDE_INDEX_HPP_

#include <iostream>
#include <tuple>
#include <utility>
#include <string>
#include <type_traits>
#include <array>

namespace Tensors {

template<char C>
struct Index {
	static constexpr char value = C;
	constexpr operator char() const {
		return C;
	}
};

template<typename T>
struct IsIndex {
	static constexpr bool value = false;
};

template<char C>
struct IsIndex<Index<C>> {
	static constexpr bool value = true;
	operator char() const {
		return value;
	}
};

template<typename ... Args>
constexpr auto extractChars() {
	std::array<char, (size_t(IsIndex<Args>::value) + ...)> chars;
	size_t i = 0;
	((IsIndex<Args>::value ? (chars[i++] = Args(), void()) : void()), ...);
	return chars;
}

// Make a getter for each argument
template<typename T>
auto make_getter(T &&val) {
	if constexpr (IsIndex<std::decay_t<T>>::value) {
		constexpr size_t pos = T::value - 'a'; // maps 'a' to 0, etc.
		return [](auto &&... args) {
			return std::get < pos > (std::forward_as_tuple(args...));
		};
	} else {
		return [v = std::forward<T>(val)](auto&&...) {
			return v;
		};
	}
}

template<typename F, typename ... Args>
auto indexBind(F f, Args ... args) {
	constexpr auto chars = extractChars<Args...>();
	constexpr size_t N = chars.size();

	// Create lambda that instantiates f<chars...>
	auto inst = [&]<std::size_t... Is>(std::index_sequence<Is...>) {
		return [=]<typename... CallArgs>(CallArgs&&... call_args) {
			return f.template operator()<chars[Is]...>(std::forward<CallArgs>(call_args)...);
		};
	}(std::make_index_sequence<N> { });

	// Build runtime argument getters
	auto getters = std::make_tuple(make_getter(std::forward<Args>(args))...);

	// Return bound callable
	return [inst, getters](auto &&... call_args) mutable {
		return std::apply([&](auto &&... gs) {
			return inst(gs(std::forward<decltype(call_args)>(call_args)...)...);
		},
		getters);
	};
}

}

#endif /* INCLUDE_INDEX_HPP_ */
