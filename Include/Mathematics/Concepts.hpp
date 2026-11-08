/*
 *	Name: Concepts
 *	Author: Pawel Mrochen
 */

#pragma once

#include <concepts>

namespace mathematics {

struct Uninitialized {};
constexpr Uninitialized UNINITIALIZED{};

struct Identity {};
//constexpr Identity IDENTITY{};

struct Normalized {};
struct Unnormalized {};

template<typename T>
concept Normalization = (std::same_as<T, Normalized> || std::same_as<T, Unnormalized>);

template<typename T>
concept Arithmetic = (std::floating_point<T> || std::integral<T>);

} // namespace mathematics
