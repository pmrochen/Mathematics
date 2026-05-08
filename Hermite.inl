/*
 *	Name: Hermite
 *	Author: Pawel Mrochen
 */

#pragma once

#include <concepts>

namespace mathematics::hermite {

template<typename V, std::floating_point T>
inline V basis(T t) noexcept
{
	T t2 = t*t;
	T t3 = t*t2;
	T h2 = T(3)*t2 - t3 - t3;
	T h4 = t3 - t2;
	return { T(1) - h2, h4 - t2 + t, h4, h2 };	// 2t^3 - 3t^2 + 1, t^3 - 2t^2 + t, t^3 - t^2, -2t^3 + 3t^2 
}

} // namespace mathematics::hermite
