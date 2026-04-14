/*
 *	Name: Bezier
 *	Author: Pawel Mrochen
 */

#pragma once

#include <concepts>
#include "Vector4.hpp"

namespace mathematics::bezier {

template<std::floating_point T>
inline Vector4<T> basis(T t)
{
	T t2 = t*t;
	T omt = T(1) - t;
	T omt2 = omt*omt;
	return Vector4<T>(omt*omt2, T(3)*t*omt2, T(3)*t2*omt, t*t2);
}

template<std::floating_point T>
inline Vector4<T> derivativeBasis(T t)
{
	T t2 = t*t;
	T omt = T(1) - t;
	T omt2 = omt*omt;
	return Vector4<T>(T(-3)*omt2, T(-6)*t*omt + T(3)*omt2, T(-3)*t2 + T(6)*t*omt, T(3)*t2);
}

} // namespace mathematics::bezier
