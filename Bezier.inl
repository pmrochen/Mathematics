/*
 *	Name: Bezier
 *	Author: Pawel Mrochen
 */

#pragma once

#include <concepts>
#include <Simd/Intrinsics.hpp>

namespace mathematics::bezier {

template<std::floating_point T>
inline T evaluate(T x0, T x1, T x2, T x3, T t) noexcept
{
	T t2 = t*t;
	T t3 = t2*t;
	T c = T(3)*(x1 - x0);
	T b = T(3)*(x2 - x1) - c;
	T a = x3 - x0 - c - b;
	return a*t3 + b*t2 + c*t + x0;
}

#if SIMD_HAS_FLOAT4

template<>
inline float evaluate(float x0, float x1, float x2, float x3, float t) noexcept
{
	float t2 = t*t;
	float t3 = t2*t;
	float c = 3.f*(x1 - x0);
	float b = 3.f*(x2 - x1) - c;
	float a = x3 - x0 - c - b;
	return simd::dot4(simd::float4(a, b, c, x0), simd::float4(t3, t2, t, 1.f));
}

#endif /* SIMD_HAS_FLOAT4 */

template<typename V, std::floating_point T>
inline V basis(T t) noexcept
{
	T t2 = t*t;
	T omt = T(1) - t;
	T omt2 = omt*omt;
	return { omt*omt2, T(3)*t*omt2, T(3)*t2*omt, t*t2 };
}

template<typename V, std::floating_point T>
inline V derivativeBasis(T t) noexcept
{
	T t2 = t*t;
	T omt = T(1) - t;
	T omt2 = omt*omt;
	return { T(-3)*omt2, T(-6)*t*omt + T(3)*omt2, T(-3)*t2 + T(6)*t*omt, T(3)*t2 };
}

} // namespace mathematics::bezier
