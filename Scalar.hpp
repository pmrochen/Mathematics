/*
 *	Name: Scalar
 *	Author: Pawel Mrochen
 */

#pragma once

#include <concepts>
#include <utility>
#include <cmath>
#include "Constants.hpp"

#ifndef MATHEMATICS_API
#if defined(MATHEMATICS_DLL) || defined(CORE_DLL)
#if defined(MATHEMATICS_EXPORTS) || defined(CORE_EXPORTS)
#define MATHEMATICS_API __declspec(dllexport)
#else
#define MATHEMATICS_API __declspec(dllimport)
#endif
#else
#define MATHEMATICS_API
#endif
#endif

#ifndef MATHEMATICS_FAST_NORMALIZE
#define MATHEMATICS_FAST_NORMALIZE 1
#endif

namespace core::mathematics {

template<typename T>
	requires (std::floating_point<T> || std::integral<T>)
inline T sign(T x)
{ 
	return (x > T(0)) ? T(1) : ((x < T(0)) ? T(-1) : T(0));
}

template<std::floating_point T>
inline T frac(T x)
{ 
	return x - std::floor(x);
}

template<std::floating_point T>
inline bool approxEquals(T a, T b)
{ 
	return (std::fabs(b - a) < Constants<T>::TOLERANCE);
}

template<std::floating_point T>
inline bool approxEquals(T a, T b, T tolerance)
{ 
	return (std::fabs(b - a) < tolerance); 
}

//template<std::floating_point T>
//inline bool isNan(T x) 
//{ 
//	return !(x == x); 
//}

template<std::floating_point T>
inline T log(T x, T b) // logarithm for a base b
{ 
	return std::log(x)/std::log(b);
}

template<std::floating_point T>
inline T log2(T x) // base-2 logarithm
{ 
	return std::log(x)/Constants<T>::LN2; 
}

template<std::floating_point T>
inline T radians(T x)
{ 
	return x*Constants<T>::DEG_TO_RAD; 
}

template<std::floating_point T>
inline T degrees(T x)
{ 
	return x*Constants<T>::RAD_TO_DEG; 
}

template<typename T>
	requires (std::floating_point<T> || std::integral<T>)
inline T step(T a, T t)
{ 
	return (t >= a) ? T(1) : T(0); 
}

template<typename T>
	requires (std::floating_point<T> || std::integral<T>)
inline T pulse(T a, T b, T t)
{ 
	return step(a, t) - step(b, t); 
}

template<std::floating_point T>
inline T boxStep(T a, T b, T t)
{
	if (t <= a)
		return T(0);
	if (t >= b)
		return T(1);
	return (t - a)/(b - a);
}

template<std::floating_point T>
inline T smoothStep(T a, T b, T t)
{
	if (t <= a)
		return T(0);
	if (t >= b)
		return T(1);
	t = (t - a)/(b - a);
	return t*t*(T(3) - T(2)*t);
}

template<std::floating_point T>
inline T smootherStep(T a, T b, T t)
{
	if (t <= a)
		return T(0);
	if (t >= b)
		return T(1);
	t = (t - a)/(b - a);
	return t*t*t*(t*(t*T(6) - T(15)) + T(10));
}

} // namespace core::mathematics
