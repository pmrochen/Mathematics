/*
 *	Name: Scalar
 *	Author: Pawel Mrochen
 */

#pragma once

#include <limits>
#include <type_traits>
#include <concepts>
#include <utility>
#include <bit>
#include <cmath>
#include "Constants.hpp"

namespace mathematics {

template<typename T>
	requires (std::floating_point<T> || std::integral<T>)
inline T abs(T x) noexcept
{
	if constexpr (std::is_floating_point_v<T>)
	{
		return std::fabs(x);
	}
	else if constexpr (std::is_signed_v<T>)
	{
#if defined(__GNUC__) || defined(__clang__)
		if constexpr (std::is_same_v<T, int>)
			return __builtin_abs(x);
		else if constexpr (std::is_same_v<T, long>)
			return __builtin_labs(x);
		else if constexpr (std::is_same_v<T, long long>)
			return __builtin_llabs(x);
#endif
		T y = x >> std::numeric_limits<T>::digits;
		return (x + y) ^ y;
	}
	else // unsigned integer
	{
		return x;
	}
}

template<typename T>
	requires (std::floating_point<T> || std::integral<T>)
inline T sign(T x) noexcept
{ 
	return (x > T(0)) ? T(1) : ((x < T(0)) ? T(-1) : T(0));
}

template<std::floating_point T>
inline T frac(T x) noexcept
{ 
	return x - std::floor(x);
}

template<typename T>
	requires (std::floating_point<T> || std::integral<T>)
inline T sqr(T x) noexcept
{
	return x*x;
}

template<std::floating_point T>
inline bool approxEquals(T a, T b) noexcept
{ 
	return (std::fabs(b - a) < Constants<T>::TOLERANCE);
}

template<std::floating_point T>
inline bool approxEquals(T a, T b, T tolerance) noexcept
{ 
	return (std::fabs(b - a) < tolerance); 
}

//template<std::floating_point T>
//inline bool isNan(T x) noexcept
//{ 
//	return !(x == x); 
//}

template<std::floating_point T>
inline T log(T x, T b) noexcept // logarithm for a base b
{ 
	return std::log(x)/std::log(b);
}

template<typename T>
	requires (std::floating_point<T> || std::integral<T>)
inline T log2(T x) noexcept // base-2 logarithm
{ 
	if constexpr (std::is_floating_point_v<T>)
	{
		return std::log(x)/Constants<T>::LN2;
	}
	else
	{
		if constexpr (std::is_signed_v<T>)
			x -= x & (x >> std::numeric_limits<T>::digits); // branchless max(x, 0)
#if defined(__GNUC__) || defined(__clang__)
		if constexpr (std::is_same_v<T, unsigned int> || std::is_same_v<T, int>)
			return T(sizeof(T)*8 - 1 - __builtin_clz((unsigned int)x));
		else if constexpr (std::is_same_v<T, unsigned long> || std::is_same_v<T, long>)
			return T(sizeof(T)*8 - 1 - __builtin_clzl((unsigned long)x));
		else if constexpr (std::is_same_v<T, unsigned long long> || std::is_same_v<T, long long>)
			return T(sizeof(T)*8 - 1 - __builtin_clzll((unsigned long long)x));
#endif
		return T(sizeof(T)*8 - 1 - std::countl_zero(std::make_unsigned_t<T>(x)));
	}
}

template<std::floating_point T>
inline T radians(T x) noexcept
{ 
	return x*Constants<T>::DEG_TO_RAD; 
}

template<std::floating_point T>
inline T degrees(T x) noexcept
{ 
	return x*Constants<T>::RAD_TO_DEG; 
}

template<std::floating_point T>
inline T lerp(T a, T b, T t) noexcept
{ 
	return a + t*(b - a);
}

template<typename T>
	requires (std::floating_point<T> || std::integral<T>)
inline T step(T a, T t) noexcept
{ 
	return (t >= a) ? T(1) : T(0); 
}

template<typename T>
	requires (std::floating_point<T> || std::integral<T>)
inline T pulse(T a, T b, T t) noexcept
{ 
	return step(a, t) - step(b, t); 
}

template<std::floating_point T>
inline T boxStep(T a, T b, T t) noexcept
{
	if (t <= a)
		return T(0);
	if (t >= b)
		return T(1);
	return (t - a)/(b - a);
}

template<std::floating_point T>
inline T smoothStep(T a, T b, T t) noexcept
{
	if (t <= a)
		return T(0);
	if (t >= b)
		return T(1);
	t = (t - a)/(b - a);
	return t*t*(T(3) - T(2)*t);
}

template<std::floating_point T>
inline T smootherStep(T a, T b, T t) noexcept
{
	if (t <= a)
		return T(0);
	if (t >= b)
		return T(1);
	t = (t - a)/(b - a);
	return t*t*t*(t*(t*T(6) - T(15)) + T(10));
}

} // namespace mathematics
