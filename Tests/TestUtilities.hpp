/*
 *	Name: TestUtilities
 *	Author: Pawel Mrochen
 */

#pragma once

#include <doctest/doctest.h>
#include <type_traits>
#include <concepts>
#include <algorithm>
#include <sstream>
#include <string>
#include <cmath>
#include <Mathematics/Mathematics.hpp>

// Every floating-point test case is instantiated for both precisions. With SIMD enabled, float exercises the
// vectorized specializations while double exercises the generic scalar code paths.
#define TEST_REAL_TYPES float, double

namespace test {

using namespace mathematics;
using namespace mathematics::templates;

// Tolerance suitable for values of magnitude ~1 after a handful of arithmetic operations.
template<std::floating_point T>
constexpr T tolerance() noexcept
{
	if constexpr (std::is_same_v<T, float>)
		return T(1e-4);
	else
		return T(1e-9);
}

template<typename T>
struct RealOf { using Type = typename T::Real; };

template<std::floating_point T>
struct RealOf<T> { using Type = T; };

template<typename T>
using Real = typename RealOf<T>::Type;

template<std::floating_point T>
bool approxEqual(T actual, T expected, T tol = tolerance<T>()) noexcept
{
	if (std::isinf(expected) || std::isinf(actual))
		return actual == expected;
	return std::fabs(actual - expected) <= tol*std::max(T(1), std::max(std::fabs(actual), std::fabs(expected)));
}

template<typename V>
	requires (!std::floating_point<V>)
bool approxEqual(const V& actual, const V& expected, Real<V> tol = tolerance<Real<V>>()) noexcept
{
	return actual.approxEquals(expected, tol);
}

template<typename V>
std::string toString(const V& value)
{
	std::ostringstream s;
	s.precision(10);
	s << value;
	return s.str();
}

} // namespace test

// Approximate comparison of scalars or any library type providing approxEquals(other, tolerance) and operator<<.
#define CHECK_APPROX(actual, expected) \
	CHECK_MESSAGE(::test::approxEqual((actual), (expected)), \
		"actual:   ", ::test::toString(actual), "\nexpected: ", ::test::toString(expected))

#define CHECK_APPROX_TOL(actual, expected, tol) \
	CHECK_MESSAGE(::test::approxEqual((actual), (expected), (tol)), \
		"actual:   ", ::test::toString(actual), "\nexpected: ", ::test::toString(expected))
