#include "TestUtilities.hpp"

using namespace mathematics;

TEST_SUITE("Scalar")
{

TEST_CASE_TEMPLATE("Constants are consistent", T, TEST_REAL_TYPES)
{
	using C = Constants<T>;
	CHECK_APPROX(C::PI, T(std::acos(-1.0)));
	CHECK_APPROX(C::HALF_PI, C::PI/T(2));
	CHECK_APPROX(C::TWO_PI, C::PI*T(2));
	CHECK_APPROX(C::INV_PI, T(1)/C::PI);
	CHECK_APPROX(C::PI_SQR, C::PI*C::PI);
	CHECK_APPROX(C::E, T(std::exp(1.0)));
	CHECK_APPROX(C::LN2, T(std::log(2.0)));
	CHECK_APPROX(C::SQRT_2, T(std::sqrt(2.0)));
	CHECK_APPROX(C::SQRT_3, T(std::sqrt(3.0)));
	CHECK_APPROX(C::DEG_TO_RAD*C::RAD_TO_DEG, T(1));
	CHECK(C::INF == std::numeric_limits<T>::infinity());
	CHECK(C::MINUS_INF == -std::numeric_limits<T>::infinity());
	CHECK(C::EPSILON == std::numeric_limits<T>::epsilon());
	CHECK(C::TOLERANCE > T(0));
}

TEST_CASE_TEMPLATE("abs", T, TEST_REAL_TYPES)
{
	CHECK(mathematics::abs(T(-2.5)) == T(2.5));
	CHECK(mathematics::abs(T(2.5)) == T(2.5));
	CHECK(mathematics::abs(T(0)) == T(0));
}

TEST_CASE_TEMPLATE("abs of integers", T, short, int, long, long long)
{
	CHECK(mathematics::abs(T(-7)) == T(7));
	CHECK(mathematics::abs(T(7)) == T(7));
	CHECK(mathematics::abs(T(0)) == T(0));
}

TEST_CASE("abs of unsigned integers")
{
	CHECK(mathematics::abs(7u) == 7u);
}

TEST_CASE_TEMPLATE("sign", T, float, double, int, long long)
{
	CHECK(sign(T(-3)) == T(-1));
	CHECK(sign(T(3)) == T(1));
	CHECK(sign(T(0)) == T(0));
}

TEST_CASE_TEMPLATE("frac", T, TEST_REAL_TYPES)
{
	CHECK_APPROX(frac(T(2.25)), T(0.25));
	CHECK_APPROX(frac(T(-2.25)), T(0.75)); // x - floor(x)
	CHECK(frac(T(3)) == T(0));
}

TEST_CASE_TEMPLATE("square", T, float, double, int)
{
	CHECK(square(T(3)) == T(9));
	CHECK(square(T(-3)) == T(9));
}

TEST_CASE_TEMPLATE("approxEquals", T, TEST_REAL_TYPES)
{
	CHECK(approxEquals(T(1), T(1) + Constants<T>::TOLERANCE*T(0.5)));
	CHECK_FALSE(approxEquals(T(1), T(1) + Constants<T>::TOLERANCE*T(2)));
	CHECK(approxEquals(T(1), T(1.05), T(0.1)));
	CHECK_FALSE(approxEquals(T(1), T(1.2), T(0.1)));
}

TEST_CASE_TEMPLATE("logarithms", T, TEST_REAL_TYPES)
{
	CHECK_APPROX(mathematics::log(T(8), T(2)), T(3));
	CHECK_APPROX(mathematics::log(T(1000), T(10)), T(3));
	CHECK_APPROX(mathematics::log2(T(8)), T(3));
	CHECK_APPROX(mathematics::log2(T(0.5)), T(-1));
}

TEST_CASE_TEMPLATE("integer log2", T, int, unsigned int, long, unsigned long, long long, unsigned long long, short)
{
	CHECK(mathematics::log2(T(1)) == T(0));
	CHECK(mathematics::log2(T(2)) == T(1));
	CHECK(mathematics::log2(T(3)) == T(1));
	CHECK(mathematics::log2(T(8)) == T(3));
	CHECK(mathematics::log2(T(1000)) == T(9));
}

TEST_CASE_TEMPLATE("radians and degrees", T, TEST_REAL_TYPES)
{
	CHECK_APPROX(radians(T(180)), Constants<T>::PI);
	CHECK_APPROX(degrees(Constants<T>::HALF_PI), T(90));
	CHECK_APPROX(degrees(radians(T(37))), T(37));
}

TEST_CASE_TEMPLATE("lerp", T, TEST_REAL_TYPES)
{
	CHECK(mathematics::lerp(T(2), T(6), T(0)) == T(2));
	CHECK(mathematics::lerp(T(2), T(6), T(1)) == T(6));
	CHECK(mathematics::lerp(T(2), T(6), T(0.25)) == T(3));
	CHECK(mathematics::lerp(T(2), T(6), T(1.5)) == T(8)); // extrapolation
}

TEST_CASE_TEMPLATE("step and pulse", T, float, double, int)
{
	CHECK(step(T(1), T(0)) == T(0));
	CHECK(step(T(1), T(1)) == T(1));
	CHECK(step(T(1), T(2)) == T(1));

	CHECK(pulse(T(1), T(3), T(0)) == T(0));
	CHECK(pulse(T(1), T(3), T(1)) == T(1));
	CHECK(pulse(T(1), T(3), T(2)) == T(1));
	CHECK(pulse(T(1), T(3), T(3)) == T(0));
}

TEST_CASE_TEMPLATE("boxStep, smoothStep and smootherStep", T, TEST_REAL_TYPES)
{
	SUBCASE("clamp outside of the edges")
	{
		CHECK(boxStep(T(1), T(3), T(0)) == T(0));
		CHECK(boxStep(T(1), T(3), T(4)) == T(1));
		CHECK(smoothStep(T(1), T(3), T(0)) == T(0));
		CHECK(smoothStep(T(1), T(3), T(4)) == T(1));
		CHECK(smootherStep(T(1), T(3), T(0)) == T(0));
		CHECK(smootherStep(T(1), T(3), T(4)) == T(1));
	}

	SUBCASE("midpoint is one half")
	{
		CHECK_APPROX(boxStep(T(1), T(3), T(2)), T(0.5));
		CHECK_APPROX(smoothStep(T(1), T(3), T(2)), T(0.5));
		CHECK_APPROX(smootherStep(T(1), T(3), T(2)), T(0.5));
	}

	SUBCASE("polynomial values")
	{
		CHECK_APPROX(boxStep(T(0), T(1), T(0.25)), T(0.25));
		CHECK_APPROX(smoothStep(T(0), T(1), T(0.25)), T(0.15625));	 // 3t^2 - 2t^3
		CHECK_APPROX(smootherStep(T(0), T(1), T(0.25)), T(0.103515625)); // 6t^5 - 15t^4 + 10t^3
	}
}

TEST_CASE("Axis and tag types")
{
	CHECK(int(Axis::X) == 0);
	CHECK(int(Axis::Y) == 1);
	CHECK(int(Axis::Z) == 2);
	CHECK(int(Axis::W) == 3);
	static_assert(Arithmetic<int> && Arithmetic<float> && !Arithmetic<Vector3>);
	static_assert(Normalization<Normalized> && Normalization<Unnormalized> && !Normalization<int>);
}

} // TEST_SUITE
