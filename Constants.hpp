/*
 *	Name: Constants
 *	Author: Pawel Mrochen
 */

#pragma once

#include <limits>

namespace core {
namespace mathematics {

//template<typename T> constexpr T /*ZERO_*/TOLERANCE = T();
//template<> constexpr float /*ZERO_*/TOLERANCE = 1e-6f;
//template<> constexpr double /*ZERO_*/TOLERANCE = 1e-15;
//template<typename T> constexpr T INFINITY = std::numeric_limits<T>::infinity();
//template<typename T> constexpr T E = T(2.71828182845904553488L);
//template<typename T> constexpr T LN2 = T(0.693147180559945309417L);
//template<typename T> constexpr T PI = T(3.14159265358979323846L);
//template<typename T> constexpr T PI_SQR = T(9.86960440108935861883L);
//template<typename T> constexpr T HALF_PI = T(1.57079632679489661923L);
//template<typename T> constexpr T TWO_PI = T(6.28318530717958647693L);
//template<typename T> constexpr T INV_PI = T(0.31830988618379069122L);
//template<typename T> constexpr T SQRT_2 = T(1.41421356237309514547L); 
//template<typename T> constexpr T SQRT_3 = T(1.73205080756887719318L);
//template<typename T> constexpr T DEG_TO_RAD = T(0.01745329251994329547L);
//template<typename T> constexpr T RAD_TO_DEG = T(57.29577951308232286465L);

template<typename T>
struct Constants;
//{
//	static constexpr T TOLERANCE/*ZERO_TOLERANCE*/ = T();
//	static constexpr T EPSILON = T();
//	static constexpr T INF = T();
//	static constexpr T MINUS_INF = T();
//	static constexpr T E = T();
//	static constexpr T LN2 = T();
//	static constexpr T PI = T();
//	static constexpr T PI_SQR = T();
//	static constexpr T HALF_PI = T();
//	static constexpr T TWO_PI = T();
//	static constexpr T INV_PI = T();
//	static constexpr T SQRT_2 = T();
//	static constexpr T SQRT_3 = T();
//	static constexpr T DEG_TO_RAD = T();
//	static constexpr T RAD_TO_DEG = T();
//
//	//static constexpr T /*zero*/tolerance() noexcept { return T(); }
//	//static constexpr T infinity() noexcept { return T(); }
//	//static constexpr T e() noexcept { return T(); }
//	//static constexpr T ln2() noexcept { return T(); }
//	//static constexpr T pi() noexcept { return T(); }
//	//static constexpr T piSquare() noexcept { return T(); }
//	//static constexpr T halfPi() noexcept { return T(); }
//	//static constexpr T twoPi() noexcept { return T(); }
//	//static constexpr T inversePi() noexcept { return T(); }
//	//static constexpr T sqrt2() noexcept { return T(); }
//	//static constexpr T sqrt3() noexcept { return T(); }
//	//static constexpr T degToRad() noexcept { return T(); }
//	//static constexpr T radToDeg() noexcept { return T(); }
//};

template<>
struct Constants<float>
{
	static constexpr float TOLERANCE/*ZERO_TOLERANCE*/ = 1e-6f;
	static constexpr float EPSILON = std::numeric_limits<float>::epsilon();
	static constexpr float INF = std::numeric_limits<float>::infinity();
	static constexpr float MINUS_INF = -std::numeric_limits<float>::infinity();
	static constexpr float E = 2.71828182845904553488f;
	static constexpr float LN2 = 0.693147180559945309417f;
	static constexpr float PI = 3.14159265358979323846f;
	static constexpr float PI_SQR = 9.86960440108935861883f;
	static constexpr float HALF_PI = 1.57079632679489661923f;
	static constexpr float TWO_PI = 6.28318530717958647693f;
	static constexpr float INV_PI = 0.31830988618379069122f;
	static constexpr float SQRT_2 = 1.41421356237309514547f;
	static constexpr float SQRT_3 = 1.73205080756887719318f;
	static constexpr float DEG_TO_RAD = 0.01745329251994329547f;
	static constexpr float RAD_TO_DEG = 57.29577951308232286465f;

	//static constexpr float /*zero*/tolerance() noexcept { return 1e-6f; }
	//static constexpr float infinity() noexcept { return std::numeric_limits<float>::infinity(); }
	//static constexpr float e() noexcept { return 2.71828182845904553488f; }
	//static constexpr float ln2() noexcept { return 0.693147180559945309417f; }
	//static constexpr float pi() noexcept { return 3.14159265358979323846f; }
	//static constexpr float piSquare() noexcept { return 9.86960440108935861883f; }
	//static constexpr float halfPi() noexcept { return 1.57079632679489661923f; }
	//static constexpr float twoPi() noexcept { return 6.28318530717958647693f; }
	//static constexpr float inversePi() noexcept { return 0.31830988618379069122f; }
	//static constexpr float sqrt2() noexcept { return 1.41421356237309514547f; }
	//static constexpr float sqrt3() noexcept { return 1.73205080756887719318f; }
	//static constexpr float degToRad() noexcept { return 0.01745329251994329547f; }
	//static constexpr float radToDeg() noexcept { return 57.29577951308232286465f; }
};

template<>
struct Constants<double>
{
	static constexpr double TOLERANCE/*ZERO_TOLERANCE*/ = 1e-15;
	static constexpr double EPSILON = std::numeric_limits<double>::epsilon();
	static constexpr double INF = std::numeric_limits<double>::infinity();
	static constexpr double MINUS_INF = -std::numeric_limits<double>::infinity();
	static constexpr double E = 2.71828182845904553488;
	static constexpr double LN2 = 0.693147180559945309417;
	static constexpr double PI = 3.14159265358979323846;
	static constexpr double PI_SQR = 9.86960440108935861883;
	static constexpr double HALF_PI = 1.57079632679489661923;
	static constexpr double TWO_PI = 6.28318530717958647693;
	static constexpr double INV_PI = 0.31830988618379069122;
	static constexpr double SQRT_2 = 1.41421356237309514547;
	static constexpr double SQRT_3 = 1.73205080756887719318;
	static constexpr double DEG_TO_RAD = 0.01745329251994329547;
	static constexpr double RAD_TO_DEG = 57.29577951308232286465;

	//static constexpr double /*zero*/tolerance() noexcept { return 1e-15; }
	//static constexpr double infinity() noexcept { return std::numeric_limits<double>::infinity(); }
	//static constexpr double e() noexcept { return 2.71828182845904553488; }
	//static constexpr double ln2() noexcept { return 0.693147180559945309417; }
	//static constexpr double pi() noexcept { return 3.14159265358979323846; }
	//static constexpr double piSquare() noexcept { return 9.86960440108935861883; }
	//static constexpr double halfPi() noexcept { return 1.57079632679489661923; }
	//static constexpr double twoPi() noexcept { return 6.28318530717958647693; }
	//static constexpr double inversePi() noexcept { return 0.31830988618379069122; }
	//static constexpr double sqrt2() noexcept { return 1.41421356237309514547; }
	//static constexpr double sqrt3() noexcept { return 1.73205080756887719318; }
	//static constexpr double degToRad() noexcept { return 0.01745329251994329547; }
	//static constexpr double radToDeg() noexcept { return 57.29577951308232286465; }
};

} // namespace mathematics
} // namespace core
