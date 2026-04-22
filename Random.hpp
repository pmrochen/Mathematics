/*
 *	Name: Random
 *	Author: Pawel Mrochen
 */

#pragma once

#include <type_traits>
#include <concepts>
#include <random>
#include <algorithm>
#include <cmath>
#include "Constants.h"
#include "Vector2.hpp"
#include "Vector3.hpp"
#include "Matrix3.hpp"
#include "Quaternion.hpp"

namespace mathematics {

template<typename T, typename U = std::mt19937>
	requires (std::floating_point<T> || std::integral<T>)
struct Random
{
	using GeneratorType = U;

	static U& getDefaultGenerator() noexcept { return defaultGenerator; }

	static T get(T minimum, T maximum) { return get(minimum, maximum, defaultGenerator); }
	static T get(T minimum, T maximum, U& generator);
	static T getRightOpen(T minimum, T maximum) requires std::floating_point<T> { return getRightOpen(minimum, maximum, defaultGenerator); }
	static T getRightOpen(T minimum, T maximum, U& generator) requires std::floating_point<T>;
	static T get01() requires std::floating_point<T> { return get01(defaultGenerator); }
	static T get01(U& generator) requires std::floating_point<T>;
	static T get01RightOpen() requires std::floating_point<T> { return get01RightOpen(defaultGenerator); }
	static T get01RightOpen(U& generator) requires std::floating_point<T>;
	static Vector2<T> getVectorOnUnitCircle() requires std::floating_point<T> { return getVectorOnUnitCircle(defaultGenerator); }
	static Vector2<T> getVectorOnUnitCircle(U& generator) requires std::floating_point<T>;
	static Vector2<T> getVectorOnUnitDisk() requires std::floating_point<T> { return getVectorOnUnitDisk(defaultGenerator); }
	static Vector2<T> getVectorOnUnitDisk(U& generator) requires std::floating_point<T>;
	static Vector3<T> getVectorInUnitSphere() requires std::floating_point<T> { return getVectorInUnitSphere(defaultGenerator); }
	static Vector3<T> getVectorInUnitSphere(U& generator) requires std::floating_point<T>;
	static Vector3<T> getVectorOnUnitSphere() requires std::floating_point<T> { return getVectorOnUnitSphere(defaultGenerator); }
	static Vector3<T> getVectorOnUnitSphere(U& generator) requires std::floating_point<T>;
	static Vector3<T> getVectorOnUnitSphericalCone(T halfAngle) requires std::floating_point<T> { return getVectorOnUnitSphericalCone(halfAngle, defaultGenerator); }
	static Vector3<T> getVectorOnUnitSphericalCone(T halfAngle, U& generator) requires std::floating_point<T>;
	static Vector3<T> getVectorOnUnitSphericalCone(const Vector3<T>& direction, T halfAngle) requires std::floating_point<T> { return getVectorOnUnitSphericalCone(direction, halfAngle, defaultGenerator); }
	static Vector3<T> getVectorOnUnitSphericalCone(const Vector3<T>& direction, T halfAngle, U& generator) requires std::floating_point<T>;
	static Matrix3<T> getRotationMatrix() requires std::floating_point<T> { return getRotationMatrix(defaultGenerator); }
	static Matrix3<T> getRotationMatrix(U& generator) requires std::floating_point<T>;
	static Quaternion<T> getQuaternion() requires std::floating_point<T> { return getQuaternion(defaultGenerator); }
	static Quaternion<T> getQuaternion(U& generator) requires std::floating_point<T>;

	static U defaultGenerator;
};

template<typename T, typename U> U Random<T, U>::defaultGenerator{};

template<typename T, typename U>
/*static*/ inline T Random<T, U>::get(T minimum, T maximum, U& generator)
{
	double double01 = (generator() - U::min())/(double)(U::max() - U::min() - 1);
	if constexpr (std::is_intergral_v<T>)
	{
		if constexpr (sizeof(T) > sizeof(long))
			return std::clamp(T(std::llround(minimum - 0.5 + ((maximum - minimum) + 1.0)*double01)), minimum, maximum);
		else
			return std::clamp(T(std::lround(minimum - 0.5 + ((maximum - minimum) + 1.0)*double01)), minimum, maximum);
	}
	else // floating point
		return T(minimum + (maximum - minimum)*double01);
}

template<typename T, typename U>
/*static*/ inline T Random<T, U>::getRightOpen(T minimum, T maximum, U& generator)
{
	double double01RightOpen = (generator() - U::min())/(double)(U::max() - U::min());
	return T(minimum + (maximum - minimum)*double01RightOpen);
}

template<typename T, typename U>
/*static*/ inline T Random<T, U>::get01(U& generator)
{
	return T((generator() - U::min())/(double)(U::max() - U::min() - 1));
}

template<typename T, typename U>
/*static*/ inline T Random<T, U>::get01RightOpen(U& generator)
{
	return T((generator() - U::min())/(double)(U::max() - U::min()));
}

template<typename T, typename U>
/*static*/ inline Vector2<T> Random<T, U>::getVectorOnUnitCircle(U& generator)
{
	T theta = getRightOpen<T>(T(0), Constants<T>::TWO_PI, generator);
	return Vector2<T>(std::cos(theta), std::sin(theta));
}

template<typename T, typename U>
/*static*/ inline Vector2<T> Random<T, U>::getVectorOnUnitDisk(U& generator)
{
	Vector2<T> v(Unintialized());
	do
	{
		v.set(get<T>(T(-1), T(1), generator), get<T>(T(-1), T(1), generator));
	} while (v.getMagnitudeSquared() > T(1));
	return v;
}

template<typename T, typename U>
/*static*/ inline Vector3<T> Random<T, U>::getVectorInUnitSphere(U& generator)
{
	Vector3<T> v(Unintialized());
	do
	{
		v.set(get<T>(T(-1), T(1), generator), get<T>(T(-1), T(1), generator), get<T>(T(-1), T(1), generator));
	} while (v.getMagnitudeSquared() > T(1));
	return v;
}

template<typename T, typename U>
/*static*/ inline Vector3<T> Random<T, U>::getVectorOnUnitSphere(U& generator)
{
	T z = get<T>(T(-1), T(1), generator);
	T phi = getRightOpen<T>(T(0), Constants<T>::TWO_PI, generator);
	T s = std::sqrt(T(1) - z*z);
	return Vector3<T>(s*std::cos(phi), s*std::sin(phi), z);
}

template<typename T, typename U>
/*static*/ inline Vector3<T> Random<T, U>::getVectorOnUnitSphericalCone(T halfAngle, U& generator)
{
	T z = get<T>(cos(halfAngle), T(1), generator);
	T phi = getRightOpen<T>(T(0), Constants<T>::TWO_PI, generator);
	T s = std::sqrt(T(1) - z*z);
	return Vector3(s*std::cos(phi), s*std::sin(phi), z);
}

template<typename T, typename U>
/*static*/ inline Vector3<T> Random<T, U>::getVectorOnUnitSphericalCone(const Vector3<T>& direction, T halfAngle, U& generator)
{
	Vector3<T> tangentX = normalize(cross((std::fabs(direction.z) < T(0.999)) ? Vector3<T>::UNIT_Z : Vector3<T>::UNIT_X, direction));
	Vector3<T> tangentY = cross(direction, tangentX);
	T theta = std::acos(get<T>(cos(halfAngle), T(1), generator));
	T phi = getRightOpen<T>(T(0), Constants<T>::TWO_PI, generator);
	return std::sin(theta)*(std::cos(phi)*tangentX + std::sin(phi)*tangentY) + std::cos(theta)*direction;
}

template<typename T, typename U>
/*static*/ inline Matrix3<T> Random<T, U>::getRotationMatrix(U& generator)
{
	//T d = std::clamp(range, T(0), T(1));
	T theta = get<T>(T(0), Constants<T>::TWO_PI/* *d*/, generator);
	T phi = get<T>(T(0), Constants<T>::TWO_PI, generator);
	T m = get<T>(T(0), T(2)/* *d*/, generator);
	T r = std::sqrt(m);
	T sp = std::sin(phi);
	T cp = std::cos(phi);
	T vx = sp*r;
	T vy = cp*r;
	T vz = std::sqrt(T(2) - m);
	T st = std::sin(theta);
	T ct = std::cos(theta);
	T sx = vx*ct - vy*st;
	T sy = vx*st + vy*ct;
	return Matrix3<T>(vx*sx - ct, vx*sy - st, vx*vz, vy*sx + st, vy*sy - ct, vy*vz, vz*sx, vz*sy, T(1) - m);
}

template<typename T, typename U>
/*static*/ inline Quaternion<T> Random<T, U>::getQuaternion(U& generator)
{
	T t1 = get<T>(T(0), Constants<T>::TWO_PI, generator);
	T t2 = get<T>(T(0), Constants<T>::TWO_PI, generator);
	T m = get01<T>(generator);
	T r1 = std::sqrt(T(1) - m);
	T r2 = std::sqrt(m);
	T s1 = std::sin(t1);
	T c1 = std::cos(t1);
	T s2 = std::sin(t2);
	T c2 = std::cos(t2);
	return Quaternion<T>(s1*r1, c1*r1, s2*r2, c2*r2);
}

template<typename T>
inline T random(T minimum, T maximum)
{
	return Random<T>::get(minimum, maximum);
}

} // namespace mathematics
