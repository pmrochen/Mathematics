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
#include "Concepts.hpp"
#include "Constants.hpp"
#include "Algebra/Vector2.hpp"
#include "Algebra/Vector3.hpp"
#include "Algebra/Matrix3.hpp"
#include "Algebra/Quaternion.hpp"

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
	static templates::Vector2<T> getVectorOnUnitCircle() requires std::floating_point<T> { return getVectorOnUnitCircle(defaultGenerator); }
	static templates::Vector2<T> getVectorOnUnitCircle(U& generator) requires std::floating_point<T>;
	static templates::Vector2<T> getVectorOnUnitDisk() requires std::floating_point<T> { return getVectorOnUnitDisk(defaultGenerator); }
	static templates::Vector2<T> getVectorOnUnitDisk(U& generator) requires std::floating_point<T>;
	static templates::Vector3<T> getVectorInUnitSphere() requires std::floating_point<T> { return getVectorInUnitSphere(defaultGenerator); }
	static templates::Vector3<T> getVectorInUnitSphere(U& generator) requires std::floating_point<T>;
	static templates::Vector3<T> getVectorOnUnitSphere() requires std::floating_point<T> { return getVectorOnUnitSphere(defaultGenerator); }
	static templates::Vector3<T> getVectorOnUnitSphere(U& generator) requires std::floating_point<T>;
	static templates::Vector3<T> getVectorOnUnitSphericalCone(T halfAngle) requires std::floating_point<T> { return getVectorOnUnitSphericalCone(halfAngle, defaultGenerator); }
	static templates::Vector3<T> getVectorOnUnitSphericalCone(T halfAngle, U& generator) requires std::floating_point<T>;
	static templates::Vector3<T> getVectorOnUnitSphericalCone(const templates::Vector3<T>& direction, T halfAngle) requires std::floating_point<T> { return getVectorOnUnitSphericalCone(direction, halfAngle, defaultGenerator); }
	static templates::Vector3<T> getVectorOnUnitSphericalCone(const templates::Vector3<T>& direction, T halfAngle, U& generator) requires std::floating_point<T>;
	static auto getRotationMatrix() requires std::floating_point<T> { return getRotationMatrix(defaultGenerator); }
	static auto getRotationMatrix(U& generator) requires std::floating_point<T>;
	static auto getQuaternion() requires std::floating_point<T> { return getQuaternion(defaultGenerator); }
	static auto getQuaternion(U& generator) requires std::floating_point<T>;

	static U defaultGenerator;
};

template<typename T, typename U> requires (std::floating_point<T> || std::integral<T>) U Random<T, U>::defaultGenerator{};

template<typename T, typename U>
	requires (std::floating_point<T> || std::integral<T>)
/*static*/ inline T Random<T, U>::get(T minimum, T maximum, U& generator)
{
	double double01 = double(generator() - U::min())/(double(U::max()) - double(U::min()));
	if constexpr (std::is_integral_v<T>)
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
	requires (std::floating_point<T> || std::integral<T>)
/*static*/ inline T Random<T, U>::getRightOpen(T minimum, T maximum, U& generator) requires std::floating_point<T>
{
	double double01RightOpen = double(generator() - U::min())/(double(U::max()) - double(U::min()) + 1.0);
	T result = T(minimum + (maximum - minimum)*double01RightOpen);
	return (result < maximum) ? result : std::nextafter(maximum, minimum); // guard against rounding up when narrowing to T
}

template<typename T, typename U>
	requires (std::floating_point<T> || std::integral<T>)
/*static*/ inline T Random<T, U>::get01(U& generator) requires std::floating_point<T>
{
	return T(double(generator() - U::min())/(double(U::max()) - double(U::min())));
}

template<typename T, typename U>
	requires (std::floating_point<T> || std::integral<T>)
/*static*/ inline T Random<T, U>::get01RightOpen(U& generator) requires std::floating_point<T>
{
	T result = T(double(generator() - U::min())/(double(U::max()) - double(U::min()) + 1.0));
	return (result < T(1)) ? result : std::nextafter(T(1), T(0)); // guard against rounding up when narrowing to T
}

template<typename T, typename U>
	requires (std::floating_point<T> || std::integral<T>)
/*static*/ inline templates::Vector2<T> Random<T, U>::getVectorOnUnitCircle(U& generator) requires std::floating_point<T>
{
	T theta = getRightOpen(T(0), Constants<T>::TWO_PI, generator);
	return templates::Vector2<T>(std::cos(theta), std::sin(theta));
}

template<typename T, typename U>
	requires (std::floating_point<T> || std::integral<T>)
/*static*/ inline templates::Vector2<T> Random<T, U>::getVectorOnUnitDisk(U& generator) requires std::floating_point<T>
{
	templates::Vector2<T> v{ Uninitialized() };
	do
	{
		v.set(get(T(-1), T(1), generator), get(T(-1), T(1), generator));
	} while (v.getMagnitudeSquared() > T(1));
	return v;
}

template<typename T, typename U>
	requires (std::floating_point<T> || std::integral<T>)
/*static*/ inline templates::Vector3<T> Random<T, U>::getVectorInUnitSphere(U& generator) requires std::floating_point<T>
{
	templates::Vector3<T> v{ Uninitialized() };
	do
	{
		v.set(get(T(-1), T(1), generator), get(T(-1), T(1), generator), get(T(-1), T(1), generator));
	} while (v.getMagnitudeSquared() > T(1));
	return v;
}

template<typename T, typename U>
	requires (std::floating_point<T> || std::integral<T>)
/*static*/ inline templates::Vector3<T> Random<T, U>::getVectorOnUnitSphere(U& generator) requires std::floating_point<T>
{
	T z = get(T(-1), T(1), generator);
	T phi = getRightOpen(T(0), Constants<T>::TWO_PI, generator);
	T s = std::sqrt(T(1) - z*z);
	return templates::Vector3<T>(s*std::cos(phi), s*std::sin(phi), z);
}

template<typename T, typename U>
	requires (std::floating_point<T> || std::integral<T>)
/*static*/ inline templates::Vector3<T> Random<T, U>::getVectorOnUnitSphericalCone(T halfAngle, U& generator) requires std::floating_point<T>
{
	T z = get(std::cos(halfAngle), T(1), generator);
	T phi = getRightOpen(T(0), Constants<T>::TWO_PI, generator);
	T s = std::sqrt(T(1) - z*z);
	return templates::Vector3<T>(s*std::cos(phi), s*std::sin(phi), z);
}

template<typename T, typename U>
	requires (std::floating_point<T> || std::integral<T>)
/*static*/ inline templates::Vector3<T> Random<T, U>::getVectorOnUnitSphericalCone(const templates::Vector3<T>& direction, T halfAngle, U& generator) requires std::floating_point<T>
{
	templates::Vector3<T> tangentX = normalize(cross((std::fabs(direction.z) < T(0.999)) ? templates::Vector3<T>::UNIT_Z : templates::Vector3<T>::UNIT_X, direction));
	templates::Vector3<T> tangentY = cross(direction, tangentX);
	T theta = std::acos(get(std::cos(halfAngle), T(1), generator));
	T phi = getRightOpen(T(0), Constants<T>::TWO_PI, generator);
	return std::sin(theta)*(std::cos(phi)*tangentX + std::sin(phi)*tangentY) + std::cos(theta)*direction;
}

template<typename T, typename U>
	requires (std::floating_point<T> || std::integral<T>)
/*static*/ inline auto Random<T, U>::getRotationMatrix(U& generator) requires std::floating_point<T>
{
	//T d = std::clamp(range, T(0), T(1));
	T theta = get(T(0), Constants<T>::TWO_PI/* *d*/, generator);
	T phi = get(T(0), Constants<T>::TWO_PI, generator);
	T m = get(T(0), T(2)/* *d*/, generator);
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
	return templates::Matrix3<T>(vx*sx - ct, vx*sy - st, vx*vz, vy*sx + st, vy*sy - ct, vy*vz, vz*sx, vz*sy, T(1) - m);
}

template<typename T, typename U>
	requires (std::floating_point<T> || std::integral<T>)
/*static*/ inline auto Random<T, U>::getQuaternion(U& generator) requires std::floating_point<T>
{
	T t1 = get(T(0), Constants<T>::TWO_PI, generator);
	T t2 = get(T(0), Constants<T>::TWO_PI, generator);
	T m = get01(generator);
	T r1 = std::sqrt(T(1) - m);
	T r2 = std::sqrt(m);
	T s1 = std::sin(t1);
	T c1 = std::cos(t1);
	T s2 = std::sin(t2);
	T c2 = std::cos(t2);
	return templates::Quaternion<T>(s1*r1, c1*r1, s2*r2, c2*r2);
}

template<typename T>
inline T random(T minimum, T maximum)
{
	return Random<T>::get(minimum, maximum);
}

} // namespace mathematics
