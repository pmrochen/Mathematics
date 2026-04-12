/*
 *	Name: Containment
 *	Author: Pawel Mrochen
 */

#pragma once

#include "Vector3.hpp"

namespace mathematics::containment {

template<typename T>
	requires std::floating_point<T>
inline bool testCylinderPoint(const Vector3<T>& center, const Vector3<T>& axis, T height, T radius, const Vector3<T>& point) noexcept
{
	Vector3<T> pd = point - (center - (height*T(0.5))*axis);
	T d = dot(pd, height*axis);
	T lengthSq = height*height;
	if ((d < T(0)) || (d > lengthSq))
		return false;

	return ((dot(pd, pd) - d*d/lengthSq) <= radius*radius);
}

template<typename T>
	requires std::floating_point<T>
inline bool testConePoint(const Vector3<T>& vertex, const Vector3<T>& axis, T height, T radius, const Vector3<T>& point) noexcept
{
	Vector3<T> diff = point - vertex;
	T d = dot(diff, axis);
	if ((d < T(0)) || (d > height))
		return false;

	T r = d*radius/height;
	return (distanceSquared(diff, d*axis) <= r*r);
}

} // namespace mathematics::containment
