/*
 *	Name: Circle2
 *	Author: Pawel Mrochen
 */

#pragma once

#include <istream>
#include <ostream>
#include <type_traits>
#include <concepts>
#include <algorithm>
#include <functional>
#include <utility>
#include <cstddef>
#include <cmath>
#include "Constants.hpp"
#include "Vector2.hpp"
#include "AxisAlignedRectangle.hpp"

namespace mathematics {
namespace templates {

template<typename T>
	requires std::floating_point<T>
struct Circle2
{
	using Real = T;
	using ConstArg = const Circle2&;
	using ConstResult = const Circle2&;

	Circle2() noexcept : center(), radius() {}
	explicit Circle2(Uninitialized) noexcept : center(Uninitialized()) {}
	Circle2(const Vector2<T>& center, T radius) noexcept : center(center), radius(radius) {}
	//explicit Circle2(const Ellipse2<T>& ellipse) noexcept;

	bool operator==(const Circle2& circle) const noexcept { return (center == circle.center) && (radius == circle.radius); }
	bool operator!=(const Circle2& circle) const noexcept { return !(*this == circle); }

	template<typename A> void serialize(A& ar) { ar(center, radius); }

	// Properties
	bool approxEquals(const Circle2& circle) const noexcept;
	bool approxEquals(const Circle2& circle, T tolerance) const noexcept;
	bool isFinite() const noexcept { return center.isFinite() && radius.isFinite(); }
	Circle2& set(const Vector2<T>& center, T radius) noexcept { this->center = center; this->radius = radius; return *this; }
	const Vector2<T>& getCenter() const noexcept { return center; }
	void setCenter(const Vector2<T>& center) noexcept { this->center = center; }
	T getRadius() const noexcept { return radius; }
	void setRadius(T radius) noexcept { this->radius = radius; }
	T getDiameter() const noexcept { return radius*T(2); }
	void setDiameter(T diameter) noexcept { radius = diameter*T(0.5); }
	T getCircumference() const noexcept { return T(2)*Constants<T>::PI*radius; }
	T getArea() const noexcept { return Constants<T>::PI*radius*radius; }

	// Circumscribed rectangle
	AxisAlignedRectangle<T> getCircumscribedRectangle() const noexcept;

	// Transformation
	Circle2& translate(const Vector2<T>& offset) { center += offset; return *this; }

	// Closest point
	//Vector2<T> getClosestPoint(const Vector2<T>& point) const noexcept; // #TODO
	T getDistanceTo(const Vector2<T>& point) const { return std::max(distance(point, center) - radius, T(0)); }
	T getSignedDistanceTo(const Vector2<T>& point) const noexcept { return distance(point, center) - radius; }

	// Containment and intersection
	bool contains(const Vector2<T>& point) const noexcept { return (distanceSquared(point, center) <= radius*radius); }
	bool intersects(const AxisAlignedRectangle<T>& rectangle) const noexcept;
	bool intersects(const Circle2& circle) const noexcept;

	Vector2<T> center;
	T radius;
};

template<typename C, typename T, typename U>
	requires std::floating_point<U>
inline std::basic_istream<C, T>& operator>>(std::basic_istream<C, T>& s, Circle2<U>& circle)
{ 
	return s >> circle.center >> std::ws >> circle.radius;
}

template<typename C, typename T, typename U>
	requires std::floating_point<U>
inline std::basic_ostream<C, T>& operator<<(std::basic_ostream<C, T>& s, const Circle2<U>& circle)
{ 
	constexpr C WS(0x20);
	return s << circle.center << WS << circle.radius;
}

template<typename T>
inline bool Circle2<T>::approxEquals(const Circle2<T>& circle) const
{
	return center.approxEquals(circle.center) && 
		(std::fabs(circle.radius - radius) < Constants<T>::TOLERANCE);
}

template<typename T>
inline bool Circle2<T>::approxEquals(const Circle2<T>& circle, T tolerance) const
{
	return center.approxEquals(circle.center, tolerance) && 
		(std::fabs(circle.radius - radius) < tolerance);
}

template<typename T>
inline AxisAlignedRectangle<T> Circle2<T>::getCircumscribedRectangle() const
{
	Vector2<T> halfDims(radius);
	return AxisAlignedRectangle<T>(center - halfDims, center + halfDims);
}

template<typename T>
inline bool Circle2<T>::intersects(const Circle2<T>& circle) const
{
	T d = circle.radius + radius;
	return (distanceSquared(circle.center, center) <= d*d);
}

} // namespace templates

#if MATHEMATICS_DOUBLE
using Circle2 = templates::Circle2<double>;
using Circle2Arg = templates::Circle2<double>::ConstArg;
using Circle2Result = templates::Circle2<double>::ConstResult;
#else
using Circle2 = templates::Circle2<float>;
using Circle2Arg = templates::Circle2<float>::ConstArg;
using Circle2Result = templates::Circle2<float>::ConstResult;
#endif

} // namespace mathematics

namespace std {

template<typename T>
struct hash;

template<typename T>
struct hash<::mathematics::templates::Circle2<T>>
{
	size_t operator()(const ::mathematics::templates::Circle2<T>& circle) const noexcept
	{
		size_t seed = hash<typename ::mathematics::templates::Vector2<T>>()(circle.center) + 0x9e3779b9;
		seed ^= hash<T>()(circle.radius) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
		return seed;
	}
};

} // namespace std

#include "Intersections.inl"

namespace mathematics::templates {

template<typename T>
inline bool Circle2<T>::intersects(const AxisAlignedRectangle<T>& rectangle) const
{
	return intersections::testAxisAlignedRectangleCircle(rectangle.minimum, rectangle.maximum, center, radius);
}

} // namespace mathematics::templates
