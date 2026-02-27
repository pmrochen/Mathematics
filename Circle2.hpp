/*
 *	Name: Circle2
 *	Author: Pawel Mrochen
 */

#pragma once

#include <istream>
#include <ostream>
#include <type_traits>
#include <concepts>
#include <utility>
#include <optional>
#include <algorithm>
#include <cstddef>
#include <cmath>
#include "Constants.hpp"
#include "Vector2.hpp"
#include "Interval.hpp"
#include "AxisAlignedRectangle.hpp"

namespace core::mathematics {
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

	// Closest point
	//Vector2<T> getClosestPoint(const Vector2<T>& point) const noexcept; // #TODO
	T getDistanceTo(const Vector2<T>& point) const { return std::max(distance(point, center) - radius, T(0)); }
	T getSignedDistanceTo(const Vector2<T>& point) const noexcept { return distance(point, center) - radius; }

	// Containment and intersection
	bool contains(const Vector2<T>& point) const noexcept;
	bool intersects(const AxisAlignedRectangle<T>& rectangle) const noexcept;
	bool intersects(const Circle2<T>& circle) const noexcept;
	std::optional<Interval<T>> findIntersection(const Line2<T>& line) const noexcept; // #TODO
	std::optional<Interval<T>> findIntersection(const Ray2<T>& ray) const noexcept; // #TODO
	std::optional<Interval<T>> findIntersection(const Segment2<T>& segment) const noexcept; // #TODO

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
	Vector2<T> h(radius, radius);
	return AxisAlignedRectangle<T>(center - h, center + h);
}

template<typename T>
inline bool Circle2<T>::contains(const Vector2<T>& point) const
{
	return (distanceSquared(point, center) <= radius*radius);
}

template<typename T>
inline bool Circle2<T>::intersects(const AxisAlignedRectangle<T>& rectangle) const
{
	T d = T(0);
	
	if (center.x < rectangle.minimum.x)
	{
		T s = center.x - rectangle.minimum.x;
		d += s*s;
	}
	else if (center.x > rectangle.maximum.x)
	{
		T s = center.x - rectangle.maximum.x;
		d += s*s;
	}

	if (center.y < rectangle.minimum.y)
	{
		T s = center.y - rectangle.minimum.y;
		d += s*s;
	}
	else if (center.y > rectangle.maximum.y)
	{
		T s = center.y - rectangle.maximum.y;
		d += s*s;
	}

	return (d <= radius*radius);
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

} // namespace core::mathematics

namespace std {

template<typename T>
struct hash;

template<typename T>
struct hash<::core::mathematics::templates::Circle2<T>>
{
	std::size_t operator()(const ::core::mathematics::templates::Circle2<T>& circle) const noexcept
	{
		std::hash<T> hasher;
		std::size_t seed = hasher(circle.center) + 0x9e3779b9;
		seed ^= hasher(circle.radius) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
		return seed;
	}
};

} // namespace std
