/*
 *	Name: Ray2
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
#include "Line2.hpp"

namespace mathematics {
namespace templates {
	
template<typename T>
	requires std::floating_point<T>
struct AxisAlignedRectangle;

template<typename T>
	requires std::floating_point<T>
struct Circle2;

template<typename T>
	requires std::floating_point<T>
struct Ray2
{
	using Real = T;
	using ConstArg = const Ray2&;
	using ConstResult = const Ray2&;

	Ray2() = default;
	explicit Ray2(Uninitialized) noexcept : origin(Uninitialized()), direction(Uninitialized()) {}
	Ray2(const Vector2<T>& origin, const Vector2<T>& direction) noexcept : origin(origin), direction(direction) {}
	explicit Ray2(const Line2<T>& line) noexcept : origin(line.origin), direction(line.direction) {}
	//explicit Ray2(const Segment2<T>& segment) noexcept;

	//Vector3<T> operator()(T t) const noexcept { return (origin + t*direction); }
	bool operator==(const Ray2& ray) const noexcept { return (origin == ray.origin) && (direction == ray.direction); }
	bool operator!=(const Ray2& ray) const noexcept { return !(*this == ray); }

	template<typename A> void serialize(A& ar) { ar(origin, direction); }

	const Line2<T>& asLine() const noexcept { return reinterpret_cast<const Line2<T>&>(*this); }

	// Properties
	bool approxEquals(const Ray2& ray) const noexcept;
	bool approxEquals(const Ray2& ray, T tolerance) const noexcept;
	//bool isFinite() const noexcept { return origin.isFinite() && direction.isFinite(); }
	Ray2& set(const Vector2<T>& origin, const Vector2<T>& direction) noexcept { this->origin = origin; this->direction = direction; return *this; }
	const Vector2<T>& getOrigin() const noexcept { return origin; }
	void setOrigin(const Vector2<T>& origin) noexcept { this->origin = origin; }
	const Vector2<T>& getDirection() const noexcept { return direction; }
	void setDirection(const Vector2<T>& direction) noexcept { this->direction = direction; }
	T getSlope() const noexcept { return direction.y/direction.x; }
	T getInclinationAngle() const noexcept { return std::atan2(direction.y, direction.x); }

	// Transformation
	Ray2& translate(const Vector2<T>& offset) noexcept { origin += offset; return *this; }
	Ray2& normalize() noexcept { direction.normalize(); return *this; }

	// Evaluation
	Vector2<T> evaluate(T t) const noexcept { return (origin + t*direction); }

	// Closest points
	Vector2<T> getClosestPoint(const Vector2<T>& point) const;											// normalized ray
	template<Normalization U> Vector2<T> getClosestPoint(const Vector2<T>& point) const;
	T getDistanceTo(const Vector2<T>& point) const { return distance(getClosestPoint(point), point); }	// normalized ray
	template<Normalization U> T getDistanceTo(const Vector2<T>& point) const { return distance(getClosestPoint<U>(point), point); }

	// Intersection
	bool intersects(const Line2<T>& line) const noexcept { return findIntersection(line).has_value(); }
	//bool intersects(const Ray2& ray) const noexcept { return findIntersection(ray).has_value(); } // #TODO
	//bool intersects(const Segment2<T>& segment) const noexcept { return findIntersection(segment).has_value(); } // #TODO
	bool intersects(const AxisAlignedRectangle& rectangle) const noexcept { return findIntersection(rectangle).has_value(); }
	bool intersects(const Circle2<T>& circle) const noexcept { return findIntersection(circle).has_value(); }
	template<Normalization U> bool intersects(const Circle2<T>& circle) const noexcept { return findIntersection<U>(circle).has_value(); }
	std::optional<T> findIntersection(const Line2<T>& line) const noexcept;
	//std::optional<T> findIntersection(const Ray2& ray) const noexcept; // #TODO
	//std::optional<T> findIntersection(const Segment2<T>& segment) const noexcept; // #TODO
	std::optional<Interval<T>> findIntersection(const AxisAlignedRectangle& rectangle) const noexcept;
	std::optional<Interval<T>> findIntersection(const Circle2<T>& circle) const noexcept;
	template<Normalization U> std::optional<Interval<T>> findIntersection(const Circle2<T>& circle) const noexcept;
	//template<ScalarOrVector2<T> U> std::optional<U> findIntersection(const Line2<T>& line) const noexcept;
	//template<ScalarOrVector2<T> U> std::optional<U> findIntersection(const Ray2& ray) const noexcept;
	//template<IntervalOrSegment2<T> U> std::optional<U> findIntersection(const AxisAlignedRectangle& rectangle) const noexcept;
	//template<IntervalOrSegment2<T> U> std::optional<U> findIntersection(const Circle2<T>& circle) const noexcept;

	Vector2<T> origin;
	Vector2<T> direction;
};

template<typename C, typename T, typename U>
	requires std::floating_point<U>
inline std::basic_istream<C, T>& operator>>(std::basic_istream<C, T>& s, Ray2<U>& ray)
{ 
	return s >> ray.origin >> std::ws >> ray.direction;
}

template<typename C, typename T, typename U>
	requires std::floating_point<U>
inline std::basic_ostream<C, T>& operator<<(std::basic_ostream<C, T>& s, const Ray2<U>& ray)
{ 
	constexpr C WS(0x20);
	return s << ray.origin << WS << ray.direction;
}

template<typename T>
inline bool Ray2<T>::approxEquals(const Ray2<T>& ray) const
{
	return origin.approxEquals(ray.origin) && direction.approxEquals(ray.direction);
}

template<typename T>
inline bool Ray2<T>::approxEquals(const Ray2<T>& ray, T tolerance) const
{
	return origin.approxEquals(ray.origin, tolerance) && direction.approxEquals(ray.direction, tolerance);
}

template<typename T>
inline Vector2<T> Ray2<T>::getClosestPoint(const Vector2<T>& point) const
{
	return std::max(dot(point - origin, direction), T(0))*direction + origin;
}

template<typename T>
template<Normalization U>
inline Vector2<T> Ray2<T>::getClosestPoint(const Vector2<T>& point) const
{
	if costexpr(std::is_same_v<U, Normalized>)
		return std::max(dot(point - origin, direction), T(0))*direction + origin;
	else
		return std::max(dot(point - origin, direction)/dot(direction, direction), T(0))*direction + origin;
}

//template<typename T>
//template<ScalarOrVector2<T> U>
//inline std::optional<U> Ray2<T>::findIntersection(const Line2<T>& line) const
//{
//	std::optional<T> result = findIntersection(line);
//	if constexpr (std::is_same_v<U, T>)
//		return result;
//	else //if constexpr (std::is_same_v<U, Vector3<T>>)
//		return result.has_value() ? std::optional<U>(line.evaluate(result.value())) : std::optional<U>();
//}

template<typename T>
	requires std::floating_point<T>
inline Ray2<T> normalize(const Ray2<T>& ray) noexcept
{
	Ray2<T> r(ray);
	r.normalize();
	return r;
}

template<typename T>
	requires std::floating_point<T>
inline Ray2<T> normalize(Ray2<T>&& ray) noexcept
{
	ray.normalize();
	return ray;
}

} // namespace templates

#if MATHEMATICS_DOUBLE
using Ray2 = templates::Ray2<double>;
using Ray2Arg = templates::Ray2<double>::ConstArg;
using Ray2Result = templates::Ray2<double>::ConstResult;
#else
using Ray2 = templates::Ray2<float>;
using Ray2Arg = templates::Ray2<float>::ConstArg;
using Ray2Result = templates::Ray2<float>::ConstResult;
#endif

} // namespace mathematics

namespace std {

template<typename T>
struct hash;

template<typename T>
struct hash<::mathematics::templates::Ray2<T>>
{
	std::size_t operator()(const ::mathematics::templates::Ray2<T>& ray) const noexcept
	{
		std::hash<T> hasher;
		std::size_t seed = hasher(ray.origin) + 0x9e3779b9;
		seed ^= hasher(ray.direction) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
		return seed;
	}
};

} // namespace std

#include "AxisAlignedRectangle.hpp"
#include "Circle2.hpp"
#include "Intersections.inl"

namespace mathematics::templates {

template<typename T>
inline std::optional<T> Ray2<T>::findIntersection(const Line2<T>& line) const
{
	return intersections::findLineRay<std::optional<T>>(line.origin, line.direction, origin, direction);
}

template<typename T>
inline std::optional<Interval<T>> Ray2<T>::findIntersection(const AxisAlignedRectangle<T>& rectangle) const
{
	std::optional<Interval<T>> result = intersections::findLineAxisAlignedRectangle<std::optional<Interval<T>>>(origin, direction, 
		rectangle.minimum, rectangle.maximum);

	if (result.has_value() && (result.value().maximum >= T(0)))
	{
		const Interval<T>& interval = result.value();
		if (interval.minimum != interval.maximum)
			return { std::in_place, std::max(interval.minimum, T(0)), interval.maximum };
		else
			return result;
	}
	else
	{
		return {};
	}
}

template<typename T>
inline std::optional<Interval<T>> Ray2<T>::findIntersection(const Circle2<T>& circle) const
{
	std::optional<Interval<T>> result = intersections::findLineNSphere<std::optional<Interval<T>>>(origin, direction, circle.center, circle.radius); 
	
	if (result.has_value() && (result.value().maximum >= T(0)))
	{
		const Interval<T>& interval = result.value();
		if (interval.minimum != interval.maximum)
			return { std::in_place, std::max(interval.minimum, T(0)), interval.maximum };
		else
			return result;
	}
	else
	{
		return {};
	}
}

template<typename T>
template<Normalization U>
inline std::optional<Interval<T>> Ray2<T>::findIntersection(const Circle2<T>& circle) const
{
	std::optional<Interval<T>> result;
	if costexpr(std::is_same_v<U, Normalized>)
		result = intersections::findNormalizedLineNSphere<std::optional<Interval<T>>>(origin, direction, circle.center, circle.radius);
	else
		result = intersections::findLineNSphere<std::optional<Interval<T>>>(origin, direction, circle.center, circle.radius);

	if (result.has_value() && (result.value().maximum >= T(0)))
	{
		const Interval<T>& interval = result.value();
		if (interval.minimum != interval.maximum)
			return { std::in_place, std::max(interval.minimum, T(0)), interval.maximum };
		else
			return result;
	}
	else
	{
		return {};
	}
}

} // namespace mathematics::templates
