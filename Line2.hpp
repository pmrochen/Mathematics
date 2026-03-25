/*
 *	Name: Line2
 *	Author: Pawel Mrochen
 */

#pragma once

#include <istream>
#include <ostream>
#include <type_traits>
#include <concepts>
#include <utility>
#include <optional>
#include <iterator>
#include <algorithm>
#include <cstddef>
#include <cmath>
#include "Constants.hpp"
#include "Vector2.hpp"
#include "Interval.hpp"

namespace mathematics {
namespace templates {

//template<typename T, typename U>
//concept ScalarOrVector2 = (std::same_as<T, U> || std::same_as<T, Vector2<U>>);

//template<typename T, typename U>
//concept IntervalOrSegment2 = (std::same_as<T, Interval<U>> || std::same_as<T, Segment2<U>>);
	
template<typename T>
	requires std::floating_point<T>
struct AxisAlignedRectangle;

template<typename T>
	requires std::floating_point<T>
struct Circle2;

template<typename T>
	requires std::floating_point<T>
struct Ray2;

template<typename T>
	requires std::floating_point<T>
struct Line2
{
	using Real = T;
	using ConstArg = const Line2&;
	using ConstResult = const Line2&;

	Line2() = default;
	explicit Line2(Uninitialized) noexcept : origin(Uninitialized()), direction(Uninitialized()) {}
	Line2(const Vector2<T>& origin, const Vector2<T>& direction) noexcept : origin(origin), direction(direction) {}
	Line2(const Vector2<T>& origin, T inclinationAngle) noexcept;
	Line2(const Ray2<T>& ray) noexcept;
	//explicit Line2(const Segment2<T>& segment) noexcept;

	//Vector2<T> operator()(T t) const noexcept { return (origin + t*direction); }
	bool operator==(const Line2& line) const noexcept { return (origin == line.origin) && (direction == line.direction); }
	bool operator!=(const Line2& line) const noexcept { return !(*this == line); }

	template<typename A> void serialize(A& ar) { ar(origin, direction); }

	const Ray2<T>& asRay() const noexcept;

	// Least-squares fit of a line
	//template<std::input_iterator I, std::sentinel_for<I> S> static Line2 computeBestFit(I first, S last); // #TODO

	// Properties
	bool approxEquals(const Line2& line) const noexcept;
	bool approxEquals(const Line2& line, T tolerance) const noexcept;
	//bool isFinite() const noexcept { return origin.isFinite() && direction.isFinite(); }
	Line2& set(const Vector2<T>& origin, const Vector2<T>& direction) noexcept { this->origin = origin; this->direction = direction; return *this; }
	const Vector2<T>& getOrigin() const noexcept { return origin; }
	void setOrigin(const Vector2<T>& origin) noexcept { this->origin = origin; }
	const Vector2<T>& getDirection() const noexcept { return direction; }
	void setDirection(const Vector2<T>& direction) noexcept { this->direction = direction; }
	T getSlope() const noexcept { return direction.y/direction.x; }
	T getInclinationAngle() const noexcept { return std::atan2(direction.y, direction.x); }
	bool isParallel(const Line2& line) const noexcept;
	bool isCoincident(const Line2& line) const noexcept;

	// Transformation
	Line2& translate(const Vector2<T>& offset) noexcept { origin += offset; return *this; }
	Line2& normalize() noexcept { direction.normalize(); return *this; }

	// Evaluation
	Vector2<T> evaluate(T t) const noexcept { return (origin + t*direction); }

	// Closest points
	Vector2<T> getClosestPoint(const Vector2<T>& point) const noexcept;											// normalized line
	template<Normalization U> Vector2<T> getClosestPoint(const Vector2<T>& point) const noexcept;
	T getDistanceTo(const Vector2<T>& point) const noexcept { return distance(getClosestPoint(point), point); }	// normalized line
	template<Normalization U> T getDistanceTo(const Vector2<T>& point) const noexcept { return distance(getClosestPoint<U>(point), point); }
	T getSignedDistanceTo(const Vector2<T>& point) const noexcept;												// normalized line
	//template<Normalization U> T getSignedDistanceTo(const Vector2<T>& point) const noexcept; // #TODO
	T getDistanceTo(const Line2& line) const noexcept { return std::fabs(getSignedDistance(line)); }			// normalized line
	//template<Normalization U> T getDistanceTo(const Line2& line) const noexcept; // #TODO
	T getSignedDistanceTo(const Line2& line) const noexcept;													// normalized line
	//template<Normalization U> T getSignedDistanceTo(const Line2& line) const noexcept; // #TODO

	// Intersection
	bool intersects(const Line2& line) const noexcept;
	bool intersects(const Ray2<T>& ray) const noexcept { return findIntersection(ray).has_value(); }
	bool intersects(const Segment2<T>& segment) const noexcept { return findIntersection(segment).has_value(); }
	bool intersects(const AxisAlignedRectangle<T>& rectangle) const noexcept { return findIntersection(rectangle).has_value(); }
	bool intersects(const Circle2<T>& circle) const noexcept;
	template<Normalization U> bool intersects(const Circle2<T>& circle) const noexcept;
	std::optional<T> findIntersection(const Line2& line) const noexcept;
	std::optional<T> findIntersection(const Ray2<T>& ray) const noexcept;
	std::optional<T> findIntersection(const Segment2<T>& segment) const noexcept;
	std::optional<Interval<T>> findIntersection(const AxisAlignedRectangle<T>& rectangle) const noexcept;
	std::optional<Interval<T>> findIntersection(const Circle2<T>& circle) const noexcept;
	template<Normalization U> std::optional<Interval<T>> findIntersection(const Circle2<T>& circle) const noexcept;
	//template<ScalarOrVector2<T> U> std::optional<U> findIntersection(const Line2& line) const noexcept;
	//template<ScalarOrVector2<T> U> std::optional<U> findIntersection(const Segment2<T>& segment) const noexcept;
	//template<IntervalOrSegment2<T> U> std::optional<U> findIntersection(const AxisAlignedRectangle<T>& rectangle) const noexcept;
	//template<IntervalOrSegment2<T> U> std::optional<U> findIntersection(const Circle2<T>& circle) const noexcept;

	Vector2<T> origin;
	Vector2<T> direction;
};

template<typename T>
inline Line2<T>::Line2(const Vector2<T>& origin, T inclinationAngle) : 
	origin(origin), 
	direction(std::cos(inclinationAngle), std::sin(inclinationAngle))
{
}

template<typename C, typename T, typename U>
	requires std::floating_point<U>
inline std::basic_istream<C, T>& operator>>(std::basic_istream<C, T>& s, Line2<U>& line)
{ 
	return s >> line.origin >> std::ws >> line.direction;
}

template<typename C, typename T, typename U>
	requires std::floating_point<U>
inline std::basic_ostream<C, T>& operator<<(std::basic_ostream<C, T>& s, const Line2<U>& line)
{ 
	constexpr C WS(0x20);
	return s << line.origin << WS << line.direction;
}

template<typename T>
inline bool Line2<T>::approxEquals(const Line2<T>& line) const
{ 
	return origin.approxEquals(line.origin) && direction.approxEquals(line.direction); 
}

template<typename T>
inline bool Line2<T>::approxEquals(const Line2<T>& line, T tolerance) const
{ 
	return origin.approxEquals(line.origin, tolerance) && direction.approxEquals(line.direction, tolerance); 
}

template<typename T>
inline bool Line2<T>::isParallel(const Line2<T>& line) const
{
	return (std::fabs(cross(direction, line.direction)) < Constants<T>::TOLERANCE);
}

template<typename T>
inline bool Line2<T>::isCoincident(const Line2<T>& line) const
{
	return (std::fabs(cross(direction, line.direction)) < Constants<T>::TOLERANCE) &&
		(std::fabs(cross(normalize(line.origin - origin), direction)) < Constants<T>::TOLERANCE);
}

template<typename T>
inline Vector2<T> Line2<T>::getClosestPoint(const Vector2<T>& point) const
{
	return dot(point - origin, direction)*direction + origin;
}

template<typename T>
template<Normalization U>
inline Vector2<T> Line2<T>::getClosestPoint(const Vector2<T>& point) const
{
	if costexpr(std::is_same_v<U, Normalized>)
		return dot(point - origin, direction)*direction + origin;
	else
		return (dot(point - origin, direction)/dot(direction, direction))*direction + origin;
}

template<typename T>
inline T Line2<T>::getSignedDistanceTo(const Vector2<T>& point) const
{
	return cross(origin - point, direction);
}

template<typename T>
inline T Line2<T>::getSignedDistanceTo(const Line2& line) const
{
	return (std::fabs(cross(direction, line.direction)) < Constants<T>::TOLERANCE) ?
		cross((origin + dot(line.origin - origin, direction)*direction) - line.origin, direction) :
		T();
}

template<typename T>
inline bool Line2<T>::intersects(const Line2& line) const
{
	return !(std::fabs(cross(direction, line.direction)) < Constants<T>::TOLERANCE) ||
		(std::fabs(cross(normalize(line.origin - origin), direction)) < Constants<T>::TOLERANCE);
}

//template<typename T>
//template<ScalarOrVector2<T> U>
//inline std::optional<U> Line2<T>::findIntersection(const Line2& line) const
//{
//	std::optional<T> result = findIntersection(line);
//	if constexpr (std::is_same_v<U, T>)
//		return result;
//	else //if constexpr (std::is_same_v<U, Vector3<T>>)
//		return result.has_value() ? std::optional<U>(line.evaluate(result.value())) : std::optional<U>();
//}
//
//template<typename T>
//template<ScalarOrVector2<T> U>
//inline std::optional<U> Line2<T>::findIntersection(const Segment2<T>& segment) const
//{
//	std::optional<T> result = findIntersection(segment);
//	if constexpr (std::is_same_v<U, T>)
//		return result;
//	else //if constexpr (std::is_same_v<U, Vector3<T>>)
//		return result.has_value() ? std::optional<U>(segment.evaluate(result.value())) : std::optional<U>();
//}

template<typename T>
	requires std::floating_point<T>
inline Line2<T> normalize(const Line2<T>& line) noexcept
{
	Line2<T> l(line);
	l.normalize();
	return l;
}

template<typename T>
	requires std::floating_point<T>
inline Line2<T> normalize(Line2<T>&& line) noexcept
{
	line.normalize();
	return line;
}

} // namespace templates

#if MATHEMATICS_DOUBLE
using Line2 = templates::Line2<double>;
using Line2Arg = templates::Line2<double>::ConstArg;
using Line2Result = templates::Line2<double>::ConstResult;
#else
using Line2 = templates::Line2<float>;
using Line2Arg = templates::Line2<float>::ConstArg;
using Line2Result = templates::Line2<float>::ConstResult;
#endif

} // namespace mathematics

namespace std {

template<typename T>
struct hash;

template<typename T>
struct hash<::mathematics::templates::Line2<T>>
{
	std::size_t operator()(const ::mathematics::templates::Line2<T>& line) const noexcept
	{
		std::hash<T> hasher;
		std::size_t seed = hasher(line.origin) + 0x9e3779b9;
		seed ^= hasher(line.direction) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
		return seed;
	}
};

} // namespace std

#include "Ray2.hpp"
#include "AxisAlignedRectangle.hpp"
#include "Circle2.hpp"
#include "Intersections.inl"

namespace mathematics::templates {

template<typename T>
inline Line2<T>::Line2(const Ray2<T>& ray) : origin(ray.origin), direction(ray.direction)
{
}

template<typename T>
inline const Ray2<T>& Line2<T>::asRay() const
{ 
	return reinterpret_cast<const Ray2<T>&>(*this);
}

template<typename T>
inline bool Line2<T>::intersects(const Circle2<T>& circle) const
{
	return intersections::testLineNSphere(origin, direction, circle.center, circle.radius);
}

template<typename T>
template<Normalization U>
inline bool Line2<T>::intersects(const Circle2<T>& circle) const
{
	if costexpr(std::is_same_v<U, Normalized>)
		return intersections::testNormalizedLineNSphere(origin, direction, circle.center, circle.radius);
	else
		return intersections::testLineNSphere(origin, direction, circle.center, circle.radius);
}

template<typename T>
inline std::optional<T> Line2<T>::findIntersection(const Line2& line) const
{
	return intersections::findLineLine<std::optional<T>>(origin, direction, line.origin, line.direction);
}

template<typename T>
inline std::optional<T> Line2<T>::findIntersection(const Ray2& ray) const
{
	return intersections::findLineRay<std::optional<T>>(origin, direction, ray.origin, ray.direction);
}

template<typename T>
inline std::optional<T> Line2<T>::findIntersection(const Segment2<T>& segment) const
{
	return intersections::findLineSegment<std::optional<T>>(origin, direction, segment.start, segment.end);
}

template<typename T>
inline std::optional<Interval<T>> Line2<T>::findIntersection(const AxisAlignedRectangle& rectangle) const
{
	return intersections::findLineAxisAlignedRectangle<std::optional<Interval<T>>>(origin, direction, rectangle.minimum, rectangle.maximum);
}

template<typename T>
inline std::optional<Interval<T>> Line2<T>::findIntersection(const Circle2<T>& circle) const
{
	return intersections::findLineNSphere<std::optional<Interval<T>>>(origin, direction, circle.center, circle.radius);
}

template<typename T>
template<Normalization U>
inline std::optional<Interval<T>> Line2<T>::findIntersection(const Circle2<T>& circle) const
{
	if costexpr(std::is_same_v<U, Normalized>)
		return intersections::findNormalizedLineNSphere<std::optional<Interval<T>>>(origin, direction, circle.center, circle.radius);
	else
		return intersections::findLineNSphere<std::optional<Interval<T>>>(origin, direction, circle.center, circle.radius);
}

} // namespace mathematics::templates
