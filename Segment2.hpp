/*
 *	Name: Segment2
 *	Author: Pawel Mrochen
 */

#pragma once

#include <istream>
#include <ostream>
#include <type_traits>
#include <concepts>
#include <utility>
#include <tuple>
#include <optional>
#include <algorithm>
#include <cstddef>
#include <cmath>
#include "Constants.hpp"
#include "Vector2.hpp"
#include "Interval.hpp"
#include "Line2.hpp"
#include "Ray2.hpp"

namespace core::mathematics {
namespace templates {
	
template<typename T>
	requires std::floating_point<T>
struct AxisAlignedRectangle;

template<typename T>
	requires std::floating_point<T>
struct Circle2;

template<typename T>
	requires std::floating_point<T>
struct Segment2
{
	using Real = T;
	using ConstArg = const Segment2&;
	using ConstResult = const Segment2&;

	Segment2() = default;
	explicit Segment2(Uninitialized) noexcept : start(Uninitialized()), end(Uninitialized()) {}
	Segment2(const Vector2<T>& start, const Vector2<T>& end) noexcept : start(start), end(end) {}
	explicit Segment2(const std::pair<Vector3<T>, Vector3<T>>& t) noexcept : start(t.first), end(t.second) {}
	explicit Segment2(const std::tuple<Vector3<T>, Vector3<T>>& t) noexcept : start(std::get<0>(t)), end(std::get<1>(t)) {}
	explicit Segment2(const Line2<T>& line) noexcept : start(line.origin), end(line.origin + line.direction) {}
	Segment2(const Line2<T>& line, const Interval<T>& interval) noexcept : start(line.evaluate(interval.minimum)), end(line.evaluate(interval.maximum)) {}
	explicit Segment2(const Ray2<T>& ray) noexcept : start(ray.origin), end(ray.origin + ray.direction) {}
	Segment2(const Ray2<T>& ray, const Interval<T>& interval) noexcept : start(ray.evaluate(interval.minimum)), end(ray.evaluate(interval.maximum)) {}

	//explicit operator std::pair<Vector3<T>, Vector3<T>>() { return { start, end }; }
	//explicit operator std::tuple<Vector3<T>, Vector3<T>>() { return { start, end }; }
	//Vector2 operator()(T t) const noexcept { return lerp(start, end, t); }
	bool operator==(const Segment2& segment) const noexcept { return (start == segment.start) && (end == segment.end); }
	bool operator!=(const Segment2& segment) const noexcept { return !(*this == segment); }

	template<typename A> void serialize(A& ar) { ar(start, end); }

	// Properties
	bool approxEquals(const Segment2& segment) const noexcept;
	bool approxEquals(const Segment2& segment, T tolerance) const noexcept;
	bool isFinite() const noexcept { return start.isFinite() && end.isFinite(); }
	Segment2& set(const Vector2<T>& start, const Vector2<T>& end) noexcept { this->start = start; this->end = end; return *this; }
	const Vector2<T>& getStart() const noexcept { return start; }
	void setStart(const Vector2<T>& start) noexcept { this->start = start; }
	const Vector2<T>& getEnd() const noexcept { return end; }
	void setEnd(const Vector2<T>& end) noexcept { this->end = end; }
	Vector2<T> getDirection() const noexcept { return normalize(end - start); }
	T getLength() const noexcept { return distance(start, end); }
	Vector2<T> getCenter() const noexcept { return lerp(start, end, T(0.5)); }
	T getSlope() const noexcept { return (end.y - start.y)/(end.x - start.x); }
	T getInclinationAngle() const noexcept { return std::atan2(end.y - start.y, end.x - start.x); }

	// Evaluation
	Vector2<T> evaluate(T t) const noexcept { return lerp(start, end, t); }

	// Closest points
	Vector2<T> getClosestPoint(const Vector2<T>& point) const;
	T getDistanceTo(const Vector2<T>& point) const { return distance(getClosestPoint(point), point); }

	// Intersection
	bool intersects(const Line2<T>& line) const noexcept { return findIntersection(line).has_value(); }
	//bool intersects(const Ray2<T>& ray) const noexcept { return findIntersection(ray).has_value(); }
	bool intersects(const Segment2& segment) const noexcept { return findIntersection(segment).has_value(); }
	bool intersects(const AxisAlignedRectangle& rectangle) const noexcept { return findIntersection(rectangle).has_value(); }
	bool intersects(const Circle2<T>& circle) const noexcept { return findIntersection(circle).has_value(); }
	std::optional<T> findIntersection(const Line2<T>& line) const noexcept;
	//std::optional<T> findIntersection(const Ray2<T>& ray) const noexcept;
	std::optional<T> findIntersection(const Segment2& segment) const;
	std::optional<Interval<T>> findIntersection(const AxisAlignedRectangle& rectangle) const noexcept;
	std::optional<Interval<T>> findIntersection(const Circle2<T>& circle) const noexcept;
	//template<ScalarOrVector2<T> U> std::optional<U> findIntersection(const Line2<T>& line) const noexcept;
	//template<ScalarOrVector2<T> U> std::optional<U> findIntersection(const Ray2<T>& ray) const noexcept;
	//template<ScalarOrVector2<T> U> std::optional<U> findIntersection(const Segment2& segment) const;
	//template<IntervalOrSegment2<T> U> std::optional<U> findIntersection(const AxisAlignedRectangle& rectangle) const noexcept;
	//template<IntervalOrSegment2<T> U> std::optional<U> findIntersection(const Circle2<T>& circle) const;

	Vector2<T> start;
	Vector2<T> end;
};

template<typename C, typename T, typename U>
	requires std::floating_point<U>
inline std::basic_istream<C, T>& operator>>(std::basic_istream<C, T>& s, Segment2<U>& segment)
{ 
	return s >> segment.start >> std::ws >> segment.end;
}

template<typename C, typename T, typename U>
	requires std::floating_point<U>
inline std::basic_ostream<C, T>& operator<<(std::basic_ostream<C, T>& s, const Segment2<U>& segment)
{ 
	constexpr C WS(0x20);
	return s << segment.start << WS << segment.end;
}

template<typename T>
inline bool Segment2<T>::approxEquals(const Segment2<T>& segment) const
{
	return start.approxEquals(segment.start) && end.approxEquals(segment.end);
}

template<typename T>
inline bool Segment2<T>::approxEquals(const Segment2<T>& segment, T tolerance) const
{
	return start.approxEquals(segment.start, tolerance) && end.approxEquals(segment.end, tolerance);
}

template<typename T>
inline Vector2<T> Segment2<T>::getClosestPoint(const Vector2<T>& point) const
{
	Vector2<T> direction = end - start;
	return std::clamp(dot(point - start, direction)/dot(direction, direction), T(0), T(1))*direction + start;
}

template<typename T>
inline std::optional<T> Segment2<T>::findIntersection(const Line2<T>& line) const
{
	Vector2<T> direction = end - start;
	T d1CrossD2 = cross(direction, line.direction);
	if (std::fabs(d1CrossD2) < Constants<T>::TOLERANCE)
	{
		if (!(std::fabs(cross(normalize(line.origin - start), direction)) < Constants<T>::TOLERANCE))
			return {};
		return std::optional<T>(T(0)); // collinear
	}

	T t = cross(line.origin - start, line.direction)/d1CrossD2;
	//T u = cross(line.origin - start, direction)/d1CrossD2;
	if ((t >= T(0)) && (t <= T(1))
		return std::optional<T>(t);
	return {};
}

template<typename T>
inline std::optional<T> Segment2<T>::findIntersection(const Segment2& segment) const
{
	Vector2<T> direction = end - start;
	T d1CrossD2 = cross(direction, segment.end - segment.start);
	if (std::fabs(d1CrossD2) < Constants<T>::TOLERANCE)
	{
		if (!(std::fabs(cross(normalize(segment.start - start), direction)) < Constants<T>::TOLERANCE))
			return {};

		T d1DotD1 = dot(direction, direction);
		T d1DotD2 = dot(direction, segment.end - segment.start);
		T t0 = dot(segment.start - start, direction)/d1DotD1;
		T t1 = t0 + d1DotD2/d1DotD1;
		if (d1DotD2 < T(0))
			std::swap(t0, t1);
		if ((t0 <= T(1)) && (t1 >= T(0)))
			return std::optional<T>(std::max(t0, T(0)));
		return {};
	}

	T t = cross(segment.start - start, segment.end - segment.start)/d1CrossD2;
	T u = cross(segment.start - start, direction)/d1CrossD2;
	if ((t >= T(0)) && (t <= T(1)) && (u >= T(0)) && (u <= T(1)))
		return std::optional<T>(t);
	return {};
}

template<typename T>
template<ScalarOrVector2<T> U>
inline std::optional<U> Segment2<T>::findIntersection(const Line2<T>& line) const
{
	std::optional<T> result = findIntersection(line);
	if constexpr (std::is_same_v<U, T>)
		return result;
	else //if constexpr (std::is_same_v<U, Vector3<T>>)
		return result.has_value() ? std::optional<U>(line.evaluate(result.value())) : std::optional<U>();
}

template<typename T>
template<ScalarOrVector2<T> U> 
inline std::optional<U> Segment2<T>::findIntersection(const Segment2& segment) const
{
	std::optional<T> result = findIntersection(segment);
	if constexpr (std::is_same_v<U, T>)
		return result;
	else //if constexpr (std::is_same_v<U, Vector3<T>>)
		return result.has_value() ? std::optional<U>(segment.evaluate(result.value())) : std::optional<U>();
}

} // namespace templates

#if MATHEMATICS_DOUBLE
using Segment2 = templates::Segment2<double>;
using Segment2Arg = templates::Segment2<double>::ConstArg;
using Segment2Result = templates::Segment2<double>::ConstResult;
#else
using Segment2 = templates::Segment2<float>;
using Segment2Arg = templates::Segment2<float>::ConstArg;
using Segment2Result = templates::Segment2<float>::ConstResult;
#endif

} // namespace core::mathematics

namespace std {

template<typename T>
struct hash;

template<typename T>
struct hash<::core::mathematics::templates::Segment2<T>>
{
	std::size_t operator()(const ::core::mathematics::templates::Segment2<T>& segment) const noexcept
	{
		std::hash<T> hasher;
		std::size_t seed = hasher(segment.start) + 0x9e3779b9;
		seed ^= hasher(segment.end) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
		return seed;
	}
};

} // namespace std

#include "AxisAlignedRectangle.hpp"
#include "Circle2.hpp"

namespace core::mathematics::templates {

template<typename T>
inline std::optional<Interval<T>> Segment2<T>::findIntersection(const AxisAlignedRectangle<T>& rectangle) const
{
	std::optional<Interval<T>> result = Line2(segment.start, segment.end - segment.start).findIntersection(rectangle);
	if (result.has_value() && (result.value().maximum >= T(0)) && (result.value().minimum <= T(1)))
	{
		const Interval<T>& interval = result.value();
		if (interval.minimum != interval.maximum)
			return { std::in_place, std::max(interval.minimum, T(0)), std::min(interval.maximum, T(1)) };
		else
			return result;
	}
	else
	{
		return {};
	}
}

template<typename T>
inline std::optional<Interval<T>> Segment2<T>::findIntersection(const Circle2<T>& circle) const
{
	std::optional<Interval<T>> result = Line2(segment.start, segment.end - segment.start).findIntersection(circle);
	if (result.has_value() && (result.value().maximum >= T(0)) && (result.value().minimum <= T(1)))
	{
		const Interval<T>& interval = result.value();
		if (interval.minimum != interval.maximum)
			return { std::in_place, std::max(interval.minimum, T(0)), std::min(interval.maximum, T(1)) };
		else
			return result;
	}
	else
	{
		return {};
	}
}

} // namespace core::mathematics::templates
