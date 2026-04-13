/*
 *	Name: LineSegment2
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
#include <tuple>
#include <optional>
#include <iterator>
#include <cstddef>
#include <cmath>
#include "Constants.hpp"
#include "Vector2.hpp"
#include "Interval.hpp"
#include "Line2.hpp"
#include "Ray2.hpp"

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
struct LineSegment2
{
	using Real = T;
	using ConstArg = const LineSegment2&;
	using ConstResult = const LineSegment2&;

	LineSegment2() = default;
	explicit LineSegment2(Uninitialized) noexcept : start(Uninitialized()), end(Uninitialized()) {}
	LineSegment2(const Vector2<T>& start, const Vector2<T>& end) noexcept : start(start), end(end) {}
	explicit LineSegment2(const std::pair<Vector3<T>, Vector3<T>>& t) noexcept : start(t.first), end(t.second) {}
	explicit LineSegment2(const std::tuple<Vector3<T>, Vector3<T>>& t) noexcept : start(std::get<0>(t)), end(std::get<1>(t)) {}
	explicit LineSegment2(const Line2<T>& line) noexcept : start(line.origin), end(line.origin + line.direction) {}
	LineSegment2(const Line2<T>& line, const Interval<T>& interval) noexcept : start(line.evaluate(interval.minimum)), end(line.evaluate(interval.maximum)) {}
	explicit LineSegment2(const Ray2<T>& ray) noexcept : start(ray.origin), end(ray.origin + ray.direction) {}
	LineSegment2(const Ray2<T>& ray, const Interval<T>& interval) noexcept : start(ray.evaluate(interval.minimum)), end(ray.evaluate(interval.maximum)) {}

	//explicit operator std::pair<Vector3<T>, Vector3<T>>() { return { start, end }; }
	//explicit operator std::tuple<Vector3<T>, Vector3<T>>() { return { start, end }; }
	//Vector2 operator()(T t) const noexcept { return lerp(start, end, t); }
	bool operator==(const LineSegment2& segment) const noexcept { return (start == segment.start) && (end == segment.end); }
	bool operator!=(const LineSegment2& segment) const noexcept { return !(*this == segment); }

	template<typename A> void serialize(A& ar) { ar(start, end); }

	// Properties
	bool approxEquals(const LineSegment2& segment) const noexcept;
	bool approxEquals(const LineSegment2& segment, T tolerance) const noexcept;
	bool isFinite() const noexcept { return start.isFinite() && end.isFinite(); }
	LineSegment2& set(const Vector2<T>& start, const Vector2<T>& end) noexcept { this->start = start; this->end = end; return *this; }
	const Vector2<T>& getStart() const noexcept { return start; }
	void setStart(const Vector2<T>& start) noexcept { this->start = start; }
	const Vector2<T>& getEnd() const noexcept { return end; }
	void setEnd(const Vector2<T>& end) noexcept { this->end = end; }
	Vector2<T> getDirection() const noexcept { return normalize(end - start); }
	T getLength() const noexcept { return distance(start, end); }
	Vector2<T> getCenter() const noexcept { return lerp(start, end, T(0.5)); }
	T getSlope() const noexcept { return (end.y - start.y)/(end.x - start.x); }
	T getInclinationAngle() const noexcept { return std::atan2(end.y - start.y, end.x - start.x); }

	// Endpoints
	template<std::output_iterator<Vector2<T>> O> O copyEndpoints(O target) const;
	std::pair<Vector2<T>, Vector2<T>> getEndpoints() const noexcept { return { start, end }; }

	// Transformation
	LineSegment2& translate(const Vector2<T>& offset) noexcept { start += offset; end += offset; return *this; }

	// Evaluation
	Vector2<T> evaluate(T t) const noexcept { return lerp(start, end, t); }

	// Closest points
	Vector2<T> getClosestPoint(const Vector2<T>& point) const;
	T getDistanceTo(const Vector2<T>& point) const { return distance(getClosestPoint(point), point); }

	// Intersection
	bool intersects(const Line2<T>& line) const noexcept { return findIntersection(line).has_value(); }
	//bool intersects(const Ray2<T>& ray) const noexcept { return findIntersection(ray).has_value(); } // #TODO
	bool intersects(const LineSegment2& segment) const noexcept { return findIntersection(segment).has_value(); }
	bool intersects(const AxisAlignedRectangle<T>& rectangle) const noexcept { return findIntersection(rectangle).has_value(); }
	bool intersects(const Circle2<T>& circle) const noexcept { return findIntersection(circle).has_value(); }
	std::optional<T> findIntersection(const Line2<T>& line) const noexcept;
	//std::optional<T> findIntersection(const Ray2<T>& ray) const noexcept; // #TODO
	std::optional<T> findIntersection(const LineSegment2& segment) const;
	std::optional<Interval<T>> findIntersection(const AxisAlignedRectangle<T>& rectangle) const noexcept;
	std::optional<Interval<T>> findIntersection(const Circle2<T>& circle) const noexcept;
	//template<ScalarOrVector2<T> U> std::optional<U> findIntersection(const Line2<T>& line) const noexcept;
	//template<ScalarOrVector2<T> U> std::optional<U> findIntersection(const Ray2<T>& ray) const noexcept;
	//template<ScalarOrVector2<T> U> std::optional<U> findIntersection(const LineSegment2& segment) const;
	//template<IntervalOrLineSegment2<T> U> std::optional<U> findIntersection(const AxisAlignedRectangle& rectangle) const noexcept;
	//template<IntervalOrLineSegment2<T> U> std::optional<U> findIntersection(const Circle2<T>& circle) const;

	Vector2<T> start;
	Vector2<T> end;
};

template<typename C, typename T, typename U>
	requires std::floating_point<U>
inline std::basic_istream<C, T>& operator>>(std::basic_istream<C, T>& s, LineSegment2<U>& segment)
{ 
	return s >> segment.start >> std::ws >> segment.end;
}

template<typename C, typename T, typename U>
	requires std::floating_point<U>
inline std::basic_ostream<C, T>& operator<<(std::basic_ostream<C, T>& s, const LineSegment2<U>& segment)
{ 
	constexpr C WS(0x20);
	return s << segment.start << WS << segment.end;
}

template<typename T>
inline bool LineSegment2<T>::approxEquals(const LineSegment2<T>& segment) const
{
	return start.approxEquals(segment.start) && end.approxEquals(segment.end);
}

template<typename T>
inline bool LineSegment2<T>::approxEquals(const LineSegment2<T>& segment, T tolerance) const
{
	return start.approxEquals(segment.start, tolerance) && end.approxEquals(segment.end, tolerance);
}

template<typename T>
template<std::output_iterator<Vector2<T>> O>
inline O LineSegment2<T>::copyEndpoints(O target) const
{
	*target++ = start;
	*target++ = end;
	return target;
}

template<typename T>
inline Vector2<T> LineSegment2<T>::getClosestPoint(const Vector2<T>& point) const
{
	Vector2<T> direction = end - start;
	return std::clamp(dot(point - start, direction)/dot(direction, direction), T(0), T(1))*direction + start;
}

//template<typename T>
//template<ScalarOrVector2<T> U>
//inline std::optional<U> LineSegment2<T>::findIntersection(const Line2<T>& line) const
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
//inline std::optional<U> LineSegment2<T>::findIntersection(const LineSegment2& segment) const
//{
//	std::optional<T> result = findIntersection(segment);
//	if constexpr (std::is_same_v<U, T>)
//		return result;
//	else //if constexpr (std::is_same_v<U, Vector3<T>>)
//		return result.has_value() ? std::optional<U>(segment.evaluate(result.value())) : std::optional<U>();
//}

} // namespace templates

#if MATHEMATICS_DOUBLE
using LineSegment2 = templates::LineSegment2<double>;
using LineSegment2Arg = templates::LineSegment2<double>::ConstArg;
using LineSegment2Result = templates::LineSegment2<double>::ConstResult;
#else
using LineSegment2 = templates::LineSegment2<float>;
using LineSegment2Arg = templates::LineSegment2<float>::ConstArg;
using LineSegment2Result = templates::LineSegment2<float>::ConstResult;
#endif

} // namespace mathematics

namespace std {

template<typename T>
struct hash;

template<typename T>
struct hash<::mathematics::templates::LineSegment2<T>>
{
	size_t operator()(const ::mathematics::templates::LineSegment2<T>& segment) const noexcept
	{
		hash<typename ::mathematics::templates::Vector2<T>> hasher;
		size_t seed = hasher(segment.start) + 0x9e3779b9;
		seed ^= hasher(segment.end) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
		return seed;
	}
};

} // namespace std

#include "AxisAlignedRectangle.hpp"
#include "Circle2.hpp"
#include "Intersections.inl"

namespace mathematics::templates {

template<typename T>
inline std::optional<T> LineSegment2<T>::findIntersection(const Line2<T>& line) const
{
	return intersections::findLineLineSegment<std::optional<T>>(line.origin, line.direction, start, end);

}

template<typename T>
inline std::optional<T> LineSegment2<T>::findIntersection(const LineSegment2& segment) const
{
	return intersections::findLineSegmentLineSegment<std::optional<T>>(start, end, segment.start, segment.end);
}

template<typename T>
inline std::optional<Interval<T>> LineSegment2<T>::findIntersection(const AxisAlignedRectangle<T>& rectangle) const
{
	std::optional<Interval<T>> result = intersections::findLineAxisAlignedRectangle<std::optional<Interval<T>>>(start, end - start, 
		rectangle.minimum, rectangle.maximum);
	
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
inline std::optional<Interval<T>> LineSegment2<T>::findIntersection(const Circle2<T>& circle) const
{
	std::optional<Interval<T>> result = intersections::findLineNSphere<std::optional<Interval<T>>>(start, end - start, circle.center, circle.radius);
	
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

} // namespace mathematics::templates
