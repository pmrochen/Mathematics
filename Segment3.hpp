/*
 *	Name: Segment3
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
#include <iterator>
#include <algorithm>
#include <cstddef>
#include <cmath>
#include "Vector3.hpp"
#include "Matrix3.hpp"
#include "AffineTransform.hpp"
#include "Interval.hpp"
#include "Line3.hpp"
#include "Ray3.hpp"

namespace mathematics {
namespace templates {
	
template<typename T>
	requires std::floating_point<T>
struct HalfSpace;

template<typename T>
	requires std::floating_point<T>
struct Plane;

template<typename T>
	requires std::floating_point<T>
struct Triangle3;

template<typename T>
	requires std::floating_point<T>
struct AxisAlignedBox;

template<typename T>
	requires std::floating_point<T>
struct OrientedBox;

template<typename T>
	requires std::floating_point<T>
struct Sphere;

template<typename T>
	requires std::floating_point<T>
struct Ellipsoid;

#if SIMD_HAS_FLOAT4

template<>
struct HalfSpace<float>;

template<>
struct Plane<float>;

#endif /* SIMD_HAS_FLOAT4 */

template<typename T>
	requires std::floating_point<T>
struct Segment3
{
	using Real = T;
	using ConstArg = const Segment3&;
	using ConstResult = const Segment3&;

	Segment3() = default;
	explicit Segment3(Uninitialized) noexcept : start(Uninitialized()), end(Uninitialized()) {}
	Segment3(const Vector3<T>& start, const Vector3<T>& end) noexcept : start(start), end(end) {}
	explicit Segment3(const std::pair<Vector3<T>, Vector3<T>>& t) noexcept : start(t.first), end(t.second) {}
	explicit Segment3(const std::tuple<Vector3<T>, Vector3<T>>& t) noexcept : start(std::get<0>(t)), end(std::get<1>(t)) {}
	explicit Segment3(const Line3<T>& line) noexcept : start(line.origin), end(line.origin + line.direction) {}
	Segment3(const Line3<T>& line, const Interval<T>& interval) noexcept : start(line.evaluate(interval.minimum)), end(line.evaluate(interval.maximum)) {}
	explicit Segment3(const Ray3<T>& ray) noexcept : start(ray.origin), end(ray.origin + ray.direction) {}
	Segment3(const Ray3<T>& ray, const Interval<T>& interval) noexcept : start(ray.evaluate(interval.minimum)), end(ray.evaluate(interval.maximum)) {}

	//explicit operator std::pair<Vector3<T>, Vector3<T>>() { return { start, end }; }
	//explicit operator std::tuple<Vector3<T>, Vector3<T>>() { return { start, end }; }
	//Vector3 operator()(T t) const noexcept { return lerp(start, end, t); }
	bool operator==(const Segment3& segment) const noexcept { return (start == segment.start) && (end == segment.end); }
	bool operator!=(const Segment3& segment) const noexcept { return !(*this == segment); }

	template<typename A> void serialize(A& ar) { ar(start, end); }

	// Properties
	bool approxEquals(const Segment3& segment) const noexcept;
	bool approxEquals(const Segment3& segment, T tolerance) const noexcept;
	bool isFinite() const noexcept { return start.isFinite() && end.isFinite(); }
	Segment3& set(const Vector3<T>& start, const Vector3<T>& end) noexcept { this->start = start; this->end = end; return *this; }
	const Vector3<T>& getStart() const noexcept { return start; }
	void setStart(const Vector3<T>& start) noexcept { this->start = start; }
	const Vector3<T>& getEnd() const noexcept { return end; }
	void setEnd(const Vector3<T>& end) noexcept { this->end = end; }
	Vector3<T> getDirection() const noexcept { return normalize(end - start); }
	T getLength() const noexcept { return distance(start, end); }
	Vector3<T> getCenter() const noexcept { return lerp(start, end, T(0.5)); }

	// Endpoints
	template<std::output_iterator<Vector3<T>> O> O copyEndpoints(O target) const;
	std::pair<Vector3<T>, Vector3<T>> getEndpoints() const noexcept { return { start, end }; }

	// Transformation
	Segment3& translate(const Vector3<T>& offset) noexcept { start += offset; end += offset; return *this; }
	Segment3& transform(const Matrix3<T>& matrix) noexcept;
	Segment3& transform(const AffineTransform<T>& transformation) noexcept;

	// Evaluation
	Vector3<T> evaluate(T t) const noexcept { return lerp(start, end, t); }

	// Closest points
	Vector3<T> getClosestPoint(const Vector3<T>& point) const;
	T getDistanceTo(const Vector3<T>& point) const { return distance(getClosestPoint(point), point); }

	// Intersection
	bool intersects(const HalfSpace<T>& halfSpace) const noexcept;
	bool intersects(const Plane<T>& plane) const noexcept { return findIntersection(plane).has_value(); }
	bool intersects(const Triangle3<T>& triangle) const noexcept { return findIntersection(triangle).has_value(); }
	bool intersects(const AxisAlignedBox& box) const noexcept { return findIntersection(box).has_value(); }
	bool intersects(const OrientedBox& box) const noexcept { return findIntersection(box).has_value(); }
	bool intersects(const Sphere<T>& sphere) const noexcept { return findIntersection(sphere).has_value(); }
	bool intersects(const Ellipsoid<T>& ellipsoid) const noexcept;
	std::optional<T> findIntersection(const Plane<T>& plane) const noexcept;
	std::optional<T> findIntersection(const Triangle3<T>& triangle) const noexcept;
	//template<ScalarOrVector3<T> U> std::optional<U> findIntersection(const Plane<T>& plane) const noexcept;
	std::optional<Interval<T>> findIntersection(const AxisAlignedBox<T>& box) const noexcept;
	std::optional<Interval<T>> findIntersection(const OrientedBox<T>& box) const noexcept;
	std::optional<Interval<T>> findIntersection(const Sphere<T>& sphere) const noexcept;
	std::optional<Interval<T>> findIntersection(const Ellipsoid<T>& ellipsoid) const noexcept;

	Vector3<T> start;
	Vector3<T> end;
};

template<typename C, typename T, typename U>
	requires std::floating_point<U>
inline std::basic_istream<C, T>& operator>>(std::basic_istream<C, T>& s, Segment3<U>& segment)
{ 
	return s >> segment.start >> std::ws >> segment.end;
}

template<typename C, typename T, typename U>
	requires std::floating_point<U>
inline std::basic_ostream<C, T>& operator<<(std::basic_ostream<C, T>& s, const Segment3<U>& segment)
{ 
	constexpr C WS(0x20);
	return s << segment.start << WS << segment.end;
}

template<typename T>
inline bool Segment3<T>::approxEquals(const Segment3<T>& segment) const
{
	return origin.approxEquals(segment.start) && direction.approxEquals(segment.end);
}

template<typename T>
inline bool Segment3<T>::approxEquals(const Segment3<T>& segment, T tolerance) const
{
	return origin.approxEquals(segment.start, tolerance) && direction.approxEquals(segment.end, tolerance);
}

template<typename T>
template<std::output_iterator<Vector2<T>> O>
inline O Segment3<T>::copyEndpoints(O target) const
{
	*target++ = start;
	*target++ = end;
	return target;
}

template<typename T>
inline Segment3<T>& Segment3<T>::transform(const Matrix3<T>& matrix)
{
	start *= matrix;
	end *= matrix;
	return *this;
}

template<typename T>
inline Segment3<T>& Segment3<T>::transform(const AffineTransform<T>& transformation)
{
	start.transform(transformation);
	end.transform(transformation);
	return *this;
}

template<typename T>
inline Vector3<T> Segment3<T>::getClosestPoint(const Vector3<T>& point) const
{
	Vector3<T> direction = end - start;
	return std::clamp(dot(point - start, direction)/dot(direction, direction), T(0), T(1))*direction + start;
}

} // namespace templates

#if MATHEMATICS_DOUBLE
using Segment3 = templates::Segment3<double>;
using Segment3Arg = templates::Segment3<double>::ConstArg;
using Segment3Result = templates::Segment3<double>::ConstResult;
#else
using Segment3 = templates::Segment3<float>;
using Segment3Arg = templates::Segment3<float>::ConstArg;
using Segment3Result = templates::Segment3<float>::ConstResult;
#endif

} // namespace mathematics

namespace std {

template<typename T>
struct hash;

template<typename T>
struct hash<::mathematics::templates::Segment3<T>>
{
	std::size_t operator()(const ::mathematics::templates::Segment3<T>& segment) const noexcept
	{
		std::hash<T> hasher;
		std::size_t seed = hasher(segment.start) + 0x9e3779b9;
		seed ^= hasher(segment.end) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
		return seed;
	}
};

} // namespace std

#include "HalfSpace.hpp"
#include "Plane.hpp"
#include "Triangle3.hpp"
#include "AxisAlignedBox.hpp"
#include "OrientedBox.hpp"
#include "Sphere.hpp"
#include "Ellipsoid.hpp"
#include "Intersections.inl"

namespace mathematics::templates {

template<typename T>
inline bool Segment3<T>::intersects(const HalfSpace<T>& halfSpace) const
{
	return halfSpace.contains(start) || halfSpace.contains(end);
}

template<typename T>
inline bool Segment3<T>::intersects(const Ellipsoid<T>& ellipsoid) const
{
	return intersections::testSegmentEllipsoid(start, end, ellipsoid.center, ellipsoid.getMatrix());
}

template<typename T>
inline std::optional<T> Segment3<T>::findIntersection(const Plane<T>& plane) const
{
	auto result = intersections::findLinePlane<std::optional<T>>(start, end - start, plane.getNormal(), plane.d);
	return (result.has_value() && (result.value() >= T(0)) && (result.value() <= T(1))) ? result : {};
}

//template<typename T>
//template<ScalarOrVector3<T> U>
//inline std::optional<U> Segment3<T>::findIntersection(const Plane<T>& plane) const
//{
//	std::optional<T> result = findIntersection(plane);
//	if constexpr (std::is_same_v<U, T>)
//		return result;
//	else //if constexpr (std::is_same_v<U, Vector3<T>>)
//		return result.has_value() ? { evaluate(result.value()) } : {};
//}

template<typename T>
inline std::optional<T> Segment3<T>::findIntersection(const Triangle3<T>& triangle) const
{
	//return intersections::findSegmentTriangle<std::optional<T>>(start, end, triangle.vertices[0], triangle.vertices[1],
	//	triangle.vertices[2]);
	auto result = intersections::findLineTriangle<std::optional<T>>(start, end - start, triangle.vertices[0], triangle.vertices[1],
		triangle.vertices[2]);
	return (result.has_value() && (result.value() >= T(0)) && (result.value() <= T(1))) ? result : {};
}

template<typename T>
inline std::optional<Interval<T>> Segment3<T>::findIntersection(const AxisAlignedBox<T>& box) const
{
	auto result = intersections::findLineAxisAlignedBox<std::optional<Interval<T>>>(start, end - start, box.minimum, box.maximum);
	
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
inline std::optional<Interval<T>> Segment3<T>::findIntersection(const OrientedBox<T>& box) const
{
	auto result = intersections::findLineOrientedBox<std::optional<Interval<T>>>(start, end - start, box.center, box.basis, box.halfDims);
	
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
inline std::optional<Interval<T>> Segment3<T>::findIntersection(const Sphere<T>& sphere) const
{
	auto result = intersections::findLineNSphere<std::optional<Interval<T>>>(start, end - start, sphere.center, sphere.radius);

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
inline std::optional<Interval<T>> Segment3<T>::findIntersection(const Ellipsoid<T>& ellipsoid) const
{
	return intersections::findSegmentEllipsoid<std::optional<Interval<T>>>(start, end, ellipsoid.center, ellipsoid.getMatrix());
}

} // namespace mathematics::templates
