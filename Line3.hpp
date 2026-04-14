/*
 *	Name: Line3
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
#include <optional>
#include <iterator>
#include <cstddef>
#include <cmath>
#include "Vector3.hpp"
#include "Matrix3.hpp"
#include "AffineTransform.hpp"
#include "Interval.hpp"

namespace mathematics {
namespace templates {

//template<typename T, typename U>
//concept ScalarOrVector3 = (std::same_as<T, U> || std::same_as<T, Vector3<U>>); // #TODO Move to Concepts.inl

//template<typename T, typename U>
//concept IntervalOrLineSegment3 = (std::same_as<T, Interval<U>> || std::same_as<T, LineSegment3<U>>); // #TODO Move to Concepts.inl

template<typename T>
	requires std::floating_point<T>
struct Ray3;

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
struct Plane<float>;

#endif /* SIMD_HAS_FLOAT4 */

template<typename T>
	requires std::floating_point<T>
struct Line3
{
	using Real = T;
	using ConstArg = const Line3&;
	using ConstResult = const Line3&;

	Line3() = default;
	explicit Line3(Uninitialized) noexcept : origin(Uninitialized()), direction(Uninitialized()) {}
	Line3(const Vector3<T>& origin, const Vector3<T>& direction) noexcept : origin(origin), direction(direction) {}
	Line3(const Ray3<T>& ray) noexcept;
	//explicit Line3(const LineSegment3<T>& segment) noexcept;

	//Vector3<T> operator()(T t) const noexcept { return (origin + t*direction); }
	bool operator==(const Line3& line) const noexcept { return (origin == line.origin) && (direction == line.direction); }
	bool operator!=(const Line3& line) const noexcept { return !(*this == line); }

	template<typename A> void serialize(A& ar) { ar(origin, direction); }

	const Ray3<T>& asRay() const noexcept;

	// Least-squares fit of a line
	//template<std::input_iterator I, std::sentinel_for<I> S> static Line3 computeBestFit(I first, S last); // #TODO

	// Properties
	bool approxEquals(const Line3& line) const noexcept;
	bool approxEquals(const Line3& line, T tolerance) const noexcept;
	//bool isFinite() const noexcept { return origin.isFinite() && direction.isFinite(); }
	Line3& set(const Vector3<T>& origin, const Vector3<T>& direction) noexcept { this->origin = origin; this->direction = direction; return *this; }
	const Vector3<T>& getOrigin() const noexcept { return origin; }
	void setOrigin(const Vector3<T>& origin) noexcept { this->origin = origin; }
	const Vector3<T>& getDirection() const noexcept { return direction; }
	void setDirection(const Vector3<T>& direction) noexcept { this->direction = direction; }

	// Transformation
	Line3& translate(const Vector3<T>& offset) noexcept { origin += offset; return *this; }
	Line3& transform(const Matrix3<T>& matrix) noexcept;
	Line3& transform(const AffineTransform<T>& transformation) noexcept;
	Line3& normalize() noexcept { direction.normalize(); return *this; }

	// Evaluation
	Vector3<T> evaluate(T t) const noexcept { return (origin + t*direction); }

	// Closest points
	Vector3<T> getClosestPoint(const Vector3<T>& point) const noexcept;											// normalized line
	template<Normalization U> Vector3<T> getClosestPoint(const Vector3<T>& point) const noexcept;
	T getDistanceTo(const Vector3<T>& point) const noexcept { return distance(getClosestPoint(point), point); }	// normalized line
	template<Normalization U> T getDistanceTo(const Vector3<T>& point) const noexcept { return distance(getClosestPoint<U>(point), point); }
	//T getDistanceTo(const Line3& line) const noexcept;

	// Intersection
	bool intersects(const Plane<T>& plane) const noexcept;
	bool intersects(const Triangle3<T>& triangle) const noexcept { return findIntersection(triangle).has_value(); }
	bool intersects(const AxisAlignedBox<T>& box) const noexcept { return findIntersection(box).has_value(); }
	bool intersects(const OrientedBox<T>& box) const noexcept { return findIntersection(box).has_value(); }
	bool intersects(const Sphere<T>& sphere) const noexcept;
	template<Normalization U> bool intersects(const Sphere<T>& sphere) const noexcept;
	bool intersects(const Ellipsoid<T>& ellipsoid) const noexcept;
	std::optional<T> findIntersection(const Plane<T>& plane) const noexcept;
	std::optional<T> findIntersection(const Triangle3<T>& triangle) const noexcept;
	//template<ScalarOrVector3<T> U> std::optional<U> findIntersection(const Plane<T>& plane) const noexcept;
	std::optional<Interval<T>> findIntersection(const AxisAlignedBox<T>& box) const noexcept;
	std::optional<Interval<T>> findIntersection(const OrientedBox<T>& box) const noexcept;
	std::optional<Interval<T>> findIntersection(const Sphere<T>& sphere) const noexcept;
	template<Normalization U> std::optional<Interval<T>> findIntersection(const Sphere<T>& sphere) const noexcept;
	std::optional<Interval<T>> findIntersection(const Ellipsoid<T>& ellipsoid) const noexcept;

	Vector3<T> origin;
	Vector3<T> direction;
};

template<typename C, typename T, typename U>
	requires std::floating_point<U>
inline std::basic_istream<C, T>& operator>>(std::basic_istream<C, T>& s, Line3<U>& line)
{ 
	return s >> line.origin >> std::ws >> line.direction;
}

template<typename C, typename T, typename U>
	requires std::floating_point<U>
inline std::basic_ostream<C, T>& operator<<(std::basic_ostream<C, T>& s, const Line3<U>& line)
{ 
	constexpr C WS(0x20);
	return s << line.origin << WS << line.direction;
}

template<typename T>
inline bool Line3<T>::approxEquals(const Line3<T>& line) const
{
	return origin.approxEquals(line.origin) && direction.approxEquals(line.direction);
}

template<typename T>
inline bool Line3<T>::approxEquals(const Line3<T>& line, T tolerance) const
{
	return origin.approxEquals(line.origin, tolerance) && direction.approxEquals(line.direction, tolerance);
}

template<typename T>
inline Line3<T>& Line3<T>::transform(const Matrix3<T>& matrix)
{
	origin *= matrix;
	direction *= matrix;
	return *this;
}

template<typename T>
inline Line3<T>& Line3<T>::transform(const AffineTransform<T>& transformation)
{
	origin.transform(transformation);
	direction *= transformation.getBasis();
	return *this;
}

template<typename T>
inline Vector3<T> Line3<T>::getClosestPoint(const Vector3<T>& point) const
{
	return dot(point - origin, direction)*direction + origin;
}

template<typename T>
template<Normalization U>
inline Vector3<T> Line3<T>::getClosestPoint(const Vector3<T>& point) const
{
	if costexpr(std::is_same_v<U, Normalized>)
		return dot(point - origin, direction)*direction + origin;
	else
		return (dot(point - origin, direction)/dot(direction, direction))*direction + origin;
}

template<typename T>
	requires std::floating_point<T>
inline Line3<T> normalize(const Line3<T>& line) noexcept
{
	Line3<T> l(line);
	l.normalize();
	return l;
}

template<typename T>
	requires std::floating_point<T>
inline Line3<T> normalize(Line3<T>&& line) noexcept
{
	line.normalize();
	return line;
}

} // namespace templates

#if MATHEMATICS_DOUBLE
using Line3 = templates::Line3<double>;
using Line3Arg = templates::Line3<double>::ConstArg;
using Line3Result = templates::Line3<double>::ConstResult;
#else
using Line3 = templates::Line3<float>;
using Line3Arg = templates::Line3<float>::ConstArg;
using Line3Result = templates::Line3<float>::ConstResult;
#endif

} // namespace mathematics

namespace std {

template<typename T>
struct hash;

template<typename T>
struct hash<::mathematics::templates::Line3<T>>
{
	size_t operator()(const ::mathematics::templates::Line3<T>& line) const noexcept
	{
		hash<typename ::mathematics::templates::Vector3<T>> hasher;
		size_t seed = hasher(line.origin) + 0x9e3779b9;
		seed ^= hasher(line.direction) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
		return seed;
	}
};

} // namespace std

#include "Ray3.hpp"
#include "Plane.hpp"
#include "Triangle3.hpp"
#include "AxisAlignedBox.hpp"
#include "OrientedBox.hpp"
#include "Sphere.hpp"
#include "Ellipsoid.hpp"
#include "Intersections.inl"

namespace mathematics::templates {

template<typename T>
inline Line3<T>::Line3(const Ray3<T>& ray) : origin(ray.origin), direction(ray.direction)
{
}

template<typename T>
inline const Ray3<T>& Line3<T>::asRay() const
{ 
	return reinterpret_cast<const Ray3<T>&>(*this);
}

template<typename T>
inline bool Line3<T>::intersects(const Plane<T>& plane) const
{
	return intersections::testLinePlane(direction, plane.getNormal());
}

template<typename T>
inline bool Line3<T>::intersects(const Sphere<T>& sphere) const
{
	return intersections::testLineNSphere(origin, direction, sphere.center, sphere.radius);
}

template<typename T>
template<Normalization U>
inline bool Line3<T>::intersects(const Sphere<T>& sphere) const
{
	if costexpr(std::is_same_v<U, Normalized>)
		return intersections::testNormalizedLineNSphere(origin, direction, sphere.center, sphere.radius);
	else
		return intersections::testLineNSphere(origin, direction, sphere.center, sphere.radius);
}

template<typename T>
inline bool Line3<T>::intersects(const Ellipsoid<T>& ellipsoid) const
{
	return intersections::testLineEllipsoid(origin, direction, ellipsoid.center, ellipsoid.getMatrix());
}

template<typename T>
inline std::optional<T> Line3<T>::findIntersection(const Plane<T>& plane) const
{
	return intersections::findLinePlane<std::optional<T>>(origin, direction, plane.getNormal(), plane.d);
}

//template<typename T>
//template<ScalarOrVector3<T> U>
//inline std::optional<U> Line3<T>::findIntersection(const Plane<T>& plane) const
//{
//	std::optional<T> result = findIntersection(plane);
//	if constexpr (std::is_same_v<U, T>)
//		return result;
//	else //if constexpr (std::is_same_v<U, Vector3<T>>)
//		return result.has_value() ? { evaluate(result.value()) } : {};
//}

template<typename T>
inline std::optional<T> Line3<T>::findIntersection(const Triangle3<T>& triangle) const
{
	return intersections::findLineTriangle<std::optional<T>>(origin, direction, triangle.vertices[0], triangle.vertices[1],
		triangle.vertices[2]);
}

template<typename T>
inline std::optional<Interval<T>> Line3<T>::findIntersection(const AxisAlignedBox& box) const
{
	return intersections::findLineAxisAlignedBox<std::optional<Interval<T>>>(origin, direction, box.minimum, box.maximum);
}

template<typename T>
inline std::optional<Interval<T>> Line3<T>::findIntersection(const OrientedBox& box) const
{
	return intersections::findLineOrientedBox<std::optional<Interval<T>>>(origin, direction, box.center, box.basis, box.halfDims);
}

template<typename T>
inline std::optional<Interval<T>> Line3<T>::findIntersection(const Sphere<T>& sphere) const
{
	return intersections::findLineNSphere<std::optional<Interval<T>>>(origin, direction, sphere.center, sphere.radius);
}

template<typename T>
template<Normalization U>
inline std::optional<Interval<T>> Line3<T>::findIntersection(const Sphere<T>& sphere) const
{
	if costexpr(std::is_same_v<U, Normalized>)
		return intersections::findNormalizedLineNSphere<std::optional<Interval<T>>>(origin, direction, sphere.center, sphere.radius);
	else
		return intersections::findLineNSphere<std::optional<Interval<T>>>(origin, direction, sphere.center, sphere.radius);
}

template<typename T>
inline std::optional<Interval<T>> Line3<T>::findIntersection(const Ellipsoid<T>& ellipsoid) const
{
	return intersections::findLineEllipsoid<std::optional<Interval<T>>>(origin, direction, ellipsoid.center, ellipsoid.getMatrix());
}

} // namespace mathematics::templates
