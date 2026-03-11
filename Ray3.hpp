/*
 *	Name: Ray3
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
#include "Vector3.hpp"
#include "Matrix3.hpp"
#include "AffineTransform.hpp"
#include "Interval.hpp"
#include "Line3.hpp"

namespace mathematics {
namespace templates {

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
struct Ray3
{
	using Real = T;
	using ConstArg = const Ray3&;
	using ConstResult = const Ray3&;

	Ray3() = default;
	explicit Ray3(Uninitialized) noexcept : origin(Uninitialized()), direction(Uninitialized()) {}
	Ray3(const Vector3<T>& origin, const Vector3<T>& direction) noexcept : origin(origin), direction(direction) {}
	explicit Ray3(const Line3<T>& line) noexcept : origin(line.origin), direction(line.direction) {}
	//explicit Ray3(const Segment3<T>& segment) noexcept;

	//Vector3<T> operator()(T t) const noexcept { return (origin + t*direction); }
	bool operator==(const Ray3& ray) const noexcept { return (origin == ray.origin) && (direction == ray.direction); }
	bool operator!=(const Ray3& ray) const noexcept { return !(*this == ray); }

	template<typename A> void serialize(A& ar) { ar(origin, direction); }

	const Line3<T>& asLine() const noexcept { return reinterpret_cast<const Line3<T>&>(*this); }

	// Properties
	bool approxEquals(const Ray3& ray) const noexcept;
	bool approxEquals(const Ray3& ray, T tolerance) const noexcept;
	//bool isFinite() const noexcept { return origin.isFinite() && direction.isFinite(); }
	Ray3& set(const Vector3<T>& origin, const Vector3<T>& direction) noexcept { this->origin = origin; this->direction = direction; return *this; }
	const Vector3<T>& getOrigin() const noexcept { return origin; }
	void setOrigin(const Vector3<T>& origin) noexcept { this->origin = origin; }
	const Vector3<T>& getDirection() const noexcept { return direction; }
	void setDirection(const Vector3<T>& direction) noexcept { this->direction = direction; }

	// Transformation
	Ray3& translate(const Vector3<T>& offset) noexcept { origin += offset; return *this; }
	Ray3& transform(const Matrix3<T>& matrix) noexcept;
	Ray3& transform(const AffineTransform<T>& transformation) noexcept;
	Ray3& normalize() noexcept { direction.normalize(); return *this; }

	// Evaluation
	Vector3<T> evaluate(T t) const noexcept { return (origin + t*direction); }

	// Closest points
	Vector3<T> getClosestPoint(const Vector3<T>& point) const;											// normalized ray
	template<Normalization U> Vector3<T> getClosestPoint(const Vector3<T>& point) const;
	T getDistanceTo(const Vector3<T>& point) const { return distance(getClosestPoint(point), point); }	// normalized ray
	template<Normalization U> T getDistanceTo(const Vector3<T>& point) const { return distance(getClosestPoint<U>(point), point); }

	// Intersection
	bool intersects(const Plane<T>& plane) const noexcept { return findIntersection(plane).has_value(); }
	bool intersects(const Triangle3<T>& triangle) const noexcept { return findIntersection(triangle).has_value(); }
	bool intersects(const AxisAlignedBox& box) const noexcept { return findIntersection(box).has_value(); }
	bool intersects(const OrientedBox& box) const noexcept { return findIntersection(box).has_value(); }
	bool intersects(const Sphere<T>& sphere) const noexcept { return findIntersection(sphere).has_value(); }
	template<Normalization U> bool intersects(const Sphere<T>& sphere) const noexcept { return findIntersection<U>(sphere).has_value(); }
	bool intersects(const Ellipsoid<T>& ellipsoid) const noexcept; // #TODO
	std::optional<T> findIntersection(const Plane<T>& plane) const noexcept;
	std::optional<T> findIntersection(const Triangle3<T>& triangle) const noexcept; // #TODO
	//template<ScalarOrVector3<T> U> std::optional<U> findIntersection(const Plane<T>& plane) const noexcept;
	std::optional<Interval<T>> findIntersection(const AxisAlignedBox& box) const noexcept;
	std::optional<Interval<T>> findIntersection(const OrientedBox& box) const noexcept;
	std::optional<Interval<T>> findIntersection(const Sphere<T>& sphere) const noexcept;
	template<Normalization U> std::optional<Interval<T>> findIntersection(const Sphere<T>& sphere) const noexcept;
	std::optional<Interval<T>> findIntersection(const Ellipsoid<T>& ellipsoid) const noexcept; // #TODO

	Vector3<T> origin;
	Vector3<T> direction;
};

template<typename C, typename T, typename U>
	requires std::floating_point<U>
inline std::basic_istream<C, T>& operator>>(std::basic_istream<C, T>& s, Ray3<U>& ray)
{ 
	return s >> ray.origin >> std::ws >> ray.direction;
}

template<typename C, typename T, typename U>
	requires std::floating_point<U>
inline std::basic_ostream<C, T>& operator<<(std::basic_ostream<C, T>& s, const Ray3<U>& ray)
{ 
	constexpr C WS(0x20);
	return s << ray.origin << WS << ray.direction;
}

template<typename T>
inline bool Ray3<T>::approxEquals(const Ray3<T>& ray) const
{
	return origin.approxEquals(ray.origin) && direction.approxEquals(ray.direction);
}

template<typename T>
inline bool Ray3<T>::approxEquals(const Ray3<T>& ray, T tolerance) const
{
	return origin.approxEquals(ray.origin, tolerance) && direction.approxEquals(ray.direction, tolerance);
}

template<typename T>
inline Ray3<T>& Ray3<T>::transform(const Matrix3<T>& matrix)
{
	origin *= matrix;
	direction *= matrix;
	return *this;
}

template<typename T>
inline Ray3<T>& Ray3<T>::transform(const AffineTransform<T>& transformation)
{
	origin.transform(transformation);
	direction *= transformation.getBasis();
	return *this;
}

template<typename T>
inline Vector3<T> Ray3<T>::getClosestPoint(const Vector3<T>& point) const
{
	return std::max(dot(point - origin, direction), T(0))*direction + origin;
}

template<typename T>
template<Normalization U>
inline Vector3<T> Ray3<T>::getClosestPoint(const Vector3<T>& point) const
{
	if costexpr(std::is_same_v<U, Normalized>)
		return std::max(dot(point - origin, direction), T(0))*direction + origin;
	else
		return std::max(dot(point - origin, direction)/dot(direction, direction), T(0))*direction + origin;
}

template<typename T>
	requires std::floating_point<T>
inline Ray3<T> normalize(const Ray3<T>& ray) noexcept
{
	Ray3<T> r(ray);
	r.normalize();
	return r;
}

template<typename T>
	requires std::floating_point<T>
inline Ray3<T> normalize(Ray3<T>&& ray) noexcept
{
	ray.normalize();
	return ray;
}

} // namespace templates

#if MATHEMATICS_DOUBLE
using Ray3 = templates::Ray3<double>;
using Ray3Arg = templates::Ray3<double>::ConstArg;
using Ray3Result = templates::Ray3<double>::ConstResult;
#else
using Ray3 = templates::Ray3<float>;
using Ray3Arg = templates::Ray3<float>::ConstArg;
using Ray3Result = templates::Ray3<float>::ConstResult;
#endif

} // namespace mathematics

namespace std {

template<typename T>
struct hash;

template<typename T>
struct hash<::mathematics::templates::Ray3<T>>
{
	std::size_t operator()(const ::mathematics::templates::Ray3<T>& ray) const noexcept
	{
		std::hash<T> hasher;
		std::size_t seed = hasher(ray.origin) + 0x9e3779b9;
		seed ^= hasher(ray.direction) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
		return seed;
	}
};

} // namespace std

#include "Plane.hpp"
#include "Triangle3.hpp"
#include "AxisAlignedBox.hpp"
#include "OrientedBox.hpp"
#include "Sphere.hpp"
#include "Ellipsoid.hpp"
#include "Intersections.inl"

namespace mathematics::templates {

template<typename T>
inline std::optional<T> Ray3<T>::findIntersection(const Plane<T>& plane) const
{
	std::optional<T> result = intersections::findLinePlane<std::optional<T>>(origin, direction, plane.getNormal(), plane.d);
	return (result.has_value() && (result.value() >= T(0))) ? result : {};
}

//template<typename T>
//template<ScalarOrVector3<T> U>
//inline std::optional<U> Ray3<T>::findIntersection(const Plane<T>& plane) const
//{
//	std::optional<T> result = findIntersection(plane);
//	if constexpr (std::is_same_v<U, T>)
//		return result;
//	else //if constexpr (std::is_same_v<U, Vector3<T>>)
//		return result.has_value() ? { evaluate(result.value()) } : {};
//}

template<typename T>
inline std::optional<Interval<T>> Ray3<T>::findIntersection(const AxisAlignedBox<T>& box) const
{
	std::optional<Interval<T>> result = intersections::findLineAxisAlignedBox<std::optional<Interval<T>>>(origin, direction,
		box.minimum, box.maximum);

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
inline std::optional<Interval<T>> Ray3<T>::findIntersection(const OrientedBox<T>& box) const
{
	//Matrix3<T> basisTranspose(transpose(box.basis));
	std::optional<Interval<T>> result = intersections::findLineAxisAlignedBox<std::optional<Interval<T>>>(box.basis*(origin - box.center)/*(origin - box.center)*basisTranspose*/,
		box.basis*direction/*direction*basisTranspose*/, -box.halfDims, box.halfDims);

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
inline std::optional<Interval<T>> Ray3<T>::findIntersection(const Sphere<T>& sphere) const
{
	std::optional<Interval<T>> result = intersections::findLineNSphere<std::optional<Interval<T>>>(origin, direction, sphere.center, sphere.radius);
	
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
inline std::optional<Interval<T>> Ray3<T>::findIntersection(const Sphere<T>& sphere) const
{
	std::optional<Interval<T>> result;
	if costexpr(std::is_same_v<U, Normalized>)
		result = intersections::findNormalizedLineNSphere<std::optional<Interval<T>>>(origin, direction, sphere.center, sphere.radius);
	else
		result = intersections::findLineNSphere<std::optional<Interval<T>>>(origin, direction, sphere.center, sphere.radius);

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
