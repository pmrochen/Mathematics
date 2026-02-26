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
#include "Line2.hpp"

namespace core::mathematics {
namespace templates {

template<typename T>
	requires std::floating_point<T>
struct Ray2
{
	using Real = T;
	using ConstArg = const Ray2&;
	using ConstResult = const Ray2&;

	Ray2() = default;
	Ray2(const Vector2<T>& origin, const Vector2<T>& direction) noexcept : origin(origin), direction(direction) {}
	explicit Ray2(const Line2<T>& line) noexcept : origin(line.origin), direction(line.direction) {}
	//explicit Ray2(const Segment2<T>& segment) noexcept;

	//Vector3<T> operator()(T t) const noexcept { return (origin + t*direction); }
	bool operator==(const Ray2& ray) const noexcept { return (origin == ray.origin) && (direction == ray.direction); }
	bool operator!=(const Ray2& ray) const noexcept { return !(*this == ray); }

	template<typename A> void serialize(A& ar) { ar(origin, direction); }

	const Line2<T>& asLine() const noexcept { return reinterpret_cast<const Line2<T>&>(*this); }

	// Properties
	bool isApproxEqual(const Ray2& ray) const noexcept;
	bool isApproxEqual(const Ray2& ray, T tolerance) const noexcept;
	//bool isFinite() const noexcept { return origin.isFinite() && direction.isFinite(); }
	void set(const Vector2<T>& origin, const Vector2<T>& direction) noexcept { this->origin = origin; this->direction = direction; }
	const Vector2<T>& getOrigin() const noexcept { return origin; }
	void setOrigin(const Vector2<T>& origin) noexcept { this->origin = origin; }
	const Vector2<T>& getDirection() const noexcept { return direction; }
	void setDirection(const Vector2<T>& direction) noexcept { this->direction = direction; }
	T getSlope() const noexcept { return direction.y/direction.x; }
	T getInclinationAngle() const noexcept { return std::atan2(direction.y, direction.x); }

	// Normalization
	void normalize() noexcept { direction.normalize(); }

	// Evaluation
	Vector2<T> evaluate(T t) const noexcept { return (origin + t*direction); }

	// Closest points
	Vector2<T> getClosestPoint(const Vector2<T>& point) const;											// normalized ray
	template<typename U> Vector2<T> getClosestPoint(const Vector2<T>& point) const;
	T getDistance(const Vector2<T>& point) const { return distance(getClosestPoint(point), point); }	// normalized ray
	template<typename U> T getDistance(const Vector2<T>& point) const { return distance(getClosestPoint<U>(point), point); }

	// Intersection
	bool testIntersection(const Line2<T>& line) const noexcept { return findIntersection(line).has_value(); }
	//bool testIntersection(const Ray2& ray) const noexcept { return findIntersection(ray).has_value(); }
	std::optional<T> findIntersection(const Line2<T>& line) const noexcept;
	//std::optional<T> findIntersection(const Ray2& ray) const noexcept;
	template<typename U> std::optional<U> findIntersection(const Line2<T>& line) const noexcept;
	//template<typename U> std::optional<U> findIntersection(const Ray2& ray) const noexcept;

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
inline bool Ray2<T>::isApproxEqual(const Ray2<T>& ray) const
{
	return origin.isApproxEqual(ray.origin) && direction.isApproxEqual(ray.direction);
}

template<typename T>
inline bool Ray2<T>::isApproxEqual(const Ray2<T>& ray, T tolerance) const
{
	return origin.isApproxEqual(ray.origin, tolerance) && direction.isApproxEqual(ray.direction, tolerance);
}

template<typename T>
inline Vector2<T> Ray2<T>::getClosestPoint(const Vector2<T>& point) const
{
	return std::max(dot(point - origin, direction), T(0))*direction + origin;
}

template<typename T>
template<typename U>
inline Vector2<T> Ray2<T>::getClosestPoint(const Vector2<T>& point) const
{
	static_assert(std::is_same_v<U, Normalized> || std::is_same_v<U, Unnormalized>);
	if costexpr(std::is_same_v<U, Normalized>)
		return std::max(dot(point - origin, direction), T(0))*direction + origin;
	else
		return std::max(dot(point - origin, direction)/dot(direction, direction), T(0))*direction + origin;
}

template<typename T>
inline std::optional<T> Ray2<T>::findIntersection(const Line2<T>& line) const
{
	T d1CrossD2 = cross(direction, line.direction);
	if (std::fabs(d1CrossD2) < Constants<T>::TOLERANCE)
	{
		if (!(std::fabs(cross(normalize(line.origin - origin), direction)) < Constants<T>::TOLERANCE))
			return {};
		return std::optional<T>(T(0)); // collinear
	}

	T t = cross(line.origin - origin, line.direction)/d1CrossD2;
	//T u = cross(line.origin - origin, direction)/d1CrossD2;
	if (t >= T(0))
		return std::optional<T>(t);
	return {};
}

template<typename T>
template<typename U> 
inline std::optional<U> Ray2<T>::findIntersection(const Line2<T>& line) const
{
	static_assert(std::is_same_v<U, T> || std::is_same_v<U, Vector3<T>>);
	std::optional<T> result = findIntersection(line);
	if constexpr (std::is_same_v<U, T>)
		return result;
	else //if constexpr (std::is_same_v<U, Vector3<T>>)
		return result.has_value() ? std::optional<U>(line.evaluate(result.value())) : std::optional<U>();
}

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

} // namespace core::mathematics

namespace std {

template<typename T>
struct hash;

template<typename T>
struct hash<::core::mathematics::templates::Ray2<T>>
{
	std::size_t operator()(const ::core::mathematics::templates::Ray2<T>& ray) const noexcept
	{
		std::hash<T> hasher;
		std::size_t seed = hasher(ray.origin) + 0x9e3779b9;
		seed ^= hasher(ray.direction) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
		return seed;
	}
};

} // namespace std
