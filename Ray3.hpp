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
#include <algorithm>
#include <cstddef>
#include <cmath>
#include "Vector3.hpp"
#include "Matrix3.hpp"
#include "AffineTransform.hpp"
#include "Line3.hpp"

namespace core::mathematics {
namespace templates {

template<typename T>
	requires std::floating_point<T>
struct Ray3
{
	using Real = T;
	using ConstArg = const Ray3&;
	using ConstResult = const Ray3&;

	Ray3() = default;
	Ray3(const Vector3<T>& origin, const Vector3<T>& direction) noexcept : origin(origin), direction(direction) {}
	explicit Ray3(const Line3<T>& line) noexcept : origin(line.origin), direction(line.direction) {}
	//explicit Ray3(const Segment3<T>& segment) noexcept;

	//Vector3<T> operator()(T t) const noexcept { return (origin + t*direction); }
	bool operator==(const Ray3& ray) const noexcept { return (origin == ray.origin) && (direction == ray.direction); }
	bool operator!=(const Ray3& ray) const noexcept { return !(*this == ray); }

	template<typename A> void serialize(A& ar) { ar(origin, direction); }

	const Line3<T>& asLine() const noexcept { return reinterpret_cast<const Line3<T>&>(*this); }

	// Properties
	bool isApproxEqual(const Ray3& ray) const noexcept;
	bool isApproxEqual(const Ray3& ray, T tolerance) const noexcept;
	//bool isFinite() const noexcept { return origin.isFinite() && direction.isFinite(); }
	void set(const Vector3<T>& origin, const Vector3<T>& direction) noexcept { this->origin = origin; this->direction = direction; }
	const Vector3<T>& getOrigin() const noexcept { return origin; }
	void setOrigin(const Vector3<T>& origin) noexcept { this->origin = origin; }
	const Vector3<T>& getDirection() const noexcept { return direction; }
	void setDirection(const Vector3<T>& direction) noexcept { this->direction = direction; }

	// Transformation
	void translate(const Vector3<T>& offset) noexcept { origin += offset; }
	void transform(const Matrix3<T>& matrix) noexcept;
	void transform(const AffineTransform<T>& transformation) noexcept;
	void normalize() noexcept { direction.normalize(); }

	// Evaluation
	Vector3<T> evaluate(T t) const noexcept { return (origin + t*direction); }

	// Closest points
	Vector3<T> getClosestPoint(const Vector3<T>& point) const;											// normalized ray
	template<typename U> Vector3<T> getClosestPoint(const Vector3<T>& point) const;
	T getDistance(const Vector3<T>& point) const { return distance(getClosestPoint(point), point); }	// normalized ray
	template<typename U> T getDistance(const Vector3<T>& point) const { return distance(getClosestPoint<U>(point), point); }

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
inline bool Ray3<T>::isApproxEqual(const Ray3<T>& ray) const
{
	return origin.isApproxEqual(ray.origin) && direction.isApproxEqual(ray.direction);
}

template<typename T>
inline bool Ray3<T>::isApproxEqual(const Ray3<T>& ray, T tolerance) const
{
	return origin.isApproxEqual(ray.origin, tolerance) && direction.isApproxEqual(ray.direction, tolerance);
}

template<typename T>
inline void Ray3<T>::transform(const Matrix3<T>& matrix)
{
	origin *= matrix;
	direction *= matrix;
}

template<typename T>
inline void Ray3<T>::transform(const AffineTransform<T>& transformation)
{
	origin.transform(transformation);
	direction *= transformation.getBasis();
}

template<typename T>
inline Vector3<T> Ray3<T>::getClosestPoint(const Vector3<T>& point) const
{
	return std::max(dot(point - origin, direction), T(0))*direction + origin;
}

template<typename T>
template<typename U>
inline Vector3<T> Ray3<T>::getClosestPoint(const Vector3<T>& point) const
{
	static_assert(std::is_same_v<U, Normalized> || std::is_same_v<U, Unnormalized>);
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

} // namespace core::mathematics

namespace std {

template<typename T>
struct hash;

template<typename T>
struct hash<::core::mathematics::templates::Ray3<T>>
{
	std::size_t operator()(const ::core::mathematics::templates::Ray3<T>& ray) const noexcept
	{
		std::hash<T> hasher;
		std::size_t seed = hasher(ray.origin) + 0x9e3779b9;
		seed ^= hasher(ray.direction) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
		return seed;
	}
};

} // namespace std
