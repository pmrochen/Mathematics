/*
 *	Name: Line3
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

namespace core::mathematics {
namespace templates {

template<typename T>
	requires std::floating_point<T>
struct Ray;

template<typename T>
	requires std::floating_point<T>
struct Line3
{
	using Real = T;
	using ConstArg = const Line3&;
	using ConstResult = const Line3&;

	Line3() = default;
	Line3(const Vector3<T>& origin, const Vector3<T>& direction) noexcept : origin(origin), direction(direction) {}
	Line3(const Ray3<T>& ray) noexcept;
	//explicit Line3(const Segment3<T>& segment) noexcept;

	//Vector3<T> operator()(T t) const noexcept { return (origin + t*direction); }
	bool operator==(const Line3& line) const noexcept { return (origin == line.origin) && (direction == line.direction); }
	bool operator!=(const Line3& line) const noexcept { return !(*this == line); }

	template<typename A> void serialize(A& ar) { ar(origin, direction); }

	const Ray3<T>& asRay() const noexcept;

	// Least-squares fit of a line
	//static Line3 computeBestFitOf(const Vector3<T>* points, std::size_t nPoints) noexcept; // #TODO
	//static Line3 computeBestFitOf(const std::vector<Vector3<T>>& points) noexcept;

	// Properties
	bool isApproxEqual(const Line3& line) const noexcept;
	bool isApproxEqual(const Line3& line, T tolerance) const noexcept;
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
	Vector3<T> getClosestPoint(const Vector3<T>& point) const noexcept;											// normalized line
	template<typename U> Vector3<T> getClosestPoint(const Vector3<T>& point) const noexcept;
	T getDistance(const Vector3<T>& point) const noexcept { return distance(getClosestPoint(point), point); }	// normalized line
	template<typename U> T getDistance(const Vector3<T>& point) const noexcept { return distance(getClosestPoint<U>(point), point); }
	//T getDistance(const Line3& line) const noexcept;

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
inline bool Line3<T>::isApproxEqual(const Line3<T>& line) const
{
	return origin.isApproxEqual(line.origin) && direction.isApproxEqual(line.direction);
}

template<typename T>
inline bool Line3<T>::isApproxEqual(const Line3<T>& line, T tolerance) const
{
	return origin.isApproxEqual(line.origin, tolerance) && direction.isApproxEqual(line.direction, tolerance);
}

template<typename T>
inline void Line3<T>::transform(const Matrix3<T>& matrix)
{
	origin *= matrix;
	direction *= matrix;
}

template<typename T>
inline void Line3<T>::transform(const AffineTransform<T>& transformation)
{
	origin.transform(transformation);
	direction *= transformation.getBasis();
}

template<typename T>
inline Vector3<T> Line3<T>::getClosestPoint(const Vector3<T>& point) const
{
	return dot(point - origin, direction)*direction + origin;
}

template<typename T>
template<typename U>
inline Vector3<T> Line3<T>::getClosestPoint(const Vector3<T>& point) const
{
	static_assert(std::is_same_v<U, Normalized> || std::is_same_v<U, Unnormalized>);
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

} // namespace core::mathematics

namespace std {

template<typename T>
struct hash;

template<typename T>
struct hash<::core::mathematics::templates::Line3<T>>
{
	std::size_t operator()(const ::core::mathematics::templates::Line3<T>& line) const noexcept
	{
		std::hash<T> hasher;
		std::size_t seed = hasher(line.origin) + 0x9e3779b9;
		seed ^= hasher(line.direction) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
		return seed;
	}
};

} // namespace std

#include "Ray3.hpp"

namespace core::mathematics::templates {

template<typename T>
inline Line3<T>::Line3(const Ray3<T>& ray) : origin(ray.origin), direction(ray.direction)
{
}

template<typename T>
inline const Ray3<T>& Line3<T>::asRay() const
{ 
	return reinterpret_cast<const Ray3<T>&>(*this);
}

} // namespace core::mathematics::templates
