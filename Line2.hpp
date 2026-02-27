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
#include <algorithm>
#include <cstddef>
#include <cmath>
#include "Constants.hpp"
#include "Vector2.hpp"

namespace core::mathematics {
namespace templates {

template<typename T, typename U>
concept ScalarOrVector2 = (std::same_as<T, U> || std::same_as<T, Vector2<U>>); // #TODO Move to Concepts.hpp

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

	// Normalization
	void normalize() noexcept { direction.normalize(); }

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
	bool intersects(const Segment2<T>& segment) const noexcept { return findIntersection(segment).has_value(); }
	std::optional<T> findIntersection(const Line2& line) const noexcept;
	std::optional<T> findIntersection(const Segment2<T>& segment) const noexcept;
	template<ScalarOrVector2<T> U> std::optional<U> findIntersection(const Line2& line) const noexcept;
	template<ScalarOrVector2<T> U> std::optional<U> findIntersection(const Segment2<T>& segment) const noexcept;

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
inline std::optional<T> Line2<T>::findIntersection(const Line2& line) const
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
	return std::optional<T>(t);
}

template<typename T>
inline std::optional<T> Line2<T>::findIntersection(const Segment2<T>& segment) const
{
	T d1CrossD2 = cross(direction, segment.end - segment.start);
	if (std::fabs(d1CrossD2) < Constants<T>::TOLERANCE)
	{
		if (!(std::fabs(cross(normalize(segment.start - origin), direction)) < Constants<T>::TOLERANCE))
			return {};
		return std::optional<T>(T(0)); // collinear
	}

	//T t = cross(segment.start - origin, segment.end - segment.start)/d1CrossD2;
	T u = cross(segment.start - origin, direction)/d1CrossD2;
	if ((u >= T(0)) && (u <= T(1))
		return std::optional<T>(/*t*/cross(segment.start - origin, segment.end - segment.start)/d1CrossD2);
	return {};
}

template<typename T>
template<ScalarOrVector2<T> U>
inline std::optional<U> Line2<T>::findIntersection(const Line2& line) const
{
	std::optional<T> result = findIntersection(line);
	if constexpr (std::is_same_v<U, T>)
		return result;
	else //if constexpr (std::is_same_v<U, Vector3<T>>)
		return result.has_value() ? std::optional<U>(line.evaluate(result.value())) : std::optional<U>();
}

template<typename T>
template<ScalarOrVector2<T> U>
inline std::optional<U> Line2<T>::findIntersection(const Segment2<T>& segment) const
{
	std::optional<T> result = findIntersection(segment);
	if constexpr (std::is_same_v<U, T>)
		return result;
	else //if constexpr (std::is_same_v<U, Vector3<T>>)
		return result.has_value() ? std::optional<U>(segment.evaluate(result.value())) : std::optional<U>();
}

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

} // namespace core::mathematics

namespace std {

template<typename T>
struct hash;

template<typename T>
struct hash<::core::mathematics::templates::Line2<T>>
{
	std::size_t operator()(const ::core::mathematics::templates::Line2<T>& line) const noexcept
	{
		std::hash<T> hasher;
		std::size_t seed = hasher(line.origin) + 0x9e3779b9;
		seed ^= hasher(line.direction) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
		return seed;
	}
};

} // namespace std

#include "Ray2.hpp"

namespace core::mathematics::templates {

template<typename T>
inline Line2<T>::Line2(const Ray2<T>& ray) : origin(ray.origin), direction(ray.direction)
{
}

template<typename T>
inline const Ray2<T>& Line2<T>::asRay() const
{ 
	return reinterpret_cast<const Ray2<T>&>(*this);
}

} // namespace core::mathematics::templates
