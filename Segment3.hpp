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
#include <algorithm>
#include <cstddef>
#include <cmath>
#include "Vector3.hpp"
#include "Matrix3.hpp"
#include "AffineTransform.hpp"
#include "Line3.hpp"
#include "Ray3.hpp"

namespace core::mathematics {
namespace templates {

template<typename T>
	requires std::floating_point<T>
struct Segment3
{
	using Real = T;
	using ConstArg = const Segment3&;
	using ConstResult = const Segment3&;

	Segment3() = default;
	Segment3(const Vector3<T>& start, const Vector3<T>& end) noexcept : start(start), end(end) {}
	explicit Segment3(const std::pair<Vector3<T>, Vector3<T>>& t) noexcept : start(t.first), end(t.second) {}
	explicit Segment3(const std::tuple<Vector3<T>, Vector3<T>>& t) noexcept : start(std::get<0>(t)), end(std::get<1>(t)) {}

	//explicit operator std::pair<Vector3<T>, Vector3<T>>() { return { start, end }; }
	//explicit operator std::tuple<Vector3<T>, Vector3<T>>() { return { start, end }; }
	//Vector3 operator()(T t) const noexcept { return lerp(start, end, t); }
	bool operator==(const Segment3& segment) const noexcept { return (start == segment.start) && (end == segment.end); }
	bool operator!=(const Segment3& segment) const noexcept { return !(*this == segment); }

	template<typename A> void serialize(A& ar) { ar(start, end); }

	// Properties
	bool isApproxEqual(const Segment3& segment) const noexcept;
	bool isApproxEqual(const Segment3& segment, T tolerance) const noexcept;
	bool isFinite() const noexcept { return start.isFinite() && end.isFinite(); }
	void set(const Vector3<T>& start, const Vector3<T>& end) noexcept { this->start = start; this->end = end; }
	const Vector3<T>& getStart() const noexcept { return start; }
	void setStart(const Vector3<T>& start) noexcept { this->start = start; }
	const Vector3<T>& getEnd() const noexcept { return end; }
	void setEnd(const Vector3<T>& end) noexcept { this->end = end; }
	Vector3<T> getDirection() const noexcept { return normalize(end - start); }
	T getLength() const noexcept { return distance(start, end); }
	Vector3<T> getCenter() const noexcept { return lerp(start, end, T(0.5)); }

	// Transformation
	void translate(const Vector3<T>& offset) noexcept { start += offset; end += offset; }
	void transform(const Matrix3<T>& matrix) noexcept;
	void transform(const AffineTransform<T>& transformation) noexcept;

	// Evaluation
	Vector3<T> evaluate(T t) const noexcept { return lerp(start, end, t); }

	// Closest points
	Vector3<T> getClosestPoint(const Vector3<T>& point) const;
	T getDistance(const Vector3<T>& point) const { return distance(getClosestPoint(point), point); }

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
inline bool Segment3<T>::isApproxEqual(const Segment3<T>& segment) const
{
	return origin.isApproxEqual(segment.start) && direction.isApproxEqual(segment.end);
}

template<typename T>
inline bool Segment3<T>::isApproxEqual(const Segment3<T>& segment, T tolerance) const
{
	return origin.isApproxEqual(segment.start, tolerance) && direction.isApproxEqual(segment.end, tolerance);
}

template<typename T>
inline void Segment3<T>::transform(const Matrix3<T>& matrix)
{
	start *= matrix;
	end *= matrix;
}

template<typename T>
inline void Segment3<T>::transform(const AffineTransform<T>& transformation)
{
	start.transform(transformation);
	end.transform(transformation);
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

} // namespace core::mathematics

namespace std {

template<typename T>
struct hash;

template<typename T>
struct hash<::core::mathematics::templates::Segment3<T>>
{
	std::size_t operator()(const ::core::mathematics::templates::Segment3<T>& segment) const noexcept
	{
		std::hash<T> hasher;
		std::size_t seed = hasher(segment.start) + 0x9e3779b9;
		seed ^= hasher(segment.end) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
		return seed;
	}
};

} // namespace std
