/*
 *	Name: LineSegment3
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
#include <vector>
#include <iterator>
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

template<typename T, typename U>
	requires (std::floating_point<T> && std::integral<U>)
class TriangleMesh;

#if SIMD_HAS_FLOAT4

template<>
struct HalfSpace<float>;

template<>
struct Plane<float>;

#endif /* SIMD_HAS_FLOAT4 */

template<typename T>
	requires std::floating_point<T>
struct LineSegment3
{
	using Real = T;
	using ConstArg = const LineSegment3&;
	using ConstResult = const LineSegment3&;
	using EndpointType = Vector3<T>;
	using PairType = std::pair<Vector3<T>, Vector3<T>>;
	using TupleType = std::tuple<Vector3<T>, Vector3<T>>;

	LineSegment3() = default;
	explicit LineSegment3(Uninitialized) noexcept : start(Uninitialized()), end(Uninitialized()) {}
	LineSegment3(const Vector3<T>& start, const Vector3<T>& end) noexcept : start(start), end(end) {}
	explicit LineSegment3(const PairType& t) noexcept : start(t.first), end(t.second) {}
	explicit LineSegment3(const TupleType& t) noexcept : start(std::get<0>(t)), end(std::get<1>(t)) {}
	explicit LineSegment3(const Line3<T>& line) noexcept : start(line.origin), end(line.origin + line.direction) {}
	LineSegment3(const Line3<T>& line, const Interval<T>& interval) noexcept : start(line.evaluate(interval.minimum)), end(line.evaluate(interval.maximum)) {}
	explicit LineSegment3(const Ray3<T>& ray) noexcept : start(ray.origin), end(ray.origin + ray.direction) {}
	LineSegment3(const Ray3<T>& ray, const Interval<T>& interval) noexcept : start(ray.evaluate(interval.minimum)), end(ray.evaluate(interval.maximum)) {}

	//Vector3 operator()(T t) const noexcept { return lerp(start, end, t); }
	bool operator==(const LineSegment3& segment) const noexcept { return (start == segment.start) && (end == segment.end); }
	bool operator!=(const LineSegment3& segment) const noexcept { return !(*this == segment); }

	template<typename A> void serialize(A& ar) { ar(start, end); }

	// Properties
	bool isZero() const noexcept { return start.isZero() && end.isZero(); }
	bool isApproxZero() const noexcept { return start.isApproxZero() && end.isApproxZero(); }
	bool approxEquals(const LineSegment3& segment) const noexcept;
	bool approxEquals(const LineSegment3& segment, T tolerance) const noexcept;
	bool isFinite() const noexcept { return start.isFinite() && end.isFinite(); }
	LineSegment3& setZero() noexcept { start.setZero(); end.setZero(); return *this; }
	LineSegment3& set(const Vector3<T>& start, const Vector3<T>& end) noexcept { this->start = start; this->end = end; return *this; }
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
	LineSegment3& translate(const Vector3<T>& offset) noexcept { start += offset; end += offset; return *this; }
	LineSegment3& transform(const Matrix3<T>& matrix) noexcept;
	LineSegment3& transform(const AffineTransform<T>& transformation) noexcept;

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
	template<std::integral U> bool intersects(const TriangleMesh<T, U>* mesh) const noexcept;
	std::optional<T> findIntersection(const Plane<T>& plane) const noexcept;
	std::optional<T> findIntersection(const Triangle3<T>& triangle) const noexcept;
	//template<ScalarOrVector3<T> U> std::optional<U> findIntersection(const Plane<T>& plane) const noexcept;
	std::optional<Interval<T>> findIntersection(const AxisAlignedBox<T>& box) const noexcept;
	std::optional<Interval<T>> findIntersection(const OrientedBox<T>& box) const noexcept;
	std::optional<Interval<T>> findIntersection(const Sphere<T>& sphere) const noexcept;
	std::optional<Interval<T>> findIntersection(const Ellipsoid<T>& ellipsoid) const noexcept;
	template<std::integral U> std::vector<T> findIntersections(const TriangleMesh<T, U>* mesh) const noexcept;

	Vector3<T> start;
	Vector3<T> end;
};

template<typename C, typename T, typename U>
	requires std::floating_point<U>
inline std::basic_istream<C, T>& operator>>(std::basic_istream<C, T>& s, LineSegment3<U>& segment)
{ 
	return s >> segment.start >> std::ws >> segment.end;
}

template<typename C, typename T, typename U>
	requires std::floating_point<U>
inline std::basic_ostream<C, T>& operator<<(std::basic_ostream<C, T>& s, const LineSegment3<U>& segment)
{ 
	constexpr C WS(0x20);
	return s << segment.start << WS << segment.end;
}

template<typename T>
inline bool LineSegment3<T>::approxEquals(const LineSegment3<T>& segment) const
{
	return start.approxEquals(segment.start) && end.approxEquals(segment.end);
}

template<typename T>
inline bool LineSegment3<T>::approxEquals(const LineSegment3<T>& segment, T tolerance) const
{
	return start.approxEquals(segment.start, tolerance) && end.approxEquals(segment.end, tolerance);
}

template<typename T>
template<std::output_iterator<Vector2<T>> O>
inline O LineSegment3<T>::copyEndpoints(O target) const
{
	*target++ = start;
	*target++ = end;
	return target;
}

template<typename T>
inline LineSegment3<T>& LineSegment3<T>::transform(const Matrix3<T>& matrix)
{
	start *= matrix;
	end *= matrix;
	return *this;
}

template<typename T>
inline LineSegment3<T>& LineSegment3<T>::transform(const AffineTransform<T>& transformation)
{
	start.transform(transformation);
	end.transform(transformation);
	return *this;
}

template<typename T>
inline Vector3<T> LineSegment3<T>::getClosestPoint(const Vector3<T>& point) const
{
	Vector3<T> direction = end - start;
	return std::clamp(dot(point - start, direction)/dot(direction, direction), T(0), T(1))*direction + start;
}

} // namespace templates

#if MATHEMATICS_DOUBLE
using LineSegment3 = templates::LineSegment3<double>;
using LineSegment3Arg = templates::LineSegment3<double>::ConstArg;
using LineSegment3Result = templates::LineSegment3<double>::ConstResult;
#else
using LineSegment3 = templates::LineSegment3<float>;
using LineSegment3Arg = templates::LineSegment3<float>::ConstArg;
using LineSegment3Result = templates::LineSegment3<float>::ConstResult;
#endif

} // namespace mathematics

namespace std {

template<typename T>
struct hash;

template<typename T>
struct hash<::mathematics::templates::LineSegment3<T>>
{
	size_t operator()(const ::mathematics::templates::LineSegment3<T>& segment) const noexcept
	{
		hash<typename ::mathematics::templates::Vector3<T>> hasher;
		size_t seed = hasher(segment.start) + 0x9e3779b9;
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
#include "TriangleMesh.hpp"
#include "Intersections.inl"

namespace mathematics::templates {

template<typename T>
inline bool LineSegment3<T>::intersects(const HalfSpace<T>& halfSpace) const
{
	return halfSpace.contains(start) || halfSpace.contains(end);
}

template<typename T>
inline bool LineSegment3<T>::intersects(const Ellipsoid<T>& ellipsoid) const
{
	return intersections::testLineSegmentEllipsoid(start, end, ellipsoid.center, ellipsoid.getMatrix());
}

template<typename T>
template<std::integral U> 
bool LineSegment3<T>::intersects(const TriangleMesh<T, U>* mesh) const
{
	if (!mesh)
		return false;

	if (mesh->indices.empty())
	{
		if (mesh->vertices.size() < 3)
			return false;

		Vector3<T> direction = end - start;
		const Vector3<T>* positions = &mesh->vertices[0];
		std::size_t count = mesh->vertices.size()/3;
		do
		{
			auto result = intersections::findLineTriangle<std::optional<T>>(start, direction, 
				positions[0], positions[1], positions[2]);
			if (result.has_value() && (result.value() >= T(0)) && (result.value() <= T(1)))
				return true;
			positions += 3;
		} while (--count);
	}
	else
	{
		if (mesh->indices.size() < 3)
			return false;

		Vector3<T> direction = end - start;
		const Vector3<T>* positions = &mesh->vertices[0];
		const U* indices = &mesh->indices[0];
		std::size_t count = mesh->indices.size()/3;
		do
		{
			auto result = intersections::findLineTriangle<std::optional<T>>(start, direction, 
				positions[indices[0]], positions[indices[1]], positions[indices[2]]);
			if (result.has_value() && (result.value() >= T(0)) && (result.value() <= T(1)))
				return true;
			indices += 3;
		} while (--count);
	}

	return false;
}

template<typename T>
inline std::optional<T> LineSegment3<T>::findIntersection(const Plane<T>& plane) const
{
	auto result = intersections::findLinePlane<std::optional<T>>(start, end - start, plane.getNormal(), plane.d);
	return (result.has_value() && (result.value() >= T(0)) && (result.value() <= T(1))) ? result : {};
}

//template<typename T>
//template<ScalarOrVector3<T> U>
//inline std::optional<U> LineSegment3<T>::findIntersection(const Plane<T>& plane) const
//{
//	std::optional<T> result = findIntersection(plane);
//	if constexpr (std::is_same_v<U, T>)
//		return result;
//	else //if constexpr (std::is_same_v<U, Vector3<T>>)
//		return result.has_value() ? { evaluate(result.value()) } : {};
//}

template<typename T>
inline std::optional<T> LineSegment3<T>::findIntersection(const Triangle3<T>& triangle) const
{
	//return intersections::findLineSegmentTriangle<std::optional<T>>(start, end, triangle.vertices[0], triangle.vertices[1],
	//	triangle.vertices[2]);
	auto result = intersections::findLineTriangle<std::optional<T>>(start, end - start, triangle.vertices[0], triangle.vertices[1],
		triangle.vertices[2]);
	return (result.has_value() && (result.value() >= T(0)) && (result.value() <= T(1))) ? result : {};
}

template<typename T>
inline std::optional<Interval<T>> LineSegment3<T>::findIntersection(const AxisAlignedBox<T>& box) const
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
inline std::optional<Interval<T>> LineSegment3<T>::findIntersection(const OrientedBox<T>& box) const
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
inline std::optional<Interval<T>> LineSegment3<T>::findIntersection(const Sphere<T>& sphere) const
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
inline std::optional<Interval<T>> LineSegment3<T>::findIntersection(const Ellipsoid<T>& ellipsoid) const
{
	return intersections::findLineSegmentEllipsoid<std::optional<Interval<T>>>(start, end, ellipsoid.center, ellipsoid.getMatrix());
}

template<typename T>
template<std::integral U> 
std::vector<T> LineSegment3<T>::findIntersections(const TriangleMesh<T, U>* mesh) const
{
	std::vector<T> intersections;
	if (!mesh)
		return intersections;

	if (mesh->indices.empty())
	{
		if (mesh->vertices.size() < 3)
			return intersections;

		Vector3<T> direction = end - start;
		const Vector3<T>* positions = &mesh->vertices[0];
		std::size_t count = mesh->vertices.size()/3;
		do
		{
			auto result = intersections::findLineTriangle<std::optional<T>>(start, direction, 
				positions[0], positions[1], positions[2]);
			if (result.has_value() && (result.value() >= T(0)) && (result.value() <= T(1)))
				intersections.push_back(result.value());
			positions += 3;
		} while (--count);
	}
	else
	{
		if (mesh->indices.size() < 3)
			return intersections;

		Vector3<T> direction = end - start;
		const Vector3<T>* positions = &mesh->vertices[0];
		const U* indices = &mesh->indices[0];
		std::size_t count = mesh->indices.size()/3;
		do
		{
			auto result = intersections::findLineTriangle<std::optional<T>>(start, direction, 
				positions[indices[0]], positions[indices[1]], positions[indices[2]]);
			if (result.has_value() && (result.value() >= T(0)) && (result.value() <= T(1)))
				intersections.push_back(result.value());
			indices += 3;
		} while (--count);
	}

	return intersections;
}

} // namespace mathematics::templates
