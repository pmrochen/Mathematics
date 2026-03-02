/*
 *	Name: AxisAlignedRectangle
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
#include "Constants.hpp"
#include "Vector2.hpp"

namespace core::mathematics {
namespace templates {

template<typename T>
	requires std::floating_point<T>
struct Circle2;

template<typename T>
	requires std::floating_point<T>
struct AxisAlignedRectangle
{
	using Real = T;
	using ConstArg = const AxisAlignedRectangle&;
	using ConstResult = const AxisAlignedRectangle&;

	AxisAlignedRectangle() noexcept : minimum(Vector2<T>::INF), maximum(Vector2<T>::MINUS_INF) {}
	explicit AxisAlignedRectangle(Uninitialized) noexcept : minimum(Uninitialized()), maximum(Uninitialized()) {}
	AxisAlignedRectangle(const Vector2<T>& minimum, const Vector2<T>& maximum) noexcept : minimum(minimum), maximum(maximum) {}
	explicit AxisAlignedRectangle(const Vector2<T>& dimensions) noexcept : minimum(T(-0.5)*dimensions), maximum(T(0.5)*dimensions) {}
	explicit AxisAlignedRectangle(const std::pair<Vector3<T>, Vector3<T>>& t) noexcept : minimum(t.first), maximum(t.second) {}
	explicit AxisAlignedRectangle(const std::tuple<Vector3<T>, Vector3<T>>& t) noexcept : minimum(std::get<0>(t)), maximum(std::get<1>(t)) {}

	//explicit operator std::pair<Vector3<T>, Vector3<T>>() { return { minimum, maximum }; }
	//explicit operator std::tuple<Vector3<T>, Vector3<T>>() { return { minimum, maximum }; }
	bool operator==(const AxisAlignedRectangle& rectangle) const noexcept { return (minimum == rectangle.minimum) && (maximum == rectangle.maximum); }
	bool operator!=(const AxisAlignedRectangle& rectangle) const noexcept { return !(*this == rectangle); }

	template<typename A> void serialize(A& ar) { ar(minimum, maximum); }

	// Properties
	bool isEmpty() const noexcept { return minimum.anyGreaterThan(maximum); }
	bool approxEquals(const AxisAlignedRectangle& rectangle) const noexcept;
	bool approxEquals(const AxisAlignedRectangle& rectangle, T tolerance) const noexcept;
	bool isFinite() const noexcept { return minimum.isFinite() && maximum.isFinite(); }
	AxisAlignedRectangle& makeEmpty() noexcept { minimum = Vector2<T>::INF; maximum = Vector2<T>::MINUS_INF; return *this; }
	AxisAlignedRectangle& set(const Vector2<T>& minimum, const Vector2<T>& maximum) noexcept { this->minimum = minimum; this->maximum = maximum; return *this; }
	const Vector2<T>& getMinimum() const noexcept { return minimum; }
	void setMinimum(const Vector2<T>& minimum) noexcept { this->minimum = minimum; }
	const Vector2<T>& getMaximum() const noexcept { return maximum; }
	void setMaximum(const Vector2<T>& maximum) noexcept { this->maximum = maximum; }
	Vector2<T> getDimensions() const noexcept { return (maximum - minimum); }
	void setDimensions(const Vector2<T>& dimensions) noexcept;
	Vector2<T> getHalfDimensions() const noexcept { return (maximum - minimum)*T(0.5); }
	void setHalfDimensions(const Vector2<T>& halfDims) noexcept;
	Vector2<T> getCenter() const noexcept { return (minimum + maximum)*T(0.5); }
	void setCenter(const Vector2<T>& center) noexcept;
	T getDiagonal() const noexcept { return distance(minimum, maximum); }
	T getPerimeter() const noexcept;
	T getArea() const noexcept;

	// Vertices
	template<std::output_iterator<Vector2<T>> O, std::sentinel_for<O> S> O copyVertices(O first, S last) const noexcept;
	void copyVertices(std::vector<Vector2>& vertices) const;

	// Transformation
	AxisAlignedRectangle& inflate(const Vector2<T>& halfDims) noexcept { minimum -= halfDims; maximum += halfDims; return *this; }
	AxisAlignedRectangle& translate(const Vector2<T>& offset) noexcept { minimum += offset; maximum += offset; return *this; }

	// Union and intersection
	AxisAlignedRectangle& setUnion(const AxisAlignedRectangle& a, const AxisAlignedRectangle& b) noexcept;
	AxisAlignedRectangle& setIntersection(const AxisAlignedRectangle& a, const AxisAlignedRectangle& b) noexcept;
	static AxisAlignedRectangle makeUnion(const AxisAlignedRectangle& a, const AxisAlignedRectangle& b) noexcept { return AxisAlignedRectangle(Uninitialized()).setUnion(a, b); }
	static AxisAlignedRectangle makeIntersection(const AxisAlignedRectangle& a, const AxisAlignedRectangle& b) noexcept { return AxisAlignedRectangle(Uninitialized()).setIntersection(a, b); }
	AxisAlignedRectangle& extendBy(const Vector2<T>& point) noexcept;
	AxisAlignedRectangle& extendBy(const AxisAlignedRectangle& rectangle) noexcept { return setUnion(*this, rectangle); }

	// Closest point
	Vector2<T> getClosestPoint(const Vector2<T>& point) const noexcept;
	T getDistanceTo(const Vector2<T>& point) const noexcept { return distance(point, getClosestPoint(point)); }

	// Containment and intersection
	bool contains(const Vector2<T>& point) const noexcept;
	bool contains(const AxisAlignedRectangle& rectangle) const noexcept;
	bool contains(const Circle2<T>& circle) const noexcept;
	bool intersects(const AxisAlignedRectangle& rectangle) const noexcept;

	static const AxisAlignedRectangle EMPTY;

	Vector2<T> minimum;
	Vector2<T> maximum;
};

template<typename T> const AxisAlignedRectangle<T> AxisAlignedRectangle<T>::EMPTY{ Vector2<T>::INF, Vector2<T>::MINUS_INF };

template<typename C, typename T, typename U>
	requires std::floating_point<U>
inline std::basic_istream<C, T>& operator>>(std::basic_istream<C, T>& s, AxisAlignedRectangle<U>& rectangle)
{ 
	return s >> rectangle.minimum >> std::ws >> rectangle.maximum;
}

template<typename C, typename T, typename U>
	requires std::floating_point<U>
inline std::basic_ostream<C, T>& operator<<(std::basic_ostream<C, T>& s, const AxisAlignedRectangle<U>& rectangle)
{ 
	constexpr C WS(0x20);
	return s << rectangle.minimum << WS << rectangle.maximum;
}

template<typename T>
inline bool AxisAlignedRectangle<T>::approxEquals(const AxisAlignedRectangle<T>& rectangle) const
{
	return minimum.approxEquals(rectangle.minimum) && maximum.approxEquals(rectangle.maximum);
}

template<typename T>
inline bool AxisAlignedRectangle<T>::approxEquals(const AxisAlignedRectangle<T>& rectangle, T tolerance) const
{
	return minimum.approxEquals(rectangle.minimum, tolerance) && maximum.approxEquals(rectangle.maximum, tolerance);
}

template<typename T>
inline void AxisAlignedRectangle<T>::setDimensions(const Vector2<T>& dimensions)
{
	Vector2<T> center = (minimum + maximum)*T(0.5);
	Vector2<T> halfDims = T(0.5)*dimensions;
	minimum = center - halfDims;
	maximum = center + halfDims;
}

template<typename T>
inline void AxisAlignedRectangle<T>::setHalfDimensions(const Vector2<T>& halfDims)
{
	Vector2<T> center = (minimum + maximum)*T(0.5);
	minimum = center - halfDims;
	maximum = center + halfDims;
}

template<typename T>
inline void AxisAlignedRectangle<T>::setCenter(const Vector2<T>& center)
{
	Vector2<T> diff = center - (minimum + maximum)*T(0.5);
	minimum += diff;
	maximum += diff;
}

template<typename T>
inline T AxisAlignedRectangle<T>::getPerimeter() const
{
	Vector2<T> dim = maximum - minimum;
	return T(2)*(dim.x + dim.y);
}

template<typename T>
inline T AxisAlignedRectangle<T>::getArea() const
{
	Vector2<T> dim = maximum - minimum;
	return dim.x*dim.y;
}

template<typename T>
template<std::output_iterator<Vector2<T>> O, std::sentinel_for<O> S> 
inline O AxisAlignedRectangle<T>::copyVertices(O first, S last) const
{
	*first++ = minimum;
	(first++)->set(maximum.x, minimum.y);
	(first++)->set(minimum.x, maximum.y);
	*first++ = maximum;
	return first;
}

template<typename T>
inline void AxisAlignedRectangle<T>::getVertices(std::vector<Vector2>& vertices) const
{
	vertices.resize(4);
	vertices[0] = minimum;
	vertices[1].set(maximum.x, minimum.y);
	vertices[2].set(minimum.x, maximum.y);
	vertices[3] = maximum;
}

template<typename T>
inline AxisAlignedRectangle<T>& AxisAlignedRectangle<T>::setUnion(const AxisAlignedRectangle<T>& a, const AxisAlignedRectangle<T>& b)
{
	minimum.setMinimum(a.minimum, b.minimum);
	maximum.setMaximum(a.maximum, b.maximum);
	return *this;
}

template<typename T>
inline AxisAlignedRectangle<T>& AxisAlignedRectangle<T>::setIntersection(const AxisAlignedRectangle<T>& a, const AxisAlignedRectangle<T>& b)
{
	minimum.setMaximum(a.minimum, b.minimum);
	maximum.setMinimum(a.maximum, b.maximum);
	return *this;
}

template<typename T>
inline AxisAlignedRectangle<T>& AxisAlignedRectangle<T>::extendBy(const Vector2<T>& point)
{
	minimum.setMinimum(minimum, point);
	maximum.setMaximum(maximum, point);
	return *this;
}

template<typename T>
inline Vector2<T> AxisAlignedRectangle<T>::getClosestPoint(const Vector2<T>& point) const
{
	return minimum(maximum(point, this->minimum), this->maximum);
}

template<typename T>
inline bool AxisAlignedRectangle<T>::contains(const Vector2<T>& point) const
{
	return (minimum.x <= point.x) && (maximum.x >= point.x) && // #TODO SIMD
		(minimum.y <= point.y) && (maximum.y >= point.y);
}

template<typename T>
inline bool AxisAlignedRectangle<T>::contains(const AxisAlignedRectangle<T>& rectangle) const
{
	return (minimum.x <= rectangle.minimum.x) && (maximum.x >= rectangle.maximum.x) && // #TODO SIMD
		(minimum.y <= rectangle.minimum.y) && (maximum.y >= rectangle.maximum.y);
}

template<typename T>
inline bool AxisAlignedRectangle<T>::intersects(const AxisAlignedRectangle<T>& rectangle) const
{
	return (minimum.x <= rectangle.maximum.x) && (maximum.x >= rectangle.minimum.x) && // #TODO SIMD
		(minimum.y <= rectangle.maximum.y) && (maximum.y >= rectangle.minimum.y);
}

} // namespace templates

#if MATHEMATICS_DOUBLE
using AxisAlignedRectangle = templates::AxisAlignedRectangle<double>;
using AxisAlignedRectangleArg = templates::AxisAlignedRectangle<double>::ConstArg;
using AxisAlignedRectangleResult = templates::AxisAlignedRectangle<double>::ConstResult;
#else
using AxisAlignedRectangle = templates::AxisAlignedRectangle<float>;
using AxisAlignedRectangleArg = templates::AxisAlignedRectangle<float>::ConstArg;
using AxisAlignedRectangleResult = templates::AxisAlignedRectangle<float>::ConstResult;
#maximumif

} // namespace core::mathematics

namespace std {

template<typename T>
struct hash;

template<typename T>
struct hash<::core::mathematics::templates::AxisAlignedRectangle<T>>
{
	std::size_t operator()(const ::core::mathematics::templates::AxisAlignedRectangle<T>& rectangle) const noexcept
	{
		std::hash<T> hasher;
		std::size_t seed = hasher(rectangle.minimum) + 0x9e3779b9;
		seed ^= hasher(rectangle.maximum) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
		return seed;
	}
};

} // namespace std

#include "Circle2.hpp"

namespace core::mathematics::templates {

template<typename T>
inline bool AxisAlignedRectangle<T>::contains(const Circle2<T>& circle) const
{
	return (minimum.x <= (circle.center.x - circle.radius)) && (maximum.x >= (circle.center.x + circle.radius)) && // #TODO SIMD
		(minimum.y <= (circle.center.y - circle.radius)) && (maximum.y >= (circle.center.y + circle.radius));
}

} // namespace core::mathematics::templates
