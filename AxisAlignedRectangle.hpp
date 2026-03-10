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

namespace mathematics {
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
	explicit AxisAlignedRectangle(const std::pair<Vector2<T>, Vector2<T>>& t) noexcept : minimum(t.first), maximum(t.second) {}
	explicit AxisAlignedRectangle(const std::tuple<Vector2<T>, Vector2<T>>& t) noexcept : minimum(std::get<0>(t)), maximum(std::get<1>(t)) {}

	//explicit operator std::pair<Vector2<T>, Vector2<T>>() { return { minimum, maximum }; }
	//explicit operator std::tuple<Vector2<T>, Vector2<T>>() { return { minimum, maximum }; }
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
	template<std::output_iterator<Vector2<T>> O> O copyVertices(O target) const noexcept;
	std::vector<Vector2> getVertices() const;

	// Circumscribed sphere
	Circle2<T> getCircumscribedCircle() const noexcept;

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
	AxisAlignedRectangle& extendBy(const Circle2<T>& circle) noexcept;

	// Closest point
	Vector2<T> getClosestPoint(const Vector2<T>& point) const noexcept { return clamp(point, minimum, maximum); }
	T getDistanceTo(const Vector2<T>& point) const noexcept { return distance(point, getClosestPoint(point)); }

	// Containment and intersection
	bool contains(const Vector2<T>& point) const noexcept;
	bool contains(const AxisAlignedRectangle& rectangle) const noexcept;
	bool contains(const Circle2<T>& circle) const noexcept;
	bool intersects(const AxisAlignedRectangle& rectangle) const noexcept;
	bool intersects(const Circle2<T>& circle) const noexcept;

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
template<std::output_iterator<Vector2<T>> O> 
inline O AxisAlignedRectangle<T>::copyVertices(O target) const
{
	*target++ = minimum;
	(target++)->set(maximum.x, minimum.y);
	(target++)->set(minimum.x, maximum.y);
	*target++ = maximum;
	return target;
}

template<typename T>
inline std::vector<Vector2<T>> AxisAlignedRectangle<T>::getVertices() const
{
	std::vector<Vector2<T>> vertices;
	vertices.resize(4);
	vertices[0] = minimum;
	vertices[1].set(maximum.x, minimum.y);
	vertices[2].set(minimum.x, maximum.y);
	vertices[3] = maximum;
	return vertices;
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
inline bool AxisAlignedRectangle<T>::contains(const Vector2<T>& point) const
{
	return minimum.allLessThanEqual(point) && maximum.allGreaterThanEqual(point);
}

template<typename T>
inline bool AxisAlignedRectangle<T>::contains(const AxisAlignedRectangle<T>& rectangle) const
{
	return minimum.allLessThanEqual(rectangle.minimum) && maximum.allGreaterThanEqual(rectangle.maximum);
}

template<typename T>
inline bool AxisAlignedRectangle<T>::intersects(const AxisAlignedRectangle<T>& rectangle) const
{
	return minimum.allLessThanEqual(rectangle.maximum) && maximum.allGreaterThanEqual(rectangle.minimum);
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

} // namespace mathematics

namespace std {

template<typename T>
struct hash;

template<typename T>
struct hash<::mathematics::templates::AxisAlignedRectangle<T>>
{
	std::size_t operator()(const ::mathematics::templates::AxisAlignedRectangle<T>& rectangle) const noexcept
	{
		std::hash<T> hasher;
		std::size_t seed = hasher(rectangle.minimum) + 0x9e3779b9;
		seed ^= hasher(rectangle.maximum) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
		return seed;
	}
};

} // namespace std

#include "Circle2.hpp"

namespace mathematics::templates {
	
template<typename T>
inline Circle2<T> AxisAlignedRectangle<T>::getCircumscribedCircle() const
{
	Vector2<T> center = (minimum + maximum)*T(0.5);
	return Circle2<T>(center, distance(center, maximum));
}

template<typename T>
inline AxisAlignedRectangle<T>& AxisAlignedRectangle<T>::extendBy(const Circle2<T>& circle)
{
	Vector2<T> radius(circle.radius);
	setUnion(*this, AxisAlignedRectangle<T>(circle.center - radius, circle.center + radius));
	return *this;
}

template<typename T>
inline bool AxisAlignedRectangle<T>::contains(const Circle2<T>& circle) const
{
	Vector2<T> radius(circle.radius);
	return minimum.allLessThanEqual(circle.center - radius) && maximum.allGreaterThanEqual(circle.center + radius);
}

template<typename T>
inline bool AxisAlignedRectangle<T>::intersects(const Circle2<T>& circle) const
{
	return circle.intersects(*this);
}

} // namespace mathematics::templates
