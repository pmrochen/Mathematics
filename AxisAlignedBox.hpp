/*
 *	Name: AxisAlignedBox
 *	Author: Pawel Mrochen
 */

#pragma once

#include <istream>
#include <ostream>
#include <type_traits>
#include <concepts>
#include <functional>
#include <utility>
#include <tuple>
#include <vector>
#include <iterator>
#include <algorithm>
#include <cstddef>
#include <cmath>
#include "Constants.hpp"
#include "Vector3.hpp"
#include "Matrix3.hpp"
#include "AffineTransform.hpp"
#include "HalfSpace.hpp"
#include "Plane.hpp"
#include "Triangle3.hpp"

namespace mathematics {
namespace templates {

template<typename T>
	requires std::floating_point<T>
struct OrientedBox;

template<typename T>
	requires std::floating_point<T>
struct Sphere;

template<typename T>
	requires std::floating_point<T>
struct SymmetricFrustum;

template<typename T>
	requires std::floating_point<T>
struct AxisAlignedBox
{
	using Real = T;
	using ConstArg = const AxisAlignedBox&;
	using ConstResult = const AxisAlignedBox&;
	using VertexType = Vector3<T>;
	using PairType = std::pair<Vector3<T>, Vector3<T>>;
	using TupleType = std::tuple<Vector3<T>, Vector3<T>>;

	AxisAlignedBox() noexcept : minimum(Vector3<T>::INF), maximum(Vector3<T>::MINUS_INF) {}
	explicit AxisAlignedBox(Uninitialized) noexcept : minimum(Uninitialized()), maximum(Uninitialized()) {}
	AxisAlignedBox(const Vector3<T>& minimum, const Vector3<T>& maximum) noexcept : minimum(minimum), maximum(maximum) {}
	explicit AxisAlignedBox(const Vector3<T>& dimensions) noexcept : minimum(T(-0.5)* dimensions), maximum(T(0.5)* dimensions) {}
	explicit AxisAlignedBox(const PairType& t) noexcept : minimum(t.first), maximum(t.second) {}
	explicit AxisAlignedBox(const TupleType& t) noexcept : minimum(std::get<0>(t)), maximum(std::get<1>(t)) {}
	explicit AxisAlignedBox(const OrientedBox<T>& box) noexcept;

	bool operator==(const AxisAlignedBox& box) const noexcept { return (minimum == box.minimum) && (maximum == box.maximum); }
	bool operator!=(const AxisAlignedBox& box) const noexcept { return !(*this == box); }

	template<typename A> void serialize(A& ar) { ar(minimum, maximum); }

	// Properties
	bool isEmpty() const noexcept { return minimum.anyGreaterThan(maximum); }
	bool approxEquals(const AxisAlignedBox& box) const noexcept;
	bool approxEquals(const AxisAlignedBox& box, T tolerance) const noexcept;
	bool isFinite() const noexcept { return minimum.isFinite() && maximum.isFinite(); }
	AxisAlignedBox& makeEmpty() noexcept { minimum = Vector3<T>::INF; maximum = Vector3<T>::MINUS_INF; return *this; }
	AxisAlignedBox& set(const Vector3<T>& minimum, const Vector3<T>& maximum) noexcept { this->minimum = minimum; this->maximum = maximum; return *this; }
	const Vector3<T>& getMinimum() const noexcept { return minimum; }
	void setMinimum(const Vector3<T>& minimum) noexcept { this->minimum = minimum; }
	const Vector3<T>& getMaximum() const noexcept { return maximum; }
	void setMaximum(const Vector3<T>& maximum) noexcept { this->maximum = maximum; }
	Vector3<T> getDimensions() const noexcept { return (maximum - minimum); }
	void setDimensions(const Vector3<T>& dimensions) noexcept;
	Vector3<T> getHalfDimensions() const noexcept { return (maximum - minimum)*T(0.5); }
	void setHalfDimensions(const Vector3<T>& halfDims) noexcept;
	Vector3<T> getCenter() const noexcept { return (minimum + maximum)*T(0.5); }
	void setCenter(const Vector3<T>& center) noexcept;
	T getDiagonal() const noexcept { return distance(minimum, maximum); }
	T getSurfaceArea() const noexcept;
	T getVolume() const noexcept;

	// Vertices
	template<std::output_iterator<Vector3<T>> O> O copyVertices(O target) const;
	std::vector<Vector3<T>> getVertices() const;

	// Primitives
	template<std::integral U> std::pair<const U*, const U*> getPrimitives(int nVerticesPerPrimitive) const noexcept; // #TODO return range
	std::size_t getPrimitiveCount(int nVerticesPerPrimitive) const noexcept;

	// Half spaces
	template<std::output_iterator<HalfSpace<T>> O> O copyHalfSpaces(O target) const;
	std::vector<HalfSpace<T>> getHalfSpaces() const;

	// Circumscribed sphere
	Sphere<T> getCircumscribedSphere() const noexcept;

	// Transformation
	AxisAlignedBox& inflate(const Vector3<T>& halfDims) noexcept { minimum -= halfDims; maximum += halfDims; return *this; }
	AxisAlignedBox& translate(const Vector3<T>& offset) noexcept { minimum += offset; maximum += offset; return *this; }
	AxisAlignedBox& transform(const Matrix3<T>& matrix, bool orthogonal = false) noexcept;
	AxisAlignedBox& transform(const AffineTransform<T>& transformation, bool orthogonal = false) noexcept;
	AxisAlignedBox& scaleAroundCenter(T factor) noexcept;

	// Union and intersection
	AxisAlignedBox& setUnion(const AxisAlignedBox& a, const AxisAlignedBox& b) noexcept;
	AxisAlignedBox& setIntersection(const AxisAlignedBox& a, const AxisAlignedBox& b) noexcept;
	static AxisAlignedBox makeUnion(const AxisAlignedBox& a, const AxisAlignedBox& b) noexcept { return AxisAlignedBox(Uninitialized()).setUnion(a, b); }
	static AxisAlignedBox makeIntersection(const AxisAlignedBox& a, const AxisAlignedBox& b) noexcept { return AxisAlignedBox(Uninitialized()).setIntersection(a, b); }
	AxisAlignedBox& extendBy(const Vector3<T>& point) noexcept;
	AxisAlignedBox& extendBy(const AxisAlignedBox& box) noexcept { return setUnion(*this, box); }
	AxisAlignedBox& extendBy(const Sphere<T>& sphere) noexcept;

	// Closest point
	Vector3<T> getClosestPoint(const Vector3<T>& point) const noexcept { return clamp(point, minimum, maximum); }
	T getDistanceTo(const Vector3<T>& point) const noexcept { return distance(point, getClosestPoint(point)); }

	// Containment and intersection
	bool contains(const Vector3<T>& point) const noexcept;
	bool contains(const AxisAlignedBox& box) const noexcept;
	bool contains(const Sphere<T>& sphere) const noexcept;
	int classify(const HalfSpace<T>& halfSpace) const noexcept; // -1 = inside, 1 = outside, 0 = partial // #TODO return enum
	bool intersects(const HalfSpace<T>& halfSpace) const noexcept;
	bool intersects(const Plane<T>& plane) const noexcept;
	bool intersects(const Triangle3<T>& triangle) const noexcept;
	bool intersects(const AxisAlignedBox& box) const noexcept;
	bool intersects(const OrientedBox<T>& box) const noexcept;
	bool intersects(const Sphere<T>& sphere) const noexcept;
	bool intersects(const SymmetricFrustum<T>& frustum) const noexcept;

	static const AxisAlignedBox EMPTY;

	Vector3<T> minimum;
	Vector3<T> maximum;
};

template<typename T> const AxisAlignedBox<T> AxisAlignedBox<T>::EMPTY{ Vector3<T>::INF, Vector3<T>::MINUS_INF };

template<typename C, typename T, typename U>
	requires std::floating_point<U>
inline std::basic_istream<C, T>& operator>>(std::basic_istream<C, T>& s, AxisAlignedBox<U>& box)
{ 
	return s >> box.minimum >> std::ws >> box.maximum;
}

template<typename C, typename T, typename U>
	requires std::floating_point<U>
inline std::basic_ostream<C, T>& operator<<(std::basic_ostream<C, T>& s, const AxisAlignedBox<U>& box)
{ 
	constexpr C WS(0x20);
	return s << box.minimum << WS << box.maximum;
}

template<typename T>
inline bool AxisAlignedBox<T>::approxEquals(const AxisAlignedBox<T>& box) const
{
	return minimum.approxEquals(box.minimum) && maximum.approxEquals(box.maximum);
}

template<typename T>
inline bool AxisAlignedBox<T>::approxEquals(const AxisAlignedBox<T>& box, T tolerance) const
{
	return minimum.approxEquals(box.minimum, tolerance) && maximum.approxEquals(box.maximum, tolerance);
}

template<typename T>
inline void AxisAlignedBox<T>::setDimensions(const Vector3<T>& dimensions)
{
	Vector3<T> center = (minimum + maximum)*T(0.5);
	Vector3<T> halfDims = T(0.5)*dimensions;
	minimum = center - halfDims;
	maximum = center + halfDims;
}

template<typename T>
inline void AxisAlignedBox<T>::setHalfDimensions(const Vector3<T>& halfDims)
{
	Vector3<T> center = (minimum + maximum)*T(0.5);
	minimum = center - halfDims;
	maximum = center + halfDims;
}

template<typename T>
inline void AxisAlignedBox<T>::setCenter(const Vector3<T>& center)
{
	Vector3<T> diff = center - (minimum + maximum)*T(0.5);
	minimum += diff;
	maximum += diff;
}

template<typename T>
inline T AxisAlignedBox<T>::getSurfaceArea() const
{
	Vector3<T> dim = maximum - minimum;
	return T(2)*(dim.x*dim.y + dim.y*dim.z + dim.z*dim.x); // T(2)*dot(dim, dim.yzx())
}

template<typename T>
inline T AxisAlignedBox<T>::getVolume() const
{
	Vector3<T> dim = maximum - minimum;
	return dim.x*dim.y*dim.z;
}

template<typename T>
template<std::output_iterator<Vector3<T>> O> 
inline O AxisAlignedBox<T>::copyVertices(O target) const
{
	*target++ = minimum;
	(target++)->set(maximum.x, minimum.y, minimum.z);
	(target++)->set(minimum.x, maximum.y, minimum.z);
	(target++)->set(maximum.x, maximum.y, minimum.z);
	(target++)->set(minimum.x, minimum.y, maximum.z);
	(target++)->set(maximum.x, minimum.y, maximum.z);
	(target++)->set(minimum.x, maximum.y, maximum.z);
	*target++ = maximum;
	return target;
}

template<typename T>
inline std::vector<Vector3<T>> AxisAlignedBox<T>::getVertices() const
{
	return { minimum, Vector3<T>(maximum.x, minimum.y, minimum.z), Vector3<T>(minimum.x, maximum.y, minimum.z),
		Vector3<T>(maximum.x, maximum.y, minimum.z), Vector3<T>(minimum.x, minimum.y, maximum.z),
		Vector3<T>(maximum.x, minimum.y, maximum.z), Vector3<T>(minimum.x, maximum.y, maximum.z), maximum };
}

template<typename T>
template<std::integral U> 
std::pair<const U*, const U*> AxisAlignedBox<T>::getPrimitives(int nVerticesPerPrimitive) const
{
	static const U edges[24] = { 0, 2, 2, 3, 3, 1, 1, 0, 3, 7, 5, 1, 6, 2, 0, 4, 5, 7, 7, 6, 6, 4, 4, 5 };
	static const U triangles[36] = { 0, 2, 1, 3, 1, 2, 1, 3, 5, 7, 5, 3, 5, 7, 4, 6, 4, 7, 4, 6, 0, 2, 0, 6, 2, 6, 3, 7, 3, 6, 4, 0, 5, 1, 5, 0 };
	static const U quads[24] = { 0, 2, 3, 1, 1, 3, 7, 5, 5, 7, 6, 4, 4, 6, 2, 0, 2, 6, 7, 3, 4, 0, 1, 5 };

	switch (nVerticesPerPrimitive)
	{
		case 2:
			return { edges, edges + 24 };
		case 3:
			return { triangles, triangles + 36 };
		case 4:
			return { quads, quads + 24 };
		default:
			return { nullptr, nullptr };
	}
}

template<typename T>
inline std::size_t AxisAlignedBox<T>::getPrimitiveCount(int nVerticesPerPrimitive) const
{
	switch (nVerticesPerPrimitive)
	{
		case 2:
			return 12;
		case 3:
			return 12;
		case 4:
			return 6;
		default:
			return 0;
	}
}

template<typename T>
template<std::output_iterator<HalfSpace<T>> O>
inline O AxisAlignedBox<T>::copyHalfSpaces(O target) const
{
	*target++ = HalfSpace<T>(Vector3<T>::MINUS_UNIT_X, Vector3<T>(minimum.x, T(0), T(0)));
	*target++ = HalfSpace<T>(Vector3<T>::UNIT_X, Vector3<T>(maximum.x, T(0), T(0)));
	*target++ = HalfSpace<T>(Vector3<T>::MINUS_UNIT_Y, Vector3<T>(T(0), minimum.y, T(0)));
	*target++ = HalfSpace<T>(Vector3<T>::UNIT_Y, Vector3<T>(T(0), maximum.y, T(0)));
	*target++ = HalfSpace<T>(Vector3<T>::MINUS_UNIT_Z, Vector3<T>(T(0), T(0), minimum.z));
	*target++ = HalfSpace<T>(Vector3<T>::UNIT_Z, Vector3<T>(T(0), T(0), maximum.z));
	return target;
}

template<typename T>
inline std::vector<HalfSpace<T>> AxisAlignedBox<T>::getHalfSpaces() const
{
	return { HalfSpace<T>(Vector3<T>::MINUS_UNIT_X, Vector3<T>(minimum.x, T(0), T(0))),
		HalfSpace<T>(Vector3<T>::UNIT_X, Vector3<T>(maximum.x, T(0), T(0))),
		HalfSpace<T>(Vector3<T>::MINUS_UNIT_Y, Vector3<T>(T(0), minimum.y, T(0))),
		HalfSpace<T>(Vector3<T>::UNIT_Y, Vector3<T>(T(0), maximum.y, T(0))),
		HalfSpace<T>(Vector3<T>::MINUS_UNIT_Z, Vector3<T>(T(0), T(0), minimum.z)),
		HalfSpace<T>(Vector3<T>::UNIT_Z, Vector3<T>(T(0), T(0), maximum.z)) };
}

template<typename T>
inline AxisAlignedBox<T>& AxisAlignedBox<T>::scaleAroundCenter(T factor)
{
	Vector3<T> center = (minimum + maximum)*T(0.5);
	minimum = (minimum - center)*factor + center;
	maximum = (maximum - center)*factor + center;
	return *this;
}

template<typename T>
inline AxisAlignedBox<T>& AxisAlignedBox<T>::setUnion(const AxisAlignedBox<T>& a, const AxisAlignedBox<T>& b)
{
	minimum.setMinimum(a.minimum, b.minimum);
	maximum.setMaximum(a.maximum, b.maximum);
	return *this;
}

template<typename T>
inline AxisAlignedBox<T>& AxisAlignedBox<T>::setIntersection(const AxisAlignedBox<T>& a, const AxisAlignedBox<T>& b)
{
	minimum.setMaximum(a.minimum, b.minimum);
	maximum.setMinimum(a.maximum, b.maximum);
	return *this;
}

template<typename T>
inline AxisAlignedBox<T>& AxisAlignedBox<T>::extendBy(const Vector3<T>& point)
{
	minimum.setMinimum(minimum, point);
	maximum.setMaximum(maximum, point);
	return *this;
}

template<typename T>
inline bool AxisAlignedBox<T>::contains(const Vector3<T>& point) const
{
	return minimum.allLessThanEqual(point) && maximum.allGreaterThanEqual(point);
}

template<typename T>
inline bool AxisAlignedBox<T>::contains(const AxisAlignedBox<T>& box) const
{
	return minimum.allLessThanEqual(box.minimum) && maximum.allGreaterThanEqual(box.maximum);
}

template<typename T>
inline bool AxisAlignedBox<T>::intersects(const AxisAlignedBox<T>& box) const
{
	return minimum.allLessThanEqual(box.maximum) && maximum.allGreaterThanEqual(box.minimum);
}

} // namespace templates

#if MATHEMATICS_DOUBLE
using AxisAlignedBox = templates::AxisAlignedBox<double>;
using AxisAlignedBoxArg = templates::AxisAlignedBox<double>::ConstArg;
using AxisAlignedBoxResult = templates::AxisAlignedBox<double>::ConstResult;
#else
using AxisAlignedBox = templates::AxisAlignedBox<float>;
using AxisAlignedBoxArg = templates::AxisAlignedBox<float>::ConstArg;
using AxisAlignedBoxResult = templates::AxisAlignedBox<float>::ConstResult;
#endif

} // namespace mathematics

namespace std {

template<typename T>
struct hash;

template<typename T>
struct hash<::mathematics::templates::AxisAlignedBox<T>>
{
	size_t operator()(const ::mathematics::templates::AxisAlignedBox<T>& box) const noexcept
	{
		hash<typename ::mathematics::templates::Vector3<T>> hasher;
		size_t seed = hasher(box.minimum) + 0x9e3779b9;
		seed ^= hasher(box.maximum) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
		return seed;
	}
};

} // namespace std

#include "OrientedBox.hpp"
#include "Sphere.hpp"
#include "SymmetricFrustum.hpp"
#include "Intersections.inl"

namespace mathematics::templates {

template<typename T>
inline AxisAlignedBox<T>::AxisAlignedBox(const OrientedBox<T>& box)
{
	Vector3<T> halfDims = abs(box.halfDims.x*box.basis[0]) + abs(box.halfDims.y*box.basis[1]) + abs(box.halfDims.z*box.basis[2]);
	minimum = box.center - halfDims;
	maximum = box.center + halfDims;
}

template<typename T>
inline Sphere<T> AxisAlignedBox<T>::getCircumscribedSphere() const
{
	Vector3<T> center = (minimum + maximum)*T(0.5);
	return { center, distance(center, maximum) };
}

template<typename T>
inline AxisAlignedBox<T>& AxisAlignedBox<T>::transform(const Matrix3<T>& matrix, bool orthogonal)
{
	*this = AxisAlignedBox<T>(OrientedBox<T>(*this, matrix, orthogonal));
	return *this;
}

template<typename T>
inline AxisAlignedBox<T>& AxisAlignedBox<T>::transform(const AffineTransform<T>& transformation, bool orthogonal)
{
	*this = AxisAlignedBox<T>(OrientedBox<T>(*this, transformation, orthogonal));
	return *this;
}

template<typename T>
inline AxisAlignedBox<T>& AxisAlignedBox<T>::extendBy(const Sphere<T>& sphere)
{
	Vector3<T> radius(sphere.radius);
	setUnion(*this, AxisAlignedBox<T>(sphere.center - radius, sphere.center + radius));
	return *this;
}

template<typename T>
inline bool AxisAlignedBox<T>::contains(const Sphere<T>& sphere) const
{
	Vector3<T> radius(sphere.radius);
	return minimum.allLessThanEqual(sphere.center - radius) && maximum.allGreaterThanEqual(sphere.center + radius);
}

template<typename T>
inline int AxisAlignedBox<T>::classify(const HalfSpace<T>& halfSpace) const
{
	return intersections::classifyAxisAlignedBoxHalfSpace(getCenter(), getHalfDimensions(), halfSpace.getNormal(), halfSpace.d);
}

template<typename T>
inline bool AxisAlignedBox<T>::intersects(const HalfSpace<T>& halfSpace) const
{
	return intersections::testAxisAlignedBoxHalfSpace(getCenter(), getHalfDimensions(), halfSpace.getNormal(), halfSpace.d);
}

template<typename T>
inline bool AxisAlignedBox<T>::intersects(const Plane<T>& plane) const
{
	return intersections::testAxisAlignedBoxPlane(getCenter(), getHalfDimensions(), plane.getNormal(), plane.d);
}

template<typename T>
inline bool AxisAlignedBox<T>::intersects(const Triangle3<T>& triangle) const
{
	return intersections::testAxisAlignedBoxTriangle(getCenter(), getHalfDimensions(), triangle.vertices[0], triangle.vertices[1], 
		triangle.vertices[2]);
}

template<typename T>
inline bool AxisAlignedBox<T>::intersects(const OrientedBox<T>& box) const
{
	return intersections::testOrientedBoxAxisAlignedBox(box.center, box.basis, box.halfDims, getCenter(), getHalfDimensions());
}

template<typename T>
inline bool AxisAlignedBox<T>::intersects(const Sphere<T>& sphere) const
{
	return intersections::testAxisAlignedBoxSphere(minimum, maximum, sphere.center, sphere.radius);
}

template<typename T>
inline bool AxisAlignedBox<T>::intersects(const SymmetricFrustum<T>& frustum) const
{
	return frustum.intersects(*this);
}

} // namespace mathematics::templates
