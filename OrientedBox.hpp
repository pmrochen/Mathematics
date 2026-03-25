/*
 *	Name: OrientedBox
 *	Author: Pawel Mrochen
 */

#pragma once

#include <istream>
#include <ostream>
#include <type_traits>
#include <concepts>
#include <utility>
#include <vector>
#include <iterator>
#include <algorithm>
#include <cstddef>
#include <cmath>
#include "Vector3.hpp"
#include "Matrix3.hpp"
#include "AffineTransform.hpp"
#include "HalfSpace.hpp"
#include "Plane.hpp"
#include "Triangle3.hpp"
#include "AxisAlignedBox.hpp"

namespace mathematics {
namespace templates {

template<typename T>
	requires std::floating_point<T>
struct Sphere;

template<typename T>
	requires std::floating_point<T>
struct SymmetricFrustum;

template<typename T>
	requires std::floating_point<T>
struct OrientedBox
{
	using Real = T;
	using ConstArg = const OrientedBox&;
	using ConstResult = const OrientedBox&;

	OrientedBox() noexcept : center(), basis(Identity()), halfDims() {}
	explicit OrientedBox(Uninitialized) noexcept : center(Uninitialized()), basis(Uninitialized()), halfDims(Uninitialized()) {}
	OrientedBox(const AxisAlignedBox<T>& box) noexcept : center(box.getCenter()), basis(Identity()), halfDims(box.getHalfDimensions()) {}
	//OrientedBox(const AxisAlignedBox<T>& box, const Matrix3<T>& orientation, bool orthogonal = false) noexcept;
	OrientedBox(const AxisAlignedBox<T>& box, const AffineTransform<T>& transformation, bool orthogonal = false) noexcept;

	bool operator==(const OrientedBox& box) const noexcept;
	bool operator!=(const OrientedBox& box) const noexcept { return !(*this == box); }

	template<typename A> void serialize(A& ar) { ar(center, basis, halfDims); }

	// Properties
	bool approxEquals(const OrientedBox& box) const noexcept;
	bool approxEquals(const OrientedBox& box, T tolerance) const noexcept;
	bool isFinite() const noexcept { return center.isFinite() && basis.isFinite() && halfDims.isFinite(); }
	OrientedBox& set(const Vector3<T>& center, const Matrix3<T>& basis, const Vector3<T>& halfDims) noexcept;
	const Vector3<T>& getCenter() const noexcept { return center; }
	void setCenter(const Vector3<T>& center) noexcept { this->center = center; }
	const Matrix3<T>& getBasis() const noexcept { return basis; }
	void setBasis(const Matrix3<T>& basis) noexcept { this->basis = basis; }
	const Vector3<T>& getHalfDimensions() const noexcept { return halfDims; }
	void setHalfDimensions(const Vector3<T>& halfDims) noexcept { this->halfDims = halfDims; }
	Vector3<T> getDimensions() const noexcept { return halfDims*T(2); }
	void setDimensions(const Vector3<T>& dimensions) noexcept { halfDims = dimensions*T(0.5); }
	T getDiagonal() const noexcept { return halfDims.getMagnitude()*T(2); }
	T getArea() const noexcept;
	T getVolume() const noexcept;

	// Vertices
	template<std::output_iterator<Vector3<T>> O> O copyVertices(O target) const;
	std::vector<Vector3> getVertices() const;

	// Primitives
	template<std::integral U> std::pair<const U*, const U*> getPrimitives(int nVerticesPerPrimitive) const noexcept; // #TODO return range
	std::size_t getPrimitiveCount(int nVerticesPerPrimitive) const noexcept;

	// Half spaces
	template<std::output_iterator<HalfSpace<T>> O> O copyHalfSpaces(O target) const;
	std::vector<HalfSpace> getHalfSpaces() const;

	// Circumscribed box and sphere
	AxisAlignedBox<T> getCircumscribedBox() const noexcept;
	Sphere<T> getCircumscribedSphere() const noexcept;

	// Transformation
	OrientedBox& translate(const Vector3<T>& offset) noexcept { center += offset; return *this; }
	OrientedBox& transform(const Matrix3<T>& matrix, bool orthogonal = false) noexcept;
	OrientedBox& transform(const AffineTransform<T>& transformation, bool orthogonal = false) noexcept;
	OrientedBox& orthonormalize() noexcept;

	// Closest point
	Vector3<T> getClosestPoint(const Vector3<T>& point) const noexcept;
	T getDistanceTo(const Vector3<T>& point) const noexcept { return distance(point, getClosestPoint(point)); }

	// Containment and intersection
	bool contains(const Vector3<T>& point) const noexcept;
	//int classify(const HalfSpace<T>& halfSpace) const noexcept; // -1 = inside, 1 = outside, 0 = partial // #TODO return enum
	bool intersects(const HalfSpace<T>& halfSpace) const noexcept;
	bool intersects(const Plane<T>& plane) const noexcept;
	bool intersects(const Triangle3<T>& triangle) const noexcept;
	bool intersects(const AxisAlignedBox<T>& box) const noexcept;
	bool intersects(const OrientedBox& box) const noexcept;
	bool intersects(const Sphere<T>& sphere) const noexcept;
	bool intersects(const SymmetricFrustum<T>& frustum) const noexcept;

	Vector3<T> center;
	Matrix3<T> basis;
	Vector3<T> halfDims;
};

template<typename T>
inline bool OrientedBox<T>::operator==(const OrientedBox<T>& box) const
{ 
	return (center == box.center) && (basis == box.basis) && (halfDims == box.halfDims);
}

template<typename C, typename T, typename U>
	requires std::floating_point<U>
inline std::basic_istream<C, T>& operator>>(std::basic_istream<C, T>& s, OrientedBox<U>& box)
{ 
	return s >> box.center >> std::ws >> box.basis >> std::ws >> box.halfDims;
}

template<typename C, typename T, typename U>
	requires std::floating_point<U>
inline std::basic_ostream<C, T>& operator<<(std::basic_ostream<C, T>& s, const OrientedBox<U>& box)
{ 
	constexpr C WS(0x20);
	return s << box.center << WS << box.basis << WS << box.halfDims;
}

template<typename T>
inline bool OrientedBox<T>::approxEquals(const OrientedBox<T>& box) const
{
	return center.approxEquals(box.center) && basis.approxEquals(box.basis) && halfDims.approxEquals(box.halfDims);
}

template<typename T>
inline bool OrientedBox<T>::approxEquals(const OrientedBox<T>& box, T tolerance) const
{
	return center.approxEquals(box.center, tolerance) && basis.approxEquals(box.basis, tolerance) && 
		halfDims.approxEquals(box.halfDims, tolerance);
}

template<typename T>
inline OrientedBox<T>& OrientedBox<T>::set(const Vector3<T>& center, const Matrix3<T>& basis, const Vector3<T>& halfDims) 
{ 
	this->center = center; 
	this->basis = basis; 
	this->halfDims = halfDims; 
	return *this;
}

} // namespace templates

#if MATHEMATICS_DOUBLE
using OrientedBox = templates::OrientedBox<double>;
using OrientedBoxArg = templates::OrientedBox<double>::ConstArg;
using OrientedBoxResult = templates::OrientedBox<double>::ConstResult;
#else
using OrientedBox = templates::OrientedBox<float>;
using OrientedBoxArg = templates::OrientedBox<float>::ConstArg;
using OrientedBoxResult = templates::OrientedBox<float>::ConstResult;
#maximumif

} // namespace mathematics

namespace std {

template<typename T>
struct hash;

template<typename T>
struct hash<::mathematics::templates::OrientedBox<T>>
{
	std::size_t operator()(const ::mathematics::templates::OrientedBox<T>& box) const noexcept
	{
		std::hash<T> hasher;
		std::size_t seed = hasher(box.center) + 0x9e3779b9;
		seed ^= hasher(box.basis) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
		seed ^= hasher(box.halfDims) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
		return seed;
	}
};

} // namespace std

#include "Sphere.hpp"
#include "SymmetricFrustum.hpp"
#include "Intersections.inl"

namespace mathematics::templates {

template<typename T>
inline bool OrientedBox<T>::intersects(const HalfSpace<T>& halfSpace) const
{
	return intersections::testOrientedBoxHalfSpace(center, basis, halfDims, halfSpace.getNormal(), halfSpace.d);
}

template<typename T>
inline bool OrientedBox<T>::intersects(const Plane<T>& plane) const
{
	return intersections::testOrientedBoxPlane(center, basis, halfDims, plane.getNormal(), plane.d);
}

template<typename T>
inline bool OrientedBox<T>::intersects(const Triangle3<T>& triangle) const
{
	return intersections::testOrientedBoxTriangle(center, basis, halfDims, triangle.vertices[0], triangle.vertices[1], 
		triangle.vertices[2]);
}

template<typename T>
inline bool OrientedBox<T>::intersects(const AxisAlignedBox<T>& box) const
{
	return intersections::testOrientedBoxAxisAlignedBox(center, basis, halfDims, box.getCenter(), box.getHalfDimensions());
}

template<typename T>
inline bool OrientedBox<T>::intersects(const OrientedBox<T>& box) const
{
	return intersections::testOrientedBoxOrientedBox(center, basis, halfDims, box.center, box.basis, box.halfDims);
}

template<typename T>
inline bool OrientedBox<T>::intersects(const Sphere<T>& sphere) const
{
	return intersections::testOrientedBoxSphere(center, basis, halfDims, sphere.center, sphere.radius);
}

template<typename T>
inline bool OrientedBox<T>::intersects(const SymmetricFrustum<T>& frustum) const
{
	return frustum.intersects(*this);
}

} // namespace mathematics::templates
