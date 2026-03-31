/*
 *	Name: Sphere
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
#include <cstddef>
#include <cmath>
#include "Constants.hpp"
#include "Vector3.hpp"
#include "HalfSpace.hpp"
#include "Plane.hpp"
#include "Triangle3.hpp"
#include "AxisAlignedBox.hpp"
#include "OrientedBox.hpp"

namespace mathematics {
namespace templates {

template<typename T>
	requires std::floating_point<T>
struct Ellipsoid;

template<typename T>
	requires std::floating_point<T>
struct Cone;

template<typename T>
	requires std::floating_point<T>
struct SymmetricFrustum;

template<typename T>
	requires std::floating_point<T>
struct Sphere
{
	using Real = T;
	using ConstArg = const Sphere&;
	using ConstResult = const Sphere&;

	Sphere() noexcept : center(), radius() {}
	explicit Sphere(Uninitialized) noexcept : center(Uninitialized()) {}
	Sphere(const Vector3<T>& center, T radius) noexcept : center(center), radius(radius) {}
	explicit Sphere(const Ellipsoid<T>& ellipsoid) noexcept;

	bool operator==(const Sphere& sphere) const noexcept { return (center == sphere.center) && (radius == sphere.radius); }
	bool operator!=(const Sphere& sphere) const noexcept { return !(*this == sphere); }

	template<typename A> void serialize(A& ar) { ar(center, radius); }

	// Properties
	bool approxEquals(const Sphere& sphere) const noexcept;
	bool approxEquals(const Sphere& sphere, T tolerance) const noexcept;
	bool isFinite() const noexcept { return center.isFinite() && radius.isFinite(); }
	Sphere& set(const Vector3<T>& center, T radius) noexcept { this->center = center; this->radius = radius; return *this; }
	const Vector3<T>& getCenter() const noexcept { return center; }
	void setCenter(const Vector3<T>& center) noexcept { this->center = center; }
	T getRadius() const noexcept { return radius; }
	void setRadius(T radius) noexcept { this->radius = radius; }
	T getDiameter() const noexcept { return radius*T(2); }
	void setDiameter(T diameter) noexcept { radius = diameter*T(0.5); }
	T getSurfaceArea() const noexcept { return T(4)*Constants<T>::PI*radius*radius; }
	T getVolume() const noexcept { return T(4)*Constants<T>::PI*radius*radius*radius/T(3); }

	// Circumscribed box
	AxisAlignedBox<T> getCircumscribedBox() const noexcept;

	// Transformation
	Sphere& translate(const Vector3<T>& offset) { center += offset; return *this; }
	Sphere& transform(const Matrix3<T>& matrix, bool orthogonal = false) noexcept;
	Sphere& transform(const AffineTransform<T>& transformation, bool orthogonal = false) noexcept;

	// Closest point
	//Vector3<T> getClosestPoint(const Vector3<T>& point) const noexcept; // #TODO
	T getDistanceTo(const Vector3<T>& point) const { return std::max(distance(point, center) - radius, T(0)); }
	T getSignedDistanceTo(const Vector3<T>& point) const noexcept { return distance(point, center) - radius; }

	// Containment and intersection
	bool contains(const Vector3<T>& point) const noexcept { return (distanceSquared(point, center) <= radius*radius); }
	bool intersects(const HalfSpace<T>& halfSpace) const noexcept;
	bool intersects(const Plane<T>& plane) const noexcept;
	bool intersects(const Triangle3<T>& triangle) const noexcept;
	bool intersects(const AxisAlignedBox<T>& box) const noexcept;
	bool intersects(const OrientedBox<T>& box) const noexcept;
	bool intersects(const Sphere& sphere) const noexcept;
	bool intersects(const Cone<T>& cone) const noexcept;
	bool intersects(const SymmetricFrustum<T>& frustum) const noexcept;

	Vector3<T> center;
	T radius;
};

template<typename C, typename T, typename U>
	requires std::floating_point<U>
inline std::basic_istream<C, T>& operator>>(std::basic_istream<C, T>& s, Sphere<U>& sphere)
{ 
	return s >> sphere.center >> std::ws >> sphere.radius;
}

template<typename C, typename T, typename U>
	requires std::floating_point<U>
inline std::basic_ostream<C, T>& operator<<(std::basic_ostream<C, T>& s, const Sphere<U>& sphere)
{ 
	constexpr C WS(0x20);
	return s << sphere.center << WS << sphere.radius;
}

template<typename T>
inline bool Sphere<T>::approxEquals(const Sphere<T>& sphere) const
{
	return center.approxEquals(sphere.center) && 
		(std::fabs(sphere.radius - radius) < Constants<T>::TOLERANCE);
}

template<typename T>
inline bool Sphere<T>::approxEquals(const Sphere<T>& sphere, T tolerance) const
{
	return center.approxEquals(sphere.center, tolerance) && 
		(std::fabs(sphere.radius - radius) < tolerance);
}

template<typename T>
inline AxisAlignedBox<T> Sphere<T>::getCircumscribedBox() const
{
	Vector3<T> halfDims(radius);
	return AxisAlignedBox<T>(center - halfDims, center + halfDims);
}

template<typename T>
inline bool Sphere<T>::intersects(const HalfSpace<T>& halfSpace) const
{
	return ((dot(halfSpace.getNormal(), center) + halfSpace.d) <= radius);
}

template<typename T>
inline bool Sphere<T>::intersects(const Plane<T>& plane) const
{
	return (std::fabs(dot(plane.getNormal(), center) + plane.d) <= radius);
}

template<typename T>
inline bool Sphere<T>::intersects(const Sphere<T>& sphere) const
{
	T d = sphere.radius + radius;
	return (distanceSquared(sphere.center, center) <= d*d);
}

} // namespace templates

#if MATHEMATICS_DOUBLE
using Sphere = templates::Sphere<double>;
using SphereArg = templates::Sphere<double>::ConstArg;
using SphereResult = templates::Sphere<double>::ConstResult;
#else
using Sphere = templates::Sphere<float>;
using SphereArg = templates::Sphere<float>::ConstArg;
using SphereResult = templates::Sphere<float>::ConstResult;
#endif

} // namespace mathematics

namespace std {

template<typename T>
struct hash;

template<typename T>
struct hash<::mathematics::templates::Sphere<T>>
{
	size_t operator()(const ::mathematics::templates::Sphere<T>& sphere) const noexcept
	{
		size_t seed = hash<typename ::mathematics::templates::Vector3<T>>()(sphere.center) + 0x9e3779b9;
		seed ^= hash<T>()(sphere.radius) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
		return seed;
	}
};

} // namespace std

#include "Ellipsoid.hpp"
#include "Cone.hpp"
#include "SymmetricFrustum.hpp"
#include "Distances.inl"
#include "Intersections.inl"

namespace mathematics::templates {

template<typename T>
inline Sphere<T>::Sphere(const Ellipsoid<T>& ellipsoid) :
	center(ellipsoid.center),
	radius(ellipsoid.radii.getMaxComponent())
{
}

template<typename T>
inline Sphere<T>& Sphere<T>::transform(const Matrix3<T>& matrix, bool orthogonal)
{
	if (!orthogonal)
		*this = Sphere<T>(Ellipsoid<T>(*this, matrix, orthogonal));
	return *this;
}

template<typename T>
inline Sphere<T>& Sphere<T>::transform(const AffineTransform<T>& transformation, bool orthogonal)
{
	if (orthogonal)
		center.transform(transformation);
	else
		*this = Sphere<T>(Ellipsoid<T>(*this, transformation, orthogonal));
	return *this;
}

template<typename T>
inline bool Sphere<T>::intersects(const Triangle3<T>& triangle) const
{
	return (distances::getPointTriangleSquared(center, triangle.vertex0, triangle.vertex1, triangle.vertex2) <= radius*radius);
}

template<typename T>
inline bool Sphere<T>::intersects(const AxisAlignedBox<T>& box) const
{
	return intersections::testAxisAlignedBoxSphere(box.minimum, box.maximum, center, radius);
}

template<typename T>
inline bool Sphere<T>::intersects(const OrientedBox<T>& box) const
{
	return intersections::testOrientedBoxSphere(box.center, box.basis, box.halfDims, center, radius);
}

template<typename T>
inline bool Sphere<T>::intersects(const Cone<T>& cone) const
{
	return intersections::testConeSphere(cone.vertex, cone.axis, cone.height, cone.radius, center, radius);
}

template<typename T>
inline bool Sphere<T>::intersects(const SymmetricFrustum<T>& frustum) const
{
	return frustum.intersects(*this);
}

} // namespace mathematics::templates
