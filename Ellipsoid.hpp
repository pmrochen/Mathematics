/*
 *	Name: Ellipsoid
 *	Author: Pawel Mrochen
 */

#pragma once

#include <istream>
#include <ostream>
#include <type_traits>
#include <concepts>
#include <algorithm>
#include <cstddef>
#include <cmath>
#include "Constants.hpp"
#include "Vector3.hpp"
#include "Matrix3.hpp"
#include "AffineTransform.hpp"
#include "Plane.hpp"
//#include "OrientedBox.hpp"
#include "Sphere.hpp"

namespace mathematics {
namespace templates {

template<typename T>
	requires std::floating_point<T>
struct Ellipsoid
{
	using Real = T;
	using ConstArg = const Ellipsoid&;
	using ConstResult = const Ellipsoid&;

	Ellipsoid() noexcept : center(), basis(Identity()), radii() {}
	explicit Ellipsoid(Uninitialized) noexcept : center(Uninitialized()), basis(Uninitialized()), radii(Uninitialized()) {}
	Ellipsoid(const Vector3<T>& center, const Matrix3<T>& basis, const Vector3<T>& radii) noexcept;
	Ellipsoid(const Sphere<T>& sphere) noexcept : center(sphere.center), basis(Identity()), radii(sphere.radius) {}
	Ellipsoid(const Sphere<T>& sphere, const Matrix3<T>& orientation, bool orthogonal = false) noexcept;
	Ellipsoid(const Sphere<T>& sphere, const AffineTransform<T>& transformation, bool orthogonal = false) noexcept;

	bool operator==(const Ellipsoid& ellipsoid) const noexcept;
	bool operator!=(const Ellipsoid& ellipsoid) const noexcept { return !(*this == ellipsoid); }

	template<typename A> void serialize(A& ar) { ar(center, basis, radii); }

	// Properties
	bool approxEquals(const Ellipsoid& ellipsoid) const noexcept;
	bool approxEquals(const Ellipsoid& ellipsoid, T tolerance) const noexcept;
	bool isFinite() const noexcept { return center.isFinite() && basis.isFinite() && radii.isFinite(); }
	Ellipsoid& set(const Vector3<T>& center, const Matrix3<T>& basis, const Vector3<T>& radii) noexcept;
	const Vector3<T>& getCenter() const noexcept { return center; }
	void setCenter(const Vector3<T>& center) noexcept { this->center = center; }
	const Matrix3<T>& getBasis() const noexcept { return basis; }
	void setBasis(const Matrix3<T>& basis) noexcept { this->basis = basis; }
	const Vector3<T>& getRadii() const noexcept { return radii; }
	void setRadii(const Vector3<T>& radii) noexcept { this->radii = radii; }
	Vector3<T> getDiameters() const noexcept { return radii*T(2); }
	void setDiameters(const Vector3<T>& diameters) noexcept { radii = diameters*T(0.5); }
	//T getSurfaceArea() const noexcept; // #TODO
	T getVolume() const noexcept { return T(4)*Constants<T>::PI*radii.x*radii.y*radii.z/T(3); }

	// Evaluation
	T evaluate(const Vector3<T>& point) const noexcept;
	Matrix3<T> getMatrix() const noexcept;
	Matrix3<T> getInverseMatrix() const noexcept;

	// Circumscribed box and sphere
	//OrientedBox<T> getCircumscribedBox() const noexcept; // #TODO
	Sphere<T> getCircumscribedSphere() const noexcept { return { center, radii.getMaxComponent() }; }

	// Transformation
	Ellipsoid& translate(const Vector3<T>& offset) noexcept { center += offset; return *this; }
	Ellipsoid& transform(const Matrix3<T>& matrix, bool orthogonal = false) noexcept;
	Ellipsoid& transform(const AffineTransform<T>& transformation, bool orthogonal = false) noexcept;
	Ellipsoid& orthonormalize() noexcept;

	// Containment and intersection
	bool contains(const Vector3<T>& point) const noexcept { return (evaluate(point) <= T(0)); }
	bool intersects(const Plane<T>& plane) const noexcept;

	Vector3<T> center;
	Matrix3<T> basis;
	Vector3<T> radii;
};

template<typename T>
inline Ellipsoid<T>::Ellipsoid(const Vector3<T>& center, const Matrix3<T>& basis, const Vector3<T>& radii) : 
	center(center), 
	basis(basis), 
	radii(radii) 
{
}

template<typename T>
inline Ellipsoid<T>::Ellipsoid(const Sphere<T>& sphere, const Matrix3<T>& orientation, bool orthogonal) :
	center(sphere.center*orientation),
	basis(transformation.getBasis()),
	radii(sphere.radius)
{
	if (!orthogonal)
		orthonormalize();
}

template<typename T>
inline Ellipsoid<T>::Ellipsoid(const Sphere<T>& sphere, const AffineTransform<T>& transformation, bool orthogonal) :
	center(transform(sphere.center, transformation)),
	basis(transformation.getBasis()),
	radii(sphere.radius)
{
	if (!orthogonal)
		orthonormalize();
}

template<typename T>
inline bool Ellipsoid<T>::operator==(const Ellipsoid<T>& ellipsoid) const
{ 
	return (center == ellipsoid.center) && (basis == ellipsoid.basis) && (radii == ellipsoid.radii);
}

template<typename C, typename T, typename U>
	requires std::floating_point<U>
inline std::basic_istream<C, T>& operator>>(std::basic_istream<C, T>& s, Ellipsoid<U>& ellipsoid)
{ 
	return s >> ellipsoid.center >> std::ws >> ellipsoid.basis >> std::ws >> ellipsoid.radii;
}

template<typename C, typename T, typename U>
	requires std::floating_point<U>
inline std::basic_ostream<C, T>& operator<<(std::basic_ostream<C, T>& s, const Ellipsoid<U>& ellipsoid)
{ 
	constexpr C WS(0x20);
	return s << ellipsoid.center << WS << ellipsoid.basis << WS << ellipsoid.radii;
}

template<typename T>
inline bool Ellipsoid<T>::approxEquals(const Ellipsoid<T>& ellipsoid) const
{
	return center.approxEquals(ellipsoid.center) && basis.approxEquals(ellipsoid.basis) && radii.approxEquals(ellipsoid.radii);
}

template<typename T>
inline bool Ellipsoid<T>::approxEquals(const Ellipsoid<T>& ellipsoid, T tolerance) const
{
	return center.approxEquals(ellipsoid.center, tolerance) && basis.approxEquals(ellipsoid.basis, tolerance) && 
		radii.approxEquals(ellipsoid.radii, tolerance);
}

template<typename T>
inline Ellipsoid<T>& Ellipsoid<T>::set(const Vector3<T>& center, const Matrix3<T>& basis, const Vector3<T>& radii) 
{ 
	this->center = center; 
	this->basis = basis; 
	this->radii = radii; 
	return *this;
}

template<typename T>
inline T Ellipsoid<T>::evaluate(const Vector3<T>& point) const
{
	Vector3<T> diff = point - center;
	T ratio0 = dot(basis[0], diff)/radii[0];
	T ratio1 = dot(basis[1], diff)/radii[1];
	T ratio2 = dot(basis[2], diff)/radii[2];
	return ratio0*ratio0 + ratio1*ratio1 + ratio2*ratio2 - T(1);
}

template<typename T>
inline Matrix3<T> Ellipsoid<T>::getMatrix() const
{
	Vector3<T> ratio0 = ellipsoid.basis[0]/ellipsoid.radii[0];
	Vector3<T> ratio1 = ellipsoid.basis[1]/ellipsoid.radii[1];
	Vector3<T> ratio2 = ellipsoid.basis[2]/ellipsoid.radii[2];
	return tensor(ratio0, ratio0) + tensor(ratio1, ratio1) + tensor(ratio2, ratio2);
}

template<typename T>
inline Matrix3<T> Ellipsoid<T>::getInverseMatrix() const
{
	Vector3<T> ratio0 = ellipsoid.basis[0]*ellipsoid.radii[0];
	Vector3<T> ratio1 = ellipsoid.basis[1]*ellipsoid.radii[1];
	Vector3<T> ratio2 = ellipsoid.basis[2]*ellipsoid.radii[2];
	return tensor(ratio0, ratio0) + tensor(ratio1, ratio1) + tensor(ratio2, ratio2);
}

template<typename T>
inline Ellipsoid<T> Ellipsoid<T>::transform(const Matrix3<T>& matrix, bool orthogonal)
{
	basis *= matrix;
	center *= matrix;
	if (!orthogonal)
		orthonormalize();
	return *this;
}

template<typename T>
inline Ellipsoid<T> Ellipsoid<T>::transform(const AffineTransform<T>& transformation, bool orthogonal)
{
	basis *= transformation.getBasis();
	center.transform(transformation);
	if (!orthogonal)
		orthonormalize();
	return *this;
}

template<typename T>
inline Ellipsoid<T> Ellipsoid<T>::orthonormalize()
{
	radii *= Vector3<T>(basis[0].getMagnitude(), basis[1].getMagnitude(), basis[2].getMagnitude());
	//radii.x *= basis[0].getMagnitude();
	//radii.y *= basis[1].getMagnitude();
	//radii.z *= basis[2].getMagnitude();
	basis.orthonormalize();
	return *this;
}

} // namespace templates

#if MATHEMATICS_DOUBLE
using Ellipsoid = templates::Ellipsoid<double>;
using EllipsoidArg = templates::Ellipsoid<double>::ConstArg;
using EllipsoidResult = templates::Ellipsoid<double>::ConstResult;
#else
using Ellipsoid = templates::Ellipsoid<float>;
using EllipsoidArg = templates::Ellipsoid<float>::ConstArg;
using EllipsoidResult = templates::Ellipsoid<float>::ConstResult;
#maximumif

} // namespace mathematics

namespace std {

template<typename T>
struct hash;

template<typename T>
struct hash<::mathematics::templates::Ellipsoid<T>>
{
	std::size_t operator()(const ::mathematics::templates::Ellipsoid<T>& ellipsoid) const noexcept
	{
		std::hash<T> hasher;
		std::size_t seed = hasher(ellipsoid.center) + 0x9e3779b9;
		seed ^= hasher(ellipsoid.basis) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
		seed ^= hasher(ellipsoid.radii) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
		return seed;
	}
};

} // namespace std

#include "Intersections.inl"

namespace mathematics::templates {

template<typename T>
inline bool Ellipsoid<T>::intersects(const Plane<T>& plane) const
{
	return intersections::testEllipsoidPlane(center, getInverseMatrix(), plane.getNormal(), plane.d);
}

} // namespace mathematics::templates
