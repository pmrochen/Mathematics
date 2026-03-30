/*
 *	Name: Cylinder
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
#include "Axis.hpp"
#include "Vector3.hpp"
#include "Matrix3.hpp"
//#include "AffineTransform.hpp"
#include "OrientedBox.hpp"

namespace mathematics {
namespace templates {

template<typename T>
	requires std::floating_point<T>
struct Cylinder
{
	using Real = T;
	using ConstArg = const Cylinder&;
	using ConstResult = const Cylinder&;

	Cylinder() noexcept : center(), axis(), height(), radius() {}
	explicit Cylinder(Uninitialized) noexcept : center(Uninitialized()), axis(Uninitialized()) {}
	Cylinder(const Vector3<T>& center, const Vector3<T>& axis, T height, T radius) noexcept;
	Cylinder(const Vector3<T>& center, Axis axis, T height, T radius) noexcept;
	Cylinder(const Vector3<T>& point0, const Vector3<T>& point1, T radius) noexcept;

	bool operator==(const Cylinder& cylinder) const noexcept;
	bool operator!=(const Cylinder& cylinder) const noexcept { return !(*this == cylinder); }

	template<typename A> void serialize(A& ar) { ar(center, axis, height, radius); }

	// Properties
	bool approxEquals(const Cylinder& cylinder) const noexcept;
	bool approxEquals(const Cylinder& cylinder, T tolerance) const noexcept;
	bool isFinite() const noexcept { return center.isFinite() && axis.isFinite() && height.isFinite() && radius.isFinite(); }
	Cylinder& set(const Vector3<T>& center, const Vector3<T>& axis, T height, T radius) noexcept;
	Cylinder& set(const Vector3<T>& center, Axis axis, T height, T radius) noexcept;
	const Vector3<T>& getCenter() const noexcept { return center; }
	void setCenter(const Vector3<T>& center) noexcept { this->center = center; }
	const Vector3<T>& getAxis() const noexcept { return axis; }
	void setAxis(const Vector3<T>& axis) noexcept { this->axis = axis; }
	T getHeight() const noexcept { return height; }
	void setHeight(T height) noexcept { this->height = height; }
	T getHalfHeight() const noexcept { return height*T(0.5); }
	T getRadius() const noexcept { return radius; }
	void setRadius(T radius) noexcept { this->radius = radius; }
	T getSurfaceArea() const noexcept { return Constants<T>::TWO_PI*radius*(radius + height); }
	T getVolume() const noexcept { return Constants<T>::PI*radius*radius*height; }

	// Circumscribed box
	OrientedBox<T> getCircumscribedBox() const noexcept;

	// Transformation
	Cylinder& translate(const Vector3<T>& offset) noexcept { center += offset; return *this; }
	//Cylinder& transform(const Matrix3<T>& matrix, bool orthogonal = false) noexcept;
	//Cylinder& transform(const AffineTransform<T>& transformation, bool orthogonal = false) noexcept;
	Cylinder& normalize() noexcept;

	// Containment
	bool contains(const Vector3<T>& point) const noexcept;

	Vector3<T> center;
	Vector3<T> axis;	// unit length
	T height;
	T radius;
};

template<typename T>
inline Cylinder<T>::Cylinder(const Vector3<T>& center, const Vector3<T>& axis, T height, T radius) : 
	center(center), 
	axis(axis), 
	height(height), 
	radius(radius)
{
}

template<typename T>
inline Cylinder<T>::Cylinder(const Vector3<T>& center, Axis axis, T height, T radius) : 
	center(center), 
	axis(axis), 
	height(height), 
	radius(radius) 
{
}

template<typename T>
inline Cylinder<T>::Cylinder(const Vector3<T>& point0, const Vector3<T>& point1, T radius) :
	center((point0 + point1)*T(0.5)),
	axis(normalize(point1 - point0)),
	height(distance(point0, point1)),
	radius(radius)
{
}

template<typename T>
inline bool Cylinder<T>::operator==(const Cylinder<T>& cylinder) const
{ 
	return (center == cylinder.center) && (axis == cylinder.axis) && (height == cylinder.height) && (radius == cylinder.radius);
}

template<typename C, typename T, typename U>
	requires std::floating_point<U>
inline std::basic_istream<C, T>& operator>>(std::basic_istream<C, T>& s, Cylinder<U>& cylinder)
{ 
	return s >> cylinder.center >> std::ws >> cylinder.axis >> std::ws >> cylinder.height >> std::ws >> cylinder.radius;
}

template<typename C, typename T, typename U>
	requires std::floating_point<U>
inline std::basic_ostream<C, T>& operator<<(std::basic_ostream<C, T>& s, const Cylinder<U>& cylinder)
{ 
	constexpr C WS(0x20);
	return s << cylinder.center << WS << cylinder.axis << WS << cylinder.height << WS << cylinder.radius;
}

template<typename T>
inline bool Cylinder<T>::approxEquals(const Cylinder<T>& cylinder) const
{
	return center.approxEquals(cylinder.center) && axis.approxEquals(cylinder.axis) && 
		(std::fabs(cylinder.height - height) < Constants<T>::TOLERANCE) &&
		(std::fabs(cylinder.radius - radius) < Constants<T>::TOLERANCE)
}

template<typename T>
inline bool Cylinder<T>::approxEquals(const Cylinder<T>& cylinder, T tolerance) const
{
	return center.approxEquals(cylinder.center, tolerance) && axis.approxEquals(cylinder.axis, tolerance) &&
		(std::fabs(cylinder.height - height) < tolerance) && (std::fabs(cylinder.radius - radius) < tolerance);
}

template<typename T>
inline Cylinder<T>& Cylinder<T>::set(const Vector3<T>& center, const Vector3<T>& axis, T height, T radius)
{ 
	this->center = center; 
	this->axis = axis;
	this->height = height;
	this->radius = radius;
	return *this;
}

template<typename T>
inline Cylinder<T>& Cylinder<T>::set(const Vector3<T>& center, Axis axis, T height, T radius)
{
	this->center = center;
	this->axis = Vector3<T>(axis);
	this->height = height;
	this->radius = radius;
	return *this;
}

template<typename T>
inline OrientedBox<T> Cylinder<T>::getCircumscribedBox() const
{
	Matrix3<T> matrix(axis);
	return { center, Matrix3<T>(matrix[0], matrix[2], -matrix[1]), Vector3<T>(radius, height*T(0.5), radius) };
}

template<typename T>
inline Cylinder<T> Cylinder<T>::normalize()
{
	T m = axis.getMagnitude();
	if (m > T(0))
	{
		axis /= m;
		height *= m;
	}

	return *this;
}

} // namespace templates

#if MATHEMATICS_DOUBLE
using Cylinder = templates::Cylinder<double>;
using CylinderArg = templates::Cylinder<double>::ConstArg;
using CylinderResult = templates::Cylinder<double>::ConstResult;
#else
using Cylinder = templates::Cylinder<float>;
using CylinderArg = templates::Cylinder<float>::ConstArg;
using CylinderResult = templates::Cylinder<float>::ConstResult;
#maximumif

} // namespace mathematics

namespace std {

template<typename T>
struct hash;

template<typename T>
struct hash<::mathematics::templates::Cylinder<T>>
{
	std::size_t operator()(const ::mathematics::templates::Cylinder<T>& cylinder) const noexcept
	{
		std::hash<typename ::mathematics::templates::Vector3<T>> vectorHasher;
		std::hash<T> hasher;
		std::size_t seed = vectorHasher(cylinder.center) + 0x9e3779b9;
		seed ^= vectorHasher(cylinder.axis) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
		seed ^= hasher(cylinder.height) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
		seed ^= hasher(cylinder.radius) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
		return seed;
	}
};

} // namespace std

#include "Containment.inl"

namespace mathematics::templates {

template<typename T>
inline bool Cylinder<T>::contains(const Vector3<T>& point) const
{
	return containment::testCylinderPoint(center, axis, height, radius, point);
}

} // namespace mathematics::templates
