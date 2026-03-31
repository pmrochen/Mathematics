/*
 *	Name: Capsule
 *	Author: Pawel Mrochen
 */

#pragma once

#include <istream>
#include <ostream>
#include <type_traits>
#include <concepts>
#include <algorithm>
#include <functional>
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
struct Capsule
{
	using Real = T;
	using ConstArg = const Capsule&;
	using ConstResult = const Capsule&;

	Capsule() noexcept : center(), axis(), height(), radius() {}
	explicit Capsule(Uninitialized) noexcept : center(Uninitialized()), axis(Uninitialized()) {}
	Capsule(const Vector3<T>& center, const Vector3<T>& axis, T height, T radius) noexcept;
	Capsule(const Vector3<T>& center, Axis axis, T height, T radius) noexcept;
	Capsule(const Vector3<T>& point0, const Vector3<T>& point1, T radius) noexcept;

	bool operator==(const Capsule& capsule) const noexcept;
	bool operator!=(const Capsule& capsule) const noexcept { return !(*this == capsule); }

	template<typename A> void serialize(A& ar) { ar(center, axis, height, radius); }

	// Properties
	bool approxEquals(const Capsule& capsule) const noexcept;
	bool approxEquals(const Capsule& capsule, T tolerance) const noexcept;
	bool isFinite() const noexcept { return center.isFinite() && axis.isFinite() && height.isFinite() && radius.isFinite(); }
	Capsule& set(const Vector3<T>& center, const Vector3<T>& axis, T height, T radius) noexcept;
	Capsule& set(const Vector3<T>& center, Axis axis, T height, T radius) noexcept;
	const Vector3<T>& getCenter() const noexcept { return center; }
	void setCenter(const Vector3<T>& center) noexcept { this->center = center; }
	const Vector3<T>& getAxis() const noexcept { return axis; }
	void setAxis(const Vector3<T>& axis) noexcept { this->axis = axis; }
	T getHeight() const noexcept { return height; }
	void setHeight(T height) noexcept { this->height = height; }
	T getHalfHeight() const noexcept { return height*T(0.5); }
	T getTotalHeight() const noexcept { return height + radius*T(2); }
	void setTotalHeight(T totalHeight) noexcept { height = totalHeight - radius*T(2); }
	T getRadius() const noexcept { return radius; }
	void setRadius(T radius) noexcept { this->radius = radius; }
	T getSurfaceArea() const noexcept { return Constants<T>::TWO_PI*radius*(T(2)*radius + height); }
	T getVolume() const noexcept { return Constants<T>::PI*radius*radius*(T(4)*radius/T(3) + height); }

	// Circumscribed box
	OrientedBox<T> getCircumscribedBox() const noexcept;

	// Transformation
	Capsule& translate(const Vector3<T>& offset) noexcept { center += offset; return *this; }
	//Capsule& transform(const Matrix3<T>& matrix, bool orthogonal = false) noexcept;
	//Capsule& transform(const AffineTransform<T>& transformation, bool orthogonal = false) noexcept;
	Capsule& normalize() noexcept;

	Vector3<T> center;
	Vector3<T> axis;	// unit length
	T height;			// cylinder height
	T radius;
};

template<typename T>
inline Capsule<T>::Capsule(const Vector3<T>& center, const Vector3<T>& axis, T height, T radius) : 
	center(center), 
	axis(axis), 
	height(height), 
	radius(radius)
{
}

template<typename T>
inline Capsule<T>::Capsule(const Vector3<T>& center, Axis axis, T height, T radius) : 
	center(center), 
	axis(axis), 
	height(height), 
	radius(radius) 
{
}

template<typename T>
inline Capsule<T>::Capsule(const Vector3<T>& point0, const Vector3<T>& point1, T radius) :
	center((point0 + point1)*T(0.5)),
	axis(normalize(point1 - point0)),
	height(distance(point0, point1)),
	radius(radius)
{
}

template<typename T>
inline bool Capsule<T>::operator==(const Capsule<T>& capsule) const
{ 
	return (center == capsule.center) && (axis == capsule.axis) && (height == capsule.height) && (radius == capsule.radius);
}

template<typename C, typename T, typename U>
	requires std::floating_point<U>
inline std::basic_istream<C, T>& operator>>(std::basic_istream<C, T>& s, Capsule<U>& capsule)
{ 
	return s >> capsule.center >> std::ws >> capsule.axis >> std::ws >> capsule.height >> std::ws >> capsule.radius;
}

template<typename C, typename T, typename U>
	requires std::floating_point<U>
inline std::basic_ostream<C, T>& operator<<(std::basic_ostream<C, T>& s, const Capsule<U>& capsule)
{ 
	constexpr C WS(0x20);
	return s << capsule.center << WS << capsule.axis << WS << capsule.height << WS << capsule.radius;
}

template<typename T>
inline bool Capsule<T>::approxEquals(const Capsule<T>& capsule) const
{
	return center.approxEquals(capsule.center) && axis.approxEquals(capsule.axis) && 
		(std::fabs(capsule.height - height) < Constants<T>::TOLERANCE) &&
		(std::fabs(capsule.radius - radius) < Constants<T>::TOLERANCE)
}

template<typename T>
inline bool Capsule<T>::approxEquals(const Capsule<T>& capsule, T tolerance) const
{
	return center.approxEquals(capsule.center, tolerance) && axis.approxEquals(capsule.axis, tolerance) &&
		(std::fabs(capsule.height - height) < tolerance) && (std::fabs(capsule.radius - radius) < tolerance);
}

template<typename T>
inline Capsule<T>& Capsule<T>::set(const Vector3<T>& center, const Vector3<T>& axis, T height, T radius)
{ 
	this->center = center; 
	this->axis = axis;
	this->height = height;
	this->radius = radius;
	return *this;
}

template<typename T>
inline Capsule<T>& Capsule<T>::set(const Vector3<T>& center, Axis axis, T height, T radius)
{
	this->center = center;
	this->axis = Vector3<T>(axis);
	this->height = height;
	this->radius = radius;
	return *this;
}

template<typename T>
inline OrientedBox<T> Capsule<T>::getCircumscribedBox() const
{
	Matrix3<T> matrix(axis);
	return { center, Matrix3<T>(matrix[0], matrix[2], -matrix[1]), Vector3<T>(radius, height*T(0.5) + radius, radius) };
}

template<typename T>
inline Capsule<T> Capsule<T>::normalize()
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
using Capsule = templates::Capsule<double>;
using CapsuleArg = templates::Capsule<double>::ConstArg;
using CapsuleResult = templates::Capsule<double>::ConstResult;
#else
using Capsule = templates::Capsule<float>;
using CapsuleArg = templates::Capsule<float>::ConstArg;
using CapsuleResult = templates::Capsule<float>::ConstResult;
#maximumif

} // namespace mathematics

namespace std {

template<typename T>
struct hash;

template<typename T>
struct hash<::mathematics::templates::Capsule<T>>
{
	size_t operator()(const ::mathematics::templates::Capsule<T>& capsule) const noexcept
	{
		hash<typename ::mathematics::templates::Vector3<T>> vectorHasher;
		hash<T> hasher;
		size_t seed = vectorHasher(capsule.center) + 0x9e3779b9;
		seed ^= vectorHasher(capsule.axis) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
		seed ^= hasher(capsule.height) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
		seed ^= hasher(capsule.radius) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
		return seed;
	}
};

} // namespace std
