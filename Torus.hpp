/*
 *	Name: Torus
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
struct Torus
{
	using Real = T;
	using ConstArg = const Torus&;
	using ConstResult = const Torus&;

	Torus() noexcept : center(), axis(), majorRadius(), minorRadius() {}
	explicit Torus(Uninitialized) noexcept : center(Uninitialized()), axis(Uninitialized()) {}
	Torus(const Vector3<T>& center, const Vector3<T>& axis, T majorRadius, T minorRadius) noexcept;
	Torus(const Vector3<T>& center, Axis axis, T majorRadius, T minorRadius) noexcept;

	bool operator==(const Torus& torus) const noexcept;
	bool operator!=(const Torus& torus) const noexcept { return !(*this == torus); }

	template<typename A> void serialize(A& ar) { ar(center, axis, majorRadius, minorRadius); }

	// Properties
	bool approxEquals(const Torus& torus) const noexcept;
	bool approxEquals(const Torus& torus, T tolerance) const noexcept;
	bool isFinite() const noexcept { return center.isFinite() && axis.isFinite() && majorRadius.isFinite() && minorRadius.isFinite(); }
	Torus& set(const Vector3<T>& center, const Vector3<T>& axis, T majorRadius, T minorRadius) noexcept;
	Torus& set(const Vector3<T>& center, Axis axis, T majorRadius, T minorRadius) noexcept;
	const Vector3<T>& getCenter() const noexcept { return center; }
	void setCenter(const Vector3<T>& center) noexcept { this->center = center; }
	const Vector3<T>& getAxis() const noexcept { return axis; }
	void setAxis(const Vector3<T>& axis) noexcept { this->axis = axis; }
	T getMajorRadius() const noexcept { return majorRadius; }
	void setMajorRadius(T majorRadius) noexcept { this->majorRadius = majorRadius; }
	T getMinorRadius() const noexcept { return minorRadius; }
	void setMinorRadius(T minorRadius) noexcept { this->minorRadius = minorRadius; }
	T getInnerRadius() const noexcept { return majorRadius - minorRadius; }
	T getOuterRadius() const noexcept { return majorRadius + minorRadius; }
	T getAspectRatio() const noexcept { return majorRadius/minorRadius; }
	T getSurfaceArea() const noexcept { return T(4)*Constants<T>::PI_SQR*majorRadius*minorRadius;  }
	T getVolume() const noexcept { return T(2)*Constants<T>::PI_SQR*majorRadius*minorRadius*minorRadius;  }

	// Circumscribed box
	OrientedBox<T> getCircumscribedBox() const noexcept;

	// Transformation
	Torus& translate(const Vector3<T>& offset) noexcept { center += offset; return *this; }
	//Torus& transform(const Matrix3<T>& matrix, bool orthogonal = false) noexcept;
	//Torus& transform(const AffineTransform<T>& transformation, bool orthogonal = false) noexcept;

	Vector3<T> center;
	Vector3<T> axis;	// unit length
	T majorRadius;		// distance from the center of the tube to the center of the torus
	T minorRadius;		// radius of the tube
};

template<typename T>
inline Torus<T>::Torus(const Vector3<T>& center, const Vector3<T>& axis, T majorRadius, T minorRadius) : 
	center(center), 
	axis(axis), 
	majorRadius(majorRadius), 
	minorRadius(minorRadius)
{
}

template<typename T>
inline Torus<T>::Torus(const Vector3<T>& center, Axis axis, T majorRadius, T minorRadius) : 
	center(center), 
	axis(axis), 
	majorRadius(majorRadius), 
	minorRadius(minorRadius) 
{
}

template<typename T>
inline bool Torus<T>::operator==(const Torus<T>& torus) const
{ 
	return (center == torus.center) && (axis == torus.axis) && (majorRadius == torus.majorRadius) && (minorRadius == torus.minorRadius);
}

template<typename C, typename T, typename U>
	requires std::floating_point<U>
inline std::basic_istream<C, T>& operator>>(std::basic_istream<C, T>& s, Torus<U>& torus)
{ 
	return s >> torus.center >> std::ws >> torus.axis >> std::ws >> torus.majorRadius >> std::ws >> torus.minorRadius;
}

template<typename C, typename T, typename U>
	requires std::floating_point<U>
inline std::basic_ostream<C, T>& operator<<(std::basic_ostream<C, T>& s, const Torus<U>& torus)
{ 
	constexpr C WS(0x20);
	return s << torus.center << WS << torus.axis << WS << torus.majorRadius << WS << torus.minorRadius;
}

template<typename T>
inline bool Torus<T>::approxEquals(const Torus<T>& torus) const
{
	return center.approxEquals(torus.center) && axis.approxEquals(torus.axis) && 
		(std::fabs(torus.majorRadius - majorRadius) < Constants<T>::TOLERANCE) &&
		(std::fabs(torus.minorRadius - minorRadius) < Constants<T>::TOLERANCE)
}

template<typename T>
inline bool Torus<T>::approxEquals(const Torus<T>& torus, T tolerance) const
{
	return center.approxEquals(torus.center, tolerance) && axis.approxEquals(torus.axis, tolerance) &&
		(std::fabs(torus.majorRadius - majorRadius) < tolerance) && (std::fabs(torus.minorRadius - minorRadius) < tolerance);
}

template<typename T>
inline Torus<T>& Torus<T>::set(const Vector3<T>& center, const Vector3<T>& axis, T majorRadius, T minorRadius)
{ 
	this->center = center; 
	this->axis = axis;
	this->majorRadius = majorRadius;
	this->minorRadius = minorRadius;
	return *this;
}

template<typename T>
inline Torus<T>& Torus<T>::set(const Vector3<T>& center, Axis axis, T majorRadius, T minorRadius)
{
	this->center = center;
	this->axis = Vector3<T>(axis);
	this->majorRadius = majorRadius;
	this->minorRadius = minorRadius;
	return *this;
}

template<typename T>
inline OrientedBox<T> Torus<T>::getCircumscribedBox() const
{
	Matrix3<T> matrix(axis);
	return { center, Matrix3<T>(matrix[0], matrix[2], -matrix[1]), 
		Vector3<T>(majorRadius + minorRadius, minorRadius, majorRadius + minorRadius) };
}

} // namespace templates

#if MATHEMATICS_DOUBLE
using Torus = templates::Torus<double>;
using TorusArg = templates::Torus<double>::ConstArg;
using TorusResult = templates::Torus<double>::ConstResult;
#else
using Torus = templates::Torus<float>;
using TorusArg = templates::Torus<float>::ConstArg;
using TorusResult = templates::Torus<float>::ConstResult;
#endif

} // namespace mathematics

namespace std {

template<typename T>
struct hash;

template<typename T>
struct hash<::mathematics::templates::Torus<T>>
{
	size_t operator()(const ::mathematics::templates::Torus<T>& torus) const noexcept
	{
		hash<typename ::mathematics::templates::Vector3<T>> vectorHasher;
		hash<T> hasher;
		size_t seed = vectorHasher(torus.center) + 0x9e3779b9;
		seed ^= vectorHasher(torus.axis) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
		seed ^= hasher(torus.majorRadius) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
		seed ^= hasher(torus.minorRadius) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
		return seed;
	}
};

} // namespace std
