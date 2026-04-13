/*
 *	Name: Quadrilateral3
 *	Author: Pawel Mrochen
 */

#pragma once

#include <stdexcept>
#include <istream>
#include <ostream>
#include <type_traits>
#include <concepts>
#include <algorithm>
#include <functional>
#include <utility>
#include <tuple>
#include <iterator>
#include <cstddef>
#include <cmath>
#include "Vector3.hpp"
#include "Matrix3.hpp"
#include "AffineTransform.hpp"

namespace mathematics {
namespace templates {

template<typename T>
	requires std::floating_point<T>
struct AxisAlignedBox;

template<typename T>
	requires std::floating_point<T>
struct Quadrilateral3
{
	using Real = T;
	using ConstArg = const Quadrilateral3&;
	using ConstResult = const Quadrilateral3&;

	Quadrilateral3() = default;
	//explicit Quadrilateral3(Uninitialized) noexcept;
	Quadrilateral3(const Vector3<T>& v0, const Vector3<T>& v1, const Vector3<T>& v2, const Vector3<T>& v3) noexcept;
	explicit Quadrilateral3(const std::tuple<Vector3<T>, Vector3<T>, Vector3<T>, Vector3<T>>& t) noexcept;

	//explicit operator std::tuple<Vector3<T>, Vector3<T>, Vector3<T>, Vector3<T>>() { return { vertices[0], vertices[1], vertices[2], vertices[3] }; }
	bool operator==(const Quadrilateral3& quad) const noexcept;
	bool operator!=(const Quadrilateral3& quad) const noexcept { return !(*this == quad); }

	template<typename A> void serialize(A& ar) { ar(vertices[0], vertices[1], vertices[2], vertices[3]); }

	// Properties
	bool approxEquals(const Quadrilateral3& quad) const noexcept;
	bool approxEquals(const Quadrilateral3& quad, T tolerance) const noexcept;
	bool isFinite() const noexcept;
	Quadrilateral3& set(const Vector3<T>& v0, const Vector3<T>& v1, const Vector3<T>& v2, const Vector3<T>& v3) noexcept;
	const Vector3<T>& getVertex(int index) const noexcept { return ((unsigned int)index < 4u) ? vertices[index] : Vector3<T>::ZERO; }
	void setVertex(int index, const Vector3<T>& vertex) noexcept; // throw (std::out_of_range)
	T getPerimeter() const noexcept;
	T getArea() const noexcept;

	// Vertices
	template<std::output_iterator<Vector3<T>> O> O copyVertices(O target) const;
	std::tuple<Vector3<T>, Vector3<T>, Vector3<T>, Vector3<T>> getVertices() const noexcept;

	// Normal
	static Vector3<T> computeNormal(const Vector3<T>& v0, const Vector3<T>& v1, const Vector3<T>& v2, const Vector3<T>& v3) noexcept;
	Vector3<T> getNormal() const noexcept;

	// Circumscribed box
	AxisAlignedBox<T> getCircumscribedBox() const noexcept;

	// Triangulation
	//std::pair<Triangle3<T>, Triangle3<T>> triangulate() noexcept; // #TODO

	// Transformation
	Quadrilateral3& translate(const Vector3<T>& offset) noexcept;
	Quadrilateral3& transform(const Matrix3<T>& matrix) noexcept;
	Quadrilateral3& transform(const AffineTransform<T>& transformation) noexcept;
	Quadrilateral3& flip() noexcept;

	Vector3<T> vertices[4];
};

template<typename T>
inline Quadrilateral3<T>::Quadrilateral3(const Vector3<T>& v0, const Vector3<T>& v1, const Vector3<T>& v2, const Vector3<T>& v3)
{ 
	vertices[0] = v0; 
	vertices[1] = v1; 
	vertices[2] = v2; 
	vertices[3] = v3;
}

template<typename T>
inline Quadrilateral3<T>::Quadrilateral3(const std::tuple<Vector3<T>, Vector3<T>, Vector3<T>, Vector3<T>>& t)
{
	vertices[0] = std::get<0>(t);
	vertices[1] = std::get<1>(t);
	vertices[2] = std::get<2>(t);
	vertices[3] = std::get<3>(t);
}

template<typename T>
inline bool Quadrilateral3<T>::operator==(const Quadrilateral3<T>& quad) const
{ 
	return (vertices[0] == quad.vertices[0]) && (vertices[1] == quad.vertices[1]) && (vertices[2] == quad.vertices[2]) &&
		(vertices[3] == quad.vertices[3]);
}

template<typename C, typename T, typename U>
	requires std::floating_point<U>
inline std::basic_istream<C, T>& operator>>(std::basic_istream<C, T>& s, Quadrilateral3<U>& quad)
{ 
	return s >> quad.vertices[0] >> std::ws >> quad.vertices[1] >> std::ws >> quad.vertices[2] >> std::ws >> quad.vertices[3];
}

template<typename C, typename T, typename U>
	requires std::floating_point<U>
inline std::basic_ostream<C, T>& operator<<(std::basic_ostream<C, T>& s, const Quadrilateral3<U>& quad)
{ 
	constexpr C WS(0x20);
	return s << quad.vertices[0] << WS << quad.vertices[1] << WS << quad.vertices[2] << WS << quad.vertices[3];
}

template<typename T>
inline bool Quadrilateral3<T>::approxEquals(const Quadrilateral3<T>& quad) const
{
	return vertices[0].approxEquals(quad.vertices[0]) &&
		vertices[1].approxEquals(quad.vertices[1]) &&
		vertices[2].approxEquals(quad.vertices[2]) &&
		vertices[3].approxEquals(quad.vertices[3]);
}

template<typename T>
inline bool Quadrilateral3<T>::approxEquals(const Quadrilateral3<T>& quad, T tolerance) const
{
	return vertices[0].approxEquals(quad.vertices[0], tolerance) &&
		vertices[1].approxEquals(quad.vertices[1], tolerance) &&
		vertices[2].approxEquals(quad.vertices[2], tolerance) &&
		vertices[3].approxEquals(quad.vertices[3], tolerance);
}

template<typename T>
inline bool Quadrilateral3<T>::isFinite() const
{ 
	return vertices[0].isFinite() && vertices[1].isFinite() && vertices[2].isFinite() && vertices[3].isFinite();
}

template<typename T>
inline Quadrilateral3<T>& Quadrilateral3<T>::set(const Vector3<T>& v0, const Vector3<T>& v1, const Vector3<T>& v2, const Vector3<T>& v3)
{
	vertices[0] = v0;
	vertices[1] = v1;
	vertices[2] = v2;
	vertices[3] = v3;
	return *this;
}

template<typename T>
inline void Quadrilateral3<T>::setVertex(int index, const Vector3<T>& vertex)
{
	if ((unsigned int)index >= 4u)
		throw std::out_of_range("Quadrilateral3::setVertex() : index");
	vertices[index] = vertex;
}

template<typename T>
inline T Quadrilateral3<T>::getPerimeter() const
{
	return distance(vertices[0], vertices[1]) + distance(vertices[1], vertices[2]) + distance(vertices[2], vertices[3]) +
		distance(vertices[3], vertices[0]);
}

template<typename T>
inline T Quadrilateral3<T>::getArea() const
{
	return cross(vertices[2] - vertices[0], vertices[3] - vertices[1]).getMagnitude()*T(0.5);
}

template<typename T>
template<std::output_iterator<Vector3<T>> O>
inline O Quadrilateral3<T>::copyVertices(O target) const
{
	*target++ = vertices[0];
	*target++ = vertices[1];
	*target++ = vertices[2];
	*target++ = vertices[3];
	return target;
}

template<typename T>
inline std::tuple<Vector3<T>, Vector3<T>, Vector3<T>, Vector3<T>> Quadrilateral3<T>::getVertices() const
{ 
	return { vertices[0], vertices[1], vertices[2], vertices[3] }; 
}

template<typename T>
/*static*/ inline Vector3<T> Quadrilateral3<T>::computeNormal(const Vector3<T>& v0, const Vector3<T>& v1, const Vector3<T>& v2, 
	const Vector3<T>& v3)
{
	return normalize(cross(v1 - v0, v2 - v0) + cross(v2 - v0, v3 - v0));
}

template<typename T>
inline Vector3<T> Quadrilateral3<T>::getNormal() const
{
	return normalize(cross(vertices[1] - vertices[0], vertices[2] - vertices[0]) + 
		cross(vertices[2] - vertices[0], vertices[3] - vertices[0]));
}

template<typename T>
inline Quadrilateral3<T>& Quadrilateral3<T>::translate(const Vector3<T>& offset)
{
	vertices[0] += offset;
	vertices[1] += offset;
	vertices[2] += offset;
	vertices[3] += offset;
	return *this;
}

template<typename T>
inline Quadrilateral3<T>& Quadrilateral3<T>::transform(const Matrix3<T>& matrix)
{
	vertices[0] *= matrix;
	vertices[1] *= matrix;
	vertices[2] *= matrix;
	vertices[3] *= matrix;
	return *this;
}

template<typename T>
inline Quadrilateral3<T>& Quadrilateral3<T>::transform(const AffineTransform<T>& transformation)
{
	vertices[0].transform(transformation);
	vertices[1].transform(transformation);
	vertices[2].transform(transformation);
	vertices[3].transform(transformation);
	return *this;
}

template<typename T>
inline Quadrilateral3<T>& Quadrilateral3<T>::flip()
{
	std::swap(vertices[0], vertices[3]);
	std::swap(vertices[1], vertices[2]);
	return *this;
}

} // namespace templates

#if MATHEMATICS_DOUBLE
using Quadrilateral3 = templates::Quadrilateral3<double>;
using Quadrilateral3Arg = templates::Quadrilateral3<double>::ConstArg;
using Quadrilateral3Result = templates::Quadrilateral3<double>::ConstResult;
#else
using Quadrilateral3 = templates::Quadrilateral3<float>;
using Quadrilateral3Arg = templates::Quadrilateral3<float>::ConstArg;
using Quadrilateral3Result = templates::Quadrilateral3<float>::ConstResult;
#endif

} // namespace mathematics

namespace std {

template<typename T>
struct hash;

template<typename T>
struct hash<::mathematics::templates::Quadrilateral3<T>>
{
	size_t operator()(const ::mathematics::templates::Quadrilateral3<T>& quad) const noexcept
	{
		hash<typename ::mathematics::templates::Vector3<T>> hasher;
		size_t seed = hasher(quad.vertices[0]) + 0x9e3779b9;
		seed ^= hasher(quad.vertices[1]) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
		seed ^= hasher(quad.vertices[2]) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
		seed ^= hasher(quad.vertices[3]) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
		return seed;
	}
};

} // namespace std

#include "AxisAlignedBox.hpp"

namespace mathematics::templates {

template<typename T>
inline AxisAlignedBox<T> Quadrilateral3<T>::getCircumscribedBox() const
{
	return { min(min(min(vertex0_, vertex1_), vertex2_), vertex3_), max(max(max(vertex0_, vertex1_), vertex2_), vertex3_) };
}

} // namespace mathematics::templates
