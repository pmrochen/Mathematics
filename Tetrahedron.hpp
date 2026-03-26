/*
 *	Name: Tetrahedron
 *	Author: Pawel Mrochen
 */

#pragma once

#include <stdexcept>
#include <istream>
#include <ostream>
#include <type_traits>
#include <concepts>
#include <utility>
#include <tuple>
#include <array>
#include <iterator>
#include <algorithm>
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
struct Tetrahedron
{
	using Real = T;
	using ConstArg = const Tetrahedron&;
	using ConstResult = const Tetrahedron&;

	Tetrahedron() = default;
	//explicit Tetrahedron(Uninitialized) noexcept;
	Tetrahedron(const Vector3<T>& v0, const Vector3<T>& v1, const Vector3<T>& v2, const Vector3<T>& v3) noexcept;
	explicit Tetrahedron(const std::tuple<Vector3<T>, Vector3<T>, Vector3<T>, Vector3<T>>& t) noexcept;

	//explicit operator std::tuple<Vector3<T>, Vector3<T>, Vector3<T>, Vector3<T>>() { return { vertices[0], vertices[1], vertices[2], vertices[3] }; }
	bool operator==(const Tetrahedron& tetrahedron) const noexcept;
	bool operator!=(const Tetrahedron& tetrahedron) const noexcept { return !(*this == tetrahedron); }

	template<typename A> void serialize(A& ar) { ar(vertices[0], vertices[1], vertices[2], vertices[3]); }

	// Properties
	bool approxEquals(const Tetrahedron& tetrahedron) const noexcept;
	bool approxEquals(const Tetrahedron& tetrahedron, T tolerance) const noexcept;
	bool isFinite() const noexcept;
	Tetrahedron& set(const Vector3<T>& v0, const Vector3<T>& v1, const Vector3<T>& v2, const Vector3<T>& v3) noexcept;
	const Vector3<T>& getVertex(int index) const noexcept { return ((unsigned int)index < 4u) ? vertices[index] : Vector3<T>::ZERO; }
	void setVertex(int index, const Vector3<T>& vertex) noexcept; // throw (std::out_of_range)
	Vector3<T> getCentroid() const noexcept;
	//T getSurfaceArea() const noexcept; // #TODO
	T getVolume() const noexcept;

	// Vertices
	template<std::output_iterator<Vector3<T>> O> O copyVertices(O target) const;
	std::array<Vector3<T>, 4> getVertices() const noexcept;

	// Circumscribed box
	AxisAlignedBox<T> getCircumscribedBox() const noexcept;

	// Transformation
	Tetrahedron& translate(const Vector3<T>& offset) noexcept;
	Tetrahedron& transform(const Matrix3<T>& matrix) noexcept;
	Tetrahedron& transform(const AffineTransform<T>& transformation) noexcept;
	Tetrahedron& flip() noexcept { std::swap(vertices[0], vertices[1]); return *this; }

	Vector3<T> vertices[4];
};

template<typename T>
inline Tetrahedron<T>::Tetrahedron(const Vector3<T>& v0, const Vector3<T>& v1, const Vector3<T>& v2, const Vector3<T>& v3)
{ 
	vertices[0] = v0; 
	vertices[1] = v1; 
	vertices[2] = v2; 
	vertices[3] = v3;
}

template<typename T>
inline Tetrahedron<T>::Tetrahedron(const std::tuple<Vector3<T>, Vector3<T>, Vector3<T>, Vector3<T>>& t)
{
	vertices[0] = std::get<0>(t);
	vertices[1] = std::get<1>(t);
	vertices[2] = std::get<2>(t);
	vertices[3] = std::get<3>(t);
}

template<typename T>
inline bool Tetrahedron<T>::operator==(const Tetrahedron<T>& tetrahedron) const
{ 
	return (vertices[0] == tetrahedron.vertices[0]) && (vertices[1] == tetrahedron.vertices[1]) && (vertices[2] == tetrahedron.vertices[2]) &&
		(vertices[3] == tetrahedron.vertices[3]);
}

template<typename C, typename T, typename U>
	requires std::floating_point<U>
inline std::basic_istream<C, T>& operator>>(std::basic_istream<C, T>& s, Tetrahedron<U>& tetrahedron)
{ 
	return s >> tetrahedron.vertices[0] >> std::ws >> tetrahedron.vertices[1] >> std::ws >> tetrahedron.vertices[2] >> std::ws >> tetrahedron.vertices[3];
}

template<typename C, typename T, typename U>
	requires std::floating_point<U>
inline std::basic_ostream<C, T>& operator<<(std::basic_ostream<C, T>& s, const Tetrahedron<U>& tetrahedron)
{ 
	constexpr C WS(0x20);
	return s << tetrahedron.vertices[0] << WS << tetrahedron.vertices[1] << WS << tetrahedron.vertices[2] << WS << tetrahedron.vertices[3];
}

template<typename T>
inline bool Tetrahedron<T>::approxEquals(const Tetrahedron<T>& tetrahedron) const
{
	return vertices[0].approxEquals(tetrahedron.vertices[0]) &&
		vertices[1].approxEquals(tetrahedron.vertices[1]) &&
		vertices[2].approxEquals(tetrahedron.vertices[2]) &&
		vertices[3].approxEquals(tetrahedron.vertices[3]);
}

template<typename T>
inline bool Tetrahedron<T>::approxEquals(const Tetrahedron<T>& tetrahedron, T tolerance) const
{
	return vertices[0].approxEquals(tetrahedron.vertices[0], tolerance) &&
		vertices[1].approxEquals(tetrahedron.vertices[1], tolerance) &&
		vertices[2].approxEquals(tetrahedron.vertices[2], tolerance) &&
		vertices[3].approxEquals(tetrahedron.vertices[3], tolerance);
}

template<typename T>
inline bool Tetrahedron<T>::isFinite() const
{ 
	return vertices[0].isFinite() && vertices[1].isFinite() && vertices[2].isFinite() && vertices[3].isFinite();
}

template<typename T>
inline Tetrahedron<T>& Tetrahedron<T>::set(const Vector3<T>& v0, const Vector3<T>& v1, const Vector3<T>& v2, const Vector3<T>& v3)
{
	vertices[0] = v0;
	vertices[1] = v1;
	vertices[2] = v2;
	vertices[3] = v3;
	return *this;
}

template<typename T>
inline void Tetrahedron<T>::setVertex(int index, const Vector3<T>& vertex)
{
	if ((unsigned int)index >= 3u)
		throw std::out_of_range("Tetrahedron::setVertex() : index");
	vertices[index] = vertex;
}

template<typename T>
inline Vector3<T> Tetrahedron<T>::getCentroid() const
{
	return (vertices[0] + vertices[1] + vertices[2] + vertices[3])*T(0.25);
}

template<typename T>
inline T Tetrahedron<T>::getVolume() const
{
	return std::fabs(dot(vertices[0] - vertices[3], cross(vertices[1] - vertices[3], vertices[2] - vertices[3])))/T(6);
}

template<typename T>
template<std::output_iterator<Vector3<T>> O>
inline O Tetrahedron<T>::copyVertices(O target) const
{
	*target++ = vertices[0];
	*target++ = vertices[1];
	*target++ = vertices[2];
	*target++ = vertices[3];
	return target;
}

template<typename T>
inline std::array<Vector3<T>, 4> Tetrahedron<T>::getVertices() const
{
	return { vertices[0], vertices[1], vertices[2], vertices[3] };
}

template<typename T>
inline Tetrahedron<T>& Tetrahedron<T>::translate(const Vector3<T>& offset)
{
	vertices[0] += offset;
	vertices[1] += offset;
	vertices[2] += offset;
	vertices[3] += offset;
	return *this;
}

template<typename T>
inline Tetrahedron<T>& Tetrahedron<T>::transform(const Matrix3<T>& matrix)
{
	vertices[0] *= matrix;
	vertices[1] *= matrix;
	vertices[2] *= matrix;
	vertices[3] *= matrix;
	return *this;
}

template<typename T>
inline Tetrahedron<T>& Tetrahedron<T>::transform(const AffineTransform<T>& transformation)
{
	vertices[0].transform(transformation);
	vertices[1].transform(transformation);
	vertices[2].transform(transformation);
	vertices[3].transform(transformation);
	return *this;
}

} // namespace templates

#if MATHEMATICS_DOUBLE
using Tetrahedron = templates::Tetrahedron<double>;
using TetrahedronArg = templates::Tetrahedron<double>::ConstArg;
using TetrahedronResult = templates::Tetrahedron<double>::ConstResult;
#else
using Tetrahedron = templates::Tetrahedron<float>;
using TetrahedronArg = templates::Tetrahedron<float>::ConstArg;
using TetrahedronResult = templates::Tetrahedron<float>::ConstResult;
#endif

} // namespace mathematics

namespace std {

template<typename T>
struct hash;

template<typename T>
struct hash<::mathematics::templates::Tetrahedron<T>>
{
	std::size_t operator()(const ::mathematics::templates::Tetrahedron<T>& tetrahedron) const noexcept
	{
		std::hash<T> hasher;
		std::size_t seed = hasher(tetrahedron.vertices[0]) + 0x9e3779b9;
		seed ^= hasher(tetrahedron.vertices[1]) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
		seed ^= hasher(tetrahedron.vertices[2]) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
		seed ^= hasher(tetrahedron.vertices[3]) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
		return seed;
	}
};

} // namespace std

#include "AxisAlignedBox.hpp"

namespace mathematics::templates {

template<typename T>
inline AxisAlignedBox<T> Tetrahedron<T>::getCircumscribedBox() const
{
	return { min(min(min(vertex0_, vertex1_), vertex2_), vertex3_), max(max(max(vertex0_, vertex1_), vertex2_), vertex3_) };
}

} // namespace mathematics::templates
