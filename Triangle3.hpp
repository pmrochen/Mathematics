/*
 *	Name: Triangle3
 *	Author: Pawel Mrochen
 */

#pragma once

#include <stdexcept>
#include <istream>
#include <ostream>
#include <limits>
#include <type_traits>
#include <concepts>
#include <algorithm>
#include <functional>
#include <utility>
#include <tuple>
#include <vector>
#include <iterator>
#include <cstddef>
#include <cmath>
#include <malloc.h>
#include "Vector3.hpp"
#include "Matrix3.hpp"
#include "AffineTransform.hpp"
#include "HalfSpace.hpp"
#include "Plane.hpp"

namespace mathematics {
namespace templates {

template<typename T>
	requires std::floating_point<T>
struct AxisAlignedBox;

template<typename T>
	requires std::floating_point<T>
struct OrientedBox;

template<typename T>
	requires std::floating_point<T>
struct Sphere;

template<typename T>
	requires std::floating_point<T>
struct Triangle3
{
	using Real = T;
	using ConstArg = const Triangle3&;
	using ConstResult = const Triangle3&;

	Triangle3() = default;
	//explicit Triangle3(Uninitialized) noexcept;
	Triangle3(const Vector3<T>& v0, const Vector3<T>& v1, const Vector3<T>& v2) noexcept;
	explicit Triangle3(const std::tuple<Vector3<T>, Vector3<T>, Vector3<T>>& t) noexcept;

	//explicit operator std::tuple<Vector3<T>, Vector3<T>, Vector3<T>>() { return { vertices[0], vertices[1], vertices[2] }; }
	//Vector3<T> operator()(T u, T v) const noexcept { return evaluate(u, v); }
	bool operator==(const Triangle3& triangle) const noexcept;
	bool operator!=(const Triangle3& triangle) const noexcept { return !(*this == triangle); }

	template<typename A> void serialize(A& ar) { ar(vertices[0], vertices[1], vertices[2]); }

	// Properties
	bool approxEquals(const Triangle3& triangle) const noexcept;
	bool approxEquals(const Triangle3& triangle, T tolerance) const noexcept;
	bool isFinite() const noexcept { return vertices[0].isFinite() && vertices[1].isFinite() && vertices[2].isFinite(); }
	Triangle3& set(const Vector3<T>& v0, const Vector3<T>& v1, const Vector3<T>& v2) noexcept;
	const Vector3<T>& getVertex(int index) const noexcept { return ((unsigned int)index < 3u) ? vertices[index] : Vector3<T>::ZERO; }
	void setVertex(int index, const Vector3<T>& vertex) noexcept; // throw (std::out_of_range)
	T getPerimeter() const noexcept;
	T getArea() const noexcept;

	// Vertices
	template<std::output_iterator<Vector3<T>> O> O copyVertices(O target) const;
	std::tuple<Vector3<T>, Vector3<T>, Vector3<T>> getVertices() const noexcept;

	// Normal
	static Vector3<T> computeNormal(const Vector3<T>& v0, const Vector3<T>& v1, const Vector3<T>& v2) noexcept;
	Vector3<T> getNormal() const noexcept;

	// Tangent space basis
	static Matrix3<T> computeTangentSpaceBasis(const Vector3<T>& v0, const Vector3<T>& v1, const Vector3<T>& v2, 
		const Vector2<T>& uv0, const Vector2<T>& uv1, const Vector2<T>& uv2, bool weightedByArea) noexcept;
	Matrix3<T> getTangentSpaceBasis(const Vector2<T>& uv0, const Vector2<T>& uv1, const Vector2<T>& uv2, bool weightedByArea) noexcept;

	// Circles
	template<typename U = /*Circle3<T>*/std::pair<Vector3<T>, T>> U getInscribedCircle() const noexcept;
	template<typename U = /*Circle3<T>*/std::pair<Vector3<T>, T>> U getCircumscribedCircle() const noexcept;

	// Circumscribed box
	AxisAlignedBox<T> getCircumscribedBox() const noexcept;

	// Circumscribed sphere
	Sphere<T> getCircumscribedSphere() const noexcept;

	// Transformation
	Triangle3& translate(const Vector3<T>& offset) noexcept;
	Triangle3& transform(const Matrix3<T>& matrix) noexcept;
	Triangle3& transform(const AffineTransform<T>& transformation) noexcept;
	Triangle3& flip() noexcept { std::swap(vertices[0], vertices[2]); }

	// Split (returns 1 or 3 triangles)
	//std::vector<Triangle3<T>> split(const Plane<T>& plane); // #TODO

	// Triangulation (returns (vertices - 2)*3 indices)
	template<std::random_access_iterator<Vector3<T>> I, std::integral U, std::output_iterator<U> O>
	static std::pair<O, bool> triangulate(I firstVertex, I lastVertex, O outIndex);
	template<std::random_access_iterator<Vector3<T>> I, std::integral U, std::random_access_iterator<U> J, std::output_iterator<U> O>
	static std::pair<O, bool> triangulate(I firstVertex, I lastVertex, J firstIndex, J lastIndex, O outIndex);

	// Evaluation (u + v + w = 1)
	Vector3<T> evaluate(T v, T w) const noexcept { return vertices[0]*(T(1) - v - w) + vertices[1]*v + vertices[2]*w; }
	Vector3<T> evaluate(T u, T v, T w) const noexcept { return vertices[0]*u + vertices[1]*v + vertices[2]*w; }
	Vector3<T> getBarycentricCoords(const Vector3<T>& point) const noexcept;

	// Closest point
	Vector3<T> getClosestPoint(const Vector3<T>& point) const noexcept;
	T getDistanceTo(const Vector3<T>& point) const noexcept;

	// Intersection
	bool intersects(const HalfSpace<T>& halfSpace) const noexcept;
	bool intersects(const Plane<T>& plane) const noexcept;
	bool intersects(const AxisAlignedBox& box) const noexcept;
	bool intersects(const OrientedBox<T>& box) const noexcept;
	bool intersects(const Sphere<T>& sphere) const noexcept;

	Vector3<T> vertices[3];
};

template<typename T>
inline Triangle3<T>::Triangle3(const Vector3<T>& v0, const Vector3<T>& v1, const Vector3<T>& v2)
{ 
	vertices[0] = v0; 
	vertices[1] = v1; 
	vertices[2] = v2; 
}

template<typename T>
inline Triangle3<T>::Triangle3(const std::tuple<Vector3<T>, Vector3<T>, Vector3<T>>& t)
{
	vertices[0] = std::get<0>(t);
	vertices[1] = std::get<1>(t);
	vertices[2] = std::get<2>(t);
}

template<typename T>
inline bool Triangle3<T>::operator==(const Triangle3<T>& triangle) const
{ 
	return (vertices[0] == triangle.vertices[0]) && (vertices[1] == triangle.vertices[1]) && (vertices[2] == triangle.vertices[2]);
}

template<typename C, typename T, typename U>
	requires std::floating_point<U>
inline std::basic_istream<C, T>& operator>>(std::basic_istream<C, T>& s, Triangle3<U>& triangle)
{ 
	return s >> triangle.vertices[0] >> std::ws >> triangle.vertices[1] >> std::ws >> triangle.vertices[2];
}

template<typename C, typename T, typename U>
	requires std::floating_point<U>
inline std::basic_ostream<C, T>& operator<<(std::basic_ostream<C, T>& s, const Triangle3<U>& triangle)
{ 
	constexpr C WS(0x20);
	return s << triangle.vertices[0] << WS << triangle.vertices[1] << WS << triangle.vertices[2];
}

template<typename T>
inline bool Triangle3<T>::approxEquals(const Triangle3<T>& triangle) const
{
	return vertices[0].approxEquals(triangle.vertices[0]) &&
		vertices[1].approxEquals(triangle.vertices[1]) &&
		vertices[2].approxEquals(triangle.vertices[2]);
}

template<typename T>
inline bool Triangle3<T>::approxEquals(const Triangle3<T>& triangle, T tolerance) const
{
	return vertices[0].approxEquals(triangle.vertices[0], tolerance) &&
		vertices[1].approxEquals(triangle.vertices[1], tolerance) &&
		vertices[2].approxEquals(triangle.vertices[2], tolerance);
}

template<typename T>
inline Triangle3<T>& Triangle3<T>::set(const Vector3<T>& v0, const Vector3<T>& v1, const Vector3<T>& v2)
{
	vertices[0] = v0;
	vertices[1] = v1;
	vertices[2] = v2;
	return *this;
}

template<typename T>
inline void Triangle3<T>::setVertex(int index, const Vector3<T>& vertex)
{
	if ((unsigned int)index >= 3u)
		throw std::out_of_range("Triangle3::setVertex() : index");
	vertices[index] = vertex;
}

template<typename T>
inline T Triangle3<T>::getPerimeter() const
{
	return distance(vertices[0], vertices[1]) + distance(vertices[1], vertices[2]) + distance(vertices[2], vertices[0]);
}

template<typename T>
inline T Triangle3<T>::getArea() const
{
	return (cross(vertices[0], vertices[1]) + cross(vertices[1], vertices[2]) + cross(vertices[2], vertices[0])).getMagnitude()*T(0.5);
}

template<typename T>
template<std::output_iterator<Vector3<T>> O>
inline O Triangle3<T>::copyVertices(O target) const
{
	*target++ = vertices[0];
	*target++ = vertices[1];
	*target++ = vertices[2];
	return target;
}

template<typename T>
inline std::tuple<Vector3<T>, Vector3<T>, Vector3<T>> Triangle3<T>::getVertices() const
{ 
	return { vertices[0], vertices[1], vertices[2] }; 
}

template<typename T>
/*static*/ inline Vector3<T> Triangle3<T>::computeNormal(const Vector3<T>& v0, const Vector3<T>& v1, const Vector3<T>& v2)
{
	return normalize(cross(v1 - v0, v2 - v0));
}

template<typename T>
inline Vector3<T> Triangle3<T>::getNormal() const
{
	return normalize(cross(vertices[1] - vertices[0], vertices[2] - vertices[0]));
}

template<typename T>
/*static*/ inline Matrix3<T> Triangle3<T>::computeTangentSpaceBasis(const Vector3<T>& v0, const Vector3<T>& v1, const Vector3<T>& v2,
	const Vector2<T>& uv0, const Vector2<T>& uv1, const Vector2<T>& uv2, bool weightedByArea)
{
	Vector3<T> e1 = v1 - v0;
	Vector3<T> e2 = v2 - v0;
	Vector3<T> normal = normalize(cross(e1, e2));

	Vector2<T> delta1 = uv1 - uv0;
	Vector2<T> delta2 = uv2 - uv0;
	T d = cross(delta1, delta2);

	Vector3<T> tangent(Uninitialized());
	Vector3<T> bitangent(Uninitialized());
	if (/*!_isnan(d) &&*/ (d != T(0)))
	{
		tangent = normalize(e1*(delta2.y/d) + e2*(-delta1.y/d));
		bitangent = normalize(e1*(-delta2.x/d) + e2*(delta1.x/d));
		if (weightedByArea)
		{
			T area = std::fabs(d)*T(0.5);
			tangent *= area; // tangent and bitangent become weighted by the uv triangle size
			bitangent *= area;
		}
	}
	else
	{
		tangent.setZero();
		bitangent.setZero();
	}

	return { tangent, bitangent, normal };
}

template<typename T>
inline Matrix3<T> Triangle3<T>::getTangentSpaceBasis(const Vector2<T>& uv0, const Vector2<T>& uv1, const Vector2<T>& uv2,
	bool weightedByArea)
{
	return computeTangentSpaceBasis(vertices[0], vertices[1], vertices[2], uv0, uv1, uv2, weightedByArea);
}

template<typename T>
template<typename U>
inline U Triangle3<T>::getInscribedCircle() const
{
	T area = getArea();
	if (area > std::numeric_limits<T>::epsilon())
	{
		T invPerim = T(1)/getPerimeter();
		return { (vertices[0]*distance(vertices[1], vertices[2]) + vertices[1]*distance(vertices[2], vertices[0]) +
			vertices[2]*distance(vertices[0], vertices[1]))*invPerim, T(2)*area*invPerim };
	}
	else
	{
		return { vertices[0], T(0) };
	}
}

template<typename T>
template<typename U>
inline U Triangle3<T>::getCircumscribedCircle() const
{
	T d1 = dot(vertices[2] - vertices[0], vertices[1] - vertices[0]);
	T d2 = dot(vertices[2] - vertices[1], vertices[0] - vertices[1]);
	T d3 = dot(vertices[0] - vertices[2], vertices[1] - vertices[2]);
	T c1 = d2*d3;
	T c2 = d3*d1;
	T c3 = d1*d2;
	T invC = T(1)/(c1 + c2 + c3);
	return { (vertices[0]*(c2 + c3) + vertices[1]*(c3 + c1) + vertices[2]*(c1 + c2))*invC*T(0.5), 
		T(0.5)*std::sqrt((d1 + d2)*(d2 + d3)*(d3 + d1)*invC) };
}

template<typename T>
inline Triangle3<T>& Triangle3<T>::translate(const Vector3<T>& offset)
{
	vertices[0] += offset;
	vertices[1] += offset;
	vertices[2] += offset;
	return *this;
}

template<typename T>
inline Triangle3<T>& Triangle3<T>::transform(const Matrix3<T>& matrix)
{
	vertices[0] *= matrix;
	vertices[1] *= matrix;
	vertices[2] *= matrix;
	return *this;
}

template<typename T>
inline Triangle3<T>& Triangle3<T>::transform(const AffineTransform<T>& transformation)
{
	vertices[0].transform(transformation);
	vertices[1].transform(transformation);
	vertices[2].transform(transformation);
	return *this;
}

template<typename T>
inline Vector3<T> Triangle3<T>::getBarycentricCoords(const Vector3<T>& point) const
{
	Vector3<T> edge1 = vertices[1] - vertices[0];
	Vector3<T> edge2 = vertices[2] - vertices[0];
	T m = T(1)/cross(edge1, edge2).getMagnitudeSquared();
	T v = std::sqrt(cross(edge2, point - vertices[2]).getMagnitudeSquared()*m);
	T w = std::sqrt(cross(edge1, point - vertices[1]).getMagnitudeSquared()*m);
	return Vector3<T>(T(1) - v - w, v, w);
}

template<typename T>
inline bool Triangle3<T>::intersects(const HalfSpace<T>& halfSpace) const
{
	return halfSpace.contains(vertices[0]) || halfSpace.contains(vertices[1]) || halfSpace.contains(vertices[2]);
}

//template<typename T>
//inline bool Triangle3<T>::intersects(const Plane<T>& plane) const // #TODO Check if all vertices are on one side
//{
//	T d0 = plane.getSignedDistanceTo(vertices[0]);
//	T d1 = plane.getSignedDistanceTo(vertices[1]);
//	T d2 = plane.getSignedDistanceTo(vertices[2]);
//	return !(((d0 < T(0)) && (d1 < T(0)) && (d2 < T(0))) || ((d0 > T(0)) && (d1 > T(0)) && (d2 > T(0))));
//}

} // namespace templates

#if MATHEMATICS_DOUBLE
using Triangle3 = templates::Triangle3<double>;
using Triangle3Arg = templates::Triangle3<double>::ConstArg;
using Triangle3Result = templates::Triangle3<double>::ConstResult;
#else
using Triangle3 = templates::Triangle3<float>;
using Triangle3Arg = templates::Triangle3<float>::ConstArg;
using Triangle3Result = templates::Triangle3<float>::ConstResult;
#endif

} // namespace mathematics

namespace std {

template<typename T>
struct hash;

template<typename T>
struct hash<::mathematics::templates::Triangle3<T>>
{
	size_t operator()(const ::mathematics::templates::Triangle3<T>& triangle) const noexcept
	{
		hash<typename ::mathematics::templates::Vector3<T>> hasher;
		size_t seed = hasher(triangle.vertices[0]) + 0x9e3779b9;
		seed ^= hasher(triangle.vertices[1]) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
		seed ^= hasher(triangle.vertices[2]) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
		return seed;
	}
};

} // namespace std

#include "AxisAlignedBox.hpp"
#include "OrientedBox.hpp"
#include "Sphere.hpp"
#include "Distances.inl"
#include "Intersections.inl"
#include "Triangulation.inl"

namespace mathematics::templates {

template<typename T>
inline AxisAlignedBox<T> Triangle3<T>::getCircumscribedBox() const
{
	return { min(min(vertex0_, vertex1_), vertex2_), max(max(vertex0_, vertex1_), vertex2_) };
}

template<typename T>
inline Sphere<T> Triangle3<T>::getCircumscribedSphere() const
{ 
	return getCircumscribedCircle<Sphere<T>>(); 
}

template<typename T>
template<std::random_access_iterator<Vector3<T>> I, std::integral U, std::output_iterator<U> O>
/*static*/ std::pair<O, bool> Triangle3<T>::triangulate(I firstVertex, I lastVertex, O outIndex)
{
	//return triangulation::triangulate3<T, I, U, O>(firstVertex, lastVertex, outIndex);

	std::ptrdiff_t nVertices = std::distance(firstVertex, lastVertex);
	if (nVertices <= 3)
	{
		for (std::ptrdiff_t i = 0; i < nVertices; i++)
			*outIndex++ = U(i);
		
		return { outIndex, false };
	}

	if (nVertices == 4) // Split along the shortest diagonal
	{
		if (distanceSquared(firstVertex[0], firstVertex[2]) <= distanceSquared(firstVertex[1], firstVertex[3]))
		{
			*outIndex++ = U(0);
			*outIndex++ = U(1);
			*outIndex++ = U(2);
			*outIndex++ = U(0);
			*outIndex++ = U(2);
			*outIndex++ = U(3);
		}
		else
		{
			*outIndex++ = U(0);
			*outIndex++ = U(1);
			*outIndex++ = U(3);
			*outIndex++ = U(1);
			*outIndex++ = U(2);
			*outIndex++ = U(3);
		}

		return { outIndex, true };
	}

	Vector3<T> e1 = firstVertex[1] - firstVertex[0];
	Vector3<T> e2 = firstVertex[nVertices - 1] - firstVertex[0];
	Matrix3<T> basis(Uninitialized());
	basis[0] = normalize(e1);
	basis[2] = normalize(cross(e1, e2));
	basis[1] = normalize(cross(basis[2], basis[0]));
	basis.transpose();

	Vector2<T>* vertices2 = (Vector2<T>*)alloca(nVertices*sizeof(Vector2<T>));
	for (std::ptrdiff_t i = 0; i < nVertices; i++)
		vertices2[i] = (firstVertex[i]*basis).xy();

	return triangulation::triangulate2(vertices2, vertices2 + nVertices, outIndex);
}

template<typename T>
template<std::random_access_iterator<Vector3<T>> I, std::integral U, std::random_access_iterator<U> J, std::output_iterator<U> O>
/*static*/ std::pair<O, bool> Triangle3<T>::triangulate(I firstVertex, I lastVertex, J firstIndex, J lastIndex, O outIndex)
{
	//return triangulation::triangulate3<T, I, U, J, O>(firstVertex, lastVertex, firstIndex, lastIndex, outIndex);

	std::ptrdiff_t nIndices = std::distance(firstIndex, lastIndex);
	if (nIndices <= 3)
	{
		for (std::ptrdiff_t i = 0; i < nIndices; i++)
			*outIndex++ = U(*firstIndex++);

		return { outIndex, false };
	}

	std::ptrdiff_t nVertices = std::distance(firstVertex, lastVertex);
	for (std::ptrdiff_t i = 0; i < nIndices; i++)
	{
		if (std::ptrdiff_t(std::size_t(firstIndex[i])) >= nVertices) // Index out of range
		{
			//throw std::runtime_error("triangulate() : index out of range");

			for (std::ptrdiff_t j = 0, n = nIndices - 2; j < n; j++)
			{
				*outIndex++ = U(firstIndex[j]);
				*outIndex++ = U(firstIndex[j + 1]);
				*outIndex++ = U(firstIndex[j + 2]);
			}

			return { outIndex, false };
		}
	}

	if (nVertices == 4) // Split along the shortest diagonal
	{
		if (distanceSquared(firstVertex[firstIndex[0]], firstVertex[firstIndex[2]]) <= distanceSquared(firstVertex[firstIndex[1]], firstVertex[firstIndex[3]]))
		{
			*outIndex++ = firstIndex[0];
			*outIndex++ = firstIndex[1];
			*outIndex++ = firstIndex[2];
			*outIndex++ = firstIndex[0];
			*outIndex++ = firstIndex[2];
			*outIndex++ = firstIndex[3];
		}
		else
		{
			*outIndex++ = firstIndex[0];
			*outIndex++ = firstIndex[1];
			*outIndex++ = firstIndex[3];
			*outIndex++ = firstIndex[1];
			*outIndex++ = firstIndex[2];
			*outIndex++ = firstIndex[3];
		}

		return { outIndex, true };
	}

	Vector3<T> e1 = firstVertex[firstIndex[1]] - firstVertex[firstIndex[0]];
	Vector3<T> e2 = firstVertex[firstIndex[nIndices - 1]] - firstVertex[firstIndex[0]];
	Matrix3<T> basis(Uninitialized());
	basis[0] = normalize(e1);
	basis[2] = normalize(cross(e1, e2));
	basis[1] = normalize(cross(basis[2], basis[0]));
	basis.transpose();

	Vector2<T>* vertices2 = (Vector2<T>*)alloca(nIndices*sizeof(Vector2<T>));
	for (std::ptrdiff_t i = 0; i < nIndices; i++)
		vertices2[i] = (firstVertex[firstIndex[i]]*basis).xy();

	auto [outEnd, result] = triangulation::triangulate2(vertices2, vertices2 + nIndices, outIndex);
	for (; outIndex != outEnd; ++outIndex)
		*outIndex = firstIndex[*outIndex];

	return { outIndex, result };
}

template<typename T>
inline Vector3<T> Triangle3<T>::getClosestPoint(const Vector3<T>& point) const
{
	Vector3<T> closestPoint(Uninitialized());
	distances::getPointTriangleSquared(point, vertices[0], vertices[1], vertices[2], &closestPoint);
	return closestPoint;
}

template<typename T>
inline T Triangle3<T>::getDistanceTo(const Vector3<T>& point) const
{
	return distances::getPointTriangle(point, vertices[0], vertices[1], vertices[2]);
}

template<typename T>
inline bool Triangle3<T>::intersects(const AxisAlignedBox<T>& box) const
{
	return intersections::testAxisAlignedBoxTriangle(box.getCenter(), box.getHalfDimensions(), vertices[0], vertices[1], vertices[2]);
}

template<typename T>
inline bool Triangle3<T>::intersects(const OrientedBox<T>& box) const
{
	return intersections::testOrientedBoxTriangle(box.center, box.basis, box.halfDims, vertices[0], vertices[1], vertices[2]);
}

template<typename T>
inline bool Triangle3<T>::intersects(const Sphere<T>& sphere) const
{
	return (distances::getPointTriangleSquared(sphere.center, vertices[0], vertices[1], vertices[2]) <= sphere.radius*sphere.radius);
}

} // namespace mathematics::templates
