/*
 *	Name: TriangleMesh
 *	Author: Pawel Mrochen
 */

#pragma once

#include <stdexcept>
#include <limits>
#include <type_traits>
#include <concepts>
#include <tuple>
#include <vector>
#include <iterator>
#include <initializer_list>
//#include <atomic>
#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <cmath>
#include "Constants.hpp"
#include "Axis.hpp"
#include "Vector3.hpp"
#include "Matrix3.hpp"
#include "AffineTransform.hpp"
#include "HalfSpace.hpp"
#include "Plane.hpp"
#include "Triangle3.hpp"
#include "AxisAlignedBox.hpp"
#include "OrientedBox.hpp"
#include "Sphere.hpp"
#include "SymmetricFrustum.hpp"
#include "ConvexPolyhedron.hpp"

namespace mathematics {
namespace templates {

template<typename T, typename U>
	requires (std::floating_point<T> && std::integral<U>)
class TriangleMesh
{
public:
	using Real = T;
	using VertexType = Vector3<T>/*std::tuple<T, T, T>*/;
	using VertexVector = std::vector<VertexType>;
	using VertexIndexType = U;
	using VertexIndexVector = std::vector<VertexIndexType>;

	TriangleMesh() noexcept : /*refCount_(),*/ vertices(), indices() {}
	template<std::input_iterator<Vector3<T>> I, std::sentinel_for<I> S> TriangleMesh(I first, S last);
	TriangleMesh(const Vector3<T>* vertices, std::size_t nVertices);
	explicit TriangleMesh(const std::vector<Vector3<T>>& vertices);
	explicit TriangleMesh(std::vector<Vector3<T>>&& vertices);
	TriangleMesh(const Vector3<T>* vertices, std::size_t nVertices, const U* indices, std::size_t nIndices);
	template<std::integral V> TriangleMesh(const Vector3<T>* vertices, std::size_t nVertices, const V* indices, std::size_t nIndices);
	TriangleMesh(const std::vector<Vector3<T>>& vertices, const std::vector<U>& indices);
	template<std::integral V> TriangleMesh(const std::vector<Vector3<T>>& vertices, const std::vector<V>& indices);
	TriangleMesh(std::vector<Vector3<T>>&& vertices, std::vector<U>&& indices);
	template<std::input_iterator<Triangle3<T>> I, std::sentinel_for<I> S> TriangleMesh(I first, S last);
	TriangleMesh(const Triangle3<T>* triangles, std::size_t nTriangles);
	explicit TriangleMesh(const std::vector<Triangle3<T>>& triangles);
	TriangleMesh(std::initializer_list<Vector3<T>> vertices) : /*refCount_(),*/ vertices(vertices), indices() {}
	TriangleMesh(const TriangleMesh& mesh);
	TriangleMesh(TriangleMesh&& mesh);
	TriangleMesh& operator=(const TriangleMesh& mesh);
	TriangleMesh& operator=(TriangleMesh&& mesh);

	// Serialization
	template<typename A> void serialize(A& ar) { ar(vertices, indices); }

	// Create
	static TriangleMesh* from(const AxisAlignedBox<T>& box);
	static TriangleMesh* from(const OrientedBox<T>& box);
	static TriangleMesh* from(const SymmetricFrustum<T>& frustum);
	static TriangleMesh* from(const ConvexPolyhedron<T>* convexPolyhedron, T tolerance = Constants<T>::TOLERANCE*T(1000));
	static TriangleMesh* makeTetrahedron(T edgeLength);
	static TriangleMesh* makeCube(const Vector3<T>& dimensions);
#if MATHEMATICS_HAS_QUICKHULL
	//template<std::input_iterator<Vector3<T>> I, std::sentinel_for<I> S> static TriangleMesh* makeConvexHull(I first, S last);
	//static TriangleMesh* makeConvexHull(const Vector3<T>* points, std::size_t nPoints); // convex hull of point cloud
	//static TriangleMesh* makeConvexHull(const std::vector<Vector3<T>>& points);
#endif

	// Clone
	TriangleMesh* clone() const { return new TriangleMesh(*this); }

	// References
	//bool hasOwner() const noexcept { return (refCount_ > 0); }
	//int getReferenceCount() const noexcept { return refCount_; }
	//void acquire() noexcept { refCount_++; }
	//void release() { if (--refCount_ <= 0) delete this; }

	// Clear
	bool isEmpty() const noexcept { return vertices.empty() && indices.empty(); }
	void clear();

	// Assign/append
	void assign(const TriangleMesh* mesh);
	void append(const TriangleMesh* mesh);

	// Integrity
	bool hasValidVertexIndices() const noexcept;

	// Vertices
	std::size_t getVertexCount() const noexcept { return vertices.size(); }
	void setVertexCount(std::size_t count) { vertices.resize(count); }
	const Vector3<T>& getVertex(std::size_t index) const noexcept;
	void setVertex(std::size_t index, const Vector3<T>& value); // throw (std::out_of_range)
	void addVertex(const Vector3<T>& value) { vertices.push_back(value); }
	template<std::input_iterator<Vector3<T>> I, std::sentinel_for<I> S> void addVertices(I first, S last);
	void addVertices(const Vector3<T>* vertices, std::size_t nVertices);
	void addVertices(const std::vector<Vector3<T>>& vertices);
	void insertVertex(std::size_t index, const Vector3<T>& value); // throw (std::out_of_range)
	void removeVertexAt(std::size_t index);
	void removeAllVertices() { vertices.resize(0); }

	// Vertex indices
	std::size_t getVertexIndexCount() const noexcept { return indices.size(); }
	void setVertexIndexCount(std::size_t count) { indices.resize(count, U()); }
	U getVertexIndex(std::size_t primitive, int vertex) const noexcept;
	U getVertexIndex(std::size_t index) const noexcept;
	void setVertexIndex(std::size_t index, U value); // throw (std::out_of_range)
	void addVertexIndex(U value) { indices.push_back(value); }
	void addVertexIndices(U v0, U v1, U v2);
	template<std::input_iterator I, std::sentinel_for<I> S> requires std::integral<std::iter_value_t<I>> void addVertexIndices(I first, S last);
	void addVertexIndices(const U* indices, std::size_t nIndices);
	void addVertexIndices(const std::vector<U>& indices);
	template<std::integral V> void addVertexIndices(const std::vector<V>& indices);
	void insertVertexIndex(std::size_t index, U value); // throw (std::out_of_range)
	void insertVertexIndices(std::size_t index, U v0, U v1, U v2); // throw (std::out_of_range)
	void removeVertexIndexAt(std::size_t index);
	bool computeVertexIndices();
	void removeAllVertexIndices() { indices.resize(0); }

	// Triangles
	std::size_t getPrimitiveCount() const noexcept { return getTriangleCount(); }
	std::size_t getTriangleCount() const noexcept;
	Triangle3<T> getTriangle(std::size_t index) const noexcept;
	template<std::output_iterator<Triangle3<T>> O> O copyTriangles(O target) const;
	std::vector<Triangle3<T>> getTriangles() const;

	// Half spaces
	template<std::output_iterator<HalfSpace<T>> O> O copyHalfSpaces(O target) const;
	std::vector<HalfSpace<T>> getHalfSpaces() const;
	template<std::output_iterator<Plane<T>> O> O copyPlanes(O target) const;
	std::vector<Plane<T>> getPlanes() const;

	// Convexity
	bool isConvex() const noexcept;

	// Surface area and volume
	Vector3<T> computeCentroid() const noexcept;
	T computeSurfaceArea() const noexcept;
	T computeVolume() const noexcept;

	// Bounding box and sphere
	AxisAlignedBox<T> computeAxisAlignedBoundingBox() const noexcept;
	//OrientedBox<T> computeOrientedBoundingBox() const noexcept;
	Sphere<T> computeBoundingSphere() const noexcept;

	// Transform
	TriangleMesh& translate(const Vector3<T>& offset) noexcept;
	TriangleMesh& transform(const Matrix3<T>& matrix) noexcept;
	TriangleMesh& transform(const AffineTransform<T>& transformation) noexcept;

	// Optimization
	bool weldVertices(T tolerance = Constants<T>::TOLERANCE); // throw (std::out_of_range)

	// Containment
	bool contains(const Vector3<T>& point) const noexcept;

	VertexVector vertices;
	VertexIndexVector indices;

private:
	//std::atomic_int refCount_;
};

template<typename T, typename U>
template<std::input_iterator<Vector3<T>> I, std::sentinel_for<I> S> 
inline TriangleMesh<T, U>::TriangleMesh(I first, S last) :
	//refCount_(), 
	vertices(first, last), 
	indices() 
{
}

template<typename T, typename U>
inline TriangleMesh<T, U>::TriangleMesh(const Vector3<T>* vertices, std::size_t nVertices) :
	//refCount_(), 
	vertices(vertices, vertices + nVertices), 
	indices() 
{
}

template<typename T, typename U>
inline TriangleMesh<T, U>::TriangleMesh(const std::vector<Vector3<T>>& vertices) : 
	//refCount_(), 
	vertices(vertices), 
	indices() 
{
}

template<typename T, typename U>
inline TriangleMesh<T, U>::TriangleMesh(std::vector<Vector3<T>>&& vertices) : 
	//refCount_(), 
	vertices(std::move(vertices)), 
	indices() 
{
}

template<typename T, typename U>
inline TriangleMesh<T, U>::TriangleMesh(const Vector3<T>* vertices, std::size_t nVertices, const U* indices, std::size_t nIndices) :
	//refCount_(), 
	vertices(vertices + nVertices), 
	indices(indices + nIndices) 
{
}

template<typename T, typename U>
template<std::integral V> 
inline TriangleMesh<T, U>::TriangleMesh(const Vector3<T>* vertices, std::size_t nVertices, const V* indices, std::size_t nIndices) :
	//refCount_(), 
	vertices(vertices + nVertices), 
	indices(indices + nIndices) 
{
}

template<typename T, typename U>
inline TriangleMesh<T, U>::TriangleMesh(const std::vector<Vector3<T>>& vertices, const std::vector<U>& indices) :
	//refCount_(), 
	vertices(vertices), 
	indices(indices) 
{
}

template<typename T, typename U>
template<std::integral V> 
inline TriangleMesh<T, U>::TriangleMesh(const std::vector<Vector3<T>>& vertices, const std::vector<V>& indices) :
	//refCount_(), 
	vertices(vertices), 
	indices(indices.begin(), indices.end()) 
{
}

template<typename T, typename U>
inline TriangleMesh<T, U>::TriangleMesh(std::vector<Vector3<T>>&& vertices, std::vector<U>&& indices) :
	//refCount_(), 
	vertices(std::move(vertices)), 
	indices(std::move(indices))
{
}

template<typename T, typename U>
template<std::input_iterator<Triangle<T>> I, std::sentinel_for<I> S> 
TriangleMesh<T, U>::TriangleMesh(I first, S last) :
	//refCount_(), 
	vertices(), 
	indices()
{
	if (first != last)
	{
		vertices.reserve(std::distance(first, last)*3);
		for (; first != last; ++first)
		{
			vertices.push_back(first->vertices[0]);
			vertices.push_back(first->vertices[1]);
			vertices.push_back(first->vertices[2]);
		}
	}
}

template<typename T, typename U>
TriangleMesh<T, U>::TriangleMesh(const Triangle<T>* triangles, std::size_t nTriangles) :
	//refCount_(), 
	vertices(), 
	indices()
{
	if (triangles && nTriangles)
	{
		vertices.reserve(nTriangles*3);
		for (size_t i = 0; i != nTriangles; i++)
		{
			vertices.push_back(triangles[i].vertices[0]);
			vertices.push_back(triangles[i].vertices[1]);
			vertices.push_back(triangles[i].vertices[2]);
		}
	}
}

template<typename T, typename U>
TriangleMesh<T, U>::TriangleMesh(const std::vector<Triangle<T>>& triangles) :
	//refCount_(), 
	vertices(), 
	indices()
{
	if (!triangles.empty())
	{
		vertices.reserve(triangles.size()*3);
		for (size_t i = 0, n = triangles.size(); i != n; i++)
		{
			vertices.push_back(triangles[i].vertices[0]);
			vertices.push_back(triangles[i].vertices[1]);
			vertices.push_back(triangles[i].vertices[2]);
		}
	}
}

template<typename T, typename U>
inline TriangleMesh<T, U>::TriangleMesh(const TriangleMesh& mesh) : 
	//refCount_(), 
	vertices(mesh.vertices), 
	indices(mesh.indices) 
{
}

template<typename T, typename U>
inline TriangleMesh<T, U>::TriangleMesh(TriangleMesh&& mesh) : 
	//refCount_(), 
	vertices(std::move(mesh.vertices)), 
	indices(std::move(mesh.indices)) 
{
}

template<typename T, typename U>
inline TriangleMesh<T, U>& TriangleMesh<T, U>::operator=(const TriangleMesh<T, U>& mesh)
{
	vertices = mesh.vertices;
	indices = mesh.indices;
	return *this;
}

template<typename T, typename U>
inline TriangleMesh<T, U>& TriangleMesh<T, U>::operator=(TriangleMesh<T, U>&& mesh)
{
	vertices = std::move(mesh.vertices);
	indices = std::move(mesh.indices);
	return *this;
}

template<typename T, typename U>
/*static*/ inline TriangleMesh<T, U>* TriangleMesh<T, U>::from(const AxisAlignedBox<T>& box)
{
	TriangleMesh<T, U>* mesh = new TriangleMesh<T, U>(box.getVertices());
	auto [first, last] = box.getPrimitives<U>(3);
	mesh->indices.assign(first, last);
	return mesh;
}

template<typename T, typename U>
/*static*/ inline TriangleMesh<T, U>* TriangleMesh<T, U>::from(const OrientedBox<T>& box)
{
	TriangleMesh<T, U>* mesh = new TriangleMesh<T, U>(box.getVertices());
	auto [first, last] = box.getPrimitives<U>(3);
	mesh->indices.assign(first, last);
	return mesh;
}

template<typename T, typename U>
/*static*/ inline TriangleMesh<T, U>* TriangleMesh<T, U>::from(const SymmetricFrustum<T>& frustum)
{
	TriangleMesh<T, U>* mesh = new TriangleMesh<T, U>(frustum.getVertices());
	auto [first, last] = frustum.getPrimitives<U>(3);
	mesh->indices.assign(first, last);
	return mesh;
}

template<typename T, typename U>
/*static*/ TriangleMesh<T, U>* TriangleMesh<T, U>::makeTetrahedron(T edgeLength)
{
	if (!(edgeLength > T(0)))
		return nullptr;

	edgeLength /= std::sqrt(T(8)/T(3));
	TriangleMesh<T, U>* mesh = new TriangleMesh<T, U>({ Vector3<T>(edgeLength*std::sqrt(T(8)/T(9)), -edgeLength/T(3), T(0)),
		Vector3<T>(-edgeLength*std::sqrt(T(2)/T(9)), -edgeLength/T(3), edgeLength*std::sqrt(T(2)/T(3))),
		Vector3<T>(-edgeLength*std::sqrt(T(2)/T(9)), -edgeLength/T(3), -edgeLength*std::sqrt(T(2)/T(3))),
		Vector3<T>(T(0), edgeLength, T(0)) });
	mesh->indices = { 0, 1, 2, 0, 3, 1, 1, 3, 2, 2, 3, 0 };
	return mesh;
}

template<typename T, typename U>
/*static*/ inline TriangleMesh<T, U>* TriangleMesh<T, U>::makeCube(const Vector3<T>& dimensions)
{
	return dimemsions.allGreaterThan(Vector3<T>::ZERO) ?
		from(AxisAlignedBox<T>(dimensions)) :
		nullptr;
}

template<typename T, typename U>
inline void TriangleMesh<T, U>::clear()
{
	vertices.resize(0);
	indices.resize(0);
}

template<typename T, typename U>
inline void TriangleMesh<T, U>::assign(const TriangleMesh<T, U>* mesh)
{
	if (mesh)
	{
		vertices.assign(mesh->vertices.begin(), mesh->vertices.end());
		indices.assign(mesh->indices.begin(), mesh->indices.end());
	}
	else
	{
		vertices.resize(0);
		indices.resize(0);
	}
}

template<typename T, typename U>
void TriangleMesh<T, U>::append(const TriangleMesh<T, U>* mesh)
{
	if (mesh)
	{
		std::size_t base = U(vertices.size());
		std::size_t begin = indices.size();

		vertices.insert(vertices.end(), mesh->vertices.begin(), mesh->vertices.end());
		indices.insert(indices.end(), mesh->indices.begin(), mesh->indices.end());

		std::size_t count = indices.size() - begin;
		if (base && count)
		{
			for (std::size_t i = begin, n = begin + count; i != n; i++)
				indices[i] += base;
		}
	}
}

template<typename T, typename U>
bool TriangleMesh::hasValidVertexIndices() const
{
	std::size_t nVertices = vertices.size();
	for (std::size_t i = 0, n = indices.size(); i != n; i++)
	{
		if (indices[i] >= nVertices)
			return false;
	}

	return true;
}

template<typename T, typename U>
inline const Vector3<T>& TriangleMesh<T, U>::getVertex(std::size_t index) const
{ 
	return (index < vertices.size()) ? Vector3<T>(vertices[index]) : Vector3<T>::ZERO; 
}

template<typename T, typename U>
inline void TriangleMesh<T, U>::setVertex(std::size_t index, const Vector3<T>& value) 
{ 
	if (index >= vertices.size()) 
		throw std::out_of_range("TriangleMesh::setVertex() : index"); 
	vertices[index] = value; 
}

template<typename T, typename U>
template<std::input_iterator<Vector3<T>> I, std::sentinel_for<I> S> 
inline void TriangleMesh<T, U>::addVertices(I first, S last)
{
	vertices.insert(vertices.end(), first, last); 
}

template<typename T, typename U>
inline void TriangleMesh<T, U>::addVertices(const Vector3<T>* vertices, std::size_t nVertices) 
{ 
	this->vertices.insert(this->vertices.end(), vertices, vertices + nVertices); 
}

template<typename T, typename U>
inline void TriangleMesh<T, U>::addVertices(const std::vector<Vector3<T>>& vertices) 
{ 
	this->vertices.insert(this->vertices.end(), vertices.begin(), vertices.end()); 
}

template<typename T, typename U>
inline void TriangleMesh<T, U>::insertVertex(std::size_t index, const Vector3<T>& value) 
{ 
	if (index > vertices.size()) 
		throw std::out_of_range("TriangleMesh::insertVertex() : index"); 
	vertices.insert(vertices.begin() + index, value); 
}

template<typename T, typename U>
inline void TriangleMesh<T, U>::removeVertexAt(std::size_t index) 
{ 
	if (index < vertices.size()) 
		vertices.erase(vertices.begin() + index); 
}

template<typename T, typename U>
inline U TriangleMesh<T, U>::getVertexIndex(std::size_t primitive, int vertex) const
{
	std::size_t index = (primitive << 1) + primitive + vertex;
	return indices.empty() ? U(index) : ((index < indices.size()) ? indices[index] : U());
}

template<typename T, typename U>
inline U TriangleMesh<T, U>::getVertexIndex(std::size_t index) const
{ 
	return (index < indices.size()) ? indices[index] : U(); 
}

template<typename T, typename U>
inline void TriangleMesh<T, U>::setVertexIndex(std::size_t index, U value) 
{ 
	if (index >= indices.size()) 
		throw std::out_of_range("TriangleMesh::setVertexIndex() : index"); 
	indices[index] = value; 
}

template<typename T, typename U>
inline void TriangleMesh<T, U>::addVertexIndices(U v0, U v1, U v2)
{
	indices.push_back(v0);
	indices.push_back(v1);
	indices.push_back(v2);
}

template<typename T, typename U>
template<std::input_iterator I, std::sentinel_for<I> S>
inline void TriangleMesh<T, U>::addVertexIndices(I first, S last)
{
	indices.insert(indices.end(), first, last); 
}

template<typename T, typename U>
inline void TriangleMesh<T, U>::addVertexIndices(const U* indices, std::size_t nIndices) 
{ 
	this->indices.insert(this->indices.end(), indices, indices + nIndices); 
}

template<typename T, typename U>
inline void TriangleMesh<T, U>::addVertexIndices(const std::vector<U>& indices) 
{ 
	this->indices.insert(this->indices.end(), indices.begin(), indices.end()); 
}

template<typename T, typename U>
template<std::integral V>
inline void TriangleMesh<T, U>::addVertexIndices(const std::vector<V>& indices) 
{ 
	this->indices.insert(this->indices.end(), indices.begin(), indices.end()); 
}

template<typename T, typename U>
inline void TriangleMesh<T, U>::insertVertexIndex(std::size_t index, U value) 
{ 
	if (index > indices.size()) 
		throw std::out_of_range("TriangleMesh::insertVertexIndex() : index"); 
	indices.insert(indices.begin() + index, value); 
}

template<typename T, typename U>
inline void TriangleMesh<T, U>::insertVertexIndices(std::size_t index, U v0, U v1, U v2)
{
	if (index > indices.size()) 
		throw std::out_of_range("TriangleMesh::insertVertexIndices() : index"); 
	indices.insert(indices.begin() + index, v0); 
	indices.insert(indices.begin() + (index + 1), v1); 
	indices.insert(indices.begin() + (index + 2), v2); 
}

template<typename T, typename U>
inline void TriangleMesh<T, U>::removeVertexIndexAt(std::size_t index) 
{ 
	if (index < indices.size()) 
		indices.erase(indices.begin() + index); 
}

template<typename T, typename U>
bool TriangleMesh<T, U>::computeVertexIndices()
{
	if (indices.empty())
	{
		indices.resize(vertices.size());
		for (std::size_t i = 0, n = vertices.size(); i != n; i++) 
			indices[i] = U(i);
		return true;
	}

	return false;
}

template<typename T, typename U>
inline std::size_t TriangleMesh<T, U>::getTriangleCount() const
{ 
	return indices.empty() ? vertices.size()/3 : indices.size()/3;
}

template<typename T, typename U>
inline Triangle3<T> TriangleMesh<T, U>::getTriangle(std::size_t index) const noexcept
{
	return { getVertex(getVertexIndex(index, 0)), getVertex(getVertexIndex(index, 1)),
		getVertex(getVertexIndex(index, 2)) };
}

template<typename T, typename U>
template<std::output_iterator<Triangle3<T>> O> 
O TriangleMesh<T, U>::copyTriangles(O target) const
{
	for (std::size_t i = 0, n = getTriangleCount(); i != n; i++)
	{
		*target++ = Triangle3<T>(getVertex(getVertexIndex(i, 0)), getVertex(getVertexIndex(i, 1)),
			getVertex(getVertexIndex(i, 2)));
	}

	return target;
}

template<typename T, typename U>
std::vector<Triangle3<T>> TriangleMesh<T, U>::getTriangles() const
{
	std::size_t n = getTriangleCount();
	std::vector<Triangle3<T>> triangles;
	triangles.reserve(n);

	for (std::size_t i = 0; i != n; i++)
	{
		triangles.emplace_back(getVertex(getVertexIndex(i, 0)), getVertex(getVertexIndex(i, 1)),
			getVertex(getVertexIndex(i, 2)));
	}

	return triangles;
}

template<typename T, typename U>
template<std::output_iterator<HalfSpace<T>> O> 
O TriangleMesh<T, U>::copyHalfSpaces(O target) const
{
	for (std::size_t i = 0, n = getTriangleCount(); i != n; i++)
	{
		*target++ = HalfSpace<T>(getVertex(getVertexIndex(i, 0)), getVertex(getVertexIndex(i, 1)),
			getVertex(getVertexIndex(i, 2)));
	}

	return target;
}

template<typename T, typename U>
std::vector<HalfSpace<T>> TriangleMesh<T, U>::getHalfSpaces() const
{
	std::size_t n = getTriangleCount();
	std::vector<HalfSpace<T>> halfSpaces;
	halfSpaces.reserve(n);

	for (std::size_t i = 0; i != n; i++)
	{
		halfSpaces.emplace_back(getVertex(getVertexIndex(i, 0)), getVertex(getVertexIndex(i, 1)),
			getVertex(getVertexIndex(i, 2)));
	}

	return halfSpaces;
}

template<typename T, typename U>
template<std::output_iterator<Plane<T>> O> 
O TriangleMesh<T, U>::copyPlanes(O target) const
{
	for (std::size_t i = 0, n = getTriangleCount(); i != n; i++)
	{
		*target++ = Plane<T>(getVertex(getVertexIndex(i, 0)), getVertex(getVertexIndex(i, 1)),
			getVertex(getVertexIndex(i, 2)));
	}

	return target;
}

template<typename T, typename U>
std::vector<Plane<T>> TriangleMesh<T, U>::getPlanes() const
{
	std::size_t n = getTriangleCount();
	std::vector<Plane<T>> planes;
	planes.reserve(n);

	for (std::size_t i = 0; i != n; i++)
	{
		planes.emplace_back(getVertex(getVertexIndex(i, 0)), getVertex(getVertexIndex(i, 1)),
			getVertex(getVertexIndex(i, 2)));
	}

	return planes;
}

template<typename T, typename U>
bool TriangleMesh<T, U>::isConvex() const
{
	std::size_t vertexCount = vertices.size();
	for (std::size_t i = 0, n = getTriangleCount(); i != n; i++)
	{
		Plane<T> p(getVertex(getVertexIndex(i, 0)), getVertex(getVertexIndex(i, 1)),
			getVertex(getVertexIndex(i, 2)));

		for (std::size_t j = 0; j != vertexCount; j++)
		{
			if (p.getSignedDistanceTo(vertices[j]) >= Constants<T>::TOLERANCE)
				return false;
		}
	}

	return true;
}

template<typename T, typename U>
Vector3<T> TriangleMesh<T, U>::computeCentroid() const
{
	if (vertices.empty())
		return Vector3::ZERO;
    
	Vector3<T> average = vertices[0];
    for (std::size_t i = 1, n = vertices.size(); i != n; i++)
        average += vertices[i];
    
	average /= T(vertices.size());
    return average;
}

template<typename T, typename U>
T TriangleMesh<T, U>::computeSurfaceArea() const
{
	T area = T(0);
	for (std::size_t i = 0, n = getTriangleCount(); i != n; i++)
	{
		const Vector3<T>& v0 = getVertex(getVertexIndex(i, 0));
		const Vector3<T>& v1 = getVertex(getVertexIndex(i, 1));
		const Vector3<T>& v2 = getVertex(getVertexIndex(i, 2));
        area += cross(v1 - v0, v2 - v0).getLength();
	}

	return area*T(0.5);
}

template<typename T, typename U>
T TriangleMesh<T, U>::computeVolume() const
{
	T vol = T(0);
	for (std::size_t i = 0, n = getTriangleCount(); i != n; i++)
	{
		const Vector3<T>& v0 = getVertex(getVertexIndex(i, 0));
		const Vector3<T>& v1 = getVertex(getVertexIndex(i, 1));
		const Vector3<T>& v2 = getVertex(getVertexIndex(i, 2));
		vol += dot(v0, cross(v1, v2));
	}

	return vol/T(6);
}

template<typename T, typename U>
AxisAlignedBox<T> TriangleMesh<T, U>::computeAxisAlignedBoundingBox() const
{
	if (vertices.empty())
		return {};
    
	Vector3<T> low = vertices[0];
	Vector3<T> high = low;
    for (std::size_t i = 1, n = vertices.size(); i != n; i++)
	{
		const Vector3<T>& vertex = vertices[i];
		low.setMinimum(low, vertex);
		high.setMaximum(high, vertex);
	}

	return { low, high };
}

template<typename T, typename U>
Sphere<T> TriangleMesh<T, U>::computeBoundingSphere() const
{
	if (vertices.empty())
		return {};

	if (vertices.size() == 1)
		return { vertex[0], T(0) };

	Vector3<T> minX = vertices[0];
	Vector3<T> maxX = minX;
	Vector3<T> minY = minX;
	Vector3<T> maxY = minX;
	Vector3<T> minZ = minX;
	Vector3<T> maxZ = minX;
    for (std::size_t i = 1, n = vertices.size(); i != n; i++)
	{
		const Vector3<T>& vertex = vertices[i];
		if (vertex.x < minX.x) 
			minX = vertex;
		if (vertex.x > maxX.x) 
			maxX = vertex;
		if (vertex.y < minY.y) 
			minY = vertex;
		if (vertex.y > maxY.y) 
			maxY = vertex;
		if (vertex.z < minZ.z) 
			minZ = vertex;
		if (vertex.z > maxZ.z) 
			maxZ = vertex;
	}

	Axis axis = Vector3<T>(distanceSquared(minX, maxX), distanceSquared(minY, maxY),
		distanceSquared(minZ, maxZ)).getMajorAxis();
	Vector3<T> dia1 = (axis == Axis::X) ? minX : ((axis == Axis::Y) ? minY : minZ);
	Vector3<T> dia2 = (axis == Axis::X) ? maxX : ((axis == Axis::Y) ? maxY : maxZ);
	Vector3<T> center = (dia1 + dia2)*T(0.5);
	T radiusSq = distanceSquared(center, dia2);
	T radius = std::sqrt(radiusSq);

    for (const Vector3<T>& vertex : vertices)
	{
		T oldToPSq = distanceSquared(center, vertex);
		if (oldToPSq > radiusSq)
		{
			T oldToP = std::sqrt(oldToPSq);
			radius = (radius + oldToP)*T(0.5);
			radiusSq = radius*radius;
			center = (radius*center + (oldToP - radius)*vertex)/oldToP;
		}
	}

	return { center, radius };
}

template<typename T, typename U>
void TriangleMesh<T, U>::translate(const Vector3<T>& offset)
{
	if (offset.isZero())
		return;

	for (Vector3<T>& v : vertices)
		v += offset;
}

template<typename T, typename U>
void TriangleMesh<T, U>::transform(const Matrix3<T>& matrix)
{
	if (vertices.empty() || matrix.isIdentity())
		return;

	for (Vector3<T>& v : vertices)
		v *= matrix;
}

template<typename T, typename U>
void TriangleMesh<T, U>::transform(const AffineTransform<T>& transformation)
{
	if (vertices.empty() || transformation.isIdentity())
		return;

	if (transformation.getBasis().isIdentity())
	{
		translate(transformation.getOrigin());
	}
	else
	{
		for (Vector3<T>& v : vertices)
			v.transform(transformation); 
	}
}

template<typename T, typename U>
bool TriangleMesh<T, U>::weldVertices(T tolerance)
{
	if (vertices.size() <= 1)
		return false;

	if (vertices.size() > (std::size_t(std::numeric_limits<U>::max()) + 1))
		throw std::out_of_range("TriangleMesh::weldVertices() : vertices.size()");

	std::vector<U> mappingTable(vertices.size());
	mappingTable[0] = U(0);
	bool remapped = false;

	for (std::size_t i = 1, n = vertices.size(); i < n; i++)
	{
		mappingTable[i] = U(i);

		for (std::ptrdiff_t j = std::ptrdiff_t(i) - 1; j >= std::ptrdiff_t(0); j--)
		{
			if (vertices[j].approxEquals(vertices[i], tolerance))
			{
				mappingTable[i] = mappingTable[j];
				remapped = true;
				break;
			}
		}
	}

	if (!remapped)
		return false;

	if (indices.empty())
		computeVertexIndices();

	U index = 0;
	for (std::size_t i = 0, n = vertices.size(); i != n; i++)
	{
		if (mappingTable[i] == U(i))
		{
			mappingTable[i] = index;
			if (index != U(i))
				vertices[index] = vertices[i];
			index++;
		}
		else
		{
			mappingTable[i] = mappingTable[mappingTable[i]];
		}
	}

	vertices.resize(index);
	for (std::size_t i = 0, n = indices.size(); i != n; i++)
		indices[i] = mappingTable[indices[i]];

	return true;
}

template<typename T, typename U>
bool TriangleMesh<T, U>::contains(const Vector3<T>& point) const
{
	for (std::size_t i = 0, n = getTriangleCount(); i != n; i++)
	{
		const Vector3& v0 = getVertex(getVertexIndex(i, 0));
		const Vector3& v1 = getVertex(getVertexIndex(i, 1));
		const Vector3& v2 = getVertex(getVertexIndex(i, 2));
		if (!HalfSpace<T>(v0, v1, v2).contains(point))
			return false;
	}

	return true;
}

} // namespace templates

#if MATHEMATICS_DOUBLE
using TriangleMesh = templates::TriangleMesh<double, std::uint32_t>;
#else
using TriangleMesh = templates::TriangleMesh<float, std::uint32_t>;
#endif

} // namespace mathematics

#include "LU.inl"

namespace mathematics::templates {

template<typename T, typename U>
/*static*/ TriangleMesh<T, U>* TriangleMesh<T, U>::from(const ConvexPolyhedron<T>* convexPolyhedron, T tolerance)
{
	if (!convexPolyhedron)
		return nullptr;

	TriangleMesh<T, U>* mesh = new TriangleMesh<T, U>;
	if (convexPolyhedron->halfSpaces.size() < 4)
		return mesh;

	std::size_t nHalfSpaces = convexPolyhedron->halfSpaces.size();
	std::vector<std::tuple<int, int, int>> vertexPlanes;
	using VerticesOnPlane = std::pair<std::vector<Vector3<T>>::iterator*, std::size_t>;
	std::vector<VerticesOnPlane> verticesOnPlanes(nHalfSpaces, { {}, {} });
	std::size_t nVerticesOnAllPlanes = 0;

	for (auto i0 = convexPolyhedron->halfSpaces.begin(); i0 != convexPolyhedron->halfSpaces.end(); ++i0)
	{ 
		auto i1 = i0; 
        for (++i1; i1 != convexPolyhedron->halfSpaces.end(); ++i1)
		{ 
			auto i2 = i1; 
            for (++i2; i2 != convexPolyhedron->halfSpaces.end(); ++i2)
			{ 
                // Test if the planes intersect at a point 
				auto d = lu::decompose(Matrix3<T>(i0->getNormal(), i1->getNormal(), i2->getNormal()));
                if (d.has_value())
				{ 
					Vector3<T> point = lu::solve(d.value(), Vector3<T>(-i0->d, -i1->d, -i2->d));
					if (point.isFinite() && convexPolyhedron->contains(point, tolerance))
					{
						auto range = std::equal_range(mesh->vertices.begin(), mesh->vertices.end(), point,
							[](const Vector3<T>& lhs, const Vector3<T>& rhs)
							{
								if (lhs.x < rhs.x) 
									return true; 
								if (lhs.x > rhs.x) 
									return false; 
								if (lhs.y < rhs.y) 
									return true; 
								if (lhs.y > rhs.y) 
									return false; 
								return (lhs.z < rhs.z); 
							});

						if (range.first == range.second)	// point not found
						{
							std::ptrdiff_t p0 = std::distance(convexPolyhedron->halfSpaces.begin(), i0);
							std::ptrdiff_t p1 = std::distance(convexPolyhedron->halfSpaces.begin(), i1);
							std::ptrdiff_t p2 = std::distance(convexPolyhedron->halfSpaces.begin(), i2);

							verticesOnPlanes[p0].second++;
							verticesOnPlanes[p1].second++;
							verticesOnPlanes[p2].second++;
							nVerticesOnAllPlanes += 3;
														
							std::ptrdiff_t vertexIndex = std::distance(mesh->vertices.begin(), range.first);
							mesh->vertices.insert(range.first, point);	// intersection of 3 planes
							vertexPlanes.insert(vertexPlanes.begin() + vertexIndex, int(p0), int(p1), int(p2));
						}
					}
				}
			}
		}
    }

	std::vector<std::vector<Vector3<T>>::iterator> iteratorBuffer(nVerticesOnAllPlanes);
	std::vector<Vector3<T>>::iterator* iteratorPtr = iteratorBuffer.data();
	for (std::size_t i = 0; i != nHalfSpaces; i++)
	{
		verticesOnPlanes[i].first = iteratorPtr;
		iteratorPtr += verticesOnPlanes[i].second;
		verticesOnPlanes[i].second = 0;
	}

	std::size_t nVertices = mesh->vertices.size();
	for (std::size_t i = 0; i != nVertices; i++)
	{
		std::ptrdiff_t p0 = std::get<0>(vertexPlanes[i]);
		std::ptrdiff_t p1 = std::get<1>(vertexPlanes[i]);
		std::ptrdiff_t p2 = std::get<2>(vertexPlanes[i]);
		verticesOnPlanes[p0].first[verticesOnPlanes[p0].second++] = mesh->vertices.begin() + i;
		verticesOnPlanes[p1].first[verticesOnPlanes[p1].second++] = mesh->vertices.begin() + i;
		verticesOnPlanes[p2].first[verticesOnPlanes[p2].second++] = mesh->vertices.begin() + i;
	}
	
	std::vector<Vector3<T>> triangulateInput;
	std::vector<U> triangulateOutput;
	for (std::size_t i = 0; i != nHalfSpaces; i++)
	{
		std::size_t nVerticesOnPlane = verticesOnPlanes[i].second;
		if (nVerticesOnPlane >= 3)
		{
			std::vector<Vector3<T>>::iterator* iVertices = verticesOnPlanes[i].first;
			Vector3<T> center = *iVertices[0];
			for (std::size_t j = 1; j < nVerticesOnPlane; j++)
				center += *iVertices[j];
			center /= T(nVerticesOnPlane);
			
			Vector3<T> normal = convexPolyhedron->halfSpaces[i].getNormal();
			auto comparer =
				[&center, &normal](const std::vector<Vector3<T>>::iterator& lhs, const std::vector<Vector3<T>>::iterator& rhs)
				{
					return dot(normal, cross(*lhs - center, *rhs - center)) > T(0);
				};

			std::sort(iVertices, iVertices + nVerticesOnPlane, comparer);
			if (comparer(iVertices[0], iVertices[nVerticesOnPlane - 1])) // correct winding of last-first vertex pair
				std::swap(iVertices[0], iVertices[nVerticesOnPlane - 1]);
				
			if (nVerticesOnPlane > 3)
			{
				triangulateInput.resize(0);
				for (std::size_t j = 0; j != nVerticesOnPlane; j++)
					triangulateInput.push_back(*iVertices[j]);
				
				triangulateOutput.resize(0);
				auto [last, result] = Triangle3<T>::triangulate(triangulateInput.begin(), triangulateInput.end(),
					std::back_inserter(triangulateOutput));

				if (result)
				{
					for (std::size_t j = 0, n = triangulateOutput.size(); j != n; j++)
						mesh->indices.emplace_back(std::distance(mesh->vertices.begin(), iVertices[triangulateOutput[j]]));
				}
			}
			else
			{
				mesh->indices.emplace_back(std::distance(mesh->vertices.begin(), iVertices[0]));
				mesh->indices.emplace_back(std::distance(mesh->vertices.begin(), iVertices[1]));
				mesh->indices.emplace_back(std::distance(mesh->vertices.begin(), iVertices[2]));
			}
		}
	}
		
	return mesh;
}

} // namespace mathematics::templates
