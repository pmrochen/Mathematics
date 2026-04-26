/*
 *	Name: TetrahedralMesh
 *	Author: Pawel Mrochen
 */

#pragma once

#include <stdexcept>
#include <limits>
#include <type_traits>
#include <concepts>
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
#include "Tetrahedron.hpp"
#include "AxisAlignedBox.hpp"
#include "Sphere.hpp"

namespace mathematics {
namespace templates {

template<typename T, typename U>
	requires (std::floating_point<T> && std::integral<U>)
class TetrahedralMesh
{
public:
	using Real = T;
	using VertexType = Vector3<T>/*std::tuple<T, T, T>*/;
	using VertexVector = std::vector<VertexType>;
	using VertexIndexType = U;
	using VertexIndexVector = std::vector<VertexIndexType>;

	TetrahedralMesh() noexcept : /*refCount_(),*/ vertices(), indices() {}
	template<std::input_iterator<Vector3<T>> I, std::sentinel_for<I> S> TetrahedralMesh(I first, S last);
	TetrahedralMesh(const Vector3<T>* vertices, std::size_t nVertices);
	explicit TetrahedralMesh(const std::vector<Vector3<T>>& vertices);
	explicit TetrahedralMesh(std::vector<Vector3<T>>&& vertices);
	TetrahedralMesh(const Vector3<T>* vertices, std::size_t nVertices, const U* indices, std::size_t nIndices);
	template<std::integral V> TetrahedralMesh(const Vector3<T>* vertices, std::size_t nVertices, const V* indices, std::size_t nIndices);
	TetrahedralMesh(const std::vector<Vector3<T>>& vertices, const std::vector<U>& indices);
	template<std::integral V> TetrahedralMesh(const std::vector<Vector3<T>>& vertices, const std::vector<V>& indices);
	TetrahedralMesh(std::vector<Vector3<T>>&& vertices, std::vector<U>&& indices);
	template<std::input_iterator<Tetrahedron<T>> I, std::sentinel_for<I> S> TetrahedralMesh(I first, S last);
	TetrahedralMesh(const Tetrahedron<T>* tetrahedrons, std::size_t nTetrahedrons);
	explicit TetrahedralMesh(const std::vector<Tetrahedron<T>>& tetrahedrons);
	TetrahedralMesh(std::initializer_list<Vector3<T>> vertices) : /*refCount_(),*/ vertices(vertices), indices() {}
	TetrahedralMesh(const TetrahedralMesh& mesh);
	TetrahedralMesh(TetrahedralMesh&& mesh);
	TetrahedralMesh& operator=(const TetrahedralMesh& mesh);
	TetrahedralMesh& operator=(TetrahedralMesh&& mesh);

	// Serialization
	template<typename A> void serialize(A& ar) { ar(vertices, indices); }

	// Clone
	TetrahedralMesh* clone() const { return new TetrahedralMesh(*this); }

	// References
	//bool hasOwner() const noexcept { return (refCount_ > 0); }
	//int getReferenceCount() const noexcept { return refCount_; }
	//void acquire() noexcept { refCount_++; }
	//void release() { if (--refCount_ <= 0) delete this; }

	// Clear
	bool isEmpty() const noexcept { return vertices.empty() && indices.empty(); }
	void clear();

	// Assign/append
	void assign(const TetrahedralMesh* mesh);
	void append(const TetrahedralMesh* mesh);

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
	void addVertexIndices(U v0, U v1, U v2, U v3);
	template<std::input_iterator I, std::sentinel_for<I> S> requires std::integral<std::iter_value_t<I>> void addVertexIndices(I first, S last);
	void addVertexIndices(const U* indices, std::size_t nIndices);
	void addVertexIndices(const std::vector<U>& indices);
	template<std::integral V> void addVertexIndices(const std::vector<V>& indices);
	void insertVertexIndex(std::size_t index, U value); // throw (std::out_of_range)
	void insertVertexIndices(std::size_t index, U v0, U v1, U v2, U v3); // throw (std::out_of_range)
	void removeVertexIndexAt(std::size_t index);
	bool computeVertexIndices();
	void removeAllVertexIndices() { indices.resize(0); }

	// Tetrahedra
	std::size_t getPrimitiveCount() const noexcept { return getTetrahedronCount(); }
	std::size_t getTetrahedronCount() const noexcept;
	Tetrahedron<T> getTetrahedron(std::size_t index) const noexcept;
	template<std::output_iterator<Tetrahedron<T>> O> O copyTetrahedrons(O target) const;
	std::vector<Tetrahedron<T>> getTetrahedrons() const;

	// Area and volume
	Vector3<T> computeCentroid() const noexcept;
	T computeVolume() const noexcept;

	// Bounding box and sphere
	AxisAlignedBox<T> computeAxisAlignedBoundingBox() const noexcept;
	Sphere<T> computeBoundingSphere() const noexcept;

	// Transform
	TetrahedralMesh& translate(const Vector3<T>& offset) noexcept;
	TetrahedralMesh& transform(const Matrix3<T>& matrix) noexcept;
	TetrahedralMesh& transform(const AffineTransform<T>& transformation) noexcept;

	// Optimization
	bool weldVertices(T tolerance = Constants<T>::TOLERANCE); // throw (std::out_of_range)

	VertexVector vertices;
	VertexIndexVector indices;
	//VertexIndexVector boundaryFaceIndices;
	//VertexIndexVector boundaryEdgeIndices;

private:
	//std::atomic_int refCount_;
};

template<typename T, typename U>
template<std::input_iterator<Vector3<T>> I, std::sentinel_for<I> S> 
inline TetrahedralMesh<T, U>::TetrahedralMesh(I first, S last) :
	//refCount_(), 
	vertices(first, last), 
	indices() 
{
}

template<typename T, typename U>
inline TetrahedralMesh<T, U>::TetrahedralMesh(const Vector3<T>* vertices, std::size_t nVertices) :
	//refCount_(), 
	vertices(vertices, vertices + nVertices), 
	indices() 
{
}

template<typename T, typename U>
inline TetrahedralMesh<T, U>::TetrahedralMesh(const std::vector<Vector3<T>>& vertices) : 
	//refCount_(), 
	vertices(vertices), 
	indices() 
{
}

template<typename T, typename U>
inline TetrahedralMesh<T, U>::TetrahedralMesh(std::vector<Vector3<T>>&& vertices) : 
	//refCount_(), 
	vertices(std::move(vertices)), 
	indices() 
{
}

template<typename T, typename U>
inline TetrahedralMesh<T, U>::TetrahedralMesh(const Vector3<T>* vertices, std::size_t nVertices, const U* indices, std::size_t nIndices) :
	//refCount_(), 
	vertices(vertices + nVertices), 
	indices(indices + nIndices) 
{
}

template<typename T, typename U>
template<std::integral V> 
inline TetrahedralMesh<T, U>::TetrahedralMesh(const Vector3<T>* vertices, std::size_t nVertices, const V* indices, std::size_t nIndices) :
	//refCount_(), 
	vertices(vertices + nVertices), 
	indices(indices + nIndices) 
{
}

template<typename T, typename U>
inline TetrahedralMesh<T, U>::TetrahedralMesh(const std::vector<Vector3<T>>& vertices, const std::vector<U>& indices) :
	//refCount_(), 
	vertices(vertices), 
	indices(indices) 
{
}

template<typename T, typename U>
template<std::integral V> 
inline TetrahedralMesh<T, U>::TetrahedralMesh(const std::vector<Vector3<T>>& vertices, const std::vector<V>& indices) :
	//refCount_(), 
	vertices(vertices), 
	indices(indices.begin(), indices.end()) 
{
}

template<typename T, typename U>
inline TetrahedralMesh<T, U>::TetrahedralMesh(std::vector<Vector3<T>>&& vertices, std::vector<U>&& indices) :
	//refCount_(), 
	vertices(std::move(vertices)), 
	indices(std::move(indices))
{
}

template<typename T, typename U>
template<std::input_iterator<Tetrahedron<T>> I, std::sentinel_for<I> S> 
TetrahedralMesh<T, U>::TetrahedralMesh(I first, S last) :
	//refCount_(), 
	vertices(), 
	indices()
{
	if (first != last)
	{
		vertices.reserve(std::distance(first, last)*4);
		for (; first != last; ++first)
		{
			vertices.push_back(first->vertices[0]);
			vertices.push_back(first->vertices[1]);
			vertices.push_back(first->vertices[2]);
			vertices.push_back(first->vertices[3]);
		}
	}
}

template<typename T, typename U>
TetrahedralMesh<T, U>::TetrahedralMesh(const Tetrahedron<T>* tetrahedrons, std::size_t nTetrahedrons) :
	//refCount_(), 
	vertices(), 
	indices()
{
	if (tetrahedrons && nTetrahedrons)
	{
		vertices.reserve(nTetrahedrons*4);
		for (size_t i = 0; i != nTetrahedrons; i++)
		{
			vertices.push_back(tetrahedrons[i].vertices[0]);
			vertices.push_back(tetrahedrons[i].vertices[1]);
			vertices.push_back(tetrahedrons[i].vertices[2]);
			vertices.push_back(tetrahedrons[i].vertices[3]);
		}
	}
}

template<typename T, typename U>
TetrahedralMesh<T, U>::TetrahedralMesh(const std::vector<Tetrahedron<T>>& tetrahedrons) :
	//refCount_(), 
	vertices(), 
	indices()
{
	if (!tetrahedrons.empty())
	{
		vertices.reserve(tetrahedrons.size()*4);
		for (size_t i = 0, n = tetrahedrons.size(); i != n; i++)
		{
			vertices.push_back(tetrahedrons[i].vertices[0]);
			vertices.push_back(tetrahedrons[i].vertices[1]);
			vertices.push_back(tetrahedrons[i].vertices[2]);
			vertices.push_back(tetrahedrons[i].vertices[3]);
		}
	}
}

template<typename T, typename U>
inline TetrahedralMesh<T, U>::TetrahedralMesh(const TetrahedralMesh& mesh) : 
	//refCount_(), 
	vertices(mesh.vertices), 
	indices(mesh.indices) 
{
}
	
template<typename T, typename U>
inline TetrahedralMesh<T, U>::TetrahedralMesh(TetrahedralMesh&& mesh) : 
	//refCount_(), 
	vertices(std::move(mesh.vertices)), 
	indices(std::move(mesh.indices)) 
{
}

template<typename T, typename U>
inline TetrahedralMesh<T, U>& TetrahedralMesh<T, U>::operator=(const TetrahedralMesh<T, U>& mesh)
{
	vertices = mesh.vertices;
	indices = mesh.indices;
	return *this;
}

template<typename T, typename U>
inline TetrahedralMesh<T, U>& TetrahedralMesh<T, U>::operator=(TetrahedralMesh<T, U>&& mesh)
{
	vertices = std::move(mesh.vertices);
	indices = std::move(mesh.indices);
	return *this;
}

template<typename T, typename U>
inline void TetrahedralMesh<T, U>::clear()
{
	vertices.resize(0);
	indices.resize(0);
}

template<typename T, typename U>
inline void TetrahedralMesh<T, U>::assign(const TetrahedralMesh<T, U>* mesh)
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
void TetrahedralMesh<T, U>::append(const TetrahedralMesh<T, U>* mesh)
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
bool TetrahedralMesh::hasValidVertexIndices() const
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
inline const Vector3<T>& TetrahedralMesh<T, U>::getVertex(std::size_t index) const
{ 
	return (index < vertices.size()) ? Vector3<T>(vertices[index]) : Vector3<T>::ZERO; 
}

template<typename T, typename U>
inline void TetrahedralMesh<T, U>::setVertex(std::size_t index, const Vector3<T>& value) 
{ 
	if (index >= vertices.size()) 
		throw std::out_of_range("TetrahedralMesh::setVertex() : index"); 
	vertices[index] = value; 
}

template<typename T, typename U>
template<std::input_iterator<Vector3<T>> I, std::sentinel_for<I> S> 
inline void TetrahedralMesh<T, U>::addVertices(I first, S last)
{
	vertices.insert(vertices.end(), first, last); 
}

template<typename T, typename U>
inline void TetrahedralMesh<T, U>::addVertices(const Vector3<T>* vertices, std::size_t nVertices) 
{ 
	this->vertices.insert(this->vertices.end(), vertices, vertices + nVertices); 
}

template<typename T, typename U>
inline void TetrahedralMesh<T, U>::addVertices(const std::vector<Vector3<T>>& vertices) 
{ 
	this->vertices.insert(this->vertices.end(), vertices.begin(), vertices.end()); 
}

template<typename T, typename U>
inline void TetrahedralMesh<T, U>::insertVertex(std::size_t index, const Vector3<T>& value) 
{ 
	if (index > vertices.size()) 
		throw std::out_of_range("TetrahedralMesh::insertVertex() : index"); 
	vertices.insert(vertices.begin() + index, value); 
}

template<typename T, typename U>
inline void TetrahedralMesh<T, U>::removeVertexAt(std::size_t index) 
{ 
	if (index < vertices.size()) 
		vertices.erase(vertices.begin() + index); 
}

template<typename T, typename U>
inline U TetrahedralMesh<T, U>::getVertexIndex(std::size_t primitive, int vertex) const
{
	std::size_t index = (primitive << 2) + vertex;
	return indices.empty() ? U(index) : ((index < indices.size()) ? indices[index] : U());
}

template<typename T, typename U>
inline U TetrahedralMesh<T, U>::getVertexIndex(std::size_t index) const
{ 
	return (index < indices.size()) ? indices[index] : U(); 
}

template<typename T, typename U>
inline void TetrahedralMesh<T, U>::setVertexIndex(std::size_t index, U value) 
{ 
	if (index >= indices.size()) 
		throw std::out_of_range("TetrahedralMesh::setVertexIndex() : index"); 
	indices[index] = value; 
}

template<typename T, typename U>
inline void TetrahedralMesh<T, U>::addVertexIndices(U v0, U v1, U v2, U v3)
{
	indices.push_back(v0);
	indices.push_back(v1);
	indices.push_back(v2);
	indices.push_back(v3);
}

template<typename T, typename U>
template<std::input_iterator I, std::sentinel_for<I> S>
inline void TetrahedralMesh<T, U>::addVertexIndices(I first, S last)
{
	indices.insert(indices.end(), first, last); 
}

template<typename T, typename U>
inline void TetrahedralMesh<T, U>::addVertexIndices(const U* indices, std::size_t nIndices) 
{ 
	this->indices.insert(this->indices.end(), indices, indices + nIndices); 
}

template<typename T, typename U>
inline void TetrahedralMesh<T, U>::addVertexIndices(const std::vector<U>& indices) 
{ 
	this->indices.insert(this->indices.end(), indices.begin(), indices.end()); 
}

template<typename T, typename U>
template<std::integral V>
inline void TetrahedralMesh<T, U>::addVertexIndices(const std::vector<V>& indices) 
{ 
	this->indices.insert(this->indices.end(), indices.begin(), indices.end()); 
}

template<typename T, typename U>
inline void TetrahedralMesh<T, U>::insertVertexIndex(std::size_t index, U value) 
{ 
	if (index > indices.size()) 
		throw std::out_of_range("TetrahedralMesh::insertVertexIndex() : index"); 
	indices.insert(indices.begin() + index, value); 
}

template<typename T, typename U>
inline void TetrahedralMesh<T, U>::insertVertexIndices(std::size_t index, U v0, U v1, U v2, U v3)
{
	if (index > indices.size()) 
		throw std::out_of_range("TetrahedralMesh::insertVertexIndices() : index"); 
	indices.insert(indices.begin() + index, v0); 
	indices.insert(indices.begin() + (index + 1), v1); 
	indices.insert(indices.begin() + (index + 2), v2); 
	indices.insert(indices.begin() + (index + 3), v3); 
}

template<typename T, typename U>
inline void TetrahedralMesh<T, U>::removeVertexIndexAt(std::size_t index) 
{ 
	if (index < indices.size()) 
		indices.erase(indices.begin() + index); 
}

template<typename T, typename U>
bool TetrahedralMesh<T, U>::computeVertexIndices()
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
inline std::size_t TetrahedralMesh<T, U>::getTetrahedronCount() const
{ 
	return indices.empty() ? (vertices.size() >> 2) : (indices.size() >> 2);
}

template<typename T, typename U>
inline Tetrahedron<T> TetrahedralMesh<T, U>::getTetrahedron(std::size_t index) const noexcept
{
	return { getVertex(getVertexIndex(index, 0)), getVertex(getVertexIndex(index, 1)),
		getVertex(getVertexIndex(index, 2)), getVertex(getVertexIndex(index, 3)) };
}

template<typename T, typename U>
template<std::output_iterator<Tetrahedron<T>> O> 
O TetrahedralMesh<T, U>::copyTetrahedrons(O target) const
{
	for (std::size_t i = 0, n = getTetrahedronCount(); i != n; i++)
	{
		*target++ = Tetrahedron<T>(getVertex(getVertexIndex(i, 0)), getVertex(getVertexIndex(i, 1)),
			getVertex(getVertexIndex(i, 2)), getVertex(getVertexIndex(i, 3)));
	}

	return target;
}

template<typename T, typename U>
std::vector<Tetrahedron<T>> TetrahedralMesh<T, U>::getTetrahedrons() const
{
	std::size_t n = getTetrahedronCount();
	std::vector<Tetrahedron<T>> tetrahedrons;
	tetrahedrons.reserve(n);

	for (std::size_t i = 0; i != n; i++)
	{
		tetrahedrons.emplace_back(getVertex(getVertexIndex(i, 0)), getVertex(getVertexIndex(i, 1)),
			getVertex(getVertexIndex(i, 2)), getVertex(getVertexIndex(i, 3)));
	}

	return tetrahedrons;
}

template<typename T, typename U>
Vector3<T> TetrahedralMesh<T, U>::computeCentroid() const
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
T TetrahedralMesh<T, U>::computeVolume() const
{
	T vol = T(0);
	for (std::size_t i = 0, n = getTetrahedronCount(); i != n; i++)
	{
		const Vector3<T>& v0 = getVertex(getVertexIndex(i, 0));
		const Vector3<T>& v1 = getVertex(getVertexIndex(i, 1));
		const Vector3<T>& v2 = getVertex(getVertexIndex(i, 2));
		const Vector3<T>& v3 = getVertex(getVertexIndex(i, 3));
		vol += std::fabs(dot(v0 - v3, cross(v1 - v3, v2 - v3)));
	}

	return vol/T(6);
}

template<typename T, typename U>
AxisAlignedBox<T> TetrahedralMesh<T, U>::computeAxisAlignedBoundingBox() const
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
Sphere<T> TetrahedralMesh<T, U>::computeBoundingSphere() const
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
void TetrahedralMesh<T, U>::translate(const Vector3<T>& offset)
{
	if (offset.isZero())
		return;

	for (Vector3<T>& v : vertices)
		v += offset;
}

template<typename T, typename U>
void TetrahedralMesh<T, U>::transform(const Matrix3<T>& matrix)
{
	if (vertices.empty() || matrix.isIdentity())
		return;

	for (Vector3<T>& v : vertices)
		v *= matrix;
}

template<typename T, typename U>
void TetrahedralMesh<T, U>::transform(const AffineTransform<T>& transformation)
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
bool TetrahedralMesh<T, U>::weldVertices(T tolerance)
{
	if (vertices.size() <= 1)
		return false;

	if (vertices.size() > (std::size_t(std::numeric_limits<U>::max()) + 1))
		throw std::out_of_range("TetrahedralMesh::weldVertices() : vertices.size()");

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

} // namespace templates

#if MATHEMATICS_DOUBLE
using TetrahedralMesh = templates::TetrahedralMesh<double, std::uint32_t>;
#else
using TetrahedralMesh = templates::TetrahedralMesh<float, std::uint32_t>;
#endif

} // namespace mathematics
