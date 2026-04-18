/*
 *	Name: ConvexPolyhedron
 *	Author: Pawel Mrochen
 */

#pragma once

#include <stdexcept>
#include <type_traits>
#include <concepts>
#include <vector>
#include <iterator>
#include <initializer_list>
#include <atomic>
#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <cmath>
#include "Constants.hpp"
#include "Vector3.hpp"
#include "Matrix3.hpp"
#include "AffineTransform.hpp"
#include "HalfSpace.hpp"
#include "Plane.hpp"
#include "AxisAlignedBox.hpp"
#include "OrientedBox.hpp"
#include "Sphere.hpp"
#include "SymmetricFrustum.hpp"

namespace mathematics {
namespace templates {

template<typename T, typename U = unsigned long long>
	requires (std::floating_point<T> && std::integral<U>)
class ConvexPolyhedron
{
public:
	using Real = T;
	using HalfSpaceType = HalfSpace<T>;
	using HalfSpaceVector = std::vector<HalfSpaceType>;
	using MaskType = U;

	ConvexPolyhedron() noexcept : refCount_(), halfSpaces() {}
	template<std::input_iterator<HalfSpace<T>> I, std::sentinel_for<I> S> ConvexPolyhedron(I first, S last);
	ConvexPolyhedron(const HalfSpace<T>* halfSpaces, std::size_t nHalfSpaces);
	explicit ConvexPolyhedron(const std::vector<HalfSpace<T>>& halfSpaces);
	explicit ConvexPolyhedron(std::vector<HalfSpace<T>>&& halfSpaces);
	template<std::input_iterator<Plane<T>> I, std::sentinel_for<I> S> ConvexPolyhedron(I first, S last);
	ConvexPolyhedron(const Plane<T>* planes, std::size_t nPlanes);
	explicit ConvexPolyhedron(const std::vector<Plane<T>>& planes);
	ConvexPolyhedron(std::initializer_list<HalfSpace<T>> halfSpaces) : refCount_(), halfSpaces(halfSpaces) {}
	ConvexPolyhedron(const ConvexPolyhedron& polyhedron);
	ConvexPolyhedron(ConvexPolyhedron&& polyhedron);
	ConvexPolyhedron& operator=(const ConvexPolyhedron& polyhedron);
	ConvexPolyhedron& operator=(ConvexPolyhedron&& polyhedron);

	// Serialization
	template<typename A> void serialize(A& ar) { ar(halfSpaces); }

	// Create
	static ConvexPolyhedron* from(const AxisAlignedBox<T>& box);
	static ConvexPolyhedron* from(const Box<T>& box);
	static ConvexPolyhedron* from(const SymmetricFrustum<T>& frustum);
	//static ConvexPolyhedron* makeConvexHull(const Vector3<T>* points, std::size_t nPoints); // convex hull of point cloud
	//static ConvexPolyhedron* makeConvexHull(const std::vector<Vector3<T>>& points);

	// Clone
	ConvexPolyhedron* clone() const { return new ConvexPolyhedron(*this); }

	// References
	bool hasOwner() const noexcept { return (refCount_ > 0); }
	int getReferenceCount() const noexcept { return refCount_; }
	void acquire() noexcept { refCount_++; }
	void release() { if (--refCount_ <= 0) delete this; }

	// Clear
	static const ConvexPolyhedron* getEmpty();
	bool isEmpty() const noexcept { return halfSpaces.empty(); }
	void clear() { halfSpaces.resize(0); }

	// Assign/append
	void assign(const ConvexPolyhedron* polyhedron);
	void append(const ConvexPolyhedron* polyhedron);

	// Half spaces
	HalfSpaceVector& getHalfSpaces() noexcept { return halfSpaces; }
	const HalfSpaceVector& getHalfSpaces() const noexcept { return halfSpaces; }
	std::size_t getHalfSpaceCount() const noexcept { return halfSpaces.size(); }
	MaskType getHalfSpaceMask() const noexcept { return (MaskType(1) << halfSpaces.size()) - MaskType(1); }
	void setHalfSpaceCount(std::size_t count) { halfSpaces.resize(count); }
	const HalfSpace<T>& getHalfSpace(std::size_t index) const noexcept;
	void setHalfSpace(std::size_t index, const HalfSpace<T>& value); // throw (std::out_of_range)
	void addHalfSpace(const HalfSpace<T>& value) { halfSpaces.push_back(value); }
	template<std::input_iterator<HalfSpace<T>> I, std::sentinel_for<I> S> void addHalfSpaces(I first, S last);
	void addHalfSpaces(const HalfSpace<T>* halfSpaces, std::size_t nHalfSpaces);
	void addHalfSpaces(const std::vector<HalfSpace<T>>& halfSpaces);
	void insertHalfSpace(std::size_t index, const HalfSpace<T>& value); // throw (std::out_of_range)
	void removeHalfSpaceAt(std::size_t index);
	void removeAllHalfSpaces() { halfSpaces.resize(0); }

	// Planes
	template<std::output_iterator<Plane<T>> O> O copyPlanes(O target) const;
	std::vector<Plane<T>> getPlanes() const { return { halfSpaces.begin(), halfSpaces.end() }; }

	// Vertices
	std::vector<Vector3<T>> computeVertices(T tolerance = Constants<T>::TOLERANCE*T(1000)) const;

	// Bounding box
	AxisAlignedBox<T> computeAxisAlignedBoundingBox(T tolerance = Constants<T>::TOLERANCE*T(1000)) const noexcept;

	// Transformation
	void translate(const Vector3<T>& offset) noexcept;
	void transform(const Matrix3<T>& matrix, bool orthogonal = false) noexcept;
	void transform(const AffineTransform<T>& transformation, bool orthogonal = false) noexcept;

	// Intersection
	bool contains(const Vector3<T>& point) const noexcept;
	bool contains(const Vector3<T>& point, T tolerance) const noexcept;
	//bool contains(const AxisAlignedBox<T>& box) const noexcept;
	//bool contains(const Sphere<T>& sphere) const noexcept;
	bool intersects(const AxisAlignedBox<T>& box) const noexcept;
	bool intersects(const AxisAlignedBox<T>& box, U& mask) const noexcept;
	bool intersects(const OrientedBox<T>& box) const noexcept;
	bool intersects(const OrientedBox<T>& box, U& mask) const noexcept;
	bool intersects(const Sphere<T>& sphere) const noexcept;

	HalfSpaceVector halfSpaces;		// normals point outwards

private:
	std::atomic_int refCount_;
};

template<typename T, typename U>
template<std::input_iterator<HalfSpace<T>> I, std::sentinel_for<I> S> 
inline ConvexPolyhedron<T, U>::ConvexPolyhedron(I first, S last) :
	refCount_(), 
	halfSpaces(first, last)
{
}

template<typename T, typename U>
inline ConvexPolyhedron<T, U>::ConvexPolyhedron(const HalfSpace<T>* halfSpaces, std::size_t nHalfSpaces) :
	refCount_(), 
	halfSpaces(halfSpaces, halfSpaces + nHalfSpaces)
{
}

template<typename T, typename U>
inline ConvexPolyhedron<T, U>::ConvexPolyhedron(const std::vector<HalfSpace<T>>& halfSpaces) : 
	refCount_(), 
	halfSpaces(halfSpaces)
{
}

template<typename T, typename U>
inline ConvexPolyhedron<T, U>::ConvexPolyhedron(std::vector<HalfSpace<T>>&& halfSpaces) : 
	refCount_(), 
	halfSpaces(std::move(halfSpaces))
{
}

template<typename T, typename U>
template<std::input_iterator<Plane<T>> I, std::sentinel_for<I> S> 
inline ConvexPolyhedron<T, U>::ConvexPolyhedron(I first, S last) :
	refCount_(), 
	halfSpaces(first, last)
{
}

template<typename T, typename U>
inline ConvexPolyhedron<T, U>::ConvexPolyhedron(const Plane<T>* planes, std::size_t nPlanes) :
	refCount_(), 
	halfSpaces(planes, planes + nPlanes)
{
}

template<typename T, typename U>
inline ConvexPolyhedron<T, U>::ConvexPolyhedron(const std::vector<Plane<T>>& planes) : 
	refCount_(), 
	halfSpaces(planes.begin(), planes.end())
{
}

template<typename T, typename U>
inline ConvexPolyhedron<T, U>::ConvexPolyhedron(const ConvexPolyhedron& polyhedron) : 
	refCount_(), 
	halfSpaces(polyhedron.halfSpaces) 
{
}

template<typename T, typename U>
inline ConvexPolyhedron<T, U>::ConvexPolyhedron(ConvexPolyhedron&& polyhedron) : 
	refCount_(), 
	halfSpaces(std::move(polyhedron.halfSpaces)) 
{
}

template<typename T, typename U>
inline ConvexPolyhedron<T, U>& ConvexPolyhedron<T, U>::operator=(const ConvexPolyhedron<T, U>& polyhedron)
{
	halfSpaces = polyhedron.halfSpaces;
	return *this;
}

template<typename T, typename U>
inline ConvexPolyhedron<T, U>& ConvexPolyhedron<T, U>::operator=(ConvexPolyhedron<T, U>&& mesh)
{
	halfSpaces = std::move(polyhedron.halfSpaces);
	return *this;
}

template<typename T, typename U>
/*static*/ inline ConvexPolyhedron<T, U>* ConvexPolyhedron<T, U>::from(const AxisAlignedBox<T>& box)
{
	return new ConvexPolyhedron<T, U>(box.getHalfSpaces());
}

template<typename T, typename U>
/*static*/ inline ConvexPolyhedron<T, U>* ConvexPolyhedron<T, U>::from(const OrientedBox<T>& box)
{
	return new ConvexPolyhedron<T, U>(box.getHalfSpaces());
}

template<typename T, typename U>
/*static*/ inline ConvexPolyhedron<T, U>* ConvexPolyhedron<T, U>::from(const SymmetricFrustum<T>& frustum)
{
	return new ConvexPolyhedron<T, U>(frustum.getHalfSpaces());
}

template<typename T, typename U>
/*static*/ inline const ConvexPolyhedron<T, U>* ConvexPolyhedron<T, U>::getEmpty()
{
	static const ConvexPolyhedron<T, U> empty;
	return &empty;
}

template<typename T, typename U>
inline void ConvexPolyhedron<T, U>::assign(const ConvexPolyhedron<T, U>* polyhedron)
{
	if (polyhedron)
		halfSpaces.assign(polyhedron->halfSpaces.begin(), polyhedron->halfSpaces.end());
	else
		halfSpaces.resize(0);
}

template<typename T, typename U>
inline void ConvexPolyhedron<T, U>::append(const ConvexPolyhedron<T, U>* polyhedron)
{
	if (polyhedron)
		halfSpaces.insert(halfSpaces.end(), polyhedron->halfSpaces.begin(), polyhedron->halfSpaces.end());
}

template<typename T, typename U>
inline const HalfSpace<T>& ConvexPolyhedron<T, U>::getHalfSpace(std::size_t index) const
{ 
	return (index < halfSpaces.size()) ? halfSpaces[index] : HalfSpace<T>::EMPTY; 
}

template<typename T, typename U>
inline void ConvexPolyhedron<T, U>::setHalfSpace(std::size_t index, const HalfSpace<T>& value) 
{ 
	if (index >= halfSpaces.size()) 
		throw std::out_of_range("ConvexPolyhedron::setHalfSpace() : index"); 
	halfSpaces[index] = value; 
}

template<typename T, typename U>
template<std::input_iterator<HalfSpace<T>> I, std::sentinel_for<I> S> 
inline void ConvexPolyhedron<T, U>::addHalfSpaces(I first, S last)
{
	halfSpaces.insert(halfSpaces.end(), first, last); 
}

template<typename T, typename U>
inline void ConvexPolyhedron<T, U>::addHalfSpaces(const HalfSpace<T>* halfSpaces, std::size_t nHalfSpaces) 
{ 
	this->halfSpaces.insert(this->halfSpaces.end(), halfSpaces, halfSpaces + nHalfSpaces); 
}

template<typename T, typename U>
inline void ConvexPolyhedron<T, U>::addHalfSpaces(const std::vector<HalfSpace<T>>& halfSpaces) 
{ 
	this->halfSpaces.insert(this->halfSpaces.end(), halfSpaces.begin(), halfSpaces.end()); 
}

template<typename T, typename U>
inline void ConvexPolyhedron<T, U>::insertHalfSpace(std::size_t index, const HalfSpace<T>& value) 
{ 
	if (index > halfSpaces.size()) 
		throw std::out_of_range("ConvexPolyhedron::insertHalfSpace() : index"); 
	halfSpaces.insert(halfSpaces.begin() + index, value); 
}

template<typename T, typename U>
inline void ConvexPolyhedron<T, U>::removeHalfSpaceAt(std::size_t index) 
{ 
	if (index < halfSpaces.size()) 
		halfSpaces.erase(halfSpaces.begin() + index); 
}

template<typename T, typename U>
template<std::output_iterator<Plane<T>> O> 
inline O ConvexPolyhedron<T, U>::copyPlanes(O target) const
{
	return std::copy(halfSpaces.begin(), halfSpaces.end(), target);
}

template<typename T, typename U>
void ConvexPolyhedron<T, U>::translate(const Vector3<T>& offset)
{
	if (!offset.isZero())
	{
		for (HalfSpace<T>& h : halfSpaces)
			h.translate(offset);
	}
}

template<typename T, typename U>
void ConvexPolyhedron<T, U>::transform(const Matrix3<T>& matrix, bool orthogonal)
{
	if (!halfSpaces.empty() && !matrix.isIdentity())
	{
		if (orthogonal)
		{
			for (HalfSpace<T>& h : halfSpaces)
				h.transform(matrix, true);
		}
		else
		{
			Matrix3<T> normalMatrix = inverseTranspose(matrix);
			for (HalfSpace<T>& h : halfSpaces)
				h.set(normalize(h.getNormal()*normalMatrix), -dot(h.getNormal(), (h.getNormal()*(-h.getConstant()))*matrix));
		}
	}
}

template<typename T, typename U>
void ConvexPolyhedron<T, U>::transform(const AffineTransform<T>& transformation, bool orthogonal)
{
	if (!halfSpaces.empty() && !transformation.isIdentity())
	{
		if (transformation.getBasis().isIdentity())
		{
			translate(transformation.getOrigin());
		}
		else if (orthogonal)
		{
			for (HalfSpace<T>& h : halfSpaces)
				h.transform(transformation, true);
		}
		else
		{
			Matrix3<T> normalMatrix = inverseTranspose(transformation.getBasis()));
			for (HalfSpace<T>& h : halfSpaces)
				h.set(normalize(h.getNormal()*normalMatrix), -dot(h.getNormal(), transform(h.getNormal()*(-h.getConstant()), transformation)));
		}
	}
}

template<typename T, typename U>
bool ConvexPolyhedron<T, U>::contains(const Vector3<T>& point) const
{
	for (const HalfSpace<T>& h : halfSpaces)
	{
		if (!h.contains(point))
			return false;
	}

	return true;
}

template<typename T, typename U>
bool ConvexPolyhedron<T, U>::contains(const Vector3<T>& point, T tolerance) const
{
	for (const HalfSpace<T>& h : halfSpaces)
	{
		if (h.getSignedDistanceTo(point) > std::max(tolerance, std::fabs(tolerance*h.getConstant())))
			return false;
	}

	return true;
}

template<typename T, typename U>
bool ConvexPolyhedron<T, U>::intersects(const AxisAlignedBox<T>& box) const
{
	for (const HalfSpace<T>& h : halfSpaces)
	{
		if (!box.intersects(h))
			return false;
	}

	return true;
}

template<typename T, typename U>
bool ConvexPolyhedron<T, U>::intersects(const AxisAlignedBox<T>& box, U& mask) const
{
	std::size_t count = halfSpaces.size();
	mask &= (U(1) << count) - U(1);

	if (mask)
	{
		U m = U(1);
		for (std::size_t i = 0; i != count; i++, m <<= 1)
		{
			if (mask & m)
			{
				int side = box.classify(halfSpaces[i]);
				if (side > 0)
					return false;
				if (side < 0)
					mask &= ~m;
			}
		}
	}

	return true;
}

template<typename T, typename U>
bool ConvexPolyhedron<T, U>::intersects(const OrientedBox<T>& box) const
{
	for (const HalfSpace<T>& h : halfSpaces)
	{
		if (!box.intersects(h))
			return false;
	}

	return true;
}

template<typename T, typename U>
bool ConvexPolyhedron<T, U>::intersects(const OrientedBox<T>& box, U& mask) const
{
	std::size_t count = halfSpaces.size();
	mask &= (U(1) << count) - U(1);

	if (mask)
	{
		U m = U(1);
		for (std::size_t i = 0; i != count; i++, m <<= 1)
		{
			if (mask & m)
			{
				int side = box.classify(halfSpaces[i]);
				if (side > 0)
					return false;
				if (side < 0)
					mask &= ~m;
			}
		}
	}

	return true;
}

template<typename T, typename U>
bool ConvexPolyhedron<T, U>::intersects(const Sphere<T>& sphere) const
{
	for (const HalfSpace<T>& h : halfSpaces)
	{
		if (!sphere.intersects(h))
			return false;
	}

	return true;
}

} // namespace templates

#if MATHEMATICS_DOUBLE
using ConvexPolyhedron = templates::ConvexPolyhedron<double>;
#else
using ConvexPolyhedron = templates::ConvexPolyhedron<float>;
#endif

} // namespace mathematics

#include "LU.inl"

namespace mathematics::templates {

template<typename T, typename U>
std::vector<Vector3<T>> ConvexPolyhedron<T, U>::computeVertices(T tolerance) const
{
	// http://www.gamedev.net/page/resources/_/technical/game-programming/shadow-caster-volumes-for-the-culling-of-potential-shadow-casters-r2330

	std::vector<Vector3<T>> vertices;
	for (auto i0 = halfSpaces.begin(), last = halfSpaces.end(); i0 != last; ++i0)
	{ 
		auto i1 = i0; 
        for (++i1; i1 != last; ++i1)
		{
			auto i2 = i1; 
            for (++i2; i2 != last; ++i2)
			{ 
                // Test if the planes intersect at a point
				auto d = lu::decompose(Matrix3<T>(i0->getNormal(), i1->getNormal(), i2->getNormal()));
                if (d.has_value())
				{ 
					Vector3<T> point = lu::solve(d.value(), Vector3<T>(-i0->d, -i1->d, -i2->d));
					if (point.isFinite() && contains(point, tolerance))
					{
						auto range = std::equal_range(vertices.begin(), vertices.end(), point, 
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
							
						if (range.first == range.second)			// point not found
							vertices.insert(range.first, point);	// intersection of 3 planes
					}
				}
			}
		}
    }

	return vertices;
}

template<typename T, typename U>
AxisAlignedBox<T> ConvexPolyhedron<T, U>::computeAxisAlignedBoundingBox(T tolerance) const
{
	// http://www.gamedev.net/page/resources/_/technical/game-programming/shadow-caster-volumes-for-the-culling-of-potential-shadow-casters-r2330

	AxisAlignedBox<T> box;
	for (auto i0 = halfSpaces.begin(), last = halfSpaces.end(); i0 != last; ++i0)
	{ 
		auto i1 = i0; 
        for (++i1; i1 != last; ++i1)
		{ 
			auto i2 = i1; 
            for (++i2; i2 != last; ++i2)
			{ 
                // Test if the planes intersect at a point 
				auto d = lu::decompose(Matrix3<T>(i0->getNormal(), i1->getNormal(), i2->getNormal()));
                if (d.has_value())
				{ 
					Vector3<T> point = lu::solve(d.value(), Vector3<T>(-i0->d, -i1->d, -i2->d));
					if (point.isFinite() && contains(point, tolerance))
					{
						box.minimum.setMinimum(box.minimum, point);
						box.maximum.setMaximum(box.maximum, point);
					}
				}
			}
		}
    }

	return box;
}

} // namespace mathematics::templates
