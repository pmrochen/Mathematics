/*
 *	Name: SymmetricFrustum
 *	Author: Pawel Mrochen
 */

#pragma once

#include <istream>
#include <ostream>
#include <limits>
#include <type_traits>
#include <concepts>
#include <algorithm>
#include <functional>
#include <utility>
#include <vector>
#include <iterator>
#include <cstddef>
#include <cmath>
#include "Scalar.hpp"
#include "Vector3.hpp"
#include "Matrix3.hpp"
#include "AffineTransform.hpp"
#include "HalfSpace.hpp"
#include "AxisAlignedBox.hpp"
#include "OrientedBox.hpp"
#include "Sphere.hpp"
#include "Cone.hpp"

namespace mathematics {
namespace templates {

template<typename T>
	requires std::floating_point<T>
struct SymmetricFrustum
{
	using Real = T;
	using ConstArg = const SymmetricFrustum&;
	using ConstResult = const SymmetricFrustum&;

	SymmetricFrustum() noexcept : origin(), basis(Identity()), halfDims(), depthRange(T()) {}
	explicit SymmetricFrustum(Uninitialized) noexcept : origin(Uninitialized()), basis(Uninitialized()), halfDims(Uninitialized()), depthRange(Uninitialized()) {}
	SymmetricFrustum(const Vector3<T>& origin, const Matrix3<T>& basis, const Vector3<T>& halfDims, const Interval<T>& depthRange) noexcept;
	SymmetricFrustum(const Vector3<T>& origin, const Matrix3<T>& basis, const Vector3<T>& halfDims, T depthMin, T depthMax) noexcept;

	bool operator==(const SymmetricFrustum& frustum) const noexcept;
	bool operator!=(const SymmetricFrustum& frustum) const noexcept { return !(*this == frustum); }

	template<typename A> void serialize(A& ar) { ar(origin, basis, halfDims, depthRange); }

	// Properties
	//bool isEmpty() const noexcept { return halfDims.anyLessThanEqual(Vector2<T>::ZERO) || (depthRange.minimum >= depthRange.maximum); }
	bool approxEquals(const SymmetricFrustum& frustum) const noexcept;
	bool approxEquals(const SymmetricFrustum& frustum, T tolerance) const noexcept;
	bool isFinite() const noexcept { return origin.isFinite() && basis.isFinite() && halfDims.isFinite() && depthRange.isFinite(); }
	SymmetricFrustum& set(const Vector3<T>& origin, const Matrix3<T>& basis, const Vector3<T>& halfDims, const Interval<T>& depthRange) noexcept;
	const Vector3<T>& getOrigin() const noexcept { return origin; }
	void setOrigin(const Vector3<T>& origin) noexcept { this->origin = origin; }
	const Matrix3<T>& getBasis() const noexcept { return basis; }
	void setBasis(const Matrix3<T>& basis) noexcept { this->basis = basis; }
	const Vector2<T>& getHalfDimensions() const noexcept { return halfDims; }
	void setHalfDimensions(const Vector2<T>& halfDims) noexcept { this->halfDims = halfDims; }
	Vector2<T> getDimensions() const noexcept { return halfDims*T(2); }
	void setDimensions(const Vector2<T>& dimensions) noexcept { halfDims = dimensions*T(0.5); }
	Vector2<T> getBaseHalfDimensions() const noexcept { return halfDims*(depthRange.maximum/depthRange.minimum); }
	void setBaseHalfDimensions(const Vector2<T>& halfDims) noexcept { this->halfDims = halfDims*(depthRange.minimum/depthRange.maximum); }
	Vector2<T> getBaseDimensions() const noexcept { return halfDims*(T(2)*depthRange.maximum/depthRange.minimum); }
	void setBaseDimensions(const Vector2<T>& dimensions) noexcept { halfDims = dimensions*(T(0.5)*depthRange.minimum/depthRange.maximum); }
	const Interval<T>& getDepthRange() const noexcept { return depthRange; }
	void setDepthRange(const Interval<T>& depthRange) noexcept { this->depthRange = depthRange; }
	T getMinDepth() const noexcept { return depthRange.minimum; }
	void setMinDepth(T depthMin) noexcept { depthRange.minimum = depthMin; }
	T getMaxDepth() const noexcept { return depthRange.maximum; }
	void setMaxDepth(T depthMax) noexcept { depthRange.maximum = depthMax; }
	T getDepthRatio() const noexcept { return depthRange.maximum/depthRange.minimum; }

	// Vertices
	template<std::output_iterator<Vector3<T>> O> O copyVertices(O target) const;
	std::vector<Vector3<T>> getVertices() const;

	// Primitives
	template<std::integral U> std::pair<const U*, const U*> getPrimitives(int nVerticesPerPrimitive) const noexcept; // #TODO return range
	std::size_t getPrimitiveCount(int nVerticesPerPrimitive) const noexcept;

	// Half spaces
	template<std::output_iterator<HalfSpace<T>> O> O copyHalfSpaces(O target) const;	// 0-left, 1-right, 2-bottom, 3-top, 4-near, 5-far
	std::vector<HalfSpace<T>> getHalfSpaces() const;									// normals point outwards
	template<typename F> bool enumerateHalfSpaces(F&& f) const noexcept;

	// Circumscribed box, sphere, cone
	OrientedBox<T> getCircumscribedBox() const noexcept;
	Sphere<T> getCircumscribedSphere() const noexcept;
	Cone<T> getCircumscribedCone() const noexcept;

	// Transformation
	SymmetricFrustum& translate(const Vector3<T>& offset) noexcept { origin += offset; return *this; }
	SymmetricFrustum& orthonormalize() noexcept;

	// Closest point
	Vector3<T> getClosestPoint(const Vector3<T>& point) const noexcept;
	T getDistanceTo(const Vector3<T>& point) const noexcept;

	// Containment and intersection
	bool contains(const Vector3<T>& point) const noexcept;
	bool intersects(const AxisAlignedBox<T>& box) const noexcept;
	bool intersects(const OrientedBox<T>& box) const noexcept;
	bool intersects(const Sphere<T>& sphere) const noexcept;

	Vector3<T> origin;
	Matrix3<T> basis;		// unit length axes
	Vector2<T> halfDims;	// half dimensions at depthRange.minimum
	Interval<T> depthRange;
};

template<typename T>
inline SymmetricFrustum<T>::SymmetricFrustum(const Vector3<T>& origin, const Matrix3<T>& basis, const Vector3<T>& halfDims, 
	const Interval<T>& depthRange) : 
	origin(origin), 
	basis(basis), 
	halfDims(halfDims),
	depthRange(depthRange)
{
}

template<typename T>
inline SymmetricFrustum<T>::SymmetricFrustum(const Vector3<T>& origin, const Matrix3<T>& basis, const Vector3<T>& halfDims, 
	T depthMin, T depthMax) : 
	origin(origin), 
	basis(basis), 
	halfDims(halfDims),
	depthRange(depthMin, depthMax)
{
}

template<typename T>
inline bool SymmetricFrustum<T>::operator==(const SymmetricFrustum<T>& frustum) const
{ 
	return (origin == frustum.origin) && (basis == frustum.basis) && (halfDims == frustum.halfDims) &&
		(depthRange == frustum.depthRange);
}

template<typename C, typename T, typename U>
	requires std::floating_point<U>
inline std::basic_istream<C, T>& operator>>(std::basic_istream<C, T>& s, SymmetricFrustum<U>& frustum)
{ 
	return s >> frustum.origin >> std::ws >> frustum.basis >> std::ws >> frustum.halfDims >> std::ws >> frustum.depthRange;
}

template<typename C, typename T, typename U>
	requires std::floating_point<U>
inline std::basic_ostream<C, T>& operator<<(std::basic_ostream<C, T>& s, const SymmetricFrustum<U>& frustum)
{ 
	constexpr C WS(0x20);
	return s << frustum.origin << WS << frustum.basis << WS << frustum.halfDims << WS << frustum.depthRange;
}

template<typename T>
inline bool SymmetricFrustum<T>::approxEquals(const SymmetricFrustum<T>& frustum) const
{
	return origin.approxEquals(frustum.origin) && basis.approxEquals(frustum.basis) && halfDims.approxEquals(frustum.halfDims) &&
		depthRange.approxEquals(frustum.depthRange);
}

template<typename T>
inline bool SymmetricFrustum<T>::approxEquals(const SymmetricFrustum<T>& frustum, T tolerance) const
{
	return origin.approxEquals(frustum.origin, tolerance) && basis.approxEquals(frustum.basis, tolerance) && 
		halfDims.approxEquals(frustum.halfDims, tolerance) && depthRange.approxEquals(frustum.depthRange, tolerance);
}

template<typename T>
inline SymmetricFrustum<T>& SymmetricFrustum<T>::set(const Vector3<T>& origin, const Matrix3<T>& basis, const Vector3<T>& halfDims,
	const Interval<T>& depthRange) 
{ 
	this->origin = origin; 
	this->basis = basis; 
	this->halfDims = halfDims;
	this->depthRange = depthRange;
	return *this;
}

template<typename T>
template<std::output_iterator<Vector3<T>> O>
inline O SymmetricFrustum<T>::copyVertices(O target) const
{
	AffineTransform<T> m(basis, origin);
	T depthRatio = depthRange.maximum/depthRange.minimum;
	*target++ = transform(Vector3<T>(-halfDims, depthRange.minimum), m);
	*target++ = transform(Vector3<T>(halfDims.x, -halfDims.y, depthRange.minimum), m);
	*target++ = transform(Vector3<T>(-halfDims.x, halfDims.y, depthRange.minimum), m);
	*target++ = transform(Vector3<T>(halfDims, depthRange.minimum), m);
	*target++ = transform(Vector3<T>(-halfDims*depthRatio, depthRange.maximum), m);
	*target++ = transform(Vector3<T>(halfDims.x*depthRatio, -halfDims.y*depthRatio, depthRange.maximum), m);
	*target++ = transform(Vector3<T>(-halfDims.x*depthRatio, halfDims.y*depthRatio, depthRange.maximum), m);
	*target++ = transform(Vector3<T>(halfDims*depthRatio, depthRange.maximum), m);
	return target;
}

template<typename T>
inline std::vector<Vector3<T>> SymmetricFrustum<T>::getVertices() const
{
	AffineTransform<T> m(basis, origin);
	T depthRatio = depthRange.maximum/depthRange.minimum;
	return { transform(Vector3<T>(-halfDims, depthRange.minimum), m),
		transform(Vector3<T>(halfDims.x, -halfDims.y, depthRange.minimum), m),
		transform(Vector3<T>(-halfDims.x, halfDims.y, depthRange.minimum), m),
		transform(Vector3<T>(halfDims, depthRange.minimum), m),
		transform(Vector3<T>(-halfDims*depthRatio, depthRange.maximum), m),
		transform(Vector3<T>(halfDims.x*depthRatio, -halfDims.y*depthRatio, depthRange.maximum), m),
		transform(Vector3<T>(-halfDims.x*depthRatio, halfDims.y*depthRatio, depthRange.maximum), m),
		transform(Vector3<T>(halfDims*depthRatio, depthRange.maximum), m) };
}

template<typename T>
template<std::integral U>
std::pair<const U*, const U*> SymmetricFrustum<T>::getPrimitives(int nVerticesPerPrimitive) const
{
	static const U edges[24] = { 0, 2, 2, 3, 3, 1, 1, 0, 3, 7, 5, 1, 6, 2, 0, 4, 5, 7, 7, 6, 6, 4, 4, 5 };
	static const U triangles[36] = { 0, 2, 1, 3, 1, 2, 1, 3, 5, 7, 5, 3, 5, 7, 4, 6, 4, 7, 4, 6, 0, 2, 0, 6, 2, 6, 3, 7, 3, 6, 4, 0, 5, 1, 5, 0 };
	static const U quads[24] = { 0, 2, 3, 1, 1, 3, 7, 5, 5, 7, 6, 4, 4, 6, 2, 0, 2, 6, 7, 3, 4, 0, 1, 5 };

	switch (nVerticesPerPrimitive)
	{
		case 2:
			return { edges, edges + 24 };
		case 3:
			return { triangles, triangles + 36 };
		case 4:
			return { quads, quads + 24 };
		default:
			return { nullptr, nullptr };
	}
}

template<typename T>
inline std::size_t SymmetricFrustum<T>::getPrimitiveCount(int nVerticesPerPrimitive) const
{
	switch (nVerticesPerPrimitive)
	{
		case 2:
			return 12;
		case 3:
			return 12;
		case 4:
			return 6;
		default:
			return 0;
	}
}

template<typename T>
template<std::output_iterator<HalfSpace<T>> O>
inline O SymmetricFrustum<T>::copyHalfSpaces(O target) const
{
	AffineTransform<T> m(basis, origin);
	Vector3<T> lb = transform(Vector3<T>(-halfDims, depthRange.minimum), m);
	Vector3<T> rb = transform(Vector3<T>(halfDims.x, -halfDims.y, depthRange.minimum), m);
	Vector3<T> lt = transform(Vector3<T>(-halfDims.x, halfDims.y, depthRange.minimum), m);
	Vector3<T> rt = transform(Vector3<T>(halfDims, depthRange.minimum), m);

	bool flip = (basis.getDeterminant() < T(0));
	*target++ = flip ? HalfSpace<T>(origin, lt, lb) : HalfSpace<T>(origin, lb, lt);
	*target++ = flip ? HalfSpace<T>(origin, rb, rt) : HalfSpace<T>(origin, rt, rb);
	*target++ = flip ? HalfSpace<T>(origin, lb, rb) : HalfSpace<T>(origin, rb, lb);
	*target++ = flip ? HalfSpace<T>(origin, rt, lt) : HalfSpace<T>(origin, lt, rt);
	*target++ = HalfSpace<T>(-basis[2], depthRange.minimum*basis[2] + origin);
	if (depthRange.maximum < std::numeric_limits<T>::max())
		*target++ = HalfSpace<T>(basis[2], depthRange.maximum*basis[2] + origin);

	return target;
}

template<typename T>
inline std::vector<HalfSpace<T>> SymmetricFrustum<T>::getHalfSpaces() const
{
	AffineTransform<T> m(basis, origin);
	Vector3<T> lb = transform(Vector3<T>(-halfDims, depthRange.minimum), m);
	Vector3<T> rb = transform(Vector3<T>(halfDims.x, -halfDims.y, depthRange.minimum), m);
	Vector3<T> lt = transform(Vector3<T>(-halfDims.x, halfDims.y, depthRange.minimum), m);
	Vector3<T> rt = transform(Vector3<T>(halfDims, depthRange.minimum), m);

	bool flip = (basis.getDeterminant() < T(0));
	if (depthRange.maximum < std::numeric_limits<T>::max())
	{
		return { flip ? HalfSpace<T>(origin, lt, lb) : HalfSpace<T>(origin, lb, lt),
			flip ? HalfSpace<T>(origin, rb, rt) : HalfSpace<T>(origin, rt, rb),
			flip ? HalfSpace<T>(origin, lb, rb) : HalfSpace<T>(origin, rb, lb),
			flip ? HalfSpace<T>(origin, rt, lt) : HalfSpace<T>(origin, lt, rt),
			HalfSpace<T>(-basis[2], depthRange.minimum*basis[2] + origin),
			HalfSpace<T>(basis[2], depthRange.maximum*basis[2] + origin) };
	}
	else
	{
		return { flip ? HalfSpace<T>(origin, lt, lb) : HalfSpace<T>(origin, lb, lt),
			flip ? HalfSpace<T>(origin, rb, rt) : HalfSpace<T>(origin, rt, rb),
			flip ? HalfSpace<T>(origin, lb, rb) : HalfSpace<T>(origin, rb, lb),
			flip ? HalfSpace<T>(origin, rt, lt) : HalfSpace<T>(origin, lt, rt),
			HalfSpace<T>(-basis[2], depthRange.minimum*basis[2] + origin) };
	}
}

template<typename T>
template<typename F> 
inline bool SymmetricFrustum<T>::enumerateHalfSpaces(F&& f) const
{
	AffineTransform<T> m(basis, origin);
	Vector3<T> lb = transform(Vector3<T>(-halfDims, depthRange.minimum), m);
	Vector3<T> rb = transform(Vector3<T>(halfDims.x, -halfDims.y, depthRange.minimum), m);
	Vector3<T> lt = transform(Vector3<T>(-halfDims.x, halfDims.y, depthRange.minimum), m);
	Vector3<T> rt = transform(Vector3<T>(halfDims, depthRange.minimum), m);

	bool flip = (basis.getDeterminant() < T(0));
	if (!f(flip ? HalfSpace<T>(origin, lt, lb) : HalfSpace<T>(origin, lb, lt)))
		return false;
	if (!f(flip ? HalfSpace<T>(origin, rb, rt) : HalfSpace<T>(origin, rt, rb)))
		return false;
	if (!f(flip ? HalfSpace<T>(origin, lb, rb) : HalfSpace<T>(origin, rb, lb)))
		return false;
	if (!f(flip ? HalfSpace<T>(origin, rt, lt) : HalfSpace<T>(origin, lt, rt)))
		return false;
	if (!f(HalfSpace<T>(-basis[2], depthRange.minimum*basis[2] + origin)))
		return false;
	if ((depthRange.maximum < std::numeric_limits<T>::max()) && !f(HalfSpace<T>(basis[2], depthRange.maximum*basis[2] + origin)))
		return false;

	return true;
}

template<typename T>
inline OrientedBox<T> SymmetricFrustum<T>::getCircumscribedBox() const
{
	Vector2<T> baseHalfDims = halfDims*(depthRange.maximum/depthRange.minimum);
	return OrientedBox<T>(origin + ((depthRange.minimum + depthRange.maximum)*T(0.5))*basis[2], basis, 
		Vector3<T>(baseHalfDims, (depthRange.maximum - depthRange.minimum)*T(0.5)));
}

template<typename T>
inline Sphere<T> SymmetricFrustum<T>::getCircumscribedSphere() const
{
	Vector2<T> baseHalfDims = halfDims*(depthRange.maximum/depthRange.minimum);
	T coneRadiusSq = baseHalfDims.getMagnitudeSquared();
	T depthMaxSq = square(depthRange.maximum);
	if (depthMaxSq > coneRadiusSq)
	{
		T sphereRadius = (coneRadiusSq + depthMaxSq)/(T(2)*depthRange.maximum);
		return Sphere<T>(origin + sphereRadius*basis[2], sphereRadius);
	}
	else
	{
		return Sphere<T>(origin + depthRange.maximum*basis[2], std::sqrt(coneRadiusSq));
	}
}

template<typename T>
inline Cone<T> SymmetricFrustum<T>::getCircumscribedCone() const
{
	Vector2<T> baseHalfDims = halfDims*(depthRange.maximum/depthRange.minimum);
	return Cone<T>(origin, basis[2], depthRange.maximum, baseHalfDims.getMagnitude());
}

template<typename T>
inline SymmetricFrustum<T> SymmetricFrustum<T>::orthonormalize()
{
	halfDims *= Vector2<T>(basis[0].getMagnitude(), basis[1].getMagnitude());
	//halfDims.x *= basis[0].getMagnitude();
	//halfDims.y *= basis[1].getMagnitude();
	depthRange.scale(basis[2].getMagnitude());
	basis.orthonormalize();
	return *this;
}

template<typename T>
inline bool SymmetricFrustum<T>::contains(const Vector3<T>& point) const
{
	return enumerateHalfSpaces([&point](const HalfSpace<T>& h) { return h.contains(point); });
}

} // namespace templates

#if MATHEMATICS_DOUBLE
using SymmetricFrustum = templates::SymmetricFrustum<double>;
using SymmetricFrustumArg = templates::SymmetricFrustum<double>::ConstArg;
using SymmetricFrustumResult = templates::SymmetricFrustum<double>::ConstResult;
#else
using SymmetricFrustum = templates::SymmetricFrustum<float>;
using SymmetricFrustumArg = templates::SymmetricFrustum<float>::ConstArg;
using SymmetricFrustumResult = templates::SymmetricFrustum<float>::ConstResult;
#endif

} // namespace mathematics

namespace std {

template<typename T>
struct hash;

template<typename T>
struct hash<::mathematics::templates::SymmetricFrustum<T>>
{
	size_t operator()(const ::mathematics::templates::SymmetricFrustum<T>& frustum) const noexcept
	{
		size_t seed = hash<typename ::mathematics::templates::Vector3<T>>()(frustum.origin) + 0x9e3779b9;
		seed ^= hash<typename ::mathematics::templates::Matrix3<T>>()(frustum.basis) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
		seed ^= hash<typename ::mathematics::templates::Vector2<T>>()(frustum.halfDims) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
		seed ^= hash<typename ::mathematics::templates::Interval<T>>()(frustum.depthRange) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
		return seed;
	}
};

} // namespace std

#include "Distances.inl"
#include "Intersections.inl"

namespace mathematics::templates {

template<typename T>
inline Vector3<T> SymmetricFrustum<T>::getClosestPoint(const Vector3<T>& point) const
{
	Vector3<T> closestPoint(Uninitialized());
	distances::getPointSymmetricFrustumSquared(point, origin, basis, halfDims, depthRange.minimum, depthRange.maximum, &closestPoint);
	return closestPoint;
}

template<typename T>
inline T SymmetricFrustum<T>::getDistanceTo(const Vector3<T>& point) const
{
	return distances::getPointSymmetricFrustum(point, origin, basis, halfDims, depthRange.minimum, depthRange.maximum);
}

template<typename T>
inline bool SymmetricFrustum<T>::intersects(const AxisAlignedBox<T>& box) const
{
	return intersections::testAxisAlignedBoxSymmetricFrustum(box.getCenter(), box.getHalfDimensions(), origin, basis, halfDims, 
		depthRange.minimum, depthRange.maximum);
}

template<typename T>
inline bool SymmetricFrustum<T>::intersects(const OrientedBox<T>& box) const
{
	return intersections::testOrientedBoxSymmetricFrustum(box.origin, box.basis, box.halfDims, origin, basis, halfDims, 
		depthRange.minimum, depthRange.maximum);
}

template<typename T>
inline bool SymmetricFrustum<T>::intersects(const Sphere<T>& sphere) const
{
	return (distances::getPointSymmetricFrustumSquared(sphere.center, origin, basis, halfDims, 
		depthRange.minimum, depthRange.maximum) <= sphere.radius*sphere.radius);
}

} // namespace mathematics::templates
