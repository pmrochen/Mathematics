/*
 *	Name: Plane
 *	Author: Pawel Mrochen
 */

#pragma once

#include <istream>
#include <ostream>
#include <limits>
#include <type_traits>
#include <concepts>
#include <utility>
#include <tuple>
#include <optional>
#include <algorithm>
#include <cstddef>
#include <cmath>
#include <Simd/Intrinsics.hpp>
#include "Constants.hpp"
#include "Vector3.hpp"
#include "Matrix3.hpp"
#include "AffineTransform.hpp"
#include "Line3.hpp"
#include "Ray3.hpp"
#include "Segment3.hpp"
#include "HalfSpace.hpp"

namespace core::mathematics {
namespace templates {

template<typename T>
	requires std::floating_point<T>
struct Plane
{
	using Real = T;
	using ConstArg = const Plane&;
	using ConstResult = const Plane&;

	constexpr Plane() noexcept : a(), b(), c(), d() {}
	explicit Plane(Uninitialized) noexcept {}
	constexpr Plane(T a, T b, T c, T d) noexcept : a(a), b(b), c(c), d(d) {}
	constexpr Plane(const Vector3<T>& normal, T constant) noexcept : x(normal.x), y(normal.y), z(normal.z), w(constant) {}
	Plane(const Vector3<T>& normal, const Vector3<T>& point) noexcept : a(normal.x), b(normal.y), c(normal.z), d(-dot(normal, point)) {}
	Plane(const Vector3<T>& p0, const Vector3<T>& p1, const Vector3<T>& p2) noexcept;
	explicit Plane(const HalfSpace<T>& h) noexcept : a(h.a), b(h.b), c(h.c), d(h.d) {}
	explicit Plane(const std::tuple<T, T, T, T>& t) noexcept : a(std::get<0>(t)), b(std::get<1>(t)), c(std::get<2>(t)), d(std::get<3>(t)) {}
	template<typename U> explicit Plane(const std::tuple<U, U, U, U>& t) noexcept : a(T(std::get<0>(t))), b(T(std::get<1>(t))), c(T(std::get<2>(t))), d(T(std::get<3>(t))) {}
	//explicit Plane(const T* p) noexcept : a(p[0]), b(p[1]), c(p[2]), d(p[3]) {}

	//explicit operator std::tuple<T, T, T, T>() noexcept { return std::tuple<T, T, T, T>(a, b, c, d); }
	//template<typename U> explicit operator std::tuple<U, U, U, U>() noexcept { return std::tuple<U, U, U, U>(U(a), U(b), U(c), U(d)); }
	//explicit operator T*() noexcept { return &a; }
	//explicit operator const T*() const noexcept { return &a; }
	//T& operator[](int i) noexcept { return (&a)[i]; }
	//const T& operator[](int i) const noexcept { return (&a)[i]; }

	bool operator==(const Plane& p) const noexcept { return (a == p.a) && (b == p.b) && (c == p.c) && (d == p.d); }
	bool operator!=(const Plane& p) const noexcept { return !(*this == p); }

	template<typename A> void load(A& ar) { ar(a, b, c, d); }
	template<typename A> void save(A& ar) const { ar(a, b, c, d); }

	template<std::size_t I> T& get() noexcept;
	template<std::size_t I> const T& get() const noexcept;

	const HalfSpace<T>& asHalfSpace() const noexcept { return reinterpret_cast<const HalfSpace<T>&>(*this); }

	// Properties
	bool isEmpty() const noexcept { return (a == T()) && (b == T()) && (c == T()) && (d == T()); }
	bool isApproxEqual(const Plane& p) const noexcept;
	bool isApproxEqual(const Plane& p, T tolerance) const noexcept;
	bool isFinite() const noexcept { return std::isfinite(a) && std::isfinite(b) && std::isfinite(c) && std::isfinite(d); }
	const Vector3<T>& getNormal() const noexcept { return reinterpret_cast<const Vector3&>(*this); }
	void setNormal(const Vector3<T>& normal) noexcept { a = normal.x; b = normal.y; c = normal.z; }
	T getConstant() const noexcept { return d; }
	void setConstant(T constant) noexcept { d = constant; }
	Plane& makeEmpty() noexcept { a = T(); b = T(); c = T(); d = T(); return *this; }
	Plane& set(const Vector3<T>& normal, T constant) noexcept { a = normal.x; b = normal.y; c = normal.z; d = constant; return *this; }
	Plane& set(T a, T b, T c, T d) noexcept { this->a = a; this->b = b; this->c = c; this->d = d; return *this; }

	// Transformation
	Plane& translate(const Vector3<T>& offset) noexcept;
	Plane& transform(const Matrix3<T>& matrix, bool orthogonal = false) noexcept;
	Plane& transform(const AffineTransform<T>& transformation, bool orthogonal = false) noexcept;
	Plane& flip() noexcept { a = -a; b = -b; c = -c; d = -d; return *this; }
	Plane& normalize() noexcept;

	// Reflection of a point on a normalized plane
	//Vector3<T> reflect(const Vector3<T>& point) const noexcept; // #TODO -> Vector3

	// Distances
	T getDistance(const Vector3<T>& point) const noexcept { return std::fabs(dot(getNormal(), point) + d); }	// normalized plane
	template<typename U> T getDistance(const Vector3<T>& point) const noexcept;
	T getSignedDistance(const Vector3<T>& point) const noexcept { return (dot(getNormal(), point) + d); }		// normalized plane
	template<typename U> T getSignedDistance(const Vector3<T>& point) const noexcept;

	// Intersection
	bool contains(const Vector3<T>& point) const noexcept;
	bool testIntersection(const Line3<T>& line) const noexcept;
	bool testIntersection(const Ray3<T>& ray) const noexcept { return findIntersection(ray).has_value(); }
	bool testIntersection(const Segment3<T>& segment) const noexcept { return findIntersection(segment).has_value(); }
	//bool testIntersection(const Triangle<T>& triangle) const noexcept; // #TODO Check if all vertices are on the same side
	std::optional<T> findIntersection(const Line3<T>& line) const noexcept;
	std::optional<T> findIntersection(const Ray3<T>& ray) const noexcept;
	std::optional<T> findIntersection(const Segment3<T>& segment) const noexcept;
	template<typename U> std::optional<U> findIntersection(const Line3<T>& line) const noexcept;
	template<typename U> std::optional<U> findIntersection(const Ray3<T>& ray) const noexcept;
	template<typename U> std::optional<U> findIntersection(const Segment3<T>& segment) const noexcept;

	//static const Plane& getEmpty() noexcept { return EMPTY; }

	// Least-squares fit of a plane
	//static Plane computeBestFit(const Vector3<T>* points, size_t nPoints); // #TODO
	//static Plane computeBestFit(const std::vector<Vector3<T>>& points);

	static const Plane EMPTY;

	T a, b, c, d;
};

template<typename T> const Plane<T> Plane<T>::EMPTY{};

#if SIMD_HAS_FLOAT4

template<>
struct Plane<float>
{
	using Real = float;
	using ConstArg = const Plane;
	using ConstResult = const Plane;

	/*constexpr*/ Plane() noexcept : abcd(simd::zero<simd::float4>()) {}
	explicit Plane(Uninitialized) noexcept {}
	/*constexpr*/ Plane(float a, float b, float c, float d) noexcept : abcd(simd::set4(a, b, c, d)) {}
	/*constexpr*/ Plane(const Vector3<float>& normal, float constant) noexcept : abcd(simd::insert<simd::W>(constant, normal)) {}
	Plane(const Vector3<float>& normal, const Vector3<float>& point) noexcept : abcd(simd::insert<simd::W>(-dot(normal, point), normal)) {}
	Plane(const Vector3<float>& p0, const Vector3<float>& p1, const Vector3<float>& p2) noexcept;
	explicit Plane(const HalfSpace<float>& h) noexcept : abcd(h.abcd) {}
	explicit Plane(const std::tuple<float, float, float, float>& t) noexcept : abcd(simd::set4(std::get<0>(t), std::get<1>(t), std::get<2>(t), std::get<3>(t))) {}
	template<typename U> explicit Plane(const std::tuple<U, U, U, U>& t) noexcept : abcd(simd::set4((float)std::get<0>(t), (float)std::get<1>(t), (float)std::get<2>(t), (float)std::get<3>(t))) {}
	//explicit Plane(const float* p) noexcept : abcd(simd::load4(p)) {}
	explicit Plane(simd::float4 p) noexcept : abcd(p) {}
	Plane(const Plane& p) noexcept : abcd(p.abcd) {}
	Plane& operator=(const Plane& p) noexcept { abcd = p.abcd; return *this; }

	operator simd::float4() const noexcept { return abcd; }
	//explicit operator std::tuple<float, float, float, float>() noexcept { return std::tuple<float, float, float, float>(a, b, c, d); }
	//template<typename U> explicit operator std::tuple<U, U, U, U>() noexcept { return std::tuple<U, U, U, U>(U(a), U(b), U(c), U(d)); }
	//explicit operator float* () noexcept { return &a; }
	//explicit operator const float* () const noexcept { return &a; }
	//float& operator[](int i) noexcept { return (&a)[i]; }
	//const float& operator[](int i) const noexcept { return (&a)[i]; }

	bool operator==(const Plane& p) const noexcept { return simd::all4(simd::equal(abcd, p)); }
	bool operator!=(const Plane& p) const noexcept { return !(*this == p); }

	template<typename A> void load(A& ar);
	template<typename A> void save(A& ar) const { ar(a, b, c, d); }

	template<std::size_t I> float& get() noexcept;
	template<std::size_t I> const float& get() const noexcept;
	template<typename U> U& get() noexcept;				// intentionally undefined
	template<typename U> const U& get() const noexcept;	// intentionally undefined
	template<> simd::float4& get() noexcept { return abcd; }
	template<> const simd::float4& get() const noexcept { return abcd; }

	const HalfSpace<float> asHalfSpace() const noexcept { return HalfSpace<float>(abcd); }

	// Properties
	bool isEmpty() const noexcept { return simd::all4(simd::equal(abcd, simd::zero<simd::float4>())); }
	bool isApproxEqual(const Plane& p) const noexcept { simd::all4(simd::lessThan(simd::abs4(simd::sub4(abcd, p)), simd::set4(Constants<float>::TOLERANCE))); }
	bool isApproxEqual(const Plane& p, float tolerance) const noexcept { simd::all4(simd::lessThan(simd::abs4(simd::sub4(abcd, p)), simd::set4(tolerance))); }
	bool isFinite() const noexcept { return simd::all4(simd::isFinite(abcd)); }
#if MATHEMATICS_SIMD_EXPAND_LAST
	const Vector3<float> getNormal() const noexcept { return Vector3<float>(simd::xyzz(abcd)); }
#else
	const Vector3<float> getNormal() const noexcept { return Vector3<float>(simd::cutoff3(abcd)); }
#endif
	void setNormal(const Vector3<float>& normal) noexcept { abcd = simd::insert3(normal, abcd); }
	float getConstant() const noexcept { return simd::extract<simd::W>(abcd); }
	void setConstant(float constant) noexcept { abcd = simd::insert<simd::W>(constant, abcd); }
	Plane& makeEmpty() noexcept { abcd = simd::zero<simd::float4>(); return *this; }
	Plane& set(const Vector3<float>& normal, float constant) noexcept { abcd = simd::insert<simd::W>(constant, normal); return *this; }
	Plane& set(float a, float b, float c, float d) noexcept { abcd = simd::set4(a, b, c, d); return *this; }

	// Transformation
	Plane& translate(const Vector3<float>& offset) noexcept;
	Plane& transform(const Matrix3<float>& matrix, bool orthogonal = false) noexcept;
	Plane& transform(const AffineTransform<float>& transformation, bool orthogonal = false) noexcept;
	Plane& flip() noexcept { abcd = simd::neg4(abcd); return *this; }
	Plane& normalize() noexcept;

	// Reflection of a point on a normalized plane
	//Vector3<float> reflect(const Vector3<float>& point) const noexcept; // #TODO -> Vector3

	// Distances
	float getDistance(const Vector3<float>& point) const noexcept { return std::fabs(dot(getNormal(), point) + d); }	// normalized plane
	template<typename U> float getDistance(const Vector3<float>& point) const noexcept;
	float getSignedDistance(const Vector3<float>& point) const noexcept { return (dot(getNormal(), point) + d); }		// normalized plane
	template<typename U> float getSignedDistance(const Vector3<float>& point) const noexcept;

	// Intersection
	bool contains(const Vector3<float>& point) const noexcept;
	bool testIntersection(const Line3<float>& line) const noexcept;
	bool testIntersection(const Ray3<float>& ray) const noexcept { return findIntersection(ray).has_value(); }
	bool testIntersection(const Segment3<float>& segment) const noexcept { return findIntersection(segment).has_value(); }
	//bool testIntersection(const Triangle<float>& triangle) const noexcept; // #TODO Check if all vertices are on the same side
	std::optional<float> findIntersection(const Line3<float>& line) const noexcept;
	std::optional<float> findIntersection(const Ray3<float>& ray) const noexcept;
	std::optional<float> findIntersection(const Segment3<float>& segment) const noexcept;
	template<typename U> std::optional<U> findIntersection(const Line3<float>& line) const noexcept;
	template<typename U> std::optional<U> findIntersection(const Ray3<float>& ray) const noexcept;
	template<typename U> std::optional<U> findIntersection(const Segment3<float>& segment) const noexcept;

	//static const Plane& getEmpty() noexcept { return EMPTY; }

	// Least-squares fit of a plane
	//static Plane computeBestFit(const Vector3<float>* points, size_t nPoints); // #TODO
	//static Plane computeBestFit(const std::vector<Vector3<float>>& points);

	static const Plane EMPTY;

	union
	{
		simd::float4 abcd;
		struct { float a, b, c, d; };
	};
};

const Plane<float> Plane<float>::EMPTY{};

#endif /* SIMD_HAS_FLOAT4 */

template<typename T>
inline Plane<T>::Plane<T>(const Vector3<T>& p0, const Vector3<T>& p1, const Vector3<T>& p2)
{
	Vector3<T> normal(cross(p1 - p0, p2 - p0));
	normal.normalize();
	set(normal, -dot(normal, p0));
}

template<typename C, typename T, typename U>
	requires std::floating_point<U>
inline std::basic_istream<C, T>& operator>>(std::basic_istream<C, T>& s, Plane<U>& p)
{ 
	return s >> p.a >> std::ws >> p.b >> std::ws >> p.c >> std::ws >> p.d; 
}

template<typename C, typename T, typename U>
	requires std::floating_point<U>
inline std::basic_ostream<C, T>& operator<<(std::basic_ostream<C, T>& s, const Plane<U>& p)
{ 
	constexpr C WS(0x20);
	return s << p.a << WS << p.b << WS << p.c << WS << p.d;
}

template<typename T>
template<std::size_t I>
inline T& Plane<T>::get()
{
	if constexpr (I == 0)
		return a;
	else if constexpr (I == 1)
		return b;
	else if constexpr (I == 2)
		return c;
	else if constexpr (I == 3)
		return d;
	static_assert(false);
}

template<typename T>
template<std::size_t I>
inline const T& Plane<T>::get() const
{
	if constexpr (I == 0)
		return a;
	else if constexpr (I == 1)
		return b;
	else if constexpr (I == 2)
		return c;
	else if constexpr (I == 3)
		return d;
	static_assert(false);
}

template<typename T>
inline bool Plane<T>::isApproxEqual(const Plane<T>& p) const
{ 
	return (std::fabs(p.a - a) < Constants<T>::TOLERANCE) && (std::fabs(p.b - b) < Constants<T>::TOLERANCE) && 
		(std::fabs(p.c - c) < Constants<T>::TOLERANCE) && (std::fabs(p.d - d) < Constants<T>::TOLERANCE); 
}

template<typename T>
inline bool Plane<T>::isApproxEqual(const Plane<T>& p, T tolerance) const
{
	return (std::fabs(p.a - a) < tolerance) && (std::fabs(p.b - b) < tolerance) && 
		(std::fabs(p.c - c) < tolerance) && (std::fabs(p.d - d) < tolerance); 
}

template<typename T>
inline Plane<T>& Plane<T>::translate(const Vector3<T>& offset)
{
	Vector3 p(-getConstant()*getNormal());
	setConstant(-dot(getNormal(), p + offset));
	return *this;
}

template<typename T>
inline Plane<T>& Plane<T>::transform(const Matrix3<T>& matrix, bool orthogonal)
{
	if (orthogonal)
	{
		Vector3 p(-getConstant()*getNormal());
		set(getNormal()*matrix, -dot(getNormal(), p*matrix);
	}
	else
	{
		Matrix3 normalMatrix;
		normalMatrix.inverseTransposeOf(matrix);
		Vector3 p(-getConstant()*getNormal());
		set(::normalize(getNormal()*normalMatrix), -dot(getNormal(), p*matrix);
	}
	return *this;
}

template<typename T>
inline Plane<T>& Plane<T>::transform(const AffineTransform<T>& transformation, bool orthogonal)
{
	if (orthogonal)
	{
		Vector3 p(-getConstant()*getNormal());
		set(getNormal()*transformation.getBasis(), -dot(getNormal(), ::transform(p, transformation));
	}
	else
	{
		Matrix3 normalMatrix;
		normalMatrix.inverseTransposeOf(transformation.getBasis());
		Vector3 p(-getConstant()*getNormal());
		set(::normalize(getNormal()*normalMatrix), -dot(getNormal(), ::transform(p, transformation));
	}
	return *this;
}

template<typename T>
inline Plane<T>& Plane<T>::normalize()
{
//#if MATHEMATICS_FAST_NORMALIZE
//	if costexpr(std::is_same_v<T, float>)
//	{
//		T m = rcpSqrtApprox(getNormal().getMagnitudeSquared());
//		if (m <= std::numeric_limits<T>::max())
//			*this *= m;
//	}
//	else
//#endif
	{
		T m = getNormal().getMagnitude();
		if (m > T(0))
		{
			m = T(1)/m;
			a *= m;
			b *= m;
			c *= m;
			d *= m;
		}
	}
	return *this;
}

//template<typename T>
//inline Vector3<T> Plane<T>::reflect(const Vector3<T>& point) const
//{ 
//	return (getNormal()*(T(-2)*(dot(getNormal(), point) + d)) + point); 
//}

template<typename T>
template<typename U>
inline T Plane<T>::getDistance(const Vector3<T>& point) const
{
	static_assert(std::is_same_v<U, Normalized> || std::is_same_v<U, Unnormalized>);
	if costexpr(std::is_same_v<U, Normalized>)
		return std::fabs(dot(getNormal(), point) + d);
	else
		return std::fabs((dot(getNormal(), point) + d)/getNormal().getMagnitude());
}

template<typename T>
template<typename U>
inline T Plane<T>::getSignedDistance(const Vector3<T>& point) const
{
	static_assert(std::is_same_v<U, Normalized> || std::is_same_v<U, Unnormalized>);
	if costexpr(std::is_same_v<U, Normalized>)
		return dot(getNormal(), point) + d;
	else
		return (dot(getNormal(), point) + d)/getNormal().getMagnitude();
}

template<typename T>
inline bool Plane<T>::contains(const Vector3<T>& point) const
{ 
	return (std::fabs(dot(getNormal(), point) + d) < Constants<T>::TOLERANCE); 
}

template<typename T>
inline bool Plane<T>::testIntersection(const Line3<T>& line) const
{
	return (std::fabs(dot(getNormal(), line.direction)) >= Constants<T>::TOLERANCE);
}

template<typename T>
inline std::optional<T> Plane<T>::findIntersection(const Line3<T>& line) const
{
	T nd = dot(getNormal(), line.direction);
	if (std::fabs(nd) < Constants<T>::TOLERANCE)
		return {};
	return std::optional<T>((-d - dot(getNormal(), line.origin))/nd);
}

template<typename T>
inline std::optional<T> Plane<T>::findIntersection(const Ray3<T>& ray) const
{
	std::optional<T> result = findIntersection(ray.asLine());
	return (result.has_value() && (result.value() >= T(0))) ? result : std::optional<T>();
}

template<typename T>
inline std::optional<T> Plane<T>::findIntersection(const Segment3<T>& segment) const
{
	std::optional<T> result = findIntersection(Line3(segment.start, segment.end - segment.start));
	return (result.has_value() && (result.value() >= T(0)) && (result.value() <= T(1))) ? result : std::optional<T>();
}

template<typename T>
template<typename U> 
inline std::optional<U> Plane<T>::findIntersection(const Line3<T>& line) const
{
	static_assert(std::is_same_v<U, T> || std::is_same_v<U, Vector3<T>>);
	std::optional<T> result = findIntersection(line);
	if constexpr (std::is_same_v<U, T>)
		return result;
	else //if constexpr (std::is_same_v<U, Vector3<T>>)
		return result.has_value() ? std::optional<U>(line.evaluate(result.value())) : std::optional<U>();
}

template<typename T>
template<typename U> 
inline std::optional<U> Plane<T>::findIntersection(const Ray3<T>& ray) const
{
	static_assert(std::is_same_v<U, T> || std::is_same_v<U, Vector3<T>>);
	std::optional<T> result = findIntersection(ray);
	if constexpr (std::is_same_v<U, T>)
		return result;
	else //if constexpr (std::is_same_v<U, Vector3<T>>)
		return result.has_value() ? std::optional<U>(ray.evaluate(result.value())) : std::optional<U>();
}

template<typename T>
template<typename U> 
inline std::optional<U> Plane<T>::findIntersection(const Segment3<T>& segment) const
{
	static_assert(std::is_same_v<U, T> || std::is_same_v<U, Vector3<T>>);
	std::optional<T> result = findIntersection(segment);
	if constexpr (std::is_same_v<U, T>)
		return result;
	else //if constexpr (std::is_same_v<U, Vector3<T>>)
		return result.has_value() ? std::optional<U>(segment.evaluate(result.value())) : std::optional<U>();
}

#if SIMD_HAS_FLOAT4

inline Plane<float>::Plane<float>(const Vector3<float>& p0, const Vector3<float>& p1, const Vector3<float>& p2)
{
	Vector3<float> normal(cross(p1 - p0, p2 - p0));
	normal.normalize();
	set(normal, -dot(normal, p0));
}

template<typename C, typename T>
inline std::basic_istream<C, T>& operator>>(std::basic_istream<C, T>& s, Plane<float>& p)
{ 
	float a, b, c, d;
	s >> a >> std::ws >> b >> std::ws >> c >> std::ws >> d;
	p.set(a, b, c, d);
	return s; 
}

template<typename A>
inline void Plane<float>::load(A& ar)
{
	float t0, t1, t2, t3;
	ar(t0, t1, t2, t3);
	set(t0, t1, t2, t3);
}

template<std::size_t I>
inline float& Plane<float>::get()
{
	if constexpr (I == 0)
		return x;
	else if constexpr (I == 1)
		return y;
	else if constexpr (I == 2)
		return z;
	else if constexpr (I == 3)
		return w;
	static_assert(false);
}

template<std::size_t I>
inline const float& Plane<float>::get() const
{
	if constexpr (I == 0)
		return x;
	else if constexpr (I == 1)
		return y;
	else if constexpr (I == 2)
		return z;
	else if constexpr (I == 3)
		return w;
	static_assert(false);
}

inline Plane<float>& Plane<float>::translate(const Vector3<float>& offset)
{
	Vector3 p(-getConstant()*getNormal());
	setConstant(-dot(getNormal(), p + offset));
	return *this;
}

inline Plane<float>& Plane<float>::transform(const Matrix3<float>& matrix, bool orthogonal)
{
	if (orthogonal)
	{
		Vector3 p(-getConstant()*getNormal());
		set(getNormal()*matrix, -dot(getNormal(), p*matrix);
	}
	else
	{
		Matrix3 normalMatrix;
		normalMatrix.inverseTransposeOf(matrix);
		Vector3 p(-getConstant()*getNormal());
		set(::normalize(getNormal()*normalMatrix), -dot(getNormal(), p*matrix);
	}
	return *this;
}

inline Plane<float>& Plane<float>::transform(const AffineTransform<float>& transformation, bool orthogonal)
{
	if (orthogonal)
	{
		Vector3 p(-getConstant()*getNormal());
		set(getNormal()*transformation.getBasis(), -dot(getNormal(), ::transform(p, transformation));
	}
	else
	{
		Matrix3 normalMatrix;
		normalMatrix.inverseTransposeOf(transformation.getBasis());
		Vector3 p(-getConstant()*getNormal());
		set(::normalize(getNormal()*normalMatrix), -dot(getNormal(), ::transform(p, transformation));
	}
	return *this;
}

inline Plane<float>& Plane<float>::normalize()
{
#if MATHEMATICS_FAST_NORMALIZE
	float m = simd::toFloat(simd::rcpSqrtApprox1(simd::dot3(abcd, abcd)));
	if (m <= std::numeric_limits<float>::max()) 
		abcd = simd::mul4(abcd, simd::set4(m));
#else
	float m = getNormal().getMagnitude();
	if (m > 0.f) 
		abcd = simd::div4(abcd, simd::set4(m));
#endif
	return *this;
}

//inline Vector3<float> Plane<float>::reflect(const Vector3<float>& point) const
//{ 
//	return (getNormal()*(-2.f*(dot(getNormal(), point) + d)) + point); 
//}

template<typename U>
inline float Plane<float>::getDistance(const Vector3<float>& point) const
{
	static_assert(std::is_same_v<U, Normalized> || std::is_same_v<U, Unnormalized>);
	if costexpr(std::is_same_v<U, Normalized>)
		return std::fabs(dot(getNormal(), point) + d);
	else
		return std::fabs((dot(getNormal(), point) + d)/getNormal().getMagnitude());
}

template<typename U>
inline float Plane<float>::getSignedDistance(const Vector3<float>& point) const
{
	static_assert(std::is_same_v<U, Normalized> || std::is_same_v<U, Unnormalized>);
	if costexpr(std::is_same_v<U, Normalized>)
		return dot(getNormal(), point) + d;
	else
		return (dot(getNormal(), point) + d)/getNormal().getMagnitude();
}

inline bool Plane<float>::contains(const Vector3<float>& point) const
{
	return (std::fabs(dot(getNormal(), point) + d) < Constants<float>::TOLERANCE);
}

inline bool Plane<float>::testIntersection(const Line3<float>& line) const
{
	return (std::fabs(dot(getNormal(), line.direction)) >= Constants<float>::TOLERANCE);
}

inline std::optional<float> Plane<float>::findIntersection(const Line3<float>& line) const
{
	Vector3<float> n = getNormal();
	float nd = dot(n, line.direction);
	if (std::fabs(nd) < Constants<float>::TOLERANCE)
		return {};
	return std::optional<float>((-d - dot(n, line.origin))/nd);
}

inline std::optional<float> Plane<float>::findIntersection(const Ray3<float>& ray) const
{
	std::optional<float> result = findIntersection(ray.asLine());
	return (result.has_value() && (result.value() >= 0.f)) ? result : std::optional<float>();
}

inline std::optional<float> Plane<float>::findIntersection(const Segment3<float>& segment) const
{
	std::optional<float> result = findIntersection(Line3(segment.start, segment.end - segment.start));
	return (result.has_value() && (result.value() >= 0.f) && (result.value() <= 1.f)) ? result : std::optional<float>();
}

template<typename U>
inline std::optional<U> Plane<float>::findIntersection(const Line3<float>& line) const
{
	static_assert(std::is_same_v<U, float> || std::is_same_v<U, Vector3<float>>);
	std::optional<float> result = findIntersection(line);
	if constexpr (std::is_same_v<U, float>)
		return result;
	else //if constexpr (std::is_same_v<U, Vector3<float>>)
		return result.has_value() ? std::optional<U>(line.evaluate(result.value())) : std::optional<U>();
}

template<typename U>
inline std::optional<U> Plane<float>::findIntersection(const Ray3<float>& ray) const
{
	static_assert(std::is_same_v<U, float> || std::is_same_v<U, Vector3<float>>);
	std::optional<float> result = findIntersection(ray);
	if constexpr (std::is_same_v<U, float>)
		return result;
	else //if constexpr (std::is_same_v<U, Vector3<float>>)
		return result.has_value() ? std::optional<U>(ray.evaluate(result.value())) : std::optional<U>();
}

template<typename U>
inline std::optional<U> Plane<float>::findIntersection(const Segment3<float>& segment) const
{
	static_assert(std::is_same_v<U, float> || std::is_same_v<U, Vector3<float>>);
	std::optional<float> result = findIntersection(segment);
	if constexpr (std::is_same_v<U, float>)
		return result;
	else //if constexpr (std::is_same_v<U, Vector3<float>>)
		return result.has_value() ? std::optional<U>(segment.evaluate(result.value())) : std::optional<U>();
}

#endif /* SIMD_HAS_FLOAT4 */

template<std::size_t I, typename T>
	requires std::floating_point<T>
inline T& get(Plane<T>& p) noexcept
{
	if constexpr (I == 0)
		return p.a;
	else if constexpr (I == 1)
		return p.b;
	else if constexpr (I == 2)
		return p.c;
	else if constexpr (I == 3)
		return p.d;
	static_assert(false);
}

template<std::size_t I, typename T>
	requires std::floating_point<T>
inline const T& get(const Plane<T>& p) noexcept
{
	if constexpr (I == 0)
		return p.a;
	else if constexpr (I == 1)
		return p.b;
	else if constexpr (I == 2)
		return p.c;
	else if constexpr (I == 3)
		return p.d;
	static_assert(false);
}

template<std::size_t I, typename T>
	requires std::floating_point<T>
inline T&& get(Plane<T>&& p) noexcept
{
	if constexpr (I == 0)
		return p.a;
	else if constexpr (I == 1)
		return p.b;
	else if constexpr (I == 2)
		return p.c;
	else if constexpr (I == 3)
		return p.d;
	static_assert(false);
}

template<std::size_t I, typename T>
	requires std::floating_point<T>
inline const T&& get(const Plane<T>&& p) noexcept
{
	if constexpr (I == 0)
		return p.a;
	else if constexpr (I == 1)
		return p.b;
	else if constexpr (I == 2)
		return p.c;
	else if constexpr (I == 3)
		return p.d;
	static_assert(false);
}

template<typename T>
	requires std::floating_point<T>
inline Plane<T> normalize(const Plane<T>& p) noexcept
{
	Plane<T> u(p);
	u.normalize();
	return u;
}

template<typename T>
#if SIMD_HAS_FLOAT4
	requires (std::floating_point<T> && !std::same_as<T, float>)
#else
	requires std::floating_point<T>
#endif
inline Plane<T> normalize(Plane<T>&& p) noexcept
{
	p.normalize();
	return p;
}

} // namespace templates

#if MATHEMATICS_DOUBLE
using Plane = templates::Plane<double>;
using PlaneArg = templates::Plane<double>::ConstArg;
using PlaneResult = templates::Plane<double>::ConstResult;
#else
using Plane = templates::Plane<float>;
using PlaneArg = templates::Plane<float>::ConstArg;
using PlaneResult = templates::Plane<float>::ConstResult;
#endif

} // namespace core::mathematics

namespace std {

template<size_t I, typename T>
struct tuple_element;

template<size_t I, typename T>
struct tuple_element<I, ::core::mathematics::templates::Plane<T>>
{
	using type = T;
};

template<typename T>
struct tuple_size;

template<typename T>
struct tuple_size<::core::mathematics::templates::Plane<T>> : integral_constant<size_t, 4> 
{
};

template<typename T>
struct hash;

template<typename T>
struct hash<::core::mathematics::templates::Plane<T>>
{
	std::size_t operator()(const ::core::mathematics::templates::Plane<T>& p) const noexcept
	{
		std::hash<T> hasher;
		std::size_t seed = hasher(p.a) + 0x9e3779b9;
		seed ^= hasher(p.b) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
		seed ^= hasher(p.c) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
		seed ^= hasher(p.d) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
		return seed;
	}
};

} // namespace std
