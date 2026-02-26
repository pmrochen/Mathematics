/*
 *	Name: HalfSpace
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
#include <algorithm>
#include <cstddef>
#include <cmath>
#include <Simd/Intrinsics.hpp>
#include "Constants.hpp"
#include "Vector3.hpp"
#include "Matrix3.hpp"
#include "AffineTransform.hpp"
#include "Segment3.hpp"
#include "Triangle3.hpp"

namespace core::mathematics {

struct Normalized
{
};

struct Unnormalized
{
};

namespace templates {

template<typename T>
	requires std::floating_point<T>
struct Plane;

template<typename T>
	requires std::floating_point<T>
struct HalfSpace
{
	using Real = T;
	using ConstArg = const HalfSpace&;
	using ConstResult = const HalfSpace&;

	constexpr HalfSpace() noexcept : a(), b(), c(), d() {}
	explicit HalfSpace(Uninitialized) noexcept {}
	constexpr HalfSpace(T a, T b, T c, T d) noexcept : a(a), b(b), c(c), d(d) {}
	constexpr HalfSpace(const Vector3<T>& normal, T constant) noexcept : x(normal.x), y(normal.y), z(normal.z), w(constant) {}
	HalfSpace(const Vector3<T>& normal, const Vector3<T>& point) noexcept : a(normal.x), b(normal.y), c(normal.z), d(-dot(normal, point)) {}
	HalfSpace(const Vector3<T>& p0, const Vector3<T>& p1, const Vector3<T>& p2) noexcept;
	HalfSpace(const Plane<T>& p) noexcept;
	explicit HalfSpace(const std::tuple<T, T, T, T>& t) noexcept : a(std::get<0>(t)), b(std::get<1>(t)), c(std::get<2>(t)), d(std::get<3>(t)) {}
	template<typename U> explicit HalfSpace(const std::tuple<U, U, U, U>& t) noexcept : a(T(std::get<0>(t))), b(T(std::get<1>(t))), c(T(std::get<2>(t))), d(T(std::get<3>(t))) {}
	//explicit HalfSpace(const T* h) noexcept : a(h[0]), b(h[1]), c(h[2]), d(h[3]) {}

	//explicit operator std::tuple<T, T, T, T>() noexcept { return std::tuple<T, T, T, T>(a, b, c, d); }
	//template<typename U> explicit operator std::tuple<U, U, U, U>() noexcept { return std::tuple<U, U, U, U>(U(a), U(b), U(c), U(d)); }
	//explicit operator T*() noexcept { return &a; }
	//explicit operator const T*() const noexcept { return &a; }
	//T& operator[](int i) noexcept { return (&a)[i]; }
	//const T& operator[](int i) const noexcept { return (&a)[i]; }

	bool operator==(const HalfSpace& h) const noexcept { return (a == h.a) && (b == h.b) && (c == h.c) && (d == h.d); }
	bool operator!=(const HalfSpace& h) const noexcept { return !(*this == h); }

	template<typename A> void load(A& ar) { ar(a, b, c, d); }
	template<typename A> void save(A& ar) const { ar(a, b, c, d); }

	template<std::size_t I> T& get() noexcept;
	template<std::size_t I> const T& get() const noexcept;

	const Plane<T>& asPlane() const noexcept;

	// Properties
	bool isEmpty() const noexcept { return (a == T()) && (b == T()) && (c == T()) && (d == T()); }
	bool isApproxEqual(const HalfSpace& h) const noexcept;
	bool isApproxEqual(const HalfSpace& h, T tolerance) const noexcept;
	bool isFinite() const noexcept { return std::isfinite(a) && std::isfinite(b) && std::isfinite(c) && std::isfinite(d); }
	const Vector3<T>& getNormal() const noexcept { return reinterpret_cast<const Vector3<T>&>(*this); }
	void setNormal(const Vector3<T>& normal) noexcept { a = normal.x; b = normal.y; c = normal.z; }
	T getConstant() const noexcept { return d; }
	void setConstant(T constant) noexcept { d = constant; }
	HalfSpace& makeEmpty() noexcept { a = T(); b = T(); c = T(); d = T(); return *this; }
	HalfSpace& set(const Vector3<T>& normal, T constant) noexcept { a = normal.x; b = normal.y; c = normal.z; d = constant; return *this; }
	HalfSpace& set(T a, T b, T c, T d) noexcept { this->a = a; this->b = b; this->c = c; this->d = d; return *this; }

	// Transformation
	HalfSpace& translate(const Vector3<T>& offset) noexcept;
	HalfSpace& transform(const Matrix3<T>& matrix, bool orthogonal = false) noexcept;
	HalfSpace& transform(const AffineTransform<T>& transformation, bool orthogonal = false) noexcept;
	HalfSpace& flip() noexcept { a = -a; b = -b; c = -c; d = -d; return *this; }
	HalfSpace& normalize() noexcept;

	// Distances
	T getDistance(const Vector3<T>& point) const { return std::max(dot(getNormal(), point) + d, 0.f); }		// normalized half-space
	template<typename U> T getDistance(const Vector3<T>& point) const;
	T getSignedDistance(const Vector3<T>& point) const noexcept { return (dot(getNormal(), point) + d); }	// normalized half-space
	template<typename U> T getSignedDistance(const Vector3<T>& point) const noexcept;

	// Intersection
	bool contains(const Vector3<T>& point) const noexcept { return (dot(getNormal(), point) <= -d); }
	bool testIntersection(const Segment3<T>& segment) const noexcept;
	bool testIntersection(const Triangle3<T>& triangle) const noexcept;

	//static const HalfSpace& getEmpty() noexcept { return EMPTY; }

	static const HalfSpace EMPTY;

	T a, b, c, d;
};

template<typename T> const HalfSpace<T> HalfSpace<T>::EMPTY{};

#if SIMD_HAS_FLOAT4

template<>
struct Plane<float>;

template<>
struct HalfSpace<float>
{
	using Real = float;
	using ConstArg = const HalfSpace;
	using ConstResult = const HalfSpace;

	/*constexpr*/ HalfSpace() noexcept : abcd(simd::zero<simd::float4>()) {}
	explicit HalfSpace(Uninitialized) noexcept {}
	/*constexpr*/ HalfSpace(float a, float b, float c, float d) noexcept : abcd(simd::set4(a, b, c, d)) {}
	/*constexpr*/ HalfSpace(const Vector3<float>& normal, float constant) noexcept : abcd(simd::insert<simd::W>(constant, normal)) {}
	HalfSpace(const Vector3<float>& normal, const Vector3<float>& point) noexcept : abcd(simd::insert<simd::W>(-dot(normal, point), normal)) {}
	HalfSpace(const Vector3<float>& p0, const Vector3<float>& p1, const Vector3<float>& p2) noexcept;
	HalfSpace(const Plane<float>& p) noexcept;
	explicit HalfSpace(const std::tuple<float, float, float, float>& t) noexcept : abcd(simd::set4(std::get<0>(t), std::get<1>(t), std::get<2>(t), std::get<3>(t))) {}
	template<typename U> explicit HalfSpace(const std::tuple<U, U, U, U>& t) noexcept : abcd(simd::set4((float)std::get<0>(t), (float)std::get<1>(t), (float)std::get<2>(t), (float)std::get<3>(t))) {}
	//explicit HalfSpace(const float* h) noexcept : abcd(simd::load4(h)) {}
	explicit HalfSpace(simd::float4 h) noexcept : abcd(h) {}
	HalfSpace(const HalfSpace& h) noexcept : abcd(h.abcd) {}
	HalfSpace& operator=(const HalfSpace& h) noexcept { abcd = h.abcd; return *this; }

	operator simd::float4() const noexcept { return abcd; }
	//explicit operator std::tuple<float, float, float, float>() noexcept { return std::tuple<float, float, float, float>(a, b, c, d); }
	//template<typename U> explicit operator std::tuple<U, U, U, U>() noexcept { return std::tuple<U, U, U, U>(U(a), U(b), U(c), U(d)); }
	//explicit operator float* () noexcept { return &a; }
	//explicit operator const float* () const noexcept { return &a; }
	//float& operator[](int i) noexcept { return (&a)[i]; }
	//const float& operator[](int i) const noexcept { return (&a)[i]; }

	bool operator==(const HalfSpace& h) const noexcept { return simd::all4(simd::equal(abcd, h)); }
	bool operator!=(const HalfSpace& h) const noexcept { return !(*this == h); }

	template<typename A> void load(A& ar);
	template<typename A> void save(A& ar) const { ar(a, b, c, d); }

	template<std::size_t I> float& get() noexcept;
	template<std::size_t I> const float& get() const noexcept;
	template<typename U> U& get() noexcept;				// intentionally undefined
	template<typename U> const U& get() const noexcept;	// intentionally undefined
	template<> simd::float4& get() noexcept { return abcd; }
	template<> const simd::float4& get() const noexcept { return abcd; }

	const Plane<float> asPlane() const noexcept;

	// Properties
	bool isEmpty() const noexcept { return simd::all4(simd::equal(abcd, simd::zero<simd::float4>())); }
	bool isApproxEqual(const HalfSpace& h) const noexcept { simd::all4(simd::lessThan(simd::abs4(simd::sub4(abcd, h)), simd::set4(Constants<float>::TOLERANCE))); }
	bool isApproxEqual(const HalfSpace& h, float tolerance) const noexcept { simd::all4(simd::lessThan(simd::abs4(simd::sub4(abcd, h)), simd::set4(tolerance))); }
	bool isFinite() const noexcept { return simd::all4(simd::isFinite(abcd)); }
#if MATHEMATICS_SIMD_EXPAND_LAST
	const Vector3<float> getNormal() const noexcept { return Vector3<float>(simd::xyzz(abcd)); }
#else
	const Vector3<float> getNormal() const noexcept { return Vector3<float>(simd::cutoff3(abcd)); }
#endif
	void setNormal(const Vector3<float>& normal) noexcept { abcd = simd::insert3(normal, abcd); }
	float getConstant() const noexcept { return simd::extract<simd::W>(abcd); }
	void setConstant(float constant) noexcept { abcd = simd::insert<simd::W>(constant, abcd); }
	HalfSpace& makeEmpty() noexcept { abcd = simd::zero<simd::float4>(); return *this; }
	HalfSpace& set(const Vector3<float>& normal, float constant) noexcept { abcd = simd::insert<simd::W>(constant, normal); return *this; }
	HalfSpace& set(float a, float b, float c, float d) noexcept { abcd = simd::set4(a, b, c, d); return *this; }

	// Transformation
	HalfSpace& translate(const Vector3<float>& offset) noexcept;
	HalfSpace& transform(const Matrix3<float>& matrix, bool orthogonal = false) noexcept;
	HalfSpace& transform(const AffineTransform<float>& transformation, bool orthogonal = false) noexcept;
	HalfSpace& flip() noexcept { abcd = simd::neg4(abcd); return *this; }
	HalfSpace& normalize() noexcept;

	// Distances
	float getDistance(const Vector3<float>& point) const { return std::max(dot(getNormal(), point) + d, 0.f); }		// normalized half-space
	template<typename U> float getDistance(const Vector3<float>& point) const;
	float getSignedDistance(const Vector3<float>& point) const noexcept { return (dot(getNormal(), point) + d); }	// normalized half-space
	template<typename U> float getSignedDistance(const Vector3<float>& point) const noexcept;

	// Intersection
	bool contains(const Vector3<float>& point) const noexcept { return (dot(getNormal(), point) <= -d); }
	bool testIntersection(const Segment3<float>& segment) const noexcept;
	bool testIntersection(const Triangle3<float>& triangle) const noexcept;

	//static const HalfSpace& getEmpty() noexcept { return EMPTY; }

	static const HalfSpace EMPTY;

	union
	{
		simd::float4 abcd;
		struct { float a, b, c, d; };
	};
};

const HalfSpace<float> HalfSpace<float>::EMPTY{};

#endif /* SIMD_HAS_FLOAT4 */

template<typename T>
inline HalfSpace<T>::HalfSpace<T>(const Vector3<T>& p0, const Vector3<T>& p1, const Vector3<T>& p2)
{
	Vector3<T> normal(cross(p1 - p0, p2 - p0));
	normal.normalize();
	set(normal, -dot(normal, p0));
}

template<typename C, typename T, typename U>
	requires std::floating_point<U>
inline std::basic_istream<C, T>& operator>>(std::basic_istream<C, T>& s, HalfSpace<U>& h)
{ 
	return s >> h.a >> std::ws >> h.b >> std::ws >> h.c >> std::ws >> h.d; 
}

template<typename C, typename T, typename U>
	requires std::floating_point<U>
inline std::basic_ostream<C, T>& operator<<(std::basic_ostream<C, T>& s, const HalfSpace<U>& h)
{ 
	constexpr C WS(0x20);
	return s << h.a << WS << h.b << WS << h.c << WS << h.d;
}

template<typename T>
template<std::size_t I>
inline T& HalfSpace<T>::get()
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
inline const T& HalfSpace<T>::get() const
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
inline bool HalfSpace<T>::isApproxEqual(const HalfSpace<T>& h) const
{ 
	return (std::fabs(h.a - a) < Constants<T>::TOLERANCE) && (std::fabs(h.b - b) < Constants<T>::TOLERANCE) && 
		(std::fabs(h.c - c) < Constants<T>::TOLERANCE) && (std::fabs(h.d - d) < Constants<T>::TOLERANCE); 
}

template<typename T>
inline bool HalfSpace<T>::isApproxEqual(const HalfSpace<T>& h, T tolerance) const
{
	return (std::fabs(h.a - a) < tolerance) && (std::fabs(h.b - b) < tolerance) && 
		(std::fabs(h.c - c) < tolerance) && (std::fabs(h.d - d) < tolerance); 
}

template<typename T>
inline HalfSpace<T>& HalfSpace<T>::translate(const Vector3<T>& offset)
{
	Vector3 p(-getConstant()*getNormal());
	setConstant(-dot(getNormal(), p + offset));
	return *this;
}

template<typename T>
inline HalfSpace<T>& HalfSpace<T>::transform(const Matrix3<T>& matrix, bool orthogonal)
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
inline HalfSpace<T>& HalfSpace<T>::transform(const AffineTransform<T>& transformation, bool orthogonal)
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
inline HalfSpace<T>& HalfSpace<T>::normalize()
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

template<typename T>
template<typename U>
inline T HalfSpace<T>::getDistance<(const Vector3<T>& point) const
{
	static_assert(std::is_same_v<U, Normalized> || std::is_same_v<U, Unnormalized>);
	if costexpr(std::is_same_v<U, Normalized>)
		return std::max(dot(getNormal(), point) + d, T(0));
	else
		return std::max((dot(getNormal(), point) + d)/getNormal().getMagnitude(), T(0));
}

template<typename T>
template<typename U> 
inline T HalfSpace<T>::getSignedDistance(const Vector3<T>& point) const
{
	static_assert(std::is_same_v<U, Normalized> || std::is_same_v<U, Unnormalized>);
	if costexpr(std::is_same_v<U, Normalized>)
		return dot(getNormal(), point) + d;
	else
		return (dot(getNormal(), point) + d)/getNormal().getMagnitude(), T(0);
}

template<typename T>
inline bool HalfSpace<T>::testIntersection(const Segment3<T>& segment) const
{
	return contains(segment.start) || contains(segment.end);
}

template<typename T>
inline bool HalfSpace<T>::testIntersection(const Triangle3<T>& triangle) const
{
	return contains(triangle.vertices[0]) || contains(triangle.vertices[1]) || contains(triangle.vertices[2]);
}

#if SIMD_HAS_FLOAT4

inline HalfSpace<float>::HalfSpace<float>(const Vector3<float>& p0, const Vector3<float>& p1, const Vector3<float>& p2)
{
	Vector3<float> normal(cross(p1 - p0, p2 - p0));
	normal.normalize();
	set(normal, -dot(normal, p0));
}

template<typename C, typename T>
inline std::basic_istream<C, T>& operator>>(std::basic_istream<C, T>& s, HalfSpace<float>& h)
{ 
	float a, b, c, d; 
	s >> a >> std::ws >> b >> std::ws >> c >> std::ws >> d;
	h.set(a, b, c, d);
	return s; 
}

template<typename A>
inline void HalfSpace<float>::load(A& ar)
{
	float t0, t1, t2, t3;
	ar(t0, t1, t2, t3);
	set(t0, t1, t2, t3);
}

template<std::size_t I>
inline float& HalfSpace<float>::get()
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
inline const float& HalfSpace<float>::get() const
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

inline HalfSpace<float>& HalfSpace<float>::translate(const Vector3<float>& offset)
{
	Vector3 p(-getConstant()*getNormal());
	setConstant(-dot(getNormal(), p + offset));
	return *this;
}

inline HalfSpace<float>& HalfSpace<float>::transform(const Matrix3<float>& matrix, bool orthogonal)
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

inline HalfSpace<float>& HalfSpace<float>::transform(const AffineTransform<float>& transformation, bool orthogonal)
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

inline HalfSpace<float>& HalfSpace<float>::normalize()
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

template<typename U>
inline float HalfSpace<float>::getDistance<(const Vector3<float>& point) const
{
	static_assert(std::is_same_v<U, Normalized> || std::is_same_v<U, Unnormalized>);
	if costexpr(std::is_same_v<U, Normalized>)
		return std::max(dot(getNormal(), point) + d, T(0));
	else
		return std::max((dot(getNormal(), point) + d)/getNormal().getMagnitude(), T(0));
}

template<typename U>
inline float HalfSpace<float>::getSignedDistance(const Vector3<float>& point) const
{
	static_assert(std::is_same_v<U, Normalized> || std::is_same_v<U, Unnormalized>);
	if costexpr(std::is_same_v<U, Normalized>)
		return dot(getNormal(), point) + d;
	else
		return (dot(getNormal(), point) + d)/getNormal().getMagnitude(), T(0);
}

inline bool HalfSpace<float>::testIntersection(const Segment3<float>& segment) const
{
	return contains(segment.start) || contains(segment.end);
}

inline bool HalfSpace<float>::testIntersection(const Triangle3<float>& triangle) const
{
	return contains(triangle.vertices[0]) || contains(triangle.vertices[1]) || contains(triangle.vertices[2]);
}

#endif /* SIMD_HAS_FLOAT4 */

template<std::size_t I, typename T>
	requires std::floating_point<T>
inline T& get(HalfSpace<T>& h) noexcept
{
	if constexpr (I == 0)
		return h.a;
	else if constexpr (I == 1)
		return h.b;
	else if constexpr (I == 2)
		return h.c;
	else if constexpr (I == 3)
		return h.d;
	static_assert(false);
}

template<std::size_t I, typename T>
	requires std::floating_point<T>
inline const T& get(const HalfSpace<T>& h) noexcept
{
	if constexpr (I == 0)
		return h.a;
	else if constexpr (I == 1)
		return h.b;
	else if constexpr (I == 2)
		return h.c;
	else if constexpr (I == 3)
		return h.d;
	static_assert(false);
}

template<std::size_t I, typename T>
	requires std::floating_point<T>
inline T&& get(HalfSpace<T>&& h) noexcept
{
	if constexpr (I == 0)
		return h.a;
	else if constexpr (I == 1)
		return h.b;
	else if constexpr (I == 2)
		return h.c;
	else if constexpr (I == 3)
		return h.d;
	static_assert(false);
}

template<std::size_t I, typename T>
	requires std::floating_point<T>
inline const T&& get(HalfSpace<T>&& h) noexcept
{
	if constexpr (I == 0)
		return h.a;
	else if constexpr (I == 1)
		return h.b;
	else if constexpr (I == 2)
		return h.c;
	else if constexpr (I == 3)
		return h.d;
	static_assert(false);
}

template<typename T>
	requires std::floating_point<T>
inline HalfSpace<T> normalize(const HalfSpace<T>& h) noexcept
{
	HalfSpace<T> u(h);
	u.normalize();
	return u;
}

template<typename T>
#if SIMD_HAS_FLOAT4
	requires (std::floating_point<T> && !std::same_as<T, float>)
#else
	requires std::floating_point<T>
#endif
inline HalfSpace<T> normalize(HalfSpace<T>&& h) noexcept
{
	h.normalize();
	return h;
}

} // namespace templates

#if MATHEMATICS_DOUBLE
using HalfSpace = templates::HalfSpace<double>;
using HalfSpaceArg = templates::HalfSpace<double>::ConstArg;
using HalfSpaceResult = templates::HalfSpace<double>::ConstResult;
#else
using HalfSpace = templates::HalfSpace<float>;
using HalfSpaceArg = templates::HalfSpace<float>::ConstArg;
using HalfSpaceResult = templates::HalfSpace<float>::ConstResult;
#endif

} // namespace core::mathematics

namespace std {

template<size_t I, typename T>
struct tuple_element;

template<size_t I, typename T>
struct tuple_element<I, ::core::mathematics::templates::HalfSpace<T>>
{
	using type = T;
};

template<typename T>
struct tuple_size;

template<typename T>
struct tuple_size<::core::mathematics::templates::HalfSpace<T>> : integral_constant<size_t, 4> 
{
};

template<typename T>
struct hash;

template<typename T>
struct hash<::core::mathematics::templates::HalfSpace<T>>
{
	std::size_t operator()(const ::core::mathematics::templates::HalfSpace<T>& h) const noexcept
	{
		std::hash<T> hasher;
		std::size_t seed = hasher(h.a) + 0x9e3779b9;
		seed ^= hasher(h.b) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
		seed ^= hasher(h.c) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
		seed ^= hasher(h.d) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
		return seed;
	}
};

} // namespace std

#include "Plane.hpp"

namespace core::mathematics::templates {

template<typename T>
inline HalfSpace<T>::HalfSpace(const Plane<T>& p) : a(p.a), b(p.b), c(p.c), d(p.d)
{
}

template<typename T>
inline const Plane<T>& HalfSpace<T>::asPlane() const 
{ 
	return reinterpret_cast<const Plane<T>&>(*this); 
}

#if SIMD_HAS_FLOAT4

inline HalfSpace<float>::HalfSpace(const Plane<float>& p) : abcd(p.abcd)
{
}

inline const Plane<float> HalfSpace<float>::asPlane() const
{ 
	return Plane<float>(abcd); 
}

#endif /* SIMD_HAS_FLOAT4 */

} // namespace core::mathematics::templates
