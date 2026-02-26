/*
 *	Name: Vector3
 *	Author: Pawel Mrochen
 */

#pragma once

#include <istream>
#include <ostream>
#include <limits>
#include <type_traits>
#include <concepts>
#include <algorithm>
#include <utility>
#include <tuple>
#include <cstddef>
#include <cmath>
#include <Simd/Intrinsics.hpp>
//#include <Tuples/Tuple3.hpp>
#include "Constants.hpp"
#include "Axis.hpp"
#include "Vector2.hpp"

namespace core::mathematics {
namespace templates {

template<typename T>
	requires std::floating_point<T>
struct Matrix3;

template<typename T>
	requires std::floating_point<T>
struct AffineTransform;

template<typename T>
	requires std::floating_point<T>
struct Quaternion;

template<typename T>
	requires (std::floating_point<T> || std::integral<T>)
struct Vector3;

template<std::floating_point T>
struct Vector3<T>
{
	using Real = T;
	using ComponentType = T;
	using ConstArg = const Vector3&;
	using ConstResult = const Vector3&;

	static constexpr int NUM_COMPONENTS = 3;

	constexpr Vector3() noexcept : x(), y(), z() {}
	explicit Vector3(Uninitialized) noexcept {}
	constexpr explicit Vector3(T scalar) noexcept : x(scalar), y(scalar), z(scalar) {}
	constexpr Vector3(T x, T y, T z) noexcept : x(x), y(y), z(z) {}
	constexpr Vector3(const Vector2<T>& v) noexcept : x(v.x), y(v.y), z() {}
	constexpr Vector3(const Vector2<T>& v, T z) noexcept : x(v.x), y(v.y), z(z) {}
	//explicit Vector3(const tuples::templates::Tuple3<T>& t) noexcept : x(t.x), y(t.y), z(t.z) {}
	//template<typename U> explicit Vector3(const tuples::templates::Tuple3<U>& t) noexcept : x(T(t.x)), y(T(t.y)), z(T(t.z)) {}
	explicit Vector3(const std::tuple<T, T, T>& t) noexcept : x(std::get<0>(t)), y(std::get<1>(t)), z(std::get<2>(t)) {}
	template<typename U> explicit Vector3(const std::tuple<U, U, U>& t) noexcept : x(T(std::get<0>(t))), y(T(std::get<1>(t))), z(T(std::get<2>(t))) {}
	explicit Vector3(const T* v) noexcept : x(v[0]), y(v[1]), z(v[2]) {}
	explicit Vector3(Axis axis) noexcept : x((axis == Axis::X) ? T(1) : T(0)), y((axis == Axis::Y) ? T(1) : T(0)), z((axis == Axis::Z) ? T(1) : T(0)) {}
	template<typename U> explicit Vector3(const Vector3<U>& v) noexcept : x(T(t.x)), y(T(t.y)), z(T(t.z)) {}

	//explicit operator tuples::templates::Tuple3<T>() noexcept { return tuples::templates::Tuple3<T>(x, y, z); }
	//template<typename U> explicit operator tuples::templates::Tuple3<U>() noexcept { return tuples::templates::Tuple3<U>(U(x), U(y), U(z)); }
	//explicit operator std::tuple<T, T, T>() { return std::tuple<T, T, T>(x, y, z); }
	//template<typename U> explicit operator std::tuple<U, U, U>() { return std::tuple<U, U, U>(U(x), U(y), U(z)); }
	explicit operator T*() noexcept { return &x; }
	explicit operator const T*() const noexcept { return &x; }
	T& operator[](int i) noexcept { return (&x)[i]; }
	const T& operator[](int i) const noexcept { return (&x)[i]; }

	Vector3 operator+() const noexcept { return *this; }
	Vector3 operator-() const noexcept { return Vector3(-x, -y, -z); }
	Vector3& operator+=(const Vector3& v) noexcept { x += v.x; y += v.y; z += v.z; return *this; }
	Vector3& operator-=(const Vector3& v) noexcept { x -= v.x; y -= v.y; z -= v.z; return *this; }
	Vector3& operator*=(const Vector3& v) noexcept { x *= v.x; y *= v.y; z *= v.z; return *this; }
	Vector3& operator*=(T f) noexcept { x *= f; y *= f; z *= f; return *this; }
	Vector3& operator*=(const Matrix3<T>& m) noexcept;
	Vector3& operator/=(const Vector3& v) noexcept { x /= v.x; y /= v.y; z /= v.z; return *this; }
	Vector3& operator/=(T f) noexcept { return operator*=(T(1)/f); }
	bool operator==(const Vector3& v) const noexcept { return (x == v.x) && (y == v.y) && (z == v.z); }
	bool operator!=(const Vector3& v) const noexcept { return !(*this == v); }
	
	template<typename A> void load(A& ar) { ar(x, y, z); }
	template<typename A> void save(A& ar) const { ar(x, y, z); }

	template<std::size_t I> T& get() noexcept;
	template<std::size_t I> const T& get() const noexcept;

	const Vector2<T>& xy/*getXY*/() const noexcept { return reinterpret_cast<const Vector2<T>&>(*this); }
	const Vector2<T> xz/*getXZ*/() const noexcept { return Vector2<T>(x, z); }
	const Vector2<T> zy/*getZY*/() const noexcept { return Vector2<T>(z, y); }
	//void setXY(const Vector2<T>& v) noexcept { x = v.x; y = v.y; }
	bool isZero() const noexcept { return (x == T()) && (y == T()) && (z == T()); }
	bool isApproxZero() const noexcept;
	bool isApproxEqual(const Vector3& v) const noexcept;
	bool isApproxEqual(const Vector3& v, T tolerance) const noexcept;
	bool allLessThan(const Vector3& v) const noexcept { return (x < v.x) && (y < v.y) && (z < v.z); }
	bool allLessThanEqual(const Vector3& v) const noexcept { return (x <= v.x) && (y <= v.y) && (z <= v.z); }
	bool allGreaterThan(const Vector3& v) const noexcept { return (x > v.x) && (y > v.y) && (z > v.z); }
	bool allGreaterThanEqual(const Vector3& v) const noexcept { return (x >= v.x) && (y >= v.y) && (z >= v.z); }
	bool anyLessThan(const Vector3& v) const noexcept { return (x < v.x) || (y < v.y) || (z < v.z); }
	bool anyLessThanEqual(const Vector3& v) const noexcept { return (x <= v.x) || (y <= v.y) || (z <= v.z); }
	bool anyGreaterThan(const Vector3& v) const noexcept { return (x > v.x) || (y > v.y) || (z > v.z); }
	bool anyGreaterThanEqual(const Vector3& v) const noexcept { return (x >= v.x) || (y >= v.y) || (z >= v.z); }
	bool isFinite() const noexcept { return std::isfinite(x) && std::isfinite(y) && std::isfinite(z); }
	T getMagnitude() const noexcept { return std::sqrt(x*x + y*y + z*z); }
	T getMagnitudeSquared() const noexcept { return (x*x + y*y + z*z); }
	void setMagnitude(T magnitude) noexcept;
	T getLength() const noexcept { return getMagnitude(); }
	T getLengthSquared() const noexcept { return getMagnitudeSquared(); }
	void setLength(T length) noexcept { setMagnitude(length); }
	Axis getMajorAxis() const noexcept;
	T getMinComponent() const { return std::min(std::min(x, y), z); }
	T getMaxComponent() const { return std::max(std::max(x, y), z); }
	Vector3& setZero/*zero*/() noexcept { x = T(); y = T(); z = T(); return *this; }
	Vector3& set(T x, T y, T z) noexcept { this->x = x; this->y = y; this->z = z; return *this; }
	Vector3& setMinimum(const Vector3& v1, const Vector3& v2);
	Vector3& setMaximum(const Vector3& v1, const Vector3& v2);
	Vector3& negate() noexcept { x = -x; y = -y; z = -z; return *this; }
	Vector3& normalize() noexcept;
	Vector3& rotate(Axis axis, T angle) noexcept;
	Vector3& rotate(const Quaternion<T>& q) noexcept;
	//Vector3& scale(const Vector3& v) noexcept { x *= v.x; y *= v.y; z *= v.z; return *this; }
	Vector3& transform(const Matrix3<T>& m) noexcept;
	Vector3& transform(const AffineTransform<T>& m) noexcept;

	//static const Vector3& getZero() noexcept { return ZERO; }
	//static const Vector3& getUnitX() noexcept { return UNIT_X; }
	//static const Vector3& getUnitY() noexcept { return UNIT_Y; }
	//static const Vector3& getUnitZ() noexcept { return UNIT_Z; }

	static const Vector3 ZERO;
	static const Vector3 UNIT_X;
	static const Vector3 UNIT_Y;
	static const Vector3 UNIT_Z; 
	static const Vector3 ONE;
	static const Vector3 TOLERANCE;
	static const Vector3 INF;
	static const Vector3 MINUS_INF;

	T x, y, z;
};

template<std::floating_point T> const Vector3<T> Vector3<T>::ZERO{};
template<std::floating_point T> const Vector3<T> Vector3<T>::UNIT_X{ T(1), T(0), T(0) };
template<std::floating_point T> const Vector3<T> Vector3<T>::UNIT_Y{ T(0), T(1), T(0) };
template<std::floating_point T> const Vector3<T> Vector3<T>::UNIT_Z{ T(0), T(0), T(1) };
template<std::floating_point T> const Vector3<T> Vector3<T>::ONE{ T(1), T(1), T(1) };
template<std::floating_point T> const Vector3<T> Vector3<T>::TOLERANCE{ Constants<T>::TOLERANCE, Constants<T>::TOLERANCE, Constants<T>::TOLERANCE };
template<std::floating_point T> const Vector3<T> Vector3<T>::INF{ std::numeric_limits<T>::infinity(), std::numeric_limits<T>::infinity(), std::numeric_limits<T>::infinity() };
template<std::floating_point T> const Vector3<T> Vector3<T>::MINUS_INF{ -std::numeric_limits<T>::infinity(), -std::numeric_limits<T>::infinity(), -std::numeric_limits<T>::infinity() };

template<std::integral T>
struct Vector3<T>
{
	using ComponentType = T;
	using ConstArg = const Vector3&;
	using ConstResult = const Vector3&;

	static constexpr int NUM_COMPONENTS = 3;

	constexpr Vector3() noexcept : x(), y(), z() {}
	explicit Vector3(Uninitialized) noexcept {}
	constexpr explicit Vector3(T scalar) noexcept : x(scalar), y(scalar), z(scalar) {}
	constexpr Vector3(T x, T y, T z) noexcept : x(x), y(y), z(z) {}
	constexpr Vector3(const Vector2<T>& v) noexcept : x(v.x), y(v.y), z() {}
	constexpr Vector3(const Vector2<T>& v, T z) noexcept : x(v.x), y(v.y), z(z) {}
	//explicit Vector3(const tuples::templates::Tuple3<T>& t) noexcept : x(t.x), y(t.y), z(t.z) {}
	//template<typename U> explicit Vector3(const tuples::templates::Tuple3<U>& t) noexcept : x(T(t.x)), y(T(t.y)), z(T(t.z)) {}
	explicit Vector3(const std::tuple<T, T, T>& t) noexcept : x(std::get<0>(t)), y(std::get<1>(t)), z(std::get<2>(t)) {}
	template<typename U> explicit Vector3(const std::tuple<U, U, U>& t) noexcept : x(T(std::get<0>(t))), y(T(std::get<1>(t))), z(T(std::get<2>(t))) {}
	explicit Vector3(const T* v) noexcept : x(v[0]), y(v[1]), z(v[2]) {}
	template<typename U> explicit Vector3(const Vector3<U>& v) noexcept : x(T(v.x)), y(T(v.y)), z(T(v.z)) {}

	//explicit operator tuples::templates::Tuple3<T>() noexcept { return tuples::templates::Tuple3<T>(x, y, z); }
	//template<typename U> explicit operator tuples::templates::Tuple3<U>() noexcept { return tuples::templates::Tuple3<U>(U(x), U(y), U(z)); }
	//explicit operator std::tuple<T, T, T>() { return std::tuple<T, T, T>(x, y, z); }
	//template<typename U> explicit operator std::tuple<U, U, U>() { return std::tuple<U, U, U>(U(x), U(y), U(z)); }
	explicit operator T* () noexcept { return &x; }
	explicit operator const T* () const noexcept { return &x; }
	T& operator[](int i) noexcept { return (&x)[i]; }
	const T& operator[](int i) const noexcept { return (&x)[i]; }

	Vector3 operator+() const noexcept { return *this; }
	Vector3 operator-() const noexcept { return Vector3(-x, -y, -z); }
	Vector3& operator+=(const Vector3& v) noexcept { x += v.x; y += v.y; z += v.z; return *this; }
	Vector3& operator-=(const Vector3& v) noexcept { x -= v.x; y -= v.y; z -= v.z; return *this; }
	Vector3& operator*=(const Vector3& v) noexcept { x *= v.x; y *= v.y; z *= v.z; return *this; }
	Vector3& operator*=(T f) noexcept { x *= f; y *= f; z *= f; return *this; }
	Vector3& operator/=(const Vector3& v) noexcept { x /= v.x; y /= v.y; z /= v.z; return *this; }
	Vector3& operator/=(T f) noexcept { x /= f; y /= f; z /= f; return *this; }
	bool operator==(const Vector3& v) const noexcept { return (x == v.x) && (y == v.y) && (z == v.z); }
	bool operator!=(const Vector3& v) const noexcept { return !(*this == v); }

	template<typename A> void serialize(A& ar) { ar(x, y, z); }

	template<std::size_t I> T& get() noexcept;
	template<std::size_t I> const T& get() const noexcept;

	const Vector2<T>& xy/*getXY*/() const noexcept { return reinterpret_cast<const Vector2<T>&>(*this); }
	Vector2<T> xz/*getXZ*/() const noexcept { return Vector2<T>(x, z); }
	Vector2<T> zy/*getZY*/() const noexcept { return Vector2<T>(z, y); }
	//void setXY(const Vector2<T>& v) noexcept { x = v.x; y = v.y; }
	bool isZero() const noexcept { return (x == T()) && (y == T()) && (z == T()); }
	bool allLessThan(const Vector3& v) const noexcept { return (x < v.x) && (y < v.y) && (z < v.z); }
	bool allLessThanEqual(const Vector3& v) const noexcept { return (x <= v.x) && (y <= v.y) && (z <= v.z); }
	bool allGreaterThan(const Vector3& v) const noexcept { return (x > v.x) && (y > v.y) && (z > v.z); }
	bool allGreaterThanEqual(const Vector3& v) const noexcept { return (x >= v.x) && (y >= v.y) && (z >= v.z); }
	bool anyLessThan(const Vector3& v) const noexcept { return (x < v.x) || (y < v.y) || (z < v.z); }
	bool anyLessThanEqual(const Vector3& v) const noexcept { return (x <= v.x) || (y <= v.y) || (z <= v.z); }
	bool anyGreaterThan(const Vector3& v) const noexcept { return (x > v.x) || (y > v.y) || (z > v.z); }
	bool anyGreaterThanEqual(const Vector3& v) const noexcept { return (x >= v.x) || (y >= v.y) || (z >= v.z); }
	T getMinComponent() const { return std::min(std::min(x, y), z); }
	T getMaxComponent() const { return std::max(std::max(x, y), z); }
	Vector3& setZero/*zero*/() noexcept { x = T(); y = T(); z = T(); return *this; }
	Vector3& set(T x, T y, T z) noexcept { this->x = x; this->y = y; this->z = z; return *this; }
	Vector3& setMinimum(const Vector3& v1, const Vector3& v2);
	Vector3& setMaximum(const Vector3& v1, const Vector3& v2);
	Vector3& negate() noexcept { x = -x; y = -y; z = -z; return *this; }

	//static const Vector3& getZero() noexcept { return ZERO; }

	static const Vector3 ZERO;

	T x, y, z;
};

template<std::integral T> const Vector3<T> Vector3<T>::ZERO{};

#if SIMD_HAS_FLOAT4

template<>
struct Matrix3<float>;

template<>
struct AffineTransform<float>;

template<>
struct Quaternion<float>;

template<>
struct Vector3<float>
{
	using Real = float;
	using ComponentType = float;
	using ConstArg = const Vector3;
	using ConstResult = const Vector3;

	static constexpr int NUM_COMPONENTS = 3;

	/*constexpr*/ Vector3() noexcept : xyz(simd::zero<simd::float4>()) {}
	explicit Vector3(Uninitialized) noexcept {}
#if MATHEMATICS_SIMD_EXPAND_LAST
	/*constexpr*/ explicit Vector3(float scalar) noexcept : xyz(simd::set4(scalar)) {}
	/*constexpr*/ Vector3(float x, float y, float z) noexcept : xyz(simd::set4(x, y, z, z)) {}
	/*constexpr*/ Vector3(const Vector2<float>& v) noexcept : xyz(simd::cutoff2(v)) {}
	/*constexpr*/ Vector3(const Vector2<float>& v, float z) noexcept : xyz(simd::pack2x2(v, simd::set2(z))) {}
	//explicit Vector3(const tuples::templates::Tuple3<float>& t) noexcept : xyz(simd::set4(t.x, t.y, t.z, t.z)) {}
	//template<typename U> explicit Vector3(const tuples::templates::Tuple3<U>& t) noexcept : Vector3((float)t.x, (float)t.y, (float)t.z) {}
	explicit Vector3(const std::tuple<float, float, float>& t) noexcept : Vector3(std::get<0>(t), std::get<1>(t), std::get<2>(t)) {}
	template<typename U> explicit Vector3(const std::tuple<U, U, U>& t) noexcept : Vector3((float)std::get<0>(t), (float)std::get<1>(t), (float)std::get<2>(t)) {}
	explicit Vector3(const float* v) noexcept : Vector3(v[0], v[1], v[2]) {}
#else
	/*constexpr*/ explicit Vector3(float scalar) noexcept : xyz(simd::set3(scalar)) {}
	/*constexpr*/ Vector3(float x, float y, float z) noexcept : xyz(simd::set3(x, y, z)) {}
	/*constexpr*/ Vector3(const Vector2<float>& v) noexcept : xyz(simd::cutoff2(v)) {}
	/*constexpr*/ Vector3(const Vector2<float>& v, float z) noexcept : xyz(simd::pack2x2(v, simd::set1(z))) {}
	//explicit Vector3(const tuples::templates::Tuple3<float>& t) noexcept : xyz(simd::set3(t.x, t.y, t.z)) {}
	//template<typename U> explicit Vector3(const tuples::templates::Tuple3<U>& t) noexcept : xyz(simd::set3((float)t.x, (float)t.y, (float)t.z)) {}
	explicit Vector3(const std::tuple<float, float, float>& t) noexcept : xyz(simd::set3(std::get<0>(t), std::get<1>(t), std::get<2>(t))) {}
	template<typename U> explicit Vector3(const std::tuple<U, U, U>& t) noexcept : xyz(simd::set3((float)std::get<0>(t), (float)std::get<1>(t), (float)std::get<2>(t))) {}
	explicit Vector3(const float* v) noexcept : xyz(simd::load3(v)) {}
#endif
	explicit Vector3(Axis axis) noexcept : Vector3((axis == Axis::X) ? 1.f : 0.f, (axis == Axis::Y) ? 1.f : 0.f, (axis == Axis::Z) ? 1.f : 0.f) {}
	explicit Vector3(simd::float4 v) noexcept : xyz(v) {}
	Vector3(const Vector3& v) noexcept : xyz(v.xyz) {}
	template<typename U> explicit Vector3(const Vector3<U>& v) noexcept : Vector3((float)v.x, (float)v.y, (float)v.z) {}
	Vector3& operator=(const Vector3& v) noexcept { xyz = v.xyz; return *this; }

	operator simd::float4() const noexcept { return xyz; }
	//explicit operator tuples::templates::Tuple3<float>() noexcept { return tuples::templates::Tuple3<float>(x, y, z); }
	//template<typename U> explicit operator tuples::templates::Tuple3<U>() noexcept { return tuples::templates::Tuple3<U>(U(x), U(y), U(z)); }
	//explicit operator std::tuple<float, float, float>() { return std::tuple<float, float, float>(x, y, z); }
	//template<typename U> explicit operator std::tuple<U, U, U>() { return std::tuple<U, U, U>(U(x), U(y), U(z)); }
	explicit operator float* () noexcept { return &x; }
	explicit operator const float* () const noexcept { return &x; }
	float& operator[](int i) noexcept { return (&x)[i]; }
	const float& operator[](int i) const noexcept { return (&x)[i]; }

	Vector3 operator+() const noexcept { return *this; }
#if MATHEMATICS_SIMD_EXPAND_LAST
	Vector3 operator-() const noexcept { return Vector3(simd::neg4(xyz)); }
#else
	Vector3 operator-() const noexcept { return Vector3(simd::neg3(xyz)); }
#endif
	Vector3& operator+=(const Vector3& v) noexcept { xyz = simd::add4(xyz, v); return *this; }
	Vector3& operator-=(const Vector3& v) noexcept { xyz = simd::sub4(xyz, v); return *this; }
	Vector3& operator*=(const Vector3& v) noexcept { xyz = simd::mul4(xyz, v); return *this; }
	Vector3& operator*=(float f) noexcept { xyz = simd::mul4(xyz, simd::set4(f)); return *this; }
	Vector3& operator*=(const Matrix3<float>& m) noexcept;
#if MATHEMATICS_SIMD_EXPAND_LAST
	Vector3& operator/=(const Vector3& v) noexcept { xyz = simd::div4(xyz, v); return *this; }
#else
	Vector3& operator/=(const Vector3& v) noexcept { xyz = simd::div3(xyz, v); return *this; }
#endif
	Vector3& operator/=(float f) noexcept { xyz = simd::div4(xyz, simd::set4(f)); return *this; }
	bool operator==(const Vector3& v) const noexcept { return simd::all3(simd::equal(xyz, v)); }
	bool operator!=(const Vector3& v) const noexcept { return !(*this == v); }

	template<typename A> void load(A& ar);
	template<typename A> void save(A& ar) const { ar(x, y, z); }

	template<std::size_t I> float& get() noexcept;
	template<std::size_t I> const float& get() const noexcept;
	template<typename U> U& get() noexcept;				// intentionally undefined
	template<typename U> const U& get() const noexcept;	// intentionally undefined
	template<> simd::float4& get() noexcept { return xyz; }
	template<> const simd::float4& get() const noexcept { return xyz; }

#if MATHEMATICS_SIMD_EXPAND_LAST
	const Vector2<float> xy/*getXY*/() const noexcept { return Vector2<float>(simd::xyyy(xyz)); }
	const Vector2<float> xz/*getXZ*/() const noexcept { return Vector2<float>(simd::xzzz(xyz)); }
	const Vector2<float> zy/*getZY*/() const noexcept { return Vector2<float>(simd::zyyy(xyz)); }
#else
	const Vector2<float> xy/*getXY*/() const noexcept { return Vector2<float>(simd::cutoff2(xyz)); }
	const Vector2<float> xz/*getXZ*/() const noexcept { return Vector2<float>(simd::xzww(xyz)); }
	const Vector2<float> zy/*getZY*/() const noexcept { return Vector2<float>(simd::zyww(xyz)); }
#endif
	//void setXY(const Vector2<float>& v) noexcept { xyz = simd::insert2(v, xyz); }
	bool isZero() const noexcept { return simd::all3(simd::equal(xyz, simd::zero<simd::float4>())); }
	bool isApproxZero() const noexcept { return simd::all3(simd::lessThan(simd::abs4(xyz), TOLERANCE)); }
	bool isApproxEqual(const Vector3& v) const noexcept { return simd::all3(simd::lessThan(simd::abs4(simd::sub4(xyz, v)), TOLERANCE)); }
	bool isApproxEqual(const Vector3& v, float tolerance) const noexcept { return simd::all3(simd::lessThan(simd::abs4(simd::sub4(xyz, v)), simd::set4(tolerance))); }
	bool allLessThan(const Vector3& v) const noexcept { return simd::all3(simd::lessThan(xyz, v)); }
	bool allLessThanEqual(const Vector3& v) const noexcept { return simd::all3(simd::lessThanEqual(xyz, v)); }
	bool allGreaterThan(const Vector3& v) const noexcept { return simd::all3(simd::greaterThan(xyz, v)); }
	bool allGreaterThanEqual(const Vector3& v) const noexcept { return simd::all3(simd::greaterThanEqual(xyz, v)); }
	bool anyLessThan(const Vector3& v) const noexcept { return simd::any3(simd::lessThan(xyz, v)); }
	bool anyLessThanEqual(const Vector3& v) const noexcept { return simd::any3(simd::lessThanEqual(xyz, v)); }
	bool anyGreaterThan(const Vector3& v) const noexcept { return simd::any3(simd::greaterThan(xyz, v)); }
	bool anyGreaterThanEqual(const Vector3& v) const noexcept { return simd::any3(simd::greaterThanEqual(xyz, v)); }
	bool isFinite() const noexcept { return simd::all3(simd::isFinite(xyz)); }
	float getMagnitude() const noexcept { return simd::toFloat(simd::sqrt1(simd::dot3(xyz, xyz))); }
	float getMagnitudeSquared() const noexcept { return simd::toFloat(simd::dot3(xyz, xyz)); }
	void setMagnitude(float magnitude) noexcept;
	float getLength() const noexcept { return getMagnitude(); }
	float getLengthSquared() const noexcept { return getMagnitudeSquared(); }
	void setLength(float length) noexcept { setMagnitude(length); }
	Axis getMajorAxis() const noexcept { (Axis)simd::asIndex(simd::equal(xyz, simd::hMax3(xyz))); }
	float getMinComponent() const noexcept { return simd::toFloat(simd::hMin3(xyz)); }
	float getMaxComponent() const noexcept { return simd::toFloat(simd::hMax3(xyz)); }
	Vector3& setZero/*zero*/() noexcept { xyz = simd::zero<simd::float4>(); return *this; }
#if MATHEMATICS_SIMD_EXPAND_LAST
	Vector3& set(float x, float y, float z) noexcept { xyz = simd::set4(x, y, z, z); return *this; }
#else
	Vector3& set(float x, float y, float z) noexcept { xyz = simd::set3(x, y, z); return *this; }
#endif
	Vector3& setMinimum(const Vector3& v1, const Vector3& v2) noexcept { xyz = simd::min4(v1, v2); return *this; }
	Vector3& setMaximum(const Vector3& v1, const Vector3& v2) noexcept { xyz = simd::max4(v1, v2); return *this; }
#if MATHEMATICS_SIMD_EXPAND_LAST
	Vector3& negate() noexcept { xyz = simd::neg4(xyz); return *this; }
#else
	Vector3& negate() noexcept { xyz = simd::neg3(xyz); return *this; }
#endif
	Vector3& normalize() noexcept;
	Vector3& rotate(Axis axis, float angle) noexcept;
	Vector3& rotate(const Quaternion<float>& q) noexcept;
	//Vector3& scale(const Vector3& v) noexcept { xyz = simd::mul4(xyz, v); return *this; }
	Vector3& transform(const Matrix3<float>& m) noexcept;
	Vector3& transform(const AffineTransform<float>& m) noexcept;

	//static const Vector3& getZero() noexcept { return ZERO; }
	//static const Vector3& getUnitX() noexcept { return UNIT_X; }
	//static const Vector3& getUnitY() noexcept { return UNIT_Y; }
	//static const Vector3& getUnitZ() noexcept { return UNIT_Z; }

	static const Vector3 ZERO;
	static const Vector3 UNIT_X;
	static const Vector3 UNIT_Y;
	static const Vector3 UNIT_Z;
	static const Vector3 ONE;
	static const Vector3 TOLERANCE;
	static const Vector3 INF;
	static const Vector3 MINUS_INF;

	union
	{
		simd::float4 xyz;
		struct { float x, y, z/*, reserved_*/; };
	};
};

const Vector3<float> Vector3<float>::ZERO{};
const Vector3<float> Vector3<float>::UNIT_X{ 1.f, 0.f, 0.f };
const Vector3<float> Vector3<float>::UNIT_Y{ 0.f, 1.f, 0.f };
const Vector3<float> Vector3<float>::UNIT_Z{ 0.f, 0.f, 1.f };
const Vector3<float> Vector3<float>::ONE{ 1.f, 1.f, 1.f };
const Vector3<float> Vector3<float>::TOLERANCE{ Constants<float>::TOLERANCE, Constants<float>::TOLERANCE, Constants<float>::TOLERANCE };
const Vector3<float> Vector3<float>::INF{ std::numeric_limits<float>::infinity(), std::numeric_limits<float>::infinity(), std::numeric_limits<float>::infinity() };
const Vector3<float> Vector3<float>::MINUS_INF{ -std::numeric_limits<float>::infinity(), -std::numeric_limits<float>::infinity(), -std::numeric_limits<float>::infinity() };

#endif /* SIMD_HAS_FLOAT4 */

template<typename T>
	requires (std::floating_point<T> || std::integral<T>)
inline Vector3<T> operator+(const Vector3<T>& v1, const Vector3<T>& v2) noexcept
{ 
	return Vector3<T>(v1.x + v2.x, v1.y + v2.y, v1.z + v2.z);
}

template<typename T>
	requires (std::floating_point<T> || std::integral<T>)
inline Vector3<T> operator-(const Vector3<T>& v1, const Vector3<T>& v2) noexcept
{ 
	return Vector3<T>(v1.x - v2.x, v1.y - v2.y, v1.z - v2.z); 
}

template<typename T>
	requires (std::floating_point<T> || std::integral<T>)
inline Vector3<T> operator*(const Vector3<T>& v1, const Vector3<T>& v2) noexcept
{ 
	return Vector3<T>(v1.x*v2.x, v1.y*v2.y, v1.z*v2.z);
}

template<typename T>
	requires (std::floating_point<T> || std::integral<T>)
inline Vector3<T> operator*(T f, const Vector3<T>& v) noexcept
{ 
	return Vector3<T>(f*v.x, f*v.y, f*v.z); 
}

template<typename T>
	requires (std::floating_point<T> || std::integral<T>)
inline Vector3<T> operator*(const Vector3<T>& v, T f) noexcept
{ 
	return Vector3<T>(v.x*f, v.y*f, v.z*f); 
}

template<typename T>
	requires (std::floating_point<T> || std::integral<T>)
inline Vector3<T> operator/(const Vector3<T>& v1, const Vector3<T>& v2) noexcept
{ 
	return Vector3<T>(v1.x/v2.x, v1.y/v2.y, v1.z/v2.z); 
}

template<typename T>
	requires (std::floating_point<T> || std::integral<T>)
inline Vector3<T> operator/(T f, const Vector3<T>& v) noexcept
{ 
	return Vector3<T>(f/v.x, f/v.y, f/v.z);
}

template<typename T>
	requires (std::floating_point<T> || std::integral<T>)
inline Vector3<T> operator/(const Vector3<T>& v, T f) noexcept
{ 
	if constexpr (std::is_floating_point_v<T>)
		return operator*(v, T(1)/f);
	else
		return Vector3<T>(v.x/f, v.y/f, v.z/f);
}

template<typename C, typename T, typename U>
	requires std::floating_point<U> || std::integral<U>
inline std::basic_istream<C, T>& operator>>(std::basic_istream<C, T>& s, Vector3<U>& v)
{ 
	return s >> v.x >> std::ws >> v.y >> std::ws >> v.z; 
}

template<typename C, typename T, typename U>
	requires std::floating_point<U> || std::integral<U>
inline std::basic_ostream<C, T>& operator<<(std::basic_ostream<C, T>& s, const Vector3<U>& v)
{ 
	constexpr C WS(0x20);
	return s << v.x << WS << v.y << WS << v.z;
}

template<std::floating_point T>
template<std::size_t I>
inline T& Vector3<T>::get()
{
	if constexpr (I == 0)
		return x;
	else if constexpr (I == 1)
		return y;
	else if constexpr (I == 2)
		return z;
	static_assert(false);
}

template<std::floating_point T>
template<std::size_t I>
inline const T& Vector3<T>::get() const
{
	if constexpr (I == 0)
		return x;
	else if constexpr (I == 1)
		return y;
	else if constexpr (I == 2)
		return z;
	static_assert(false);
}

template<std::floating_point T>
inline bool Vector3<T>::isApproxZero() const
{ 
	return (std::fabs(x) < Constants<T>::TOLERANCE) && (std::fabs(y) < Constants<T>::TOLERANCE) && 
		(std::fabs(z) < Constants<T>::TOLERANCE);
}

template<std::floating_point T>
inline bool Vector3<T>::isApproxEqual(const Vector3<T>& v) const
{ 
	return (std::fabs(v.x - x) < Constants<T>::TOLERANCE) && (std::fabs(v.y - y) < Constants<T>::TOLERANCE) &&
		(std::fabs(v.z - z) < Constants<T>::TOLERANCE);
}

template<std::floating_point T>
inline bool Vector3<T>::isApproxEqual(const Vector3<T>& v, T tolerance) const
{ 
	return (std::fabs(v.x - x) < tolerance) && (std::fabs(v.y - y) < tolerance) &&
		(std::fabs(v.z - z) < tolerance);
}

template<std::floating_point T>
inline void Vector3<T>::setMagnitude(T magnitude) 
{ 
	T m = getMagnitude(); 
	if (m > T(0)) 
		*this *= magnitude/m;
}

template<std::floating_point T>
inline Axis Vector3<T>::getMajorAxis() const
{ 
	Axis axis = Axis::X; 
	if (std::fabs(y) > std::fabs(x)) 
		axis = Axis::Y;
	if (std::fabs(z) > std::fabs((&x)[(int)axis])) 
		axis = Axis::Z; 
	return axis;
}

template<std::floating_point T>
inline Vector3<T>& Vector3<T>::setMinimum(const Vector3<T>& v1, const Vector3<T>& v2)
{ 
	x = std::min(v1.x, v2.x); 
	y = std::min(v1.y, v2.y);
	z = std::min(v1.z, v2.z); 
	return *this;
}

template<std::floating_point T>
inline Vector3<T>& Vector3<T>::setMaximum(const Vector3<T>& v1, const Vector3<T>& v2)
{ 
	x = std::max(v1.x, v2.x); 
	y = std::max(v1.y, v2.y);
	z = std::max(v1.z, v2.z); 
	return *this;
}

template<std::floating_point T>
inline Vector3<T>& Vector3<T>::normalize()
{
//#if MATHEMATICS_FAST_NORMALIZE
//	if costexpr(std::is_same_v<T, float>)
//	{
//		T m = rcpSqrtApprox(getMagnitudeSquared());
//		if (m <= std::numeric_limits<T>::max()) 
//			*this *= m;
//	}
//	else
//#endif
	{
		T m = getMagnitude();
		if (m > T(0))
			*this /= m;
	}
	return *this;
}

template<std::floating_point T>
inline Vector3<T>& Vector3<T>::rotate(Axis axis, T angle)
{
	if (angle != T(0))
	{
		T sinAngle = std::sin(angle);
		T cosAngle = std::cos(angle);

		switch (axis)
		{
			case Axis::X:
			{
				T j = y, k = z;
				y = j*cosAngle - k*sinAngle;
				z = k*cosAngle + j*sinAngle;
			} break;

			case Axis::Y:
			{
				T i = x, k = z;
				x = i*cosAngle + k*sinAngle;
				z = k*cosAngle - i*sinAngle;
			} break;

			case Axis::Z:
			{
				T i = x, j = y;
				x = i*cosAngle - j*sinAngle;
				y = j*cosAngle + i*sinAngle;
			} break;
		}
	}

	return *this;
}

template<std::integral T>
template<std::size_t I>
inline T& Vector3<T>::get()
{
	if constexpr (I == 0)
		return x;
	else if constexpr (I == 1)
		return y;
	else if constexpr (I == 2)
		return z;
	static_assert(false);
}

template<std::integral T>
template<std::size_t I>
inline const T& Vector3<T>::get() const
{
	if constexpr (I == 0)
		return x;
	else if constexpr (I == 1)
		return y;
	else if constexpr (I == 2)
		return z;
	static_assert(false);
}

template<std::integral T>
inline Vector3<T>& Vector3<T>::setMinimum(const Vector3<T>& v1, const Vector3<T>& v2)
{
	x = std::min(v1.x, v2.x);
	y = std::min(v1.y, v2.y);
	z = std::min(v1.z, v2.z);
	return *this;
}

template<std::integral T>
inline Vector3<T>& Vector3<T>::setMaximum(const Vector3<T>& v1, const Vector3<T>& v2)
{
	x = std::max(v1.x, v2.x);
	y = std::max(v1.y, v2.y);
	z = std::max(v1.z, v2.z);
	return *this;
}

#if SIMD_HAS_FLOAT4

template<>
inline Vector3<float> operator+(const Vector3<float>& v1, const Vector3<float>& v2) noexcept 
{ 
	return Vector3<float>(simd::add4(v1, v2)); 
}

template<>
inline Vector3<float> operator-(const Vector3<float>& v1, const Vector3<float>& v2) noexcept 
{ 
	return Vector3<float>(simd::sub4(v1, v2));
}

template<>
inline Vector3<float> operator*(const Vector3<float>& v1, const Vector3<float>& v2) noexcept 
{ 
	return Vector3<float>(simd::mul4(v1, v2)); 
}

template<>
inline Vector3<float> operator*(float f, const Vector3<float>& v) noexcept 
{ 
	return Vector3<float>(simd::mul4(simd::set4(f), v)); 
}

template<>
inline Vector3<float> operator*(const Vector3<float>& v, float f) noexcept 
{ 
	return Vector3<float>(simd::mul4(v, simd::set4(f))); 
}

template<>
inline Vector3<float> operator/(const Vector3<float>& v1, const Vector3<float>& v2) noexcept 
{ 
#if IMAGING_SIMD_EXPAND_LAST
	return Vector3<float>(simd::div4(v1, v2));
#else
	return Vector3<float>(simd::div3(v1, v2));
#endif
}

template<>
inline Vector3<float> operator/(float f, const Vector3<float>& v) noexcept 
{ 
#if IMAGING_SIMD_EXPAND_LAST
	return Vector3<float>(simd::div4(simd::set4(f), v));
#else
	return Vector3<float>(simd::div3(simd::set3(f), v));
#endif
}

template<>
inline Vector3<float> operator/(const Vector3<float>& v, float f) noexcept 
{ 
	return Vector3<float>(simd::div4(v, simd::set4(f))); 
}

template<typename C, typename T>
inline std::basic_istream<C, T>& operator>>(std::basic_istream<C, T>& s, Vector3<float>& v)
{ 
	float x, y, z; 
	s >> x >> std::ws >> y >> std::ws >> z; 
	v.set(x, y, z); 
	return s; 
}

template<typename A>
inline void Vector3<float>::load(A& ar)
{
	float t0, t1, t2;
	ar(t0, t1, t2);
	set(t0, t1, t2);
}

template<std::size_t I>
inline float& Vector3<float>::get()
{
	if constexpr (I == 0)
		return x;
	else if constexpr (I == 1)
		return y;
	else if constexpr (I == 2)
		return z;
	static_assert(false);
}

template<std::size_t I>
inline const float& Vector3<float>::get() const
{
	if constexpr (I == 0)
		return x;
	else if constexpr (I == 1)
		return y;
	else if constexpr (I == 2)
		return z;
	static_assert(false);
}

inline void Vector3<float>::setMagnitude(float magnitude)
{ 
	float m = getMagnitude(); 
	if (m > 0.f) 
		*this *= magnitude/m;
}

inline Vector3<float>& Vector3<float>::normalize()
{
#if MATHEMATICS_FAST_NORMALIZE
	float m = simd::toFloat(simd::rcpSqrtApprox1(simd::dot3(xyz, xyz)));
	if (m <= std::numeric_limits<float>::max()) 
		xyz = simd::mul4(xyz, simd::set4(m));
#else
	float m = getMagnitude(); 
	if (m > 0.f) 
		xyz = simd::div4(xyz, simd::set4(m));
#endif
	return *this;
}

inline Vector3<float>& Vector3<float>::rotate(Axis axis, float angle)
{
	if (angle != 0.f)
	{
		float sinAngle = std::sin(angle);
		float cosAngle = std::cos(angle);

		switch (axis)
		{
			case Axis::X:
			{
				float j = y, k = z;
				set(x, j*cosAngle - k*sinAngle, k*cosAngle + j*sinAngle); // #TODO SIMD
			} break;

			case Axis::Y:
			{
				float i = x, k = z;
				set(i*cosAngle + k*sinAngle, y, k*cosAngle - i*sinAngle); // #TODO SIMD
			} break;

			case Axis::Z:
			{
				float i = x, j = y;
				set(i*cosAngle - j*sinAngle, j*cosAngle + i*sinAngle, z); // #TODO SIMD
			} break;
		}
	}

	return *this;
}

#endif /* SIMD_HAS_FLOAT4 */

template<std::size_t I, typename T>
	requires (std::floating_point<T> || std::integral<T>)
inline T& get(Vector3<T>& v) noexcept
{
	if constexpr (I == 0)
		return v.x;
	else if constexpr (I == 1)
		return v.y;
	else if constexpr (I == 2)
		return v.z;
	static_assert(false);
}

template<std::size_t I, typename T>
	requires (std::floating_point<T> || std::integral<T>)
inline const T& get(const Vector3<T>& v) noexcept
{
	if constexpr (I == 0)
		return v.x;
	else if constexpr (I == 1)
		return v.y;
	else if constexpr (I == 2)
		return v.z;
	static_assert(false);
}

template<std::size_t I, typename T>
	requires (std::floating_point<T> || std::integral<T>)
inline T&& get(Vector3<T>&& v) noexcept
{
	if constexpr (I == 0)
		return v.x;
	else if constexpr (I == 1)
		return v.y;
	else if constexpr (I == 2)
		return v.z;
	static_assert(false);
}

template<std::size_t I, typename T>
	requires (std::floating_point<T> || std::integral<T>)
inline const T&& get(const Vector3<T>&& v) noexcept
{
	if constexpr (I == 0)
		return v.x;
	else if constexpr (I == 1)
		return v.y;
	else if constexpr (I == 2)
		return v.z;
	static_assert(false);
}

template<typename T>
	requires std::floating_point<T>
inline T dot(const Vector3<T>& v1, const Vector3<T>& v2) noexcept
{
	return v1.x*v2.x + v1.y*v2.y + v1.z*v2.z;
}

template<typename T>
	requires std::floating_point<T>
inline Vector3<T> cross(const Vector3<T>& v1, const Vector3<T>& v2) noexcept
{
	return Vector3(v1.y*v2.z - v1.z*v2.y, v1.z*v2.x - v1.x*v2.z, v1.x*v2.y - v1.y*v2.x);
}

template<typename T>
	requires std::floating_point<T>
inline T distance(const Vector3<T>& v1, const Vector3<T>& v2) noexcept
{
	T x = v2.x - v1.x;
	T y = v2.y - v1.y;
	T z = v2.z - v1.z;
	return std::sqrt(x*x + y*y + z*z);
}

template<typename T>
	requires std::floating_point<T>
inline T distanceSquared(const Vector3<T>& v1, const Vector3<T>& v2) noexcept
{
	T x = v2.x - v1.x;
	T y = v2.y - v1.y;
	T z = v2.z - v1.z;
	return x*x + y*y + z*z;
}

template<typename T>
	requires std::floating_point<T>
inline T length(const Vector3<T>& v) noexcept
{
	return std::sqrt(v.x*v.x + v.y*v.y + v.z*v.z);
}

template<typename T>
	requires std::floating_point<T>
inline T lengthSquared(const Vector3<T>& v) noexcept
{
	return v.x*v.x + v.y*v.y + v.z*v.z;
}

template<typename T>
	requires (std::floating_point<T> || std::integral<T>)
inline Vector3<T> minimum(const Vector3<T>& v1, const Vector3<T>& v2)
{
	return Vector3<T>(std::min(v1.x, v2.x), std::min(v1.y, v2.y), std::min(v1.z, v2.z));
}

template<typename T>
	requires (std::floating_point<T> || std::integral<T>)
inline Vector3<T> maximum(const Vector3<T>& v1, const Vector3<T>& v2)
{
	return Vector3<T>(std::max(v1.x, v2.x), std::max(v1.y, v2.y), std::max(v1.z, v2.z));
}

template<typename T>
	requires std::floating_point<T>
inline Vector3<T> lerp(const Vector3<T>& v1, const Vector3<T>& v2, T t) noexcept
{
	return Vector3<T>(v1.x + t*(v2.x - v1.x), v1.y + t*(v2.y - v1.y), v1.z + t*(v2.z - v1.z));
}

template<typename T>
	requires std::floating_point<T>
inline Vector3<T> slerp(const Vector3<T>& v1, const Vector3<T>& v2, T t)
{
	T dp = v1.x*v2.x + v1.y*v2.y + v1.z*v2.z;
	if ((T(1) - dp) < Constants<T>::TOLERANCE)
	{
		Vector3<T> c(v1.x + t*(v2.x - v1.x), v1.y + t*(v2.y - v1.y), v1.z + t*(v2.z - v1.z));
		c.normalize();
		return c;
	}

	dp = std::clamp(dp, T(-1), T(1));
	T theta = std::acos(dp)*t;
	T st = std::sin(theta);
	T ct = std::cos(theta);
	Vector3<T> c(v2.x - v1.x*dp, v2.y - v1.y*dp, v2.z - v1.z*dp);
	c.normalize();
	return Vector3<T>(v1.x*ct + c.x*st, v1.y*ct + c.y*st, v1.z*ct + c.z*st);
}

#if SIMD_HAS_FLOAT4

template<>
inline float dot(const Vector3<float>& v1, const Vector3<float>& v2) noexcept
{
	return simd::toFloat(simd::dot3(v1, v2));
}

template<>
inline Vector3<float> cross(const Vector3<float>& v1, const Vector3<float>& v2) noexcept
{
	return Vector3<float>(simd::sub4(simd::mul4(simd::yzxw(v1), simd::zxyw(v2)),
		simd::mul4(simd::zxyw(v1), simd::yzxw(v2))));
}

template<>
inline float distance(const Vector3<float>& v1, const Vector3<float>& v2) noexcept
{
	auto v = simd::sub4(v2, v1);
	return simd::toFloat(simd::sqrt1(simd::dot3(v, v)));
}

template<>
inline float distanceSquared(const Vector3<float>& v1, const Vector3<float>& v2) noexcept
{
	auto v = simd::sub4(v2, v1);
	return simd::toFloat(simd::dot3(v, v));
}

template<>
inline float length(const Vector3<float>& v) noexcept
{
	return simd::toFloat(simd::sqrt1(simd::dot3(v, v)));
}

template<>
inline float lengthSquared(const Vector3<float>& v) noexcept
{
	return simd::toFloat(simd::dot3(v, v));
}

template<>
inline Vector3<float> minimum(const Vector3<float>& v1, const Vector3<float>& v2)
{
	return Vector3<float>(simd::min4(v1, v2));
}

template<>
inline Vector3<float> maximum(const Vector3<float>& v1, const Vector3<float>& v2)
{
	return Vector3<float>(simd::max4(v1, v2));
}

template<>
inline Vector3<float> lerp(const Vector3<float>& v1, const Vector3<float>& v2, float t) noexcept
{
	return Vector3<float>(simd::mulAdd4(simd::set4(t), simd::sub4(v2, v1), v1));
}

template<>
inline Vector3<float> slerp(const Vector3<float>& v1, const Vector3<float>& v2, float t)
{
	float dp = simd::toFloat(simd::dot3(v1, v2));
	if ((1.f - dp) < Constants<float>::TOLERANCE)
	{
		Vector3<float> c(simd::mulAdd4(simd::set4(t), simd::sub4(v2, v1), v1));
		c.normalize();
		return c;
	}

	dp = std::clamp(dp, -1.f, 1.f);
	float theta = std::acos(dp)*t;
	float st = std::sin(theta);
	float ct = std::cos(theta);
	Vector3<float> c(simd::sub4(v2, simd::mul4(v1, simd::set4(dp))));
	c.normalize();
	return Vector3<float>(simd::mulAdd4(v1, simd::set4(ct), simd::mul4(c, simd::set4(st))));
}

#endif /* SIMD_HAS_FLOAT4 */

template<typename T>
	requires std::floating_point<T>
inline Vector3<T> perpendicular(const Vector3<T>& v) noexcept // #TODO SIMD
{
	int min = 0, a = 1, b = 2;

	if (std::fabs(v.y) < std::fabs(v.x)) // #TODO Use getMinComponent()
	{
		a = min;
		min = 1;
	}

	if (std::fabs(v.z) < std::fabs(v[min]))
	{
		b = min;
		min = 2;
	}

	Vector3<T> u/*(Uninitialized())*/;
	//u[min] = T(0);
	u[a] = v[b];
	u[b] = -v[a];
	return u*(v.getMagnitude()/u.getMagnitude()); // #FIXME division by zero
}

template<typename T>
	requires std::floating_point<T>
inline Vector3<T> reflect(const Vector3<T>& i, const Vector3<T>& n) noexcept
{
	return i - T(2)*dot(i, n)*n;
}

template<typename T>
	requires std::floating_point<T>
inline Vector3<T> refract(const Vector3<T>& i, const Vector3<T>& n, T etaRatio) noexcept
{
	T cosI = dot(-i, n);
	T cosT2 = T(1) - etaRatio*etaRatio*(T(1) - cosI*cosI);
	return (cosT2 <= T(0)) ? Vector3<T>() : (etaRatio*i + ((etaRatio*cosI - std::sqrt(std::fabs(cosT2)))*n));
}

template<typename T>
	requires std::floating_point<T>
inline Vector3<T> normalize(const Vector3<T>& v) noexcept
{
	Vector3<T> u(v);
	u.normalize();
	return u;
}

template<typename T>
#if SIMD_HAS_FLOAT4
	requires (std::floating_point<T> && !std::same_as<T, float>)
#else
	requires std::floating_point<T>
#endif
inline Vector3<T> normalize(Vector3<T>&& v) noexcept
{
	v.normalize();
	return v;
}

template<typename T>
	requires std::floating_point<T>
inline Vector3<T> project(const Vector3<T>& v1, const Vector3<T>& v2) noexcept
{
	T m = v2.getMagnitudeSquared();
	return (m > T(0)) ? (dot(v1, v2)/m)*v2 : Vector4<T>::ZERO;
}

template<typename T>
	requires std::floating_point<T>
inline T angle(const Vector3<T>& v1, const Vector3<T>& v2)
{
#if MATHEMATICS_FAST_NORMALIZE
	Vector3<T> u(v1);
	Vector3<T> v(v2);
	u.normalize();
	v.normalize();
	return std::acos(std::clamp(dot(u, v), T(-1), T(1)));
#else
	T q = v1.getMagnitude()*v2.getMagnitude();
	if (q > T(0))
		return std::acos(std::clamp(dot(v1, v2)/q, T(-1), T(1)));
	else
		return (dot(v1, v2) >= T(0)) ? T(0) : Constants<T>::PI;
#endif
}

} // namespace templates

using Vector3f = templates::Vector3<float>;
using Vector3d = templates::Vector3<double>;

#if MATHEMATICS_DOUBLE
using Vector3 = templates::Vector3<double>;
using Vector3Arg = templates::Vector3<double>::ConstArg;
using Vector3Result = templates::Vector3<double>::ConstResult;
#else
using Vector3 = templates::Vector3<float>;
using Vector3Arg = templates::Vector3<float>::ConstArg;
using Vector3Result = templates::Vector3<float>::ConstResult;
#endif

using IntVector3 = templates::Vector3<int>;
using IntVector3Arg = templates::Vector3<int>::ConstArg;
using IntVector3Result = templates::Vector3<int>::ConstResult;

} // namespace core::mathematics

namespace std {

template<size_t I, typename T>
struct tuple_element;

template<size_t I, typename T>
struct tuple_element<I, ::core::mathematics::templates::Vector3<T>>
{
	using type = T;
};

template<typename T>
struct tuple_size;

template<typename T>
struct tuple_size<::core::mathematics::templates::Vector3<T>> : integral_constant<size_t, 3> 
{
};

template<typename T>
struct hash;

template<typename T>
struct hash<::core::mathematics::templates::Vector3<T>>
{
	std::size_t operator()(const ::core::mathematics::templates::Vector3<T>& v) const noexcept
	{
		std::hash<T> hasher;
		std::size_t seed = hasher(v.x) + 0x9e3779b9;
		seed ^= hasher(v.y) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
		seed ^= hasher(v.z) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
		return seed;
	}
};

} // namespace std

#include "Matrix3.hpp"
#include "AffineTransform.hpp"
#include "Quaternion.hpp"

namespace core::mathematics::templates {

template<typename T>
inline Vector3<T>& Vector3<T>::operator*=(const Matrix3<T>& m)
{
	set(x*m.m00 + y*m.m10 + z*m.m20, x*m.m01 + y*m.m11 + z*m.m21, x*m.m02 + y*m.m12 + z*m.m22);
	return *this;
}

template<typename T>
	requires std::floating_point<T>
inline Vector3<T> operator*(const Vector3<T>& v, const Matrix3<T>& m) noexcept
{
	return Vector3<T>(v.x*m.m00 + v.y*m.m10 + v.z*m.m20, v.x*m.m01 + v.y*m.m11 + v.z*m.m21, v.x*m.m02 + v.y*m.m12 + v.z*m.m22);
}

//template<typename T>
//	requires std::floating_point<T>
//inline Vector3<T> operator*(const Matrix3<T>& m, const Vector3<T>& v) noexcept; // valid for column vectors only

template<typename T>
inline Vector3<T>& Vector3<T>::rotate(const Quaternion<T>& q)
{
	T x1 = q.y*z - q.z*y;
	T y1 = q.z*x - q.x*z;
	T z1 = q.x*y - q.y*x;
	T x2 = q.w*x1 + q.y*z1 - q.z*y1;
	T y2 = q.w*y1 + q.z*x1 - q.x*z1;
	T z2 = q.w*z1 + q.x*y1 - q.y*x1;
	x += T(2)*x2;
	y += T(2)*y2;
	z += T(2)*z2;
	return *this;
}

template<typename T>
inline Vector3<T>& Vector3<T>::transform(const Matrix3<T>& m)
{
	*this *= m;
	return *this;
}

template<typename T>
inline Vector3<T>& Vector3<T>::transform(const AffineTransform<T>& m)
{
	set(x*m.m00 + y*m.m10 + z*m.m20 + m.m30, x*m.m01 + y*m.m11 + z*m.m21 + m.m31, x*m.m02 + y*m.m12 + z*m.m22 + m.m32);
	return *this;
}

template<typename T>
	requires std::floating_point<T>
inline Vector3<T> rotate(const Vector3<T>& v, const Quaternion<T>& q) noexcept
{
	T x1 = q.y*v.z - q.z*v.y;
	T y1 = q.z*v.x - q.x*v.z;
	T z1 = q.x*v.y - q.y*v.x;
	T x2 = q.w*x1 + q.y*z1 - q.z*y1;
	T y2 = q.w*y1 + q.z*x1 - q.x*z1;
	T z2 = q.w*z1 + q.x*y1 - q.y*x1;
	return Vector3<T>(v.x + T(2)*x2, v.y + T(2)*y2, v.z + T(2)*z2);
}

template<typename T>
	requires std::floating_point<T>
inline Vector3<T> transform(const Vector3<T>& v, const Matrix3<T>& m) noexcept
{
	return v*m;
}

template<typename T>
	requires std::floating_point<T>
inline Vector3<T> transform(const Vector3<T>& v, const AffineTransform<T>& m) noexcept
{
	return Vector3<T>(v.x*m.m00 + v.y*m.m10 + v.z*m.m20 + m.m30, v.x*m.m01 + v.y*m.m11 + v.z*m.m21 + m.m31, 
		v.x*m.m02 + v.y*m.m12 + v.z*m.m22 + m.m32);
}

template<typename T>
	requires std::floating_point<T>
inline Matrix3<T> tensor(const Vector3<T>& v1, const Vector3<T>& v2) noexcept
{
	return Matrix3<T>(v1.x*v2.x, v1.x*v2.y, v1.x*v2.z,
		v1.y*v2.x, v1.y*v2.y, v1.y*v2.z,
		v1.z*v2.x, v1.z*v2.y, v1.z*v2.z);
}

#if SIMD_HAS_FLOAT4

inline Vector3<float>& Vector3<float>::operator*=(const Matrix3<float>& m)
{
	auto t = simd::mulAdd4(simd::xxxx(xyz), m.row0, simd::mul4(simd::yyyy(xyz), m.row1));
	xyz = simd::add4(t, simd::mul4(simd::zzzz(xyz), m.row2));
	return *this;
}

template<>
inline Vector3<float> operator*(const Vector3<float>& v, const Matrix3<float>& m) noexcept
{
	auto t = simd::mulAdd4(simd::xxxx(v), m.row0, simd::mul4(simd::yyyy(v), m.row1));
	return Vector3<float>(simd::add4(t, simd::mul4(simd::zzzz(v), m.row2)));
}

//template<>
//inline Vector3<float> operator*(const Matrix3<float>& m, const Vector3<float>& v) noexcept; // valid for column vectors only

inline Vector3<float>& Vector3<float>::rotate(const Quaternion<float>& q)
{
	auto qyzx = simd::yzxx(q);
	auto qzxy = simd::zxyy(q);
	auto t1 = simd::sub4(simd::mul4(qyzx, simd::zxyy(xyz)), simd::mul4(qzxy, simd::yzxx(xyz)));
	auto t2 = simd::mulAdd4(simd::wwww(q), t1, simd::sub4(simd::mul4(qyzx, simd::zxyy(t1)), simd::mul4(qzxy, simd::yzxx(t1))));
	static const simd::float4 two = simd::set3(2.0f);
	xyz = simd::mulAdd4(two, t2, xyz);
	return *this;
}

inline Vector3<float>& Vector3<float>::transform(const Matrix3<float>& m)
{
	*this *= m;
	return *this;
}

inline Vector3<float>& Vector3<float>::transform(const AffineTransform<float>& m)
{
	auto t = simd::mulAdd4(simd::xxxx(xyz), m.row0, simd::mul4(simd::yyyy(xyz), m.row1));
	t = simd::add4(t, simd::mul4(simd::zzzz(xyz), m.row2));
	xyz = simd::add4(t, m.row3);
	return *this;
}

template<>
inline Vector3<float> rotate(const Vector3<float>& v, const Quaternion<float>& q) noexcept
{
	auto qyzx = simd::yzxx(q);
	auto qzxy = simd::zxyy(q);
	auto t1 = simd::sub4(simd::mul4(qyzx, simd::zxyy(v)), simd::mul4(qzxy, simd::yzxx(v)));
	auto t2 = simd::mulAdd4(simd::wwww(q), t1, simd::sub4(simd::mul4(qyzx, simd::zxyy(t1)), simd::mul4(qzxy, simd::yzxx(t1))));
	static const simd::float4 two = simd::set3(2.0f);
	return Vector3<float>(simd::mulAdd4(two, t2, v));
}

template<>
inline Vector3<float> transform(const Vector3<float>& v, const Matrix3<float>& m) noexcept
{
	return v*m;
}

template<>
inline Vector3<float> transform(const Vector3<float>& v, const AffineTransform<float>& m) noexcept
{
	auto t = simd::mulAdd4(simd::xxxx(v), m.row0, simd::mul4(simd::yyyy(v), m.row1));
	t = simd::add4(t, simd::mul4(simd::zzzz(v), m.row2));
	return Vector3<float>(simd::add4(t, m.row3));
}

template<>
inline Matrix3<float> tensor(const Vector3<float>& v1, const Vector3<float>& v2) noexcept
{
	return Matrix3<float>(simd::mul4(simd::xxxx(v1), v2),
		simd::mul4(simd::yyyy(v1), v2),
		simd::mul4(simd::zzzz(v1), v2));
}

#endif /* SIMD_HAS_FLOAT4 */

} // namespace core::mathematics::templates
