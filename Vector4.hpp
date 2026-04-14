/*
 *	Name: Vector4
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
#include <tuple>
#include <cstddef>
#include <cmath>
//#include <Tuples/Tuple4.hpp>
#include <Simd/Intrinsics.hpp>
#include <Simd/Functional.hpp>
#include "Constants.hpp"
#include "Scalar.hpp"
#include "Axis.hpp"
#include "Vector2.hpp"
#include "Vector3.hpp"

namespace mathematics {
namespace templates {

template<typename T>
	requires std::floating_point<T>
struct Matrix4;

template<typename T>
	requires (std::floating_point<T> || std::integral<T>)
struct Vector4;

template<std::floating_point T>
struct Vector4<T>
{
	using Real = T;
	using ComponentType = T;
	using ConstArg = const Vector4&;
	using ConstResult = const Vector4&;
	using TupleType = std::tuple<T, T, T, T>;
	template<Arithmetic U> OtherTupleType = std::tuple<U, U, U, U>;

	static constexpr int NUM_COMPONENTS = 4;

	constexpr Vector4() noexcept : x(), y(), z(), w() {}
	explicit Vector4(Uninitialized) noexcept {}
	constexpr explicit Vector4(T scalar) noexcept : x(scalar), y(scalar), z(scalar), w(scalar) {}
	constexpr Vector4(T x, T y, T z, T w) noexcept : x(x), y(y), z(z), w(w) {}
	constexpr Vector4(const Vector2<T>& v) noexcept : x(v.x), y(v.y), z(), w(T(1)) {}
	constexpr Vector4(const Vector2<T>& v, T z, T w) noexcept : x(v.x), y(v.y), z(z), w(w) {}
	constexpr Vector4(const Vector2<T>& xy, const Vector2<T>& zw) noexcept : x(xy.x), y(xy.y), z(zw.x), w(zw.y) {}
	constexpr Vector4(const Vector3<T>& v) noexcept : x(v.x), y(v.y), z(v.z), w(T(1)) {}
	constexpr Vector4(const Vector3<T>& v, T w) noexcept : x(v.x), y(v.y), z(v.z), w(w) {}
	//explicit Vector4(const Tuple4<T>& t) noexcept : x(t.x), y(t.y), z(t.z), w(t.w) {}
	//template<Arithmetic U> explicit Vector4(const Tuple4<U>& t) noexcept : x(T(t.x)), y(T(t.y)), z(T(t.z)), w(T(t.w)) {}
	explicit Vector4(const TupleType& t) noexcept : x(std::get<0>(t)), y(std::get<1>(t)), z(std::get<2>(t)), w(std::get<3>(t)) {}
	template<Arithmetic U> explicit Vector4(const OtherTupleType<U>& t) noexcept : x(T(std::get<0>(t))), y(T(std::get<1>(t))), z(T(std::get<2>(t))), w(T(std::get<3>(t))) {}
	explicit Vector4(const T* v) noexcept : x(v[0]), y(v[1]), z(v[2]), w(v[3]) {}
	explicit Vector4(Axis axis) noexcept : x((axis == Axis::X) ? T(1) : T(0)), y((axis == Axis::Y) ? T(1) : T(0)), z((axis == Axis::Z) ? T(1) : T(0)), w((axis == Axis::W) ? T(1) : T(0)) {}
	template<Arithmetic U> explicit Vector4(const Vector4<U>& v) noexcept : x(T(t.x)), y(T(t.y)), z(T(t.z)), w(T(t.w)) {}

	//explicit operator Tuple4<T>() noexcept { return Tuple4<T>(x, y, z, w); }
	//template<Arithmetic U> explicit operator Tuple4<U>() noexcept { return Tuple4<U>(U(x), U(y), U(z), U(w)); }
	//explicit operator TupleType() { return TupleType(x, y, z, w); }
	//template<Arithmetic U> explicit operator OtherTupleType<U>() { return OtherTupleType<U>(U(x), U(y), U(z), U(w)); }
	explicit operator T*() noexcept { return &x; }
	explicit operator const T*() const noexcept { return &x; }
	T& operator[](int i) noexcept { return (&x)[i]; }
	const T& operator[](int i) const noexcept { return (&x)[i]; }

	Vector4 operator+() const noexcept { return *this; }
	Vector4 operator-() const noexcept { return Vector4(-x, -y, -z, -w); }
	Vector4& operator+=(const Vector4& v) noexcept { x += v.x; y += v.y; z += v.z; w += v.w; return *this; }
	Vector4& operator-=(const Vector4& v) noexcept { x -= v.x; y -= v.y; z -= v.z; w -= v.w; return *this; }
	Vector4& operator*=(const Vector4& v) noexcept { x *= v.x; y *= v.y; z *= v.z; w *= v.w; return *this; }
	Vector4& operator*=(T f) noexcept { x *= f; y *= f; z *= f; w *= f; return *this; }
	Vector4& operator*=(const Matrix4<T>& m) noexcept;
	Vector4& operator/=(const Vector4& v) noexcept { x /= v.x; y /= v.y; z /= v.z; w /= v.w; return *this; }
	Vector4& operator/=(T f) noexcept { return operator*=(T(1)/f); }
	bool operator==(const Vector4& v) const noexcept { return (x == v.x) && (y == v.y) && (z == v.z) && (w == v.w); }
	bool operator!=(const Vector4& v) const noexcept { return !(*this == v); }
	
	template<typename A> void load(A& ar) { ar(x, y, z, w); }
	template<typename A> void save(A& ar) const { ar(x, y, z, w); }

	template<std::size_t I> T& get() noexcept;
	template<std::size_t I> const T& get() const noexcept;

	const Vector2<T>& xy/*getXY*/() const noexcept { return reinterpret_cast<const Vector2<T>&>(*this); }
	const Vector3<T>& xyz/*getXYZ*/() const noexcept { return reinterpret_cast<const Vector3<T>&>(*this); }
	//void setXY(const Vector2<T>& v) noexcept { x = v.x; y = v.y; }
	//void setXYZ(const Vector3<T>& v) noexcept { x = v.x; y = v.y; z = v.z; }
	bool isZero() const noexcept { return (x == T()) && (y == T()) && (z == T()) && (w == T()); }
	bool isApproxZero() const noexcept;
	bool approxEquals(const Vector4& v) const noexcept;
	bool approxEquals(const Vector4& v, T tolerance) const noexcept;
	bool allLessThan(const Vector4& v) const noexcept { return (x < v.x) && (y < v.y) && (z < v.z) && (w < v.w); }
	bool allLessThanEqual(const Vector4& v) const noexcept { return (x <= v.x) && (y <= v.y) && (z <= v.z) && (w <= v.w); }
	bool allGreaterThan(const Vector4& v) const noexcept { return (x > v.x) && (y > v.y) && (z > v.z) && (w > v.w); }
	bool allGreaterThanEqual(const Vector4& v) const noexcept { return (x >= v.x) && (y >= v.y) && (z >= v.z) && (w >= v.w); }
	bool anyLessThan(const Vector4& v) const noexcept { return (x < v.x) || (y < v.y) || (z < v.z) || (w < v.w); }
	bool anyLessThanEqual(const Vector4& v) const noexcept { return (x <= v.x) || (y <= v.y) || (z <= v.z) || (w <= v.w); }
	bool anyGreaterThan(const Vector4& v) const noexcept { return (x > v.x) || (y > v.y) || (z > v.z) || (w > v.w); }
	bool anyGreaterThanEqual(const Vector4& v) const noexcept { return (x >= v.x) || (y >= v.y) || (z >= v.z) || (w >= v.w); }
	bool isFinite() const noexcept { return std::isfinite(x) && std::isfinite(y) && std::isfinite(z) && std::isfinite(w); }
	T getMagnitude() const noexcept { return std::sqrt(x*x + y*y + z*z + w*w); }
	T getMagnitudeSquared() const noexcept { return (x*x + y*y + z*z + w*w); }
	void setMagnitude(T magnitude) noexcept;
	T getLength() const noexcept { return getMagnitude(); }
	T getLengthSquared() const noexcept { return getMagnitudeSquared(); }
	void setLength(T length) noexcept { setMagnitude(length); }
	T getMinComponent() const { return std::min(std::min(std::min(x, y), z), w); }
	T getMaxComponent() const { return std::max(std::max(std::max(x, y), z), w); }
	Vector4& setZero() noexcept { x = T(); y = T(); z = T(); w = T(); return *this; }
	Vector4& set(T x, T y, T z, T w) noexcept { this->x = x; this->y = y; this->z = z; this->w = w; return *this; }
	Vector4& setMinimum(const Vector4& v1, const Vector4& v2);
	Vector4& setMaximum(const Vector4& v1, const Vector4& v2);
	Vector4& negate() noexcept { x = -x; y = -y; z = -z; w = -w; return *this; }
	Vector4& normalize() noexcept;
	Vector4& transform(const Matrix4<T>& m) noexcept;

	static const Vector4 ZERO;
	static const Vector4 UNIT_X;
	static const Vector4 UNIT_Y;
	static const Vector4 UNIT_Z; 
	static const Vector4 UNIT_W; 
	static const Vector4 MINUS_UNIT_X;
	static const Vector4 MINUS_UNIT_Y;
	static const Vector4 MINUS_UNIT_Z;
	static const Vector4 MINUS_UNIT_W;
	static const Vector4 ONE;
	static const Vector4 TOLERANCE;
	static const Vector4 INF;
	static const Vector4 MINUS_INF;

	T x, y, z, w;
};

template<std::floating_point T> const Vector4<T> Vector4<T>::ZERO{};
template<std::floating_point T> const Vector4<T> Vector4<T>::UNIT_X{ T(1), T(0), T(0), T(0) };
template<std::floating_point T> const Vector4<T> Vector4<T>::UNIT_Y{ T(0), T(1), T(0), T(0) };
template<std::floating_point T> const Vector4<T> Vector4<T>::UNIT_Z{ T(0), T(0), T(1), T(0) };
template<std::floating_point T> const Vector4<T> Vector4<T>::UNIT_W{ T(0), T(0), T(0), T(1) };
template<std::floating_point T> const Vector4<T> Vector4<T>::MINUS_UNIT_X{ T(-1), T(0), T(0), T(0) };
template<std::floating_point T> const Vector4<T> Vector4<T>::MINUS_UNIT_Y{ T(0), T(-1), T(0), T(0) };
template<std::floating_point T> const Vector4<T> Vector4<T>::MINUS_UNIT_Z{ T(0), T(0), T(-1), T(0) };
template<std::floating_point T> const Vector4<T> Vector4<T>::MINUS_UNIT_W{ T(0), T(0), T(0), T(-1) };
template<std::floating_point T> const Vector4<T> Vector4<T>::ONE{ T(1), T(1), T(1), T(1) };
template<std::floating_point T> const Vector4<T> Vector4<T>::TOLERANCE{ Constants<T>::TOLERANCE, Constants<T>::TOLERANCE, Constants<T>::TOLERANCE, Constants<T>::TOLERANCE };
template<std::floating_point T> const Vector4<T> Vector4<T>::INF{ std::numeric_limits<T>::infinity(), std::numeric_limits<T>::infinity(), std::numeric_limits<T>::infinity(), std::numeric_limits<T>::infinity() };
template<std::floating_point T> const Vector4<T> Vector4<T>::MINUS_INF{ -std::numeric_limits<T>::infinity(), -std::numeric_limits<T>::infinity(), -std::numeric_limits<T>::infinity(), -std::numeric_limits<T>::infinity() };

template<std::integral T>
struct Vector4<T>
{
	using ComponentType = T;
	using ConstArg = const Vector4&;
	using ConstResult = const Vector4&;
	using TupleType = std::tuple<T, T, T, T>;
	template<Arithmetic U> OtherTupleType = std::tuple<U, U, U, U>;

	static constexpr int NUM_COMPONENTS = 4;

	constexpr Vector4() noexcept : x(), y(), z(), w() {}
	explicit Vector4(Uninitialized) noexcept {}
	constexpr explicit Vector4(T scalar) noexcept : x(scalar), y(scalar), z(scalar), w(scalar) {}
	constexpr Vector4(T x, T y, T z, T w) noexcept : x(x), y(y), z(z), w(w) {}
	constexpr Vector4(const Vector2<T>& v) noexcept : x(v.x), y(v.y), z(), w(std::numeric_limits<T>::max()) {}
	constexpr Vector4(const Vector2<T>& v, T z, T w) noexcept : x(v.x), y(v.y), z(z), w(w) {}
	constexpr Vector4(const Vector2<T>& xy, const Vector2<T>& zw) noexcept : x(xy.x), y(xy.y), z(zw.x), w(zw.y) {}
	constexpr Vector4(const Vector3<T>& v) noexcept : x(v.x), y(v.y), z(v.z), w(std::numeric_limits<T>::max()) {}
	constexpr Vector4(const Vector3<T>& v, T w) noexcept : x(v.x), y(v.y), z(v.z), w(w) {}
	//explicit Vector4(const Tuple4<T>& t) noexcept : x(t.x), y(t.y), z(t.z), w(t.w) {}
	//template<Arithmetic U> explicit Vector4(const Tuple4<U>& t) noexcept : x(T(t.x)), y(T(t.y)), z(T(t.z)), w(T(t.w)) {}
	explicit Vector4(const TupleType& t) noexcept : x(std::get<0>(t)), y(std::get<1>(t)), z(std::get<2>(t)), w(std::get<3>(t)) {}
	template<Arithmetic U> explicit Vector4(const OtherTupleType<U>& t) noexcept : x(T(std::get<0>(t))), y(T(std::get<1>(t))), z(T(std::get<2>(t))), w(T(std::get<3>(t))) {}
	explicit Vector4(const T* v) noexcept : x(v[0]), y(v[1]), z(v[2]), w(v[3]) {}
	template<Arithmetic U> explicit Vector4(const Vector4<U>& v) noexcept : x(T(v.x)), y(T(v.y)), z(T(v.z)), w(T(v.w)) {}

	//explicit operator Tuple4<T>() noexcept { return Tuple4<T>(x, y, z, w); }
	//template<Arithmetic U> explicit operator Tuple4<U>() noexcept { return Tuple4<U>(U(x), U(y), U(z), U(w)); }
	//explicit operator TupleType() { return TupleType(x, y, z, w); }
	//template<Arithmetic U> explicit operator OtherTupleType<U>() { return OtherTupleType<U>(U(x), U(y), U(z), U(w)); }
	explicit operator T* () noexcept { return &x; }
	explicit operator const T* () const noexcept { return &x; }
	T& operator[](int i) noexcept { return (&x)[i]; }
	const T& operator[](int i) const noexcept { return (&x)[i]; }

	Vector4 operator+() const noexcept { return *this; }
	Vector4 operator-() const noexcept { return Vector4(-x, -y, -z, -w); }
	Vector4& operator+=(const Vector4& v) noexcept { x += v.x; y += v.y; z += v.z; w += v.w; return *this; }
	Vector4& operator-=(const Vector4& v) noexcept { x -= v.x; y -= v.y; z -= v.z; w -= v.w; return *this; }
	Vector4& operator*=(const Vector4& v) noexcept { x *= v.x; y *= v.y; z *= v.z; w *= v.w; return *this; }
	Vector4& operator*=(T f) noexcept { x *= f; y *= f; z *= f; w *= f; return *this; }
	Vector4& operator/=(const Vector4& v) noexcept { x /= v.x; y /= v.y; z /= v.z; w /= v.w; return *this; }
	Vector4& operator/=(T f) noexcept { x /= f; y /= f; z /= f; w /= f; return *this; }
	bool operator==(const Vector4& v) const noexcept { return (x == v.x) && (y == v.y) && (z == v.z) && (w == v.w); }
	bool operator!=(const Vector4& v) const noexcept { return !(*this == v); }

	template<typename A> void serialize(A& ar) { ar(x, y, z, w); }

	template<std::size_t I> T& get() noexcept;
	template<std::size_t I> const T& get() const noexcept;

	const Vector2<T>& xy/*getXY*/() const noexcept { return reinterpret_cast<const Vector2<T>&>(*this); }
	const Vector3<T>& xyz/*getXYZ*/() const noexcept { return reinterpret_cast<const Vector3<T>&>(*this); }
	//void setXY(const Vector2<T>& v) noexcept { x = v.x; y = v.y; }
	//void setXYZ(const Vector3<T>& v) noexcept { x = v.x; y = v.y; z = v.z; }
	bool isZero() const noexcept { return (x == T()) && (y == T()) && (z == T()) && (w == T()); }
	bool allLessThan(const Vector4& v) const noexcept { return (x < v.x) && (y < v.y) && (z < v.z) && (w < v.w); }
	bool allLessThanEqual(const Vector4& v) const noexcept { return (x <= v.x) && (y <= v.y) && (z <= v.z) && (w <= v.w); }
	bool allGreaterThan(const Vector4& v) const noexcept { return (x > v.x) && (y > v.y) && (z > v.z) && (w > v.w); }
	bool allGreaterThanEqual(const Vector4& v) const noexcept { return (x >= v.x) && (y >= v.y) && (z >= v.z) && (w >= v.w); }
	bool anyLessThan(const Vector4& v) const noexcept { return (x < v.x) || (y < v.y) || (z < v.z) || (w < v.w); }
	bool anyLessThanEqual(const Vector4& v) const noexcept { return (x <= v.x) || (y <= v.y) || (z <= v.z) || (w <= v.w); }
	bool anyGreaterThan(const Vector4& v) const noexcept { return (x > v.x) || (y > v.y) || (z > v.z) || (w > v.w); }
	bool anyGreaterThanEqual(const Vector4& v) const noexcept { return (x >= v.x) || (y >= v.y) || (z >= v.z) || (w >= v.w); }
	T getMinComponent() const { return std::min(std::min(std::min(x, y), z), w); }
	T getMaxComponent() const { return std::max(std::max(std::max(x, y), z), w); }
	Vector4& setZero() noexcept { x = T(); y = T(); z = T(); w = T(); return *this; }
	Vector4& set(T x, T y, T z, T w) noexcept { this->x = x; this->y = y; this->z = z; this->w = w; return *this; }
	Vector4& setMinimum(const Vector4& v1, const Vector4& v2);
	Vector4& setMaximum(const Vector4& v1, const Vector4& v2);
	Vector4& negate() noexcept { x = -x; y = -y; z = -z; w = -w; return *this; }

	static const Vector4 ZERO;

	T x, y, z, w;
};

template<std::integral T> const Vector4<T> Vector4<T>::ZERO{};

#if SIMD_HAS_FLOAT4

template<>
struct Matrix4<float>;

template<>
struct Vector4<float>
{
	using Real = float;
	using ComponentType = float;
	using ConstArg = const Vector4;
	using ConstResult = const Vector4;
	using TupleType = std::tuple<float, float, float, float>;
	template<Arithmetic U> OtherTupleType = std::tuple<U, U, U, U>;
	using SimdType = simd::float4;

	static constexpr int NUM_COMPONENTS = 4;

	/*constexpr*/ Vector4() noexcept : xyzw(simd::zero<simd::float4>()) {}
	explicit Vector4(Uninitialized) noexcept {}
	/*constexpr*/ explicit Vector4(float scalar) noexcept : xyzw(simd::set4(scalar)) {}
	/*constexpr*/ Vector4(float x, float y, float z, float w) noexcept : xyzw(simd::set4(x, y, z, w)) {}
	/*constexpr*/ Vector4(const Vector2<float>& v) noexcept : xyzw(simd::insert2(v, UNIT_W)) {}
	/*constexpr*/ Vector4(const Vector2<float>& v, float z, float w) noexcept : xyzw(simd::pack2x2(v, simd::set2(z, w))) {}
	/*constexpr*/ Vector4(const Vector2<float>& xy, const Vector2<float>& zw) noexcept : xyzw(simd::pack2x2(xy, zw)) {}
	/*constexpr*/ Vector4(const Vector3<float>& v) noexcept : xyzw(simd::insert3(v, UNIT_W)) {}
	/*constexpr*/ Vector4(const Vector3<float>& v, float w) noexcept : xyzw(simd::insert<simd::W>(w, v)) {}
	//explicit Vector4(const Tuple4<float>& t) noexcept : xyzw(simd::set4(t.x, t.y, t.z, t.w)) {}
	//template<Arithmetic U> explicit Vector4(const Tuple4<U>& t) noexcept : xyzw(simd::set4((float)t.x, (float)t.y, (float)t.z, (float)t.w)) {}
	explicit Vector4(const TupleType& t) noexcept : xyzw(simd::set4(std::get<0>(t), std::get<1>(t), std::get<2>(t), std::get<3>(t))) {}
	template<Arithmetic U> explicit Vector4(const OtherTupleType<U>& t) noexcept : xyzw(simd::set4((float)std::get<0>(t), (float)std::get<1>(t), (float)std::get<2>(t), (float)std::get<3>(t))) {}
	explicit Vector4(const float* v) noexcept : xyzw(simd::load4(v)) {}
	explicit Vector4(Axis axis) noexcept : Vector4((axis == Axis::X) ? 1.f : 0.f, (axis == Axis::Y) ? 1.f : 0.f, (axis == Axis::Z) ? 1.f : 0.f, (axis == Axis::W) ? 1.f : 0.f)) {}
	explicit Vector4(simd::float4 v) noexcept : xyzw(v) {}
	Vector4(const Vector4& v) noexcept : xyzw(v.xyzw) {}
	template<Arithmetic U> explicit Vector4(const Vector4<U>& v) noexcept : xyzw(simd::set4((float)v.x, (float)v.y, (float)v.z, (float)v.w)) {}
	Vector4& operator=(const Vector4& v) noexcept { xyzw = v.xyzw; return *this; }

	operator simd::float4() const noexcept { return xyzw; }
	//explicit operator Tuple4<float>() noexcept { return Tuple4<float>(x, y, z, w); }
	//template<Arithmetic U> explicit operator Tuple4<U>() noexcept { return Tuple4<U>(U(x), U(y), U(z), U(w)); }
	//explicit operator TupleType() { return TupleType(x, y, z, w); }
	//template<Arithmetic U> explicit operator OtherTupleType<U>() { return OtherTupleType<U>(U(x), U(y), U(z), U(w)); }
	explicit operator float* () noexcept { return &x; }
	explicit operator const float* () const noexcept { return &x; }
	float& operator[](int i) noexcept { return (&x)[i]; }
	const float& operator[](int i) const noexcept { return (&x)[i]; }

	Vector4 operator+() const noexcept { return *this; }
	Vector4 operator-() const noexcept { return Vector4(simd::neg4(xyzw)); }
	Vector4& operator+=(const Vector4& v) noexcept { xyzw = simd::add4(xyzw, v); return *this; }
	Vector4& operator-=(const Vector4& v) noexcept { xyzw = simd::sub4(xyzw, v); return *this; }
	Vector4& operator*=(const Vector4& v) noexcept { xyzw = simd::mul4(xyzw, v); return *this; }
	Vector4& operator*=(float f) noexcept { xyzw = simd::mul4(xyzw, simd::set4(f)); return *this; }
	Vector4& operator*=(const Matrix4<float>& m) noexcept;
	Vector4& operator/=(const Vector4& v) noexcept { xyzw = simd::div4(xyzw, v); return *this; }
	Vector4& operator/=(float f) noexcept { xyzw = simd::div4(xyzw, simd::set4(f)); return *this; }
	bool operator==(const Vector4& v) const noexcept { return simd::all4(simd::equal(xyzw, v)); }
	bool operator!=(const Vector4& v) const noexcept { return !(*this == v); }

	template<typename A> void load(A& ar);
	template<typename A> void save(A& ar) const { ar(x, y, z, w); }

	template<std::size_t I> float& get() noexcept;
	template<std::size_t I> const float& get() const noexcept;
	template<typename U> U& get() noexcept;				// intentionally undefined
	template<typename U> const U& get() const noexcept;	// intentionally undefined
	template<> simd::float4& get() noexcept { return xyzw; }
	template<> const simd::float4& get() const noexcept { return xyzw; }

#if MATHEMATICS_SIMD_EXPAND_LAST
	const Vector2<float> xy/*getXY*/() const noexcept { return Vector2<float>(simd::xyyy(xyzw)); }
	const Vector3<float> xyz/*getXYZ*/() const noexcept { return Vector3<float>(simd::xyzz(xyzw)); }
#else
	const Vector2<float> xy/*getXY*/() const noexcept { return Vector2<float>(simd::cutoff2(xyzw)); }
	const Vector3<float> xyz/*getXYZ*/() const noexcept { return Vector3<float>(simd::cutoff3(xyzw)); }
#endif
	//void setXY(const Vector2<float>& v) noexcept { xyzw = simd::insert2(v, xyzw); }
	//void setXYZ(const Vector3<float>& v) noexcept { xyzw = simd::insert3(v, xyzw); }
	bool isZero() const noexcept { return simd::all4(simd::equal(xyzw, simd::zero<simd::float4>())); }
	bool isApproxZero() const noexcept { return simd::all4(simd::lessThan(simd::abs4(xyzw), TOLERANCE)); }
	bool approxEquals(const Vector4& v) const noexcept { return simd::all4(simd::lessThan(simd::abs4(simd::sub4(xyzw, v)), TOLERANCE)); }
	bool approxEquals(const Vector4& v, float tolerance) const noexcept { return simd::all4(simd::lessThan(simd::abs4(simd::sub4(xyzw, v)), simd::set4(tolerance))); }
	bool allLessThan(const Vector4& v) const noexcept { return simd::all4(simd::lessThan(xyzw, v)); }
	bool allLessThanEqual(const Vector4& v) const noexcept { return simd::all4(simd::lessThanEqual(xyzw, v)); }
	bool allGreaterThan(const Vector4& v) const noexcept { return simd::all4(simd::greaterThan(xyzw, v)); }
	bool allGreaterThanEqual(const Vector4& v) const noexcept { return simd::all4(simd::greaterThanEqual(xyzw, v)); }
	bool anyLessThan(const Vector4& v) const noexcept { return simd::any4(simd::lessThan(xyzw, v)); }
	bool anyLessThanEqual(const Vector4& v) const noexcept { return simd::any4(simd::lessThanEqual(xyzw, v)); }
	bool anyGreaterThan(const Vector4& v) const noexcept { return simd::any4(simd::greaterThan(xyzw, v)); }
	bool anyGreaterThanEqual(const Vector4& v) const noexcept { return simd::any4(simd::greaterThanEqual(xyzw, v)); }
	bool isFinite() const noexcept { return simd::all4(simd::isFinite(xyzw)); }
	float getMagnitude() const noexcept { return simd::toFloat(simd::sqrt1(simd::dot4(xyzw, xyzw))); }
	float getMagnitudeSquared() const noexcept { return simd::toFloat(simd::dot4(xyzw, xyzw)); }
	void setMagnitude(float magnitude) noexcept;
	float getLength() const noexcept { return getMagnitude(); }
	float getLengthSquared() const noexcept { return getMagnitudeSquared(); }
	void setLength(float length) noexcept { setMagnitude(length); }
	float getMinComponent() const noexcept { return simd::toFloat(simd::hMin4(xyzw)); }
	float getMaxComponent() const noexcept { return simd::toFloat(simd::hMax4(xyzw)); }
	Vector4& setZero() noexcept { xyzw = simd::zero<simd::float4>(); return *this; }
	Vector4& set(float x, float y, float z, float w) noexcept { xyzw = simd::set4(x, y, z, w); return *this; }
	Vector4& setMinimum(const Vector4& v1, const Vector4& v2) noexcept { xyzw = simd::min4(v1, v2); return *this; }
	Vector4& setMaximum(const Vector4& v1, const Vector4& v2) noexcept { xyzw = simd::max4(v1, v2); return *this; }
	Vector4& negate() noexcept { xyzw = simd::neg4(xyzw); return *this; }
	Vector4& normalize() noexcept;
	Vector4& transform(const Matrix4<float>& m) noexcept;

	static const Vector4 ZERO;
	static const Vector4 UNIT_X;
	static const Vector4 UNIT_Y;
	static const Vector4 UNIT_Z;
	static const Vector4 UNIT_W;
	static const Vector4 MINUS_UNIT_X;
	static const Vector4 MINUS_UNIT_Y;
	static const Vector4 MINUS_UNIT_Z;
	static const Vector4 MINUS_UNIT_W;
	static const Vector4 ONE;
	static const Vector4 TOLERANCE;
	static const Vector4 INF;
	static const Vector4 MINUS_INF;

	union
	{
		simd::float4 xyzw;
		struct { float x, y, z, w; };
	};
};

const Vector4<float> Vector4<float>::ZERO{};
const Vector4<float> Vector4<float>::UNIT_X{ 1.f, 0.f, 0.f, 0.f };
const Vector4<float> Vector4<float>::UNIT_Y{ 0.f, 1.f, 0.f, 0.f };
const Vector4<float> Vector4<float>::UNIT_Z{ 0.f, 0.f, 1.f, 0.f };
const Vector4<float> Vector4<float>::UNIT_W{ 0.f, 0.f, 0.f, 1.f };
const Vector4<float> Vector4<float>::MINUS_UNIT_X{ -1.f, 0.f, 0.f, 0.f };
const Vector4<float> Vector4<float>::MINUS_UNIT_Y{ 0.f, -1.f, 0.f, 0.f };
const Vector4<float> Vector4<float>::MINUS_UNIT_Z{ 0.f, 0.f, -1.f, 0.f };
const Vector4<float> Vector4<float>::MINUS_UNIT_W{ 0.f, 0.f, 0.f, -1.f };
const Vector4<float> Vector4<float>::ONE{ 1.f, 1.f, 1.f, 1.f };
const Vector4<float> Vector4<float>::TOLERANCE{ Constants<float>::TOLERANCE, Constants<float>::TOLERANCE, Constants<float>::TOLERANCE, Constants<float>::TOLERANCE };
const Vector4<float> Vector4<float>::INF{ std::numeric_limits<float>::infinity(), std::numeric_limits<float>::infinity(), std::numeric_limits<float>::infinity(), std::numeric_limits<float>::infinity() };
const Vector4<float> Vector4<float>::MINUS_INF{ -std::numeric_limits<float>::infinity(), -std::numeric_limits<float>::infinity(), -std::numeric_limits<float>::infinity(), -std::numeric_limits<float>::infinity() };

#endif /* SIMD_HAS_FLOAT4 */

template<typename T>
	requires (std::floating_point<T> || std::integral<T>)
inline Vector4<T> operator+(const Vector4<T>& v1, const Vector4<T>& v2) noexcept
{ 
	return Vector4<T>(v1.x + v2.x, v1.y + v2.y, v1.z + v2.z, v1.w + v2.w); 
}

template<typename T>
	requires (std::floating_point<T> || std::integral<T>)
inline Vector4<T> operator-(const Vector4<T>& v1, const Vector4<T>& v2) noexcept
{ 
	return Vector4<T>(v1.x - v2.x, v1.y - v2.y, v1.z - v2.z, v1.w - v2.w); 
}

template<typename T>
	requires (std::floating_point<T> || std::integral<T>)
inline Vector4<T> operator*(const Vector4<T>& v1, const Vector4<T>& v2) noexcept
{ 
	return Vector4<T>(v1.x*v2.x, v1.y*v2.y, v1.z*v2.z, v1.w*v2.w); 
}

template<typename T>
	requires (std::floating_point<T> || std::integral<T>)
inline Vector4<T> operator*(T f, const Vector4<T>& v) noexcept
{ 
	return Vector4<T>(f*v.x, f*v.y, f*v.z, f*v.w); 
}

template<typename T>
	requires (std::floating_point<T> || std::integral<T>)
inline Vector4<T> operator*(const Vector4<T>& v, T f) noexcept
{ 
	return Vector4<T>(v.x*f, v.y*f, v.z*f, v.w*f); 
}

template<typename T>
	requires (std::floating_point<T> || std::integral<T>)
inline Vector4<T> operator/(const Vector4<T>& v1, const Vector4<T>& v2) noexcept
{ 
	return Vector4<T>(v1.x/v2.x, v1.y/v2.y, v1.z/v2.z, v1.w/v2.w); 
}

template<typename T>
	requires (std::floating_point<T> || std::integral<T>)
inline Vector4<T> operator/(T f, const Vector4<T>& v) noexcept
{ 
	return Vector4<T>(f/v.x, f/v.y, f/v.z, f/v.w); 
}

template<typename T>
	requires (std::floating_point<T> || std::integral<T>)
inline Vector4<T> operator/(const Vector4<T>& v, T f) noexcept
{ 
	if constexpr (std::is_floating_point_v<T>)
		return operator*(v, T(1)/f);
	else
		return Vector4<T>(v.x/f, v.y/f, v.z/f, v.w/f);
}

template<typename C, typename T, typename U>
	requires std::floating_point<U> || std::integral<U>
inline std::basic_istream<C, T>& operator>>(std::basic_istream<C, T>& s, Vector4<U>& v)
{ 
	return s >> v.x >> std::ws >> v.y >> std::ws >> v.z >> std::ws >> v.w; 
}

template<typename C, typename T, typename U>
	requires std::floating_point<U> || std::integral<U>
inline std::basic_ostream<C, T>& operator<<(std::basic_ostream<C, T>& s, const Vector4<U>& v)
{ 
	constexpr C WS(0x20);
	return s << v.x << WS << v.y << WS << v.z << WS << v.w;
}

template<std::floating_point T>
template<std::size_t I>
inline T& Vector4<T>::get()
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

template<std::floating_point T>
template<std::size_t I>
inline const T& Vector4<T>::get() const
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

template<std::floating_point T>
inline bool Vector4<T>::isApproxZero() const
{ 
	return (std::fabs(x) < Constants<T>::TOLERANCE) && (std::fabs(y) < Constants<T>::TOLERANCE) && 
		(std::fabs(z) < Constants<T>::TOLERANCE) && (std::fabs(w) < Constants<T>::TOLERANCE);
}

template<std::floating_point T>
inline bool Vector4<T>::approxEquals(const Vector4<T>& v) const
{ 
	return (std::fabs(v.x - x) < Constants<T>::TOLERANCE) && (std::fabs(v.y - y) < Constants<T>::TOLERANCE) && 
		(std::fabs(v.z - z) < Constants<T>::TOLERANCE) && (std::fabs(v.w - w) < Constants<T>::TOLERANCE); 
}

template<std::floating_point T>
inline bool Vector4<T>::approxEquals(const Vector4<T>& v, T tolerance) const
{
	return (std::fabs(v.x - x) < tolerance) && (std::fabs(v.y - y) < tolerance) && 
		(std::fabs(v.z - z) < tolerance) && (std::fabs(v.w - w) < tolerance); 
}

template<std::floating_point T>
inline void Vector4<T>::setMagnitude(T magnitude) 
{ 
	T m = getMagnitude(); 
	if (m > T(0)) 
		*this *= magnitude/m;
}

template<std::floating_point T>
inline Vector4<T>& Vector4<T>::setMinimum(const Vector4<T>& v1, const Vector4<T>& v2)
{ 
	x = std::min(v1.x, v2.x); 
	y = std::min(v1.y, v2.y);
	z = std::min(v1.z, v2.z); 
	w = std::min(v1.w, v2.w); 
	return *this; 
}

template<std::floating_point T>
inline Vector4<T>& Vector4<T>::setMaximum(const Vector4<T>& v1, const Vector4<T>& v2)
{ 
	x = std::max(v1.x, v2.x); 
	y = std::max(v1.y, v2.y);
	z = std::max(v1.z, v2.z); 
	w = std::max(v1.w, v2.w); 
	return *this; 
}

template<std::floating_point T>
inline Vector4<T>& Vector4<T>::normalize()
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

template<std::integral T>
template<std::size_t I>
inline T& Vector4<T>::get()
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

template<std::integral T>
template<std::size_t I>
inline const T& Vector4<T>::get() const
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

template<std::integral T>
inline Vector4<T>& Vector4<T>::setMinimum(const Vector4<T>& v1, const Vector4<T>& v2)
{
	x = std::min(v1.x, v2.x);
	y = std::min(v1.y, v2.y);
	z = std::min(v1.z, v2.z);
	w = std::min(v1.w, v2.w);
	return *this;
}

template<std::integral T>
inline Vector4<T>& Vector4<T>::setMaximum(const Vector4<T>& v1, const Vector4<T>& v2)
{
	x = std::max(v1.x, v2.x);
	y = std::max(v1.y, v2.y);
	z = std::max(v1.z, v2.z);
	w = std::max(v1.w, v2.w);
	return *this;
}

#if SIMD_HAS_FLOAT4

template<>
inline Vector4<float> operator+(const Vector4<float>& v1, const Vector4<float>& v2) noexcept 
{ 
	return Vector4<float>(simd::add4(v1, v2)); 
}

template<>
inline Vector4<float> operator-(const Vector4<float>& v1, const Vector4<float>& v2) noexcept 
{ 
	return Vector4<float>(simd::sub4(v1, v2));
}

template<>
inline Vector4<float> operator*(const Vector4<float>& v1, const Vector4<float>& v2) noexcept 
{ 
	return Vector4<float>(simd::mul4(v1, v2)); 
}

template<>
inline Vector4<float> operator*(float f, const Vector4<float>& v) noexcept 
{ 
	return Vector4<float>(simd::mul4(simd::set4(f), v)); 
}

template<>
inline Vector4<float> operator*(const Vector4<float>& v, float f) noexcept 
{ 
	return Vector4<float>(simd::mul4(v, simd::set4(f))); 
}

template<>
inline Vector4<float> operator/(const Vector4<float>& v1, const Vector4<float>& v2) noexcept 
{ 
	return Vector4<float>(simd::div4(v1, v2)); 
}

template<>
inline Vector4<float> operator/(float f, const Vector4<float>& v) noexcept 
{ 
	return Vector4<float>(simd::div4(simd::set4(f), v)); 
}

template<>
inline Vector4<float> operator/(const Vector4<float>& v, float f) noexcept 
{ 
	return Vector4<float>(simd::div4(v, simd::set4(f))); 
}

template<typename C, typename T>
inline std::basic_istream<C, T>& operator>>(std::basic_istream<C, T>& s, Vector4<float>& v)
{ 
	float x, y, z, w; 
	s >> x >> std::ws >> y >> std::ws >> z >> std::ws >> w; 
	v.set(x, y, z, w); 
	return s; 
}

template<typename A>
inline void Vector4<float>::load(A& ar)
{
	float t0, t1, t2, t3;
	ar(t0, t1, t2, t3);
	set(t0, t1, t2, t3);
}

template<std::size_t I>
inline float& Vector4<float>::get()
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
inline const float& Vector4<float>::get() const
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

inline void Vector4<float>::setMagnitude(float magnitude)
{ 
	float m = getMagnitude();
	if (m > 0.f) 
		*this *= magnitude/m;
}

inline Vector4<float>& Vector4<float>::normalize()
{
#if MATHEMATICS_FAST_NORMALIZE
	float m = simd::toFloat(simd::rcpSqrtApprox1(simd::dot4(xyzw, xyzw)));
	if (m <= std::numeric_limits<float>::max())
		xyzw = simd::mul4(xyzw, simd::set4(m));
#else
	float m = getMagnitude(); 
	if (m > 0.f) 
		xyzw = simd::div4(xyzw, simd::set4(m));
#endif
	return *this;
}

#endif /* SIMD_HAS_FLOAT4 */

template<std::size_t I, typename T>
	requires (std::floating_point<T> || std::integral<T>)
inline T& get(Vector4<T>& v) noexcept
{
	if constexpr (I == 0)
		return v.x;
	else if constexpr (I == 1)
		return v.y;
	else if constexpr (I == 2)
		return v.z;
	else if constexpr (I == 3)
		return v.w;
	static_assert(false);
}

template<std::size_t I, typename T>
	requires (std::floating_point<T> || std::integral<T>)
inline const T& get(const Vector4<T>& v) noexcept
{
	if constexpr (I == 0)
		return v.x;
	else if constexpr (I == 1)
		return v.y;
	else if constexpr (I == 2)
		return v.z;
	else if constexpr (I == 3)
		return v.w;
	static_assert(false);
}

template<std::size_t I, typename T>
	requires (std::floating_point<T> || std::integral<T>)
inline T&& get(Vector4<T>&& v) noexcept
{
	if constexpr (I == 0)
		return v.x;
	else if constexpr (I == 1)
		return v.y;
	else if constexpr (I == 2)
		return v.z;
	else if constexpr (I == 3)
		return v.w;
	static_assert(false);
}

template<std::size_t I, typename T>
	requires (std::floating_point<T> || std::integral<T>)
inline const T&& get(const Vector4<T>&& v) noexcept
{
	if constexpr (I == 0)
		return v.x;
	else if constexpr (I == 1)
		return v.y;
	else if constexpr (I == 2)
		return v.z;
	else if constexpr (I == 3)
		return v.w;
	static_assert(false);
}

template<typename T>
	requires (std::floating_point<T> || std::integral<T>)
inline Vector4<T> abs(const Vector4<T>& v) noexcept
{
	if constexpr (std::is_floating_point_v<T>)
		return Vector4<T>(std::fabs(v.x), std::fabs(v.y), std::fabs(v.z), std::fabs(v.w));
	else // integral
		return Vector4<T>(abs(v.x), abs(v.y), abs(v.z), abs(v.w));
}

template<typename T>
#if SIMD_HAS_FLOAT4
	requires ((std::floating_point<T> || std::integral<T>) && !std::same_as<T, float>)
#else
	requires (std::floating_point<T> || std::integral<T>)
#endif
inline Vector4<T> abs(Vector4<T>&& v) noexcept
{
	if constexpr (std::is_floating_point_v<T>)
	{
		v.x = std::fabs(v.x);
		v.y = std::fabs(v.y);
		v.z = std::fabs(v.z);
		v.w = std::fabs(v.w);
	}
	else // integral
	{
		v.x = abs(v.x);
		v.y = abs(v.y);
		v.z = abs(v.z);
		v.w = abs(v.w);
	}
	return v;
}

template<typename T>
	requires (std::floating_point<T> || std::integral<T>)
inline T sum(const Vector4<T>& v) noexcept
{
	return v.x + v.y + v.z + v.w;
}

template<typename T>
	requires std::floating_point<T>
inline T dot(const Vector4<T>& v1, const Vector4<T>& v2) noexcept
{
	return v1.x*v2.x + v1.y*v2.y + v1.z*v2.z + v1.w*v2.w;
}

template<typename T>
	requires std::floating_point<T>
inline Vector4<T> cross(const Vector4<T>& v1, const Vector4<T>& v2, const Vector4<T>& v3) noexcept
{
	T d1 = v2.z*v3.w - v2.w*v3.z;
	T d2 = v2.y*v3.w - v2.w*v3.y;
	T d3 = v2.y*v3.z - v2.z*v3.y;
	T d4 = v2.x*v3.w - v2.w*v3.x;
	T d5 = v2.x*v3.z - v2.z*v3.x;
	T d6 = v2.x*v3.y - v2.y*v3.x;
	return Vector4<T>(-v1.y*d1 + v1.z*d2 - v1.w*d3, v1.x*d1 - v1.z*d4 + v1.w*d5, -v1.x*d2 + v1.y*d4 - v1.w*d6,
		v1.x*d3 - v1.y*d5 + v1.z*d6);
}

template<typename T>
	requires std::floating_point<T>
inline T distance(const Vector4<T>& v1, const Vector4<T>& v2) noexcept
{
	T x = v2.x - v1.x;
	T y = v2.y - v1.y;
	T z = v2.z - v1.z;
	T w = v2.w - v1.w;
	return std::sqrt(x*x + y*y + z*z + w*w);
}

template<typename T>
	requires std::floating_point<T>
inline T distanceSquared(const Vector4<T>& v1, const Vector4<T>& v2) noexcept
{
	T x = v2.x - v1.x;
	T y = v2.y - v1.y;
	T z = v2.z - v1.z;
	T w = v2.w - v1.w;
	return x*x + y*y + z*z + w*w;
}

template<typename T>
	requires std::floating_point<T>
inline T length(const Vector4<T>& v) noexcept
{
	return std::sqrt(v.x*v.x + v.y*v.y + v.z*v.z + v.w*v.w);
}

template<typename T>
	requires std::floating_point<T>
inline T lengthSquared(const Vector4<T>& v) noexcept
{
	return v.x*v.x + v.y*v.y + v.z*v.z + v.w*v.w;
}

template<typename T>
	requires (std::floating_point<T> || std::integral<T>)
inline Vector4<T> min(const Vector4<T>& v1, const Vector4<T>& v2)
{
	return Vector4<T>(std::min(v1.x, v2.x), std::min(v1.y, v2.y), std::min(v1.z, v2.z), std::min(v1.w, v2.w));
}

template<typename T>
	requires (std::floating_point<T> || std::integral<T>)
inline Vector4<T> max(const Vector4<T>& v1, const Vector4<T>& v2)
{
	return Vector4<T>(std::max(v1.x, v2.x), std::max(v1.y, v2.y), std::max(v1.z, v2.z), std::max(v1.w, v2.w));
}

template<typename T>
	requires (std::floating_point<T> || std::integral<T>)
inline Vector4<T> clamp(const Vector4<T>& v, const Vector4<T>& low, const Vector4<T>& high)
{
	return Vector4<T>(std::clamp(v.x, low.x, high.x), std::clamp(v.y, low.y, high.y), std::clamp(v.z, low.z, high.z), 
		std::clamp(v.w, low.w, high.w));
}

template<typename T>
	requires std::floating_point<T>
inline Vector4<T> lerp(const Vector4<T>& v1, const Vector4<T>& v2, T t) noexcept
{
	return Vector4<T>(v1.x + t*(v2.x - v1.x), v1.y + t*(v2.y - v1.y), v1.z + t*(v2.z - v1.z), v1.w + t*(v2.w - v1.w));
}

template<typename T>
	requires std::floating_point<T>
inline Vector4<T> slerp(const Vector4<T>& v1, const Vector4<T>& v2, T t)
{
	T dp = v1.x*v2.x + v1.y*v2.y + v1.z*v2.z + v1.w*v2.w;
	if ((T(1) - dp) < Constants<T>::TOLERANCE)
	{
		Vector4<T> c(v1.x + t*(v2.x - v1.x), v1.y + t*(v2.y - v1.y), v1.z + t*(v2.z - v1.z), v1.w + t*(v2.w - v1.w));
		c.normalize();
		return c;
	}

	dp = std::clamp(dp, T(-1), T(1));
	T theta = std::acos(dp)*t;
	T st = std::sin(theta);
	T ct = std::cos(theta);
	Vector4<T> c(v2.x - v1.x*dp, v2.y - v1.y*dp, v2.z - v1.z*dp, v2.w - v1.w*dp);
	c.normalize();
	return Vector4<T>(v1.x*ct + c.x*st, v1.y*ct + c.y*st, v1.z*ct + c.z*st, v1.w*ct + c.w*st);
}

#if SIMD_HAS_FLOAT4

template<>
inline Vector4<float> abs(const Vector4<float>& v) noexcept
{
	return Vector4<float>(simd::abs4(v));
}

template<>
inline T sum(const Vector4<float>& v) noexcept
{
	return simd::toFloat(simd::hAdd4(v));
}

template<>
inline float dot(const Vector4<float>& v1, const Vector4<float>& v2) noexcept
{
	return simd::toFloat(simd::dot4(v1, v2));
}

template<>
inline Vector4<float> cross(const Vector4<float>& v1, const Vector4<float>& v2, const Vector4<float>& v3) noexcept // #TODO SIMD
{
	float d1 = v2.z*v3.w - v2.w*v3.z;
	float d2 = v2.y*v3.w - v2.w*v3.y;
	float d3 = v2.y*v3.z - v2.z*v3.y;
	float d4 = v2.x*v3.w - v2.w*v3.x;
	float d5 = v2.x*v3.z - v2.z*v3.x;
	float d6 = v2.x*v3.y - v2.y*v3.x;
	return Vector4<float>(-v1.y*d1 + v1.z*d2 - v1.w*d3, v1.x*d1 - v1.z*d4 + v1.w*d5, -v1.x*d2 + v1.y*d4 - v1.w*d6,
		v1.x*d3 - v1.y*d5 + v1.z*d6);
}

template<>
inline float distance(const Vector4<float>& v1, const Vector4<float>& v2) noexcept
{
	auto v = simd::sub4(v2, v1);
	return simd::toFloat(simd::sqrt1(simd::dot4(v, v)));
}

template<>
inline float distanceSquared(const Vector4<float>& v1, const Vector4<float>& v2) noexcept
{
	auto v = simd::sub4(v2, v1);
	return simd::toFloat(simd::dot4(v, v));
}

template<>
inline float length(const Vector4<float>& v) noexcept
{
	return simd::toFloat(simd::sqrt1(simd::dot4(v, v)));
}

template<>
inline float lengthSquared(const Vector4<float>& v) noexcept
{
	return simd::toFloat(simd::dot4(v, v));
}

template<>
inline Vector4<float> min(const Vector4<float>& v1, const Vector4<float>& v2)
{
	return Vector4<float>(simd::min4(v1, v2));
}

template<>
inline Vector4<float> max(const Vector4<float>& v1, const Vector4<float>& v2)
{
	return Vector4<float>(simd::max4(v1, v2));
}

template<>
inline Vector4<float> clamp(const Vector4<float>& v, const Vector4<float>& low, const Vector4<float>& high)
{
	return Vector4<float>(simd::min4(simd::max4(v, low), high));
}

template<>
inline Vector4<float> lerp(const Vector4<float>& v1, const Vector4<float>& v2, float t) noexcept
{
	return Vector4<float>(simd::mulAdd4(simd::set4(t), simd::sub4(v2, v1), v1));
}

template<>
inline Vector4<float> slerp(const Vector4<float>& v1, const Vector4<float>& v2, float t)
{
	float dp = simd::toFloat(simd::dot4(v1, v2));
	if ((1.f - dp) < Constants<float>::TOLERANCE)
	{
		Vector4<float> c(simd::mulAdd4(simd::set4(t), simd::sub4(v2, v1), v1));
		c.normalize();
		return c;
	}

	dp = std::clamp(dp, -1.f, 1.f);
	float theta = std::acos(dp)*t;
	float st = std::sin(theta);
	float ct = std::cos(theta);
	Vector4<float> c(simd::sub4(v2, simd::mul4(v1, simd::set4(dp))));
	c.normalize();
	return Vector4<float>(simd::mulAdd4(v1, simd::set4(ct), simd::mul4(c, simd::set4(st))));
}

#endif /* SIMD_HAS_FLOAT4 */

template<typename T>
	requires std::floating_point<T>
inline Vector4<T> normalize(const Vector4<T>& v) noexcept
{
	Vector4<T> u(v);
	u.normalize();
	return u;
}

template<typename T>
#if SIMD_HAS_FLOAT4
	requires (std::floating_point<T> && !std::same_as<T, float>)
#else
	requires std::floating_point<T>
#endif
inline Vector4<T> normalize(Vector4<T>&& v) noexcept
{
	v.normalize();
	return v;
}

//template<typename T>
//	requires std::floating_point<T>
//inline Vector4<T> project(const Vector4<T>& v1, const Vector4<T>& v2) noexcept
//{
//	T m = v2.getMagnitudeSquared();
//	return (m > T(0)) ? (dot(v1, v2)/m)*v2 : Vector4<T>::ZERO;
//}

//template<typename T>
//	requires std::floating_point<T>
//inline T angle(const Vector4<T>& v1, const Vector4<T>& v2)
//{
//#if MATHEMATICS_FAST_NORMALIZE
//	Vector4<T> u(v1);
//	Vector4<T> v(v2);
//	u.normalize();
//	v.normalize();
//	return std::acos(std::clamp(dot(u, v), T(-1), T(1)));
//#else
//	T q = v1.getMagnitude()*v2.getMagnitude();
//	if (q > T(0))
//		return std::acos(std::clamp(dot(v1, v2)/q, T(-1), T(1)));
//	else
//		return (dot(v1, v2) >= T(0)) ? T(0) : Constants<T>::PI;
//#endif
//}

} // namespace templates

using Vector4f = templates::Vector4<float>;
using Vector4d = templates::Vector4<double>;

#if MATHEMATICS_DOUBLE
using Vector4 = templates::Vector4<double>;
using Vector4Arg = templates::Vector4<double>::ConstArg;
using Vector4Result = templates::Vector4<double>::ConstResult;
#else
using Vector4 = templates::Vector4<float>;
using Vector4Arg = templates::Vector4<float>::ConstArg;
using Vector4Result = templates::Vector4<float>::ConstResult;
#endif

using IntVector4 = templates::Vector4<int>;
using IntVector4Arg = templates::Vector4<int>::ConstArg;
using IntVector4Result = templates::Vector4<int>::ConstResult;

} // namespace mathematics

namespace std {

template<size_t I, typename T>
struct tuple_element;

template<size_t I, typename T>
struct tuple_element<I, ::mathematics::templates::Vector4<T>>
{
	using type = T;
};

template<typename T>
struct tuple_size;

template<typename T>
struct tuple_size<::mathematics::templates::Vector4<T>> : integral_constant<size_t, 4> 
{
};

template<typename T>
struct hash;

template<typename T>
struct hash<::mathematics::templates::Vector4<T>>
{
	size_t operator()(const ::mathematics::templates::Vector4<T>& v) const noexcept
	{
		hash<T> hasher;
		size_t seed = hasher(v.x) + 0x9e3779b9;
		seed ^= hasher(v.y) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
		seed ^= hasher(v.z) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
		seed ^= hasher(v.w) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
		return seed;
	}
};

#if SIMD_HAS_FLOAT4

template<>
struct hash<::mathematics::templates::Vector4<float>>
{
	size_t operator()(const ::mathematics::templates::Vector4<float>& v) const noexcept
	{
		return hash<typename ::simd::float4>()(v.xyzw);
	}
};

#endif /* SIMD_HAS_FLOAT4 */

} // namespace std

#include "Matrix4.hpp"

namespace mathematics::templates {

template<typename T>
inline Vector4<T>& Vector4<T>::operator*=(const Matrix4<T>& m)
{
	set(x*m.m00 + y*m.m10 + z*m.m20 + w*m.m30, x*m.m01 + y*m.m11 + z*m.m21 + w*m.m31, 
		x*m.m02 + y*m.m12 + z*m.m22 + w*m.m32, x*m.m03 + y*m.m13 + z*m.m23 + w*m.m33);
	return *this;
}

template<typename T>
	requires std::floating_point<T>
inline Vector4<T> operator*(const Vector4<T>& v, const Matrix4<T>& m) noexcept
{
	return Vector4<T>(v.x*m.m00 + v.y*m.m10 + v.z*m.m20 + v.w*m.m30, v.x*m.m01 + v.y*m.m11 + v.z*m.m21 + v.w*m.m31,
		v.x*m.m02 + v.y*m.m12 + v.z*m.m22 + v.w*m.m32, v.x*m.m03 + v.y*m.m13 + v.z*m.m23 + v.w*m.m33);
}

template<typename T>
	requires std::floating_point<T>
inline Vector4<T> operator*(const Matrix4<T>& m, const Vector4<T>& v) noexcept
{
	return Vector4<T>(m.m00*v.x + m.m01*v.y + m.m02*v.z + m.m03*v.w, m.m10*v.x + m.m11*v.y + m.m12*v.z + m.m13*v.w,
		m.m20*v.x + m.m21*v.y + m.m22*v.z + m.m23*v.w, m.m30*v.x + m.m31*v.y + m.m32*v.z + m.m33*v.w);
}

template<typename T>
inline Vector4& Vector4<T>::transform(const Matrix4<T>& m)
{
	*this *= m;
	return *this;
}

template<typename T>
	requires std::floating_point<T>
inline Vector4<T> transform(const Vector4<T>& v, const Matrix4<T>& m) noexcept
{
	return v*m;
}

template<typename T>
	requires std::floating_point<T>
inline Matrix4<T> tensor(const Vector4<T>& v1, const Vector4<T>& v2) noexcept
{
	return Matrix4<T>(v1.x*v2.x, v1.x*v2.y, v1.x*v2.z, v1.x*v2.w,
		v1.y*v2.x, v1.y*v2.y, v1.y*v2.z, v1.y*v2.w,
		v1.z*v2.x, v1.z*v2.y, v1.z*v2.z, v1.z*v2.w,
		v1.w*v2.x, v1.w*v2.y, v1.w*v2.z, v1.w*v2.w);
}

#if SIMD_HAS_FLOAT4

inline Vector4<float>& Vector4<float>::operator*=(const Matrix4<float>& m)
{
	auto t = simd::mulAdd4(simd::xxxx(xyzw), m.row0, simd::mul4(simd::yyyy(xyzw), m.row1));
	t = simd::add4(t, simd::mul4(simd::zzzz(xyzw), m.row2));
	xyzw = simd::add4(t, simd::mul4(simd::wwww(xyzw), m.row3));
	return *this;
}

template<>
inline Vector4<float> operator*(const Vector4<float>& v, const Matrix4<float>& m) noexcept
{
	auto t = simd::mulAdd4(simd::xxxx(v), m.row0, simd::mul4(simd::yyyy(v), m.row1));
	t = simd::add4(t, simd::mul4(simd::zzzz(v), m.row2));
	return Vector4<float>(simd::add4(t, simd::mul4(simd::wwww(v), m.row3)));
}

template<>
inline Vector4<float> operator*(const Matrix4<float>& m, const Vector4<float>& v) noexcept
{
	return Vector4<float>(simd::set4(simd::dot4(m.row0, v), simd::dot4(m.row1, v), simd::dot4(m.row2, v), simd::dot4(m.row3, v)));
}

inline Vector4<float>& Vector4<float>::transform(const Matrix4<float>& m)
{
	*this *= m;
	return *this;
}

template<>
inline Vector4<float> transform(const Vector4<float>& v, const Matrix4<float>& m) noexcept
{
	return v*m;
}

template<>
inline Matrix4<float> tensor(const Vector4<float>& v1, const Vector4<float>& v2) noexcept
{
	return Matrix4<float>(simd::mul4(simd::xxxx(v1), v2),
		simd::mul4(simd::yyyy(v1), v2),
		simd::mul4(simd::zzzz(v1), v2),
		simd::mul4(simd::wwww(v1), v2));
}

#endif /* SIMD_HAS_FLOAT4 */

} // namespace mathematics::templates
