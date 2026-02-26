/*
 *	Name: Vector2
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
//#include <Tuples/Tuple2.hpp>
#include "Constants.hpp"
#include "Axis.hpp"

namespace core::mathematics {

struct Uninitialized
{
};

constexpr Uninitialized UNINITIALIZED{};

struct Normalized
{
};

struct Unnormalized
{
};

template<typename T>
concept NormalizationType = (std::same_as<T, Normalized> || std::same_as<T, Unnormalized>); // #TODO Move to Concepts.hpp

namespace templates {

template<typename T>
	requires std::floating_point<T>
struct Matrix2;

template<typename T>
	requires (std::floating_point<T> || std::integral<T>)
struct Vector2;

template<std::floating_point T>
struct Vector2<T>
{
	using Real = T;
	using ComponentType = T;
	using ConstArg = const Vector2&;
	using ConstResult = const Vector2&;

	static constexpr int NUM_COMPONENTS = 2;

	constexpr Vector2() noexcept : x(), y() {}
	explicit Vector2(Uninitialized) noexcept {}
	constexpr explicit Vector2(T scalar) noexcept : x(scalar), y(scalar) {}
	constexpr Vector2(T x, T y) noexcept : x(x), y(y) {}
	//explicit Vector2(const tuples::templates::Tuple2<T>& t) noexcept : x(t.x), y(t.y) {}
	//template<typename U> explicit Vector2(const tuples::templates::Tuple2<U>& t) noexcept : x(T(t.x)), y(T(t.y)) {}
	explicit Vector2(const std::pair<T, T>& t) noexcept : x(t.first), y(t.second) {}
	template<typename U> explicit Vector2(const std::pair<U, U>& t) noexcept : x(T(t.first)), y(T(t.second)) {}
	explicit Vector2(const std::tuple<T, T>& t) noexcept : x(std::get<0>(t)), y(std::get<1>(t)) {}
	template<typename U> explicit Vector2(const std::tuple<U, U>& t) noexcept : x(T(std::get<0>(t))), y(T(std::get<1>(t))) {}
	explicit Vector2(const T* v) noexcept : x(v[0]), y(v[1]) {}
	explicit Vector2(Axis axis) noexcept : x((axis == Axis::X) ? T(1) : T(0)), y((axis == Axis::Y) ? T(1) : T(0)) {}
	template<typename U> explicit Vector2(const Vector2<U>& v) noexcept : x(T(t.x)), y(T(t.y)) {}

	//explicit operator tuples::templates::Tuple2<T>() noexcept { return tuples::templates::Tuple2<T>(x, y); }
	//template<typename U> explicit operator tuples::templates::Tuple2<U>() noexcept { return tuples::templates::Tuple2<U>(U(x), U(y)); }
	//explicit operator std::pair<T, T>() { return std::pair<T, T>(x, y); }
	//template<typename U> explicit operator std::pair<U, U>() { return std::pair<U, U>(U(x), U(y)); }
	//explicit operator std::tuple<T, T>() { return std::tuple<T, T>(x, y); }
	//template<typename U> explicit operator std::tuple<U, U>() { return std::tuple<U, U>(U(x), U(y)); }
	explicit operator T*() noexcept { return &x; }
	explicit operator const T*() const noexcept { return &x; }
	T& operator[](int i) noexcept { return (&x)[i]; }
	const T& operator[](int i) const noexcept { return (&x)[i]; }

	Vector2 operator+() const noexcept { return *this; }
	Vector2 operator-() const noexcept { return Vector2(-x, -y); }
	Vector2& operator+=(const Vector2& v) noexcept { x += v.x; y += v.y; return *this; }
	Vector2& operator-=(const Vector2& v) noexcept { x -= v.x; y -= v.y; return *this; }
	Vector2& operator*=(const Vector2& v) noexcept { x *= v.x; y *= v.y; return *this; }
	Vector2& operator*=(T f) noexcept { x *= f; y *= f; return *this; }
	Vector2& operator*=(const Matrix2<T>& m) noexcept;
	Vector2& operator/=(const Vector2& v) noexcept { x /= v.x; y /= v.y; return *this; }
	Vector2& operator/=(T f) noexcept { return operator*=(T(1)/f); }
	bool operator==(const Vector2& v) const noexcept { return (x == v.x) && (y == v.y); }
	bool operator!=(const Vector2& v) const noexcept { return !(*this == v); }
	
	template<typename A> void load(A& ar) { ar(x, y); }
	template<typename A> void save(A& ar) const { ar(x, y); }

	template<std::size_t I> T& get() noexcept;
	template<std::size_t I> const T& get() const noexcept;

	bool isZero() const noexcept { return (x == T()) && (y == T()); }
	bool isApproxZero() const noexcept;
	bool isApproxEqual(const Vector2& v) const noexcept;
	bool isApproxEqual(const Vector2& v, T tolerance) const noexcept;
	bool allLessThan(const Vector2& v) const noexcept { return (x < v.x) && (y < v.y); }
	bool allLessThanEqual(const Vector2& v) const noexcept { return (x <= v.x) && (y <= v.y); }
	bool allGreaterThan(const Vector2& v) const noexcept { return (x > v.x) && (y > v.y); }
	bool allGreaterThanEqual(const Vector2& v) const noexcept { return (x >= v.x) && (y >= v.y); }
	bool anyLessThan(const Vector2& v) const noexcept { return (x < v.x) || (y < v.y); }
	bool anyLessThanEqual(const Vector2& v) const noexcept { return (x <= v.x) || (y <= v.y); }
	bool anyGreaterThan(const Vector2& v) const noexcept { return (x > v.x) || (y > v.y); }
	bool anyGreaterThanEqual(const Vector2& v) const noexcept { return (x >= v.x) || (y >= v.y); }
	bool isFinite() const noexcept { return std::isfinite(x) && std::isfinite(y); }
	T getMagnitude() const noexcept { return std::sqrt(x*x + y*y); }
	T getMagnitudeSquared() const noexcept { return (x*x + y*y); }
	void setMagnitude(T magnitude) noexcept;
	T getLength() const noexcept { return getMagnitude(); }
	T getLengthSquared() const noexcept { return getMagnitudeSquared(); }
	void setLength(T length) noexcept { setMagnitude(length); }
	Axis getMajorAxis() const noexcept { return (std::fabs(y) > std::fabs(x)) ? Axis::Y : Axis::X; }
	T getMinComponent() const { return std::min(x, y); }
	T getMaxComponent() const { return std::max(x, y); }
	Vector2& setZero/*zero*/() noexcept { x = T(); y = T(); return *this; }
	Vector2& set(T x, T y) noexcept { this->x = x; this->y = y; return *this; }
	Vector2& setMinimum(const Vector2& v1, const Vector2& v2);
	Vector2& setMaximum(const Vector2& v1, const Vector2& v2);
	Vector2& negate() noexcept { x = -x; y = -y; return *this; }
	Vector2& normalize() noexcept;
	Vector2& rotate(T angle) noexcept;
	//Vector2& scale(const Vector2& v) noexcept { x *= v.x; y *= v.y; return *this; }
	Vector2& transform(const Matrix2<T>& m) noexcept;

	//static const Vector2& getZero() noexcept { return ZERO; }
	//static const Vector2& getUnitX() noexcept { return UNIT_X; }
	//static const Vector2& getUnitY() noexcept { return UNIT_Y; }

	static const Vector2 ZERO;
	static const Vector2 UNIT_X;
	static const Vector2 UNIT_Y;
	static const Vector2 ONE;
	static const Vector2 TOLERANCE;
	static const Vector2 INF;
	static const Vector2 MINUS_INF;

	T x, y;
};

template<std::floating_point T> const Vector2<T> Vector2<T>::ZERO{};
template<std::floating_point T> const Vector2<T> Vector2<T>::UNIT_X{ T(1), T(0) };
template<std::floating_point T> const Vector2<T> Vector2<T>::UNIT_Y{ T(0), T(1) };
template<std::floating_point T> const Vector2<T> Vector2<T>::ONE{ T(1), T(1) };
template<std::floating_point T> const Vector2<T> Vector2<T>::TOLERANCE{ Constants<T>::TOLERANCE, Constants<T>::TOLERANCE };
template<std::floating_point T> const Vector2<T> Vector2<T>::INF{ std::numeric_limits<T>::infinity(), std::numeric_limits<T>::infinity() };
template<std::floating_point T> const Vector2<T> Vector2<T>::MINUS_INF{ -std::numeric_limits<T>::infinity(), -std::numeric_limits<T>::infinity() };

template<std::integral T>
struct Vector2<T>
{
	using ComponentType = T;
	using ConstArg = const Vector2&;
	using ConstResult = const Vector2&;

	static constexpr int NUM_COMPONENTS = 2;

	constexpr Vector2() noexcept : x(), y() {}
	explicit Vector2(Uninitialized) noexcept {}
	constexpr explicit Vector2(T scalar) noexcept : x(scalar), y(scalar) {}
	constexpr Vector2(T x, T y) noexcept : x(x), y(y) {}
	//explicit Vector2(const tuples::templates::Tuple2<T>& t) noexcept : x(t.x), y(t.y) {}
	//template<typename U> explicit Vector2(const tuples::templates::Tuple2<U>& t) noexcept : x(T(t.x)), y(T(t.y)) {}
	explicit Vector2(const std::pair<T, T>& t) noexcept : x(t.first), y(t.second) {}
	template<typename U> explicit Vector2(const std::pair<U, U>& t) noexcept : x(T(t.first)), y(T(t.second)) {}
	explicit Vector2(const std::tuple<T, T>& t) noexcept : x(std::get<0>(t)), y(std::get<1>(t)) {}
	template<typename U> explicit Vector2(const std::tuple<U, U>& t) noexcept : x(T(std::get<0>(t))), y(T(std::get<1>(t))) {}
	explicit Vector2(const T* v) noexcept : x(v[0]), y(v[1]) {}
	template<typename U> explicit Vector2(const Vector2<U>& v) noexcept : x(T(v.x)), y(T(v.y)) {}

	//explicit operator tuples::templates::Tuple2<T>() noexcept { return tuples::templates::Tuple2<T>(x, y); }
	//template<typename U> explicit operator tuples::templates::Tuple2<U>() noexcept { return tuples::templates::Tuple2<U>(U(x), U(y)); }
	//explicit operator std::pair<T, T>() { return std::pair<T, T>(x, y); }
	//template<typename U> explicit operator std::pair<U, U>() { return std::pair<U, U>(U(x), U(y)); }
	//explicit operator std::tuple<T, T>() { return std::tuple<T, T>(x, y); }
	//template<typename U> explicit operator std::tuple<U, U>() { return std::tuple<U, U>(U(x), U(y)); }
	explicit operator T* () noexcept { return &x; }
	explicit operator const T* () const noexcept { return &x; }
	T& operator[](int i) noexcept { return (&x)[i]; }
	const T& operator[](int i) const noexcept { return (&x)[i]; }

	Vector2 operator+() const noexcept { return *this; }
	Vector2 operator-() const noexcept { return Vector2(-x, -y); }
	Vector2& operator+=(const Vector2& v) noexcept { x += v.x; y += v.y; return *this; }
	Vector2& operator-=(const Vector2& v) noexcept { x -= v.x; y -= v.y; return *this; }
	Vector2& operator*=(const Vector2& v) noexcept { x *= v.x; y *= v.y; return *this; }
	Vector2& operator*=(T f) noexcept { x *= f; y *= f; return *this; }
	Vector2& operator/=(const Vector2& v) noexcept { x /= v.x; y /= v.y; return *this; }
	Vector2& operator/=(T f) noexcept { x /= f; y /= f; return *this; }
	bool operator==(const Vector2& v) const noexcept { return (x == v.x) && (y == v.y); }
	bool operator!=(const Vector2& v) const noexcept { return !(*this == v); }

	template<typename A> void serialize(A& ar) { ar(x, y); }

	template<std::size_t I> T& get() noexcept;
	template<std::size_t I> const T& get() const noexcept;

	bool isZero() const noexcept { return (x == T()) && (y == T()); }
	bool allLessThan(const Vector2& v) const noexcept { return (x < v.x) && (y < v.y); }
	bool allLessThanEqual(const Vector2& v) const noexcept { return (x <= v.x) && (y <= v.y); }
	bool allGreaterThan(const Vector2& v) const noexcept { return (x > v.x) && (y > v.y); }
	bool allGreaterThanEqual(const Vector2& v) const noexcept { return (x >= v.x) && (y >= v.y); }
	bool anyLessThan(const Vector2& v) const noexcept { return (x < v.x) || (y < v.y); }
	bool anyLessThanEqual(const Vector2& v) const noexcept { return (x <= v.x) || (y <= v.y); }
	bool anyGreaterThan(const Vector2& v) const noexcept { return (x > v.x) || (y > v.y); }
	bool anyGreaterThanEqual(const Vector2& v) const noexcept { return (x >= v.x) || (y >= v.y); }
	T getMinComponent() const { return std::min(x, y); }
	T getMaxComponent() const { return std::max(x, y); }
	Vector2& setZero/*zero*/() noexcept { x = T(); y = T(); return *this; }
	Vector2& set(T x, T y) noexcept { this->x = x; this->y = y; return *this; }
	Vector2& setMinimum(const Vector2& v1, const Vector2& v2);
	Vector2& setMaximum(const Vector2& v1, const Vector2& v2);
	Vector2& negate() noexcept { x = -x; y = -y; return *this; }

	//static const Vector2& getZero() noexcept { return ZERO; }

	static const Vector2 ZERO;

	T x, y;
};

template<std::integral T> const Vector2<T> Vector2<T>::ZERO{};

#if SIMD_HAS_FLOAT4

template<>
struct Matrix2<float>;

template<>
struct Vector2<float>
{
	using Real = float;
	using ComponentType = float;
	using ConstArg = const Vector2;
	using ConstResult = const Vector2;

	static constexpr int NUM_COMPONENTS = 2;

	/*constexpr*/ Vector2() noexcept : xy(simd::zero<simd::float4>()) {}
	explicit Vector2(Uninitialized) noexcept {}
#if MATHEMATICS_SIMD_EXPAND_LAST
	/*constexpr*/ explicit Vector2(float scalar) noexcept : xy(simd::set4(scalar)) {}
	/*constexpr*/ Vector2(float x, float y) noexcept : xy(simd::set4(x, y, y, y)) {}
	//explicit Vector2(const tuples::templates::Tuple2<float>& t) noexcept : xy(simd::set4(t.x, t.y, t.y, t.y)) {}
	//template<typename U> explicit Vector2(const tuples::templates::Tuple2<U>& t) : Vector2((float)t.x, (float)t.y) {}
	explicit Vector2(const std::pair<float, float>& t) noexcept : xy(simd::set4(t.first, t.second, t.second, t.second)) {}
	template<typename U> explicit Vector2(const std::pair<U, U>& t) noexcept : Vector2((float)t.first, (float)t.second) {}
	explicit Vector2(const std::tuple<float, float>& t) noexcept : Vector2(std::get<0>(t), std::get<1>(t)) {}
	template<typename U> explicit Vector2(const std::tuple<U, U>& t) noexcept : Vector2((float)std::get<0>(t), (float)std::get<1>(t)) {}
	explicit Vector2(const float* v) noexcept : Vector2(v[0], v[1]) {}
#else
	/*constexpr*/ explicit Vector2(float scalar) noexcept : xy(simd::set2(scalar)) {}
	/*constexpr*/ Vector2(float x, float y) noexcept : xy(simd::set2(x, y)) {}
	//explicit Vector2(const tuples::templates::Tuple2<float>& t) noexcept : xy(simd::set2(t.x, t.y)) {}
	//template<typename U> explicit Vector2(const tuples::templates::Tuple2<U>& t) noexcept : xy(simd::set2((float)t.x, (float)t.y)) {}
	explicit Vector2(const std::pair<float, float>& t) noexcept : xy(simd::set2(t.first, t.second)) {}
	template<typename U> explicit Vector2(const std::pair<U, U>& t) noexcept : xy(simd::set2((float)t.first, (float)t.second)) {}
	explicit Vector2(const std::tuple<float, float>& t) noexcept : xy(simd::set2(std::get<0>(t), std::get<1>(t))) {}
	template<typename U> explicit Vector2(const std::tuple<U, U>& t) noexcept : xy(simd::set2((float)std::get<0>(t), (float)std::get<1>(t))) {}
	explicit Vector2(const float* v) noexcept : xy(simd::load2(v)) {}
#endif
	explicit Vector2(Axis axis) noexcept : Vector2((axis == Axis::X) ? 1.f : 0.f, (axis == Axis::Y) ? 1.f : 0.f) {}
	explicit Vector2(simd::float4 v) noexcept : xy(v) {}
	Vector2(const Vector2& v) noexcept : xy(v.xy) {}
	template<typename U> explicit Vector2(const Vector2<U>& v) noexcept : Vector2((float)v.x, (float)v.y)) {}
	Vector2& operator=(const Vector2& v) noexcept { xy = v.xy; return *this; }

	operator simd::float4() const noexcept { return xy; }
	//explicit operator tuples::templates::Tuple2<float>() noexcept { return tuples::templates::Tuple2<float>(x, y); }
	//template<typename U> explicit operator tuples::templates::Tuple2<U>() noexcept { return tuples::templates::Tuple2<U>(U(x), U(y)); }
	//explicit operator std::pair<float, float>() { return std::pair<float, float>(x, y); }
	//template<typename U> explicit operator std::pair<U, U>() { return std::pair<U, U>(U(x), U(y)); }
	//explicit operator std::tuple<float, float>() { return std::tuple<float, float>(x, y); }
	//template<typename U> explicit operator std::tuple<U, U>() { return std::tuple<U, U>(U(x), U(y)); }
	explicit operator float* () noexcept { return &x; }
	explicit operator const float* () const noexcept { return &x; }
	float& operator[](int i) noexcept { return (&x)[i]; }
	const float& operator[](int i) const noexcept { return (&x)[i]; }

	Vector2 operator+() const noexcept { return *this; }
#if MATHEMATICS_SIMD_EXPAND_LAST
	Vector2 operator-() const noexcept { return Vector2(simd::neg4(xy)); }
#else
	Vector2 operator-() const noexcept { return Vector2(simd::neg2(xy)); }
#endif
	Vector2& operator+=(const Vector2& v) noexcept { xy = simd::add4(xy, v); return *this; }
	Vector2& operator-=(const Vector2& v) noexcept { xy = simd::sub4(xy, v); return *this; }
	Vector2& operator*=(const Vector2& v) noexcept { xy = simd::mul4(xy, v); return *this; }
	Vector2& operator*=(float f) noexcept { xy = simd::mul4(xy, simd::set4(f)); return *this; }
	Vector2& operator*=(const Matrix2<float>& m) noexcept;
#if MATHEMATICS_SIMD_EXPAND_LAST
	Vector2& operator/=(const Vector2& v) noexcept { xy = simd::div4(xy, v); return *this; }
#else
	Vector2& operator/=(const Vector2& v) noexcept { xy = simd::div2(xy, v); return *this; }
#endif
	Vector2& operator/=(float f) noexcept { xy = simd::div4(xy, simd::set4(f)); return *this; }
	bool operator==(const Vector2& v) const noexcept { return simd::all2(simd::equal(xy, v)); }
	bool operator!=(const Vector2& v) const noexcept { return !(*this == v); }

	template<typename A> void load(A& ar);
	template<typename A> void save(A& ar) const { ar(x, y); }

	template<std::size_t I> float& get() noexcept;
	template<std::size_t I> const float& get() const noexcept;
	template<typename U> U& get() noexcept;				// intentionally undefined
	template<typename U> const U& get() const noexcept;	// intentionally undefined
	template<> simd::float4& get() noexcept { return xy; }
	template<> const simd::float4& get() const noexcept { return xy; }

	bool isZero() const noexcept { return simd::all2(simd::equal(xy, simd::zero<simd::float4>())); }
	bool isApproxZero() const noexcept { return simd::all2(simd::lessThan(simd::abs4(xy), TOLERANCE)); }
	bool isApproxEqual(const Vector2& v) const noexcept { return simd::all2(simd::lessThan(simd::abs4(simd::sub4(xy, v)), TOLERANCE)); }
	bool isApproxEqual(const Vector2& v, float tolerance) const noexcept { return simd::all2(simd::lessThan(simd::abs4(simd::sub4(xy, v)), simd::set4(tolerance))); }
	bool allLessThan(const Vector2& v) const noexcept { return simd::all2(simd::lessThan(xy, v)); }
	bool allLessThanEqual(const Vector2& v) const noexcept { return simd::all2(simd::lessThanEqual(xy, v)); }
	bool allGreaterThan(const Vector2& v) const noexcept { return simd::all2(simd::greaterThan(xy, v)); }
	bool allGreaterThanEqual(const Vector2& v) const noexcept { return simd::all2(simd::greaterThanEqual(xy, v)); }
	bool anyLessThan(const Vector2& v) const noexcept { return simd::any2(simd::lessThan(xy, v)); }
	bool anyLessThanEqual(const Vector2& v) const noexcept { return simd::any2(simd::lessThanEqual(xy, v)); }
	bool anyGreaterThan(const Vector2& v) const noexcept { return simd::any2(simd::greaterThan(xy, v)); }
	bool anyGreaterThanEqual(const Vector2& v) const noexcept { return simd::any2(simd::greaterThanEqual(xy, v)); }
	bool isFinite() const noexcept { return simd::all2(simd::isFinite(xy)); }
	float getMagnitude() const noexcept { return simd::toFloat(simd::sqrt1(simd::dot2(xy, xy))); }
	float getMagnitudeSquared() const noexcept { return simd::toFloat(simd::dot2(xy, xy)); }
	void setMagnitude(float magnitude) noexcept;
	float getLength() const noexcept { return getMagnitude(); }
	float getLengthSquared() const noexcept { return getMagnitudeSquared(); }
	void setLength(float length) noexcept { setMagnitude(length); }
	Axis getMajorAxis() const noexcept { (Axis)simd::asIndex(simd::equal(xy, simd::hMax2(xy))); }
	float getMinComponent() const noexcept { return simd::toFloat(simd::hMin2(xy)); }
	float getMaxComponent() const noexcept { return simd::toFloat(simd::hMax2(xy)); }
	Vector2& setZero/*zero*/() noexcept { xy = simd::zero<simd::float4>(); return *this; }
#if MATHEMATICS_SIMD_EXPAND_LAST
	Vector2& set(float x, float y) noexcept { xy = simd::set4(x, y, y, y); return *this; }
#else
	Vector2& set(float x, float y) noexcept { xy = simd::set2(x, y); return *this; }
#endif
	Vector2& setMinimum(const Vector2& v1, const Vector2& v2) noexcept { xy = simd::min4(v1, v2); return *this; }
	Vector2& setMaximum(const Vector2& v1, const Vector2& v2) noexcept { xy = simd::max4(v1, v2); return *this; }
#if MATHEMATICS_SIMD_EXPAND_LAST
	Vector2& negate() noexcept { xy = simd::neg4(xy); return *this; }
#else
	Vector2& negate() noexcept { xy = simd::neg2(xy); return *this; }
#endif
	Vector2& normalize() noexcept;
	Vector2& rotate(float angle) noexcept;
	//Vector2& scale(const Vector2& v) noexcept { xy = simd::mul4(xy, v); return *this; }
	Vector2& transform(const Matrix2<float>& m) noexcept;

	//static const Vector2& getZero() noexcept { return ZERO; }
	//static const Vector2& getUnitX() noexcept { return UNIT_X; }
	//static const Vector2& getUnitY() noexcept { return UNIT_Y; }

	static const Vector2 ZERO;
	static const Vector2 UNIT_X;
	static const Vector2 UNIT_Y;
	static const Vector2 ONE;
	static const Vector2 TOLERANCE;
	static const Vector2 INF;
	static const Vector2 MINUS_INF;

	union
	{
		simd::float4 xy;
		struct { float x, y; };
	};
};

const Vector2<float> Vector2<float>::ZERO{};
const Vector2<float> Vector2<float>::UNIT_X{ 1.f, 0.f };
const Vector2<float> Vector2<float>::UNIT_Y{ 0.f, 1.f };
const Vector2<float> Vector2<float>::ONE{ 1.f, 1.f };
const Vector2<float> Vector2<float>::TOLERANCE{ Constants<float>::TOLERANCE, Constants<float>::TOLERANCE };
const Vector2<float> Vector2<float>::INF{ std::numeric_limits<float>::infinity(), std::numeric_limits<float>::infinity() };
const Vector2<float> Vector2<float>::MINUS_INF{ -std::numeric_limits<float>::infinity(), -std::numeric_limits<float>::infinity() };

#endif /* SIMD_HAS_FLOAT4 */

template<typename T>
	requires (std::floating_point<T> || std::integral<T>)
inline Vector2<T> operator+(const Vector2<T>& v1, const Vector2<T>& v2) noexcept
{ 
	return Vector2<T>(v1.x + v2.x, v1.y + v2.y);
}

template<typename T>
	requires (std::floating_point<T> || std::integral<T>)
inline Vector2<T> operator-(const Vector2<T>& v1, const Vector2<T>& v2) noexcept
{ 
	return Vector2<T>(v1.x - v2.x, v1.y - v2.y);
}

template<typename T>
	requires (std::floating_point<T> || std::integral<T>)
inline Vector2<T> operator*(const Vector2<T>& v1, const Vector2<T>& v2) noexcept
{ 
	return Vector2<T>(v1.x*v2.x, v1.y*v2.y);
}

template<typename T>
	requires (std::floating_point<T> || std::integral<T>)
inline Vector2<T> operator*(T f, const Vector2<T>& v) noexcept
{ 
	return Vector2<T>(f*v.x, f*v.y);
}

template<typename T>
	requires (std::floating_point<T> || std::integral<T>)
inline Vector2<T> operator*(const Vector2<T>& v, T f) noexcept
{ 
	return Vector2<T>(v.x*f, v.y*f);
}

template<typename T>
	requires (std::floating_point<T> || std::integral<T>)
inline Vector2<T> operator/(const Vector2<T>& v1, const Vector2<T>& v2) noexcept
{ 
	return Vector2<T>(v1.x/v2.x, v1.y/v2.y);
}

template<typename T>
	requires (std::floating_point<T> || std::integral<T>)
inline Vector2<T> operator/(T f, const Vector2<T>& v) noexcept
{ 
	return Vector2<T>(f/v.x, f/v.y);
}

template<typename T>
	requires (std::floating_point<T> || std::integral<T>)
inline Vector2<T> operator/(const Vector2<T>& v, T f) noexcept
{ 
	if constexpr (std::is_floating_point_v<T>)
		return operator*(v, T(1)/f);
	else
		return Vector2<T>(v.x/f, v.y/f);
}

template<typename C, typename T, typename U>
	requires std::floating_point<U> || std::integral<U>
inline std::basic_istream<C, T>& operator>>(std::basic_istream<C, T>& s, Vector2<U>& v)
{ 
	return s >> v.x >> std::ws >> v.y; 
}

template<typename C, typename T, typename U>
	requires std::floating_point<U> || std::integral<U>
inline std::basic_ostream<C, T>& operator<<(std::basic_ostream<C, T>& s, const Vector2<U>& v)
{ 
	constexpr C WS(0x20);
	return s << v.x << WS << v.y;
}

template<std::floating_point T>
template<std::size_t I>
inline T& Vector2<T>::get()
{
	if constexpr (I == 0)
		return x;
	else if constexpr (I == 1)
		return y;
	static_assert(false);
}

template<std::floating_point T>
template<std::size_t I>
inline const T& Vector2<T>::get() const
{
	if constexpr (I == 0)
		return x;
	else if constexpr (I == 1)
		return y;
	static_assert(false);
}

template<std::floating_point T>
inline bool Vector2<T>::isApproxZero() const
{ 
	return (std::fabs(x) < Constants<T>::TOLERANCE) && (std::fabs(y) < Constants<T>::TOLERANCE); 
}

template<std::floating_point T>
inline bool Vector2<T>::isApproxEqual(const Vector2<T>& v) const
{ 
	return (std::fabs(v.x - x) < Constants<T>::TOLERANCE) && (std::fabs(v.y - y) < Constants<T>::TOLERANCE); 
}

template<std::floating_point T>
inline bool Vector2<T>::isApproxEqual(const Vector2<T>& v, T tolerance) const
{ 
	return (std::fabs(v.x - x) < tolerance) && (std::fabs(v.y - y) < tolerance); 
}

template<std::floating_point T>
inline void Vector2<T>::setMagnitude(T magnitude) 
{ 
	T m = getMagnitude(); 
	if (m > T(0)) 
		*this *= magnitude/m;
}

template<std::floating_point T>
inline Vector2<T>& Vector2<T>::setMinimum(const Vector2<T>& v1, const Vector2<T>& v2)
{ 
	x = std::min(v1.x, v2.x); 
	y = std::min(v1.y, v2.y); 
	return *this; 
}

template<std::floating_point T>
inline Vector2<T>& Vector2<T>::setMaximum(const Vector2<T>& v1, const Vector2<T>& v2)
{ 
	x = std::max(v1.x, v2.x); 
	y = std::max(v1.y, v2.y); 
	return *this; 
}

template<std::floating_point T>
inline Vector2<T>& Vector2<T>::normalize()
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
inline Vector2<T>& Vector2<T>::rotate(T angle)
{
	if (angle != T(0))
	{
		T sinAngle = std::sin(angle);
		T cosAngle = std::cos(angle);
		set(x*cosAngle - y*sinAngle, y*cosAngle + x*sinAngle);
	}

	return *this;
}

template<std::integral T>
template<std::size_t I>
inline T& Vector2<T>::get()
{
	if constexpr (I == 0)
		return x;
	else if constexpr (I == 1)
		return y;
	static_assert(false);
}

template<std::integral T>
template<std::size_t I>
inline const T& Vector2<T>::get() const
{
	if constexpr (I == 0)
		return x;
	else if constexpr (I == 1)
		return y;
	static_assert(false);
}

template<std::integral T>
inline Vector2<T>& Vector2<T>::setMinimum(const Vector2<T>& v1, const Vector2<T>& v2)
{
	x = std::min(v1.x, v2.x);
	y = std::min(v1.y, v2.y);
	return *this;
}

template<std::integral T>
inline Vector2<T>& Vector2<T>::setMaximum(const Vector2<T>& v1, const Vector2<T>& v2)
{
	x = std::max(v1.x, v2.x);
	y = std::max(v1.y, v2.y);
	return *this;
}

#if SIMD_HAS_FLOAT4

template<>
inline Vector2<float> operator+(const Vector2<float>& v1, const Vector2<float>& v2) noexcept 
{ 
	return Vector2<float>(simd::add4(v1, v2)); 
}

template<>
inline Vector2<float> operator-(const Vector2<float>& v1, const Vector2<float>& v2) noexcept 
{ 
	return Vector2<float>(simd::sub4(v1, v2)); 
}

template<>
inline Vector2<float> operator*(const Vector2<float>& v1, const Vector2<float>& v2) noexcept 
{ 
	return Vector2<float>(simd::mul4(v1, v2)); 
}

template<>
inline Vector2<float> operator*(float f, const Vector2<float>& v) noexcept 
{ 
	return Vector2<float>(simd::mul4(simd::set4(f), v)); 
}

template<>
inline Vector2<float> operator*(const Vector2<float>& v, float f) noexcept 
{ 
	return Vector2<float>(simd::mul4(v, simd::set4(f))); 
}

template<>
inline Vector2<float> operator/(const Vector2<float>& v1, const Vector2<float>& v2) noexcept 
{ 
#if IMAGING_SIMD_EXPAND_LAST
	return Vector2<float>(simd::div4(v1, v2));
#else
	return Vector2<float>(simd::div2(v1, v2));
#endif
}

template<>
inline Vector2<float> operator/(float f, const Vector2<float>& v) noexcept 
{ 
#if IMAGING_SIMD_EXPAND_LAST
	return Vector2<float>(simd::div4(simd::set4(f), v));
#else
	return Vector2<float>(simd::div2(simd::set2(f), v));
#endif
}

template<>
inline Vector2<float> operator/(const Vector2<float>& v, float f) noexcept 
{ 
	return Vector2<float>(simd::div4(v, simd::set4(f))); 
}

template<typename C, typename T>
inline std::basic_istream<C, T>& operator>>(std::basic_istream<C, T>& s, Vector2<float>& v)
{ 
	float x, y; 
	s >> x >> std::ws >> y; 
	v.set(x, y); 
	return s; 
}

template<typename A>
inline void Vector2<float>::load(A& ar)
{
	float t0, t1;
	ar(t0, t1);
	set(t0, t1);
}

template<std::size_t I>
inline float& Vector2<float>::get()
{
	if constexpr (I == 0)
		return x;
	else if constexpr (I == 1)
		return y;
	static_assert(false);
}

template<std::size_t I>
inline const float& Vector2<float>::get() const
{
	if constexpr (I == 0)
		return x;
	else if constexpr (I == 1)
		return y;
	static_assert(false);
}

inline void Vector2<float>::setMagnitude(float magnitude)
{ 
	float m = getMagnitude();
	if (m > 0.f) 
		*this *= magnitude/m;
}

inline Vector2<float>& Vector2<float>::normalize()
{
#if MATHEMATICS_FAST_NORMALIZE
	float m = simd::toFloat(simd::rcpSqrtApprox1(simd::dot2(xy, xy)));
	if (m <= std::numeric_limits<float>::max()) 
		xy = simd::mul4(xy, simd::set4(m));
#else
	float m = getMagnitude(); 
	if (m > 0.f) 
		xy = simd::div4(xy, simd::set4(m));
#endif
	return *this;
}

inline Vector2<float>& Vector2<float>::rotate(float angle)
{
	if (angle != 0.f)
	{
		auto sinAngle = simd::set2(std::sin(angle));
		auto cosAngle = simd::set2(std::cos(angle));
		xy = simd::subAdd4(simd::mul4(xy, cosAngle), simd::mul4(simd::yxxx(xy), sinAngle));
	}

	return *this;
}

#endif /* SIMD_HAS_FLOAT4 */

template<std::size_t I, typename T>
	requires (std::floating_point<T> || std::integral<T>)
inline T& get(Vector2<T>& v) noexcept
{
	if constexpr (I == 0)
		return v.x;
	else if constexpr (I == 1)
		return v.y;
	static_assert(false);
}

template<std::size_t I, typename T>
	requires (std::floating_point<T> || std::integral<T>)
inline const T& get(const Vector2<T>& v) noexcept
{
	if constexpr (I == 0)
		return v.x;
	else if constexpr (I == 1)
		return v.y;
	static_assert(false);
}

template<std::size_t I, typename T>
	requires (std::floating_point<T> || std::integral<T>)
inline T&& get(Vector2<T>&& v) noexcept
{
	if constexpr (I == 0)
		return v.x;
	else if constexpr (I == 1)
		return v.y;
	static_assert(false);
}

template<std::size_t I, typename T>
	requires (std::floating_point<T> || std::integral<T>)
inline const T&& get(const Vector2<T>&& v) noexcept
{
	if constexpr (I == 0)
		return v.x;
	else if constexpr (I == 1)
		return v.y;
	static_assert(false);
}

template<typename T>
	requires std::floating_point<T>
inline T dot(const Vector2<T>& v1, const Vector2<T>& v2) noexcept
{
	return v1.x*v2.x + v1.y*v2.y;
}

template<typename T>
	requires std::floating_point<T>
inline T cross(const Vector2<T>& v1, const Vector2<T>& v2) noexcept
{
	return v1.x*v2.y - v1.y*v2.x;
}

template<typename T>
	requires std::floating_point<T>
inline T distance(const Vector2<T>& v1, const Vector2<T>& v2) noexcept
{
	T x = v2.x - v1.x;
	T y = v2.y - v1.y;
	return std::sqrt(x*x + y*y);
}

template<typename T>
	requires std::floating_point<T>
inline T distanceSquared(const Vector2<T>& v1, const Vector2<T>& v2) noexcept
{
	T x = v2.x - v1.x;
	T y = v2.y - v1.y;
	return x*x + y*y;
}

template<typename T>
	requires std::floating_point<T>
inline T length(const Vector2<T>& v) noexcept
{
	return std::sqrt(v.x*v.x + v.y*v.y);
}

template<typename T>
	requires std::floating_point<T>
inline T lengthSquared(const Vector2<T>& v) noexcept
{
	return v.x*v.x + v.y*v.y;
}

template<typename T>
	requires (std::floating_point<T> || std::integral<T>)
inline Vector2<T> minimum(const Vector2<T>& v1, const Vector2<T>& v2)
{
	return Vector2<T>(std::min(v1.x, v2.x), std::min(v1.y, v2.y));
}

template<typename T>
	requires (std::floating_point<T> || std::integral<T>)
inline Vector2<T> maximum(const Vector2<T>& v1, const Vector2<T>& v2)
{
	return Vector2<T>(std::max(v1.x, v2.x), std::max(v1.y, v2.y));
}

template<typename T>
	requires std::floating_point<T>
inline Vector2<T> lerp(const Vector2<T>& v1, const Vector2<T>& v2, T t) noexcept
{
	return Vector2<T>(v1.x + t*(v2.x - v1.x), v1.y + t*(v2.y - v1.y));
}

template<typename T>
	requires std::floating_point<T>
inline Vector2<T> slerp(const Vector2<T>& v1, const Vector2<T>& v2, T t)
{
	T dp = v1.x*v2.x + v1.y*v2.y;
	if ((T(1) - dp) < Constants<T>::TOLERANCE)
	{
		Vector2<T> c(v1.x + t*(v2.x - v1.x), v1.y + t*(v2.y - v1.y));
		c.normalize();
		return c;
	}

	dp = std::clamp(dp, T(-1), T(1));
	T theta = std::acos(dp)*t;
	T st = std::sin(theta);
	T ct = std::cos(theta);
	Vector2<T> c(v2.x - v1.x*dp, v2.y - v1.y*dp);
	c.normalize();
	return Vector2<T>(v1.x*ct + c.x*st, v1.y*ct + c.y*st);
}

template<typename T>
	requires std::floating_point<T>
inline Vector2<T> perpendicular(const Vector2<T>& v) noexcept
{
	return Vector2<T>(-v.y, v.x);
}

#if SIMD_HAS_FLOAT4

template<>
inline float dot(const Vector2<float>& v1, const Vector2<float>& v2) noexcept
{
	return simd::toFloat(simd::dot2(v1, v2));
}

template<>
inline float cross(const Vector2<float>& v1, const Vector2<float>& v2) noexcept
{
	return simd::toFloat(simd::dot2(v1, simd::yxzw(simd::neg1(v2))));
}

template<>
inline float distance(const Vector2<float>& v1, const Vector2<float>& v2) noexcept
{
	auto v = simd::sub4(v2, v1);
	return simd::toFloat(simd::sqrt1(simd::dot2(v, v)));
}

template<>
inline float distanceSquared(const Vector2<float>& v1, const Vector2<float>& v2) noexcept
{
	auto v = simd::sub4(v2, v1);
	return simd::toFloat(simd::dot2(v, v));
}

template<>
inline float length(const Vector2<float>& v) noexcept
{
	return simd::toFloat(simd::sqrt1(simd::dot2(v, v)));
}

template<>
inline float lengthSquared(const Vector2<float>& v) noexcept
{
	return simd::toFloat(simd::dot2(v, v));
}

template<>
inline Vector2<float> minimum(const Vector2<float>& v1, const Vector2<float>& v2)
{
	return Vector2<float>(simd::min4(v1, v2));
}

template<>
inline Vector2<float> maximum(const Vector2<float>& v1, const Vector2<float>& v2)
{
	return Vector2<float>(simd::max4(v1, v2));
}

template<>
inline Vector2<float> lerp(const Vector2<float>& v1, const Vector2<float>& v2, float t) noexcept
{
	return Vector2<float>(simd::mulAdd4(simd::set4(t), simd::sub4(v2, v1), v1));
}

template<>
inline Vector2<float> slerp(const Vector2<float>& v1, const Vector2<float>& v2, float t)
{
	float dp = simd::toFloat(simd::dot2(v1, v2));
	if ((1.f - dp) < Constants<float>::TOLERANCE)
	{
		Vector2<float> c(simd::mulAdd4(simd::set4(t), simd::sub4(v2, v1), v1));
		c.normalize();
		return c;
	}

	dp = std::clamp(dp, -1.f, 1.f);
	float theta = std::acos(dp)*t;
	float st = std::sin(theta);
	float ct = std::cos(theta);
	Vector2<float> c(simd::sub4(v2, simd::mul4(v1, simd::set4(dp))));
	c.normalize();
	return Vector2<float>(simd::mulAdd4(v1, simd::set4(ct), simd::mul4(c, simd::set4(st))));
}

template<>
inline Vector2<float> perpendicular(const Vector2<float>& v) noexcept
{
#if MATHEMATICS_SIMD_EXPAND_LAST
	return Vector2<float>(simd::neg1(simd::yxxx(v)));
#else
	return Vector2<float>(simd::neg1(simd::yxzw(v)));
#endif
}

#endif /* SIMD_HAS_FLOAT4 */

template<typename T>
	requires std::floating_point<T>
inline Vector2<T> normalize(const Vector2<T>& v) noexcept
{
	Vector2<T> u(v);
	u.normalize();
	return u;
}

template<typename T>
#if SIMD_HAS_FLOAT4
	requires (std::floating_point<T> && !std::same_as<T, float>)
#else
	requires std::floating_point<T>
#endif
inline Vector2<T> normalize(Vector2<T>&& v) noexcept
{
	v.normalize();
	return v;
}

template<typename T>
	requires std::floating_point<T>
inline Vector2<T> project(const Vector2<T>& v1, const Vector2<T>& v2) noexcept
{
	T m = v2.getMagnitudeSquared();
	return (m > T(0)) ? (dot(v1, v2)/m)*v2 : Vector4<T>::ZERO;
}

template<typename T>
	requires std::floating_point<T>
inline T angle(const Vector2<T>& v1, const Vector2<T>& v2)
{
#if MATHEMATICS_FAST_NORMALIZE
	Vector2<T> u(v1);
	Vector2<T> v(v2);
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

using Vector2f = templates::Vector2<float>;
using Vector2d = templates::Vector2<double>;

#if MATHEMATICS_DOUBLE
using Vector2 = templates::Vector2<double>;
using Vector2Arg = templates::Vector2<double>::ConstArg;
using Vector2Result = templates::Vector2<double>::ConstResult;
#else
using Vector2 = templates::Vector2<float>;
using Vector2Arg = templates::Vector2<float>::ConstArg;
using Vector2Result = templates::Vector2<float>::ConstResult;
#endif

using IntVector2 = templates::Vector2<int>;
using IntVector2Arg = templates::Vector2<int>::ConstArg;
using IntVector2Result = templates::Vector2<int>::ConstResult;

} // namespace core::mathematics

namespace std {

template<size_t I, typename T>
struct tuple_element;

template<size_t I, typename T>
struct tuple_element<I, ::core::mathematics::templates::Vector2<T>>
{
	using type = T;
};

template<typename T>
struct tuple_size;

template<typename T>
struct tuple_size<::core::mathematics::templates::Vector2<T>> : integral_constant<size_t, 2> 
{
};

template<typename T>
struct hash;

template<typename T>
struct hash<::core::mathematics::templates::Vector2<T>>
{
	std::size_t operator()(const ::core::mathematics::templates::Vector2<T>& v) const noexcept
	{
		std::hash<T> hasher;
		std::size_t seed = hasher(v.x) + 0x9e3779b9;
		seed ^= hasher(v.y) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
		return seed;
	}
};

} // namespace std

#include "Matrix2.hpp"

namespace core::mathematics::templates {

template<typename T>
inline Vector2<T>& Vector2<T>::operator*=(const Matrix2<T>& m)
{
	set(x*m.m00 + y*m.m10, x*m.m01 + y*m.m11);
	return *this;
}

template<typename T>
	requires std::floating_point<T>
inline Vector2<T> operator*(const Vector2<T>& v, const Matrix2<T>& m) noexcept
{
	return Vector2<T>(v.x*m.m00 + v.y*m.m10, v.x*m.m01 + v.y*m.m11);
}

//template<typename T>
//	requires std::floating_point<T>
//inline Vector2<T> operator*(const Matrix2<T>& m, const Vector2<T>& v) noexcept; // valid for column vectors only

template<typename T>
inline Vector2<T>& Vector2<T>::transform(const Matrix2<T>& m)
{
	*this *= m;
	return *this;
}

template<typename T>
	requires std::floating_point<T>
inline Vector2<T> transform(const Vector2<T>& v, const Matrix2<T>& m) noexcept
{
	return v*m;
}

#if SIMD_HAS_FLOAT4

inline Vector2<float>& Vector2<float>::operator*=(const Matrix2<float>& m)
{
	auto t = simd::mul4(simd::xxyy(xy), simd::pack2x2(m.row0, m.row1));
	t = simd::add4(t, simd::zwzw(t));
#if MATHEMATICS_SIMD_EXPAND_LAST
	xy = simd::xyyy(t);
#else
	xy = simd::cutoff2(t);
#endif
	return *this;
}

template<>
inline Vector2<float> operator*(const Vector2<float>& v, const Matrix2<float>& m) noexcept
{
	auto t = simd::mul4(simd::xxyy(v), simd::pack2x2(m.row0, m.row1));
	t = simd::add4(t, simd::zwzw(t));
#if MATHEMATICS_SIMD_EXPAND_LAST
	return Vector2<float>(simd::xyyy(t));
#else
	return Vector2<float>(simd::cutoff2(t));
#endif
}

//template<>
//inline Vector2<float> operator*(const Matrix2<float>& m, const Vector2<float>& v) noexcept; // valid for column vectors only

inline Vector2<float>& Vector2<float>::transform(const Matrix2<float>& m)
{
	*this *= m;
	return *this;
}

template<>
inline Vector2<float> transform(const Vector2<float>& v, const Matrix2<float>& m) noexcept
{
	return v*m;
}

#endif /* SIMD_HAS_FLOAT4 */

} // namespace core::mathematics::templates
