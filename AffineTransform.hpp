/*
 *	Name: AffineTransform
 *	Author: Pawel Mrochen
 */

#pragma once

#include <istream>
#include <ostream>
#include <limits>
#include <type_traits>
#include <concepts>
#include <algorithm>
#include <tuple>
#include <cstddef>
#include <cmath>
#include <Simd/Intrinsics.hpp>
#include "Constants.hpp"
#include "Axis.hpp"
#include "Scalar.hpp"
#include "Vector3.hpp"
#include "Vector4.hpp"
#include "Matrix2.hpp"
#include "Matrix3.hpp"

namespace core::mathematics {

//struct Orthogonal
//{
//};

namespace templates {

template<typename T>
	requires std::floating_point<T>
struct Quaternion;

template<typename T>
	requires std::floating_point<T>
struct YawPitchRoll;

template<typename T>
	requires std::floating_point<T>
struct Euler;

template<typename T>
	requires std::floating_point<T>
struct Plane;

template<typename T>
	requires std::floating_point<T>
struct AffineTransform
{
	using Real = T;
	using ComponentType = T;
	using ConstArg = const AffineTransform&;
	using ConstResult = const AffineTransform&;

	//static constexpr int NUM_COMPONENTS = 12;

	constexpr AffineTransform() noexcept : m00(), m01(), m02(), m10(), m11(), m12(), m20(), m21(), m22(), m30(), m31(), m32() {}
	explicit AffineTransform(Uninitialized) noexcept {}
	explicit AffineTransform(Identity) noexcept : m00(1), m01(), m02(), m10(), m11(1), m12(), m20(), m21(), m22(1), m30(), m31(), m32() {}
	constexpr AffineTransform(T m00, T m01, T m02, T m10, T m11, T m12, T m20, T m21, T m22, T m30, T m31, T m32) noexcept;
	constexpr AffineTransform(const Vector3<T>& row0, const Vector3<T>& row1, const Vector3<T>& row2, const Vector3<T>& row3) noexcept;
	constexpr explicit AffineTransform(const std::tuple<Vector3<T>, Vector3<T>, Vector3<T>>& t) noexcept;
	constexpr AffineTransform(const Matrix2<T>& m) noexcept;
	constexpr AffineTransform(const Matrix3<T>& m) noexcept;
	constexpr AffineTransform(const Matrix3<T>& r, const Vector3<T>& t) noexcept;
	AffineTransform(const Vector3<T>& forward, const Vector3<T>& origin) noexcept;
	AffineTransform(const Vector3<T>& up, const Vector3<T>& forward, const Vector3<T>& origin) noexcept;
	explicit AffineTransform(const T* m) noexcept;

	//explicit operator T*() noexcept { return &m00; }
	//explicit operator const T*() const noexcept { return &m00; }
	Vector3<T>& operator[](int i) noexcept { return reinterpret_cast<Vector3<T>*>(&m00)[i]; }
	const Vector3<T>& operator[](int i) const noexcept { return reinterpret_cast<const Vector3<T>*>(&m00)[i]; }

	bool operator==(const AffineTransform& m) const noexcept;
	bool operator!=(const AffineTransform& m) const noexcept { return !(*this == m); }

	template<typename A> void load(A& ar) { ar(m00, m01, m02, m10, m11, m12, m20, m21, m22, m30, m31, m32); }
	template<typename A> void save(A& ar) const { ar(m00, m01, m02, m10, m11, m12, m20, m21, m22, m30, m31, m32); }

	// #TODO rename get...() to make...()
	static AffineTransform getTranslation(const Vector3<T>& v) noexcept { return AffineTransform(Uninitialized()).setTranslation(v); }
	static AffineTransform getScaling(const Vector3<T>& v) noexcept { return AffineTransform(Uninitialized()).setScaling(v); }
	static AffineTransform getScaling(T f) noexcept { return AffineTransform(Uninitialized()).setScaling(f); }
	static AffineTransform getRotation(Axis axis, T angle) noexcept { return AffineTransform(Uninitialized()).setRotation(axis, angle); }
	static AffineTransform getRotation(const Vector3<T>& axis, T angle) noexcept { return AffineTransform(Uninitialized()).setRotation(axis, angle); }
	static AffineTransform getRotation(const YawPitchRoll<T>& r) noexcept { return AffineTransform(Uninitialized()).setRotation(r); }
	static AffineTransform getRotation(const Euler<T>& e) noexcept { return AffineTransform(Uninitialized()).setRotation(e); }
	static AffineTransform getRotation(const Quaternion<T>& q) noexcept { return AffineTransform(Uninitialized()).setRotation(q); }
	static AffineTransform getRotationTranslation(Axis axis, T angle, const Vector3<T>& t) noexcept { return AffineTransform(Uninitialized()).setRotationTranslation(axis, angle, t); }
	static AffineTransform getRotationTranslation(const Vector3<T>& axis, T angle, const Vector3<T>& t) noexcept { return AffineTransform(Uninitialized()).setRotationTranslation(axis, angle, t); }
	static AffineTransform getRotationTranslation(const YawPitchRoll<T>& r, const Vector3<T>& t) noexcept { return AffineTransform(Uninitialized()).setRotationTranslation(r, t); }
	static AffineTransform getRotationTranslation(const Euler<T>& e, const Vector3<T>& t) noexcept { return AffineTransform(Uninitialized()).setRotationTranslation(e, t); }
	static AffineTransform getRotationTranslation(const Quaternion<T>& q, const Vector3<T>& t) noexcept { return AffineTransform(Uninitialized()).setRotationTranslation(q, t); }
	static AffineTransform getScalingRotation(const Vector3<T>& s, Axis axis, T angle) noexcept { return AffineTransform(Uninitialized()).setScalingRotation(s, axis, angle); }
	static AffineTransform getScalingRotation(const Vector3<T>& s, const Vector3<T>& axis, T angle) noexcept { return AffineTransform(Uninitialized()).setScalingRotation(s, axis, angle); }
	static AffineTransform getScalingRotation(const Vector3<T>& s, const YawPitchRoll<T>& r) noexcept { return AffineTransform(Uninitialized()).setScalingRotation(s, r); }
	static AffineTransform getScalingRotation(const Vector3<T>& s, const Euler<T>& e) noexcept { return AffineTransform(Uninitialized()).setScalingRotation(s, e); }
	static AffineTransform getScalingRotation(const Vector3<T>& s, const Quaternion<T>& q) noexcept { return AffineTransform(Uninitialized()).setScalingRotation(s, q); }
	static AffineTransform getScalingTranslation(const Vector3<T>& s, const Vector3<T>& t) noexcept { return AffineTransform(Uninitialized()).setScalingTranslation(s, t); }
	static AffineTransform getScalingRotationTranslation(const Vector3<T>& s, Axis axis, T angle, const Vector3<T>& t) noexcept { return AffineTransform(Uninitialized()).setScalingRotationTranslation(s, axis, angle, t); }
	static AffineTransform getScalingRotationTranslation(const Vector3<T>& s, const Vector3<T>& axis, T angle, const Vector3<T>& t) noexcept { return AffineTransform(Uninitialized()).setScalingRotationTranslation(s, axis, angle, t); }
	static AffineTransform getScalingRotationTranslation(const Vector3<T>& s, const YawPitchRoll<T>& r, const Vector3<T>& t) noexcept { return AffineTransform(Uninitialized()).setScalingRotationTranslation(s, r, t); }
	static AffineTransform getScalingRotationTranslation(const Vector3<T>& s, const Euler<T>& e, const Vector3<T>& t) noexcept { return AffineTransform(Uninitialized()).setScalingRotationTranslation(s, e, t); }
	static AffineTransform getScalingRotationTranslation(const Vector3<T>& s, const Quaternion<T>& q, const Vector3<T>& t) noexcept { return AffineTransform(Uninitialized()).setScalingRotationTranslation(s, q, t); }
	static AffineTransform getShearing(T xy, T xz, T yx, T yz, T zx, T zy) noexcept { return AffineTransform(Uninitialized()).setShearing(xy, xz, yx, yz, zx, zy); }
	static AffineTransform getReflection(const Plane<T>& plane) noexcept { return AffineTransform(Uninitialized()).setReflection(plane); }

	//bool isZero() const noexcept;
	//bool isApproxZero() const noexcept;
	bool isIdentity() const noexcept;
	bool isApproxIdentity() const noexcept;
	bool isApproxEqual(const AffineTransform& m) const noexcept;
	bool isApproxEqual(const AffineTransform& m, T tolerance) const noexcept;
	bool isApproxOrthogonal() const noexcept { return getBasis().isApproxOrthogonal(); }
	bool hasApproxUniformScaling() const noexcept { return getBasis().hasApproxUniformScaling(); }
	bool isFinite() const noexcept;
	const Vector3<T>& getRow(int i) const noexcept { return reinterpret_cast<const Vector3<T>*>(&m00)[i]; }
	Vector4<T> getColumn(int i) const noexcept { return Vector4<T>((&m00)[i], (&m10)[i], (&m20)[i], (&m30)[i]); }
	Vector3<T> getScaleComponent() const noexcept { return getBasis().getScaleComponent(); }
	//Vector3<T>& getOrigin() noexcept { return *reinterpret_cast<Vector3<T>*>(&m30); }
	const Vector3<T>& getOrigin() const noexcept { return *reinterpret_cast<const Vector3<T>*>(&m30); }
	void setOrigin(const Vector3<T>& t) noexcept { *reinterpret_cast<Vector3<T>*>(&m30) = t; }
	//Matrix3<T>& getBasis() noexcept { return reinterpret_cast<Matrix3<T>&>(*this); }
	const Matrix3<T>& getBasis() const noexcept { return reinterpret_cast<const Matrix3<T>&>(*this); }
	void setBasis(const Matrix3<T>& m) noexcept { reinterpret_cast<Matrix3<T>&>(*this) = m; }
	AffineTransform& setZero/*zero*/() noexcept;
	AffineTransform& setIdentity/*makeIdentity*/() noexcept;
	AffineTransform& set(const Matrix3<T>& r, const Vector3<T>& t) noexcept;
	AffineTransform& set(const Vector3<T>& row0, const Vector3<T>& row1, const Vector3<T>& row2, const Vector3<T>& row3) noexcept;
	AffineTransform& set(T m00, T m01, T m02, T m10, T m11, T m12, T m20, T m21, T m22, T m30, T m31, T m32) noexcept;
	AffineTransform& setTranslation/*makeTranslation*/(const Vector3<T>& v) noexcept;
	AffineTransform& setScaling/*makeScaling*/(const Vector3<T>& v) noexcept;
	AffineTransform& setScaling/*makeScaling*/(T f) noexcept;
	AffineTransform& setRotation(Axis axis, T angle) noexcept;
	AffineTransform& setRotation(const Vector3<T>& axis, T angle) noexcept;
	AffineTransform& setRotation(const YawPitchRoll<T>& r) noexcept;
	AffineTransform& setRotation(const Euler<T>& e) noexcept;
	AffineTransform& setRotation(const Quaternion<T>& q) noexcept;
	AffineTransform& setRotationTranslation(Axis axis, T angle, const Vector3<T>& t) noexcept { return setRotation(axis, angle).setOrigin(t); }
	AffineTransform& setRotationTranslation(const Vector3<T>& axis, T angle, const Vector3<T>& t) noexcept { return setRotation(axis, angle).setOrigin(t); }
	AffineTransform& setRotationTranslation(const YawPitchRoll<T>& r, const Vector3<T>& t) noexcept { return setRotation(r).setOrigin(t); }
	AffineTransform& setRotationTranslation(const Euler<T>& e, const Vector3<T>& t) noexcept { return setRotation(e).setOrigin(t); }
	AffineTransform& setRotationTranslation(const Quaternion<T>& q, const Vector3<T>& t) noexcept { return setRotation(q).setOrigin(t); }
	AffineTransform& setScalingRotation(const Vector3<T>& s, Axis axis, T angle) noexcept { return setRotation(axis, angle).preScale(s); }
	AffineTransform& setScalingRotation(const Vector3<T>& s, const Vector3<T>& axis, T angle) noexcept { return setRotation(axis, angle).preScale(s); }
	AffineTransform& setScalingRotation(const Vector3<T>& s, const YawPitchRoll<T>& r) noexcept { return setRotation(r).preScale(s); }
	AffineTransform& setScalingRotation(const Vector3<T>& s, const Euler<T>& e) noexcept { return setRotation(e).preScale(s); }
	AffineTransform& setScalingRotation(const Vector3<T>& s, const Quaternion<T>& q) noexcept { return setRotation(q).preScale(s); }
	AffineTransform& setScalingTranslation(const Vector3<T>& s, const Vector3<T>& t) noexcept;
	AffineTransform& setScalingRotationTranslation(const Vector3<T>& s, Axis axis, T angle, const Vector3<T>& t) noexcept { return setRotation(axis, angle).preScale(s).setOrigin(t); }
	AffineTransform& setScalingRotationTranslation(const Vector3<T>& s, const Vector3<T>& axis, T angle, const Vector3<T>& t) noexcept { return setRotation(axis, angle).preScale(s).setOrigin(t); }
	AffineTransform& setScalingRotationTranslation(const Vector3<T>& s, const YawPitchRoll<T>& r, const Vector3<T>& t) noexcept { return setRotation(r).preScale(s).setOrigin(t); }
	AffineTransform& setScalingRotationTranslation(const Vector3<T>& s, const Euler<T>& e, const Vector3<T>& t) noexcept { return setRotation(e).preScale(s).setOrigin(t); }
	AffineTransform& setScalingRotationTranslation(const Vector3<T>& s, const Quaternion<T>& q, const Vector3<T>& t) noexcept { return setRotation(q).preScale(s).setOrigin(t); }
	AffineTransform& setShearing(T xy, T xz, T yx, T yz, T zx, T zy) noexcept;
	AffineTransform& setReflection(const Plane<T>& plane) noexcept;
	AffineTransform& setInverse/*inverseOf*/(const AffineTransform& m) noexcept;
	AffineTransform& setInverseOrthogonal/*inverseOfOrthogonal*/(const AffineTransform& m) noexcept;
	//template<typename U = void> AffineTransform& setInverse(const AffineTransform& m) noexcept;
	AffineTransform& preConcatenate(const AffineTransform& m) noexcept;
	AffineTransform& concatenate(const AffineTransform& m) noexcept;
	AffineTransform& translate(const Vector3<T>& v) noexcept { m30 += v.x; m31 += v.y; m32 += v.z; return *this; }
	AffineTransform& preScale(const Vector3<T>& v) noexcept;
	AffineTransform& preScale(T f) noexcept;
	AffineTransform& scale(const Vector3<T>& v) noexcept;
	AffineTransform& scale(T f) noexcept;
	//AffineTransform& rotate(Axis axis, T angle) noexcept { concatenate(AffineTransform(Uninitialized()).setRotation(axis, angle)); return *this; }
	//AffineTransform& rotate(const Vector3<T>& axis, T angle) noexcept { concatenate(AffineTransform(Uninitialized()).setRotation(axis, angle)); return *this; }
	//AffineTransform& rotate(const YawPitchRoll<T>& r) noexcept { concatenate(AffineTransform(Uninitialized()).setRotation(r)); return *this; }
	//AffineTransform& rotate(const Euler<T>& e) noexcept { concatenate(AffineTransform(Uninitialized()).setRotation(e)); return *this; }
	//AffineTransform& rotate(const Quaternion<T>& q) noexcept { concatenate(AffineTransform(Uninitialized()).setRotation(q)); return *this; }
	//AffineTransform& shear(T xy, T xz, T yx, T yz, T zx, T zy) noexcept { concatenate(AffineTransform(Uninitialized()).setShearing(xy, xz, yx, yz, zx, zy)); return *this; }
	AffineTransform& invert() noexcept;
	AffineTransform& invertOrthogonal() noexcept;
	//template<typename U = void> AffineTransform& invert() noexcept;
	AffineTransform& orthonormalize() noexcept { reinterpret_cast<Matrix3<T>&>(*this).orthonormalize(); return *this; }

	//static const AffineTransform& getZero() noexcept { return ZERO; }
	//static const AffineTransform& getIdentity() noexcept { return IDENTITY; }

	static const AffineTransform ZERO;
	static const AffineTransform IDENTITY;

	T m00, m01, m02;
	T m10, m11, m12;
	T m20, m21, m22;
	T m30, m31, m32;
};

template<typename T> const AffineTransform<T> AffineTransform<T>::ZERO{};
template<typename T> const AffineTransform<T> AffineTransform<T>::IDENTITY{ T(1), T(0), T(0), T(0), T(1), T(0), T(0), T(0), T(1), T(0), T(0), T(0) };

#if SIMD_HAS_FLOAT4

template<>
struct Quaternion<float>;

template<>
struct YawPitchRoll<float>;

template<>
struct Euler<float>;

template<>
struct Plane<float>;

template<>
struct AffineTransform<float>
{
	using Real = float;
	using ComponentType = float;
	using ConstArg = const AffineTransform;
	using ConstResult = const AffineTransform;

	//static constexpr int NUM_COMPONENTS = 12;

	/*constexpr*/ AffineTransform() noexcept;
	explicit AffineTransform(Uninitialized) noexcept {}
	explicit AffineTransform(Identity) noexcept;
	/*constexpr*/ AffineTransform(float m00, float m01, float m02, float m10, float m11, float m12, float m20, float m21, float m22, float m30, float m31, float m32) noexcept;
	/*constexpr*/ AffineTransform(const Vector3<float>& row0, const Vector3<float>& row1, const Vector3<float>& row2, const Vector3<float>& row3) noexcept;
	/*constexpr*/ explicit AffineTransform(const std::tuple<Vector3<float>, Vector3<float>, Vector3<float>>& t) noexcept;
	/*constexpr*/ AffineTransform(const Matrix2<float>& m) noexcept;
	/*constexpr*/ AffineTransform(const Matrix3<float>& m) noexcept;
	/*constexpr*/ AffineTransform(const Matrix3<float>& r, const Vector3<float>& t) noexcept;
	AffineTransform(const Vector3<float>& forward, const Vector3<float>& origin) noexcept;
	AffineTransform(const Vector3<float>& up, const Vector3<float>& forward, const Vector3<float>& origin) noexcept;
	explicit AffineTransform(const float* m) noexcept;
	explicit AffineTransform(const simd::float4* m) noexcept : row0(m[0]), row1(m[1]), row2(m[2]), row3(m[3]) {}
	AffineTransform(simd::float4 row0, simd::float4 row1, simd::float4 row2, simd::float4 row3) noexcept;
	explicit AffineTransform(const std::tuple<simd::float4, simd::float4, simd::float4>& t) noexcept;
	AffineTransform(const AffineTransform& m) noexcept : row0(m.row0), row1(m.row1), row2(m.row2), row3(m.row3) {}
	AffineTransform& operator=(const AffineTransform& m) noexcept { row0 = m.row0; row1 = m.row1; row2 = m.row2; row3 = m.row3; return *this; }

	//explicit operator float*() noexcept { return &m00; }
	//explicit operator const float*() const noexcept { return &m00; }
	explicit operator simd::float4*() noexcept { return &row0; }
	explicit operator const simd::float4*() const noexcept { return &row0; }
	Vector3<float>& operator[](int i) noexcept { return reinterpret_cast<Vector3<float>&>((&row0)[i]); }
	const Vector3<float>& operator[](int i) const noexcept { return reinterpret_cast<const Vector3<float>&>((&row0)[i]); }

	bool operator==(const AffineTransform& m) const noexcept;
	bool operator!=(const AffineTransform& m) const noexcept { return !(*this == m); }

	template<typename A> void load(A& ar);
	template<typename A> void save(A& ar) const { ar(m00, m01, m02, m10, m11, m12, m20, m21, m22, m30, m31, m32); }

	// #TODO rename get...() to make...()
	static AffineTransform getTranslation(const Vector3<float>& v) noexcept { return AffineTransform(Uninitialized()).setTranslation(v); }
	static AffineTransform getScaling(const Vector3<float>& v) noexcept { return AffineTransform(Uninitialized()).setScaling(v); }
	static AffineTransform getScaling(float f) noexcept { return AffineTransform(Uninitialized()).setScaling(f); }
	static AffineTransform getRotation(Axis axis, float angle) noexcept { return AffineTransform(Uninitialized()).setRotation(axis, angle); }
	static AffineTransform getRotation(const Vector3<float>& axis, float angle) noexcept { return AffineTransform(Uninitialized()).setRotation(axis, angle); }
	static AffineTransform getRotation(const YawPitchRoll<float>& r) noexcept { return AffineTransform(Uninitialized()).setRotation(r); }
	static AffineTransform getRotation(const Euler<float>& e) noexcept { return AffineTransform(Uninitialized()).setRotation(e); }
	static AffineTransform getRotation(const Quaternion<float>& q) noexcept { return AffineTransform(Uninitialized()).setRotation(q); }
	static AffineTransform getRotationTranslation(Axis axis, float angle, const Vector3<float>& t) noexcept { return AffineTransform(Uninitialized()).setRotationTranslation(axis, angle, t); }
	static AffineTransform getRotationTranslation(const Vector3<float>& axis, float angle, const Vector3<float>& t) noexcept { return AffineTransform(Uninitialized()).setRotationTranslation(axis, angle, t); }
	static AffineTransform getRotationTranslation(const YawPitchRoll<float>& r, const Vector3<float>& t) noexcept { return AffineTransform(Uninitialized()).setRotationTranslation(r, t); }
	static AffineTransform getRotationTranslation(const Euler<float>& e, const Vector3<float>& t) noexcept { return AffineTransform(Uninitialized()).setRotationTranslation(e, t); }
	static AffineTransform getRotationTranslation(const Quaternion<float>& q, const Vector3<float>& t) noexcept { return AffineTransform(Uninitialized()).setRotationTranslation(q, t); }
	static AffineTransform getScalingRotation(const Vector3<float>& s, Axis axis, float angle) noexcept { return AffineTransform(Uninitialized()).setScalingRotation(s, axis, angle); }
	static AffineTransform getScalingRotation(const Vector3<float>& s, const Vector3<float>& axis, float angle) noexcept { return AffineTransform(Uninitialized()).setScalingRotation(s, axis, angle); }
	static AffineTransform getScalingRotation(const Vector3<float>& s, const YawPitchRoll<float>& r) noexcept { return AffineTransform(Uninitialized()).setScalingRotation(s, r); }
	static AffineTransform getScalingRotation(const Vector3<float>& s, const Euler<float>& e) noexcept { return AffineTransform(Uninitialized()).setScalingRotation(s, e); }
	static AffineTransform getScalingRotation(const Vector3<float>& s, const Quaternion<float>& q) noexcept { return AffineTransform(Uninitialized()).setScalingRotation(s, q); }
	static AffineTransform getScalingTranslation(const Vector3<float>& s, const Vector3<float>& t) noexcept { return AffineTransform(Uninitialized()).setScalingTranslation(s, t); }
	static AffineTransform getScalingRotationTranslation(const Vector3<float>& s, Axis axis, float angle, const Vector3<float>& t) noexcept { return AffineTransform(Uninitialized()).setScalingRotationTranslation(s, axis, angle, t); }
	static AffineTransform getScalingRotationTranslation(const Vector3<float>& s, const Vector3<float>& axis, float angle, const Vector3<float>& t) noexcept { return AffineTransform(Uninitialized()).setScalingRotationTranslation(s, axis, angle, t); }
	static AffineTransform getScalingRotationTranslation(const Vector3<float>& s, const YawPitchRoll<float>& r, const Vector3<float>& t) noexcept { return AffineTransform(Uninitialized()).setScalingRotationTranslation(s, r, t); }
	static AffineTransform getScalingRotationTranslation(const Vector3<float>& s, const Euler<float>& e, const Vector3<float>& t) noexcept { return AffineTransform(Uninitialized()).setScalingRotationTranslation(s, e, t); }
	static AffineTransform getScalingRotationTranslation(const Vector3<float>& s, const Quaternion<float>& q, const Vector3<float>& t) noexcept { return AffineTransform(Uninitialized()).setScalingRotationTranslation(s, q, t); }
	static AffineTransform getShearing(float xy, float xz, float yx, float yz, float zx, float zy) noexcept { return AffineTransform(Uninitialized()).setShearing(xy, xz, yx, yz, zx, zy); }
	static AffineTransform getReflection(const Plane<float>& plane) noexcept { return AffineTransform(Uninitialized()).setReflection(plane); }

	//bool isZero() const noexcept;
	//bool isApproxZero() const noexcept;
	bool isIdentity() const noexcept;
	bool isApproxIdentity() const noexcept;
	bool isApproxEqual(const AffineTransform& m) const noexcept;
	bool isApproxEqual(const AffineTransform& m, float tolerance) const noexcept;
	bool isApproxOrthogonal() const noexcept { return getBasis().isApproxOrthogonal(); }
	bool hasApproxUniformScaling() const noexcept { return getBasis().hasApproxUniformScaling(); }
	bool isFinite() const noexcept;
	const Vector3<float>& getRow(int i) const noexcept { return reinterpret_cast<const Vector3<float>&>((&row0)[i]); }
	Vector4<float> getColumn(int i) const noexcept; // #TODO
	Vector3<float> getScaleComponent() const noexcept { return getBasis().getScaleComponent(); }
	//Vector3<float>& getOrigin() noexcept { return reinterpret_cast<Vector3<float>&>(row3); }
	const Vector3<float>& getOrigin() const noexcept { return reinterpret_cast<const Vector3<float>&>(row3); }
	void setOrigin(const Vector3<float>& t) noexcept { row3 = t; }
	//Matrix3<float>& getBasis() noexcept { return reinterpret_cast<Matrix3<float>&>(*this); }
	const Matrix3<float>& getBasis() const noexcept { return reinterpret_cast<const Matrix3<float>&>(*this); }
	void setBasis(const Matrix3<float>& m) noexcept { row0 = m.row0; row1 = m.row1; row2 = m.row2; }
	AffineTransform& setZero/*zero*/() noexcept;
	AffineTransform& setIdentity/*makeIdentity*/() noexcept;
	AffineTransform& set(const Matrix3<float>& r, const Vector3<float>& t) noexcept;
	AffineTransform& set(const Vector3<float>& row0, const Vector3<float>& row1, const Vector3<float>& row2, const Vector3<float>& row3) noexcept;
	AffineTransform& set(simd::float4 row0, simd::float4 row1, simd::float4 row2, simd::float4 row3) noexcept;
	AffineTransform& set(float m00, float m01, float m02, float m10, float m11, float m12, float m20, float m21, float m22, float m30, float m31, float m32) noexcept;
	AffineTransform& setTranslation/*makeTranslation*/(const Vector3<float>& v) noexcept;
	AffineTransform& setScaling/*makeScaling*/(const Vector3<float>& v) noexcept;
	AffineTransform& setScaling/*makeScaling*/(float f) noexcept;
	AffineTransform& setRotation(Axis axis, float angle) noexcept;
	AffineTransform& setRotation(const Vector3<float>& axis, float angle) noexcept;
	AffineTransform& setRotation(const YawPitchRoll<float>& r) noexcept;
	AffineTransform& setRotation(const Euler<float>& e) noexcept;
	AffineTransform& setRotation(const Quaternion<float>& q) noexcept;
	AffineTransform& setRotationTranslation(Axis axis, float angle, const Vector3<float>& t) noexcept { return setRotation(axis, angle).setOrigin(t); }
	AffineTransform& setRotationTranslation(const Vector3<float>& axis, float angle, const Vector3<float>& t) noexcept { return setRotation(axis, angle).setOrigin(t); }
	AffineTransform& setRotationTranslation(const YawPitchRoll<float>& r, const Vector3<float>& t) noexcept { return setRotation(r).setOrigin(t); }
	AffineTransform& setRotationTranslation(const Euler<float>& e, const Vector3<float>& t) noexcept { return setRotation(e).setOrigin(t); }
	AffineTransform& setRotationTranslation(const Quaternion<float>& q, const Vector3<float>& t) noexcept { return setRotation(q).setOrigin(t); }
	AffineTransform& setScalingRotation(const Vector3<float>& s, Axis axis, float angle) noexcept { return setRotation(axis, angle).preScale(s); }
	AffineTransform& setScalingRotation(const Vector3<float>& s, const Vector3<float>& axis, float angle) noexcept { return setRotation(axis, angle).preScale(s); }
	AffineTransform& setScalingRotation(const Vector3<float>& s, const YawPitchRoll<float>& r) noexcept { return setRotation(r).preScale(s); }
	AffineTransform& setScalingRotation(const Vector3<float>& s, const Euler<float>& e) noexcept { return setRotation(e).preScale(s); }
	AffineTransform& setScalingRotation(const Vector3<float>& s, const Quaternion<float>& q) noexcept { return setRotation(q).preScale(s); }
	AffineTransform& setScalingTranslation(const Vector3<float>& s, const Vector3<float>& t) noexcept;
	AffineTransform& setScalingRotationTranslation(const Vector3<float>& s, Axis axis, float angle, const Vector3<float>& t) noexcept { return setRotation(axis, angle).preScale(s).setOrigin(t); }
	AffineTransform& setScalingRotationTranslation(const Vector3<float>& s, const Vector3<float>& axis, float angle, const Vector3<float>& t) noexcept { return setRotation(axis, angle).preScale(s).setOrigin(t); }
	AffineTransform& setScalingRotationTranslation(const Vector3<float>& s, const YawPitchRoll<float>& r, const Vector3<float>& t) noexcept { return setRotation(r).preScale(s).setOrigin(t); }
	AffineTransform& setScalingRotationTranslation(const Vector3<float>& s, const Euler<float>& e, const Vector3<float>& t) noexcept { return setRotation(e).preScale(s).setOrigin(t); }
	AffineTransform& setScalingRotationTranslation(const Vector3<float>& s, const Quaternion<float>& q, const Vector3<float>& t) noexcept { return setRotation(q).preScale(s).setOrigin(t); }
	AffineTransform& setShearing(float xy, float xz, float yx, float yz, float zx, float zy) noexcept;
	AffineTransform& setReflection(const Plane<float>& plane) noexcept;
	AffineTransform& setInverse/*inverseOf*/(const AffineTransform& m) noexcept;
	AffineTransform& setInverseOrthogonal/*inverseOfOrthogonal*/(const AffineTransform& m) noexcept;
	//template<typename U = void> AffineTransform& setInverse(const AffineTransform& m) noexcept;
	AffineTransform& preConcatenate(const AffineTransform& m) noexcept;
	AffineTransform& concatenate(const AffineTransform& m) noexcept;
	AffineTransform& translate(const Vector3<float>& v) noexcept { row3 = simd::add4(row3, v); return *this; }
	AffineTransform& preScale(const Vector3<float>& v) noexcept;
	AffineTransform& preScale(float f) noexcept;
	AffineTransform& scale(const Vector3<float>& v) noexcept;
	AffineTransform& scale(float f) noexcept;
	//AffineTransform& rotate(Axis axis, float angle) noexcept { concatenate(AffineTransform(Uninitialized()).setRotation(axis, angle)); return *this; }
	//AffineTransform& rotate(const Vector3<float>& axis, float angle) noexcept { concatenate(AffineTransform(Uninitialized()).setRotation(axis, angle)); return *this; }
	//AffineTransform& rotate(const YawPitchRoll<float>& r) noexcept { concatenate(AffineTransform(Uninitialized()).setRotation(r)); return *this; }
	//AffineTransform& rotate(const Euler<float>& e) noexcept { concatenate(AffineTransform(Uninitialized()).setRotation(e)); return *this; }
	//AffineTransform& rotate(const Quaternion<float>& q) noexcept { concatenate(AffineTransform(Uninitialized()).setRotation(q)); return *this; }
	//AffineTransform& shear(float xy, float xz, float yx, float yz, float zx, float zy) noexcept { concatenate(AffineTransform(Uninitialized()).setShearing(xy, xz, yx, yz, zx, zy)); return *this; }
	AffineTransform& invert() noexcept;
	AffineTransform& invertOrthogonal() noexcept;
	//template<typename U = void> AffineTransform& invert() noexcept;
	AffineTransform& orthonormalize() noexcept { reinterpret_cast<Matrix3<float>&>(*this).orthonormalize(); return *this; }

	//static const AffineTransform& getZero() noexcept { return ZERO; }
	//static const AffineTransform& getIdentity() noexcept { return IDENTITY; }

	static const AffineTransform ZERO;
	static const AffineTransform IDENTITY;

	union
	{
		simd::float4 row0;
		struct { float m00, m01, m02; };
	};

	union
	{
		simd::float4 row1;
		struct { float m10, m11, m12; };
	};

	union
	{
		simd::float4 row2;
		struct { float m20, m21, m22; };
	};

	union
	{
		simd::float4 row3;
		struct { float m30, m31, m32; };
	};
};

const AffineTransform<float> AffineTransform<float>::ZERO{};
const AffineTransform<float> AffineTransform<float>::IDENTITY{ 1.f, 0.f, 0.f, 0.f, 1.f, 0.f, 0.f, 0.f, 1.f, 0.f, 0.f, 0.f };

#endif /* SIMD_HAS_FLOAT4 */

template<typename T>
inline AffineTransform<T>::AffineTransform(T m00, T m01, T m02, T m10, T m11, T m12, T m20, T m21, T m22, T m30, T m31, T m32) :
	m00(m00), m01(m01), m02(m02), m10(m10), m11(m11), m12(m12), m20(m20), m21(m21), m22(m22), m30(m30), m31(m31), m32(m32)
{
}

template<typename T>
inline AffineTransform<T>::AffineTransform(const Vector3<T>& row0, const Vector3<T>& row1, const Vector3<T>& row2, const Vector3<T>& row3) :
	m00(row0.x), m01(row0.y), m02(row0.z), m10(row1.x), m11(row1.y), m12(row1.z), m20(row2.x), m21(row2.y), m22(row2.z),
	m30(row3.x), m31(row3.y), m32(row3.z)
{
}

template<typename T>
inline AffineTransform<T>::AffineTransform(const std::tuple<Vector3<T>, Vector3<T>, Vector3<T>, Vector3<T>>& t) :
	AffineTransform(std::get<0>(t), std::get<1>(t), std::get<2>(t), std::get<3>(t))
{
}

template<typename T>
inline AffineTransform<T>::AffineTransform(const Matrix2<T>& m) :
	m00(m.m00), m01(m.m01), m02(), m10(m.m10), m11(m.m11), m12(), m20(), m21(), m22(1), m30(), m31(), m32()
{
}

template<typename T>
inline AffineTransform<T>::AffineTransform(const Matrix3<T>& m) :
	m00(m.m00), m01(m.m01), m02(m.m02), m10(m.m10), m11(m.m11), m12(m.m12), m20(m.m20), m21(m.m21), m22(m.m22), m30(), m31(), m32()
{
}

template<typename T>
inline AffineTransform<T>::AffineTransform(const Matrix3<T>& r, const Vector3<T>& t) :
	m00(r.m00), m01(r.m01), m02(r.m02), m10(r.m10), m11(r.m11), m12(r.m12), m20(r.m20), m21(r.m21), m22(r.m22), 
	m30(t.x), m31(t.y), m32(t.z)
{
}

template<typename T>
inline AffineTransform<T>::AffineTransform(const Vector3<T>& forward, const Vector3<T>& origin) /*:
	AffineTransform<T>(Matrix3<T>(forward), origin)*/
{
	T m = forward.getMagnitude();
	if (m > T(0))
	{
		Vector3<T> at = forward/m;
		if ((T(1) - std::fabs(at.y)) >= Constants<T>::TOLERANCE)
		{
			Vector3<T> right(at.z, T(0), -at.x);
			right.normalize();
			Vector3<T> up = cross(at, right);
			set(right, up, at, origin);
		}
		else
		{
			Vector3<T> up((at.y > T(0)) ? -Vector3<T>::UNIT_Z : Vector3<T>::UNIT_Z);
			Vector3<T> right = cross(up, at);
			right.normalize();
			up = cross(at, right);
			set(right, up, at, origin);
		}
	}
	else
	{
		set(Matrix3<T>::IDENTITY, origin);
	}
}

template<typename T>
inline AffineTransform<T>::AffineTransform(const Vector3<T>& up, const Vector3<T>& forward, const Vector3<T>& origin) /*:
	AffineTransform<T>(Matrix3<T>(up, forward), origin)*/
{
	T m = forward.getMagnitude();
	if ((m > T(0)) && (up.getMagnitude() > T(0)))
	{
		Vector3<T> at = forward/m;
		Vector3<T> right = cross(up, at);
		right.normalize();
		set(right, cross(at, right), at, origin);
	}
	else
	{
		set(Matrix3<T>::IDENTITY, origin);
	}
}

template<typename T>
inline AffineTransform<T>::AffineTransform(const T* m) :
	m00(m[0]), m01(m[1]), m02(m[2]), m10(m[3]), m11(m[4]), m12(m[5]), m20(m[6]), m21(m[7]), m22(m[8]), m30(m[9]), m31(m[10]), m32(m[11])
{
}

template<typename T>
inline bool AffineTransform<T>::operator==(const AffineTransform<T>& m) const
{
	return (m00 == m.m00) && (m01 == m.m01) && (m02 == m.m02) &&
		(m10 == m.m10) && (m11 == m.m11) && (m12 == m.m12) &&
		(m20 == m.m20) && (m21 == m.m21) && (m22 == m.m22) &&
		(m30 == m.m30) && (m31 == m.m31) && (m32 == m.m32);
}

template<typename C, typename T, typename U>
	requires std::floating_point<U>
inline std::basic_istream<C, T>& operator>>(std::basic_istream<C, T>& s, AffineTransform<U>& m)
{
	return s >> m.m00 >> std::ws >> m.m01 >> std::ws >> m.m02 >> std::ws >>
		m.m10 >> std::ws >> m.m11 >> std::ws >> m.m12 >> std::ws >>
		m.m20 >> std::ws >> m.m21 >> std::ws >> m.m22 >> std::ws >>
		m.m30 >> std::ws >> m.m31 >> std::ws >> m.m32;
}

template<typename C, typename T, typename U>
	requires std::floating_point<U>
inline std::basic_ostream<C, T>& operator<<(std::basic_ostream<C, T>& s, const AffineTransform<U>& m)
{
	constexpr C WS(0x20);
	return s << m.m00 << WS << m.m01 << WS << m.m02 << WS <<
		m.m10 << WS << m.m11 << WS << m.m12 << WS <<
		m.m20 << WS << m.m21 << WS << m.m22 << WS <<
		m.m30 << WS << m.m31 << WS << m.m32;
}

template<typename T>
inline bool AffineTransform<T>::isIdentity() const
{
	return (m00 == T(1)) && (m01 == T()) && (m02 == T()) &&
		(m10 == T()) && (m11 == T(1)) && (m12 == T()) &&
		(m20 == T()) && (m21 == T()) && (m22 == T(1)) &&
		(m30 == T()) && (m31 == T()) && (m32 == T());
}

template<typename T>
inline bool AffineTransform<T>::isApproxIdentity() const
{
	return (std::fabs(T(1) - m00) < Constants<T>::TOLERANCE) &&
		(std::fabs(m01) < Constants<T>::TOLERANCE) &&
		(std::fabs(m02) < Constants<T>::TOLERANCE) &&
		(std::fabs(m10) < Constants<T>::TOLERANCE) &&
		(std::fabs(T(1) - m11) < Constants<T>::TOLERANCE) &&
		(std::fabs(m12) < Constants<T>::TOLERANCE) &&
		(std::fabs(m20) < Constants<T>::TOLERANCE) &&
		(std::fabs(m21) < Constants<T>::TOLERANCE) &&
		(std::fabs(T(1) - m22) < Constants<T>::TOLERANCE) &&
		(std::fabs(m30) < Constants<T>::TOLERANCE) &&
		(std::fabs(m31) < Constants<T>::TOLERANCE) &&
		(std::fabs(m32) < Constants<T>::TOLERANCE);
}

template<typename T>
inline bool AffineTransform<T>::isApproxEqual(const AffineTransform<T>& m) const
{
	return (std::fabs(m.m00 - m00) < Constants<T>::TOLERANCE) &&
		(std::fabs(m.m01 - m01) < Constants<T>::TOLERANCE) &&
		(std::fabs(m.m02 - m02) < Constants<T>::TOLERANCE) &&
		(std::fabs(m.m10 - m10) < Constants<T>::TOLERANCE) &&
		(std::fabs(m.m11 - m11) < Constants<T>::TOLERANCE) &&
		(std::fabs(m.m12 - m12) < Constants<T>::TOLERANCE) &&
		(std::fabs(m.m20 - m20) < Constants<T>::TOLERANCE) &&
		(std::fabs(m.m21 - m21) < Constants<T>::TOLERANCE) &&
		(std::fabs(m.m22 - m22) < Constants<T>::TOLERANCE) &&
		(std::fabs(m.m30 - m30) < Constants<T>::TOLERANCE) &&
		(std::fabs(m.m31 - m31) < Constants<T>::TOLERANCE) &&
		(std::fabs(m.m32 - m32) < Constants<T>::TOLERANCE);
}

template<typename T>
inline bool AffineTransform<T>::isApproxEqual(const AffineTransform<T>& m, T tolerance) const
{
	return (std::fabs(m.m00 - m00) < tolerance) &&
		(std::fabs(m.m01 - m01) < tolerance) &&
		(std::fabs(m.m02 - m02) < tolerance) &&
		(std::fabs(m.m10 - m10) < tolerance) &&
		(std::fabs(m.m11 - m11) < tolerance) &&
		(std::fabs(m.m12 - m12) < tolerance) &&
		(std::fabs(m.m20 - m20) < tolerance) &&
		(std::fabs(m.m21 - m21) < tolerance) &&
		(std::fabs(m.m22 - m22) < tolerance) &&
		(std::fabs(m.m30 - m30) < tolerance) &&
		(std::fabs(m.m31 - m31) < tolerance) &&
		(std::fabs(m.m32 - m32) < tolerance);
}

template<typename T>
inline bool AffineTransform<T>::isFinite() const
{
	return std::isfinite(m00) && std::isfinite(m01) && std::isfinite(m02) &&
		std::isfinite(m10) && std::isfinite(m11) && std::isfinite(m12) &&
		std::isfinite(m20) && std::isfinite(m21) && std::isfinite(m22) &&
		std::isfinite(m30) && std::isfinite(m31) && std::isfinite(m32);
}

template<typename T>
inline AffineTransform<T>& AffineTransform<T>::setZero()
{
	m00 = T(); m01 = T(); m02 = T();
	m10 = T(); m11 = T(); m12 = T();
	m20 = T(); m21 = T(); m22 = T();
	m30 = T(); m31 = T(); m32 = T();
	return *this;
}

template<typename T>
inline AffineTransform<T>& AffineTransform<T>::setIdentity()
{
	m00 = T(1); m01 = T(); m02 = T();
	m10 = T(); m11 = T(1); m12 = T();
	m20 = T(); m21 = T(); m22 = T(1);
	m30 = T(); m31 = T(); m32 = T();
	return *this;
}

template<typename T>
inline AffineTransform<T>& AffineTransform<T>::set(const Matrix3<T>& r, const Vector3<T>& t)
{
	m00 = r.m00; m01 = r.m01; m02 = r.m02; 
	m10 = r.m10; m11 = r.m11; m12 = r.m12;
	m20 = r.m20; m21 = r.m21; m22 = r.m22; 
	m30 = t.x; m31 = t.y; m32 = t.z; 
	return *this;
}

template<typename T>
inline AffineTransform<T>& AffineTransform<T>::set(const Vector3<T>& row0, const Vector3<T>& row1, const Vector3<T>& row2, 
	const Vector3<T>& row3)
{
	m00 = row0.x; m01 = row0.y; m02 = row0.z; 
	m10 = row1.x; m11 = row1.y; m12 = row1.z; 
	m20 = row2.x; m21 = row2.y; m22 = row2.z;
	m30 = row3.x; m31 = row3.y; m32 = row3.z;
	return *this;
}

template<typename T>
inline AffineTransform<T>& AffineTransform<T>::set(T m00, T m01, T m02, T m10, T m11, T m12, T m20, T m21, T m22, T m30, T m31, T m32)
{
	this->m00 = m00; this->m01 = m01; this->m02 = m02; 
	this->m10 = m10; this->m11 = m11; this->m12 = m12; 
	this->m20 = m20; this->m21 = m21; this->m22 = m22; 
	this->m30 = m30; this->m31 = m31; this->m32 = m32; 
	return *this;
}

template<typename T>
inline AffineTransform<T>& AffineTransform<T>::setTranslation(const Vector3<T>& v)
{
	m00 = T(1); m01 = T(); m02 = T();
	m10 = T(); m11 = T(1); m12 = T();
	m20 = T(); m21 = T(); m22 = T(1);
	m30 = v.x; m31 = v.y; m32 = v.z;
	return *this;
}

template<typename T>
inline AffineTransform<T>& AffineTransform<T>::setScaling(const Vector3<T>& v)
{
	m00 = v.x; m01 = T(); m02 = T();
	m10 = T(); m11 = v.y; m12 = T();
	m20 = T(); m21 = T(); m22 = v.z;
	m30 = T(); m31 = T(); m32 = T();
	return *this;
}

template<typename T>
inline AffineTransform<T>& AffineTransform<T>::setScaling(T f)
{
	m00 = f; m01 = T(); m02 = T();
	m10 = T(); m11 = f; m12 = T();
	m20 = T(); m21 = T(); m22 = f;
	m30 = T(); m31 = T(); m32 = T();
	return *this;
}

template<typename T>
inline AffineTransform<T>& AffineTransform<T>::setRotation(Axis axis, T angle)
{
	if (angle != T(0))
	{
		T sinAngle = std::sin(angle);
		T cosAngle = std::cos(angle);

		switch (axis)
		{
			case Axis::X:
				return set(T(1), T(), T(), T(), cosAngle, sinAngle, T(), -sinAngle, cosAngle, T(), T(), T());
			case Axis::Y:
				return set(cosAngle, T(), -sinAngle, T(), T(1), T(), sinAngle, T(), cosAngle, T(), T(), T());
			case Axis::Z:
				return set(cosAngle, sinAngle, T(), -sinAngle, cosAngle, T(), T(), T(), T(1), T(), T(), T());
			default:
				return setZero();
		}
	}

	return setIdentity();
}

template<typename T>
inline AffineTransform<T>& AffineTransform<T>::setRotation(const Vector3<T>& axis, T angle)
{
	T m = axis.getMagnitude();
	if ((m > T(0)) && (angle != T(0)))
	{
		Vector3<T> abc = axis/m;
		T sinAngle = std::sin(angle);
		T cosAngle = std::cos(angle);
		T mCos = T(1) - cosAngle;
		T aa = abc.x*abc.x, bb = abc.y*abc.y, cc = abc.z*abc.z;
		T abm = abc.x*abc.y*mCos, acm = abc.x*abc.z*mCos, bcm = abc.y*abc.z*mCos;
		T aSin = abc.x*sinAngle, bSin = abc.y*sinAngle, cSin = abc.z*sinAngle;

		return set(aa*mCos + cosAngle, abm + cSin, acm - bSin,
			abm - cSin, bb*mCos + cosAngle, bcm + aSin,
			acm + bSin, bcm - aSin, cc*mCos + cosAngle,
			T(), T(), T());
	}

	return setIdentity();
}

template<typename T>
inline AffineTransform<T>& AffineTransform<T>::setScalingTranslation(const Vector3<T>& s, const Vector3<T>& t)
{
	m00 = s.x; m01 = T(); m02 = T();
	m10 = T(); m11 = s.y; m12 = T();
	m20 = T(); m21 = T(); m22 = s.z;
	m30 = t.x; m31 = t.y; m32 = t.z;
	return *this;
}

template<typename T>
inline AffineTransform<T>& AffineTransform<T>::setShearing(T xy, T xz, T yx, T yz, T zx, T zy)
{
	m00 = T(1); m01 = xy; m02 = xz;
	m10 = yx; m11 = T(1); m12 = yz;
	m20 = zx; m21 = zy; m22 = T(1);
	m30 = T(); m31 = T(); m32 = T();
	return *this;
}

template<typename T>
inline AffineTransform<T>& AffineTransform<T>::setInverse(const AffineTransform<T>& m)
{
	T a = m.m11*m.m22 - m.m12*m.m21;
	T b = m.m12*m.m20 - m.m10*m.m22;
	T c = m.m10*m.m21 - m.m11*m.m20;
	T d = m.m02*m.m21 - m.m01*m.m22;
	T e = m.m00*m.m22 - m.m02*m.m20;
	T f = m.m01*m.m20 - m.m00*m.m21;
	T g = m.m01*m.m12 - m.m02*m.m11;
	T h = m.m02*m.m10 - m.m00*m.m12;
	T i = m.m00*m.m11 - m.m01*m.m10;
	T n = T(1)/(m.m00*a + m.m01*b + m.m02*c);
	m00 = a*n; m01 = d*n; m02 = g*n;
	m10 = b*n; m11 = e*n; m12 = h*n;
	m20 = c*n; m21 = f*n; m22 = i*n;
	setOrigin(-(m.getOrigin()*getBasis()));
	return *this;
}

template<typename T>
inline AffineTransform<T>& AffineTransform<T>::setInverseOrthogonal(const AffineTransform<T>& m)
{
	m00 = m.m00; m01 = m.m10; m02 = m.m20;
	m10 = m.m01; m11 = m.m11; m12 = m.m21;
	m20 = m.m02; m21 = m.m12; m22 = m.m22;
	setOrigin(-(m.getOrigin()*getBasis()));
	return *this;
}

template<typename T>
inline AffineTransform<T>& AffineTransform<T>::preConcatenate(const AffineTransform<T>& m)
{
	set(m.m00*m00 + m.m01*m10 + m.m02*m20,
		m.m00*m01 + m.m01*m11 + m.m02*m21,
		m.m00*m02 + m.m01*m12 + m.m02*m22,
		m.m10*m00 + m.m11*m10 + m.m12*m20,
		m.m10*m01 + m.m11*m11 + m.m12*m21,
		m.m10*m02 + m.m11*m12 + m.m12*m22,
		m.m20*m00 + m.m21*m10 + m.m22*m20,
		m.m20*m01 + m.m21*m11 + m.m22*m21,
		m.m20*m02 + m.m21*m12 + m.m22*m22,
		m.m30*m00 + m.m31*m10 + m.m32*m20 + m30,
		m.m30*m01 + m.m31*m11 + m.m32*m21 + m31,
		m.m30*m02 + m.m31*m12 + m.m32*m22 + m32);
	return *this;
}

template<typename T>
inline AffineTransform<T>& AffineTransform<T>::concatenate(const AffineTransform<T>& m)
{
	set(m00*m.m00 + m01*m.m10 + m02*m.m20,
		m00*m.m01 + m01*m.m11 + m02*m.m21,
		m00*m.m02 + m01*m.m12 + m02*m.m22,
		m10*m.m00 + m11*m.m10 + m12*m.m20,
		m10*m.m01 + m11*m.m11 + m12*m.m21,
		m10*m.m02 + m11*m.m12 + m12*m.m22,
		m20*m.m00 + m21*m.m10 + m22*m.m20,
		m20*m.m01 + m21*m.m11 + m22*m.m21,
		m20*m.m02 + m21*m.m12 + m22*m.m22,
		m30*m.m00 + m31*m.m10 + m32*m.m20 + m.m30,
		m30*m.m01 + m31*m.m11 + m32*m.m21 + m.m31,
		m30*m.m02 + m31*m.m12 + m32*m.m22 + m.m32);
	return *this;
}

template<typename T>
inline AffineTransform<T>& AffineTransform<T>::preScale(const Vector3<T>& v)
{
	m00 *= v.x; m01 *= v.x; m02 *= v.x;
	m10 *= v.y; m11 *= v.y; m12 *= v.y;
	m20 *= v.z; m21 *= v.z; m22 *= v.z;
	return *this;
}

template<typename T>
inline AffineTransform<T>& AffineTransform<T>::preScale(T f)
{
	m00 *= f; m01 *= f; m02 *= f;
	m10 *= f; m11 *= f; m12 *= f;
	m20 *= f; m21 *= f; m22 *= f;
	return *this;
}

template<typename T>
inline AffineTransform<T>& AffineTransform<T>::scale(const Vector3<T>& v)
{
	m00 *= v.x; m01 *= v.y; m02 *= v.z;
	m10 *= v.x; m11 *= v.y; m12 *= v.z;
	m20 *= v.x; m21 *= v.y; m22 *= v.z;
	m30 *= v.x; m31 *= v.y; m32 *= v.z;
	return *this;
}

template<typename T>
inline AffineTransform<T>& AffineTransform<T>::scale(T f)
{
	m00 *= f; m01 *= f; m02 *= f;
	m10 *= f; m11 *= f; m12 *= f;
	m20 *= f; m21 *= f; m22 *= f;
	m30 *= f; m31 *= f; m32 *= f;
	return *this;
}

template<typename T>
inline AffineTransform<T>& AffineTransform<T>::invert()
{
	T a = m11*m22 - m12*m21;
	T b = m12*m20 - m10*m22;
	T c = m10*m21 - m11*m20;
	T d = m02*m21 - m01*m22;
	T e = m00*m22 - m02*m20;
	T f = m01*m20 - m00*m21;
	T g = m01*m12 - m02*m11;
	T h = m02*m10 - m00*m12;
	T i = m00*m11 - m01*m10;
	T n = T(1)/(m00*a + m01*b + m02*c);
	m00 = a*n; m01 = d*n; m02 = g*n;
	m10 = b*n; m11 = e*n; m12 = h*n;
	m20 = c*n; m21 = f*n; m22 = i*n;
	setOrigin(-(getOrigin()*getBasis()));
	return *this;
}

template<typename T>
inline AffineTransform<T>& AffineTransform<T>::invertOrthogonal()
{
	T t = m01; m01 = m10; m10 = t;
	t = m02; m02 = m20; m20 = t;
	t = m12; m12 = m21; m21 = t;
	setOrigin(-(getOrigin()*getBasis()));
	return *this;
}

#if SIMD_HAS_FLOAT4

inline AffineTransform<float>::AffineTransform() : 
	row0(simd::zero<simd::float4>()), 
	row1(simd::zero<simd::float4>()), 
	row2(simd::zero<simd::float4>()), 
	row3(simd::zero<simd::float4>()) 
{
}

inline AffineTransform<float>::AffineTransform(Identity) :
	row0(Vector3<float>::UNIT_X),
	row1(Vector3<float>::UNIT_Y),
	row2(Vector3<float>::UNIT_Z),
	row3(simd::zero<simd::float4>())
{
}

inline AffineTransform<float>::AffineTransform(float m00, float m01, float m02, float m10, float m11, float m12, 
	float m20, float m21, float m22, float m30, float m31, float m32) :
#if MATHEMATICS_SIMD_EXPAND_LAST
	row0(simd::set4(m00, m01, m02, m02)),
	row1(simd::set4(m10, m11, m12, m12)),
	row2(simd::set4(m20, m21, m22, m22)),
	row3(simd::set4(m30, m31, m32, m32))
#else
	row0(simd::set3(m00, m01, m02)),
	row1(simd::set3(m10, m11, m12)),
	row2(simd::set3(m20, m21, m22)),
	row3(simd::set3(m30, m31, m32))
#endif
{
}

inline AffineTransform<float>::AffineTransform(const Vector3<float>& row0, const Vector3<float>& row1, const Vector3<float>& row2,
	const Vector3<float>& row3) :
	row0(row0),
	row1(row1),
	row2(row2),
	row3(row3)
{
}

inline AffineTransform<float>::AffineTransform(const std::tuple<Vector3<float>, Vector3<float>, Vector3<float>, Vector3<float>>& t) :
	AffineTransform(std::get<0>(t), std::get<1>(t), std::get<2>(t), std::get<3>(t))
{
}

inline AffineTransform<float>::AffineTransform(const Matrix2<float>& m) :
	row0(simd::cutoff2(m.row0)),
	row1(simd::cutoff2(m.row1)),
	row2(Vector3<float>::UNIT_Z),
	row3(),
{
}

inline AffineTransform<float>::AffineTransform(const Matrix3<float>& m) :
	row0(m.row0),
	row1(m.row1),
	row2(m.row2),
	row3(),
{
}

inline AffineTransform<float>::AffineTransform(const Matrix3<float>& r, const Vector3<float>& t) :
	row0(r.row0),
	row1(r.row1),
	row2(r.row2),
	row3(t)
{
}

inline AffineTransform<float>::AffineTransform(const Vector3<float>& forward, const Vector3<float>& origin) /*:
	AffineTransform<float>(Matrix3<float>(forward), origin)*/
{
	float m = forward.getMagnitude();
	if (m > 0.f)
	{
		Vector3<float> at = forward/m;
		if ((1.f - std::fabs(at.y)) >= Constants<float>::TOLERANCE)
		{
			Vector3<float> right(at.z, 0.f, -at.x);
			right.normalize();
			Vector3<float> up = cross(at, right);
			set(right, up, at, origin);
		}
		else
		{
			Vector3<float> up((at.y > 0.f) ? -Vector3<float>::UNIT_Z : Vector3<float>::UNIT_Z);
			Vector3<float> right = cross(up, at);
			right.normalize();
			up = cross(at, right);
			set(right, up, at, origin);
		}
	}
	else
	{
		set(Matrix3<float>::IDENTITY, origin);
	}
}

inline AffineTransform<float>::AffineTransform(const Vector3<float>& up, const Vector3<float>& forward, const Vector3<float>& origin) /*:
	AffineTransform<float>(Matrix3<float>(up, forward), origin)*/
{
	float m = forward.getMagnitude();
	if ((m > 0.f) && (up.getMagnitude() > 0.f))
	{
		Vector3<float> at = forward/m;
		Vector3<float> right = cross(up, at);
		right.normalize();
		set(right, cross(at, right), at, origin);
	}
	else
	{
		set(Matrix3<float>::IDENTITY, origin);
	}
}

#if MATHEMATICS_SIMD_EXPAND_LAST
inline AffineTransform<float>::AffineTransform(const float* m) :
	AffineTransform(m[0], m[1], m[2], m[3], m[4], m[5], m[6], m[7], m[8], m[9], m[10], m[11])
{
}
#else
inline AffineTransform<float>::AffineTransform(const float* m)
{
	simd::unpack4x3(m, row0, row1, row2, row3);
}
#endif

inline AffineTransform<float>::AffineTransform(simd::float4 row0, simd::float4 row1, simd::float4 row2, simd::float4 row3) : 
	row0(row0), 
	row1(row1), 
	row2(row2), 
	row3(row3) 
{
}

inline AffineTransform<float>::AffineTransform(const std::tuple<simd::float4, simd::float4, simd::float4>& t) :
	row0(std::get<0>(t)), 
	row1(std::get<1>(t)), 
	row2(std::get<2>(t)), 
	row3(std::get<3>(t))
{
}

inline bool AffineTransform<float>::operator==(const AffineTransform<float>& m) const
{
	return simd::all3(simd::equal(row0, m.row0)) && 
		simd::all3(simd::equal(row1, m.row1)) && 
		simd::all3(simd::equal(row2, m.row2)) &&
		simd::all3(simd::equal(row3, m.row3));
}

template<typename C, typename T>
inline std::basic_istream<C, T>& operator>>(std::basic_istream<C, T>& s, AffineTransform<float>& m)
{
	float m00, m01, m02;
	float m10, m11, m12;
	float m20, m21, m22;
	float m30, m31, m32;
	s >> m00 >> std::ws >> m01 >> std::ws >> m02 >> std::ws >>
		m10 >> std::ws >> m11 >> std::ws >> m12 >> std::ws >>
		m20 >> std::ws >> m21 >> std::ws >> m22 >> std::ws >>
		m30 >> std::ws >> m31 >> std::ws >> m32;
	m.set(m00, m01, m02, m10, m11, m12, m20, m21, m22, m30, m31, m32);
	return s;
}

template<typename A>
inline void AffineTransform<float>::load(A& ar)
{
	float t00, t01, t02;
	float t10, t11, t12;
	float t20, t21, t22;
	float t30, t31, t32;
	ar(t00, t01, t02, t10, t11, t12, t20, t21, t22, t30, t31, t32);
	set(t00, t01, t02, t10, t11, t12, t20, t21, t22, t30, t31, t32);
}

inline bool AffineTransform<float>::isIdentity() const
{
	return simd::all3(simd::equal(row0, Vector3<float>::UNIT_X)) &&
		simd::all3(simd::equal(row1, Vector3<float>::UNIT_Y)) &&
		simd::all3(simd::equal(row2, Vector3<float>::UNIT_Z)) &&
		simd::all3(simd::equal(row3, simd::zero<simd::float4>()));
}

inline bool AffineTransform<float>::isApproxIdentity() const
{
	return simd::all3(simd::lessThan(simd::abs4(simd::sub4(row0, Vector3<float>::UNIT_X)), Vector3<float>::TOLERANCE)) &&
		simd::all3(simd::lessThan(simd::abs4(simd::sub4(row1, Vector3<float>::UNIT_Y)), Vector3<float>::TOLERANCE)) &&
		simd::all3(simd::lessThan(simd::abs4(simd::sub4(row2, Vector3<float>::UNIT_Z)), Vector3<float>::TOLERANCE)) &&
		simd::all3(simd::lessThan(simd::abs4(row3), Vector3<float>::TOLERANCE));
}

inline bool AffineTransform<float>::isApproxEqual(const AffineTransform& m) const
{
	return simd::all3(simd::lessThan(simd::abs4(simd::sub4(row0, m.row0)), Vector3<float>::TOLERANCE)) &&
		simd::all3(simd::lessThan(simd::abs4(simd::sub4(row1, m.row1)), Vector3<float>::TOLERANCE)) &&
		simd::all3(simd::lessThan(simd::abs4(simd::sub4(row2, m.row2)), Vector3<float>::TOLERANCE)) &&
		simd::all3(simd::lessThan(simd::abs4(simd::sub4(row3, m.row3)), Vector3<float>::TOLERANCE));
}

inline bool AffineTransform<float>::isApproxEqual(const AffineTransform& m, float tolerance) const
{
	auto t = simd::set4(tolerance);
	return simd::all3(simd::lessThan(simd::abs4(simd::sub4(row0, m.row0)), t)) &&
		simd::all3(simd::lessThan(simd::abs4(simd::sub4(row1, m.row1)), t)) &&
		simd::all3(simd::lessThan(simd::abs4(simd::sub4(row2, m.row2)), t)) &&
		simd::all3(simd::lessThan(simd::abs4(simd::sub4(row3, m.row3)), t));
}

inline bool AffineTransform<float>::isFinite() const
{
	return simd::all3(simd::isFinite(row0)) &&
		simd::all3(simd::isFinite(row1)) &&
		simd::all3(simd::isFinite(row2)) &&
		simd::all3(simd::isFinite(row3));
}

inline AffineTransform<float>& AffineTransform<float>::setZero()
{
	row0 = simd::zero<simd::float4>();
	row1 = simd::zero<simd::float4>();
	row2 = simd::zero<simd::float4>();
	row3 = simd::zero<simd::float4>();
	return *this;
}

inline AffineTransform<float>& AffineTransform<float>::setIdentity()
{
	row0 = Vector3<float>::UNIT_X;
	row1 = Vector3<float>::UNIT_Y;
	row2 = Vector3<float>::UNIT_Z;
	row3 = simd::zero<simd::float4>();
	return *this;
}

inline AffineTransform<float>& AffineTransform<float>::set(const Matrix3<float>& r, const Vector3<float>& t)
{
	row0 = r.row0;
	row1 = r.row1;
	row2 = r.row2;
	row3 = t;
	return *this;
}

inline AffineTransform<float>& AffineTransform<float>::set(const Vector3<float>& row0, const Vector3<float>& row1, 
	const Vector3<float>& row2, const Vector3<float>& row3)
{
	this->row0 = row0;
	this->row1 = row1;
	this->row2 = row2;
	this->row3 = row3;
	return *this;
}

inline AffineTransform<float>& AffineTransform<float>::set(simd::float4 row0, simd::float4 row1, simd::float4 row2, simd::float4 row3)
{
	this->row0 = row0;
	this->row1 = row1;
	this->row2 = row2;
	this->row3 = row3;
	return *this;
}

inline AffineTransform<float>& AffineTransform<float>::set(float m00, float m01, float m02, float m10, float m11, float m12,
	float m20, float m21, float m22, float m30, float m31, float m32)
{
#if MATHEMATICS_SIMD_EXPAND_LAST
	row0 = simd::set4(m00, m01, m02, m02);
	row1 = simd::set4(m10, m11, m12, m12);
	row2 = simd::set4(m20, m21, m22, m22);
	row3 = simd::set4(m30, m31, m32, m32);
#else
	row0 = simd::set3(m00, m01, m02);
	row1 = simd::set3(m10, m11, m12);
	row2 = simd::set3(m20, m21, m22);
	row3 = simd::set3(m30, m31, m32);
#endif
	return *this;
}

inline AffineTransform<float>& AffineTransform<float>::setTranslation(const Vector3<float>& v)
{
	row0 = Vector3<float>::UNIT_X;
	row1 = Vector3<float>::UNIT_Y;
	row2 = Vector3<float>::UNIT_Z;
	row3 = v;
	return *this;
}

inline AffineTransform<float>& AffineTransform<float>::setScaling(const Vector3<float>& v)
{
	row0 = simd::cutoff1(v);
	row1 = simd::and4(v, simd::constant4i<simd::float4, 0, -1, 0, 0>());
#if MATHEMATICS_SIMD_EXPAND_LAST
	row2 = simd::and4(v, simd::constant4i<simd::float4, 0, 0, -1, -1>());
#else
	row2 = simd::and4(v, simd::constant4i<simd::float4, 0, 0, -1, 0>());
#endif
	row3 = simd::zero<simd::float4>();
	return *this;
}

inline AffineTransform<float>& AffineTransform<float>::setScaling(float f)
{
	const auto zero = simd::zero<simd::float4>();
	row0 = simd::set1(f);
	row1 = simd::insert<1>(f, zero);
#if MATHEMATICS_SIMD_EXPAND_LAST
	row2 = simd::set4(0.f, 0.f, f, f);
#else
	row2 = simd::insert<2>(f, zero);
#endif
	row3 = zero;
	return *this;
}

inline AffineTransform<float>& AffineTransform<float>::setRotation(Axis axis, float angle)
{
	if (angle != 0.f)
	{
		float sinAngle = std::sin(angle);
		float cosAngle = std::cos(angle);

		switch (axis)
		{
			case Axis::X:
				return set(1.f, 0.f, 0.f, 0.f, cosAngle, sinAngle, 0.f, -sinAngle, cosAngle, 0.f, 0.f, 0.f);
			case Axis::Y:
				return set(cosAngle, 0.f, -sinAngle, 0.f, 1.f, 0.f, sinAngle, 0.f, cosAngle, 0.f, 0.f, 0.f);
			case Axis::Z:
				return set(cosAngle, sinAngle, 0.f, -sinAngle, cosAngle, 0.f, 0.f, 0.f, 1.f, 0.f, 0.f, 0.f);
			default:
				return setZero();
		}
	}

	return setIdentity();
}

inline AffineTransform<float>& AffineTransform<float>::setRotation(const Vector3<float>& axis, float angle)
{
	float m = axis.getMagnitude();
	if ((m > 0.f) && (angle != 0.f))
	{
		Vector3<float> abc = axis/m;
		float sinAngle = std::sin(angle);
		float cosAngle = std::cos(angle);
		float mCos = 1.f - cosAngle;
		float aa = abc.x*abc.x, bb = abc.y*abc.y, cc = abc.z*abc.z;
		float abm = abc.x*abc.y*mCos, acm = abc.x*abc.z*mCos, bcm = abc.y*abc.z*mCos;
		float aSin = abc.x*sinAngle, bSin = abc.y*sinAngle, cSin = abc.z*sinAngle;

		return set(aa*mCos + cosAngle, abm + cSin, acm - bSin,
			abm - cSin, bb*mCos + cosAngle, bcm + aSin,
			acm + bSin, bcm - aSin, cc*mCos + cosAngle,
			0.f, 0.f, 0.f);
	}

	return setIdentity();
}

inline AffineTransform<float>& AffineTransform<float>::setScalingTranslation(const Vector3<float>& s, const Vector3<float>& t)
{
	row0 = simd::cutoff1(s);
	row1 = simd::and4(s, simd::constant4i<simd::float4, 0, -1, 0, 0>());
#if MATHEMATICS_SIMD_EXPAND_LAST
	row2 = simd::and4(s, simd::constant4i<simd::float4, 0, 0, -1, -1>());
#else
	row2 = simd::and4(s, simd::constant4i<simd::float4, 0, 0, -1, 0>());
#endif
	row3 = t;
	return *this;
}

inline AffineTransform<float>& AffineTransform<float>::setShearing(float xy, float xz, float yx, float yz, float zx, float zy)
{
#if MATHEMATICS_SIMD_EXPAND_LAST
	row0 = simd::set4(1.f, xy, xz, xz);
	row1 = simd::set4(yx, 1.f, yz, yz);
	row2 = simd::set4(zx, zy, 1.f, 1.f);
#else
	row0 = simd::set3(1.f, xy, xz);
	row1 = simd::set3(yx, 1.f, yz);
	row2 = simd::set3(zx, zy, 1.f);
#endif
	row3 = simd::zero<simd::float4>();
	return *this;
}

inline AffineTransform<float>& AffineTransform<float>::setInverse(const AffineTransform<float>& m)
{
	float a = m.m11*m.m22 - m.m12*m.m21; // #TODO SIMD
	float b = m.m12*m.m20 - m.m10*m.m22;
	float c = m.m10*m.m21 - m.m11*m.m20;
	float d = m.m02*m.m21 - m.m01*m.m22;
	float e = m.m00*m.m22 - m.m02*m.m20;
	float f = m.m01*m.m20 - m.m00*m.m21;
	float g = m.m01*m.m12 - m.m02*m.m11;
	float h = m.m02*m.m10 - m.m00*m.m12;
	float i = m.m00*m.m11 - m.m01*m.m10;
	auto n = simd::set4(1.f/(m.m00*a + m.m01*b + m.m02*c));
#if MATHEMATICS_SIMD_EXPAND_LAST
	row0 = simd::mul4(simd::set4(a, d, g, g), n);
	row1 = simd::mul4(simd::set4(b, e, h, h), n);
	row2 = simd::mul4(simd::set4(c, f, i, i), n);
#else
	row0 = simd::mul4(simd::set3(a, d, g), n);
	row1 = simd::mul4(simd::set3(b, e, h), n);
	row2 = simd::mul4(simd::set3(c, f, i), n);
#endif
	setOrigin(-(m.getOrigin()*getBasis()));
	return *this;
}

inline AffineTransform<float>& AffineTransform<float>::setInverseOrthogonal(const AffineTransform<float>& m)
{
#if MATHEMATICS_SIMD_EXPAND_LAST
	auto [r0, r1, r2] = simd::transpose3x3(m.row0, m.row1, m.row2);
	row0 = simd::xyzz(r0);
	row1 = simd::xyzz(r1);
	row2 = simd::xyzz(r2);
#else
	//simd::transpose3x3(m.row0, m.row1, m.row2, row0, row1, row2);
	std::tie(row0, row1, row2) = simd::transpose3x3(m.row0, m.row1, m.row2);
#endif
	setOrigin(-(m.getOrigin()*getBasis()));
	return *this;
}

inline AffineTransform<float>& AffineTransform<float>::preConcatenate(const AffineTransform<float>& m)
{
	auto r0 = simd::mul4(simd::broadcast<0>(m.row0), row0);
	r0 = simd::add4(r0, simd::mul4(simd::broadcast<1>(m.row0), row1));
	r0 = simd::add4(r0, simd::mul4(simd::broadcast<2>(m.row0), row2));
	auto r1 = simd::mul4(simd::broadcast<0>(m.row1), row0);
	r1 = simd::add4(r1, simd::mul4(simd::broadcast<1>(m.row1), row1));
	r1 = simd::add4(r1, simd::mul4(simd::broadcast<2>(m.row1), row2));
	auto r2 = simd::mul4(simd::broadcast<0>(m.row2), row0);
	r2 = simd::add4(r2, simd::mul4(simd::broadcast<1>(m.row2), row1));
	r2 = simd::add4(r2, simd::mul4(simd::broadcast<2>(m.row2), row2));
	auto r3 = simd::mul4(simd::broadcast<0>(m.row3), row0);
	r3 = simd::add4(r3, simd::mul4(simd::broadcast<1>(m.row3), row1));
	r3 = simd::add4(r3, simd::mul4(simd::broadcast<2>(m.row3), row2));
	r3 = simd::add4(r3, row3);
	row0 = r0;
	row1 = r1;
	row2 = r2;
	row3 = r3;
	return *this;
}

inline AffineTransform<float>& AffineTransform<float>::concatenate(const AffineTransform<float>& m)
{
	auto r0 = simd::mul4(simd::broadcast<0>(row0), m.row0);
	r0 = simd::add4(r0, simd::mul4(simd::broadcast<1>(row0), m.row1));
	r0 = simd::add4(r0, simd::mul4(simd::broadcast<2>(row0), m.row2));
	auto r1 = simd::mul4(simd::broadcast<0>(row1), m.row0);
	r1 = simd::add4(r1, simd::mul4(simd::broadcast<1>(row1), m.row1));
	r1 = simd::add4(r1, simd::mul4(simd::broadcast<2>(row1), m.row2));
	auto r2 = simd::mul4(simd::broadcast<0>(row2), m.row0);
	r2 = simd::add4(r2, simd::mul4(simd::broadcast<1>(row2), m.row1));
	r2 = simd::add4(r2, simd::mul4(simd::broadcast<2>(row2), m.row2));
	auto r3 = simd::mul4(simd::broadcast<0>(row3), m.row0);
	r3 = simd::add4(r3, simd::mul4(simd::broadcast<1>(row3), m.row1));
	r3 = simd::add4(r3, simd::mul4(simd::broadcast<2>(row3), m.row2));
	r3 = simd::add4(r3, m.row3);
	row0 = r0;
	row1 = r1;
	row2 = r2;
	row3 = r3;
	return *this;
}

inline AffineTransform<float>& AffineTransform<float>::preScale(const Vector3<T>& v)
{
	row0 = simd::mul4(row0, simd::xxxx(v));
	row1 = simd::mul4(row1, simd::yyyy(v));
	row2 = simd::mul4(row2, simd::zzzz(v));
	return *this;
}

inline AffineTransform<float>& AffineTransform<float>::preScale(float f)
{
	auto t = simd::set4(f);
	row0 = simd::mul4(row0, t);
	row1 = simd::mul4(row1, t);
	row2 = simd::mul4(row2, t);
	return *this;
}

inline AffineTransform<float>& AffineTransform<float>::scale(const Vector3<float>& v)
{
	row0 = simd::mul4(row0, v);
	row1 = simd::mul4(row1, v);
	row2 = simd::mul4(row2, v);
	row3 = simd::mul4(row3, v);
	return *this;
}

inline AffineTransform<float>& AffineTransform<float>::scale(float f)
{
	auto t = simd::set4(f);
	row0 = simd::mul4(row0, t);
	row1 = simd::mul4(row1, t);
	row2 = simd::mul4(row2, t);
	row3 = simd::mul4(row3, t);
	return *this;
}

inline AffineTransform<float>& AffineTransform<float>::invert()
{
	float a = m11*m22 - m12*m21; // #TODO SIMD
	float b = m12*m20 - m10*m22;
	float c = m10*m21 - m11*m20;
	float d = m02*m21 - m01*m22;
	float e = m00*m22 - m02*m20;
	float f = m01*m20 - m00*m21;
	float g = m01*m12 - m02*m11;
	float h = m02*m10 - m00*m12;
	float i = m00*m11 - m01*m10;
	auto n = simd::set4(1.f/(m00*a + m01*b + m02*c));
#if MATHEMATICS_SIMD_EXPAND_LAST
	row0 = simd::mul4(simd::set4(a, d, g, g), n);
	row1 = simd::mul4(simd::set4(b, e, h, h), n);
	row2 = simd::mul4(simd::set4(c, f, i, i), n);
#else
	row0 = simd::mul4(simd::set3(a, d, g), n);
	row1 = simd::mul4(simd::set3(b, e, h), n);
	row2 = simd::mul4(simd::set3(c, f, i), n);
#endif
	setOrigin(-(getOrigin()*getBasis()));
	return *this;
}

inline AffineTransform<float>& AffineTransform<float>::invertOrthogonal()
{
#if MATHEMATICS_SIMD_EXPAND_LAST
	auto [r0, r1, r2] = simd::transpose3x3(row0, row1, row2);
	row0 = simd::xyzz(r0);
	row1 = simd::xyzz(r1);
	row2 = simd::xyzz(r2);
#else
	//simd::transpose3x3(row0, row1, row2, row0, row1, row2);
	std::tie(row0, row1, row2) = simd::transpose3x3(row0, row1, row2);
#endif
	setOrigin(-(getOrigin()*getBasis()));
	return *this;
}

#endif /* SIMD_HAS_FLOAT4 */

template<typename T>
	requires std::floating_point<T>
inline AffineTransform<T> concatenate(const AffineTransform<T>& m1, const AffineTransform<T>& m2) noexcept
{
	return AffineTransform<T>(m1.m00*m2.m00 + m1.m01*m2.m10 + m1.m02*m2.m20,
		m1.m00*m2.m01 + m1.m01*m2.m11 + m1.m02*m2.m21,
		m1.m00*m2.m02 + m1.m01*m2.m12 + m1.m02*m2.m22,
		m1.m10*m2.m00 + m1.m11*m2.m10 + m1.m12*m2.m20,
		m1.m10*m2.m01 + m1.m11*m2.m11 + m1.m12*m2.m21,
		m1.m10*m2.m02 + m1.m11*m2.m12 + m1.m12*m2.m22,
		m1.m20*m2.m00 + m1.m21*m2.m10 + m1.m22*m2.m20,
		m1.m20*m2.m01 + m1.m21*m2.m11 + m1.m22*m2.m21,
		m1.m20*m2.m02 + m1.m21*m2.m12 + m1.m22*m2.m22,
		m1.m30*m2.m00 + m1.m31*m2.m10 + m1.m32*m2.m20 + m2.m30,
		m1.m30*m2.m01 + m1.m31*m2.m11 + m1.m32*m2.m21 + m2.m31,
		m1.m30*m2.m02 + m1.m31*m2.m12 + m1.m32*m2.m22 + m2.m32);
}

template<typename T>
	requires std::floating_point<T>
inline AffineTransform<T> concatenate(const AffineTransform<T>& m1, const AffineTransform<T>& m2, const AffineTransform<T>& m3) noexcept
{
	return concatenate(concatenate(m1, m2), m3);
}

template<typename T>
	requires std::floating_point<T>
inline AffineTransform<T> concatenate(const AffineTransform<T>& m1, const AffineTransform<T>& m2, const AffineTransform<T>& m3,
	const AffineTransform<T>& m4) noexcept
{
	return concatenate(concatenate(concatenate(m1, m2), m3), m4);
}

template<typename T>
	requires std::floating_point<T>
inline AffineTransform<T> inverse(const AffineTransform<T>& m) noexcept
{
	return AffineTransform<T>(Uninitialized()).setInverse(m);
}

template<typename T>
	requires std::floating_point<T>
inline AffineTransform<T> inverse(AffineTransform<T>&& m) noexcept
{
	m.invert();
	return m;
}

template<typename T>
	requires std::floating_point<T>
inline AffineTransform<T> inverseOrthogonal(const AffineTransform<T>& m) noexcept
{
	AffineTransform<T> n(Uninitialized());
	n.m00 = m.m00; n.m01 = m.m10; n.m02 = m.m20;
	n.m10 = m.m01; n.m11 = m.m11; n.m12 = m.m21;
	n.m20 = m.m02; n.m21 = m.m12; n.m22 = m.m22;
	n.setOrigin(-(m.getOrigin()*n.getBasis()));
	return n;
}

template<typename T>
	requires std::floating_point<T>
inline AffineTransform<T> inverseOrthogonal(AffineTransform<T>&& m) noexcept
{
	T t = m.m01; m.m01 = m.m10; m.m10 = t;
	t = m.m02; m.m02 = m.m20; m.m20 = t;
	t = m.m12; m.m12 = m.m21; m.m21 = t;
	m.setOrigin(-(m.getOrigin()*m.getBasis()));
	return m;
}

#if SIMD_HAS_FLOAT4

template<>
inline AffineTransform<float> concatenate(const AffineTransform<float>& m1, const AffineTransform<float>& m2) noexcept
{
	auto r0 = simd::mul4(simd::broadcast<0>(m1.row0), m2.row0);
	r0 = simd::add4(r0, simd::mul4(simd::broadcast<1>(m1.row0), m2.row1));
	r0 = simd::add4(r0, simd::mul4(simd::broadcast<2>(m1.row0), m2.row2));
	auto r1 = simd::mul4(simd::broadcast<0>(m1.row1), m2.row0);
	r1 = simd::add4(r1, simd::mul4(simd::broadcast<1>(m1.row1), m2.row1));
	r1 = simd::add4(r1, simd::mul4(simd::broadcast<2>(m1.row1), m2.row2));
	auto r2 = simd::mul4(simd::broadcast<0>(m1.row2), m2.row0);
	r2 = simd::add4(r2, simd::mul4(simd::broadcast<1>(m1.row2), m2.row1));
	r2 = simd::add4(r2, simd::mul4(simd::broadcast<2>(m1.row2), m2.row2));
	auto r3 = simd::mul4(simd::broadcast<0>(m1.row3), m2.row0);
	r3 = simd::add4(r3, simd::mul4(simd::broadcast<1>(m1.row3), m2.row1));
	r3 = simd::add4(r3, simd::mul4(simd::broadcast<2>(m1.row3), m2.row2));
	r3 = simd::add4(r3, m2.row3);
	return AffineTransform<float>(r0, r1, r2, r3);
}

template<>
inline AffineTransform<float> inverseOrthogonal(const AffineTransform<float>& m) noexcept
{
	AffineTransform<T> n(Uninitialized());
#if MATHEMATICS_SIMD_EXPAND_LAST
	auto [r0, r1, r2] = simd::transpose3x3(m.row0, m.row1, m.row2);
	n.row0 = simd::xyzz(r0);
	n.row1 = simd::xyzz(r1);
	n.row2 = simd::xyzz(r2);
#else
	//simd::transpose3x3(m.row0, m.row1, m.row2, n.row0, n.row1, n.row2); 
	std::tie(n.row0, n.row1, n.row2) = simd::transpose3x3(m.row0, m.row1, m.row2);
#endif
	n.setOrigin(-(m.getOrigin()*n.getBasis()));
	return n;
}

template<>
inline AffineTransform<float> inverseOrthogonal(AffineTransform<float>&& m) noexcept
{
#if MATHEMATICS_SIMD_EXPAND_LAST
	auto [r0, r1, r2] = simd::transpose3x3(m.row0, m.row1, m.row2);
	m.row0 = simd::xyzz(r0);
	m.row1 = simd::xyzz(r1);
	m.row2 = simd::xyzz(r2);
#else
	//simd::transpose3x3(m.row0, m.row1, m.row2, m.row0, m.row1, m.row2); 
	std::tie(m.row0, m.row1, m.row2) = simd::transpose3x3(m.row0, m.row1, m.row2);
#endif
	m.setOrigin(-(m.getOrigin()*m.getBasis()));
	return m;
}

#endif /* SIMD_HAS_FLOAT4 */

} // namespace templates

//using AffineTransformf = templates::AffineTransform<float>;
//using AffineTransformd = templates::AffineTransform<double>;

#if MATHEMATICS_DOUBLE
using AffineTransform = templates::AffineTransform<double>;
using AffineTransformArg = templates::AffineTransform<double>::ConstArg;
using AffineTransformResult = templates::AffineTransform<double>::ConstResult;
#else
using AffineTransform = templates::AffineTransform<float>;
using AffineTransformArg = templates::AffineTransform<float>::ConstArg;
using AffineTransformResult = templates::AffineTransform<float>::ConstResult;
#endif

} // namespace core::mathematics

#include "Quaternion.hpp"
#include "YawPitchRoll.hpp"
#include "Euler.hpp"
#include "Plane.hpp"

namespace core::mathematics::templates {

template<typename T>
inline AffineTransform<T>& AffineTransform<T>::setRotation(const YawPitchRoll<T>& r)
{
	if (!r.isZero())
	{
		T sy = std::sin(r.yaw), cy = std::cos(r.yaw);
		T sp = std::sin(r.pitch), cp = std::cos(r.pitch);
		T sr = std::sin(r.roll), cr = std::cos(r.roll);
		T ss = sy*sp, cs = cy*sp;

		return set(cr*cy + sr*ss, sr*cp, sr*cs - cr*sy,
			cr*ss - sr*cy, cr*cp, cr*cs + sr*sy,
			sy*cp, -sp, cy*cp,
			T(), T(), T());
	}

	return setIdentity();
}

template<typename T>
inline AffineTransform<T>& AffineTransform<T>::setRotation(const Euler<T>& e)
{
	if (!e.isZero() && (e.order != EulerOrder::UNSPECIFIED))
	{
		static const int safe[] = { 0, 1, 2, 0 };
		static const int next[] = { 1, 2, 0, 1 };
		unsigned int o = (unsigned int)e.order;
		int f = o & 1; o >>= 1;
		int s = o & 1; o >>= 1;
		int n = o & 1; o >>= 1;
		int i = safe[o & 3];
		int j = next[i + n];
		int k = next[i + 1 - n];
		T ti = e[f ? (s ? i : k) : i];
		T tj = e[j];
		T th = e[f ? i : (s ? i : k)];

		if (f)
		{
			T t = ti;
			ti = th;
			th = t;
		}

		if (n)
		{
			ti = -ti;
			tj = -tj;
			th = -th;
		}

		T si = std::sin(ti), ci = std::cos(ti);
		T sj = std::sin(tj), cj = std::cos(tj);
		T sh = std::sin(th), ch = std::cos(th);
		T cc = ci*ch, cs = ci*sh, sc = si*ch, ss = si*sh;

		T m[3][3];
		if (s)
		{
			m[i][i] = cj; m[j][i] = sj*si; m[k][i] = sj*ci;
			m[i][j] = sj*sh; m[j][j] = -cj*ss+cc; m[k][j] = -cj*cs-sc;
			m[i][k] = -sj*ch; m[j][k] = cj*sc+cs; m[k][k] = cj*cc-ss;
		}
		else
		{
			m[i][i] = cj*ch; m[j][i] = sj*sc-cs; m[k][i] = sj*cc+ss;
			m[i][j] = cj*sh; m[j][j] = sj*ss+cc; m[k][j] = sj*cs-sc;
			m[i][k] = -sj; m[j][k] = cj*si; m[k][k] = cj*ci;
		}

		return set(m[0][0], m[0][1], m[0][2],
			m[1][0], m[1][1], m[1][2],
			m[2][0], m[2][1], m[2][2],
			T(), T(), T());
	}

	return setIdentity();
}

template<typename T>
inline AffineTransform<T>& AffineTransform<T>::setRotation(const Quaternion<T>& q)
{
	if (!q.isIdentity())
	{
		T nq = q.getNorm();
		T s = (nq > T(0)) ? T(2)/nq : T(2);
		T xs = q.x*s, ys = q.y*s, zs = q.z*s;
		T wx = q.w*xs, wy = q.w*ys, wz = q.w*zs;
		T xx = q.x*xs, xy = q.x*ys, xz = q.x*zs;
		T yy = q.y*ys, yz = q.y*zs, zz = q.z*zs;

		return set(T(1) - (yy + zz), xy + wz, xz - wy,
			xy - wz, T(1) - (xx + zz), yz + wx,
			xz + wy, yz - wx, T(1) - (xx + yy),
			T(), T(), T());
	}

	return setIdentity();
}

template<typename T>
inline AffineTransform<T>& AffineTransform<T>::setReflection(const Plane<T>& plane)
{
	Plane<T> p(plane);
	p.normalize();
	m00 = T(-2)*p.a*p.a + T(1); m01 = T(-2)*p.b*p.a; m02 = T(-2)*p.c*p.a;
	m10 = T(-2)*p.a*p.b; m11 = T(-2)*p.b*p.b + T(1); m12 = T(-2)*p.c*p.b;
	m20 = T(-2)*p.a*p.c; m21 = T(-2)*p.b*p.c; m22 = T(-2)*p.c*p.c + T(1);
	m30 = T(-2)*p.a*p.d; m31 = T(-2)*p.b*p.d; m32 = T(-2)*p.c*p.d;
	return *this;
}

#if SIMD_HAS_FLOAT4

inline AffineTransform<float>& AffineTransform<float>::setRotation(const YawPitchRoll<float>& r)
{
	if (!r.isZero())
	{
		float sy = std::sin(r.yaw), cy = std::cos(r.yaw);
		float sp = std::sin(r.pitch), cp = std::cos(r.pitch);
		float sr = std::sin(r.roll), cr = std::cos(r.roll);
		float ss = sy*sp, cs = cy*sp;

		return set(cr*cy + sr*ss, sr*cp, sr*cs - cr*sy,
			cr*ss - sr*cy, cr*cp, cr*cs + sr*sy,
			sy*cp, -sp, cy*cp,
			0.f, 0.f, 0.f);
	}

	return setIdentity();
}

inline AffineTransform<float>& AffineTransform<float>::setRotation(const Euler<float>& e)
{
	if (!e.isZero() && (e.order != EulerOrder::UNSPECIFIED))
	{
		static const int safe[] = { 0, 1, 2, 0 };
		static const int next[] = { 1, 2, 0, 1 };
		unsigned int o = (unsigned int)e.order;
		int f = o & 1; o >>= 1;
		int s = o & 1; o >>= 1;
		int n = o & 1; o >>= 1;
		int i = safe[o & 3];
		int j = next[i + n];
		int k = next[i + 1 - n];
		float ti = e[f ? (s ? i : k) : i];
		float tj = e[j];
		float th = e[f ? i : (s ? i : k)];

		if (f)
		{
			float t = ti;
			ti = th;
			th = t;
		}

		if (n)
		{
			ti = -ti;
			tj = -tj;
			th = -th;
		}

		float si = std::sin(ti), ci = std::cos(ti);
		float sj = std::sin(tj), cj = std::cos(tj);
		float sh = std::sin(th), ch = std::cos(th);
		float cc = ci*ch, cs = ci*sh, sc = si*ch, ss = si*sh;

		float m[3][3];
		if (s)
		{
			m[i][i] = cj; m[j][i] = sj*si; m[k][i] = sj*ci;
			m[i][j] = sj*sh; m[j][j] = -cj*ss+cc; m[k][j] = -cj*cs-sc;
			m[i][k] = -sj*ch; m[j][k] = cj*sc+cs; m[k][k] = cj*cc-ss;
		}
		else
		{
			m[i][i] = cj*ch; m[j][i] = sj*sc-cs; m[k][i] = sj*cc+ss;
			m[i][j] = cj*sh; m[j][j] = sj*ss+cc; m[k][j] = sj*cs-sc;
			m[i][k] = -sj; m[j][k] = cj*si; m[k][k] = cj*ci;
		}

		return set(m[0][0], m[0][1], m[0][2],
			m[1][0], m[1][1], m[1][2],
			m[2][0], m[2][1], m[2][2],
			0.f, 0.f, 0.f);
	}

	return setIdentity();
}

inline AffineTransform<float>& AffineTransform<float>::setRotation(const Quaternion<float>& q)
{
	if (!q.isIdentity())
	{
		float nq = q.getNorm();
		float s = (nq > 0.f) ? 2.f/nq : 2.f;
		float xs = q.x*s, ys = q.y*s, zs = q.z*s; // #TODO SIMD
		float wx = q.w*xs, wy = q.w*ys, wz = q.w*zs;
		float xx = q.x*xs, xy = q.x*ys, xz = q.x*zs;
		float yy = q.y*ys, yz = q.y*zs, zz = q.z*zs;

		return set(1.f - (yy + zz), xy + wz, xz - wy,
			xy - wz, 1.f - (xx + zz), yz + wx,
			xz + wy, yz - wx, 1.f - (xx + yy),
			0.f, 0.f, 0.f);
	}

	return setIdentity();
}

inline AffineTransform<float>& AffineTransform<float>::setReflection(const Plane<float>& plane)
{
	Plane<float> p(plane);
	p.normalize();
	static const auto mTwo = simd::set4(-2.f);
#if MATHEMATICS_SIMD_EXPAND_LAST
	auto n = simd::xyzz(p);
#else
	auto n = simd::cutoff3(p);
#endif
	auto t = simd::mul4(mTwo, p);
	row0 = simd::mulAdd4(simd::xxxx(t), n, Vector3<float>::UNIT_X);
	row1 = simd::mulAdd4(simd::yyyy(t), n, Vector3<float>::UNIT_Y);
	row2 = simd::mulAdd4(simd::zzzz(t), n, Vector3<float>::UNIT_Z);
	row3 = simd::mul4(simd::wwww(t), n);
	return *this;
}

#endif /* SIMD_HAS_FLOAT4 */

} // namespace core::mathematics::templates
