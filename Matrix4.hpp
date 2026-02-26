/*
 *	Name: Matrix4
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
#include "Vector3.hpp"
#include "Vector4.hpp"
#include "Matrix2.hpp"
#include "Matrix3.hpp"
#include "AffineTransform.hpp"

namespace core::mathematics {
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
struct Matrix4
{
	using Real = T;
	using ComponentType = T;
	using ConstArg = const Matrix4&;
	using ConstResult = const Matrix4&;

	//static constexpr int NUM_COMPONENTS = 16;

	constexpr Matrix4() noexcept : m00(), m01(), m02(), m03(), m10(), m11(), m12(), m13(), m20(), m21(), m22(), m23(), m30(), m31(), m32(), m33() {}
	explicit Matrix4(Uninitialized) noexcept {}
	explicit Matrix4(Identity) noexcept : m00(1), m01(), m02(), m03(), m10(), m11(1), m12(), m13(), m20(), m21(), m22(1), m23(), m30(), m31(), m32(), m33(1) {}
	constexpr Matrix4(T m00, T m01, T m02, T m03, T m10, T m11, T m12, T m13, T m20, T m21, T m22, T m23, T m30, T m31, T m32, T m33) noexcept;
	constexpr Matrix4(const Vector4<T>& row0, const Vector4<T>& row1, const Vector4<T>& row2, const Vector4<T>& row3) noexcept;
	constexpr explicit Matrix4(const std::tuple<Vector4<T>, Vector4<T>, Vector4<T>, Vector4<T>>& t) noexcept;
	constexpr Matrix4(const Matrix2<T>& m) noexcept;
	constexpr Matrix4(const Matrix3<T>& m) noexcept;
	constexpr Matrix4(const AffineTransform<T>& m) noexcept;
	explicit Matrix4(const T* m) noexcept;

	explicit operator T*() noexcept { return &m00; }
	explicit operator const T*() const noexcept { return &m00; }
	Vector4<T>& operator[](int i) noexcept { return reinterpret_cast<Vector4<T>*>(&m00)[i]; }
	const Vector4<T>& operator[](int i) const noexcept { return reinterpret_cast<const Vector4<T>*>(&m00)[i]; }

	Matrix4 operator+() const noexcept { return *this; }
	Matrix4 operator-() const noexcept;
	Matrix4& operator+=(const Matrix4& m) noexcept;
	Matrix4& operator-=(const Matrix4& m) noexcept;
	Matrix4& operator*=(T f) noexcept;
	Matrix4& operator/=(T f) noexcept { return operator*=(T(1)/f); }
	Matrix4& operator*=(const Matrix4& m) noexcept;
	bool operator==(const Matrix4& m) const noexcept;
	bool operator!=(const Matrix4& m) const noexcept { return !(*this == m); }

	template<typename A> void load(A& ar) { ar(m00, m01, m02, m03, m10, m11, m12, m13, m20, m21, m22, m23, m30, m31, m32, m33); }
	template<typename A> void save(A& ar) const { ar(m00, m01, m02, m03, m10, m11, m12, m13, m20, m21, m22, m23, m30, m31, m32, m33); }

	// #TODO rename get...() to make...()
	static Matrix4 getTranslation(const Vector3<T>& v) noexcept { return Matrix4(Uninitialized()).setTranslation(v); }
	static Matrix4 getScaling(const Vector3<T>& v) noexcept { return Matrix4(Uninitialized()).setScaling(v); }
	static Matrix4 getRotation(Axis axis, T angle) noexcept { return Matrix4(Uninitialized()).setRotation(axis, angle); }
	static Matrix4 getRotation(const Vector3<T>& axis, T angle) noexcept { return Matrix4(Uninitialized()).setRotation(axis, angle); }
	static Matrix4 getRotation(const YawPitchRoll<T>& r) noexcept { return Matrix4(Uninitialized()).setRotation(r); }
	static Matrix4 getRotation(const Euler<T>& e) noexcept { return Matrix4(Uninitialized()).setRotation(e); }
	static Matrix4 getRotation(const Quaternion<T>& q) noexcept { return Matrix4(Uninitialized()).setRotation(q); }
	static Matrix4 getRotationTranslation(Axis axis, T angle, const Vector3<T>& t) noexcept { return Matrix4(Uninitialized()).setRotationTranslation(axis, angle, t); }
	static Matrix4 getRotationTranslation(const Vector3<T>& axis, T angle, const Vector3<T>& t) noexcept { return Matrix4(Uninitialized()).setRotationTranslation(axis, angle, t); }
	static Matrix4 getRotationTranslation(const YawPitchRoll<T>& r, const Vector3<T>& t) noexcept { return Matrix4(Uninitialized()).setRotationTranslation(r, t); }
	static Matrix4 getRotationTranslation(const Euler<T>& e, const Vector3<T>& t) noexcept { return Matrix4(Uninitialized()).setRotationTranslation(e, t); }
	static Matrix4 getRotationTranslation(const Quaternion<T>& q, const Vector3<T>& t) noexcept { return Matrix4(Uninitialized()).setRotationTranslation(q, t); }
	static Matrix4 getScalingRotation(const Vector3<T>& s, Axis axis, T angle) noexcept { return Matrix4(Uninitialized()).setScalingRotation(s, axis, angle); }
	static Matrix4 getScalingRotation(const Vector3<T>& s, const Vector3<T>& axis, T angle) noexcept { return Matrix4(Uninitialized()).setScalingRotation(s, axis, angle); }
	static Matrix4 getScalingRotation(const Vector3<T>& s, const YawPitchRoll<T>& r) noexcept { return Matrix4(Uninitialized()).setScalingRotation(s, r); }
	static Matrix4 getScalingRotation(const Vector3<T>& s, const Euler<T>& e) noexcept { return Matrix4(Uninitialized()).setScalingRotation(s, e); }
	static Matrix4 getScalingRotation(const Vector3<T>& s, const Quaternion<T>& q) noexcept { return Matrix4(Uninitialized()).setScalingRotation(s, q); }
	static Matrix4 getScalingTranslation(const Vector3<T>& s, const Vector3<T>& t) noexcept { return Matrix4(Uninitialized()).setScalingTranslation(s, t); }
	static Matrix4 getScalingRotationTranslation(const Vector3<T>& s, Axis axis, T angle, const Vector3<T>& t) noexcept { return Matrix4(Uninitialized()).setScalingRotationTranslation(s, axis, angle, t); }
	static Matrix4 getScalingRotationTranslation(const Vector3<T>& s, const Vector3<T>& axis, T angle, const Vector3<T>& t) noexcept { return Matrix4(Uninitialized()).setScalingRotationTranslation(s, axis, angle, t); }
	static Matrix4 getScalingRotationTranslation(const Vector3<T>& s, const YawPitchRoll<T>& r, const Vector3<T>& t) noexcept { return Matrix4(Uninitialized()).setScalingRotationTranslation(s, r, t); }
	static Matrix4 getScalingRotationTranslation(const Vector3<T>& s, const Euler<T>& e, const Vector3<T>& t) noexcept { return Matrix4(Uninitialized()).setScalingRotationTranslation(s, e, t); }
	static Matrix4 getScalingRotationTranslation(const Vector3<T>& s, const Quaternion<T>& q, const Vector3<T>& t) noexcept { return Matrix4(Uninitialized()).setScalingRotationTranslation(s, q, t); }
	static Matrix4 getShearing(T xy, T xz, T yx, T yz, T zx, T zy) noexcept { return Matrix4(Uninitialized()).setShearing(xy, xz, yx, yz, zx, zy); }
	static Matrix4 getReflection(const Plane<T>& plane) noexcept { return Matrix4(Uninitialized()).setReflection(plane); }

	bool isZero() const noexcept;
	bool isApproxZero() const noexcept;
	bool isIdentity() const noexcept;
	bool isApproxIdentity() const noexcept;
	bool isApproxEqual(const Matrix4& m) const noexcept;
	bool isApproxEqual(const Matrix4& m, T tolerance) const noexcept;
	bool isFinite() const noexcept;
	const Vector4<T>& getRow(int i) const noexcept { return reinterpret_cast<const Vector4<T>*>(&m00)[i]; }
	Vector4<T> getColumn(int i) const noexcept { return Vector4<T>((&m00)[i], (&m10)[i], (&m20)[i], (&m30)[i]); }
	T getTrace() const noexcept { return (m00 + m11 + m22 + m33); }
	T getDeterminant() const noexcept;
	bool isSingular() const noexcept { return (getDeterminant() == T(0)); }
	Matrix4& setZero/*zero*/() noexcept;
	Matrix4& setIdentity/*makeIdentity*/() noexcept;
	Matrix4& set(const Vector4<T>& row0, const Vector4<T>& row1, const Vector4<T>& row2, const Vector4<T>& row3) noexcept;
	Matrix4& set(T m00, T m01, T m02, T m03, T m10, T m11, T m12, T m13, T m20, T m21, T m22, T m23, T m30, T m31, T m32, T m33) noexcept;
	Matrix4& setTranslation/*makeTranslation*/(const Vector3<T>& v) noexcept;
	Matrix4& setScaling/*makeScaling*/(const Vector3<T>& v) noexcept;
	Matrix4& setRotation(Axis axis, T angle) noexcept;
	Matrix4& setRotation(const Vector3<T>& axis, T angle) noexcept;
	Matrix4& setRotation(const YawPitchRoll<T>& r) noexcept;
	Matrix4& setRotation(const Euler<T>& e) noexcept;
	Matrix4& setRotation(const Quaternion<T>& q) noexcept;
	Matrix4& setRotationTranslation(Axis axis, T angle, const Vector3<T>& t) noexcept;
	Matrix4& setRotationTranslation(const Vector3<T>& axis, T angle, const Vector3<T>& t) noexcept;
	Matrix4& setRotationTranslation(const YawPitchRoll<T>& r, const Vector3<T>& t) noexcept;
	Matrix4& setRotationTranslation(const Euler<T>& e, const Vector3<T>& t) noexcept;
	Matrix4& setRotationTranslation(const Quaternion<T>& q, const Vector3<T>& t) noexcept;
	Matrix4& setScalingRotation(const Vector3<T>& s, Axis axis, T angle) noexcept;
	Matrix4& setScalingRotation(const Vector3<T>& s, const Vector3<T>& axis, T angle) noexcept;
	Matrix4& setScalingRotation(const Vector3<T>& s, const YawPitchRoll<T>& r) noexcept;
	Matrix4& setScalingRotation(const Vector3<T>& s, const Euler<T>& e) noexcept;
	Matrix4& setScalingRotation(const Vector3<T>& s, const Quaternion<T>& q) noexcept;
	Matrix4& setScalingTranslation(const Vector3<T>& s, const Vector3<T>& t) noexcept;
	Matrix4& setScalingRotationTranslation(const Vector3<T>& s, Axis axis, T angle, const Vector3<T>& t) noexcept;
	Matrix4& setScalingRotationTranslation(const Vector3<T>& s, const Vector3<T>& axis, T angle, const Vector3<T>& t) noexcept;
	Matrix4& setScalingRotationTranslation(const Vector3<T>& s, const YawPitchRoll<T>& r, const Vector3<T>& t) noexcept;
	Matrix4& setScalingRotationTranslation(const Vector3<T>& s, const Euler<T>& e, const Vector3<T>& t) noexcept;
	Matrix4& setScalingRotationTranslation(const Vector3<T>& s, const Quaternion<T>& q, const Vector3<T>& t) noexcept;
	Matrix4& setShearing(T xy, T xz, T yx, T yz, T zx, T zy) noexcept;
	Matrix4& setReflection(const Plane<T>& plane) noexcept;
	Matrix4& setTranspose/*transposeOf*/(const Matrix4& m) noexcept;
	Matrix4& setInverse/*inverseOf*/(const Matrix4& m) noexcept;
	Matrix4& setInverseTranspose/*inverseTransposeOf*/(const Matrix4& m) noexcept { return setInverse(m).transpose(); }
	Matrix4& preConcatenate(const Matrix4& m) noexcept;
	Matrix4& concatenate(const Matrix4& m) noexcept { return operator*=(m); }
	Matrix4& translate(const Vector3<T>& v) noexcept;
	Matrix4& preScale(const Vector3<T>& v) noexcept;
	Matrix4& scale(const Vector3<T>& v) noexcept;
	//Matrix4& rotate(Axis axis, T angle) noexcept { concatenate(Matrix4(Uninitialized()).setRotation(axis, angle)); return *this; }
	//Matrix4& rotate(const Vector3<T>& axis, float angle) noexcept { concatenate(Matrix4(Uninitialized()).setRotation(axis, angle)); return *this; }
	//Matrix4& rotate(const YawPitchRoll<T>& r) noexcept { concatenate(Matrix4(Uninitialized()).setRotation(r)); return *this; }
	//Matrix4& rotate(const Euler<T>& e) noexcept { concatenate(Matrix4(Uninitialized()).setRotation(e)); return *this; }
	//Matrix4& rotate(const Quaternion<T>& q) noexcept { concatenate(Matrix4(Uninitialized()).setRotation(q)); return *this; }
	//Matrix4& shear(T xy, T xz, T yx, T yz, T zx, T zy) noexcept { concatenate(Matrix4(Uninitialized()).setShearing(xy, xz, yx, yz, zx, zy)); return *this; }
	Matrix4& negate() noexcept;
	Matrix4& transpose() noexcept;
	Matrix4& invert() noexcept { setInverse(*this); }

	//static const Matrix4& getZero() noexcept { return ZERO; }
	//static const Matrix4& getIdentity() noexcept { return IDENTITY; }

	static const Matrix4 ZERO;
	static const Matrix4 IDENTITY;

	T m00, m01, m02, m03;
	T m10, m11, m12, m13;
	T m20, m21, m22, m23;
	T m30, m31, m32, m33;
};

template<typename T> const Matrix4<T> Matrix4<T>::ZERO{};
template<typename T> const Matrix4<T> Matrix4<T>::IDENTITY{ T(1), T(0), T(0), T(0), T(0), T(1), T(0), T(0), T(0), T(0), T(1), T(0), T(0), T(0), T(0), T(1) };

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
struct Matrix4<float>
{
	using Real = float;
	using ComponentType = float;
	using ConstArg = const Matrix4;
	using ConstResult = const Matrix4;

	//static constexpr int NUM_COMPONENTS = 16;

	/*constexpr*/ Matrix4() noexcept;
	explicit Matrix4(Uninitialized) noexcept {}
	explicit Matrix4(Identity) noexcept;
	/*constexpr*/ Matrix4(float m00, float m01, float m02, float m03, float m10, float m11, float m12, float m13, float m20, float m21, float m22, float m23, float m30, float m31, float m32, float m33) noexcept;
	/*constexpr*/ Matrix4(const Vector4<float>& row0, const Vector4<float>& row1, const Vector4<float>& row2, const Vector4<float>& row3) noexcept;
	/*constexpr*/ explicit Matrix4(const std::tuple<Vector4<float>, Vector4<float>, Vector4<float>, Vector4<float>>& t) noexcept;
	/*constexpr*/ Matrix4(const Matrix2<float>& m) noexcept;
	/*constexpr*/ Matrix4(const Matrix3<float>& m) noexcept;
	/*constexpr*/ Matrix4(const AffineTransform<float>& m) noexcept;
	explicit Matrix4(const float* m) noexcept;
	explicit Matrix4(const simd::float4* m) noexcept : row0(m[0]), row1(m[1]), row2(m[2]), row3(m[3]) {}
	Matrix4(simd::float4 row0, simd::float4 row1, simd::float4 row2, simd::float4 row3) noexcept;
	explicit Matrix4(const std::tuple<simd::float4, simd::float4, simd::float4, simd::float4>& t) noexcept;
	Matrix4(const Matrix4& m) noexcept : row0(m.row0), row1(m.row1), row2(m.row2), row3(m.row3) {}
	Matrix4& operator=(const Matrix4& m) noexcept { row0 = m.row0; row1 = m.row1; row2 = m.row2; row3 = m.row3; return *this; }

	explicit operator float*() noexcept { return &m00; }
	explicit operator const float*() const noexcept { return &m00; }
	explicit operator simd::float4*() noexcept { return &row0; }
	explicit operator const simd::float4*() const noexcept { return &row0; }
	Vector4<float>& operator[](int i) noexcept { return reinterpret_cast<Vector4<float>&>((&row0)[i]); }
	const Vector4<float>& operator[](int i) const noexcept { return reinterpret_cast<const Vector4<float>&>((&row0)[i]); }

	Matrix4 operator+() const noexcept { return *this; }
	Matrix4 operator-() const noexcept;
	Matrix4& operator+=(const Matrix4& m) noexcept;
	Matrix4& operator-=(const Matrix4& m) noexcept;
	Matrix4& operator*=(float f) noexcept;
	Matrix4& operator/=(float f) noexcept { return operator*=(1.f/f); }
	Matrix4& operator*=(const Matrix4& m) noexcept;
	bool operator==(const Matrix4& m) const noexcept;
	bool operator!=(const Matrix4& m) const noexcept { return !(*this == m); }

	template<typename A> void load(A& ar);
	template<typename A> void save(A& ar) const { ar(m00, m01, m02, m03, m10, m11, m12, m13, m20, m21, m22, m23, m30, m31, m32, m33); }

	// #TODO rename get...() to make...()
	static Matrix4 getTranslation(const Vector3<float>& v) noexcept { return Matrix4(Uninitialized()).setTranslation(v); }
	static Matrix4 getScaling(const Vector3<float>& v) noexcept { return Matrix4(Uninitialized()).setScaling(v); }
	static Matrix4 getRotation(Axis axis, float angle) noexcept { return Matrix4(Uninitialized()).setRotation(axis, angle); }
	static Matrix4 getRotation(const Vector3<float>& axis, float angle) noexcept { return Matrix4(Uninitialized()).setRotation(axis, angle); }
	static Matrix4 getRotation(const YawPitchRoll<float>& r) noexcept { return Matrix4(Uninitialized()).setRotation(r); }
	static Matrix4 getRotation(const Euler<float>& e) noexcept { return Matrix4(Uninitialized()).setRotation(e); }
	static Matrix4 getRotation(const Quaternion<float>& q) noexcept { return Matrix4(Uninitialized()).setRotation(q); }
	static Matrix4 getRotationTranslation(Axis axis, float angle, const Vector3<float>& t) noexcept { return Matrix4(Uninitialized()).setRotationTranslation(axis, angle, t); }
	static Matrix4 getRotationTranslation(const Vector3<float>& axis, float angle, const Vector3<float>& t) noexcept { return Matrix4(Uninitialized()).setRotationTranslation(axis, angle, t); }
	static Matrix4 getRotationTranslation(const YawPitchRoll<float>& r, const Vector3<float>& t) noexcept { return Matrix4(Uninitialized()).setRotationTranslation(r, t); }
	static Matrix4 getRotationTranslation(const Euler<float>& e, const Vector3<float>& t) noexcept { return Matrix4(Uninitialized()).setRotationTranslation(e, t); }
	static Matrix4 getRotationTranslation(const Quaternion<float>& q, const Vector3<float>& t) noexcept { return Matrix4(Uninitialized()).setRotationTranslation(q, t); }
	static Matrix4 getScalingRotation(const Vector3<float>& s, Axis axis, float angle) noexcept { return Matrix4(Uninitialized()).setScalingRotation(s, axis, angle); }
	static Matrix4 getScalingRotation(const Vector3<float>& s, const Vector3<float>& axis, float angle) noexcept { return Matrix4(Uninitialized()).setScalingRotation(s, axis, angle); }
	static Matrix4 getScalingRotation(const Vector3<float>& s, const YawPitchRoll<float>& r) noexcept { return Matrix4(Uninitialized()).setScalingRotation(s, r); }
	static Matrix4 getScalingRotation(const Vector3<float>& s, const Euler<float>& e) noexcept { return Matrix4(Uninitialized()).setScalingRotation(s, e); }
	static Matrix4 getScalingRotation(const Vector3<float>& s, const Quaternion<float>& q) noexcept { return Matrix4(Uninitialized()).setScalingRotation(s, q); }
	static Matrix4 getScalingTranslation(const Vector3<float>& s, const Vector3<float>& t) noexcept { return Matrix4(Uninitialized()).setScalingTranslation(s, t); }
	static Matrix4 getScalingRotationTranslation(const Vector3<float>& s, Axis axis, float angle, const Vector3<float>& t) noexcept { return Matrix4(Uninitialized()).setScalingRotationTranslation(s, axis, angle, t); }
	static Matrix4 getScalingRotationTranslation(const Vector3<float>& s, const Vector3<float>& axis, float angle, const Vector3<float>& t) noexcept { return Matrix4(Uninitialized()).setScalingRotationTranslation(s, axis, angle, t); }
	static Matrix4 getScalingRotationTranslation(const Vector3<float>& s, const YawPitchRoll<float>& r, const Vector3<float>& t) noexcept { return Matrix4(Uninitialized()).setScalingRotationTranslation(s, r, t); }
	static Matrix4 getScalingRotationTranslation(const Vector3<float>& s, const Euler<float>& e, const Vector3<float>& t) noexcept { return Matrix4(Uninitialized()).setScalingRotationTranslation(s, e, t); }
	static Matrix4 getScalingRotationTranslation(const Vector3<float>& s, const Quaternion<float>& q, const Vector3<float>& t) noexcept { return Matrix4(Uninitialized()).setScalingRotationTranslation(s, q, t); }
	static Matrix4 getShearing(float xy, float xz, float yx, float yz, float zx, float zy) noexcept { return Matrix4(Uninitialized()).setShearing(xy, xz, yx, yz, zx, zy); }
	static Matrix4 getReflection(const Plane<float>& plane) noexcept { return Matrix4(Uninitialized()).setReflection(plane); }

	bool isZero() const noexcept;
	bool isApproxZero() const noexcept;
	bool isIdentity() const noexcept;
	bool isApproxIdentity() const noexcept;
	bool isApproxEqual(const Matrix4& m) const noexcept;
	bool isApproxEqual(const Matrix4& m, float tolerance) const noexcept;
	bool isFinite() const noexcept;
	const Vector4<float>& getRow(int i) const noexcept { return reinterpret_cast<const Vector4<float>&>((&row0)[i]); }
	Vector4<float> getColumn(int i) const noexcept; // #TODO
	float getTrace() const noexcept { return (m00 + m11 + m22 + m33); }
	float getDeterminant() const noexcept;
	bool isSingular() const noexcept { return (getDeterminant() == 0.f); }
	Matrix4& setZero/*zero*/() noexcept;
	Matrix4& setIdentity/*makeIdentity*/() noexcept;
	Matrix4& set(const Vector4<float>& row0, const Vector4<float>& row1, const Vector4<float>& row2, const Vector4<float>& row3) noexcept;
	Matrix4& set(simd::float4 row0, simd::float4 row1, simd::float4 row2, simd::float4 row3) noexcept;
	Matrix4& set(float m00, float m01, float m02, float m03, float m10, float m11, float m12, float m13, float m20, float m21, float m22, float m23, float m30, float m31, float m32, float m33) noexcept;
	Matrix4& setTranslation/*makeTranslation*/(const Vector3<float>& v) noexcept;
	Matrix4& setScaling/*makeScaling*/(const Vector3<float>& v) noexcept;
	Matrix4& setRotation(Axis axis, float angle) noexcept;
	Matrix4& setRotation(const Vector3<float>& axis, float angle) noexcept;
	Matrix4& setRotation(const YawPitchRoll<float>& r) noexcept;
	Matrix4& setRotation(const Euler<float>& e) noexcept;
	Matrix4& setRotation(const Quaternion<float>& q) noexcept;
	Matrix4& setRotationTranslation(Axis axis, float angle, const Vector3<float>& t) noexcept;
	Matrix4& setRotationTranslation(const Vector3<float>& axis, float angle, const Vector3<float>& t) noexcept;
	Matrix4& setRotationTranslation(const YawPitchRoll<float>& r, const Vector3<float>& t) noexcept;
	Matrix4& setRotationTranslation(const Euler<float>& e, const Vector3<float>& t) noexcept;
	Matrix4& setRotationTranslation(const Quaternion<float>& q, const Vector3<float>& t) noexcept;
	Matrix4& setScalingRotation(const Vector3<float>& s, Axis axis, float angle) noexcept { return setRotation(axis, angle).preScale(s); }
	Matrix4& setScalingRotation(const Vector3<float>& s, const Vector3<float>& axis, float angle) noexcept { return setRotation(axis, angle).preScale(s); }
	Matrix4& setScalingRotation(const Vector3<float>& s, const YawPitchRoll<float>& r) noexcept { return setRotation(r).preScale(s); }
	Matrix4& setScalingRotation(const Vector3<float>& s, const Euler<float>& e) noexcept { return setRotation(e).preScale(s); }
	Matrix4& setScalingRotation(const Vector3<float>& s, const Quaternion<float>& q) noexcept { return setRotation(q).preScale(s); }
	Matrix4& setScalingTranslation(const Vector3<float>& s, const Vector3<float>& t) noexcept;
	Matrix4& setScalingRotationTranslation(const Vector3<float>& s, Axis axis, float angle, const Vector3<float>& t) noexcept;
	Matrix4& setScalingRotationTranslation(const Vector3<float>& s, const Vector3<float>& axis, float angle, const Vector3<float>& t) noexcept;
	Matrix4& setScalingRotationTranslation(const Vector3<float>& s, const YawPitchRoll<float>& r, const Vector3<float>& t) noexcept;
	Matrix4& setScalingRotationTranslation(const Vector3<float>& s, const Euler<float>& e, const Vector3<float>& t) noexcept;
	Matrix4& setScalingRotationTranslation(const Vector3<float>& s, const Quaternion<float>& q, const Vector3<float>& t) noexcept;
	Matrix4& setShearing(float xy, float xz, float yx, float yz, float zx, float zy) noexcept;
	Matrix4& setReflection(const Plane<float>& plane) noexcept;
	Matrix4& setTranspose/*transposeOf*/(const Matrix4& m) noexcept;
	Matrix4& setInverse/*inverseOf*/(const Matrix4& m) noexcept;
	Matrix4& setInverseTranspose/*inverseTransposeOf*/(const Matrix4& m) noexcept { return setInverse(m).transpose(); }
	Matrix4& preConcatenate(const Matrix4& m) noexcept;
	Matrix4& concatenate(const Matrix4& m) noexcept { return operator*=(m); }
	Matrix4& translate(const Vector3<float>& v) noexcept;
	Matrix4& preScale(const Vector3<float>& v) noexcept;
	Matrix4& scale(const Vector3<float>& v) noexcept;
	//Matrix4& rotate(Axis axis, float angle) noexcept { concatenate(Matrix4(Uninitialized()).setRotation(axis, angle)); return *this; }
	//Matrix4& rotate(const Vector3<float>& axis, float angle) noexcept { concatenate(Matrix4(Uninitialized()).setRotation(axis, angle)); return *this; }
	//Matrix4& rotate(const YawPitchRoll<float>& r) noexcept { concatenate(Matrix4(Uninitialized()).setRotation(r)); return *this; }
	//Matrix4& rotate(const Euler<float>& e) noexcept { concatenate(Matrix4(Uninitialized()).setRotation(e)); return *this; }
	//Matrix4& rotate(const Quaternion<float>& q) noexcept { concatenate(Matrix4(Uninitialized()).setRotation(q)); return *this; }
	//Matrix4& shear(float xy, float xz, float yx, float yz, float zx, float zy) noexcept { concatenate(Matrix4(Uninitialized()).setShearing(xy, xz, yx, yz, zx, zy)); return *this; }
	Matrix4& negate() noexcept;
	Matrix4& transpose() noexcept;
	Matrix4& invert() noexcept { setInverse(*this); }

	//static const Matrix4& getZero() noexcept { return ZERO; }
	//static const Matrix4& getIdentity() noexcept { return IDENTITY; }

	static const Matrix4 ZERO;
	static const Matrix4 IDENTITY;

	union
	{
		simd::float4 row0;
		struct { float m00, m01, m02, m03; };
	};

	union
	{
		simd::float4 row1;
		struct { float m10, m11, m12, m13; };
	};

	union
	{
		simd::float4 row2;
		struct { float m20, m21, m22, m23; };
	};

	union
	{
		simd::float4 row3;
		struct { float m30, m31, m32, m33; };
	};
};

const Matrix4<float> Matrix4<float>::ZERO{};
const Matrix4<float> Matrix4<float>::IDENTITY{ 1.f, 0.f, 0.f, 0.f, 0.f, 1.f, 0.f, 0.f, 0.f, 0.f, 1.f, 0.f, 0.f, 0.f, 0.f, 1.f };

#endif /* SIMD_HAS_FLOAT4 */

template<typename T>
inline Matrix4<T>::Matrix4(T m00, T m01, T m02, T m03, T m10, T m11, T m12, T m13, T m20, T m21, T m22, T m23, T m30, T m31, T m32, T m33) :
	m00(m00), m01(m01), m02(m02), m03(m03), m10(m10), m11(m11), m12(m12), m13(m13), m20(m20), m21(m21), m22(m22), m23(m23), 
	m30(m30), m31(m31), m32(m32), m33(m33) 
{
}

template<typename T>
inline Matrix4<T>::Matrix4(const Vector4<T>& row0, const Vector4<T>& row1, const Vector4<T>& row2, const Vector4<T>& row3) :
	m00(row0.x), m01(row0.y), m02(row0.z), m03(row0.w), m10(row1.x), m11(row1.y), m12(row1.z), m13(row1.w),
	m20(row2.x), m21(row2.y), m22(row2.z), m23(row2.w), m30(row3.x), m31(row3.y), m32(row3.z), m33(row3.w) 
{
}

template<typename T>
inline Matrix4<T>::Matrix4(const std::tuple<Vector3<T>, Vector3<T>, Vector3<T>, Vector3<T>>& t) :
	Matrix4(std::get<0>(t), std::get<1>(t), std::get<2>(t), std::get<3>(t))
{
}

template<typename T>
inline Matrix4<T>::Matrix4(const Matrix2<T>& m) : 
	m00(m.m00), m01(m.m01), m02(), m03(), m10(m.m10), m11(m.m11), m12(), m13(), m20(), m21(), m22(1), m23(), 
	m30(), m31(), m32(), m33(1) 
{
}

template<typename T>
inline Matrix4<T>::Matrix4(const Matrix3<T>& m) : 
	m00(m.m00), m01(m.m01), m02(m.m02), m03(), m10(m.m10), m11(m.m11), m12(m.m12), m13(), m20(m.m20), m21(m.m21), m22(m.m22), m23(), 
	m30(), m31(), m32(), m33(1) 
{
}

template<typename T>
inline Matrix4<T>::Matrix4(const AffineTransform<T>& m) : 
	m00(m.m00), m01(m.m01), m02(m.m02), m03(), m10(m.m10), m11(m.m11), m12(m.m12), m13(),
	m20(m.m20), m21(m.m21), m22(m.m22), m23(), m30(m.x), m31(m.y), m32(m.z), m33(1) 
{
}

template<typename T>
inline Matrix4<T>::Matrix4(const T* m) :
	m00(m[0]), m01(m[1]), m02(m[2]), m03(m[3]), m10(m[4]), m11(m[5]), m12(m[6]), m13(m[7]), 
	m20(m[8]), m21(m[9]), m22(m[10]), m23(m[11]), m30(m[12]), m31(m[13]), m32(m[14]), m33(m[15])
{
}

template<typename T>
inline Matrix4<T> Matrix4<T>::operator-() const
{
	return Matrix4<T>(-m00, -m01, -m02, -m03,
		-m10, -m11, -m12, -m13,
		-m20, -m21, -m22, -m23,
		-m30, -m31, -m32, -m33);
}

template<typename T>
inline Matrix4<T>& Matrix4<T>::operator+=(const Matrix4<T>& m)
{
	m00 += m.m00; m01 += m.m01; m02 += m.m02; m03 += m.m03;
	m10 += m.m10; m11 += m.m11; m12 += m.m12; m13 += m.m13;
	m20 += m.m20; m21 += m.m21; m22 += m.m22; m23 += m.m23;
	m30 += m.m30; m31 += m.m31; m32 += m.m32; m33 += m.m33;
	return *this;
}

template<typename T>
inline Matrix4<T>& Matrix4<T>::operator-=(const Matrix4<T>& m)
{
	m00 -= m.m00; m01 -= m.m01; m02 -= m.m02; m03 -= m.m03;
	m10 -= m.m10; m11 -= m.m11; m12 -= m.m12; m13 -= m.m13;
	m20 -= m.m20; m21 -= m.m21; m22 -= m.m22; m23 -= m.m23;
	m30 -= m.m30; m31 -= m.m31; m32 -= m.m32; m33 -= m.m33;
	return *this;
}

template<typename T>
inline Matrix4<T>& Matrix4<T>::operator*=(T f)
{
	m00 *= f; m01 *= f; m02 *= f; m03 *= f;
	m10 *= f; m11 *= f; m12 *= f; m13 *= f;
	m20 *= f; m21 *= f; m22 *= f; m23 *= f;
	m30 *= f; m31 *= f; m32 *= f; m33 *= f;
	return *this;
}

template<typename T>
inline Matrix4<T>& Matrix4<T>::operator*=(const Matrix4<T>& m)
{
	set(m00*m.m00 + m01*m.m10 + m02*m.m20 + m03*m.m30,
		m00*m.m01 + m01*m.m11 + m02*m.m21 + m03*m.m31,
		m00*m.m02 + m01*m.m12 + m02*m.m22 + m03*m.m32,
		m00*m.m03 + m01*m.m13 + m02*m.m23 + m03*m.m33,
		m10*m.m00 + m11*m.m10 + m12*m.m20 + m13*m.m30,
		m10*m.m01 + m11*m.m11 + m12*m.m21 + m13*m.m31,
		m10*m.m02 + m11*m.m12 + m12*m.m22 + m13*m.m32,
		m10*m.m03 + m11*m.m13 + m12*m.m23 + m13*m.m33,
		m20*m.m00 + m21*m.m10 + m22*m.m20 + m23*m.m30,
		m20*m.m01 + m21*m.m11 + m22*m.m21 + m23*m.m31,
		m20*m.m02 + m21*m.m12 + m22*m.m22 + m23*m.m32,
		m20*m.m03 + m21*m.m13 + m22*m.m23 + m23*m.m33,
		m30*m.m00 + m31*m.m10 + m32*m.m20 + m33*m.m30,
		m30*m.m01 + m31*m.m11 + m32*m.m21 + m33*m.m31,
		m30*m.m02 + m31*m.m12 + m32*m.m22 + m33*m.m32,
		m30*m.m03 + m31*m.m13 + m32*m.m23 + m33*m.m33);
	return *this;
}

template<typename T>
	requires std::floating_point<T>
inline Matrix4<T> operator+(const Matrix4<T>& m1, const Matrix4<T>& m2) noexcept
{
	return Matrix4<T>(m1.m00 + m2.m00, m1.m01 + m2.m01, m1.m02 + m2.m02, m1.m03 + m2.m03,
		m1.m10 + m2.m10, m1.m11 + m2.m11, m1.m12 + m2.m12, m1.m13 + m2.m13,
		m1.m20 + m2.m20, m1.m21 + m2.m21, m1.m22 + m2.m22, m1.m23 + m2.m23,
		m1.m30 + m2.m30, m1.m31 + m2.m31, m1.m32 + m2.m32, m1.m33 + m2.m33);
}

template<typename T>
	requires std::floating_point<T>
inline Matrix4<T> operator-(const Matrix4<T>& m1, const Matrix4<T>& m2) noexcept
{
	return Matrix4<T>(m1.m00 - m2.m00, m1.m01 - m2.m01, m1.m02 - m2.m02, m1.m03 - m2.m03,
		m1.m10 - m2.m10, m1.m11 - m2.m11, m1.m12 - m2.m12, m1.m13 - m2.m13,
		m1.m20 - m2.m20, m1.m21 - m2.m21, m1.m22 - m2.m22, m1.m23 - m2.m23,
		m1.m30 - m2.m30, m1.m31 - m2.m31, m1.m32 - m2.m32, m1.m33 - m2.m33);
}

template<typename T>
	requires std::floating_point<T>
inline Matrix4<T> operator*(T f, const Matrix4<T>& m) noexcept
{
	return Matrix4<T>(f*m.m00, f*m.m01, f*m.m02, f*m.m03,
		f*m.m10, f*m.m11, f*m.m12, f*m.m13,
		f*m.m20, f*m.m21, f*m.m22, f*m.m23,
		f*m.m30, f*m.m31, f*m.m32, f*m.m33);
}

template<typename T>
	requires std::floating_point<T>
inline Matrix4<T> operator*(const Matrix4<T>& m, T f) noexcept
{
	return Matrix4<T>(m.m00*f, m.m01*f, m.m02*f, m.m03*f,
		m.m10*f, m.m11*f, m.m12*f, m.m13*f,
		m.m20*f, m.m21*f, m.m22*f, m.m23*f,
		m.m30*f, m.m31*f, m.m32*f, m.m33*f);
}

template<typename T>
	requires std::floating_point<T>
inline Matrix4<T> operator*(const Matrix4<T>& m1, const Matrix4<T>& m2) noexcept
{
	return Matrix4<T>(m1.m00*m2.m00 + m1.m01*m2.m10 + m1.m02*m2.m20 + m1.m03*m2.m30,
		m1.m00*m2.m01 + m1.m01*m2.m11 + m1.m02*m2.m21 + m1.m03*m2.m31,
		m1.m00*m2.m02 + m1.m01*m2.m12 + m1.m02*m2.m22 + m1.m03*m2.m32,
		m1.m00*m2.m03 + m1.m01*m2.m13 + m1.m02*m2.m23 + m1.m03*m2.m33,
		m1.m10*m2.m00 + m1.m11*m2.m10 + m1.m12*m2.m20 + m1.m13*m2.m30,
		m1.m10*m2.m01 + m1.m11*m2.m11 + m1.m12*m2.m21 + m1.m13*m2.m31,
		m1.m10*m2.m02 + m1.m11*m2.m12 + m1.m12*m2.m22 + m1.m13*m2.m32,
		m1.m10*m2.m03 + m1.m11*m2.m13 + m1.m12*m2.m23 + m1.m13*m2.m33,
		m1.m20*m2.m00 + m1.m21*m2.m10 + m1.m22*m2.m20 + m1.m23*m2.m30,
		m1.m20*m2.m01 + m1.m21*m2.m11 + m1.m22*m2.m21 + m1.m23*m2.m31,
		m1.m20*m2.m02 + m1.m21*m2.m12 + m1.m22*m2.m22 + m1.m23*m2.m32,
		m1.m20*m2.m03 + m1.m21*m2.m13 + m1.m22*m2.m23 + m1.m23*m2.m33,
		m1.m30*m2.m00 + m1.m31*m2.m10 + m1.m32*m2.m20 + m1.m33*m2.m30,
		m1.m30*m2.m01 + m1.m31*m2.m11 + m1.m32*m2.m21 + m1.m33*m2.m31,
		m1.m30*m2.m02 + m1.m31*m2.m12 + m1.m32*m2.m22 + m1.m33*m2.m32,
		m1.m30*m2.m03 + m1.m31*m2.m13 + m1.m32*m2.m23 + m1.m33*m2.m33);
}

template<typename T>
	requires std::floating_point<T>
inline Matrix4<T> operator/(const Matrix4<T>& m, T f) noexcept
{ 
	return operator*(m, T(1)/f); 
}

template<typename T>
inline bool Matrix4<T>::operator==(const Matrix4<T>& m) const
{
	return (m00 == m.m00) && (m01 == m.m01) && (m02 == m.m02) && (m03 == m.m03) &&
		(m10 == m.m10) && (m11 == m.m11) && (m12 == m.m12) && (m13 == m.m13) &&
		(m20 == m.m20) && (m21 == m.m21) && (m22 == m.m22) && (m23 == m.m23) &&
		(m30 == m.m30) && (m31 == m.m31) && (m32 == m.m32) && (m33 == m.m33);
}

template<typename C, typename T, typename U>
	requires std::floating_point<U>
inline std::basic_istream<C, T>& operator>>(std::basic_istream<C, T>& s, Matrix4<U>& m)
{
	return s >> m.m00 >> std::ws >> m.m01 >> std::ws >> m.m02 >> std::ws >> m.m03 >> std::ws >>
		m.m10 >> std::ws >> m.m11 >> std::ws >> m.m12 >> std::ws >> m.m13 >> std::ws >>
		m.m20 >> std::ws >> m.m21 >> std::ws >> m.m22 >> std::ws >> m.m23 >> std::ws >>
		m.m30 >> std::ws >> m.m31 >> std::ws >> m.m32 >> std::ws >> m.m33;
}

template<typename C, typename T, typename U>
	requires std::floating_point<U>
inline std::basic_ostream<C, T>& operator<<(std::basic_ostream<C, T>& s, const Matrix4<U>& m)
{
	constexpr C WS(0x20);
	return s << m.m00 << WS << m.m01 << WS << m.m02 << WS << m.m03 << WS <<
		m.m10 << WS << m.m11 << WS << m.m12 << WS << m.m13 << WS <<
		m.m20 << WS << m.m21 << WS << m.m22 << WS << m.m23 << WS <<
		m.m30 << WS << m.m31 << WS << m.m32 << WS << m.m33;
}

template<typename T>
inline bool Matrix4<T>::isZero() const
{
	return (m00 == T()) && (m01 == T()) && (m02 == T()) && (m03 == T()) &&
		(m10 == T()) && (m11 == T()) && (m12 == T()) && (m13 == T()) &&
		(m20 == T()) && (m21 == T()) && (m22 == T()) && (m23 == T()) &&
		(m30 == T()) && (m31 == T()) && (m32 == T()) && (m33 == T());
}

template<typename T>
inline bool Matrix4<T>::isApproxZero() const
{
	return (std::fabs(m00) < Constants<T>::TOLERANCE) &&
		(std::fabs(m01) < Constants<T>::TOLERANCE) &&
		(std::fabs(m02) < Constants<T>::TOLERANCE) &&
		(std::fabs(m03) < Constants<T>::TOLERANCE) &&
		(std::fabs(m10) < Constants<T>::TOLERANCE) &&
		(std::fabs(m11) < Constants<T>::TOLERANCE) &&
		(std::fabs(m12) < Constants<T>::TOLERANCE) &&
		(std::fabs(m13) < Constants<T>::TOLERANCE) &&
		(std::fabs(m20) < Constants<T>::TOLERANCE) &&
		(std::fabs(m21) < Constants<T>::TOLERANCE) &&
		(std::fabs(m22) < Constants<T>::TOLERANCE) &&
		(std::fabs(m23) < Constants<T>::TOLERANCE) &&
		(std::fabs(m30) < Constants<T>::TOLERANCE) &&
		(std::fabs(m31) < Constants<T>::TOLERANCE) &&
		(std::fabs(m32) < Constants<T>::TOLERANCE) &&
		(std::fabs(m33) < Constants<T>::TOLERANCE);
}

template<typename T>
inline bool Matrix4<T>::isIdentity() const
{
	return (m00 == T(1)) && (m01 == T()) && (m02 == T()) && (m03 == T()) &&
		(m10 == T()) && (m11 == T(1)) && (m12 == T()) && (m13 == T()) &&
		(m20 == T()) && (m21 == T()) && (m22 == T(1)) && (m23 == T()) &&
		(m30 == T()) && (m31 == T()) && (m32 == T()) && (m33 == T(1));
}

template<typename T>
inline bool Matrix4<T>::isApproxIdentity() const
{
	return (std::fabs(T(1) - m00) < Constants<T>::TOLERANCE) &&
		(std::fabs(m01) < Constants<T>::TOLERANCE) &&
		(std::fabs(m02) < Constants<T>::TOLERANCE) &&
		(std::fabs(m03) < Constants<T>::TOLERANCE) &&
		(std::fabs(m10) < Constants<T>::TOLERANCE) &&
		(std::fabs(T(1) - m11) < Constants<T>::TOLERANCE) &&
		(std::fabs(m12) < Constants<T>::TOLERANCE) &&
		(std::fabs(m13) < Constants<T>::TOLERANCE) &&
		(std::fabs(m20) < Constants<T>::TOLERANCE) &&
		(std::fabs(m21) < Constants<T>::TOLERANCE) &&
		(std::fabs(T(1) - m22) < Constants<T>::TOLERANCE) &&
		(std::fabs(m23) < Constants<T>::TOLERANCE) &&
		(std::fabs(m30) < Constants<T>::TOLERANCE) &&
		(std::fabs(m31) < Constants<T>::TOLERANCE) &&
		(std::fabs(m32) < Constants<T>::TOLERANCE) &&
		(std::fabs(T(1) - m33) < Constants<T>::TOLERANCE);
}

template<typename T>
inline bool Matrix4<T>::isApproxEqual(const Matrix4<T>& m) const
{
	return (std::fabs(m.m00 - m00) < Constants<T>::TOLERANCE) &&
		(std::fabs(m.m01 - m01) < Constants<T>::TOLERANCE) &&
		(std::fabs(m.m02 - m02) < Constants<T>::TOLERANCE) &&
		(std::fabs(m.m03 - m03) < Constants<T>::TOLERANCE) &&
		(std::fabs(m.m10 - m10) < Constants<T>::TOLERANCE) &&
		(std::fabs(m.m11 - m11) < Constants<T>::TOLERANCE) &&
		(std::fabs(m.m12 - m12) < Constants<T>::TOLERANCE) &&
		(std::fabs(m.m13 - m13) < Constants<T>::TOLERANCE) &&
		(std::fabs(m.m20 - m20) < Constants<T>::TOLERANCE) &&
		(std::fabs(m.m21 - m21) < Constants<T>::TOLERANCE) &&
		(std::fabs(m.m22 - m22) < Constants<T>::TOLERANCE) &&
		(std::fabs(m.m23 - m23) < Constants<T>::TOLERANCE) &&
		(std::fabs(m.m30 - m30) < Constants<T>::TOLERANCE) &&
		(std::fabs(m.m31 - m31) < Constants<T>::TOLERANCE) &&
		(std::fabs(m.m32 - m32) < Constants<T>::TOLERANCE) &&
		(std::fabs(m.m33 - m33) < Constants<T>::TOLERANCE);
}

template<typename T>
inline bool Matrix4<T>::isApproxEqual(const Matrix4<T>& m, T tolerance) const
{
	return (std::fabs(m.m00 - m00) < tolerance) &&
		(std::fabs(m.m01 - m01) < tolerance) &&
		(std::fabs(m.m02 - m02) < tolerance) &&
		(std::fabs(m.m03 - m03) < tolerance) &&
		(std::fabs(m.m10 - m10) < tolerance) &&
		(std::fabs(m.m11 - m11) < tolerance) &&
		(std::fabs(m.m12 - m12) < tolerance) &&
		(std::fabs(m.m13 - m13) < tolerance) &&
		(std::fabs(m.m20 - m20) < tolerance) &&
		(std::fabs(m.m21 - m21) < tolerance) &&
		(std::fabs(m.m22 - m22) < tolerance) &&
		(std::fabs(m.m23 - m23) < tolerance) &&
		(std::fabs(m.m30 - m30) < tolerance) &&
		(std::fabs(m.m31 - m31) < tolerance) &&
		(std::fabs(m.m32 - m32) < tolerance) &&
		(std::fabs(m.m33 - m33) < tolerance);
}

template<typename T>
inline bool Matrix4<T>::isFinite() const
{
	return std::isfinite(m00) && std::isfinite(m01) && std::isfinite(m02) && std::isfinite(m03) &&
		std::isfinite(m10) && std::isfinite(m11) && std::isfinite(m12) && std::isfinite(m13) &&
		std::isfinite(m20) && std::isfinite(m21) && std::isfinite(m22) && std::isfinite(m23) &&
		std::isfinite(m30) && std::isfinite(m31) && std::isfinite(m32) && std::isfinite(m33);
}

template<typename T>
inline T Matrix4<T>::getDeterminant() const
{
	auto det3 = [](T m00, T m01, T m02, T m10, T m11, T m12, T m20, T m21, T m22) -> T
		{
			return m00*(m11*m22 - m12*m21) - m01*(m10*m22 - m12*m20) + m02*(m10*m21 - m11*m20);
		};

	T q = det3(m11, m21, m31, m12, m22, m32, m13, m23, m33);
	T r = det3(m10, m20, m30, m12, m22, m32, m13, m23, m33);
	T s = det3(m10, m20, m30, m11, m21, m31, m13, m23, m33);
	T t = det3(m10, m20, m30, m11, m21, m31, m12, m22, m32);
	return (m00*q - m01*r + m02*s - m03*t);
}

template<typename T>
inline Matrix4<T>& Matrix4<T>::setZero()
{
	m00 = T(); m01 = T(); m02 = T(); m03 = T();
	m10 = T(); m11 = T(); m12 = T(); m13 = T();
	m20 = T(); m21 = T(); m22 = T(); m23 = T();
	m30 = T(); m31 = T(); m32 = T(); m33 = T();
	return *this;
}

template<typename T>
inline Matrix4<T>& Matrix4<T>::setIdentity()
{
	m00 = T(1); m01 = T(); m02 = T(); m03 = T();
	m10 = T(); m11 = T(1); m12 = T(); m13 = T();
	m20 = T(); m21 = T(); m22 = T(1); m23 = T();
	m30 = T(); m31 = T(); m32 = T(); m33 = T(1);
	return *this;
}

template<typename T>
inline Matrix4<T>& Matrix4<T>::set(const Vector4<T>& row0, const Vector4<T>& row1, const Vector4<T>& row2, const Vector4<T>& row3)
{
	m00 = row0.x; m01 = row0.y; m02 = row0.z; m03 = row0.w; 
	m10 = row1.x; m11 = row1.y; m12 = row1.z; m13 = row1.w; 
	m20 = row2.x; m21 = row2.y; m22 = row2.z; m23 = row2.w; 
	m30 = row3.x; m31 = row3.y; m32 = row3.z; m33 = row3.w; 
	return *this;
}

template<typename T>
inline Matrix4<T>& Matrix4<T>::set(T m00, T m01, T m02, T m03, T m10, T m11, T m12, T m13, T m20, T m21, T m22, T m23, 
	T m30, T m31, T m32, T m33)
{
	this->m00 = m00; this->m01 = m01; this->m02 = m02; this->m03 = m03;
	this->m10 = m10; this->m11 = m11; this->m12 = m12; this->m13 = m13; 
	this->m20 = m20; this->m21 = m21; this->m22 = m22; this->m23 = m23; 
	this->m30 = m30; this->m31 = m31; this->m32 = m32; this->m33 = m33; 
	return *this;
}

template<typename T>
inline Matrix4<T>& Matrix4<T>::setTranslation(const Vector3<T>& v)
{
	m00 = T(1); m01 = T(); m02 = T(); m03 = T();
	m10 = T(); m11 = T(1); m12 = T(); m13 = T();
	m20 = T(); m21 = T(); m22 = T(1); m23 = T();
	m30 = v.x; m31 = v.y; m32 = v.z; m33 = T(1);
	return *this;
}

template<typename T>
inline Matrix4<T>& Matrix4<T>::setScaling(const Vector3<T>& v)
{
	m00 = v.x; m01 = T(); m02 = T(); m03 = T();
	m10 = T(); m11 = v.y; m12 = T(); m13 = T();
	m20 = T(); m21 = T(); m22 = v.z; m23 = T();
	m30 = T(); m31 = T(); m32 = T(); m33 = T(1);
	return *this;
}

template<typename T>
inline Matrix4<T>& Matrix4<T>::setRotation(Axis axis, T angle)
{
	if (angle != T(0))
	{
		T sinAngle = std::sin(angle);
		T cosAngle = std::cos(angle);

		switch (axis)
		{
			case Axis::X:
				return set(T(1), T(), T(), T(), T(), cosAngle, sinAngle, T(), T(), -sinAngle, cosAngle, T(), T(), T(), T(), T(1));
			case Axis::Y:
				return set(cosAngle, T(), -sinAngle, T(), T(), T(1), T(), T(), sinAngle, T(), cosAngle, T(), T(), T(), T(), T(1));
			case Axis::Z:
				return set(cosAngle, sinAngle, T(), T(), -sinAngle, cosAngle, T(), T(), T(), T(), T(1), T(), T(), T(), T(), T(1));
			default:
				return setZero();
		}
	}

	return setIdentity();
}

template<typename T>
inline Matrix4<T>& Matrix4<T>::setRotation(const Vector3<T>& axis, T angle)
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

		return set(aa*mCos + cosAngle, abm + cSin, acm - bSin, T(),
			abm - cSin, bb*mCos + cosAngle, bcm + aSin, T(),
			acm + bSin, bcm - aSin, cc*mCos + cosAngle, T(),
			T(), T(), T(), T(1));
	}

	return setIdentity();
}

template<typename T>
inline Matrix4<T>& Matrix4<T>::setRotationTranslation(Axis axis, T angle, const Vector3<T>& t)
{
	setRotation(axis, angle);
	m30 = t.x; m31 = t.y; m32 = t.z;
	return *this;
}

template<typename T>
inline Matrix4<T>& Matrix4<T>::setRotationTranslation(const Vector3<T>& axis, T angle, const Vector3<T>& t)
{
	setRotation(axis, angle);
	m30 = t.x; m31 = t.y; m32 = t.z;
	return *this;
}

template<typename T>
inline Matrix4<T>& Matrix4<T>::setScalingRotation(const Vector3<T>& s, Axis axis, T angle)
{
	setRotation(axis, angle);
	m00 *= s.x; m01 *= s.x; m02 *= s.x;
	m10 *= s.y; m11 *= s.y; m12 *= s.y;
	m20 *= s.z; m21 *= s.z; m22 *= s.z;
	return *this;
}

template<typename T>
inline Matrix4<T>& Matrix4<T>::setScalingRotation(const Vector3<T>& s, const Vector3<T>& axis, T angle)
{
	setRotation(axis, angle);
	m00 *= s.x; m01 *= s.x; m02 *= s.x;
	m10 *= s.y; m11 *= s.y; m12 *= s.y;
	m20 *= s.z; m21 *= s.z; m22 *= s.z;
	return *this;
}

template<typename T>
inline Matrix4<T>& Matrix4<T>::setScalingTranslation(const Vector3<T>& s, const Vector3<T>& t)
{
	m00 = s.x; m01 = T(); m02 = T(); m03 = T();
	m10 = T(); m11 = s.y; m12 = T(); m13 = T();
	m20 = T(); m21 = T(); m22 = s.z; m23 = T();
	m30 = t.x; m31 = t.y; m32 = t.z; m33 = T(1);
	return *this;
}

template<typename T>
inline Matrix4<T>& Matrix4<T>::setScalingRotationTranslation(const Vector3<T>& s, Axis axis, T angle, const Vector3<T>& t)
{
	setRotation(axis, angle);
	m00 *= s.x; m01 *= s.x; m02 *= s.x;
	m10 *= s.y; m11 *= s.y; m12 *= s.y;
	m20 *= s.z; m21 *= s.z; m22 *= s.z;
	m30 = t.x; m31 = t.y; m32 = t.z;
	return *this;
}

template<typename T>
inline Matrix4<T>& Matrix4<T>::setScalingRotationTranslation(const Vector3<T>& s, const Vector3<T>& axis, T angle, const Vector3<T>& t)
{
	setRotation(axis, angle);
	m00 *= s.x; m01 *= s.x; m02 *= s.x;
	m10 *= s.y; m11 *= s.y; m12 *= s.y;
	m20 *= s.z; m21 *= s.z; m22 *= s.z;
	m30 = t.x; m31 = t.y; m32 = t.z;
	return *this;
}

template<typename T>
inline Matrix4<T>& Matrix4<T>::setShearing(T xy, T xz, T yx, T yz, T zx, T zy)
{
	m00 = T(1); m01 = xy; m02 = xz; m03 = T();
	m10 = yx; m11 = T(1); m12 = yz; m13 = T();
	m20 = zx; m21 = zy; m22 = T(1); m23 = T();
	m30 = T(); m31 = T(); m32 = T(); m33 = T(1);
	return *this;
}

template<typename T>
inline Matrix4<T>& Matrix4<T>::setTranspose(const Matrix4<T>& m)
{
	m00 = m.m00; m01 = m.m10; m02 = m.m20; m03 = m.m30;
	m10 = m.m01; m11 = m.m11; m12 = m.m21; m13 = m.m31;
	m20 = m.m02; m21 = m.m12; m22 = m.m22; m23 = m.m32;
	m30 = m.m03; m31 = m.m13; m32 = m.m23; m33 = m.m33;
	return *this;
}

template<typename T>
inline Matrix4<T>& Matrix4<T>::setInverse(const Matrix4<T>& m)
{
	T di = m.m00;
	T d = T(1)/di;
	T n00 = d;
	T n10 = m.m10*(-d);
	T n20 = m.m20*(-d);
	T n30 = m.m30*(-d);
	T n01 = m.m01*d;
	T n02 = m.m02*d;
	T n03 = m.m03*d;
	T n11 = m.m11 + n10*n01*di;
	T n12 = m.m12 + n10*n02*di;
	T n13 = m.m13 + n10*n03*di;
	T n21 = m.m21 + n20*n01*di;
	T n22 = m.m22 + n20*n02*di;
	T n23 = m.m23 + n20*n03*di;
	T n31 = m.m31 + n30*n01*di;
	T n32 = m.m32 + n30*n02*di;
	T n33 = m.m33 + n30*n03*di;
	di = n11;
	d = T(1)/di;
	n11 = d;
	n01 *= -d;
	n21 *= -d;
	n31 *= -d;
	n10 *= d;
	n12 *= d;
	n13 *= d;
	n00 += n01*n10*di;
	n02 += n01*n12*di;
	n03 += n01*n13*di;
	n20 += n21*n10*di;
	n22 += n21*n12*di;
	n23 += n21*n13*di;
	n30 += n31*n10*di;
	n32 += n31*n12*di;
	n33 += n31*n13*di;
	di = n22;
	d = T(1)/di;
	n22 = d;
	n02 *= -d;
	n12 *= -d;
	n32 *= -d;
	n20 *= d;
	n21 *= d;
	n23 *= d;
	n00 += n02*n20*di;
	n01 += n02*n21*di;
	n03 += n02*n23*di;
	n10 += n12*n20*di;
	n11 += n12*n21*di;
	n13 += n12*n23*di;
	n30 += n32*n20*di;
	n31 += n32*n21*di;
	n33 += n32*n23*di;
	di = n33;
	d = T(1)/di;
	n33 = d;
	n03 *= -d;
	n13 *= -d;
	n23 *= -d;
	n30 *= d;
	n31 *= d;
	n32 *= d;
	m00 = n00 + n03*n30*di;
	m01 = n01 + n03*n31*di;
	m02 = n02 + n03*n32*di;
	m03 = n03;
	m10 = n10 + n13*n30*di;
	m11 = n11 + n13*n31*di;
	m12 = n12 + n13*n32*di;
	m13 = n13;
	m20 = n20 + n23*n30*di;
	m21 = n21 + n23*n31*di;
	m22 = n22 + n23*n32*di;
	m23 = n23;
	m30 = n30;
	m31 = n31;
	m32 = n32;
	m33 = n33;
	return *this;
}

template<typename T>
inline Matrix4<T>& Matrix4<T>::preConcatenate(const Matrix4<T>& m)
{
	set(m.m00*m00 + m.m01*m10 + m.m02*m20 + m.m03*m30,
		m.m00*m01 + m.m01*m11 + m.m02*m21 + m.m03*m31,
		m.m00*m02 + m.m01*m12 + m.m02*m22 + m.m03*m32,
		m.m00*m03 + m.m01*m13 + m.m02*m23 + m.m03*m33,
		m.m10*m00 + m.m11*m10 + m.m12*m20 + m.m13*m30,
		m.m10*m01 + m.m11*m11 + m.m12*m21 + m.m13*m31,
		m.m10*m02 + m.m11*m12 + m.m12*m22 + m.m13*m32,
		m.m10*m03 + m.m11*m13 + m.m12*m23 + m.m13*m33,
		m.m20*m00 + m.m21*m10 + m.m22*m20 + m.m23*m30,
		m.m20*m01 + m.m21*m11 + m.m22*m21 + m.m23*m31,
		m.m20*m02 + m.m21*m12 + m.m22*m22 + m.m23*m32,
		m.m20*m03 + m.m21*m13 + m.m22*m23 + m.m23*m33,
		m.m30*m00 + m.m31*m10 + m.m32*m20 + m.m33*m30,
		m.m30*m01 + m.m31*m11 + m.m32*m21 + m.m33*m31,
		m.m30*m02 + m.m31*m12 + m.m32*m22 + m.m33*m32,
		m.m30*m03 + m.m31*m13 + m.m32*m23 + m.m33*m33);
	return *this;
}

template<typename T>
inline Matrix4<T>& Matrix4<T>::preScale(const Vector3<T>& v)
{
	m00 *= v.x; m01 *= v.x; m02 *= v.x; m03 *= v.x;
	m10 *= v.y; m11 *= v.y; m12 *= v.y; m13 *= v.y;
	m20 *= v.z; m21 *= v.z; m22 *= v.z; m23 *= v.z;
	return *this;
}

template<typename T>
inline Matrix4<T>& Matrix4<T>::scale(const Vector3<T>& v)
{
	m00 *= v.x; m01 *= v.y; m02 *= v.z;
	m10 *= v.x; m11 *= v.y; m12 *= v.z;
	m20 *= v.x; m21 *= v.y; m22 *= v.z;
	m30 *= v.x; m31 *= v.y; m32 *= v.z;
	return *this;
}

template<typename T>
inline Matrix4<T>& Matrix4<T>::negate()
{
	m00 = -m00; m01 = -m01; m02 = -m02; m03 = -m03;
	m10 = -m10; m11 = -m11; m12 = -m12; m13 = -m13;
	m20 = -m20; m21 = -m21; m22 = -m22; m23 = -m23;
	m30 = -m30; m31 = -m31; m32 = -m32; m33 = -m33;
	return *this;
}

template<typename T>
inline Matrix4<T>& Matrix4<T>::transpose()
{
	T t = m01; m01 = m10; m10 = t;
	t = m02; m02 = m20; m20 = t;
	t = m03; m03 = m30; m30 = t;
	t = m12; m12 = m21; m21 = t;
	t = m13; m13 = m31; m31 = t;
	t = m23; m23 = m32; m32 = t;
	return *this;
}

#if SIMD_HAS_FLOAT4

inline Matrix4<float>::Matrix4() : 
	row0(simd::zero<simd::float4>()), 
	row1(simd::zero<simd::float4>()), 
	row2(simd::zero<simd::float4>()), 
	row3(simd::zero<simd::float4>()) 
{
}

inline Matrix4<float>::Matrix4(Identity) :
	row0(Vector4<float>::UNIT_X),
	row1(Vector4<float>::UNIT_Y),
	row2(Vector4<float>::UNIT_Z),
	row3(Vector4<float>::UNIT_W)
{
}

inline Matrix4<float>::Matrix4(float m00, float m01, float m02, float m03, float m10, float m11, float m12, float m13,
	float m20, float m21, float m22, float m23, float m30, float m31, float m32, float m33) :
	row0(simd::set4(m00, m01, m02, m03)),
	row1(simd::set4(m10, m11, m12, m13)),
	row2(simd::set4(m20, m21, m22, m23)),
	row3(simd::set4(m30, m31, m32, m33))
{
}

inline Matrix4<float>::Matrix4(const Vector4<float>& row0, const Vector4<float>& row1, const Vector4<float>& row2, 
	const Vector4<float>& row3) :
	row0(row0),
	row1(row1),
	row2(row2),
	row3(row3)
{
}

inline Matrix4<float>::Matrix4(const std::tuple<Vector3<float>, Vector3<float>, Vector3<float>, Vector3<float>>& t) :
	Matrix4(std::get<0>(t), std::get<1>(t), std::get<2>(t), std::get<3>(t))
{
}

inline Matrix4<float>::Matrix4(const Matrix2<float>& m) : 
	row0(simd::cutoff2(m.row0)),
	row1(simd::cutoff2(m.row1)),
	row2(Vector4<float>::UNIT_Z),
	row3(Vector4<float>::UNIT_W)
{
}

inline Matrix4<float>::Matrix4(const Matrix3<float>& m) :
	row0(simd::cutoff3(m.row0)),
	row1(simd::cutoff3(m.row1)),
	row2(simd::cutoff3(m.row2)),
	row3(Vector4<float>::UNIT_W)
{
}

inline Matrix4<float>::Matrix4(const AffineTransform<float>& m) :
	row0(simd::cutoff3(m.row0)),
	row1(simd::cutoff3(m.row1)),
	row2(simd::cutoff3(m.row2)),
	row3(simd::insert3(m.row3, Vector4<float>::UNIT_W))
{
}

inline Matrix4<float>::Matrix4(const float* m) : 
	row0(simd::load4(m)),
	row1(simd::load4(m + 4)),
	row2(simd::load4(m + 8)),
	row3(simd::load4(m + 12))
{
}

inline Matrix4<float>::Matrix4(simd::float4 row0, simd::float4 row1, simd::float4 row2, simd::float4 row3) :
	row0(row0), 
	row1(row1), 
	row2(row2), 
	row3(row3) 
{
}

inline Matrix4<float>::Matrix4(const std::tuple<simd::float4, simd::float4, simd::float4, simd::float4>& t) :
	row0(std::get<0>(t)), 
	row1(std::get<1>(t)), 
	row2(std::get<2>(t)), 
	row3(std::get<3>(t))
{
}

inline Matrix4<float> Matrix4<float>::operator-() const
{
	return Matrix4<float>(simd::neg4(row0), simd::neg4(row1), simd::neg4(row2), simd::neg4(row3));
}

inline Matrix4<float>& Matrix4<float>::operator+=(const Matrix4<float>& m)
{
	row0 = simd::add4(row0, m.row0);
	row1 = simd::add4(row1, m.row1);
	row2 = simd::add4(row2, m.row2);
	row3 = simd::add4(row3, m.row3);
	return *this;
}

inline Matrix4<float>& Matrix4<float>::operator-=(const Matrix4<float>& m)
{
	row0 = simd::sub4(row0, m.row0);
	row1 = simd::sub4(row1, m.row1);
	row2 = simd::sub4(row2, m.row2);
	row3 = simd::sub4(row3, m.row3);
	return *this;
}

inline Matrix4<float>& Matrix4<float>::operator*=(float f)
{
	auto t = simd::set4(f);
	row0 = simd::mul4(row0, t);
	row1 = simd::mul4(row1, t);
	row2 = simd::mul4(row2, t);
	row3 = simd::mul4(row3, t);
	return *this;
}

inline Matrix4<float>& Matrix4<float>::operator*=(const Matrix4<float>& m)
{
	auto r0 = simd::mul4(simd::broadcast<0>(row0), m.row0);
	r0 = simd::add4(r0, simd::mul4(simd::broadcast<1>(row0), m.row1));
	r0 = simd::add4(r0, simd::mul4(simd::broadcast<2>(row0), m.row2));
	r0 = simd::add4(r0, simd::mul4(simd::broadcast<3>(row0), m.row3));
	auto r1 = simd::mul4(simd::broadcast<0>(row1), m.row0);
	r1 = simd::add4(r1, simd::mul4(simd::broadcast<1>(row1), m.row1));
	r1 = simd::add4(r1, simd::mul4(simd::broadcast<2>(row1), m.row2));
	r1 = simd::add4(r1, simd::mul4(simd::broadcast<3>(row1), m.row3));
	auto r2 = simd::mul4(simd::broadcast<0>(row2), m.row0);
	r2 = simd::add4(r2, simd::mul4(simd::broadcast<1>(row2), m.row1));
	r2 = simd::add4(r2, simd::mul4(simd::broadcast<2>(row2), m.row2));
	r2 = simd::add4(r2, simd::mul4(simd::broadcast<3>(row2), m.row3));
	auto r3 = simd::mul4(simd::broadcast<0>(row3), m.row0);
	r3 = simd::add4(r3, simd::mul4(simd::broadcast<1>(row3), m.row1));
	r3 = simd::add4(r3, simd::mul4(simd::broadcast<2>(row3), m.row2));
	r3 = simd::add4(r3, simd::mul4(simd::broadcast<3>(row3), m.row3));
	row0 = r0;
	row1 = r1;
	row2 = r2;
	row3 = r3;
	return *this;
}

template<>
inline Matrix4<float> operator+(const Matrix4<float>& m1, const Matrix4<float>& m2) noexcept
{
	return Matrix4<float>(simd::add4(m1.row0, m2.row0),
		simd::add4(m1.row1, m2.row1),
		simd::add4(m1.row2, m2.row2),
		simd::add4(m1.row3, m2.row3));
}

template<>
inline Matrix4<float> operator-(const Matrix4<float>& m1, const Matrix4<float>& m2) noexcept
{
	return Matrix4<float>(simd::sub4(m1.row0, m2.row0),
		simd::sub4(m1.row1, m2.row1),
		simd::sub4(m1.row2, m2.row2),
		simd::sub4(m1.row3, m2.row3));
}

template<>
inline Matrix4<float> operator*(float f, const Matrix4<float>& m) noexcept
{
	auto t = simd::set4(f);
	return Matrix4<float>(simd::mul4(t, m.row0),
		simd::mul4(t, m.row1),
		simd::mul4(t, m.row2),
		simd::mul4(t, m.row3));
}

template<>
inline Matrix4<float> operator*(const Matrix4<float>& m, float f) noexcept
{
	auto t = simd::set4(f);
	return Matrix4<float>(simd::mul4(m.row0, t),
		simd::mul4(m.row1, t),
		simd::mul4(m.row2, t),
		simd::mul4(m.row3, t));
}

template<>
inline Matrix4<float> operator*(const Matrix4<float>& m1, const Matrix4<float>& m2) noexcept
{
	auto r0 = simd::mul4(simd::broadcast<0>(m1.row0), m2.row0);
	r0 = simd::add4(r0, simd::mul4(simd::broadcast<1>(m1.row0), m2.row1));
	r0 = simd::add4(r0, simd::mul4(simd::broadcast<2>(m1.row0), m2.row2));
	r0 = simd::add4(r0, simd::mul4(simd::broadcast<3>(m1.row0), m2.row3));
	auto r1 = simd::mul4(simd::broadcast<0>(m1.row1), m2.row0);
	r1 = simd::add4(r1, simd::mul4(simd::broadcast<1>(m1.row1), m2.row1));
	r1 = simd::add4(r1, simd::mul4(simd::broadcast<2>(m1.row1), m2.row2));
	r1 = simd::add4(r1, simd::mul4(simd::broadcast<3>(m1.row1), m2.row3));
	auto r2 = simd::mul4(simd::broadcast<0>(m1.row2), m2.row0);
	r2 = simd::add4(r2, simd::mul4(simd::broadcast<1>(m1.row2), m2.row1));
	r2 = simd::add4(r2, simd::mul4(simd::broadcast<2>(m1.row2), m2.row2));
	r2 = simd::add4(r2, simd::mul4(simd::broadcast<3>(m1.row2), m2.row3));
	auto r3 = simd::mul4(simd::broadcast<0>(m1.row3), m2.row0);
	r3 = simd::add4(r3, simd::mul4(simd::broadcast<1>(m1.row3), m2.row1));
	r3 = simd::add4(r3, simd::mul4(simd::broadcast<2>(m1.row3), m2.row2));
	r3 = simd::add4(r3, simd::mul4(simd::broadcast<3>(m1.row3), m2.row3));
	return Matrix4<float>(r0, r1, r2, r3);
}

template<>
inline Matrix4<float> operator/(const Matrix4<float>& m, float f) noexcept 
{ 
	return operator*(m, 1.f/f); 
}

inline bool Matrix4<float>::operator==(const Matrix4<float>& m) const
{
	return simd::all4(simd::equal(row0, m.row0)) &&
		simd::all4(simd::equal(row1, m.row1)) &&
		simd::all4(simd::equal(row2, m.row2)) &&
		simd::all4(simd::equal(row3, m.row3));
}

template<typename C, typename T>
inline std::basic_istream<C, T>& operator>>(std::basic_istream<C, T>& s, Matrix4<float>& m)
{
	float m00, m01, m02, m03;
	float m10, m11, m12, m13;
	float m20, m21, m22, m23;
	float m30, m31, m32, m33;
	s >> m00 >> std::ws >> m01 >> std::ws >> m02 >> std::ws >> m03 >> std::ws >>
		m10 >> std::ws >> m11 >> std::ws >> m12 >> std::ws >> m13 >> std::ws >>
		m20 >> std::ws >> m21 >> std::ws >> m22 >> std::ws >> m23 >> std::ws >>
		m30 >> std::ws >> m31 >> std::ws >> m32 >> std::ws >> m33;
	m.set(m00, m01, m02, m03, m10, m11, m12, m13, m20, m21, m22, m23, m30, m31, m32, m33);
	return s;
}

template<typename A>
inline void Matrix4<float>::load(A& ar)
{
	float t00, t01, t02, t03;
	float t10, t11, t12, t13;
	float t20, t21, t22, t23;
	float t30, t31, t32, t33;
	ar(t00, t01, t02, t03, t10, t11, t12, t13, t20, t21, t22, t23, t30, t31, t32, t33);
	set(t00, t01, t02, t03, t10, t11, t12, t13, t20, t21, t22, t23, t30, t31, t32, t33);
}

inline bool Matrix4<float>::isZero() const
{
	const auto zero = simd::zero<simd::float4>();
	return simd::all4(simd::equal(row0, zero)) && 
		simd::all4(simd::equal(row1, zero)) && 
		simd::all4(simd::equal(row2, zero)) && 
		simd::all4(simd::equal(row3, zero));
}

inline bool Matrix4<float>::isApproxZero() const
{
	return simd::all4(simd::lessThan(simd::abs4(row0), Vector4<float>::TOLERANCE)) &&
		simd::all4(simd::lessThan(simd::abs4(row1), Vector4<float>::TOLERANCE)) &&
		simd::all4(simd::lessThan(simd::abs4(row2), Vector4<float>::TOLERANCE)) &&
		simd::all4(simd::lessThan(simd::abs4(row3), Vector4<float>::TOLERANCE));
}

inline bool Matrix4<float>::isIdentity() const
{
	return simd::all4(simd::equal(row0, Vector4<float>::UNIT_X)) && 
		simd::all4(simd::equal(row1, Vector4<float>::UNIT_Y)) &&
		simd::all4(simd::equal(row2, Vector4<float>::UNIT_Z)) &&
		simd::all4(simd::equal(row3, Vector4<float>::UNIT_W));
}

inline bool Matrix4<float>::isApproxIdentity() const
{
	return simd::all4(simd::lessThan(simd::abs4(simd::sub4(row0, Vector4<float>::UNIT_X)), Vector4<float>::TOLERANCE)) &&
		simd::all4(simd::lessThan(simd::abs4(simd::sub4(row1, Vector4<float>::UNIT_Y)), Vector4<float>::TOLERANCE)) &&
		simd::all4(simd::lessThan(simd::abs4(simd::sub4(row2, Vector4<float>::UNIT_Z)), Vector4<float>::TOLERANCE)) &&
		simd::all4(simd::lessThan(simd::abs4(simd::sub4(row3, Vector4<float>::UNIT_W)), Vector4<float>::TOLERANCE));
}

inline bool Matrix4<float>::isApproxEqual(const Matrix4& m) const
{
	return simd::all4(simd::lessThan(simd::abs4(simd::sub4(row0, m.row0)), Vector4<float>::TOLERANCE)) &&
		simd::all4(simd::lessThan(simd::abs4(simd::sub4(row1, m.row1)), Vector4<float>::TOLERANCE)) &&
		simd::all4(simd::lessThan(simd::abs4(simd::sub4(row2, m.row2)), Vector4<float>::TOLERANCE)) &&
		simd::all4(simd::lessThan(simd::abs4(simd::sub4(row3, m.row3)), Vector4<float>::TOLERANCE));
}

inline bool Matrix4<float>::isApproxEqual(const Matrix4& m, float tolerance) const
{
	auto t = simd::set4(tolerance);
	return simd::all4(simd::lessThan(simd::abs4(simd::sub4(row0, m.row0)), t)) &&
		simd::all4(simd::lessThan(simd::abs4(simd::sub4(row1, m.row1)), t)) &&
		simd::all4(simd::lessThan(simd::abs4(simd::sub4(row2, m.row2)), t)) &&
		simd::all4(simd::lessThan(simd::abs4(simd::sub4(row3, m.row3)), t));
}

inline bool Matrix4<float>::isFinite() const
{
	return simd::all4(simd::isFinite(row0)) &&
		simd::all4(simd::isFinite(row1)) &&
		simd::all4(simd::isFinite(row2)) &&
		simd::all4(simd::isFinite(row3));
}

inline float Matrix4<float>::getDeterminant() const // #TODO SIMD
{
	auto det3 = [](float m00, float m01, float m02, float m10, float m11, float m12, float m20, float m21, float m22) -> float
		{
			return m00*(m11*m22 - m12*m21) - m01*(m10*m22 - m12*m20) + m02*(m10*m21 - m11*m20);
		};

	float q = det3(m11, m21, m31, m12, m22, m32, m13, m23, m33);
	float r = det3(m10, m20, m30, m12, m22, m32, m13, m23, m33);
	float s = det3(m10, m20, m30, m11, m21, m31, m13, m23, m33);
	float t = det3(m10, m20, m30, m11, m21, m31, m12, m22, m32);
	return (m00*q - m01*r + m02*s - m03*t);
}

inline Matrix4<float>& Matrix4<float>::setZero()
{
	row0 = simd::zero<simd::float4>();
	row1 = simd::zero<simd::float4>();
	row2 = simd::zero<simd::float4>();
	row3 = simd::zero<simd::float4>();
	return *this;
}

inline Matrix4<float>& Matrix4<float>::setIdentity()
{
	row0 = Vector4<float>::UNIT_X;
	row1 = Vector4<float>::UNIT_Y;
	row2 = Vector4<float>::UNIT_Z;
	row3 = Vector4<float>::UNIT_W;
	return *this;
}

inline Matrix4<float>& Matrix4<float>::set(const Vector4<float>& row0, const Vector4<float>& row1, const Vector4<float>& row2,
	const Vector4<float>& row3)
{
	this->row0 = row0;
	this->row1 = row1;
	this->row2 = row2;
	this->row3 = row3;
	return *this;
}

inline Matrix4<float>& Matrix4<float>::set(simd::float4 row0, simd::float4 row1, simd::float4 row2, simd::float4 row3)
{
	this->row0 = row0;
	this->row1 = row1;
	this->row2 = row2;
	this->row3 = row3;
	return *this;
}

inline Matrix4<float>& Matrix4<float>::set(float m00, float m01, float m02, float m03, float m10, float m11, float m12, float m13,
	float m20, float m21, float m22, float m23, float m30, float m31, float m32, float m33)
{
	row0 = simd::set4(m00, m01, m02, m03);
	row1 = simd::set4(m10, m11, m12, m13);
	row2 = simd::set4(m20, m21, m22, m23);
	row3 = simd::set4(m30, m31, m32, m33);
	return *this;
}

inline Matrix4<float>& Matrix4<float>::setTranslation(const Vector3<float>& v)
{
	row0 = Vector4<float>::UNIT_X;
	row1 = Vector4<float>::UNIT_Y;
	row2 = Vector4<float>::UNIT_Z;
	row3 = simd::insert3(v, Vector4<float>::UNIT_W);
	return *this;
}

inline Matrix4<float>& Matrix4<float>::setScaling(const Vector3<float>& v)
{
	row0 = simd::cutoff1(v);
	row1 = simd::and4(v, simd::constant4i<simd::float4, 0, -1, 0, 0>());
	row2 = simd::and4(v, simd::constant4i<simd::float4, 0, 0, -1, 0>());
	row3 = Vector4<float>::UNIT_W;
	return *this;
}

inline Matrix4<float>& Matrix4<float>::setRotation(Axis axis, float angle)
{
	if (angle != 0.f)
	{
		float sinAngle = std::sin(angle);
		float cosAngle = std::cos(angle);

		switch (axis)
		{
			case Axis::X:
				return set(1.f, 0.f, 0.f, 0.f, 0.f, cosAngle, sinAngle, 0.f, 0.f, -sinAngle, cosAngle, 0.f, 0.f, 0.f, 0.f, 1.f);
			case Axis::Y:
				return set(cosAngle, 0.f, -sinAngle, 0.f, 0.f, 1.f, 0.f, 0.f, sinAngle, 0.f, cosAngle, 0.f, 0.f, 0.f, 0.f, 1.f);
			case Axis::Z:
				return set(cosAngle, sinAngle, 0.f, 0.f, -sinAngle, cosAngle, 0.f, 0.f, 0.f, 0.f, 1.f, 0.f, 0.f, 0.f, 0.f, 1.f);
			default:
				return setZero();
		}
	}

	return setIdentity();
}

inline Matrix4<float>& Matrix4<float>::setRotation(const Vector3<float>& axis, float angle)
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

		return set(aa*mCos + cosAngle, abm + cSin, acm - bSin, 0.f,
			abm - cSin, bb*mCos + cosAngle, bcm + aSin, 0.f,
			acm + bSin, bcm - aSin, cc*mCos + cosAngle, 0.f,
			0.f, 0.f, 0.f, 1.f);
	}

	return setIdentity();
}

inline Matrix4<float>& Matrix4<float>::setRotationTranslation(Axis axis, float angle, const Vector3<float>& t)
{
	setRotation(axis, angle);
	row3 = simd::insert3(t, Vector4<float>::UNIT_W);
	return *this;
}

inline Matrix4<float>& Matrix4<float>::setRotationTranslation(const Vector3<float>& axis, float angle, const Vector3<float>& t)
{
	setRotation(axis, angle);
	row3 = simd::insert3(t, Vector4<float>::UNIT_W);
	return *this;
}

inline Matrix4<float>& Matrix4<float>::setScalingTranslation(const Vector3<float>& s, const Vector3<float>& t)
{
	row0 = simd::cutoff1(s);
	row1 = simd::and4(s, simd::constant4i<simd::float4, 0, -1, 0, 0>());
	row2 = simd::and4(s, simd::constant4i<simd::float4, 0, 0, -1, 0>());
	row3 = simd::insert3(t, Vector4<float>::UNIT_W);
	return *this;
}

inline Matrix4<float>& Matrix4<float>::setScalingRotationTranslation(const Vector3<float>& s, Axis axis, float angle,
	const Vector3<float>& t)
{
	setRotation(axis, angle);
	preScale(s);
	row3 = simd::insert3(t, Vector4<float>::UNIT_W);
	return *this;
}

inline Matrix4<float>& Matrix4<float>::setScalingRotationTranslation(const Vector3<float>& s, const Vector3<float>& axis, float angle,
	const Vector3<float>& t)
{
	setRotation(axis, angle);
	preScale(s);
	row3 = simd::insert3(t, Vector4<float>::UNIT_W);
	return *this;
}

inline Matrix4<float>& Matrix4<float>::setShearing(float xy, float xz, float yx, float yz, float zx, float zy)
{
	row0 = simd::set3(1.f, xy, xz);
	row1 = simd::set3(yx, 1.f, yz);
	row2 = simd::set3(zx, zy, 1.f);
	row3 = Vector4<float>::UNIT_W;
	return *this;
}

inline Matrix4<float>& Matrix4<float>::setTranspose(const Matrix4<float>& m)
{
	//simd::transpose4x4(m.row0, m.row1, m.row2, m.row3, row0, row1, row2, row3);
	std::tie(row0, row1, row2, row3) = simd::transpose4x4(m.row0, m.row1, m.row2, m.row3);
	return *this;
}

inline Matrix4<float>& Matrix4<float>::setInverse(const Matrix4<float>& m)
{
	auto r2 = m.row2;
	auto t = simd::mul4(r2, m.row3);
	t = simd::swizzle<1, 0, 3, 2>(t);
	auto minor0 = simd::mul4(m.row1, t);
	auto minor1 = simd::mul4(m.row0, t);
	t = simd::swizzle<2, 3, 0, 1>(t);
	minor0 = simd::sub4(simd::mul4(m.row1, t), minor0);
	minor1 = simd::sub4(simd::mul4(m.row0, t), minor1);
	minor1 = simd::swizzle<2, 3, 0, 1>(minor1);
	t = simd::mul4(m.row1, r2);
	t = simd::swizzle<1, 0, 3, 2>(t);
	minor0 = simd::add4(simd::mul4(m.row3, t), minor0);
	auto minor3 = simd::mul4(m.row0, t);
	t = simd::swizzle<2, 3, 0, 1>(t);
	minor0 = simd::sub4(minor0, simd::mul4(m.row3, t));
	minor3 = simd::sub4(simd::mul4(m.row0, t), minor3);
	minor3 = simd::swizzle<2, 3, 0, 1>(minor3);
	t = simd::mul4(simd::swizzle<2, 3, 0, 1>(m.row1), m.row3);
	t = simd::swizzle<1, 0, 3, 2>(t);
	r2 = simd::swizzle<2, 3, 0, 1>(r2);
	minor0 = simd::add4(simd::mul4(r2, t), minor0);
	auto minor2 = simd::mul4(m.row0, t);
	t = simd::swizzle<2, 3, 0, 1>(t);
	minor0 = simd::sub4(minor0, simd::mul4(r2, t));
	minor2 = simd::sub4(simd::mul4(m.row0, t), minor2);
	minor2 = simd::swizzle<2, 3, 0, 1>(minor2);
	t = simd::mul4(m.row0, m.row1);
	t = simd::swizzle<1, 0, 3, 2>(t);
	minor2 = simd::add4(simd::mul4(m.row3, t), minor2);
	minor3 = simd::sub4(simd::mul4(r2, t), minor3);
	t = simd::swizzle<2, 3, 0, 1>(t);
	minor2 = simd::sub4(simd::mul4(m.row3, t), minor2);
	minor3 = simd::sub4(minor3, simd::mul4(r2, t));
	t = simd::mul4(m.row0, m.row3);
	t = simd::swizzle<1, 0, 3, 2>(t);
	minor1 = simd::sub4(minor1, simd::mul4(r2, t));
	minor2 = simd::add4(simd::mul4(m.row1, t), minor2);
	t = simd::swizzle<2, 3, 0, 1>(t);
	minor1 = simd::add4(simd::mul4(r2, t), minor1);
	minor2 = simd::sub4(minor2, simd::mul4(m.row1, t));
	t = simd::mul4(m.row0, r2);
	t = simd::swizzle<1, 0, 3, 2>(t);
	minor1 = simd::add4(simd::mul4(m.row3, t), minor1);
	minor3 = simd::sub4(minor3, simd::mul4(m.row1, t));
	t = simd::swizzle<2, 3, 0, 1>(t);
	minor1 = simd::sub4(minor1, simd::mul4(m.row3, t));
	minor3 = simd::add4(simd::mul4(m.row1, t), minor3);
	auto det = simd::mul4(m.row0, minor0);
	det = simd::add4(simd::swizzle<2, 3, 0, 1>(det), det);
	det = simd::add4(simd::swizzle<1, 0, 3, 2>(det), det);
	//static const auto one = simd::set4(1.f);
	auto n = simd::div4(/*one*/Vector4<float>::ONE, det);
	row0 = simd::mul4(minor0, n);
	row1 = simd::mul4(minor1, n);
	row2 = simd::mul4(minor2, n);
	row3 = simd::mul4(minor3, n);
	return *this;
}

inline Matrix4<float>& Matrix4<float>::preConcatenate(const Matrix4<float>& m)
{
	auto r0 = simd::mul4(simd::broadcast<0>(m.row0), row0);
	r0 = simd::add4(r0, simd::mul4(simd::broadcast<1>(m.row0), row1));
	r0 = simd::add4(r0, simd::mul4(simd::broadcast<2>(m.row0), row2));
	r0 = simd::add4(r0, simd::mul4(simd::broadcast<3>(m.row0), row3));
	auto r1 = simd::mul4(simd::broadcast<0>(m.row1), row0);
	r1 = simd::add4(r1, simd::mul4(simd::broadcast<1>(m.row1), row1));
	r1 = simd::add4(r1, simd::mul4(simd::broadcast<2>(m.row1), row2));
	r1 = simd::add4(r1, simd::mul4(simd::broadcast<3>(m.row1), row3));
	auto r2 = simd::mul4(simd::broadcast<0>(m.row2), row0);
	r2 = simd::add4(r2, simd::mul4(simd::broadcast<1>(m.row2), row1));
	r2 = simd::add4(r2, simd::mul4(simd::broadcast<2>(m.row2), row2));
	r2 = simd::add4(r2, simd::mul4(simd::broadcast<3>(m.row2), row3));
	auto r3 = simd::mul4(simd::broadcast<0>(m.row3), row0);
	r3 = simd::add4(r3, simd::mul4(simd::broadcast<1>(m.row3), row1));
	r3 = simd::add4(r3, simd::mul4(simd::broadcast<2>(m.row3), row2));
	r3 = simd::add4(r3, simd::mul4(simd::broadcast<3>(m.row3), row3));
	row0 = r0;
	row1 = r1;
	row2 = r2;
	row3 = r3;
	return *this;
}

inline Matrix4<float>& Matrix4<float>::preScale(const Vector3<float>& v)
{
	row0 = simd::mul4(row0, simd::xxxx(v));
	row1 = simd::mul4(row1, simd::yyyy(v));
	row2 = simd::mul4(row2, simd::zzzz(v));
	return *this;
}

inline Matrix4<float>& Matrix4<float>::scale(const Vector3<float>& v)
{
	auto t = simd::insert3(v, Vector4<float>::UNIT_W);
	row0 = simd::mul4(row0, t);
	row1 = simd::mul4(row1, t);
	row2 = simd::mul4(row2, t);
	row3 = simd::mul4(row3, t);
	return *this;
}

inline Matrix4<float>& Matrix4<float>::negate()
{
	row0 = simd::neg4(row0);
	row1 = simd::neg4(row1);
	row2 = simd::neg4(row2);
	row3 = simd::neg4(row3);
	return *this;
}

inline Matrix4<float>& Matrix4<float>::transpose()
{
	//simd::transpose4x4(row0, row1, row2, row3, row0, row1, row2, row3);
	std::tie(row0, row1, row2, row3) = simd::transpose4x4(row0, row1, row2, row3);
	return *this;
}

#endif /* SIMD_HAS_FLOAT4 */

template<typename T>
	requires std::floating_point<T>
inline Matrix4<T> concatenate(const Matrix4<T>& m1, const Matrix4<T>& m2) noexcept
{
	return m1*m2;
}

template<typename T>
	requires std::floating_point<T>
inline Matrix4<T> concatenate(const Matrix4<T>& m1, const Matrix4<T>& m2, const Matrix4<T>& m3) noexcept
{
	return concatenate(concatenate(m1, m2), m3);
}

template<typename T>
	requires std::floating_point<T>
inline Matrix4<T> concatenate(const Matrix4<T>& m1, const Matrix4<T>& m2, const Matrix4<T>& m3, const Matrix4<T>& m4) noexcept
{
	return concatenate(concatenate(concatenate(m1, m2), m3), m4);
}

template<typename T>
	requires std::floating_point<T>
inline Matrix4<T> transpose(const Matrix4<T>& m) noexcept
{
	return Matrix4<T>(m.m00, m.m10, m.m20, m.m30,
		m.m01, m.m11, m.m21, m.m31,
		m.m02, m.m12, m.m22, m.m32,
		m.m03, m.m13, m.m23, m.m33);
}

template<typename T>
	requires std::floating_point<T>
inline Matrix4<T> transpose(Matrix4<T>&& m) noexcept
{
	m.transpose();
	return m;
}

template<typename T>
	requires std::floating_point<T>
inline Matrix4<T> inverse(const Matrix4<T>& m) noexcept
{
	return Matrix4<T>(Uninitialized()).setInverse(m);
}

template<typename T>
	requires std::floating_point<T>
inline Matrix4<T> inverse(Matrix4<T>&& m) noexcept
{
	m.invert();
	return m;
}

template<typename T>
	requires std::floating_point<T>
inline Matrix4<T> adjoint(const Matrix4<T>& m) noexcept
{
	auto det3 = [](T m00, T m01, T m02, T m10, T m11, T m12, T m20, T m21, T m22) -> T 
		{ 
			return m00*(m11*m22 - m12*m21) - m01*(m10*m22 - m12*m20) + m02*(m10*m21 - m11*m20);
		};

	return Matrix4<T>(det3(m.m11, m.m21, m.m31, m.m12, m.m22, m.m32, m.m13, m.m23, m.m33),
		-det3(m.m01, m.m21, m.m31, m.m02, m.m22, m.m32, m.m03, m.m23, m.m33),
		det3(m.m01, m.m11, m.m31, m.m02, m.m12, m.m32, m.m03, m.m13, m.m33),
		-det3(m.m01, m.m11, m.m21, m.m02, m.m12, m.m22, m.m03, m.m13, m.m23),
		-det3(m.m10, m.m20, m.m30, m.m12, m.m22, m.m32, m.m13, m.m23, m.m33),
		det3(m.m00, m.m20, m.m30, m.m02, m.m22, m.m32, m.m03, m.m23, m.m33),
		-det3(m.m00, m.m10, m.m30, m.m02, m.m12, m.m32, m.m03, m.m13, m.m33),
		det3(m.m00, m.m10, m.m20, m.m02, m.m12, m.m22, m.m03, m.m13, m.m23),
		det3(m.m10, m.m20, m.m30, m.m11, m.m21, m.m31, m.m13, m.m23, m.m33),
		-det3(m.m00, m.m20, m.m30, m.m01, m.m21, m.m31, m.m03, m.m23, m.m33),
		det3(m.m00, m.m10, m.m30, m.m01, m.m11, m.m31, m.m03, m.m13, m.m33),
		-det3(m.m00, m.m10, m.m20, m.m01, m.m11, m.m21, m.m03, m.m13, m.m23),
		-det3(m.m10, m.m20, m.m30, m.m11, m.m21, m.m31, m.m12, m.m22, m.m32),
		det3(m.m00, m.m20, m.m30, m.m01, m.m21, m.m31, m.m02, m.m22, m.m32),
		-det3(m.m00, m.m10, m.m30, m.m01, m.m11, m.m31, m.m02, m.m12, m.m32),
		det3(m.m00, m.m10, m.m20, m.m01, m.m11, m.m21, m.m02, m.m12, m.m22));
}

#if SIMD_HAS_FLOAT4

template<>
inline Matrix4<float> transpose(const Matrix4<float>& m) noexcept
{
	Matrix4<float> n(Uninitialized());
	//simd::transpose4x4(m.row0, m.row1, m.row2, m.row3, n.row0, n.row1, n.row2, n.row3);
	std::tie(n.row0, n.row1, n.row2, n.row3) = simd::transpose4x4(m.row0, m.row1, m.row2, m.row3);
	return n;
}

template<>
inline Matrix4<float> transpose(Matrix4<float>&& m) noexcept
{
	//simd::transpose4x4(m.row0, m.row1, m.row2, m.row3, m.row0, m.row1, m.row2, m.row3);
	std::tie(m.row0, m.row1, m.row2, m.row3) = simd::transpose4x4(m.row0, m.row1, m.row2, m.row3);
	return m;
}

#endif /* SIMD_HAS_FLOAT4 */

} // namespace templates

using Matrix4f = templates::Matrix4<float>;
using Matrix4d = templates::Matrix4<double>;

#if MATHEMATICS_DOUBLE
using Matrix4 = templates::Matrix4<double>;
using Matrix4Arg = templates::Matrix4<double>::ConstArg;
using Matrix4Result = templates::Matrix4<double>::ConstResult;
#else
using Matrix4 = templates::Matrix4<float>;
using Matrix4Arg = templates::Matrix4<float>::ConstArg;
using Matrix4Result = templates::Matrix4<float>::ConstResult;
#endif

} // namespace core::mathematics

#include "Quaternion.hpp"
#include "YawPitchRoll.hpp"
#include "Euler.hpp"
#include "Plane.hpp"

namespace core::mathematics::templates {

template<typename T>
inline Matrix4<T>& Matrix4<T>::setRotation(const YawPitchRoll<T>& r)
{
	if (!r.isZero())
	{
		T sy = std::sin(r.yaw), cy = std::cos(r.yaw);
		T sp = std::sin(r.pitch), cp = std::cos(r.pitch);
		T sr = std::sin(r.roll), cr = std::cos(r.roll);
		T ss = sy*sp, cs = cy*sp;

		return set(cr*cy + sr*ss, sr*cp, sr*cs - cr*sy, T(),
			cr*ss - sr*cy, cr*cp, cr*cs + sr*sy, T(),
			sy*cp, -sp, cy*cp, T(),
			T(), T(), T(), T(1));
	}

	return setIdentity();
}

template<typename T>
inline Matrix4<T>& Matrix4<T>::setRotation(const Euler<T>& e)
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

		return set(m[0][0], m[0][1], m[0][2], T(),
			m[1][0], m[1][1], m[1][2], T(),
			m[2][0], m[2][1], m[2][2], T(), 
			T(), T(), T(), T(1));
	}

	return setIdentity();
}

template<typename T>
inline Matrix4<T>& Matrix4<T>::setRotation(const Quaternion<T>& q)
{
	if (!q.isIdentity())
	{
		T nq = q.getNorm();
		T s = (nq > T(0)) ? T(2)/nq : T(2);
		T xs = q.x*s, ys = q.y*s, zs = q.z*s;
		T wx = q.w*xs, wy = q.w*ys, wz = q.w*zs;
		T xx = q.x*xs, xy = q.x*ys, xz = q.x*zs;
		T yy = q.y*ys, yz = q.y*zs, zz = q.z*zs;

		return set(T(1) - (yy + zz), xy + wz, xz - wy, T(),
			xy - wz, T(1) - (xx + zz), yz + wx, T(),
			xz + wy, yz - wx, T(1) - (xx + yy), T(),
			T(), T(), T(), T(1));
	}

	return setIdentity();
}

template<typename T>
inline Matrix4<T>& Matrix4<T>::setRotationTranslation(const YawPitchRoll<T>& r, const Vector3<T>& t)
{
	setRotation(r);
	m30 = t.x; m31 = t.y; m32 = t.z;
	return *this;
}

template<typename T>
inline Matrix4<T>& Matrix4<T>::setRotationTranslation(const Euler<T>& e, const Vector3<T>& t)
{
	setRotation(e);
	m30 = t.x; m31 = t.y; m32 = t.z;
	return *this;
}

template<typename T>
inline Matrix4<T>& Matrix4<T>::setRotationTranslation(const Quaternion<T>& q, const Vector3<T>& t)
{
	setRotation(q);
	m30 = t.x; m31 = t.y; m32 = t.z;
	return *this;
}

template<typename T>
inline Matrix4<T>& Matrix4<T>::setScalingRotation(const Vector3<T>& s, const YawPitchRoll<T>& r)
{
	setRotation(r);
	m00 *= s.x; m01 *= s.x; m02 *= s.x;
	m10 *= s.y; m11 *= s.y; m12 *= s.y;
	m20 *= s.z; m21 *= s.z; m22 *= s.z;
	return *this;
}

template<typename T>
inline Matrix4<T>& Matrix4<T>::setScalingRotation(const Vector3<T>& s, const Euler<T>& e)
{
	setRotation(e);
	m00 *= s.x; m01 *= s.x; m02 *= s.x;
	m10 *= s.y; m11 *= s.y; m12 *= s.y;
	m20 *= s.z; m21 *= s.z; m22 *= s.z;
	return *this;
}

template<typename T>
inline Matrix4<T>& Matrix4<T>::setScalingRotation(const Vector3<T>& s, const Quaternion<T>& q)
{
	setRotation(q);
	m00 *= s.x; m01 *= s.x; m02 *= s.x;
	m10 *= s.y; m11 *= s.y; m12 *= s.y;
	m20 *= s.z; m21 *= s.z; m22 *= s.z;
	return *this;
}

template<typename T>
inline Matrix4<T>& Matrix4<T>::setScalingRotationTranslation(const Vector3<T>& s, const YawPitchRoll<T>& r, const Vector3<T>& t)
{
	setRotation(r);
	m00 *= s.x; m01 *= s.x; m02 *= s.x;
	m10 *= s.y; m11 *= s.y; m12 *= s.y;
	m20 *= s.z; m21 *= s.z; m22 *= s.z;
	m30 = t.x; m31 = t.y; m32 = t.z;
	return *this;
}

template<typename T>
inline Matrix4<T>& Matrix4<T>::setScalingRotationTranslation(const Vector3<T>& s, const Euler<T>& e, const Vector3<T>& t)
{
	setRotation(e);
	m00 *= s.x; m01 *= s.x; m02 *= s.x;
	m10 *= s.y; m11 *= s.y; m12 *= s.y;
	m20 *= s.z; m21 *= s.z; m22 *= s.z;
	m30 = t.x; m31 = t.y; m32 = t.z;
	return *this;
}

template<typename T>
inline Matrix4<T>& Matrix4<T>::setScalingRotationTranslation(const Vector3<T>& s, const Quaternion<T>& q, const Vector3<T>& t)
{
	setRotation(q);
	m00 *= s.x; m01 *= s.x; m02 *= s.x;
	m10 *= s.y; m11 *= s.y; m12 *= s.y;
	m20 *= s.z; m21 *= s.z; m22 *= s.z;
	m30 = t.x; m31 = t.y; m32 = t.z;
	return *this;
}

template<typename T>
inline Matrix4<T>& Matrix4<T>::setReflection(const Plane<T>& plane)
{
	Plane<T> p(plane);
	p.normalize();
	m00 = T(-2)*p.a*p.a + T(1); m01 = T(-2)*p.b*p.a; m02 = T(-2)*p.c*p.a; m03 = T();
	m10 = T(-2)*p.a*p.b; m11 = T(-2)*p.b*p.b + T(1); m12 = T(-2)*p.c*p.b; m13 = T();
	m20 = T(-2)*p.a*p.c; m21 = T(-2)*p.b*p.c; m22 = T(-2)*p.c*p.c + T(1); m23 = T();
	m30 = T(-2)*p.a*p.d; m31 = T(-2)*p.b*p.d; m32 = T(-2)*p.c*p.d; m33 = T(1);
	return *this;
}

#if SIMD_HAS_FLOAT4

inline Matrix4<float>& Matrix4<float>::setRotation(const YawPitchRoll<float>& r)
{
	if (!r.isZero())
	{
		float sy = std::sin(r.yaw), cy = std::cos(r.yaw);
		float sp = std::sin(r.pitch), cp = std::cos(r.pitch);
		float sr = std::sin(r.roll), cr = std::cos(r.roll);
		float ss = sy*sp, cs = cy*sp;

		return set(simd::set3(cr*cy + sr*ss, sr*cp, sr*cs - cr*sy),
			simd::set3(cr*ss - sr*cy, cr*cp, cr*cs + sr*sy),
			simd::set3(sy*cp, -sp, cy*cp),
			(simd::float4)Vector4<float>::UNIT_W);
	}

	return setIdentity();
}

inline Matrix4<float>& Matrix4<float>::setRotation(const Euler<float>& e)
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

		return set(simd::set3(m[0][0], m[0][1], m[0][2]),
			simd::set3(m[1][0], m[1][1], m[1][2]),
			simd::set3(m[2][0], m[2][1], m[2][2]),
			(simd::float4)Vector4<float>::UNIT_W);
	}

	return setIdentity();
}

inline Matrix4<float>& Matrix4<float>::setRotation(const Quaternion<float>& q)
{
	if (!q.isIdentity())
	{
		float nq = q.getNorm();
		float s = (nq > 0.f) ? 2.f/nq : 2.f;
		float xs = q.x*s, ys = q.y*s, zs = q.z*s; // #TODO SIMD
		float wx = q.w*xs, wy = q.w*ys, wz = q.w*zs;
		float xx = q.x*xs, xy = q.x*ys, xz = q.x*zs;
		float yy = q.y*ys, yz = q.y*zs, zz = q.z*zs;

		return set(simd::set3(1.f - (yy + zz), xy + wz, xz - wy),
			simd::set3(xy - wz, 1.f - (xx + zz), yz + wx),
			simd::set3(xz + wy, yz - wx, 1.f - (xx + yy)),
			(simd::float4)Vector4<float>::UNIT_W);
	}

	return setIdentity();
}

inline Matrix4<float>& Matrix4<float>::setRotationTranslation(const YawPitchRoll<float>& r, const Vector3<float>& t)
{
	setRotation(r);
	row3 = simd::insert3(t, Vector4<float>::UNIT_W);
	return *this;
}

inline Matrix4<float>& Matrix4<float>::setRotationTranslation(const Euler<float>& e, const Vector3<float>& t)
{
	setRotation(e);
	row3 = simd::insert3(t, Vector4<float>::UNIT_W);
	return *this;
}

inline Matrix4<float>& Matrix4<float>::setRotationTranslation(const Quaternion<float>& q, const Vector3<float>& t)
{
	setRotation(q);
	row3 = simd::insert3(t, Vector4<float>::UNIT_W);
	return *this;
}

inline Matrix4<float>& Matrix4<float>::setScalingRotationTranslation(const Vector3<float>& s, const YawPitchRoll<float>& r, 
	const Vector3<float>& t)
{
	setRotation(r);
	preScale(s);
	row3 = simd::insert3(t, Vector4<float>::UNIT_W);
	return *this;
}

inline Matrix4<float>& Matrix4<float>::setScalingRotationTranslation(const Vector3<float>& s, const Euler<float>& e, 
	const Vector3<float>& t)
{
	setRotation(e);
	preScale(s);
	row3 = simd::insert3(t, Vector4<float>::UNIT_W);
	return *this;
}

inline Matrix4<float>& Matrix4<float>::setScalingRotationTranslation(const Vector3<float>& s, const Quaternion<float>& q, 
	const Vector3<float>& t)
{
	setRotation(q);
	preScale(s);
	row3 = simd::insert3(t, Vector4<float>::UNIT_W);
	return *this;
}

inline Matrix4<float>& Matrix4<float>::setReflection(const Plane<float>& plane)
{
	Plane<float> p(plane);
	p.normalize();
	static const auto mTwo = simd::set4(-2.f);
	auto n = simd::cutoff3(p);
	auto t = simd::mul4(mTwo, p);
	row0 = simd::mulAdd4(simd::xxxx(t), n, Vector4<float>::UNIT_X);
	row1 = simd::mulAdd4(simd::yyyy(t), n, Vector4<float>::UNIT_Y);
	row2 = simd::mulAdd4(simd::zzzz(t), n, Vector4<float>::UNIT_Z);
	row3 = simd::mulAdd4(simd::wwww(t), n, Vector4<float>::UNIT_W);
	return *this;
}

#endif /* SIMD_HAS_FLOAT4 */

} // namespace core::mathematics::templates
