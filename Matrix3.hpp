/*
 *	Name: Matrix3
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
#include "Matrix2.hpp"

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
struct Matrix3
{
	using Real = T;
	using ComponentType = T;
	using ConstArg = const Matrix3&;
	using ConstResult = const Matrix3&;

	//static constexpr int NUM_COMPONENTS = 9;

	constexpr Matrix3() noexcept : m00(), m01(), m02(), m10(), m11(), m12(), m20(), m21(), m22() {}
	explicit Matrix3(Uninitialized) noexcept {}
	explicit Matrix3(Identity) noexcept : m00(1), m01(), m02(), m10(), m11(1), m12(), m20(), m21(), m22(1) {}
	constexpr Matrix3(T m00, T m01, T m02, T m10, T m11, T m12, T m20, T m21, T m22) noexcept;
	constexpr Matrix3(const Vector3<T>& row0, const Vector3<T>& row1, const Vector3<T>& row2) noexcept;
	constexpr explicit Matrix3(const std::tuple<Vector3<T>, Vector3<T>, Vector3<T>>& t) noexcept;
	constexpr Matrix3(const Matrix2<T>& m) noexcept;
	explicit Matrix3(const YawPitchRoll<T>& r) noexcept { setRotation(r); }
	explicit Matrix3(const Euler<T>& e) noexcept { setRotation(e); }
	explicit Matrix3(const Quaternion<T>& q) noexcept; // the input quaternion has to be normalized
	explicit Matrix3(const Vector3<T>& forward) noexcept;
	Matrix3(const Vector3<T>& up, const Vector3<T>& forward) noexcept;
	explicit Matrix3(const T* m) noexcept;

	//explicit operator T*() noexcept { return &m00; }
	//explicit operator const T*() const noexcept { return &m00; }
	Vector3<T>& operator[](int i) noexcept { return reinterpret_cast<Vector3<T>*>(&m00)[i]; }
	const Vector3<T>& operator[](int i) const noexcept { return reinterpret_cast<const Vector3<T>*>(&m00)[i]; }

	Matrix3 operator+() const noexcept { return *this; }
	Matrix3 operator-() const noexcept { return Matrix3(-m00, -m01, -m02, -m10, -m11, -m12, -m20, -m21, -m22); }
	Matrix3& operator+=(const Matrix3& m) noexcept;
	Matrix3& operator-=(const Matrix3& m) noexcept;
	Matrix3& operator*=(T f) noexcept;
	Matrix3& operator/=(T f) noexcept { return operator*=(T(1)/f); }
	Matrix3& operator*=(const Matrix3& m) noexcept;
	bool operator==(const Matrix3& m) const noexcept;
	bool operator!=(const Matrix3& m) const noexcept { return !(*this == m); }

	template<typename A> void load(A& ar) { ar(m00, m01, m02, m10, m11, m12, m20, m21, m22); }
	template<typename A> void save(A& ar) const { ar(m00, m01, m02, m10, m11, m12, m20, m21, m22); }

	// #TODO rename get...() to make...()
	static Matrix3 getScaling(const Vector3<T>& v) noexcept { return Matrix3(Uninitialized()).setScaling(v); }
	static Matrix3 getScaling(T f) noexcept { return Matrix3(Uninitialized()).setScaling(f); }
	static Matrix3 getRotation(Axis axis, T angle) noexcept { return Matrix3(Uninitialized()).setRotation(axis, angle); }
	static Matrix3 getRotation(const Vector3<T>& axis, T angle) noexcept { return Matrix3(Uninitialized()).setRotation(axis, angle); }
	static Matrix3 getRotation(const YawPitchRoll<T>& r) noexcept { return Matrix3(Uninitialized()).setRotation(r); }
	static Matrix3 getRotation(const Euler<T>& e) noexcept { return Matrix3(Uninitialized()).setRotation(e); }
	static Matrix3 getRotation(const Quaternion<T>& q) noexcept { return Matrix3(Uninitialized()).setRotation(q); }
	static Matrix3 getScalingRotation(const Vector3<T>& s, Axis axis, T angle) noexcept { return Matrix3(Uninitialized()).setScalingRotation(s, axis, angle); }
	static Matrix3 getScalingRotation(const Vector3<T>& s, const Vector3<T>& axis, T angle) noexcept { return Matrix3(Uninitialized()).setScalingRotation(s, axis, angle); }
	static Matrix3 getScalingRotation(const Vector3<T>& s, const YawPitchRoll<T>& r) noexcept { return Matrix3(Uninitialized()).setScalingRotation(s, r); }
	static Matrix3 getScalingRotation(const Vector3<T>& s, const Euler<T>& e) noexcept { return Matrix3(Uninitialized()).setScalingRotation(s, e); }
	static Matrix3 getScalingRotation(const Vector3<T>& s, const Quaternion<T>& q) noexcept { return Matrix3(Uninitialized()).setScalingRotation(s, q); }
	static Matrix3 getShearing(T xy, T xz, T yx, T yz, T zx, T zy) noexcept { return Matrix3(Uninitialized()).setShearing(xy, xz, yx, yz, zx, zy); }

	bool isZero() const noexcept;
	bool isApproxZero() const noexcept;
	bool isIdentity() const noexcept;
	bool isApproxIdentity() const noexcept;
	bool isApproxEqual(const Matrix3& m) const noexcept;
	bool isApproxEqual(const Matrix3& m, T tolerance) const noexcept;
	bool isApproxOrthogonal() const noexcept;
	bool hasApproxUniformScaling() const noexcept;
	bool isFinite() const noexcept;
	const Vector3<T>& getRow(int i) const noexcept { return reinterpret_cast<const Vector3<T>*>(&m00)[i]; }
	Vector3<T> getColumn(int i) const noexcept { return Vector3<T>((&m00)[i], (&m10)[i], (&m20)[i]); }
	Vector3<T> getScaleComponent() const noexcept { return Vector3<T>((*this)[0].getMagnitude(), (*this)[1].getMagnitude(), (*this)[2].getMagnitude()); }
	T getTrace() const noexcept { return (m00 + m11 + m22); }
	T getDeterminant() const noexcept;
	bool isSingular() const noexcept { return (getDeterminant() == T(0)); }
	Matrix3& setZero/*zero*/() noexcept;
	Matrix3& setIdentity/*makeIdentity*/() noexcept;
	Matrix3& set(const Vector3<T>& row0, const Vector3<T>& row1, const Vector3<T>& row2) noexcept;
	Matrix3& set(T m00, T m01, T m02, T m10, T m11, T m12, T m20, T m21, T m22) noexcept;
	Matrix3& setScaling/*makeScaling*/(const Vector3<T>& v) noexcept;
	Matrix3& setScaling/*makeScaling*/(T f) noexcept;
	Matrix3& setRotation(Axis axis, T angle) noexcept;
	Matrix3& setRotation(const Vector3<T>& axis, T angle) noexcept;
	Matrix3& setRotation(const YawPitchRoll<T>& r) noexcept;
	Matrix3& setRotation(const Euler<T>& e) noexcept;
	Matrix3& setRotation(const Quaternion<T>& q) noexcept;
	Matrix3& setScalingRotation(const Vector3<T>& s, Axis axis, T angle) noexcept { return setRotation(axis, angle).preScale(s); }
	Matrix3& setScalingRotation(const Vector3<T>& s, const Vector3<T>& axis, T angle) noexcept { return setRotation(axis, angle).preScale(s); }
	Matrix3& setScalingRotation(const Vector3<T>& s, const YawPitchRoll<T>& r) noexcept { return setRotation(r).preScale(s); }
	Matrix3& setScalingRotation(const Vector3<T>& s, const Euler<T>& e) noexcept { return setRotation(e).preScale(s); }
	Matrix3& setScalingRotation(const Vector3<T>& s, const Quaternion<T>& q) noexcept { return setRotation(q).preScale(s); }
	Matrix3& setShearing(T xy, T xz, T yx, T yz, T zx, T zy) noexcept;
	Matrix3& setTranspose/*transposeOf*/(const Matrix3& m) noexcept;
	Matrix3& setInverse/*inverseOf*/(const Matrix3& m) noexcept;
	Matrix3& setInverseTranspose/*inverseTransposeOf*/(const Matrix3& m) noexcept;
	Matrix3& preConcatenate(const Matrix3& m) noexcept;
	Matrix3& concatenate(const Matrix3& m) noexcept { return operator*=(m); }
	Matrix3& preScale(const Vector3<T>& v) noexcept;
	Matrix3& preScale(T f) noexcept { return operator*=(f); }
	Matrix3& scale(const Vector3<T>& v) noexcept;
	Matrix3& scale(T f) noexcept { return operator*=(f); }
	//Matrix3& rotate(Axis axis, T angle) noexcept { concatenate(Matrix3(Uninitialized()).setRotation(axis, angle)); return *this; }
	//Matrix3& rotate(const Vector3<T>& axis, T angle) noexcept { concatenate(Matrix3(Uninitialized()).setRotation(axis, angle)); return *this; }
	//Matrix3& rotate(const YawPitchRoll<T>& r) noexcept { concatenate(Matrix3(Uninitialized()).setRotation(r)); return *this; }
	//Matrix3& rotate(const Euler<T>& e) noexcept { concatenate(Matrix3(Uninitialized()).setRotation(e)); return *this; }
	//Matrix3& rotate(const Quaternion<T>& q) noexcept { concatenate(Matrix3(Uninitialized()).setRotation(q)); return *this; }
	//Matrix3& shear(T xy, T xz, T yx, T yz, T zx, T zy) noexcept { concatenate(Matrix3(Uninitialized()).setShearing(xy, xz, yx, yz, zx, zy)); return *this; }
	Matrix3& negate() noexcept;
	Matrix3& transpose() noexcept;
	Matrix3& invert() noexcept;
	Matrix3& orthonormalize() noexcept;

	//static const Matrix3& getZero() noexcept { return ZERO; }
	//static const Matrix3& getIdentity() noexcept { return IDENTITY; }

	static const Matrix3 ZERO;
	static const Matrix3 IDENTITY;

	T m00, m01, m02;
	T m10, m11, m12;
	T m20, m21, m22;
};

template<typename T> const Matrix3<T> Matrix3<T>::ZERO{};
template<typename T> const Matrix3<T> Matrix3<T>::IDENTITY{ T(1), T(0), T(0), T(0), T(1), T(0), T(0), T(0), T(1) };

#if SIMD_HAS_FLOAT4

template<>
struct Quaternion<float>;

template<>
struct YawPitchRoll<float>;

template<>
struct Euler<float>;

template<>
struct Matrix3<float>
{
	using Real = float;
	using ComponentType = float;
	using ConstArg = const Matrix3;
	using ConstResult = const Matrix3;

	//static constexpr int NUM_COMPONENTS = 9;

	/*constexpr*/ Matrix3() noexcept;
	explicit Matrix3(Uninitialized) noexcept {}
	explicit Matrix3(Identity) noexcept;
	/*constexpr*/ Matrix3(float m00, float m01, float m02, float m10, float m11, float m12, float m20, float m21, float m22) noexcept;
	/*constexpr*/ Matrix3(const Vector3<float>& row0, const Vector3<float>& row1, const Vector3<float>& row2) noexcept;
	/*constexpr*/ explicit Matrix3(const std::tuple<Vector3<float>, Vector3<float>, Vector3<float>>& t) noexcept;
	/*constexpr*/ Matrix3(const Matrix2<float>& m) noexcept;
	explicit Matrix3(const YawPitchRoll<float>& r) noexcept { setRotation(r); }
	explicit Matrix3(const Euler<float>& e) noexcept { setRotation(e); }
	explicit Matrix3(const Quaternion<float>& q) noexcept; // the input quaternion has to be normalized
	explicit Matrix3(const Vector3<float>& forward) noexcept;
	Matrix3(const Vector3<float>& up, const Vector3<float>& forward) noexcept;
	explicit Matrix3(const float* m) noexcept;
	explicit Matrix3(const simd::float4* m) noexcept : row0(m[0]), row1(m[1]), row2(m[2]) {}
	Matrix3(simd::float4 row0, simd::float4 row1, simd::float4 row2) noexcept : row0(row0), row1(row1), row2(row2) {}
	explicit Matrix3(const std::tuple<simd::float4, simd::float4, simd::float4>& t) noexcept;
	Matrix3(const Matrix3& m) noexcept : row0(m.row0), row1(m.row1), row2(m.row2) {}
	Matrix3& operator=(const Matrix3& m) noexcept { row0 = m.row0; row1 = m.row1; row2 = m.row2; return *this; }

	//explicit operator float*() noexcept { return &m00; }
	//explicit operator const float*() const noexcept { return &m00; }
	explicit operator simd::float4*() noexcept { return &row0; }
	explicit operator const simd::float4*() const noexcept { return &row0; }
	Vector3<float>& operator[](int i) noexcept { return reinterpret_cast<Vector3<float>&>((&row0)[i]); }
	const Vector3<float>& operator[](int i) const noexcept { return reinterpret_cast<const Vector3<float>&>((&row0)[i]); }

	Matrix3 operator+() const noexcept { return *this; }
	Matrix3 operator-() const noexcept;
	Matrix3& operator+=(const Matrix3& m) noexcept;
	Matrix3& operator-=(const Matrix3& m) noexcept;
	Matrix3& operator*=(float f) noexcept;
	Matrix3& operator/=(float f) noexcept { return operator*=(1.f/f); }
	Matrix3& operator*=(const Matrix3& m) noexcept;
	bool operator==(const Matrix3& m) const noexcept;
	bool operator!=(const Matrix3& m) const noexcept { return !(*this == m); }

	template<typename A> void load(A& ar);
	template<typename A> void save(A& ar) const { ar(m00, m01, m02, m10, m11, m12, m20, m21, m22); }

	// #TODO rename get...() to make...()
	static Matrix3 getScaling(const Vector3<float>& v) noexcept { return Matrix3(Uninitialized()).setScaling(v); }
	static Matrix3 getScaling(float f) noexcept { return Matrix3(Uninitialized()).setScaling(f); }
	static Matrix3 getRotation(Axis axis, float angle) noexcept { return Matrix3(Uninitialized()).setRotation(axis, angle); }
	static Matrix3 getRotation(const Vector3<float>& axis, float angle) noexcept { return Matrix3(Uninitialized()).setRotation(axis, angle); }
	static Matrix3 getRotation(const YawPitchRoll<float>& r) noexcept { return Matrix3(Uninitialized()).setRotation(r); }
	static Matrix3 getRotation(const Euler<float>& e) noexcept { return Matrix3(Uninitialized()).setRotation(e); }
	static Matrix3 getRotation(const Quaternion<float>& q) noexcept { return Matrix3(Uninitialized()).setRotation(q); }
	static Matrix3 getScalingRotation(const Vector3<float>& s, Axis axis, float angle) noexcept { return Matrix3(Uninitialized()).setScalingRotation(s, axis, angle); }
	static Matrix3 getScalingRotation(const Vector3<float>& s, const Vector3<float>& axis, float angle) noexcept { return Matrix3(Uninitialized()).setScalingRotation(s, axis, angle); }
	static Matrix3 getScalingRotation(const Vector3<float>& s, const YawPitchRoll<float>& r) noexcept { return Matrix3(Uninitialized()).setScalingRotation(s, r); }
	static Matrix3 getScalingRotation(const Vector3<float>& s, const Euler<float>& e) noexcept { return Matrix3(Uninitialized()).setScalingRotation(s, e); }
	static Matrix3 getScalingRotation(const Vector3<float>& s, const Quaternion<float>& q) noexcept { return Matrix3(Uninitialized()).setScalingRotation(s, q); }
	static Matrix3 getShearing(float xy, float xz, float yx, float yz, float zx, float zy) noexcept { return Matrix3(Uninitialized()).setShearing(xy, xz, yx, yz, zx, zy); }

	bool isZero() const noexcept;
	bool isApproxZero() const noexcept;
	bool isIdentity() const noexcept;
	bool isApproxIdentity() const noexcept;
	bool isApproxEqual(const Matrix3& m) const noexcept;
	bool isApproxEqual(const Matrix3& m, float tolerance) const noexcept;
	bool isApproxOrthogonal() const noexcept;
	bool hasApproxUniformScaling() const noexcept;
	bool isFinite() const noexcept;
	const Vector3<float>& getRow(int i) const noexcept { return reinterpret_cast<const Vector3<float>&>((&row0)[i]); }
	Vector3<float> getColumn(int i) const noexcept; // #TODO
	Vector3<float> getScaleComponent() const noexcept { return Vector3<float>((*this)[0].getMagnitude(), (*this)[1].getMagnitude(), (*this)[2].getMagnitude()); }
	float getTrace() const noexcept { return (m00 + m11 + m22); }
	float getDeterminant() const noexcept;
	bool isSingular() const noexcept { return (getDeterminant() == 0.f); }
	Matrix3& setZero/*zero*/() noexcept;
	Matrix3& setIdentity/*makeIdentity*/() noexcept;
	Matrix3& set(const Vector3<float>& row0, const Vector3<float>& row1, const Vector3<float>& row2) noexcept;
	Matrix3& set(simd::float4 row0, simd::float4 row1, simd::float4 row2) noexcept;
	Matrix3& set(float m00, float m01, float m02, float m10, float m11, float m12, float m20, float m21, float m22) noexcept;
	Matrix3& setScaling/*makeScaling*/(const Vector3<float>& v) noexcept;
	Matrix3& setScaling/*makeScaling*/(float f) noexcept;
	Matrix3& setRotation(Axis axis, float angle) noexcept;
	Matrix3& setRotation(const Vector3<float>& axis, float angle) noexcept;
	Matrix3& setRotation(const YawPitchRoll<float>& r) noexcept;
	Matrix3& setRotation(const Euler<float>& e) noexcept;
	Matrix3& setRotation(const Quaternion<float>& q) noexcept;
	Matrix3& setScalingRotation(const Vector3<float>& s, Axis axis, float angle) noexcept { return setRotation(axis, angle).preScale(s); }
	Matrix3& setScalingRotation(const Vector3<float>& s, const Vector3<float>& axis, float angle) noexcept { return setRotation(axis, angle).preScale(s); }
	Matrix3& setScalingRotation(const Vector3<float>& s, const YawPitchRoll<float>& r) noexcept { return setRotation(r).preScale(s); }
	Matrix3& setScalingRotation(const Vector3<float>& s, const Euler<float>& e) noexcept { return setRotation(e).preScale(s); }
	Matrix3& setScalingRotation(const Vector3<float>& s, const Quaternion<float>& q) noexcept { return setRotation(q).preScale(s); }
	Matrix3& setShearing(float xy, float xz, float yx, float yz, float zx, float zy) noexcept;
	Matrix3& setTranspose/*transposeOf*/(const Matrix3& m) noexcept;
	Matrix3& setInverse/*inverseOf*/(const Matrix3& m) noexcept;
	Matrix3& setInverseTranspose/*inverseTransposeOf*/(const Matrix3& m) noexcept;
	Matrix3& preConcatenate(const Matrix3& m) noexcept;
	Matrix3& concatenate(const Matrix3& m) noexcept { return operator*=(m); }
	Matrix3& preScale(const Vector3<float>& v) noexcept;
	Matrix3& preScale(float f) noexcept { return operator*=(f); }
	Matrix3& scale(const Vector3<float>& v) noexcept;
	Matrix3& scale(float f) noexcept { return operator*=(f); }
	//Matrix3& rotate(Axis axis, float angle) noexcept { concatenate(Matrix3(Uninitialized()).setRotation(axis, angle)); return *this; }
	//Matrix3& rotate(const Vector3<float>& axis, float angle) noexcept { concatenate(Matrix3(Uninitialized()).setRotation(axis, angle)); return *this; }
	//Matrix3& rotate(const YawPitchRoll<float>& r) noexcept { concatenate(Matrix3(Uninitialized()).setRotation(r)); return *this; }
	//Matrix3& rotate(const Euler<float>& e) noexcept { concatenate(Matrix3(Uninitialized()).setRotation(e)); return *this; }
	//Matrix3& rotate(const Quaternion<float>& q) noexcept { concatenate(Matrix3(Uninitialized()).setRotation(q)); return *this; }
	//Matrix3& shear(float xy, float xz, float yx, float yz, float zx, float zy) noexcept { concatenate(Matrix3(Uninitialized()).setShearing(xy, xz, yx, yz, zx, zy)); return *this; }
	Matrix3& negate() noexcept;
	Matrix3& transpose() noexcept;
	Matrix3& invert() noexcept;
	Matrix3& orthonormalize() noexcept;

	//static const Matrix3& getZero() noexcept { return ZERO; }
	//static const Matrix3& getIdentity() noexcept { return IDENTITY; }

	static const Matrix3 ZERO;
	static const Matrix3 IDENTITY;

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
};

const Matrix3<float> Matrix3<float>::ZERO{};
const Matrix3<float> Matrix3<float>::IDENTITY{ 1.f, 0.f, 0.f, 0.f, 1.f, 0.f, 0.f, 0.f, 1.f };

#endif /* SIMD_HAS_FLOAT4 */

template<typename T>
inline Matrix3<T>::Matrix3(T m00, T m01, T m02, T m10, T m11, T m12, T m20, T m21, T m22) :
	m00(m00), m01(m01), m02(m02), m10(m10), m11(m11), m12(m12), m20(m20), m21(m21), m22(m22)
{
}

template<typename T>
inline Matrix3<T>::Matrix3(const Vector3<T>& row0, const Vector3<T>& row1, const Vector3<T>& row2) :
	m00(row0.x), m01(row0.y), m02(row0.z), m10(row1.x), m11(row1.y), m12(row1.z), m20(row2.x), m21(row2.y), m22(row2.z)
{
}

template<typename T>
inline Matrix3<T>::Matrix3(const std::tuple<Vector3<T>, Vector3<T>, Vector3<T>>& t) :
	Matrix3(std::get<0>(t), std::get<1>(t), std::get<2>(t))
{
}

template<typename T>
inline Matrix3<T>::Matrix3(const Matrix2<T>& m) : 
	m00(m.m00), m01(m.m01), m02(), m10(m.m10), m11(m.m11), m12(), m20(), m21(), m22(1) 
{
}

template<typename T>
inline Matrix3<T>::Matrix3(const Vector3<T>& forward)
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
			set(right, up, at);
		}
		else
		{
			Vector3<T> up((at.y > T(0)) ? -Vector3<T>::UNIT_Z : Vector3<T>::UNIT_Z);
			Vector3<T> right = cross(up, at);
			right.normalize();
			up = cross(at, right);
			set(right, up, at);
		}
	}
	else
	{
		setIdentity();
	}
}

template<typename T>
inline Matrix3<T>::Matrix3(const Vector3<T>& up, const Vector3<T>& forward)
{
	T m = forward.getMagnitude();
	if ((m > T(0)) && (up.getMagnitude() > T(0)))
	{
		Vector3<T> at = forward/m;
		Vector3<T> right = cross(up, at);
		right.normalize();
		set(right, cross(at, right), at);
	}
	else
	{
		setIdentity();
	}
}

template<typename T>
inline Matrix3<T>::Matrix3(const T* m) : 
	m00(m[0]), m01(m[1]), m02(m[2]), m10(m[3]), m11(m[4]), m12(m[5]), m20(m[6]), m21(m[7]), m22(m[8]) 
{
}

template<typename T>
inline Matrix3<T>& Matrix3<T>::operator+=(const Matrix3<T>& m)
{
	m00 += m.m00; m01 += m.m01; m02 += m.m02;
	m10 += m.m10; m11 += m.m11; m12 += m.m12;
	m20 += m.m20; m21 += m.m21; m22 += m.m22;
	return *this;
}

template<typename T>
inline Matrix3<T>& Matrix3<T>::operator-=(const Matrix3<T>& m)
{
	m00 -= m.m00; m01 -= m.m01; m02 -= m.m02;
	m10 -= m.m10; m11 -= m.m11; m12 -= m.m12;
	m20 -= m.m20; m21 -= m.m21; m22 -= m.m22;
	return *this;
}

template<typename T>
inline Matrix3<T>& Matrix3<T>::operator*=(T f)
{
	m00 *= f; m01 *= f; m02 *= f;
	m10 *= f; m11 *= f; m12 *= f;
	m20 *= f; m21 *= f; m22 *= f;
	return *this;
}

template<typename T>
inline Matrix3<T>& Matrix3<T>::operator*=(const Matrix3<T>& m)
{
	set(m00*m.m00 + m01*m.m10 + m02*m.m20,
		m00*m.m01 + m01*m.m11 + m02*m.m21,
		m00*m.m02 + m01*m.m12 + m02*m.m22,
		m10*m.m00 + m11*m.m10 + m12*m.m20,
		m10*m.m01 + m11*m.m11 + m12*m.m21,
		m10*m.m02 + m11*m.m12 + m12*m.m22,
		m20*m.m00 + m21*m.m10 + m22*m.m20,
		m20*m.m01 + m21*m.m11 + m22*m.m21,
		m20*m.m02 + m21*m.m12 + m22*m.m22);
	return *this;
}

template<typename T>
	requires std::floating_point<T>
inline Matrix3<T> operator+(const Matrix3<T>& m1, const Matrix3<T>& m2) noexcept
{
	return Matrix3<T>(m1.m00 + m2.m00, m1.m01 + m2.m01, m1.m02 + m2.m02,
		m1.m10 + m2.m10, m1.m11 + m2.m11, m1.m12 + m2.m12,
		m1.m20 + m2.m20, m1.m21 + m2.m21, m1.m22 + m2.m22);
}

template<typename T>
	requires std::floating_point<T>
inline Matrix3<T> operator-(const Matrix3<T>& m1, const Matrix3<T>& m2) noexcept
{
	return Matrix3<T>(m1.m00 - m2.m00, m1.m01 - m2.m01, m1.m02 - m2.m02,
		m1.m10 - m2.m10, m1.m11 - m2.m11, m1.m12 - m2.m12,
		m1.m20 - m2.m20, m1.m21 - m2.m21, m1.m22 - m2.m22);
}

template<typename T>
	requires std::floating_point<T>
inline Matrix3<T> operator*(T f, const Matrix3<T>& m) noexcept
{
	return Matrix3<T>(f*m.m00, f*m.m01, f*m.m02,
		f*m.m10, f*m.m11, f*m.m12,
		f*m.m20, f*m.m21, f*m.m22);
}

template<typename T>
	requires std::floating_point<T>
inline Matrix3<T> operator*(const Matrix3<T>& m, T f) noexcept
{
	return Matrix3<T>(m.m00*f, m.m01*f, m.m02*f,
		m.m10*f, m.m11*f, m.m12*f,
		m.m20*f, m.m21*f, m.m22*f);
}

template<typename T>
	requires std::floating_point<T>
inline Matrix3<T> operator*(const Matrix3<T>& m1, const Matrix3<T>& m2) noexcept
{
	return Matrix3<T>(m1.m00*m2.m00 + m1.m01*m2.m10 + m1.m02*m2.m20,
		m1.m00*m2.m01 + m1.m01*m2.m11 + m1.m02*m2.m21,
		m1.m00*m2.m02 + m1.m01*m2.m12 + m1.m02*m2.m22,
		m1.m10*m2.m00 + m1.m11*m2.m10 + m1.m12*m2.m20,
		m1.m10*m2.m01 + m1.m11*m2.m11 + m1.m12*m2.m21,
		m1.m10*m2.m02 + m1.m11*m2.m12 + m1.m12*m2.m22,
		m1.m20*m2.m00 + m1.m21*m2.m10 + m1.m22*m2.m20,
		m1.m20*m2.m01 + m1.m21*m2.m11 + m1.m22*m2.m21,
		m1.m20*m2.m02 + m1.m21*m2.m12 + m1.m22*m2.m22);
}

template<typename T>
	requires std::floating_point<T>
inline Matrix3<T> operator/(const Matrix3<T>& m, T f) noexcept
{ 
	return operator*(m, T(1)/f); 
}

template<typename T>
inline bool Matrix3<T>::operator==(const Matrix3<T>& m) const
{
	return (m00 == m.m00) && (m01 == m.m01) && (m02 == m.m02) &&
		(m10 == m.m10) && (m11 == m.m11) && (m12 == m.m12) &&
		(m20 == m.m20) && (m21 == m.m21) && (m22 == m.m22);
}

template<typename C, typename T, typename U>
	requires std::floating_point<U>
inline std::basic_istream<C, T>& operator>>(std::basic_istream<C, T>& s, Matrix3<U>& m)
{
	return s >> m.m00 >> std::ws >> m.m01 >> std::ws >> m.m02 >> std::ws >>
		m.m10 >> std::ws >> m.m11 >> std::ws >> m.m12 >> std::ws >>
		m.m20 >> std::ws >> m.m21 >> std::ws >> m.m22;
}

template<typename C, typename T, typename U>
	requires std::floating_point<U>
inline std::basic_ostream<C, T>& operator<<(std::basic_ostream<C, T>& s, const Matrix3<U>& m)
{
	constexpr C WS(0x20);
	return s << m.m00 << WS << m.m01 << WS << m.m02 << WS <<
		m.m10 << WS << m.m11 << WS << m.m12 << WS <<
		m.m20 << WS << m.m21 << WS << m.m22;
}

template<typename T>
inline bool Matrix3<T>::isZero() const
{
	return (m00 == T()) && (m01 == T()) && (m02 == T()) &&
		(m10 == T()) && (m11 == T()) && (m12 == T()) &&
		(m20 == T()) && (m21 == T()) && (m22 == T());
}

template<typename T>
inline bool Matrix3<T>::isApproxZero() const
{
	return (std::fabs(m00) < Constants<T>::TOLERANCE) &&
		(std::fabs(m01) < Constants<T>::TOLERANCE) &&
		(std::fabs(m02) < Constants<T>::TOLERANCE) &&
		(std::fabs(m10) < Constants<T>::TOLERANCE) &&
		(std::fabs(m11) < Constants<T>::TOLERANCE) &&
		(std::fabs(m12) < Constants<T>::TOLERANCE) &&
		(std::fabs(m20) < Constants<T>::TOLERANCE) &&
		(std::fabs(m21) < Constants<T>::TOLERANCE) &&
		(std::fabs(m22) < Constants<T>::TOLERANCE);
}

template<typename T>
inline bool Matrix3<T>::isIdentity() const
{
	return (m00 == T(1)) && (m01 == T()) && (m02 == T()) &&
		(m10 == T()) && (m11 == T(1)) && (m12 == T()) &&
		(m20 == T()) && (m21 == T()) && (m22 == T(1));
}

template<typename T>
inline bool Matrix3<T>::isApproxIdentity() const
{
	return (std::fabs(T(1) - m00) < Constants<T>::TOLERANCE) &&
		(std::fabs(m01) < Constants<T>::TOLERANCE) &&
		(std::fabs(m02) < Constants<T>::TOLERANCE) &&
		(std::fabs(m10) < Constants<T>::TOLERANCE) &&
		(std::fabs(T(1) - m11) < Constants<T>::TOLERANCE) &&
		(std::fabs(m12) < Constants<T>::TOLERANCE) &&
		(std::fabs(m20) < Constants<T>::TOLERANCE) &&
		(std::fabs(m21) < Constants<T>::TOLERANCE) &&
		(std::fabs(T(1) - m22) < Constants<T>::TOLERANCE);
}

template<typename T>
inline bool Matrix3<T>::isApproxEqual(const Matrix3<T>& m) const
{
	return (std::fabs(m.m00 - m00) < Constants<T>::TOLERANCE) &&
		(std::fabs(m.m01 - m01) < Constants<T>::TOLERANCE) &&
		(std::fabs(m.m02 - m02) < Constants<T>::TOLERANCE) &&
		(std::fabs(m.m10 - m10) < Constants<T>::TOLERANCE) &&
		(std::fabs(m.m11 - m11) < Constants<T>::TOLERANCE) &&
		(std::fabs(m.m12 - m12) < Constants<T>::TOLERANCE) &&
		(std::fabs(m.m20 - m20) < Constants<T>::TOLERANCE) &&
		(std::fabs(m.m21 - m21) < Constants<T>::TOLERANCE) &&
		(std::fabs(m.m22 - m22) < Constants<T>::TOLERANCE);
}

template<typename T>
inline bool Matrix3<T>::isApproxEqual(const Matrix3<T>& m, T tolerance) const
{
	return (std::fabs(m.m00 - m00) < tolerance) &&
		(std::fabs(m.m01 - m01) < tolerance) &&
		(std::fabs(m.m02 - m02) < tolerance) &&
		(std::fabs(m.m10 - m10) < tolerance) &&
		(std::fabs(m.m11 - m11) < tolerance) &&
		(std::fabs(m.m12 - m12) < tolerance) &&
		(std::fabs(m.m20 - m20) < tolerance) &&
		(std::fabs(m.m21 - m21) < tolerance) &&
		(std::fabs(m.m22 - m22) < tolerance);
}

template<typename T>
inline bool Matrix3<T>::isApproxOrthogonal() const
{
	Matrix3<T> m(*this);
	m *= Matrix3<T>(Uninitialized()).setTranspose(*this);
	return m.isApproxIdentity();
}

template<typename T>
inline bool Matrix3<T>::hasApproxUniformScaling() const
{
	Matrix3<T> m(*this);
	m *= Matrix3<T>(Uninitialized()).setTranspose(*this);
	return (std::fabs(m.m01) < Constants<T>::TOLERANCE) && (std::fabs(m.m02) < Constants<T>::TOLERANCE) &&
		(std::fabs(m.m10) < Constants<T>::TOLERANCE) && (std::fabs(m.m12) < Constants<T>::TOLERANCE) &&
		(std::fabs(m.m20) < Constants<T>::TOLERANCE) && (std::fabs(m.m21) < Constants<T>::TOLERANCE) &&
		(std::fabs(m.m11 - m.m00) < Constants<T>::TOLERANCE) && (std::fabs(m.m22 - m.m11) < Constants<T>::TOLERANCE) &&
		(std::fabs(m.m00 - m.m22) < Constants<T>::TOLERANCE);
}

template<typename T>
inline bool Matrix3<T>::isFinite() const
{
	return std::isfinite(m00) && std::isfinite(m01) && std::isfinite(m02) &&
		std::isfinite(m10) && std::isfinite(m11) && std::isfinite(m12) &&
		std::isfinite(m20) && std::isfinite(m21) && std::isfinite(m22);
}

template<typename T>
inline T Matrix3<T>::getDeterminant() const
{
	return m00*(m11*m22 - m12*m21) - m01*(m10*m22 - m12*m20) + m02*(m10*m21 - m11*m20);
}

template<typename T>
inline Matrix3<T>& Matrix3<T>::setZero()
{
	m00 = T(); m01 = T(); m02 = T();
	m10 = T(); m11 = T(); m12 = T();
	m20 = T(); m21 = T(); m22 = T();
	return *this;
}

template<typename T>
inline Matrix3<T>& Matrix3<T>::setIdentity()
{
	m00 = T(1); m01 = T(); m02 = T();
	m10 = T(); m11 = T(1); m12 = T();
	m20 = T(); m21 = T(); m22 = T(1);
	return *this;
}

template<typename T>
inline Matrix3<T>& Matrix3<T>::set(const Vector3<T>& row0, const Vector3<T>& row1, const Vector3<T>& row2)
{
	m00 = row0.x; m01 = row0.y; m02 = row0.z;
	m10 = row1.x; m11 = row1.y; m12 = row1.z;
	m20 = row2.x; m21 = row2.y; m22 = row2.z; 
	return *this;
}

template<typename T>
inline Matrix3<T>& Matrix3<T>::set(T m00, T m01, T m02, T m10, T m11, T m12, T m20, T m21, T m22)
{
	this->m00 = m00; this->m01 = m01; this->m02 = m02; 
	this->m10 = m10; this->m11 = m11; this->m12 = m12; 
	this->m20 = m20; this->m21 = m21; this->m22 = m22; 
	return *this;
}

template<typename T>
inline Matrix3<T>& Matrix3<T>::setScaling(const Vector3<T>& v)
{
	m00 = v.x; m01 = T(); m02 = T();
	m10 = T(); m11 = v.y; m12 = T();
	m20 = T(); m21 = T(); m22 = v.z;
	return *this;
}

template<typename T>
inline Matrix3<T>& Matrix3<T>::setScaling(T f)
{
	m00 = f; m01 = T(); m02 = T();
	m10 = T(); m11 = f; m12 = T();
	m20 = T(); m21 = T(); m22 = f;
	return *this;
}

template<typename T>
inline Matrix3<T>& Matrix3<T>::setRotation(Axis axis, T angle)
{
	if (angle != T(0))
	{
		T sinAngle = std::sin(angle);
		T cosAngle = std::cos(angle);

		switch (axis)
		{
			case Axis::X:
				return set(T(1), T(), T(), T(), cosAngle, sinAngle, T(), -sinAngle, cosAngle);
			case Axis::Y:
				return set(cosAngle, T(), -sinAngle, T(), T(1), T(), sinAngle, T(), cosAngle);
			case Axis::Z:
				return set(cosAngle, sinAngle, T(), -sinAngle, cosAngle, T(), T(), T(), T(1));
			default:
				return setZero();
		}
	}

	return setIdentity();
}

template<typename T>
inline Matrix3<T>& Matrix3<T>::setRotation(const Vector3<T>& axis, T angle)
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
			acm + bSin, bcm - aSin, cc*mCos + cosAngle);
	}

	return setIdentity();
}

template<typename T>
inline Matrix3<T>& Matrix3<T>::setShearing(T xy, T xz, T yx, T yz, T zx, T zy)
{
	m00 = T(1); m01 = xy; m02 = xz;
	m10 = yx; m11 = T(1); m12 = yz;
	m20 = zx; m21 = zy; m22 = T(1);
	return *this;
}

template<typename T>
inline Matrix3<T>& Matrix3<T>::setTranspose(const Matrix3<T>& m)
{
	m00 = m.m00; m01 = m.m10; m02 = m.m20;
	m10 = m.m01; m11 = m.m11; m12 = m.m21;
	m20 = m.m02; m21 = m.m12; m22 = m.m22;
	return *this;
}

template<typename T>
inline Matrix3<T>& Matrix3<T>::setInverse(const Matrix3<T>& m)
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
	return *this;
}

template<typename T>
inline Matrix3<T>& Matrix3<T>::setInverseTranspose(const Matrix3<T>& m)
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
	m00 = a*n; m01 = b*n; m02 = c*n;
	m10 = d*n; m11 = e*n; m12 = f*n;
	m20 = g*n; m21 = h*n; m22 = i*n;
	return *this;
}

template<typename T>
inline Matrix3<T>& Matrix3<T>::preConcatenate(const Matrix3<T>& m)
{
	set(m.m00*m00 + m.m01*m10 + m.m02*m20,
		m.m00*m01 + m.m01*m11 + m.m02*m21,
		m.m00*m02 + m.m01*m12 + m.m02*m22,
		m.m10*m00 + m.m11*m10 + m.m12*m20,
		m.m10*m01 + m.m11*m11 + m.m12*m21,
		m.m10*m02 + m.m11*m12 + m.m12*m22,
		m.m20*m00 + m.m21*m10 + m.m22*m20,
		m.m20*m01 + m.m21*m11 + m.m22*m21,
		m.m20*m02 + m.m21*m12 + m.m22*m22);
	return *this;
}

template<typename T>
inline Matrix3<T>& Matrix3<T>::preScale(const Vector3<T>& v)
{
	m00 *= v.x; m01 *= v.x; m02 *= v.x;
	m10 *= v.y; m11 *= v.y; m12 *= v.y;
	m20 *= v.z; m21 *= v.z; m22 *= v.z;
	return *this;
}

template<typename T>
inline Matrix3<T>& Matrix3<T>::scale(const Vector3<T>& v)
{
	m00 *= v.x; m01 *= v.y; m02 *= v.z;
	m10 *= v.x; m11 *= v.y; m12 *= v.z;
	m20 *= v.x; m21 *= v.y; m22 *= v.z;
	return *this;
}

template<typename T>
inline Matrix3<T>& Matrix3<T>::negate()
{
	m00 = -m00; m01 = -m01; m02 = -m02;
	m10 = -m10; m11 = -m11; m12 = -m12;
	m20 = -m20; m21 = -m21; m22 = -m22;
	return *this;
}

template<typename T>
inline Matrix3<T>& Matrix3<T>::transpose()
{
	T t = m01; m01 = m10; m10 = t;
	t = m02; m02 = m20; m20 = t;
	t = m12; m12 = m21; m21 = t;
	return *this;
}

template<typename T>
inline Matrix3<T>& Matrix3<T>::invert()
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
	return *this;
}

template<typename T>
inline Matrix3<T>& Matrix3<T>::orthonormalize()
{
	Vector3<T>& v1 = (*this)[0];
	Vector3<T>& v2 = (*this)[1];
	Vector3<T>& v3 = (*this)[2];
	v2 -= project(v2, v1);
	v3 -= project(v3, v1);
	v3 -= project(v3, v2);
	v1.normalize();
	v2.normalize();
	v3.normalize();
	return *this;
}

#if SIMD_HAS_FLOAT4

inline Matrix3<float>::Matrix3() : 
	row0(simd::zero<simd::float4>()), 
	row1(simd::zero<simd::float4>()), 
	row2(simd::zero<simd::float4>()) 
{
}

inline Matrix3<float>::Matrix3(Identity) : 
	row0(Vector3<float>::UNIT_X), 
	row1(Vector3<float>::UNIT_Y), 
	row2(Vector3<float>::UNIT_Z) 
{
}

inline Matrix3<float>::Matrix3(float m00, float m01, float m02, float m10, float m11, float m12, float m20, float m21, float m22) :
#if MATHEMATICS_SIMD_EXPAND_LAST
	row0(simd::set4(m00, m01, m02, m02)),
	row1(simd::set4(m10, m11, m12, m12)),
	row2(simd::set4(m20, m21, m22, m22))
#else
	row0(simd::set3(m00, m01, m02)),
	row1(simd::set3(m10, m11, m12)),
	row2(simd::set3(m20, m21, m22))
#endif
{
}

inline Matrix3<float>::Matrix3(const Vector3<float>& row0, const Vector3<T>& row1, const Vector3<T>& row2) :
	row0(row0),
	row1(row1),
	row2(row2)
{
}

inline Matrix3<float>::Matrix3(const std::tuple<Vector3<float>, Vector3<float>, Vector3<float>>& t) :
	Matrix3(std::get<0>(t), std::get<1>(t), std::get<2>(t))
{
}

inline Matrix3<float>::Matrix3(const Matrix2<float>& m) :
	row0(simd::cutoff2(m.row0)),
	row1(simd::cutoff2(m.row1)),
	row2(Vector3<float>::UNIT_Z)
{
}

inline Matrix3<float>::Matrix3(const Vector3<float>& forward)
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
			set(right, up, at);
		}
		else
		{
			Vector3<float> up((at.y > 0.f) ? -Vector3<float>::UNIT_Z : Vector3<float>::UNIT_Z);
			Vector3<float> right = cross(up, at);
			right.normalize();
			up = cross(at, right);
			set(right, up, at);
		}
	}
	else
	{
		setIdentity();
	}
}

inline Matrix3<float>::Matrix3(const Vector3<float>& up, const Vector3<float>& forward)
{
	float m = forward.getMagnitude();
	if ((m > 0.f) && (up.getMagnitude() > 0.f))
	{
		Vector3<float> at = forward/m;
		Vector3<float> right = cross(up, at);
		right.normalize();
		set(right, cross(at, right), at);
	}
	else
	{
		setIdentity();
	}
}

#if MATHEMATICS_SIMD_EXPAND_LAST
inline Matrix3<float>::Matrix3(const float* m) :
	Matrix3(m[0], m[1], m[2], m[3], m[4], m[5], m[6], m[7], m[8])
{
}
#else
inline Matrix3<float>::Matrix3(const float* m)
{
	simd::unpack3x3(m, row0, row1, row2);
}
#endif

inline Matrix3<float>::Matrix3(const std::tuple<simd::float4, simd::float4, simd::float4>& t) :
	row0(std::get<0>(t)), 
	row1(std::get<1>(t)), 
	row2(std::get<2>(t))
{
}

inline Matrix3<float> Matrix3<float>::operator-() const
{
#if MATHEMATICS_SIMD_EXPAND_LAST
	return Matrix3<float>(simd::neg4(row0), simd::neg4(row1), simd::neg4(row2));
#else
	return Matrix3<float>(simd::neg3(row0), simd::neg3(row1), simd::neg3(row2));
#endif
}

inline Matrix3<float>& Matrix3<float>::operator+=(const Matrix3<float>& m)
{
	row0 = simd::add4(row0, m.row0);
	row1 = simd::add4(row1, m.row1);
	row2 = simd::add4(row2, m.row2);
	return *this;
}

inline Matrix3<float>& Matrix3<float>::operator-=(const Matrix3<float>& m)
{
	row0 = simd::sub4(row0, m.row0);
	row1 = simd::sub4(row1, m.row1);
	row2 = simd::sub4(row2, m.row2);
	return *this;
}

inline Matrix3<float>& Matrix3<float>::operator*=(float f)
{
	auto t = simd::set4(f);
	row0 = simd::mul4(row0, t);
	row1 = simd::mul4(row1, t);
	row2 = simd::mul4(row2, t);
	return *this;
}

inline Matrix3<float>& Matrix3<float>::operator*=(const Matrix3<float>& m)
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
	row0 = r0;
	row1 = r1;
	row2 = r2;
	return *this;
}

template<>
inline Matrix3<float> operator+(const Matrix3<float>& m1, const Matrix3<float>& m2) noexcept
{
	return Matrix3<float>(simd::add4(m1.row0, m2.row0),
		simd::add4(m1.row1, m2.row1),
		simd::add4(m1.row2, m2.row2));
}

template<>
inline Matrix3<float> operator-(const Matrix3<float>& m1, const Matrix3<float>& m2) noexcept
{
	return Matrix3<float>(simd::sub4(m1.row0, m2.row0),
		simd::sub4(m1.row1, m2.row1),
		simd::sub4(m1.row2, m2.row2));
}

template<>
inline Matrix3<float> operator*(float f, const Matrix3<float>& m) noexcept
{
	auto t = simd::set4(f);
	return Matrix3<float>(simd::mul4(t, m.row0),
		simd::mul4(t, m.row1),
		simd::mul4(t, m.row2));
}

template<>
inline Matrix3<float> operator*(const Matrix3<float>& m, float f) noexcept
{
	auto t = simd::set4(f);
	return Matrix3<float>(simd::mul4(m.row0, t),
		simd::mul4(m.row1, t),
		simd::mul4(m.row2, t));
}

template<>
inline Matrix3<float> operator*(const Matrix3<float>& m1, const Matrix3<float>& m2) noexcept
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
	return Matrix3<float>(r0, r1, r2);
}

template<>
inline Matrix3<float> operator/(const Matrix3<float>& m, float f) noexcept 
{ 
	return operator*(m, 1.f/f); 
}

inline bool Matrix3<float>::operator==(const Matrix3<float>& m) const
{
	return simd::all3(simd::equal(row0, m.row0)) &&
		simd::all3(simd::equal(row1, m.row1)) &&
		simd::all3(simd::equal(row2, m.row2));
}

template<typename C, typename T>
inline std::basic_istream<C, T>& operator>>(std::basic_istream<C, T>& s, Matrix3<float>& m)
{
	float m00, m01, m02;
	float m10, m11, m12;
	float m20, m21, m22;
	s >> m00 >> std::ws >> m01 >> std::ws >> m02 >> std::ws >>
		m10 >> std::ws >> m11 >> std::ws >> m12 >> std::ws >>
		m20 >> std::ws >> m21 >> std::ws >> m22;
	m.set(m00, m01, m02, m10, m11, m12, m20, m21, m22);
	return s;
}

template<typename A>
inline void Matrix3<float>::load(A& ar)
{
	float t00, t01, t02;
	float t10, t11, t12;
	float t20, t21, t22;
	ar(t00, t01, t02, t10, t11, t12, t20, t21, t22);
	set(t00, t01, t02, t10, t11, t12, t20, t21, t22);
}

inline bool Matrix3<float>::isZero() const
{
	const auto zero = simd::zero<simd::float4>();
	return simd::all3(simd::equal(row0, zero)) && 
		simd::all3(simd::equal(row1, zero)) && 
		simd::all3(simd::equal(row2, zero));
}

inline bool Matrix3<float>::isApproxZero() const
{
	return simd::all3(simd::lessThan(simd::abs4(row0), Vector3<float>::TOLERANCE)) &&
		simd::all3(simd::lessThan(simd::abs4(row1), Vector3<float>::TOLERANCE)) &&
		simd::all3(simd::lessThan(simd::abs4(row2), Vector3<float>::TOLERANCE));
}

inline bool Matrix3<float>::isIdentity() const
{
	return simd::all3(simd::equal(row0, Vector3<float>::UNIT_X)) &&
		simd::all3(simd::equal(row1, Vector3<float>::UNIT_Y)) &&
		simd::all3(simd::equal(row2, Vector3<float>::UNIT_Z));
}

inline bool Matrix3<float>::isApproxIdentity() const
{
	return simd::all3(simd::lessThan(simd::abs4(simd::sub4(row0, Vector3<float>::UNIT_X)), Vector3<float>::TOLERANCE)) &&
		simd::all3(simd::lessThan(simd::abs4(simd::sub4(row1, Vector3<float>::UNIT_Y)), Vector3<float>::TOLERANCE)) &&
		simd::all3(simd::lessThan(simd::abs4(simd::sub4(row2, Vector3<float>::UNIT_Z)), Vector3<float>::TOLERANCE));
}

inline bool Matrix3<float>::isApproxEqual(const Matrix3& m) const
{
	return simd::all3(simd::lessThan(simd::abs4(simd::sub4(row0, m.row0)), Vector3<float>::TOLERANCE)) &&
		simd::all3(simd::lessThan(simd::abs4(simd::sub4(row1, m.row1)), Vector3<float>::TOLERANCE)) &&
		simd::all3(simd::lessThan(simd::abs4(simd::sub4(row2, m.row2)), Vector3<float>::TOLERANCE));
}

inline bool Matrix3<float>::isApproxEqual(const Matrix3& m, float tolerance) const
{
	auto t = simd::set4(tolerance);
	return simd::all3(simd::lessThan(simd::abs4(simd::sub4(row0, m.row0)), t)) &&
		simd::all3(simd::lessThan(simd::abs4(simd::sub4(row1, m.row1)), t)) &&
		simd::all3(simd::lessThan(simd::abs4(simd::sub4(row2, m.row2)), t));
}

inline bool Matrix3<float>::isApproxOrthogonal() const
{
	Matrix3<float> m(*this);
	m *= Matrix3<float>(Uninitialized()).setTranspose(*this);
	return m.isApproxIdentity();
}

inline bool Matrix3<float>::hasApproxUniformScaling() const
{
	Matrix3<float> m(*this);
	m *= Matrix3<float>(Uninitialized()).setTranspose(*this);
	return simd::all3(simd::lessThan(simd::abs4(simd::set3(m.m01, m.m02, m.m10)), Vector3<float>::TOLERANCE)) &&
		simd::all3(simd::lessThan(simd::abs4(simd::set3(m.m12, m.m20, m.m21)), Vector3<float>::TOLERANCE)) &&
		simd::all3(simd::lessThan(simd::abs4(simd::sub4(simd::set3(m.m11, m.m22, m.m00), simd::set3(m.m00, m.m11, m.m22))), Vector3<float>::TOLERANCE));
}

inline bool Matrix3<float>::isFinite() const
{
	return simd::all3(simd::isFinite(row0)) &&
		simd::all3(simd::isFinite(row1)) &&
		simd::all3(simd::isFinite(row2));
}

inline float Matrix3<float>::getDeterminant() const // #TODO SIMD
{
	return m00*(m11*m22 - m12*m21) - m01*(m10*m22 - m12*m20) + m02*(m10*m21 - m11*m20);
}

inline Matrix3<float>& Matrix3<float>::setZero()
{
	row0 = simd::zero<simd::float4>();
	row1 = simd::zero<simd::float4>();
	row2 = simd::zero<simd::float4>();
	return *this;
}

inline Matrix3<float>& Matrix3<float>::setIdentity()
{
	row0 = Vector3<float>::UNIT_X;
	row1 = Vector3<float>::UNIT_Y;
	row2 = Vector3<float>::UNIT_Z;
	return *this;
}

inline Matrix3<float>& Matrix3<float>::set(const Vector3<float>& row0, const Vector3<float>& row1, const Vector3<float>& row2)
{
	this->row0 = row0;
	this->row1 = row1;
	this->row2 = row2;
	return *this;
}

inline Matrix3<float>& Matrix3<float>::set(simd::float4 row0, simd::float4 row1, simd::float4 row2)
{
	this->row0 = row0;
	this->row1 = row1;
	this->row2 = row2;
	return *this;
}

inline Matrix3<float>& Matrix3<float>::set(float m00, float m01, float m02, float m10, float m11, float m12,
	float m20, float m21, float m22)
{
#if MATHEMATICS_SIMD_EXPAND_LAST
	row0 = simd::set4(m00, m01, m02, m02);
	row1 = simd::set4(m10, m11, m12, m12);
	row2 = simd::set4(m20, m21, m22, m22);
#else
	row0 = simd::set3(m00, m01, m02);
	row1 = simd::set3(m10, m11, m12);
	row2 = simd::set3(m20, m21, m22);
#endif
	return *this;
}

inline Matrix3<float>& Matrix3<float>::setScaling(const Vector3<float>& v)
{
	row0 = simd::cutoff1(v);
	row1 = simd::and4(v, simd::constant4i<simd::float4, 0, -1, 0, 0>());
#if MATHEMATICS_SIMD_EXPAND_LAST
	row2 = simd::and4(v, simd::constant4i<simd::float4, 0, 0, -1, -1>());
#else
	row2 = simd::and4(v, simd::constant4i<simd::float4, 0, 0, -1, 0>());
#endif
	return *this;
}

inline Matrix3<float>& Matrix3<float>::setScaling(float f)
{
	const auto zero = simd::zero<simd::float4>();
	row0 = simd::set1(f);
	row1 = simd::insert<1>(f, zero);
#if MATHEMATICS_SIMD_EXPAND_LAST
	row2 = simd::set4(0.f, 0.f, f, f);
#else
	row2 = simd::insert<2>(f, zero);
#endif
	return *this;
}

inline Matrix3<float>& Matrix3<float>::setRotation(Axis axis, float angle)
{
	if (angle != 0.f)
	{
		float sinAngle = std::sin(angle);
		float cosAngle = std::cos(angle);

		switch (axis)
		{
			case Axis::X:
				return set(1.f, 0.f, 0.f, 0.f, cosAngle, sinAngle, 0.f, -sinAngle, cosAngle);
			case Axis::Y:
				return set(cosAngle, 0.f, -sinAngle, 0.f, 1.f, 0.f, sinAngle, 0.f, cosAngle);
			case Axis::Z:
				return set(cosAngle, sinAngle, 0.f, -sinAngle, cosAngle, 0.f, 0.f, 0.f, 1.f);
			default:
				return setZero();
		}
	}

	return setIdentity();
}

inline Matrix3<float>& Matrix3<float>::setRotation(const Vector3<float>& axis, float angle)
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
			acm + bSin, bcm - aSin, cc*mCos + cosAngle);
	}

	return setIdentity();
}

inline Matrix3<float>& Matrix3<float>::setShearing(float xy, float xz, float yx, float yz, float zx, float zy)
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
	return *this;
}

inline Matrix3<float>& Matrix3<float>::setTranspose(const Matrix3<float>& m)
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
	return *this;
}

inline Matrix3<float>& Matrix3<float>::setInverse(const Matrix3<float>& m)
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
	return *this;
}

inline Matrix3<float>& Matrix3<float>::setInverseTranspose(const Matrix3<float>& m)
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
	row0 = simd::mul4(simd::set4(a, b, c, c), n);
	row1 = simd::mul4(simd::set4(d, e, f, f), n);
	row2 = simd::mul4(simd::set4(g, h, i, i), n);
#else
	row0 = simd::mul4(simd::set3(a, b, c), n);
	row1 = simd::mul4(simd::set3(d, e, f), n);
	row2 = simd::mul4(simd::set3(g, h, i), n);
#endif
	return *this;
}

inline Matrix3<float>& Matrix3<float>::preConcatenate(const Matrix3<float>& m)
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
	row0 = r0;
	row1 = r1;
	row2 = r2;
	return *this;
}

inline Matrix3<float>& Matrix3<float>::preScale(const Vector3<float>& v)
{
	row0 = simd::mul4(row0, simd::xxxx(v));
	row1 = simd::mul4(row1, simd::yyyy(v));
	row2 = simd::mul4(row2, simd::zzzz(v));
	return *this;
}

inline Matrix3<float>& Matrix3<float>::scale(const Vector3<float>& v)
{
	row0 = simd::mul4(row0, v);
	row1 = simd::mul4(row1, v);
	row2 = simd::mul4(row2, v);
	return *this;
}

inline Matrix3<float>& Matrix3<float>::negate()
{
#if MATHEMATICS_SIMD_EXPAND_LAST
	row0 = simd::neg4(row0);
	row1 = simd::neg4(row1);
	row2 = simd::neg4(row2);
#else
	row0 = simd::neg3(row0);
	row1 = simd::neg3(row1);
	row2 = simd::neg3(row2);
#endif
	return *this;
}

inline Matrix3<float>& Matrix3<float>::transpose()
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
	return *this;
}

inline Matrix3<float>& Matrix3<float>::invert()
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
	return *this;
}

inline Matrix3<float>& Matrix3<float>::orthonormalize()
{
	Vector3<float>& v1 = (*this)[0];
	Vector3<float>& v2 = (*this)[1];
	Vector3<float>& v3 = (*this)[2];
	v2 -= project(v2, v1);
	v3 -= project(v3, v1);
	v3 -= project(v3, v2);
	v1.normalize();
	v2.normalize();
	v3.normalize();
	return *this;
}

#endif /* SIMD_HAS_FLOAT4 */

template<typename T>
	requires std::floating_point<T>
inline Matrix3<T> concatenate(const Matrix3<T>& m1, const Matrix3<T>& m2) noexcept
{
	return m1*m2;
}

template<typename T>
	requires std::floating_point<T>
inline Matrix3<T> concatenate(const Matrix3<T>& m1, const Matrix3<T>& m2, const Matrix3<T>& m3) noexcept
{
	return concatenate(concatenate(m1, m2), m3);
}

template<typename T>
	requires std::floating_point<T>
inline Matrix3<T> concatenate(const Matrix3<T>& m1, const Matrix3<T>& m2, const Matrix3<T>& m3, const Matrix3<T>& m4) noexcept
{
	return concatenate(concatenate(concatenate(m1, m2), m3), m4);
}

template<typename T>
	requires std::floating_point<T>
inline Matrix3<T> transpose(const Matrix3<T>& m) noexcept
{
	return Matrix3<T>(m.m00, m.m10, m.m20,
		m.m01, m.m11, m.m21,
		m.m02, m.m12, m.m22);
}

template<typename T>
	requires std::floating_point<T>
inline Matrix3<T> transpose(Matrix3<T>&& m) noexcept
{
	m.transpose();
	return m;
}

template<typename T>
	requires std::floating_point<T>
inline Matrix3<T> inverse(const Matrix3<T>& m) noexcept
{
	return Matrix3<T>(Uninitialized()).setInverse(m);
}

template<typename T>
	requires std::floating_point<T>
inline Matrix3<T> inverse(Matrix3<T>&& m) noexcept
{
	m.invert();
	return m;
}

template<typename T>
	requires std::floating_point<T>
inline Matrix3<T> adjoint(const Matrix3<T>& m) noexcept
{
	auto det2 = [](T m00, T m01, T m10, T m11) -> T { return m00*m11 - m01*m10; };
	return Matrix3<T>(det2(m.m11, m.m12, m.m21, m.m22), -det2(m.m01, m.m02, m.m21, m.m22), det2(m.m01, m.m02, m.m11, m.m12),
		-det2(m.m10, m.m12, m.m20, m.m22), det2(m.m00, m.m02, m.m20, m.m22), -det2(m.m00, m.m02, m.m10, m.m12),
		det2(m.m10, m.m11, m.m20, m.m21), -det2(m.m00, m.m01, m.m20, m.m21), det2(m.m00, m.m01, m.m10, m.m11));
}

#if SIMD_HAS_FLOAT4

template<>
inline Matrix3<float> transpose(const Matrix3<float>& m) noexcept
{
	Matrix3<float> n(Uninitialized());
#if MATHEMATICS_SIMD_EXPAND_LAST
	auto [r0, r1, r2] = simd::transpose3x3(m.row0, m.row1, m.row2);
	n.row0 = simd::xyzz(r0);
	n.row1 = simd::xyzz(r1);
	n.row2 = simd::xyzz(r2);
#else
	//simd::transpose3x3(m.row0, m.row1, m.row2, n.row0, n.row1, n.row2);
	std::tie(n.row0, n.row1, n.row2) = simd::transpose3x3(m.row0, m.row1, m.row2);
#endif
	return n;
}

template<>
inline Matrix3<float> transpose(Matrix3<float>&& m) noexcept
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
	return m;
}

#endif /* SIMD_HAS_FLOAT4 */

} // namespace templates

using Matrix3f = templates::Matrix3<float>;
using Matrix3d = templates::Matrix3<double>;

#if MATHEMATICS_DOUBLE
using Matrix3 = templates::Matrix3<double>;
using Matrix3Arg = templates::Matrix3<double>::ConstArg;
using Matrix3Result = templates::Matrix3<double>::ConstResult;
#else
using Matrix3 = templates::Matrix3<float>;
using Matrix3Arg = templates::Matrix3<float>::ConstArg;
using Matrix3Result = templates::Matrix3<float>::ConstResult;
#endif

} // namespace core::mathematics

#include "Quaternion.hpp"
#include "YawPitchRoll.hpp"
#include "Euler.hpp"

namespace core::mathematics::templates {

template<typename T>
inline Matrix3<T>::Matrix3(const Quaternion<T>& q)
{
	if (!q.isIdentity())
	{
		//T nq = q.getNorm();
		const T s = /*(nq > T(0)) ? T(2)/nq :*/ T(2);
		T xs = q.x*s, ys = q.y*s, zs = q.z*s;
		T wx = q.w*xs, wy = q.w*ys, wz = q.w*zs;
		T xx = q.x*xs, xy = q.x*ys, xz = q.x*zs;
		T yy = q.y*ys, yz = q.y*zs, zz = q.z*zs;

		set(T(1) - (yy + zz), xy + wz, xz - wy,
			xy - wz, T(1) - (xx + zz), yz + wx,
			xz + wy, yz - wx, T(1) - (xx + yy));
	}
	else
		setIdentity();
}

template<typename T>
inline Matrix3<T>& Matrix3<T>::setRotation(const YawPitchRoll<T>& r)
{
	if (!r.isZero())
	{
		T sy = std::sin(r.yaw), cy = std::cos(r.yaw);
		T sp = std::sin(r.pitch), cp = std::cos(r.pitch);
		T sr = std::sin(r.roll), cr = std::cos(r.roll);
		T ss = sy*sp, cs = cy*sp;

		return set(cr*cy + sr*ss, sr*cp, sr*cs - cr*sy,
			cr*ss - sr*cy, cr*cp, cr*cs + sr*sy,
			sy*cp, -sp, cy*cp);
	}

	return setIdentity();
}

template<typename T>
inline Matrix3<T>& Matrix3<T>::setRotation(const Euler<T>& e)
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
			m[2][0], m[2][1], m[2][2]);
	}

	return setIdentity();
}

template<typename T>
inline Matrix3<T>& Matrix3<T>::setRotation(const Quaternion<T>& q)
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
			xz + wy, yz - wx, T(1) - (xx + yy));
	}

	return setIdentity();
}

#if SIMD_HAS_FLOAT4

inline Matrix3<float>::Matrix3(const Quaternion<float>& q)
{
	if (!q.isIdentity())
	{
		//float nq = q.getNorm();
		const float s = /*(nq > 0.f) ? 2.f/nq :*/ 2.f;
		float xs = q.x*s, ys = q.y*s, zs = q.z*s; // #TODO SIMD
		float wx = q.w*xs, wy = q.w*ys, wz = q.w*zs;
		float xx = q.x*xs, xy = q.x*ys, xz = q.x*zs;
		float yy = q.y*ys, yz = q.y*zs, zz = q.z*zs;

		set(1.f - (yy + zz), xy + wz, xz - wy,
			xy - wz, 1.f - (xx + zz), yz + wx,
			xz + wy, yz - wx, 1.f - (xx + yy));
	}
	else
		setIdentity();
}

inline Matrix3<float>& Matrix3<float>::setRotation(const YawPitchRoll<float>& r)
{
	if (!r.isZero())
	{
		float sy = std::sin(r.yaw), cy = std::cos(r.yaw);
		float sp = std::sin(r.pitch), cp = std::cos(r.pitch);
		float sr = std::sin(r.roll), cr = std::cos(r.roll);
		float ss = sy*sp, cs = cy*sp;

		return set(cr*cy + sr*ss, sr*cp, sr*cs - cr*sy,
			cr*ss - sr*cy, cr*cp, cr*cs + sr*sy,
			sy*cp, -sp, cy*cp);
	}

	return setIdentity();
}

inline Matrix3<float>& Matrix3<float>::setRotation(const Euler<float>& e)
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
			m[2][0], m[2][1], m[2][2]);
	}

	return setIdentity();
}

inline Matrix3<float>& Matrix3<float>::setRotation(const Quaternion<float>& q)
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
			xz + wy, yz - wx, 1.f - (xx + yy));
	}

	return setIdentity();
}

#endif /* SIMD_HAS_FLOAT4 */

} // namespace core::mathematics::templates
