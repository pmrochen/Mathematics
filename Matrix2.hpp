/*
 *	Name: Matrix2
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
#include "Constants.hpp"
#include "Vector2.hpp"
#include "Vector3.hpp"
#include "Vector4.hpp"

namespace core::mathematics {
	
struct Identity
{
};

//constexpr Identity IDENTITY{};

namespace templates {

template<typename T>
	requires std::floating_point<T>
struct Matrix2
{
	using Real = T;
	using ComponentType = T;
	using ConstArg = const Matrix2&;
	using ConstResult = const Matrix2&;

	//static constexpr int NUM_COMPONENTS = 4;

	constexpr Matrix2() noexcept : m00(), m01(), m10(), m11() {}
	explicit Matrix2(Uninitialized) noexcept {}
	explicit Matrix2(Identity) noexcept : m00(1), m01(), m10(), m11(1) {}
	constexpr Matrix2(T m00, T m01, T m10, T m11) noexcept : m00(m00), m01(m01), m10(m10), m11(m11) {}
	constexpr Matrix2(const Vector2<T>& row0, const Vector2<T>& row1) noexcept : m00(row0.x), m01(row0.y), m10(row1.x), m11(row1.y) {}
	constexpr explicit Matrix2(const std::pair<Vector2<T>, Vector2<T>>& t) noexcept : Matrix2(t.first, t.second) {}
	constexpr explicit Matrix2(const std::tuple<Vector2<T>, Vector2<T>>& t) noexcept : Matrix2(std::get<0>(t), std::get<1>(t)) {}
	explicit Matrix2(const T* m) noexcept : m00(m[0]), m01(m[1]), m10(m[2]), m11(m[3]) {}

	//explicit operator T*() noexcept { return &m00; }
	//explicit operator const T*() const noexcept { return &m00; }
	Vector2<T>& operator[](int i) noexcept { return reinterpret_cast<Vector2<T>*>(&m00)[i]; }
	const Vector2<T>& operator[](int i) const noexcept { return reinterpret_cast<const Vector2<T>*>(&m00)[i]; }

	Matrix2 operator+() const noexcept { return *this; }
	Matrix2 operator-() const noexcept { return Matrix2(-m00, -m01, -m10, -m11); }
	Matrix2& operator+=(const Matrix2& m) noexcept;
	Matrix2& operator-=(const Matrix2& m) noexcept;
	Matrix2& operator*=(T f) noexcept;
	Matrix2& operator/=(T f) noexcept { return operator*=(T(1)/f); }
	Matrix2& operator*=(const Matrix2& m) noexcept;
	bool operator==(const Matrix2& m) const noexcept { return (m00 == m.m00) && (m01 == m.m01) && (m10 == m.m10) && (m11 == m.m11); }
	bool operator!=(const Matrix2& m) const noexcept { return !(*this == m); }

	template<typename A> void load(A& ar) { ar(m00, m01, m10, m11); }
	template<typename A> void save(A& ar) const { ar(m00, m01, m10, m11); }

	// #TODO rename get...() to make...()
	static Matrix2 getScaling(const Vector2<T>& v) noexcept { return Matrix2(Uninitialized()).setScaling(v); }
	static Matrix2 getScaling(T f) noexcept { return Matrix2(Uninitialized()).setScaling(f); }
	static Matrix2 getRotation(T angle) noexcept { return Matrix2(Uninitialized()).setRotation(angle); }
	static Matrix2 getScalingRotation(const Vector2<T>& s, T angle) noexcept { return Matrix2(Uninitialized()).setScalingRotation(s, angle); }
	static Matrix2 getShearing(T xy, T yx) noexcept { return Matrix2(Uninitialized()).setShearing(xy, yx); }

	bool isZero() const noexcept { return (m00 == T()) && (m01 == T()) && (m10 == T()) && (m11 == T()); }
	bool isApproxZero() const noexcept;
	bool isIdentity() const noexcept { return (m00 == T(1)) && (m01 == T()) && (m10 == T()) && (m11 == T(1)); }
	bool isApproxIdentity() const noexcept;
	bool isApproxEqual(const Matrix2& m) const noexcept;
	bool isApproxEqual(const Matrix2& m, T tolerance) const noexcept;
	bool isApproxOrthogonal() const noexcept;
	bool hasApproxUniformScaling() const noexcept;
	bool isFinite() const noexcept { return std::isfinite(m00) && std::isfinite(m01) && std::isfinite(m10) && std::isfinite(m11); }
	const Vector2<T>& getRow(int i) const noexcept { return reinterpret_cast<const Vector2<T>*>(&m00)[i]; }
	Vector2<T> getColumn(int i) const noexcept { return Vector2<T>((&m00)[i], (&m10)[i]); }
	Vector2<T> getScaleComponent() const noexcept { return Vector2<T>((*this)[0].getMagnitude(), (*this)[1].getMagnitude()); }
	T getTrace() const noexcept { return (m00 + m11); }
	T getDeterminant() const noexcept { return (m00*m11 - m01*m10); }
	bool isSingular() const noexcept { return (getDeterminant() == T(0)); }
	Matrix2& setZero/*zero*/() noexcept { m00 = T(); m01 = T(); m10 = T(); m11 = T(); return *this; }
	Matrix2& setIdentity/*makeIdentity*/() noexcept { m00 = T(1); m01 = T(); m10 = T(); m11 = T(1); return *this; }
	Matrix2& set(const Vector2<T>& row0, const Vector2<T>& row1) noexcept;
	Matrix2& set(T m00, T m01, T m10, T m11) noexcept;
	Matrix2& setScaling/*makeScaling*/(const Vector2<T>& v) noexcept;
	Matrix2& setScaling/*makeScaling*/(T f) noexcept;
	Matrix2& setRotation(T angle) noexcept;
	Matrix2& setScalingRotation(const Vector2<T>& s, T angle) noexcept { return setRotation(angle).preScale(s); }
	Matrix2& setShearing(T xy, T yx) noexcept;
	Matrix2& setTranspose/*transposeOf*/(const Matrix2& m) noexcept;
	Matrix2& setInverse/*inverseOf*/(const Matrix2& m) noexcept;
	Matrix2& setInverseTranspose/*inverseTransposeOf*/(const Matrix2& m) noexcept;
	Matrix2& preConcatenate(const Matrix2& m) noexcept;
	Matrix2& concatenate(const Matrix2& m) noexcept { return operator*=(m); }
	Matrix2& preScale(const Vector2<T>& v) noexcept;
	Matrix2& preScale(T f) noexcept { return operator*=(f); }
	Matrix2& scale(const Vector2<T>& v) noexcept;
	Matrix2& scale(T f) noexcept { return operator*=(f); }
	//Matrix2& rotate(T angle) noexcept { concatenate(Matrix2(Uninitialized()).setRotation(angle)); return *this; }
	//Matrix2& shear(T xy, T yx) noexcept { concatenate(Matrix2(Uninitialized()).setShearing(xy, yx)); return *this; }
	Matrix2& negate() noexcept;
	Matrix2& transpose() noexcept;
	Matrix2& invert() noexcept;
	Matrix2& orthonormalize() noexcept;

	//static const Matrix2& getZero() noexcept { return ZERO; }
	//static const Matrix2& getIdentity() noexcept { return IDENTITY; }

	static const Matrix2 ZERO;
	static const Matrix2 IDENTITY;

	T m00, m01;
	T m10, m11;
};

template<typename T> const Matrix2<T> Matrix2<T>::ZERO{};
template<typename T> const Matrix2<T> Matrix2<T>::IDENTITY{ T(1), T(0), T(0), T(1) };

#if SIMD_HAS_FLOAT4

template<>
struct Matrix2<float>
{
	using Real = float;
	using ComponentType = float;
	using ConstArg = const Matrix2;
	using ConstResult = const Matrix2;

	//static constexpr int NUM_COMPONENTS = 4;

	/*constexpr*/ Matrix2() noexcept : row0(simd::zero<simd::float4>()), row1(simd::zero<simd::float4>()) {}
	explicit Matrix2(Uninitialized) noexcept {}
	explicit Matrix2(Identity) noexcept : row0(Vector2<float>::UNIT_X), row1(Vector2<float>::UNIT_Y) {}
	/*constexpr*/ Matrix2(float m00, float m01, float m10, float m11) noexcept;
	/*constexpr*/ Matrix2(const Vector2<float>& row0, const Vector2<float>& row1) noexcept : row0(row0), row1(row1) {}
	/*constexpr*/ explicit Matrix2(const std::pair<Vector2<float>, Vector2<float>>& t) noexcept : Matrix2(t.first, t.second) {}
	/*constexpr*/ explicit Matrix2(const std::tuple<Vector2<float>, Vector2<float>>& t) noexcept : Matrix2(std::get<0>(t), std::get<1>(t)) {}
	explicit Matrix2(const float* m) noexcept;
	explicit Matrix2(const simd::float4* m) noexcept : row0(m[0]), row1(m[1]) {}
	Matrix2(simd::float4 row0, simd::float4 row1) noexcept : row0(row0), row1(row1) {}
	explicit Matrix2(simd::float4 m) noexcept;
	explicit Matrix2(const std::pair<simd::float4, simd::float4>& t) noexcept : row0(t.first), row1(t.second) {}
	explicit Matrix2(const std::tuple<simd::float4, simd::float4>& t) noexcept : row0(std::get<0>(t)), row1(std::get<1>(t)) {}
	Matrix2(const Matrix2& m) noexcept : row0(m.row0), row1(m.row1) {}
	Matrix2& operator=(const Matrix2& m) noexcept { row0 = m.row0; row1 = m.row1; return *this; }

	//explicit operator float*() noexcept { return &m00; }
	//explicit operator const float*() const noexcept { return &m00; }
	explicit operator simd::float4*() noexcept { return &row0; }
	explicit operator const simd::float4*() const noexcept { return &row0; }
	explicit operator simd::float4() const noexcept { return simd::pack2x2(row0, row1); }
	Vector2<float>& operator[](int i) noexcept { return reinterpret_cast<Vector2<float>&>((&row0)[i]); }
	const Vector2<float>& operator[](int i) const noexcept { return reinterpret_cast<const Vector2<float>&>((&row0)[i]); }

	Matrix2 operator+() const noexcept { return *this; }
	Matrix2 operator-() const noexcept;
	Matrix2& operator+=(const Matrix2& m) noexcept;
	Matrix2& operator-=(const Matrix2& m) noexcept;
	Matrix2& operator*=(float f) noexcept;
	Matrix2& operator/=(float f) noexcept { return operator*=(1.f/f); }
	Matrix2& operator*=(const Matrix2& m) noexcept;
	bool operator==(const Matrix2& m) const noexcept;
	bool operator!=(const Matrix2& m) const noexcept { return !(*this == m); }

	template<typename A> void load(A& ar);
	template<typename A> void save(A& ar) const { ar(m00, m01, m10, m11); }

	// #TODO rename get...() to make...()
	static Matrix2 getScaling(const Vector2<float>& v) noexcept { return Matrix2(Uninitialized()).setScaling(v); }
	static Matrix2 getScaling(float f) noexcept { return Matrix2(Uninitialized()).setScaling(f); }
	static Matrix2 getRotation(float angle) noexcept { return Matrix2(Uninitialized()).setRotation(angle); }
	static Matrix2 getScalingRotation(const Vector2<float>& s, float angle) noexcept { return Matrix2(Uninitialized()).setScalingRotation(s, angle); }
	static Matrix2 getShearing(float xy, float yx) noexcept { return Matrix2(Uninitialized()).setShearing(xy, yx); }

	bool isZero() const noexcept;
	bool isApproxZero() const noexcept;
	bool isIdentity() const noexcept;
	bool isApproxIdentity() const noexcept;
	bool isApproxEqual(const Matrix2& m) const noexcept;
	bool isApproxEqual(const Matrix2& m, float tolerance) const noexcept;
	bool isApproxOrthogonal() const noexcept;
	bool hasApproxUniformScaling() const noexcept;
	bool isFinite() const noexcept { return simd::all4(simd::isFinite(simd::pack2x2(row0, row1))); }
	const Vector2<float>& getRow(int i) const noexcept { return reinterpret_cast<const Vector2<float>&>((&row0)[i]); }
	Vector2<float> getColumn(int i) const noexcept; // #TODO
	Vector2<float> getScaleComponent() const noexcept { return Vector2<float>((*this)[0].getMagnitude(), (*this)[1].getMagnitude()); }
	float getTrace() const noexcept { return (m00 + m11); }
	float getDeterminant() const noexcept { return (m00*m11 - m01*m10); }
	bool isSingular() const noexcept { return (getDeterminant() == 0.f); }
	Matrix2& setZero/*zero*/() noexcept { row0 = simd::zero<simd::float4>(); row1 = simd::zero<simd::float4>(); return *this; }
	Matrix2& setIdentity/*makeIdentity*/() noexcept { row0 = Vector2<float>::UNIT_X; row1 = Vector2<float>::UNIT_Y; return *this; }
	Matrix2& set(const Vector2<float>& row0, const Vector2<float>& row1) noexcept;
	Matrix2& set(simd::float4 row0, simd::float4 row1) noexcept;
	//Matrix2& set(simd::float4 packed) noexcept;
	Matrix2& set(float m00, float m01, float m10, float m11) noexcept;
	Matrix2& setScaling/*makeScaling*/(const Vector2<float>& v) noexcept;
	Matrix2& setScaling/*makeScaling*/(float f) noexcept;
	Matrix2& setRotation(float angle) noexcept;
	Matrix2& setScalingRotation(const Vector2<float>& s, float angle) noexcept { return setRotation(angle).preScale(s); }
	Matrix2& setShearing(float xy, float yx) noexcept;
	Matrix2& setTranspose/*transposeOf*/(const Matrix2& m) noexcept;
	Matrix2& setInverse/*inverseOf*/(const Matrix2& m) noexcept;
	Matrix2& setInverseTranspose/*inverseTransposeOf*/(const Matrix2& m) noexcept;
	Matrix2& preConcatenate(const Matrix2& m) noexcept;
	Matrix2& concatenate(const Matrix2& m) noexcept { return operator*=(m); }
	Matrix2& preScale(const Vector2<float>& v) noexcept;
	Matrix2& preScale(float f) noexcept { return operator*=(f); }
	Matrix2& scale(const Vector2<float>& v) noexcept;
	Matrix2& scale(float f) noexcept { return operator*=(f); }
	//Matrix2& rotate(float angle) noexcept { concatenate(Matrix2(Uninitialized()).setRotation(angle)); return *this; }
	//Matrix2& shear(float xy, float yx) noexcept { concatenate(Matrix2(Uninitialized()).setShearing(xy, yx)); return *this; }
	Matrix2& negate() noexcept;
	Matrix2& transpose() noexcept;
	Matrix2& invert() noexcept;
	Matrix2& orthonormalize() noexcept;

	//static const Matrix2& getZero() noexcept { return ZERO; }
	//static const Matrix2& getIdentity() noexcept { return IDENTITY; }

	static const Matrix2 ZERO;
	static const Matrix2 IDENTITY;

	union
	{
		simd::float4 row0;
		struct { float m00, m01; };
	};

	union
	{
		simd::float4 row1;
		struct { float m10, m11; };
	};
};

const Matrix2<float> Matrix2<float>::ZERO{};
const Matrix2<float> Matrix2<float>::IDENTITY{ 1.f, 0.f, 0.f, 1.f };

#endif /* SIMD_HAS_FLOAT4 */

template<typename T>
inline Matrix2<T>& Matrix2<T>::operator+=(const Matrix2<T>& m)
{
	m00 += m.m00; m01 += m.m01;
	m10 += m.m10; m11 += m.m11;
	return *this;
}

template<typename T>
inline Matrix2<T>& Matrix2<T>::operator-=(const Matrix2<T>& m)
{
	m00 -= m.m00; m01 -= m.m01;
	m10 -= m.m10; m11 -= m.m11;
	return *this;
}

template<typename T>
inline Matrix2<T>& Matrix2<T>::operator*=(T f)
{
	m00 *= f; m01 *= f;
	m10 *= f; m11 *= f;
	return *this;
}

template<typename T>
inline Matrix2<T>& Matrix2<T>::operator*=(const Matrix2<T>& m)
{
	set(m00*m.m00 + m01*m.m10, m00*m.m01 + m01*m.m11,
		m10*m.m00 + m11*m.m10, m10*m.m01 + m11*m.m11);
	return *this;
}

template<typename T>
	requires std::floating_point<T>
inline Matrix2<T> operator+(const Matrix2<T>& m1, const Matrix2<T>& m2) noexcept
{
	return Matrix2<T>(m1.m00 + m2.m00, m1.m01 + m2.m01, m1.m10 + m2.m10, m1.m11 + m2.m11);
}

template<typename T>
	requires std::floating_point<T>
inline Matrix2<T> operator-(const Matrix2<T>& m1, const Matrix2<T>& m2) noexcept
{
	return Matrix2<T>(m1.m00 - m2.m00, m1.m01 - m2.m01, m1.m10 - m2.m10, m1.m11 - m2.m11);
}

template<typename T>
	requires std::floating_point<T>
inline Matrix2<T> operator*(T f, const Matrix2<T>& m) noexcept
{
	return Matrix2<T>(f*m.m00, f*m.m01, f*m.m10, f*m.m11);
}

template<typename T>
	requires std::floating_point<T>
inline Matrix2<T> operator*(const Matrix2<T>& m, T f) noexcept
{
	return Matrix2<T>(m.m00*f, m.m01*f, m.m10*f, m.m11*f);
}

template<typename T>
	requires std::floating_point<T>
inline Matrix2<T> operator*(const Matrix2<T>& m1, const Matrix2<T>& m2) noexcept
{
	return Matrix2<T>(m1.m00*m2.m00 + m1.m01*m2.m10, m1.m00*m2.m01 + m1.m01*m2.m11,
		m1.m10*m2.m00 + m1.m11*m2.m10, m1.m10*m2.m01 + m1.m11*m2.m11);
}

template<typename T>
	requires std::floating_point<T>
inline Matrix2 operator/(const Matrix2& m, T f) noexcept
{ 
	return operator*(m, T(1)/f); 
}

template<typename C, typename T, typename U>
	requires std::floating_point<U>
inline std::basic_istream<C, T>& operator>>(std::basic_istream<C, T>& s, Matrix2<U>& m)
{
	return s >> m.m00 >> std::ws >> m.m01 >> std::ws >> m.m10 >> std::ws >> m.m11;
}

template<typename C, typename T, typename U>
	requires std::floating_point<U>
inline std::basic_ostream<C, T>& operator<<(std::basic_ostream<C, T>& s, const Matrix2<U>& m)
{
	constexpr C WS(0x20);
	return s << m.m00 << WS << m.m01 << WS << m.m10 << WS << m.m11;
}

template<typename T>
inline bool Matrix2<T>::isApproxZero() const
{
	return (std::fabs(m00) < Constants<T>::TOLERANCE) &&
		(std::fabs(m01) < Constants<T>::TOLERANCE) &&
		(std::fabs(m10) < Constants<T>::TOLERANCE) &&
		(std::fabs(m11) < Constants<T>::TOLERANCE);
}

template<typename T>
inline bool Matrix2<T>::isApproxIdentity() const
{
	return (std::fabs(T(1) - m00) < Constants<T>::TOLERANCE) &&
		(std::fabs(m01) < Constants<T>::TOLERANCE) &&
		(std::fabs(m10) < Constants<T>::TOLERANCE) &&
		(std::fabs(T(1) - m11) < Constants<T>::TOLERANCE);
}

template<typename T>
inline bool Matrix2<T>::isApproxEqual(const Matrix2<T>& m) const
{
	return (std::fabs(m.m00 - m00) < Constants<T>::TOLERANCE) &&
		(std::fabs(m.m01 - m01) < Constants<T>::TOLERANCE) &&
		(std::fabs(m.m10 - m10) < Constants<T>::TOLERANCE) &&
		(std::fabs(m.m11 - m11) < Constants<T>::TOLERANCE);
}

template<typename T>
inline bool Matrix2<T>::isApproxEqual(const Matrix2<T>& m, T tolerance) const
{
	return (std::fabs(m.m00 - m00) < tolerance) &&
		(std::fabs(m.m01 - m01) < tolerance) &&
		(std::fabs(m.m10 - m10) < tolerance) &&
		(std::fabs(m.m11 - m11) < tolerance);
}

template<typename T>
inline bool Matrix2<T>::isApproxOrthogonal() const
{
	Matrix2<T> m(*this);
	m *= Matrix2<T>(Uninitialized()).setTranspose(*this);
	return m.isApproxIdentity();
}

template<typename T>
inline bool Matrix2<T>::hasApproxUniformScaling() const
{
	Matrix2<T> m(*this);
	m *= Matrix2<T>(Uninitialized()).setTranspose(*this);
	return (std::fabs(m.m01) < Constants<T>::TOLERANCE) && (std::fabs(m.m10) < Constants<T>::TOLERANCE) && 
		(std::fabs(m.m11 - m.m00) < Constants<T>::TOLERANCE);
}

template<typename T>
inline Matrix2<T>& Matrix2<T>::set(const Vector2<T>& row0, const Vector2<T>& row1)
{ 
	m00 = row0.x; m01 = row0.y; 
	m10 = row1.x; m11 = row1.y; 
	return *this; 
}

template<typename T>
inline Matrix2<T>& Matrix2<T>::set(T m00, T m01, T m10, T m11)
{ 
	this->m00 = m00; this->m01 = m01; 
	this->m10 = m10; this->m11 = m11; 
	return *this; 
}

template<typename T>
inline Matrix2<T>& Matrix2<T>::setScaling(const Vector2<T>& v)
{
	m00 = v.x; m01 = T();
	m10 = T(); m11 = v.y;
	return *this;
}

template<typename T>
inline Matrix2<T>& Matrix2<T>::setScaling(T f)
{
	m00 = f; m01 = T();
	m10 = T(); m11 = f;
	return *this;
}

template<typename T>
inline Matrix2<T>& Matrix2<T>::setRotation(T angle)
{
	if (angle != T(0))
	{
		T sinAngle = std::sin(angle);
		T cosAngle = std::cos(angle);
		return set(cosAngle, sinAngle, -sinAngle, cosAngle);
	}

	return setIdentity();
}

template<typename T>
inline Matrix2<T>& Matrix2<T>::setShearing(T xy, T yx)
{
	m00 = T(1); m01 = xy;
	m10 = yx; m11 = T(1);
	return *this;
}

template<typename T>
inline Matrix2<T>& Matrix2<T>::setTranspose(const Matrix2<T>& m)
{
	m00 = m.m00; m01 = m.m10;
	m10 = m.m01; m11 = m.m11;
	return *this;
}

template<typename T>
inline Matrix2<T>& Matrix2<T>::setInverse(const Matrix2<T>& m)
{
	T a = m.m00, b = m.m01;
	T c = m.m10, d = m.m11;
	T n = T(1)/(a*d - b*c);
	m00 = d*n; m01 = -b*n;
	m10 = -c*n; m11 = a*n;
	return *this;
}

template<typename T>
inline Matrix2<T>& Matrix2<T>::setInverseTranspose(const Matrix2<T>& m)
{
	T a = m.m00, b = m.m01;
	T c = m.m10, d = m.m11;
	T n = T(1)/(a*d - b*c);
	m00 = d*n; m01 = -c*n;
	m10 = -b*n; m11 = a*n;
	return *this;
}

template<typename T>
inline Matrix2<T>& Matrix2<T>::preConcatenate(const Matrix2<T>& m)
{
	set(m.m00*m00 + m.m01*m10, m.m00*m01 + m.m01*m11,
		m.m10*m00 + m.m11*m10, m.m10*m01 + m.m11*m11);
	return *this;
}

template<typename T>
inline Matrix2<T>& Matrix2<T>::preScale(const Vector2<T>& v)
{
	m00 *= v.x; m01 *= v.x;
	m10 *= v.y; m11 *= v.y;
	return *this;
}

template<typename T>
inline Matrix2<T>& Matrix2<T>::scale(const Vector2<T>& v)
{
	m00 *= v.x; m01 *= v.y;
	m10 *= v.x; m11 *= v.y;
	return *this;
}

template<typename T>
inline Matrix2<T>& Matrix2<T>::negate()
{
	m00 = -m00; m01 = -m01;
	m10 = -m10; m11 = -m11;
	return *this;
}

template<typename T>
inline Matrix2<T>& Matrix2<T>::transpose()
{
	T t = m01; m01 = m10; m10 = t;
	return *this;
}

template<typename T>
inline Matrix2<T>& Matrix2<T>::invert()
{
	T a = m00, b = m01;
	T c = m10, d = m11;
	T n = T(1)/(a*d - b*c);
	m00 = d*n; m01 = -b*n;
	m10 = -c*n; m11 = a*n;
	return *this;
}

template<typename T>
inline Matrix2<T>& Matrix2<T>::orthonormalize()
{
	Vector2<T>& v1 = (*this)[0];
	Vector2<T>& v2 = (*this)[1];
	v2 -= project(v2, v1);
	v1.normalize();
	v2.normalize();
	return *this;
}

#if SIMD_HAS_FLOAT4

inline Matrix2<float>::Matrix2(float m00, float m01, float m10, float m11) : 
#if MATHEMATICS_SIMD_EXPAND_LAST
	row0(simd::set4(m00, m01, m01, m01)),
	row1(simd::set4(m10, m11, m11, m11))
#else
	row0(simd::set2(m00, m01)),
	row1(simd::set2(m10, m11))
#endif
{
}

#if MATHEMATICS_SIMD_EXPAND_LAST
inline Matrix2<float>::Matrix2(const float* m) :
	Matrix2(m[0], m[1], m[2], m[3])
{
}
#else
inline Matrix2<float>::Matrix2(const float* m)
{
	//simd::unpack2x2(m, row0, row1);
	std::tie(row0, row1) = simd::unpack2x2(m);
}
#endif

inline Matrix2<float>::Matrix2(simd::float4 m)
{ 
#if MATHEMATICS_SIMD_EXPAND_LAST
	auto [r0, r1] = simd::unpack2x2(m);
	row0 = simd::xyyy(r0);
	row1 = simd::xyyy(r1);
#else
	//simd::unpack2x2(m, row0, row1);
	std::tie(row0, row1) = simd::unpack2x2(m);
#endif
}

inline Matrix2<float> Matrix2<float>::operator-() const
{
#if MATHEMATICS_SIMD_EXPAND_LAST
	return Matrix2<float>(simd::neg4(row0), simd::neg4(row1));
#else
	return Matrix2<float>(simd::neg2(row0), simd::neg2(row1));
#endif
}

inline Matrix2<float>& Matrix2<float>::operator+=(const Matrix2<float>& m)
{
	row0 = simd::add4(row0, m.row0);
	row1 = simd::add4(row1, m.row1);
	return *this;
}

inline Matrix2<float>& Matrix2<float>::operator-=(const Matrix2<float>& m)
{
	row0 = simd::sub4(row0, m.row0);
	row1 = simd::sub4(row1, m.row1);
	return *this;
}

inline Matrix2<float>& Matrix2<float>::operator*=(float f)
{
	auto t = simd::set4(f);
	row0 = simd::mul4(row0, t);
	row1 = simd::mul4(row1, t);
	return *this;
}

inline Matrix2<float>& Matrix2<float>::operator*=(const Matrix2<float>& m)
{
	auto a = simd::pack2x2(row0, row1);
	auto b = simd::pack2x2(m.row0, m.row1);
	auto a1 = simd::swizzle<0, 0, 2, 2>(a);
	auto b1 = simd::swizzle<0, 1, 0, 1>(b);
	auto a2 = simd::swizzle<1, 1, 3, 3>(a);
	auto b2 = simd::swizzle<2, 3, 2, 3>(b);
	auto n = simd::add4(simd::mul4(a1, b1), simd::mul4(a2, b2));
#if MATHEMATICS_SIMD_EXPAND_LAST
	auto [r0, r1] = simd::unpack2x2(n);
	row0 = simd::xyyy(r0);
	row1 = simd::xyyy(r1);
#else
	//simd::unpack2x2(n, row0, row1);
	std::tie(row0, row1) = simd::unpack2x2(n);
#endif
	return *this;
}

template<>
inline Matrix2<float> operator+(const Matrix2<float>& m1, const Matrix2<float>& m2) noexcept
{
	return Matrix2<float>(simd::add4(m1.row0, m2.row0),
		simd::add4(m1.row1, m2.row1));
}

template<>
inline Matrix2<float> operator-(const Matrix2<float>& m1, const Matrix2<float>& m2) noexcept
{
	return Matrix2<float>(simd::sub4(m1.row0, m2.row0),
		simd::sub4(m1.row1, m2.row1));
}

template<>
inline Matrix2<float> operator*(float f, const Matrix2<float>& m) noexcept
{
	auto t = simd::set4(f);
	return Matrix2<float>(simd::mul4(t, m.row0),
		simd::mul4(t, m.row1));
}

template<>
inline Matrix2<float> operator*(const Matrix2<float>& m, float f) noexcept
{
	auto t = simd::set4(f);
	return Matrix2<float>(simd::mul4(m.row0, t),
		simd::mul4(m.row1, t));
}

template<>
inline Matrix2<float> operator*(const Matrix2<float>& m1, const Matrix2<float>& m2) noexcept
{
	auto a = simd::pack2x2(m1.row0, m1.row1);
	auto b = simd::pack2x2(m2.row0, m2.row1);
	auto a1 = simd::swizzle<0, 0, 2, 2>(a);
	auto b1 = simd::swizzle<0, 1, 0, 1>(b);
	auto a2 = simd::swizzle<1, 1, 3, 3>(a);
	auto b2 = simd::swizzle<2, 3, 2, 3>(b);
	auto n = simd::add4(simd::mul4(a1, b1), simd::mul4(a2, b2));
	Matrix2<float> m(Uninitialized());
#if MATHEMATICS_SIMD_EXPAND_LAST
	auto [r0, r1] = simd::unpack2x2(n);
	m.row0 = simd::xyyy(r0);
	m.row1 = simd::xyyy(r1);
#else
	//simd::unpack2x2(n, m.row0, m.row1);
	std::tie(m.row0, m.row1) = simd::unpack2x2(n);
#endif
	return m;
}

template<>
inline Matrix2<float> operator/(const Matrix2<float>& m, float f) noexcept 
{ 
	return operator*(m, 1.f/f); 
}

inline bool Matrix2<float>::operator==(const Matrix2<float>& m) const
{
	return simd::all4(simd::equal(simd::pack2x2(row0, row1), simd::pack2x2(m.row0, m.row1)));
}

template<typename C, typename T>
inline std::basic_istream<C, T>& operator>>(std::basic_istream<C, T>& s, Matrix2<float>& m)
{
	float m00, m01;
	float m10, m11;
	s >> m00 >> std::ws >> m01 >> std::ws >> m10 >> std::ws >> m11;
	m.set(m00, m01, m10, m11);
	return s;
}

template<typename A>
inline void Matrix2<float>::load(A& ar)
{
	float t00, t01;
	float t10, t11;
	ar(t00, t01, t10, t11);
	set(t00, t01, t10, t11);
}

inline bool Matrix2<float>::isZero() const
{ 
	return simd::all4(simd::equal(simd::pack2x2(row0, row1), simd::zero<simd::float4>())); 
}

inline bool Matrix2<float>::isApproxZero() const
{ 
	simd::all4(simd::lessThan(simd::abs4(simd::pack2x2(row0, row1)), Vector4<float>::TOLERANCE)); 
}

inline bool Matrix2<float>::isIdentity() const
{ 
	return simd::all4(simd::equal(simd::pack2x2(row0, row1), simd::constant4<simd::float4, 1, 0, 0, 1>())); 
}

inline bool Matrix2<float>::isApproxIdentity() const
{
	return simd::all4(simd::lessThan(simd::abs4(simd::sub4(simd::pack2x2(row0, row1), simd::constant4<simd::float4, 1, 0, 0, 1>())), 
		Vector4<float>::TOLERANCE));
}

inline bool Matrix2<float>::isApproxEqual(const Matrix2& m) const
{
	return simd::all4(simd::lessThan(simd::abs4(simd::sub4(simd::pack2x2(row0, row1), simd::pack2x2(m.row0, m.row1))), 
		Vector4<float>::TOLERANCE));
}

inline bool Matrix2<float>::isApproxEqual(const Matrix2& m, float tolerance) const
{
	return simd::all4(simd::lessThan(simd::abs4(simd::sub4(simd::pack2x2(row0, row1), simd::pack2x2(m.row0, m.row1))), 
		simd::set4(tolerance)));
}

inline bool Matrix2<float>::isApproxOrthogonal() const
{
	Matrix2<float> m(*this);
	m *= Matrix2<float>(Uninitialized()).setTranspose(*this);
	return m.isApproxIdentity();
}

inline bool Matrix2<float>::hasApproxUniformScaling() const
{
	Matrix2<float> m(*this);
	m *= Matrix2<float>(Uninitialized()).setTranspose(*this);
	return simd::all3(simd::lessThan(simd::abs4(simd::set3(m.m01, m.m10, m.m11 - m.m00)), Vector3<float>::TOLERANCE));
}

inline Matrix2<float>& Matrix2<float>::set(const Vector2<float>& row0, const Vector2<float>& row1)
{ 
	this->row0 = row0; 
	this->row1 = row1; 
	return *this; 
}

inline Matrix2<float>& Matrix2<float>::set(simd::float4 row0, simd::float4 row1)
{
	this->row0 = row0;
	this->row1 = row1;
	return *this;
}

inline Matrix2<float>& Matrix2<float>::set(float m00, float m01, float m10, float m11)
{ 
#if MATHEMATICS_SIMD_EXPAND_LAST
	row0 = simd::set4(m00, m01, m01, m01);
	row1 = simd::set4(m10, m11, m11, m11);
#else
	row0 = simd::set2(m00, m01);
	row1 = simd::set2(m10, m11);
#endif
	return *this;
}

inline Matrix2<float>& Matrix2<float>::setScaling(const Vector2<T>& v)
{
	row0 = simd::cutoff1(v);
#if MATHEMATICS_SIMD_EXPAND_LAST
	row1 = simd::and4(v, simd::constant4i<simd::float4, 0, -1, -1, -1>());
#else
	row1 = simd::and4(v, simd::constant4i<simd::float4, 0, -1, 0, 0>());
#endif
	return *this;
}

inline Matrix2<float>& Matrix2<float>::setScaling(float f)
{
	row0 = simd::set1(f);
#if MATHEMATICS_SIMD_EXPAND_LAST
	row1 = simd::set4(0.f, f, f, f);
#else
	row1 = simd::insert<1>(f, simd::zero<simd::float4>());
#endif
	return *this;
}

inline Matrix2<float>& Matrix2<float>::setRotation(float angle)
{
	if (angle != 0.f)
	{
		T sinAngle = std::sin(angle);
		T cosAngle = std::cos(angle);
		return set(cosAngle, sinAngle, -sinAngle, cosAngle);
	}

	return setIdentity();
}

inline Matrix2<float>& Matrix2<float>::setShearing(float xy, float yx)
{
#if MATHEMATICS_SIMD_EXPAND_LAST
	row0 = simd::set4(1.f, xy, xy, xy);
	row1 = simd::set4(yx, 1.f, 1.f, 1.f);
#else
	row0 = simd::set2(1.f, xy);
	row1 = simd::set2(yx, 1.f);
#endif
	return *this;
}

inline Matrix2<float>& Matrix2<float>::setTranspose(const Matrix2<float>& m)
{
#if MATHEMATICS_SIMD_EXPAND_LAST
	auto [r0, r1] = simd::transpose2x2(m.row0, m.row1);
	row0 = simd::xyyy(r0);
	row1 = simd::xyyy(r1);
#else
	//simd::transpose2x2(m.row0, m.row1, row0, row1);
	std::tie(row0, row1) = simd::transpose2x2(m.row0, m.row1);
#endif
	return *this;
}

inline Matrix2<float>& Matrix2<float>::setInverse(const Matrix2<float>& m)
{
	float a = m.m00, b = m.m01;
	float c = m.m10, d = m.m11;
	auto n = simd::set4(1.f/(a*d - b*c));
#if MATHEMATICS_SIMD_EXPAND_LAST
	b = -b;
	row0 = simd::mul4(simd::set4(d, b, b, b), n);
	row1 = simd::mul4(simd::set4(-c, a, a, a), n);
#else
	row0 = simd::mul4(simd::set2(d, -b), n);
	row1 = simd::mul4(simd::set2(-c, a), n);
#endif
	return *this;
}

inline Matrix2<float>& Matrix2<float>::setInverseTranspose(const Matrix2<float>& m)
{
	float a = m.m00, b = m.m01;
	float c = m.m10, d = m.m11;
	auto n = simd::set4(1.f/(a*d - b*c));
#if MATHEMATICS_SIMD_EXPAND_LAST
	c = -c;
	row0 = simd::mul4(simd::set4(d, c, c, c), n);
	row1 = simd::mul4(simd::set4(-b, a, a, a), n);
#else
	row0 = simd::mul4(simd::set2(d, -c), n);
	row1 = simd::mul4(simd::set2(-b, a), n);
#endif
	return *this;
}

inline Matrix2<float>& Matrix2<float>::preConcatenate(const Matrix2<float>& m)
{
	auto a = simd::pack2x2(m.row0, m.row1);
	auto b = simd::pack2x2(row0, row1);
	auto a1 = simd::swizzle<0, 0, 2, 2>(a);
	auto b1 = simd::swizzle<0, 1, 0, 1>(b);
	auto a2 = simd::swizzle<1, 1, 3, 3>(a);
	auto b2 = simd::swizzle<2, 3, 2, 3>(b);
	auto n = simd::add4(simd::mul4(a1, b1), simd::mul4(a2, b2));
#if MATHEMATICS_SIMD_EXPAND_LAST
	auto [r0, r1] = simd::unpack2x2(n);
	row0 = simd::xyyy(r0);
	row1 = simd::xyyy(r1);
#else
	//simd::unpack2x2(n, row0, row1);
	std::tie(row0, row1) = simd::unpack2x2(n);
#endif
	return *this;
}

inline Matrix2<float>& Matrix2<float>::preScale(const Vector2<float>& v)
{
	row0 = simd::mul4(row0, simd::xxxx(v));
	row1 = simd::mul4(row1, simd::yyyy(v));
	return *this;
}

inline Matrix2<float>& Matrix2<float>::scale(const Vector2<float>& v)
{
	row0 = simd::mul4(row0, v);
	row1 = simd::mul4(row1, v);
	return *this;
}

inline Matrix2<float>& Matrix2<float>::negate()
{
#if MATHEMATICS_SIMD_EXPAND_LAST
	row0 = simd::neg4(row0);
	row1 = simd::neg4(row1);
#else
	row0 = simd::neg2(row0);
	row1 = simd::neg2(row1);
#endif
	return *this;
}

inline Matrix2<float>& Matrix2<float>::transpose()
{
#if MATHEMATICS_SIMD_EXPAND_LAST
	auto [r0, r1] = simd::transpose2x2(row0, row1);
	row0 = simd::xyyy(r0);
	row1 = simd::xyyy(r1);
#else
	//simd::transpose2x2(row0, row1, row0, row1);
	std::tie(row0, row1) = simd::transpose2x2(row0, row1);
#endif
	return *this;
}

inline Matrix2<float>& Matrix2<float>::invert()
{
	float a = m00, b = m01;
	float c = m10, d = m11;
	auto n = simd::set4(1.f/(a*d - b*c));
#if MATHEMATICS_SIMD_EXPAND_LAST
	b = -b;
	row0 = simd::mul4(simd::set4(d, b, b, b), n);
	row1 = simd::mul4(simd::set4(-c, a, a, a), n);
#else
	row0 = simd::mul4(simd::set2(d, -b), n);
	row1 = simd::mul4(simd::set2(-c, a), n);
#endif
	return *this;
}

inline Matrix2<float>& Matrix2<float>::orthonormalize()
{
	Vector2<float>& v1 = (*this)[0];
	Vector2<float>& v2 = (*this)[1];
	v2 -= project(v2, v1);
	v1.normalize();
	v2.normalize();
	return *this;
}

#endif /* SIMD_HAS_FLOAT4 */

template<typename T>
	requires std::floating_point<T>
inline Matrix2<T> concatenate(const Matrix2<T>& m1, const Matrix2<T>& m2) noexcept
{
	return m1*m2;
}

template<typename T>
	requires std::floating_point<T>
inline Matrix2<T> concatenate(const Matrix2<T>& m1, const Matrix2<T>& m2, const Matrix2<T>& m3) noexcept
{
	return concatenate(concatenate(m1, m2), m3);
}

template<typename T>
	requires std::floating_point<T>
inline Matrix2<T> concatenate(const Matrix2<T>& m1, const Matrix2<T>& m2, const Matrix2<T>& m3, const Matrix2<T>& m4) noexcept
{
	return concatenate(concatenate(concatenate(m1, m2), m3), m4);
}

template<typename T>
	requires std::floating_point<T>
inline Matrix2<T> transpose(const Matrix2<T>& m) noexcept
{
	return Matrix2<T>(m.m00, m.m10, m.m01, m.m11);
}

template<typename T>
	requires std::floating_point<T>
inline Matrix2<T> transpose(Matrix2<T>&& m) noexcept
{
	m.transpose();
	return m;
}

template<typename T>
	requires std::floating_point<T>
inline Matrix2<T> inverse(const Matrix2<T>& m) noexcept
{
	return Matrix2<T>(Uninitialized()).setInverse(m);
}

template<typename T>
	requires std::floating_point<T>
inline Matrix2<T> inverse(Matrix2<T>&& m) noexcept
{
	m.invert();
	return m;
}

template<typename T>
	requires std::floating_point<T>
inline Matrix2<T> adjoint(const Matrix2<T>& m) noexcept
{
	return Matrix2<T>(m.m11, -m.m10, -m.m01, m.m00);
}

#if SIMD_HAS_FLOAT4

template<>
inline Matrix2<float> transpose(const Matrix2<float>& m) noexcept
{
	Matrix2<float> n(Uninitialized());
#if MATHEMATICS_SIMD_EXPAND_LAST
	auto [r0, r1] = simd::transpose2x2(m.row0, m.row1);
	n.row0 = simd::xyyy(r0);
	n.row1 = simd::xyyy(r1);
#else
	//simd::transpose2x2(m.row0, m.row1, n.row0, n.row1);
	std::tie(n.row0, n.row1) = simd::transpose2x2(m.row0, m.row1);
#endif
	return n;
}

template<>
inline Matrix2<float> transpose(Matrix2<float>&& m) noexcept
{
#if MATHEMATICS_SIMD_EXPAND_LAST
	auto [r0, r1] = simd::transpose2x2(m.row0, m.row1);
	m.row0 = simd::xyyy(r0);
	m.row1 = simd::xyyy(r1);
#else
	//simd::transpose2x2(m.row0, m.row1, m.row0, m.row1);
	std::tie(m.row0, m.row1) = simd::transpose2x2(m.row0, m.row1);
#endif
	return m;
}

#endif /* SIMD_HAS_FLOAT4 */

} // namespace templates

using Matrix2f = templates::Matrix2<float>;
using Matrix2d = templates::Matrix2<double>;

#if MATHEMATICS_DOUBLE
using Matrix2 = templates::Matrix2<double>;
using Matrix2Arg = templates::Matrix2<double>::ConstArg;
using Matrix2Result = templates::Matrix2<double>::ConstResult;
#else
using Matrix2 = templates::Matrix2<float>;
using Matrix2Arg = templates::Matrix2<float>::ConstArg;
using Matrix2Result = templates::Matrix2<float>::ConstResult;
#endif

} // namespace core::mathematics
