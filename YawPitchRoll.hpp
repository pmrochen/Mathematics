/*
 *	Name: YawPitchRoll
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
#include "Constants.hpp"
#include "Vector3.hpp"
#include "Matrix3.hpp"
#include "Quaternion.hpp"

namespace mathematics {
namespace templates {

template<typename T>
	requires std::floating_point<T>
struct Euler;

template<typename T>
	requires std::floating_point<T>
struct YawPitchRoll
{
	using Real = T;
	using TupleType = std::tuple<T, T, T>;

	constexpr YawPitchRoll() noexcept : yaw(), pitch(), roll() {}
	explicit YawPitchRoll(Uninitialized) noexcept {}
	constexpr YawPitchRoll(T yaw, T pitch, T roll) noexcept : yaw(yaw), pitch(pitch), roll(roll) {}
	explicit YawPitchRoll(const Euler<T>& e) noexcept;
	explicit YawPitchRoll(const Quaternion<T>& q) noexcept;
	explicit YawPitchRoll(const Matrix3<T>& m) noexcept;
	//explicit YawPitchRoll(const Vector3<T>& direction) noexcept; // -> fromForward
	//YawPitchRoll(const Vector3<T>& direction, T roll) noexcept;
	explicit YawPitchRoll(const TupleType& t) noexcept : yaw(std::get<0>(t)), pitch(std::get<1>(t)), roll(std::get<2>(t)) {}
	//YawPitchRoll(const T* r) noexcept : yaw(r[0]), pitch(r[1]), roll(r[2]) {}

	//explicit operator TupleType() { return { yaw, pitch, roll }; }
	//explicit operator T*() noexcept { return &yaw; }
	//explicit operator const T*() const noexcept { return &yaw; }

	YawPitchRoll operator+() const noexcept { return *this; }
	YawPitchRoll operator-() const noexcept { return YawPitchRoll(-yaw, -pitch, -roll); }
	YawPitchRoll& operator+=(const YawPitchRoll& r) noexcept { yaw += r.yaw; pitch += r.pitch; roll += r.roll; return *this; }
	YawPitchRoll& operator-=(const YawPitchRoll& r) noexcept { yaw -= r.yaw; pitch -= r.pitch; roll -= r.roll; return *this; }
	YawPitchRoll& operator*=(T f) noexcept { yaw *= f; pitch *= f; roll *= f; return *this; }
	YawPitchRoll& operator/=(T f) noexcept { return operator*=(T(1)/f); }
	bool operator==(const YawPitchRoll& r) const noexcept { return (yaw == r.yaw) && (pitch == r.pitch) && (roll == r.roll); }
	bool operator!=(const YawPitchRoll& r) const noexcept { return !(*this == r); }

	template<typename A> void serialize(A& ar) { ar(yaw, pitch, roll); }

	template<std::size_t I> T& get() noexcept;
	template<std::size_t I> const T& get() const noexcept;

    static YawPitchRoll fromForward(const Vector3<T>& forward, T roll = T()) noexcept;
	Vector3 toForward/*toDirection*/() const noexcept;

    bool isZero() const noexcept { return (yaw == T()) && (pitch == T()) && (roll == T()); }
	bool isApproxZero() const noexcept;
	bool approxEquals(const YawPitchRoll& r) const noexcept;
	bool approxEquals(const YawPitchRoll& r, T tolerance) const noexcept;
	bool isFinite() const noexcept { return std::isfinite(yaw) && std::isfinite(pitch) && std::isfinite(roll); }
	YawPitchRoll& setZero() noexcept { yaw = T(); pitch = T(); roll = T(); return *this; }
	YawPitchRoll& set(T yaw, T pitch, T roll) noexcept { this->yaw = yaw; this->pitch = pitch; this->roll = roll; return *this; }
	YawPitchRoll& negate() noexcept { yaw = -yaw; pitch = -pitch; roll = -roll; return *this; }
	//YawPitchRoll& clamp(T low, T high);

	static const YawPitchRoll ZERO;

    T yaw, pitch, roll;
};

template<typename T> const YawPitchRoll<T> YawPitchRoll<T>::ZERO{};

template<typename T>
inline YawPitchRoll<T>::YawPitchRoll(const Quaternion<T>& q) : 
	YawPitchRoll(Matrix3<T>::makeRotation(q)) // #FIXME Don't call converting constructor when q is identity
{
	if (q.isIdentity())
		setZero();
}

template<typename T>
	requires std::floating_point<T>
inline YawPitchRoll<T> operator+(const YawPitchRoll<T>& r1, const YawPitchRoll<T>& r2) noexcept
{
	return YawPitchRoll<T>(r1.yaw + r2.yaw, r1.pitch + r2.pitch, r1.roll + r2.roll);
}

template<typename T>
	requires std::floating_point<T>
inline YawPitchRoll<T> operator-(const YawPitchRoll<T>& r1, const YawPitchRoll<T>& r2) noexcept
{
	return YawPitchRoll<T>(r1.yaw - r2.yaw, r1.pitch - r2.pitch, r1.roll - r2.roll);
}

template<typename T>
	requires std::floating_point<T>
inline YawPitchRoll<T> operator*(T f, const YawPitchRoll<T>& r) noexcept
{ 
	return YawPitchRoll(f*r.yaw, f*r.pitch, f*r.roll); 
}

template<typename T>
	requires std::floating_point<T>
inline YawPitchRoll<T> operator*(const YawPitchRoll<T>& r, T f) noexcept
{ 
	return YawPitchRoll(r.yaw*f, r.pitch*f, r.roll*f); 
}

template<typename T>
	requires std::floating_point<T>
inline YawPitchRoll<T> operator/(T f, const YawPitchRoll<T>& r) noexcept
{ 
	return YawPitchRoll(f/r.yaw, f/r.pitch, f/r.roll); 
}

template<typename T>
	requires std::floating_point<T>
inline YawPitchRoll<T> operator/(const YawPitchRoll<T>& r, T f) noexcept
{ 
	return operator*(r, T(1)/f); 
}

template<typename C, typename T, typename U>
	requires std::floating_point<U>
inline std::basic_istream<C, T>& operator>>(std::basic_istream<C, T>& s, YawPitchRoll<U>& r)
{ 
	return s >> r.yaw >> std::ws >> r.pitch >> std::ws >> r.roll; 
}

template<typename C, typename T, typename U>
	requires std::floating_point<U>
inline std::basic_ostream<C, T>& operator<<(std::basic_ostream<C, T>& s, const YawPitchRoll<U>& r)
{ 
	constexpr C WS(0x20);
	return s << r.yaw << WS << r.pitch << WS << r.roll;
}

template<typename T>
template<std::size_t I>
inline T& YawPitchRoll<T>::get()
{
	if constexpr (I == 0)
		return yaw;
	else if constexpr (I == 1)
		return pitch;
	else if constexpr (I == 2)
		return roll;
	static_assert(false);
}

template<typename T>
template<std::size_t I>
inline const T& YawPitchRoll<T>::get() const
{
	if constexpr (I == 0)
		return yaw;
	else if constexpr (I == 1)
		return pitch;
	else if constexpr (I == 2)
		return roll;
	static_assert(false);
}

template<typename T>
inline YawPitchRoll<T> YawPitchRoll<T>::fromForward(const Vector3<T>& forward, T roll)
{
    T h = T(0);
    T p = T(0);
    T m = direction.getMagnitude();
    if (m > T(0))
    {
        Vector3<T> direction = forward/m;
        p = std::asin(std::clamp(-direction.y, T(-1), T(1)));
        if ((T(1) - std::fabs(direction.y)) >= Constants<T>::TOLERANCE)
        {
            T cosP = std::cos(p);
            h = std::acos(std::clamp(direction.z/cosP, T(-1), T(1)));
            if (direction.x < T(0)) 
                h = -h; // h = Constants<T>::TWO_PI - h;
        }
	}

    return YawPitchRoll<T>(h, p, roll);
}

template<typename T>
inline Vector3<T> YawPitchRoll<T>::toForward() const
{
	T sh = std::sin(yaw);
    T ch = std::cos(yaw);
    T sp = std::sin(pitch);
    T cp = std::cos(pitch);
	return Vector3<T>(sh*cp, -sp, ch*cp);
}

template<typename T>
inline bool YawPitchRoll<T>::isApproxZero() const
{ 
	return (std::fabs(yaw) < Constants<T>::TOLERANCE) && (std::fabs(pitch) < Constants<T>::TOLERANCE) && 
		(std::fabs(roll) < Constants<T>::TOLERANCE);
}

template<typename T>
inline bool YawPitchRoll<T>::approxEquals(const YawPitchRoll<T>& r) const
{ 
	return (std::fabs(r.yaw - yaw) < Constants<T>::TOLERANCE) && (std::fabs(r.pitch - pitch) < Constants<T>::TOLERANCE) &&
		(std::fabs(r.roll - roll) < Constants<T>::TOLERANCE);
}

template<typename T>
inline bool YawPitchRoll<T>::approxEquals(const YawPitchRoll<T>& r, T tolerance) const
{ 
	return (std::fabs(r.yaw - yaw) < tolerance) && (std::fabs(r.pitch - pitch) < tolerance) &&
		(std::fabs(r.roll - roll) < tolerance);
}

//template<typename T>
//inline YawPitchRoll<T>& YawPitchRoll<T>::clamp(T low, T high)
//{
//	yaw = std::clamp(yaw, low, high);
//	pitch = std::clamp(pitch, low, high);
//	roll = std::clamp(roll, low, high);
//	return *this;
//}

template<std::size_t I, typename T>
	requires std::floating_point<T>
inline T& get(YawPitchRoll<T>& r) noexcept
{
	if constexpr (I == 0)
		return r.yaw;
	else if constexpr (I == 1)
		return r.pitch;
	else if constexpr (I == 2)
		return r.roll;
	static_assert(false);
}

template<std::size_t I, typename T>
	requires std::floating_point<T>
inline const T& get(const YawPitchRoll<T>& r) noexcept
{
	if constexpr (I == 0)
		return r.yaw;
	else if constexpr (I == 1)
		return r.pitch;
	else if constexpr (I == 2)
		return r.roll;
	static_assert(false);
}

template<std::size_t I, typename T>
	requires std::floating_point<T>
inline T&& get(YawPitchRoll<T>&& r) noexcept
{
	if constexpr (I == 0)
		return r.yaw;
	else if constexpr (I == 1)
		return r.pitch;
	else if constexpr (I == 2)
		return r.roll;
	static_assert(false);
}

template<std::size_t I, typename T>
	requires std::floating_point<T>
inline const T&& get(const YawPitchRoll<T>&& r) noexcept
{
	if constexpr (I == 0)
		return r.yaw;
	else if constexpr (I == 1)
		return r.pitch;
	else if constexpr (I == 2)
		return r.roll;
	static_assert(false);
}

} // namespace templates

#if MATHEMATICS_DOUBLE
using YawPitchRoll = templates::YawPitchRoll<double>;
#else
using YawPitchRoll = templates::YawPitchRoll<float>;
#endif

} // namespace mathematics

namespace std {

template<size_t I, typename T>
struct tuple_element;

template<size_t I, typename T>
struct tuple_element<I, ::mathematics::templates::YawPitchRoll<T>>
{
	using type = T;
};

template<typename T>
struct tuple_size;

template<typename T>
struct tuple_size<::mathematics::templates::YawPitchRoll<T>> : integral_constant<size_t, 3> 
{
};

template<typename T>
struct hash;

template<typename T>
struct hash<::mathematics::templates::YawPitchRoll<T>>
{
	size_t operator()(const ::mathematics::templates::YawPitchRoll<T>& r) const noexcept
	{
		hash<T> hasher;
		size_t seed = hasher(r.yaw) + 0x9e3779b9;
		seed ^= hasher(r.pitch) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
		seed ^= hasher(r.roll) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
		return seed;
	}
};

} // namespace std

#include "Euler.hpp"

namespace mathematics::templates {

template<typename T>
inline YawPitchRoll<T>::YawPitchRoll(const Euler<T>& e)
{
	if (e.order == EulerOrder::ZXY)
	{
        yaw = e.y;
		pitch = e.x;
		roll = e.z;
	}
	else if (e.order != EulerOrder::UNSPECIFIED)
	{
		Euler<T> euler(Matrix3<T>(e), EulerOrder::ZXY);
		yaw = euler.y;
		pitch = euler.x;
		roll = euler.z;
	}
	else
	{
		yaw = T();
		pitch = T();
		roll = T();
	}
}

template<typename T>
inline YawPitchRoll<T>::YawPitchRoll(const Matrix3<T>& m) : YawPitchRoll<T>(Euler<T>(m, EulerOrder::ZXY))
{
}

} // namespace mathematics::templates
