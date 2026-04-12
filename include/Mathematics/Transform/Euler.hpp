/*
 *	Name: Euler
 *	Author: Pawel Mrochen
 */

#pragma once

#include <stdexcept>
#include <istream>
#include <ostream>
#include <limits>
#include <type_traits>
#include <concepts>
#include <algorithm>
#include <functional>
#include <utility>
#include <tuple>
#include <string>
#include <cstddef>
#include <cctype>
#include <cmath>
#include "Constants.hpp"
#include "Vector3.hpp"
#include "Matrix3.hpp"
#include "Quaternion.hpp"
#include "YawPitchRoll.hpp"

namespace mathematics::templates {

template<typename T>
	requires std::floating_point<T>
struct Euler;

} // namespace mathematics::templates

namespace std {

template<size_t I, typename T>
struct tuple_element;

template<size_t I, typename T>
struct tuple_element<I, ::mathematics::templates::Euler<T>>
{
	using type = T;
};

template<typename T>
struct tuple_element<3, ::mathematics::templates::Euler<T>>
{
	using type = EulerOrder;
};

template<typename T>
struct tuple_size;

template<typename T>
struct tuple_size<::mathematics::templates::Euler<T>> : integral_constant<size_t, 4> 
{
};

} // namespace std

namespace mathematics {
namespace templates {

template<typename T>
	requires std::floating_point<T>
struct Matrix3;

template<typename T>
	requires std::floating_point<T>
struct Euler
{
	using Real = T;

	constexpr Euler() noexcept : x(), y(), z(), order() {}
	explicit Euler(Uninitialized) noexcept {}
	constexpr explicit Euler(EulerOrder order) noexcept : x(), y(), z(), order(order) {}
	constexpr Euler(T x, T y, T z, EulerOrder order) noexcept : x(x), y(y), z(z), order(order) {}
	constexpr Euler(const Euler& e, EulerOrder order) noexcept : x(e.x), y(e.y), z(e.z), order(order) {}
	constexpr Euler(const YawPitchRoll<T>& r) noexcept : x(r.pitch), y(r.yaw), z(r.roll), order(EulerOrder::ZXY) {}
	Euler(const Quaternion<T>& q, EulerOrder order) noexcept;
	Euler(const Matrix3<T>& m, EulerOrder order) noexcept;
	explicit Euler(const std::tuple<T, T, T, EulerOrder>& t) noexcept : x(std::get<0>(t)), y(std::get<1>(t)), z(std::get<2>(t)), order(std::get<3>(t)) {}
	template<Arithmetic U> explicit Euler(const std::tuple<U, U, U, EulerOrder>& t) noexcept : x(T(std::get<0>(t))), y(T(std::get<1>(t))), z(T(std::get<2>(t))), order(std::get<3>(t)) {}
	Euler(const T* e, EulerOrder order) noexcept : x(e[0]), y(e[1]), z(e[2]), order(order) {}

	//explicit operator std::tuple<T, T, T, EulerOrder>() { return std::tuple<T, T, T, EulerOrder>(x, y, z, order); }
	//template<Arithmetic U> explicit operator std::tuple<U, U, U, EulerOrder>() { return std::tuple<U, U, U, EulerOrder>(U(x), U(y), U(z), order); }
	explicit operator T*() noexcept { return &x; }
	explicit operator const T*() const noexcept { return &x; }
	T& operator[](int i) noexcept { return (&x)[i]; }
	const T& operator[](int i) const noexcept { return (&x)[i]; }

	Euler operator+() const noexcept { return *this; }
	Euler operator-() const noexcept { return Euler(-x, -y, -z, order); }
	Euler& operator+=(const Euler& e); // throw (std::invalid_argument);
	Euler& operator-=(const Euler& e); // throw (std::invalid_argument);
	Euler& operator*=(T f) noexcept { x *= f; y *= f; z *= f; return *this; }
	Euler& operator/=(T f) noexcept { return operator*=(T(1)/f); }
	bool operator==(const Euler& e) const noexcept { return (order == e.order) && (x == e.x) && (y == e.y) && (z == e.z); }
	bool operator!=(const Euler& e) const noexcept { return !(*this == e); }

	template<typename A> void serialize(A& ar) { ar(x, y, z, order); }

	template<std::size_t I> typename std::tuple_element<I, Euler>::type& get() noexcept;
	template<std::size_t I> const typename std::tuple_element<I, Euler>::type& get() const noexcept;

	bool isZero() const noexcept { return (x == T()) && (y == T()) && (z == T()); }
	bool isApproxZero() const noexcept;
	bool approxEquals(const Euler& e) const noexcept;
	bool approxEquals(const Euler& e, T tolerance) const noexcept;
	bool isFinite() const noexcept { return std::isfinite(x) && std::isfinite(y) && std::isfinite(z); }
	Euler& setZero() noexcept { x = T(); y = T(); z = T(); return *this; }
	Euler& setZero(EulerOrder order) noexcept { x = T(); y = T(); z = T(); this->order = order; return *this; }
	Euler& set(T x, T y, T z) noexcept { this->x = x; this->y = y; this->z = z; return *this; }
	Euler& set(T x, T y, T z, EulerOrder order) noexcept { this->x = x; this->y = y; this->z = z; this->order = order; return *this; }
	Euler& negate() noexcept { x = -x; y = -y; z = -z; return *this; }
	//Euler& clamp(T low, T high);

	static const Euler ZERO;
	static const Euler ZERO_XYZ;
	static const Euler ZERO_XZY;
	static const Euler ZERO_YZX;
	static const Euler ZERO_YXZ;
	static const Euler ZERO_ZXY;
	static const Euler ZERO_ZYX;

    T x, y, z;
	EulerOrder order;
};

template<typename T> const Euler<T> Euler<T>::ZERO{};
template<typename T> const Euler<T> Euler<T>::ZERO_XYZ{ T(), T(), T(), EulerOrder::XYZ };
template<typename T> const Euler<T> Euler<T>::ZERO_XZY{ T(), T(), T(), EulerOrder::XZY };
template<typename T> const Euler<T> Euler<T>::ZERO_YZX{ T(), T(), T(), EulerOrder::YZX };
template<typename T> const Euler<T> Euler<T>::ZERO_YXZ{ T(), T(), T(), EulerOrder::YXZ };
template<typename T> const Euler<T> Euler<T>::ZERO_ZXY{ T(), T(), T(), EulerOrder::ZXY };
template<typename T> const Euler<T> Euler<T>::ZERO_ZYX{ T(), T(), T(), EulerOrder::ZYX };

template<typename T>
inline Euler<T>::Euler(const Quaternion<T>& q, EulerOrder order) : 
	Euler(Matrix3<T>::makeRotation(q), order) // #FIXME Don't call converting constructor when q is identity
{
	if (q.isIdentity())
		setZero(order);
}

template<typename T>
inline Euler<T>::Euler(const Matrix3<T>& m, EulerOrder order)
{
	if (order != EulerOrder::UNSPECIFIED)
	{
		// Source: http://www.realtimerendering.com/resources/GraphicsGems/gemsiv/euler_angle/

		static const int safe[] = { 0, 1, 2, 0 };
		static const int next[] = { 1, 2, 0, 1 };
		unsigned int o = (unsigned int)order;
		int f = o & 1; o >>= 1;
		int s = o & 1; o >>= 1;
		int n = o & 1; o >>= 1;
		int i = safe[o & 3];
		int j = next[i + n];
		int k = next[i + 1 - n];

		T ti, tj, th;
		if (s)
		{
			T sy = std::sqrt(m[j][i]*m[j][i] + m[k][i]*m[k][i]);
			if (sy > T(16)*Constants<T>::EPSILON/*Constants<T>::TOLERANCE*/)
			{
				ti = std::atan2(m[j][i], m[k][i]);
				tj = std::atan2(sy, m[i][i]);
				th = std::atan2(m[i][j], -m[i][k]);
			}
			else
			{
				ti = std::atan2(-m[k][j], m[j][j]);
				tj = std::atan2(sy, m[i][i]);
				th = T(0);
			}
		}
		else
		{
			T cy = std::sqrt(m[i][i]*m[i][i] + m[i][j]*m[i][j]);
			if (cy > T(16)*Constants<T>::EPSILON/*Constants<T>::TOLERANCE*/)
			{
				ti = std::atan2(m[j][k], m[k][k]);
				tj = std::atan2(-m[i][k], cy);
				th = std::atan2(m[i][j], m[i][i]);
			}
			else
			{
				ti = std::atan2(-m[k][j], m[j][j]);
				tj = std::atan2(-m[i][k], cy);
				th = T(0);
			}
		}

		if (n)
		{
			ti = -ti;
			tj = -tj;
			th = -th;
		}

		if (f)
		{
			T t = ti;
			ti = th;
			th = t;
		}

		T e[3];
		e[f ? (s ? i : k) : i] = ti;
		e[j] = tj;
		e[f ? i : (s ? i : k)] = th;
		set(e[0], e[1], e[2], order);
	}
	else
		setZero(order);
}

template<typename T>
inline Euler<T>& Euler<T>::operator+=(const Euler<T>& e)
{
	if (order != e.order) 
		throw std::invalid_argument("Euler<T>::operator+=(const Euler<T>&) : order");
	x += e.x; 
    y += e.y; 
    z += e.z;
	return *this;
}

template<typename T>
inline Euler<T>& Euler<T>::operator-=(const Euler<T>& e)
{
	if (order != e.order) 
		throw std::invalid_argument("Euler<T>::operator-=(const Euler<T>&) : order");
	x -= e.x; 
    y -= e.y; 
    z -= e.z;
	return *this;
}

template<typename T>
	requires std::floating_point<T>
inline Euler<T> operator+(const Euler<T>& e1, const Euler<T>& e2) // throw (std::invalid_argument);
{
	if (e1.order != e2.order)
		throw std::invalid_argument("operator+(const Euler<T>&, const Euler<T>&) : order");
	return Euler<T>(e1.x + e2.x, e1.y + e2.y, e1.z + e2.z, e1.order);
}

template<typename T>
	requires std::floating_point<T>
inline Euler<T> operator-(const Euler<T>& e1, const Euler<T>& e2) // throw (std::invalid_argument);
{
	if (e1.order != e2.order) 
		throw std::invalid_argument("operator-(const Euler<T>&, const Euler<T>&) : order");
	return Euler<T>(e1.x - e2.x, e1.y - e2.y, e1.z - e2.z, e1.order);
}

template<typename T>
	requires std::floating_point<T>
inline Euler<T> operator*(T f, const Euler<T>& e) noexcept
{ 
	return Euler(f*e.x, f*e.y, f*e.z, e.order); 
}

template<typename T>
	requires std::floating_point<T>
inline Euler<T> operator*(const Euler<T>& e, T f) noexcept
{ 
	return Euler(e.x*f, e.y*f, e.z*f, e.order); 
}

template<typename T>
	requires std::floating_point<T>
inline Euler<T> operator/(T f, const Euler<T>& e) noexcept
{ 
	return Euler(f/e.x, f/e.y, f/e.z, e.order); 
}

template<typename T>
	requires std::floating_point<T>
inline Euler<T> operator/(const Euler<T>& e, T f) noexcept
{ 
	return operator*(e, T(1)/f); 
}

template<typename C, typename T, typename U>
	requires std::floating_point<U>
inline std::basic_istream<C, T>& operator>>(std::basic_istream<C, T>& s, Euler<U>& e)
{
	std::basic_string<C, T> order;
	s >> e.x >> std::ws >> e.y >> std::ws >> e.z >> std::ws >> order;
	if (order == "XYZ")
		e.order = EulerOrder::XYZ;
	else if (order == "XZY")
		e.order = EulerOrder::XZY;
	else if (order == "YZX")
		e.order = EulerOrder::YZX;
	else if (order == "YXZ")
		e.order = EulerOrder::YXZ;
	else if (order == "ZXY")
		e.order = EulerOrder::ZXY;
	else if (order == "ZYX")
		e.order = EulerOrder::ZYX;
	else if (std::all_of(order.begin(), order.end(), [](const auto& c) { return std::isdigit(c) || (c == C('-')); }))
		e.order = (EulerOrder)std::stoi(order);
    else
        e.order = EulerOrder::UNSPECIFIED;
	return s;
}

template<typename C, typename T, typename U>
	requires std::floating_point<U>
inline std::basic_ostream<C, T>& operator<<(std::basic_ostream<C, T>& s, const Euler<U>& e)
{
	const C* order = "Unspecified";
	switch (e.order)
	{
		case EulerOrder::XYZ: 
			order = "XYZ"; 
			break;
		case EulerOrder::XZY: 
			order = "XZY"; 
			break;
		case EulerOrder::YZX: 
			order = "YZX"; 
			break;
		case EulerOrder::YXZ: 
			order = "YXZ"; 
			break;
		case EulerOrder::ZXY: 
			order = "ZXY"; 
			break;
		case EulerOrder::ZYX: 
			order = "ZYX"; 
			break;
	}
	constexpr C WS(0x20);
	return s << e.x << WS << e.y << WS << e.z << WS << order;
}

template<typename T>
template<std::size_t I>
inline typename std::tuple_element<I, Euler<T>>::type& Euler<T>::get()
{
	if constexpr (I == 0)
		return x;
	else if constexpr (I == 1)
		return y;
	else if constexpr (I == 2)
		return z;
	else if constexpr (I == 3)
		return order;
	static_assert(false);
}

template<typename T>
template<std::size_t I>
inline const typename std::tuple_element<I, Euler<T>>::type& Euler<T>::get() const
{
	if constexpr (I == 0)
		return x;
	else if constexpr (I == 1)
		return y;
	else if constexpr (I == 2)
		return z;
	else if constexpr (I == 3)
		return order;
	static_assert(false);
}

template<typename T>
inline bool Euler<T>::isApproxZero() const
{ 
	return (std::fabs(x) < Constants<T>::TOLERANCE) && (std::fabs(y) < Constants<T>::TOLERANCE) && 
		(std::fabs(z) < Constants<T>::TOLERANCE);
}

template<typename T>
inline bool Euler<T>::approxEquals(const Euler<T>& e) const
{ 
	return (order == e.order) && (std::fabs(e.x - x) < Constants<T>::TOLERANCE) && 
        (std::fabs(e.y - y) < Constants<T>::TOLERANCE) && (std::fabs(e.z - z) < Constants<T>::TOLERANCE);
}

template<typename T>
inline bool Euler<T>::approxEquals(const Euler<T>& e, T tolerance) const
{ 
	return (order == e.order) && (std::fabs(e.x - x) < tolerance) && 
        (std::fabs(e.y - y) < tolerance) && (std::fabs(e.z - z) < tolerance);
}

//template<typename T>
//inline Euler<T>& Euler<T>::clamp(T low, T high)
//{
//	x = std::clamp(x, low, high);
//	y = std::clamp(y, low, high);
//	z = std::clamp(z, low, high);
//	return *this;
//}

template<std::size_t I, typename T>
	requires std::floating_point<T>
inline typename std::tuple_element<I, Euler<T>>::type& get(Euler<T>& e) noexcept
{
	if constexpr (I == 0)
		return e.x;
	else if constexpr (I == 1)
		return e.y;
	else if constexpr (I == 2)
		return e.z;
	else if constexpr (I == 3)
		return e.order;
	static_assert(false);
}

template<std::size_t I, typename T>
	requires std::floating_point<T>
inline const typename std::tuple_element<I, Euler<T>>::type& get(const Euler<T>& e) noexcept
{
	if constexpr (I == 0)
		return e.x;
	else if constexpr (I == 1)
		return e.y;
	else if constexpr (I == 2)
		return e.z;
	else if constexpr (I == 3)
		return e.order;
	static_assert(false);
}

template<std::size_t I, typename T>
	requires std::floating_point<T>
inline typename std::tuple_element<I, Euler<T>>::type&& get(Euler<T>&& e) noexcept
{
	if constexpr (I == 0)
		return e.x;
	else if constexpr (I == 1)
		return e.y;
	else if constexpr (I == 2)
		return e.z;
	else if constexpr (I == 3)
		return e.order;
	static_assert(false);
}

template<std::size_t I, typename T>
	requires std::floating_point<T>
inline const typename std::tuple_element<I, Euler<T>>::type&& get(const Euler<T>&& e) noexcept
{
	if constexpr (I == 0)
		return e.x;
	else if constexpr (I == 1)
		return e.y;
	else if constexpr (I == 2)
		return e.z;
	else if constexpr (I == 3)
		return e.order;
	static_assert(false);
}

} // namespace templates

#if MATHEMATICS_DOUBLE
using Euler = templates::Euler<double>;
#else
using Euler = templates::Euler<float>;
#endif

} // namespace mathematics

namespace std {

template<typename T>
struct hash;

template<typename T>
struct hash<::mathematics::templates::Euler<T>>
{
	size_t operator()(const ::mathematics::templates::Euler<T>& e) const noexcept
	{
		hash<T> hasher;
		size_t seed = hasher(e.x) + 0x9e3779b9;
		seed ^= hasher(e.y) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
		seed ^= hasher(e.z) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
		seed ^= std::hash<int/*typename ::mathematics::EulerOrder*/>()((int)e.order) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
		return seed;
	}
};

} // namespace std
