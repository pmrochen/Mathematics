/*
 *	Name: BezierCurve2
 *	Author: Pawel Mrochen
 */

#pragma once

#include <stdexcept>
#include <istream>
#include <ostream>
#include <type_traits>
#include <concepts>
#include <algorithm>
#include <functional>
#include <utility>
#include <tuple>
#include <optional>
#include <iterator>
#include <cstddef>
#include <cmath>
#include "Constants.hpp"
#include "Vector2.hpp"
#include "Matrix2.hpp"

namespace mathematics {
namespace templates {

template<typename T>
	requires std::floating_point<T>
struct BezierCurve2
{
	using Real = T;
	using ComponentType = T;
	using ConstArg = const BezierCurve2&;
	using ConstResult = const BezierCurve2&;
	using PointType = Vector2<T>;
	using TupleType = std::tuple<Vector2<T>, Vector2<T>, Vector2<T>, Vector2<T>>;

	static constexpr int NUM_COMPONENTS = 2;

	constexpr BezierCurve2() = default;
	explicit BezierCurve2(Uninitialized) noexcept;
	constexpr BezierCurve2(const Vector2<T>& p0, const Vector2<T>& p1, const Vector2<T>& p2, const Vector2<T>& p3) noexcept;
	explicit BezierCurve2(const TupleType& t) noexcept : points{ std::get<0>(t), std::get<1>(t), std::get<2>(t), std::get<3>(t) } {}
	explicit BezierCurve2(const Vector2<T>* p) noexcept : points{ p[0], p[1], p[2], p[3] } {}

	bool operator==(const BezierCurve2& curve) const noexcept;
	bool operator!=(const BezierCurve2& curve) const noexcept { return !(*this == curve); }

	template<typename A> void serialize(A& ar) { ar(points[0], points[1], points[2], points[3]); }

	// Properties
	bool isZero() const noexcept { return points[0].isZero() && points[1].isZero() && points[2].isZero() && points[3].isZero(); }
	bool isApproxZero() const noexcept;
	bool approxEquals(const BezierCurve2& curve) const noexcept;
	bool approxEquals(const BezierCurve2& curve, T tolerance) const noexcept;
	bool isFinite() const noexcept;
	BezierCurve2& setZero() noexcept;
	BezierCurve2& set(const Vector2<T>& p0, const Vector2<T>& p1, const Vector2<T>& p2, const Vector2<T>& p3) noexcept;

	// Control points
	const Vector2<T>& getControlPoint(int index) const noexcept { return ((unsigned int)index < 4u) ? points[index] : Vector2<T>::ZERO; }
	void setControlPoint(int index, const Vector2<T>& point); // throw (std::out_of_range)
	template<std::output_iterator<Vector2<T>> O> O copyControlPoints(O target) const;
	TupleType getControlPoints() const noexcept;

	// Transformation
	BezierCurve2& translate(const Vector2<T>& offset) noexcept;
	BezierCurve2& transform(const Matrix2<T>& matrix) noexcept;

	// Evaluation
	Vector2<T> evaluate(T t) const noexcept;
	template<std::size_t/*int*/ I> T evaluate(T t) const noexcept;
	Vector2<T> calculateDerivative(T t) const noexcept;
	Vector2<T> calculateTangent(T t) const noexcept { return calculateDerivative(t).normalize(); }
	T calculateSpeed(T t) const noexcept { return calculateDerivative(t).getMagnitude(); }
	template<int Order = 8> T calculateLength(T t = T(1)) const noexcept;
	template<int Order = 8, int Interations = 32> T calculateTime(T s, T tolerance = Constants<T>::TOLERANCE) const noexcept; // inverse mapping of s = Length(t) given by t = Length^{-1}(s)

	static const BezierCurve2 ZERO;

	Vector2<T> points[4];
};

template<typename T> const BezierCurve<T> BezierCurve<T>::ZERO{};

template<typename T>
inline BezierCurve2<T>::BezierCurve2(Uninitialized) : 
	points{ { Uninitialized() }, { Uninitialized() }, { Uninitialized() }, { Uninitialized() } } 
{
}

template<typename T>
inline BezierCurve2<T>::BezierCurve2(const Vector2<T>& p0, const Vector2<T>& p1, const Vector2<T>& p2, const Vector2<T>& p3) :
	points{ p0, p1, p2, p3 }
{ 
}

template<typename T>
inline bool BezierCurve2<T>::operator==(const BezierCurve2<T>& curve) const
{ 
	return (points[0] == curve.points[0]) && (points[1] == curve.points[1]) && (points[2] == curve.points[2]) &&
		(points[3] == curve.points[3]);
}

template<typename C, typename T, typename U>
	requires std::floating_point<U>
inline std::basic_istream<C, T>& operator>>(std::basic_istream<C, T>& s, BezierCurve2<U>& curve)
{ 
	return s >> curve.points[0] >> std::ws >> curve.points[1] >> std::ws >> curve.points[2] >> std::ws >> curve.points[3];
}

template<typename C, typename T, typename U>
	requires std::floating_point<U>
inline std::basic_ostream<C, T>& operator<<(std::basic_ostream<C, T>& s, const BezierCurve2<U>& curve)
{ 
	constexpr C WS(0x20);
	return s << curve.points[0] << WS << curve.points[1] << WS << curve.points[2] << WS << curve.points[3];
}

template<typename T>
inline bool BezierCurve2<T>::isApproxZero() const
{
	return points[0].isApproxZero() && points[1].isApproxZero() && points[2].isApproxZero() && points[3].isApproxZero();
}

template<typename T>
inline bool BezierCurve2<T>::approxEquals(const BezierCurve2<T>& curve) const
{
	return points[0].approxEquals(curve.points[0]) && points[1].approxEquals(curve.points[1]) && 
		points[2].approxEquals(curve.points[2]) && points[3].approxEquals(curve.points[3]);
}

template<typename T>
inline bool BezierCurve2<T>::approxEquals(const BezierCurve2<T>& curve, T tolerance) const
{
	return points[0].approxEquals(curve.points[0], tolerance) && points[1].approxEquals(curve.points[1], tolerance) &&
		points[2].approxEquals(curve.points[2], tolerance) && points[3].approxEquals(curve.points[3], tolerance);
}

template<typename T>
inline bool BezierCurve2<T>::isFinite() const
{ 
	return points[0].isFinite() && points[1].isFinite() && points[2].isFinite() && points[3].isFinite();
}

template<typename T>
inline BezierCurve2<T>& BezierCurve2<T>::setZero()
{
	points[0].setZero();
	points[1].setZero();
	points[2].setZero();
	points[3].setZero();
	return *this;
}

template<typename T>
inline BezierCurve2<T>& BezierCurve2<T>::set(const Vector2<T>& p0, const Vector2<T>& p1, const Vector2<T>& p2, const Vector2<T>& p3)
{
	points[0] = p0;
	points[1] = p1;
	points[2] = p2;
	points[3] = p3;
	return *this;
}

template<typename T>
inline void BezierCurve2<T>::setControlPoint(int index, const Vector2<T>& point)
{
	if ((unsigned int)index >= 4u)
		throw std::out_of_range("BezierCurve2::setControlPoint() : index");
	points[index] = point;
}

template<typename T>
template<std::output_iterator<Vector2<T>> O>
inline O BezierCurve2<T>::copyControlPoints(O target) const
{
	*target++ = points[0];
	*target++ = points[1];
	*target++ = points[2];
	*target++ = points[3];
	return target;
}

template<typename T>
inline typename BezierCurve2<T>::TupleType BezierCurve2<T>::getControlPoints() const
{ 
	return { points[0], points[1], points[2], points[3] }; 
}

template<typename T>
inline BezierCurve2<T>& BezierCurve2<T>::translate(const Vector2<T>& offset)
{
	points[0] += offset;
	points[1] += offset;
	points[2] += offset;
	points[3] += offset;
	return *this;
}

template<typename T>
inline BezierCurve2<T>& BezierCurve2<T>::transform(const Matrix2<T>& matrix)
{
	points[0] *= matrix;
	points[1] *= matrix;
	points[2] *= matrix;
	points[3] *= matrix;
	return *this;
}

} // namespace templates

#if MATHEMATICS_DOUBLE
using BezierCurve2 = templates::BezierCurve2<double>;
using BezierCurve2Arg = templates::BezierCurve2<double>::ConstArg;
using BezierCurve2Result = templates::BezierCurve2<double>::ConstResult;
#else
using BezierCurve2 = templates::BezierCurve2<float>;
using BezierCurve2Arg = templates::BezierCurve2<float>::ConstArg;
using BezierCurve2Result = templates::BezierCurve2<float>::ConstResult;
#endif

} // namespace mathematics

namespace std {

template<typename T>
struct hash;

template<typename T>
struct hash<::mathematics::templates::BezierCurve2<T>>
{
	size_t operator()(const ::mathematics::templates::BezierCurve2<T>& curve) const noexcept
	{
		hash<typename ::mathematics::templates::Vector2<T>> hasher;
		size_t seed = hasher(curve.points[0]) + 0x9e3779b9;
		seed ^= hasher(curve.points[1]) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
		seed ^= hasher(curve.points[2]) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
		seed ^= hasher(curve.points[3]) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
		return seed;
	}
};

} // namespace std

#include "Bezier.inl"
#include "Romberg.inl"

namespace mathematics::templates {

template<typename T>
inline Vector2<T> Bezier2<T>::evaluate(T t) const // #TODO SIMD
{
	Vector4<T> bt = bezier::basis(t);
	return Vector2<T>(points[0].x*bt[0] + points[1].x*bt[1] + points[2].x*bt[2] + points[3].x*bt[3],
		points[0].y*bt[0] + points[1].y*bt[1] + points[2].y*bt[2] + points[3].y*bt[3]);
}

template<typename T>
template<std::size_t/*int*/ I>
inline T Bezier2<T>::evaluate(T t) const // #TODO SIMD
{
	Vector4<T> bt = bezier::basis(t);
	return points[0].get<I>()*bt[0] + points[1].get<I>()*bt[1] + points[2].get<I>()*bt[2] + points[3].get<I>()*bt[3];
}

template<typename T>
inline Vector2<T> Bezier2<T>::calculateDerivative(T t) const // #TODO SIMD
{
	Vector4 dbt = bezier::derivativeBasis(t);
	return Vector2<T>(points[0].x*dbt[0] + points[1].x*dbt[1] + points[2].x*dbt[2] + points[3].x*dbt[3],
		points[0].y*dbt[0] + points[1].y*dbt[1] + points[2].y*dbt[2] + points[3].y*dbt[3]);
}

template<typename T>
template<int Order>
inline T Bezier2<T>::calculateLength(T t) const
{
	if (t <= T(0)) 
		return T(0);
	
	return romberg::estimate<Order>(T(0), std::min(t, T(1)), [this](T t) { return calculateSpeed(t); });
}

template<typename T>
template<int Order, int Iterations>
inline std::optional<T> Bezier2<T>::calculateTime(T s, T tolerance) const
{
	if (s <= T(0)) 
		return { T(0) };

	T totalLen = calculateLength(T(1));
	if (s >= totalLen) 
		return { T(1) };

	T time = s/totalLen;
    for (int i = 0; i < Iterations; i++)
    {
        T difference = calculateLength<Order>(time) - s;
        if (std::fabs(difference) < tolerance)
			return { time };

        time -= difference/calculateSpeed(time);
    }
    
	return {}; // Newton's method failed. If this happens, increase iterations or tolerance or integration accuracy.
}

} // namespace mathematics::templates
