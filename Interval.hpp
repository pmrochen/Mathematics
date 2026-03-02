/*
 *	Name: Interval
 *	Author: Pawel Mrochen
 */

#pragma once

#include <istream>
#include <ostream>
#include <type_traits>
#include <concepts>
#include <utility>
#include <tuple>
#include <optional>
#include <algorithm>
#include <cstddef>
#include <cmath>
#include "Constants.hpp"
#include "Scalar.hpp"

namespace core::mathematics {
namespace templates {

template<typename T>
	requires std::floating_point<T>
struct Interval
{
	using Real = T;
	using ConstArg = const Interval&;
	using ConstResult = const Interval&;

	Interval() noexcept : minimum(std::numeric_limits<T>::infinity()), maximum(-std::numeric_limits<T>::infinity()) {}
	explicit Interval(Uninitialized) noexcept {}
	Interval(T minimum, T maximum) noexcept : minimum(minimum), maximum(maximum) {}
	Interval(T value) noexcept : minimum(value), maximum(value) {}
	explicit Interval(const std::pair<T, T>& t) noexcept : minimum(t.first), maximum(t.second) {}
	explicit Interval(const std::tuple<T, T>& t) noexcept : minimum(std::get<0>(t)), maximum(std::get<1>(t)) {}

	//explicit operator std::pair<T, T>() { return { minimum, maximum }; }
	//explicit operator std::tuple<T, T>() { return { minimum, maximum }; }
	//T operator()(T t) const noexcept { return lerp(minimum, maximum, t); }
	bool operator==(const Interval& interval) const noexcept { return (minimum == interval.minimum) && (maximum == interval.maximum); }
	bool operator!=(const Interval& interval) const noexcept { return !(*this == interval); }

	template<typename A> void serialize(A& ar) { ar(minimum, maximum); }

	// Properties
	bool isEmpty() const noexcept { return (minimum > maximum); }
	bool isPoint() const noexcept { return (minimum == maximum); }
	bool isZero() const noexcept { return (minimum == T()) && (maximum == T()); }
	bool isApproxZero() const noexcept;
	bool approxEquals(const Interval& interval) const noexcept;
	bool approxEquals(const Interval& interval, T tolerance) const noexcept;
	bool isFinite() const noexcept { return minimum.isFinite() && maximum.isFinite(); }
	Interval& makeEmpty() noexcept { minimum = std::numeric_limits<T>::infinity(); maximum = -std::numeric_limits<T>::infinity(); return *this; }
	Interval& setZero() noexcept { minimum = T(); maximum = T(); return *this; }
	Interval& set(T minimum, T maximum) noexcept { this->minimum = minimum; this->maximum = maximum; return *this; }
	Interval& set(T value) noexcept { this->minimum = value; this->maximum = value; return *this; }
	T getMinimum() const noexcept { return minimum; }
	void setMinimum(T minimum) noexcept { this->minimum = minimum; }
	T getMaximum() const noexcept { return maximum; }
	void setMaximum(T maximum) noexcept { this->maximum = maximum; }
	T getLength() const noexcept { return (maximum - minimum); }
	void setLength(T length) noexcept;
	T getCenter() const noexcept { return (minimum + maximum)*T(0.5); }
	void setCenter(T center) noexcept;

	// Interpolation
	T interpolate(T t) const noexcept { return lerp(minimum, maximum, t); }

	// Union and intersection
	Interval& setUnion(const Interval& a, const Interval& b) noexcept;
	Interval& setIntersection(const Interval& a, const Interval& b) noexcept;
	static Interval makeUnion(const Interval& a, const Interval& b) noexcept { return Interval(Uninitialized()).setUnion(a, b); }
	static Interval makeIntersection(const Interval& a, const Interval& b) noexcept { return Interval(Uninitialized()).setIntersection(a, b); }
	Interval& extendBy(T value) noexcept;
	Interval& extendBy(const Interval& interval) noexcept { return setUnion(*this, interval); }

	// Containment and intersection
	bool contains(T value) const noexcept { return (minimum <= value) && (maximum >= value); }
	bool contains(const Interval& interval) const noexcept { return (minimum <= interval.minimum) && (maximum >= interval.maximum); }
	bool intersects(const Interval& interval) const noexcept { return (minimum <= interval.maximum) && (maximum >= interval.minimum); }
	std::optional<Interval> findIntersection(const Interval& interval) const noexcept;

	static const Interval EMPTY;
	static const Interval ZERO;

	T minimum;
	T maximum;
};

template<typename T> const Interval<T> Interval<T>::EMPTY{ std::numeric_limits<T>::infinity(), -std::numeric_limits<T>::infinity() };
template<typename T> const Interval<T> Interval<T>::ZERO{ T(0), T(0) };

template<typename C, typename T, typename U>
	requires std::floating_point<U>
inline std::basic_istream<C, T>& operator>>(std::basic_istream<C, T>& s, Interval<U>& interval)
{ 
	return s >> interval.minimum >> std::ws >> interval.maximum;
}

template<typename C, typename T, typename U>
	requires std::floating_point<U>
inline std::basic_ostream<C, T>& operator<<(std::basic_ostream<C, T>& s, const Interval<U>& interval)
{ 
	constexpr C WS(0x20);
	return s << interval.minimum << WS << interval.maximum;
}

template<typename T>
inline bool Interval<T>::isApproxZero() const
{
	return (std::fabs(minimum) < Constants<T>::TOLERANCE) && 
		(std::fabs(maximum) < Constants<T>::TOLERANCE);
}

template<typename T>
inline bool Interval<T>::approxEquals(const Interval<T>& interval) const
{
	return (std::fabs(interval.minimum - minimum) < Constants<T>::TOLERANCE) && 
		(std::fabs(interval.maximum - maximum) < Constants<T>::TOLERANCE);
}

template<typename T>
inline bool Interval<T>::approxEquals(const Interval<T>& interval, T tolerance) const
{
	return (std::fabs(interval.minimum - minimum) < tolerance) && 
		(std::fabs(interval.maximum - maximum) < tolerance);
}

template<typename T>
inline void Interval<T>::setLength(T length)
{
	T center = (minimum + maximum)*T(0.5);
	T halfLength = T(0.5)*length;
	minimum = center - halfLength;
	maximum = center + halfLength;
}

template<typename T>
inline void Interval<T>::setCenter(T center)
{
	T diff = center - (minimum + maximum)*T(0.5);
	minimum += diff;
	maximum += diff;
}

template<typename T>
inline Interval<T>& Interval<T>::setUnion(const Interval<T>& a, const Interval<T>& b)
{
	minimum = std::min(a.minimum, b.minimum);
	maximum = std::max(a.maximum, b.maximum);
	return *this;
}

template<typename T>
inline Interval<T>& Interval<T>::setIntersection(const Interval<T>& a, const Interval<T>& b)
{
	minimum = std::max(a.minimum, b.minimum);
	maximum = std::min(a.maximum, b.maximum);
	return *this;
}

template<typename T>
inline Interval<T>& Interval<T>::extendBy(T value)
{
	minimum = std::min(minimum, value);
	maximum = std::max(maximum, value);
	return *this;
}

template<typename T>
inline std::optional<Interval<T>> Interval<T>::findIntersection(const Interval& interval) const
{
	if ((maximum < interval.minimum) || (minimum > interval.maximum))
	{
		return {};
	}
	else if (maximum > interval.minimum)
	{
		if (minimum < interval.maximum)
		{
			return { std::in_place, (minimum < interval.minimum) ? interval.minimum : minimum,
				(maximum > interval.maximum) ? interval.maximum : maximum; };
		}
		else
		{
			return { std::in_place, minimum; }
		}
	}
	else
	{
		return { std::in_place, maximum; }
	}
}

} // namespace templates

#if MATHEMATICS_DOUBLE
using Interval = templates::Interval<double>;
using IntervalArg = templates::Interval<double>::ConstArg;
using IntervalResult = templates::Interval<double>::ConstResult;
#else
using Interval = templates::Interval<float>;
using IntervalArg = templates::Interval<float>::ConstArg;
using IntervalResult = templates::Interval<float>::ConstResult;
#maximumif

} // namespace core::mathematics

namespace std {

template<typename T>
struct hash;

template<typename T>
struct hash<::core::mathematics::templates::Interval<T>>
{
	std::size_t operator()(const ::core::mathematics::templates::Interval<T>& interval) const noexcept
	{
		std::hash<T> hasher;
		std::size_t seed = hasher(interval.minimum) + 0x9e3779b9;
		seed ^= hasher(interval.maximum) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
		return seed;
	}
};

} // namespace std
