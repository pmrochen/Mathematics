/*
 *	Name: Intersections
 *	Author: Pawel Mrochen
 */

#pragma once

#include <limits>
#include <type_traits>
#include <concepts>
#include <optional>
#include <algorithm>
#include <cmath>
#include "Simd/Intrinsics.hpp"
#include "Constants.hpp"
#include "Scalar.hpp"
#include "Vector2.hpp"
#include "Vector3.hpp"
#include "Vector4.hpp"
#include "Matrix3.hpp"
#include "Distances.inl"

namespace mathematics::intersections {
namespace detail {

//template <template <typename...> class, template<typename...> class>
//struct IsSameTemplate : public std::false_type {};
//
//template <template <typename...> class T>
//struct IsSameTemplate<T, T> : public std::true_type {};

template <typename T>
struct IsOptional : public std::false_type {};

template <typename T>
struct IsOptional<std::optional<T>> : public std::true_type {};

template <typename T>
	requires (std::floating_point<T> /*|| std::integral<T>*/)
struct FloatIntTypeConverter;

template<>
struct FloatIntTypeConverter<float>
{
	using Type = int;
};

template<>
struct FloatIntTypeConverter<double>
{
	using Type = long long;
};

//template<>
//struct FloatIntTypeConverter<int>
//{
//	using Type = float;
//};
//
//template<>
//struct FloatIntTypeConverter<long long>
//{
//	using Type = double;
//};

//template<typename T>
//	requires std::floating_point<T>
//struct IntersectionLookup
//{
//	static int get(int i) noexcept
//	{
//		static const int lut[] = { 0, 1, 2, -1, 0, 1, 5, -1, 0, 4, 2, -1, 0, 4, 5, -1, 3, 1, 2, -1, 3, 1, 5, -1, 3, 4, 2, -1, 3, 4, 5, -1 };
//		return lut[i];
//	};
//};
//
//#if SIMD_HAS_FLOAT4
//
//template<>
//struct IntersectionLookup<float>
//{
//	static int get(int i) noexcept
//	{
//		static const int lut[] = { 0, 1, 2, -1, 0, 1, 6, -1, 0, 5, 2, -1, 0, 5, 6, -1, 4, 1, 2, -1, 4, 1, 6, -1, 4, 5, 2, -1, 4, 5, 6, -1 };
//		return lut[i];
//	};
//};
//
//#endif /* SIMD_HAS_FLOAT4 */

template<typename O, typename T>
	requires std::floating_point<T>
inline O infinity() noexcept
{
	if constexpr (IsOptional<O>::value)
		return {};
	else
		return { std::numeric_limits<T>::infinity() };
}

template<typename O, typename T>
	requires std::floating_point<T>
inline O interval(T value) noexcept
{
	if constexpr (IsOptional<O>::value)
		return { std::in_place, value };
	else
		return { value };
}

template<typename O, typename T>
	requires std::floating_point<T>
inline O interval(T minimum, T maximum) noexcept
{
	if constexpr (IsOptional<O>::value)
		return { std::in_place, minimum, maximum };
	else
		return { minimum, maximum };
}

template<typename O, typename T>
	requires std::floating_point<T>
inline O intervalChecked(T minimum, T maximum) noexcept
{
	if constexpr (detail::IsOptional<O>::value)
	{
		if (minimum > maximum)
			return {};
		else
			return { std::in_place, minimum, maximum };
	}
	else
		return { minimum, maximum };
}

template<typename T>
	requires std::floating_point<T>
inline bool testAxisX01(T a, T b, T fa, T fb, const Vector3<T>& halfDims, const Vector3<T>& v0, const Vector3<T>& v1,
	const Vector3<T>& v2) noexcept
{
	T p0 = a*v0.y - b*v0.z;
	T p2 = a*v2.y - b*v2.z;
	T minimum, maximum;
	if (p0 < p2)
	{
		minimum = p0;
		maximum = p2;
	}
	else
	{
		minimum = p2;
		maximum = p0;
	}

	T rad = fa*halfDims.y + fb*halfDims.z;
	return !((minimum > rad) || (maximum < -rad));
}

template<typename T>
	requires std::floating_point<T>
inline bool testAxisX2(T a, T b, T fa, T fb, const Vector3<T>& halfDims, const Vector3<T>& v0, const Vector3<T>& v1,
	const Vector3<T>& v2) noexcept
{
	T p0 = a*v0.y - b*v0.z;
	T p1 = a*v1.y - b*v1.z;
	T minimum, maximum;
	if (p0 < p1)
	{
		minimum = p0;
		maximum = p1;
	}
	else
	{
		minimum = p1;
		maximum = p0;
	}

	T rad = fa*halfDims.y + fb*halfDims.z;
	return !((minimum > rad) || (maximum < -rad));
}

template<typename T>
	requires std::floating_point<T>
inline bool testAxisY02(T a, T b, T fa, T fb, const Vector3<T>& halfDims, const Vector3<T>& v0, const Vector3<T>& v1,
	const Vector3<T>& v2) noexcept
{
	T p0 = -a*v0.x + b*v0.z;
	T p2 = -a*v2.x + b*v2.z;
	T minimum, maximum;
	if (p0 < p2)
	{
		minimum = p0;
		maximum = p2;
	}
	else
	{
		minimum = p2;
		maximum = p0;
	}

	T rad = fa*halfDims.x + fb*halfDims.z;
	return !((minimum > rad) || (maximum < -rad));
}

template<typename T>
	requires std::floating_point<T>
inline bool testAxisY1(T a, T b, T fa, T fb, const Vector3<T>& halfDims, const Vector3<T>& v0, const Vector3<T>& v1,
	const Vector3<T>& v2) noexcept
{
	T p0 = -a*v0.x + b*v0.z;
	T p1 = -a*v1.x + b*v1.z;
	T minimum, maximum;
	if (p0 < p1)
	{
		minimum = p0;
		maximum = p1;
	}
	else
	{
		minimum = p1;
		maximum = p0;
	}

	T rad = fa*halfDims.x + fb*halfDims.z;
	return !((minimum > rad) || (maximum < -rad));
}

template<typename T>
	requires std::floating_point<T>
inline bool testAxisZ12(T a, T b, T fa, T fb, const Vector3<T>& halfDims, const Vector3<T>& v0, const Vector3<T>& v1,
	const Vector3<T>& v2) noexcept
{
	T p1 = a*v1.x - b*v1.y;
	T p2 = a*v2.x - b*v2.y;
	T minimum, maximum;
	if (p2 < p1)
	{
		minimum = p2;
		maximum = p1;
	}
	else
	{
		minimum = p1;
		maximum = p2;
	}

	T rad = fa*halfDims.x + fb*halfDims.y;
	return !((minimum > rad) || (maximum < -rad));
}

template<typename T>
	requires std::floating_point<T>
inline bool testAxisZ0(T a, T b, T fa, T fb, const Vector3<T>& halfDims, const Vector3<T>& v0, const Vector3<T>& v1,
	const Vector3<T>& v2) noexcept
{
	T p0 = a*v0.x - b*v0.y;
	T p1 = a*v1.x - b*v1.y;
	T minimum, maximum;
	if (p0 < p1)
	{
		minimum = p0;
		maximum = p1;
	}
	else
	{
		minimum = p1;
		maximum = p0;
	}

	T rad = fa*halfDims.x + fb*halfDims.y;
	return !((minimum > rad) || (maximum < -rad));
}

} // namespace detail

template<typename T>
	requires std::floating_point<T>
inline bool testLineLine(const Vector2<T>& originA, const Vector2<T>& directionA, const Vector2<T>& originB,
	const Vector2<T>& directionB) noexcept
{
	return !(std::fabs(cross(directionA, directionB)) < Constants<T>::TOLERANCE) ||
		(std::fabs(cross(normalize(originB - originA), directionA)) < Constants<T>::TOLERANCE);
}

template<typename O, typename T>
	requires std::floating_point<T>
inline O findLineLine(const Vector2<T>& originA, const Vector2<T>& directionA, const Vector2<T>& originB, 
	const Vector2<T>& directionB) noexcept
{
	// http://www.geometrictools.com/

	Vector2<T> diff = originB - originA;
	T d1CrossD2 = cross(directionA, directionB);
	if (std::fabs(d1CrossD2) < Constants<T>::TOLERANCE)
	{
		if (!(std::fabs(cross(normalize(diff), directionA)) < Constants<T>::TOLERANCE))
			return detail::infinity<O, T>();
		return { T(0) }; // collinear
	}

	return { cross(diff, directionB)/d1CrossD2 };
}

template<typename O, typename T>
	requires std::floating_point<T>
inline O findLineRay(const Vector2<T>& originA, const Vector2<T>& directionA, const Vector2<T>& originB,
	const Vector2<T>& directionB) noexcept
{
	// http://www.geometrictools.com/

	Vector2<T> diff = originA - originB;
	T d2CrossD1 = cross(directionB, directionA);
	if (std::fabs(d2CrossD1) < Constants<T>::TOLERANCE)
	{
		if (!(std::fabs(cross(normalize(diff), directionB)) < Constants<T>::TOLERANCE))
			return detail::infinity<O, T>();
		return { T(0) }; // collinear
	}

	T t = cross(diff, directionA)/d2CrossD1;
	if (t >= T(0))
		return { t };
	return detail::infinity<O, T>();
}

template<typename O, typename T>
	requires std::floating_point<T>
inline O findLineSegment(const Vector2<T>& origin, const Vector2<T>& direction, const Vector2<T>& start, const Vector2<T>& end) noexcept
{
	// http://www.geometrictools.com/

	Vector2<T> directionB = end - start;
	Vector2<T> diff = origin - start;
	T d2CrossD1 = cross(directionB, direction);
	if (std::fabs(d2CrossD1) < Constants<T>::TOLERANCE)
	{
		if (!(std::fabs(cross(normalize(diff), directionB)) < Constants<T>::TOLERANCE))
			return detail::infinity<O, T>();
		return { T(0) }; // collinear
	}

	T t = cross(diff, direction)/d2CrossD1;
	if ((t >= T(0)) && (t <= T(1))
		return { t };
	return detail::infinity<O, T>();
}

template<typename O, typename T>
	requires std::floating_point<T>
inline O findSegmentSegment(const Vector2<T>& startA, const Vector2<T>& endA, const Vector2<T>& startB, const Vector2<T>& endB) noexcept
{
	// http://www.geometrictools.com/

	Vector2<T> directionA = endA - startA;
	Vector2<T> directionB = endB - startB;
	Vector2<T> diff = startB - startA;
	T d1CrossD2 = cross(directionA, directionB);
	if (std::fabs(d1CrossD2) < Constants<T>::TOLERANCE)
	{
		if (!(std::fabs(cross(normalize(diff), directionA)) < Constants<T>::TOLERANCE))
			return detail::infinity<O, T>();

		T d1DotD1 = dot(directionA, directionA);
		T d1DotD2 = dot(directionA, directionB);
		T t0 = dot(diff, directionA)/d1DotD1;
		T t1 = t0 + d1DotD2/d1DotD1;
		if (d1DotD2 < T(0))
			std::swap(t0, t1);
		if ((t0 <= T(1)) && (t1 >= T(0)))
			return std::optional<T>(std::max(t0, T(0)));
		return detail::infinity<O, T>();
	}

	T t = cross(diff, directionB)/d1CrossD2;
	T u = cross(diff, directionA)/d1CrossD2;
	if ((t >= T(0)) && (t <= T(1)) && (u >= T(0)) && (u <= T(1)))
		return std::optional<T>(t);
	return detail::infinity<O, T>();
}

template<typename T>
	requires std::floating_point<T>
inline bool testLinePlane(const Vector3<T>& direction, const Vector3<T>& normal) noexcept
{
	return (std::fabs(dot(normal, direction)) >= Constants<T>::TOLERANCE);
}

template<typename O, typename T>
	requires std::floating_point<T>
inline O findLinePlane(const Vector3<T>& origin, const Vector3<T>& direction, const Vector3<T>& normal, T constant) noexcept
{
	T nDotD = dot(normal, direction);
	if (std::fabs(nDotD) < Constants<T>::TOLERANCE)
		return detail::infinity<O, T>();
	return { (-plane.d - dot(normal, origin))/nDotD };
}

template<typename O, typename T>
	requires std::floating_point<T>
inline O findLineTriangle(const Vector3<T>& origin, const Vector3<T>& direction, const Vector3<T>& vertex0, const Vector3<T>& vertex1,
	const Vector3<T>& vertex2) noexcept
{
	// http://fileadmin.cs.lth.se/cs/Personal/Tomas_Akenine-Moller/code/

	Vector3<T> edge1 = vertex1 - vertex0;
	Vector3<T> edge2 = vertex2 - vertex0;
	Vector3<T> pVec = cross(direction, edge2);
	T det = dot(edge1, pVec);
	if (std::fabs(det) < Constants<T>::TOLERANCE)
		return detail::infinity<O, T>();

	T invDet = T(1)/det;
	Vector3<T> tVec = origin - vertex0;
	T u = dot(tVec, pVec)*invDet;
	if ((u < T(0)) || (u > T(1)))
		return detail::infinity<O, T>();

	Vector3<T> qVec = cross(tVec, edge1);
	T v = dot(direction, qVec)*invDet;
	if ((v < T(0)) || ((u + v) > T(1)))
		return detail::infinity<O, T>();

	return { dot(edge2, qVec)*invDet };
}

//template<typename O, typename T>
//	requires std::floating_point<T>
//inline O findRayTriangle(const Vector3<T>& origin, const Vector3<T>& direction, const Vector3<T>& vertex0, const Vector3<T>& vertex1,
//	const Vector3<T>& vertex2) noexcept
//{
//	// http://fileadmin.cs.lth.se/cs/Personal/Tomas_Akenine-Moller/code/
//
//	Vector3<T> edge1 = vertex1 - vertex0;
//	Vector3<T> edge2 = vertex2 - vertex0;
//	Vector3<T> pVec = cross(ray.direction, edge2);
//	T det = dot(edge1, pVec);
//	if (std::fabs(det) < Constants<T>::TOLERANCE)
//		return detail::infinity<O, T>();
//
//	T invDet = T(1)/det;
//	Vector3<T> tVec = origin - vertex0;
//	T u = dot(tVec, pVec)*invDet;
//	if ((u < T(0)) || (u > T(1)))
//		return detail::infinity<O, T>();
//
//	Vector3<T> qVec = cross(tVec, edge1);
//	T v = dot(direction, qVec)*invDet;
//	if ((v < T(0)) || ((u + v) > T(1)))
//		return detail::infinity<O, T>();
//
//	T t = dot(edge2, qVec)*invDet;
//	if (t >= T(0))
//		return { t };
//	else
//		return detail::infinity<O, T>();
//}
//
//template<typename O, typename T>
//	requires std::floating_point<T>
//inline O findSegmentTriangle(const Vector3<T>& start, const Vector3<T>& end, const Vector3<T>& vertex0, const Vector3<T>& vertex1,
//	const Vector3<T>& vertex2) noexcept
//{
//	// http://fileadmin.cs.lth.se/cs/Personal/Tomas_Akenine-Moller/code/
//
//	Vector3<T> edge1 = vertex1 - vertex0;
//	Vector3<T> edge2 = vertex2 - vertex0;
//	Vector3<T> direction = segment.end - segment.start;
//	Vector3<T> pVec = cross(direction, edge2);
//	T det = dot(edge1, pVec);
//	if (std::fabs(det) < Constants<T>::TOLERANCE)
//		return detail::infinity<O, T>();
//
//	T invDet = T(1)/det;
//	Vector3<T> tVec = start - vertex0;
//	T u = dot(tVec, pVec)*invDet;
//	if ((u < T(0)) || (u > T(1)))
//		return detail::infinity<O, T>();
//
//	Vector3<T> qVec = cross(tVec, edge1);
//	T v = dot(direction, qVec)*invDet;
//	if ((v < T(0)) || ((u + v) > T(1)))
//		return detail::infinity<O, T>();
//
//	T t = dot(edge2, qVec)*invDet;
//	if ((t >= T(0)) && (t <= T(1)))
//		return { t };
//	else
//		return detail::infinity<O, T>();
//}

template<typename O, typename T>
	requires std::floating_point<T>
inline O findLineAxisAlignedRectangle(const Vector2<T>& origin, const Vector2<T>& direction, const Vector2<T>& minimum,
	const Vector2<T>& maximum) noexcept
{
	Vector2<T> invDir = T(1)/direction;

	T tMin, tMax;
	if (invDir/*direction*/.x >= T(0))
	{
		tMin = (minimum.x - origin.x)*invDir.x;
		tMax = (maximum.x - origin.x)*invDir.x;
	}
	else
	{
		tMin = (maximum.x - origin.x)*invDir.x;
		tMax = (minimum.x - origin.x)*invDir.x;
	}

	T tyMin, tyMax;
	if (invDir/*direction*/.y >= T(0))
	{
		tyMin = (minimum.y - origin.y)*invDir.y;
		tyMax = (maximum.y - origin.y)*invDir.y;
	}
	else
	{
		tyMin = (maximum.y - origin.y)*invDir.y;
		tyMax = (minimum.y - origin.y)*invDir.y;
	}

	if ((tMin > tyMax) || (tyMin > tMax))
		return {};

	if (tyMin > tMin)
		tMin = tyMin;
	if (tyMax < tMax)
		tMax = tyMax;

	return detail::intervalChecked<O>(tMin, tMax);
}

template<typename O, typename T>
	requires std::floating_point<T>
inline O findLineAxisAlignedBox(const Vector3<T>& origin, const Vector3<T>& direction, const Vector3<T>& minimum,
	const Vector3<T>& maximum) noexcept
{
	Vector3<T> invDir = T(1)/direction;

	T tMin, tMax;
	if (invDir/*direction*/.x >= T(0))
	{
		tMin = (minimum.x - origin.x)*invDir.x;
		tMax = (maximum.x - origin.x)*invDir.x;
	}
	else
	{
		tMin = (maximum.x - origin.x)*invDir.x;
		tMax = (minimum.x - origin.x)*invDir.x;
	}

	T tyMin, tyMax;
	if (invDir/*direction*/.y >= T(0))
	{
		tyMin = (minimum.y - origin.y)*invDir.y;
		tyMax = (maximum.y - origin.y)*invDir.y;
	}
	else
	{
		tyMin = (maximum.y - origin.y)*invDir.y;
		tyMax = (minimum.y - origin.y)*invDir.y;
	}

	if ((tMin > tyMax) || (tyMin > tMax))
		return {};

	if (tyMin > tMin)
		tMin = tyMin;
	if (tyMax < tMax)
		tMax = tyMax;

	T tzMin, tzMax;
	if (invDir/*direction*/.z >= T(0))
	{
		tzMin = (minimum.z - origin.z)*invDir.z;
		tzMax = (maximum.z - origin.z)*invDir.z;
	}
	else
	{
		tzMin = (maximum.z - origin.z)*invDir.z;
		tzMax = (minimum.z - origin.z)*invDir.z;
	}

	if ((tMin > tzMax) || (tzMin > tMax))
		return {};

	if (tzMin > tMin)
		tMin = tzMin;
	if (tzMax < tMax)
		tMax = tzMax;

	return detail::intervalChecked<O>(tMin, tMax);
}

#if SIMD_HAS_FLOAT4

//template<typename O>
//inline O findLineAxisAlignedRectangle(const Vector2<float>& origin, const Vector2<float>& direction, 
//	const Vector2<float>& minimum, const Vector2<float>& maximum) noexcept // #TODO
//{
//}

template<typename O>
inline O findLineAxisAlignedBox(const Vector3<float>& origin, const Vector3<float>& direction,
	const Vector3<float>& minimum, const Vector3<float>& maximum) noexcept
{
	auto invDir = simd::div4(Vector4<float>::ONE, direction);
	auto l1 = simd::mul4(simd::sub4(minimum, origin), invDir);
	auto l2 = simd::mul4(simd::sub4(maximum, origin), invDir);

	// The order we use for those min/max is vital to filter out NaNs that happens
	// when an invDir is +/- inf and (minimum - origin) is 0. inf*0 = NaN
	auto filteredL1a = simd::min4(l1, Vector4<float>::INF);
	auto filteredL2a = simd::min4(l2, Vector4<float>::INF);
	auto filteredL1b = simd::max4(l1, Vector4<float>::MINUS_INF);
	auto filteredL2b = simd::max4(l2, Vector4<float>::MINUS_INF);
	auto lMax = simd::max4(filteredL1a, filteredL2a);
	auto lMin = simd::min4(filteredL1b, filteredL2b);
	auto lMax0 = simd::swizzle<1, 2, 3, 0>(lMax);
	auto lMin0 = simd::swizzle<1, 2, 3, 0>(lMin);
	lMax = simd::min4/*min1*/(lMax, lMax0);
	lMin = simd::max4/*max1*/(lMin, lMin0);
	auto lMax1 = simd::swizzle<2, 3, 2, 3>(lMax);
	auto lMin1 = simd::swizzle<2, 3, 2, 3>(lMin);
	lMax = simd::min4/*min1*/(lMax, lMax1);
	lMin = simd::max4/*max1*/(lMin, lMin1);

	return detail::intervalChecked<O>(simd::toFloat(lMin), simd::toFloat(lMax));
}

#endif /* SIMD_HAS_FLOAT4 */

template<typename O, typename T>
	requires std::floating_point<T>
inline O findLineOrientedBox(const Vector3<T>& origin, const Vector3<T>& direction, const Vector3<T>& center,
	const Matrix3<T>& basis, const Vector3<T>& halfDims) noexcept
{
	//Matrix3<T> basisTranspose(transpose(basis));
	return findLineAxisAlignedBox<O, T>(basis*(origin - center)/*(origin - center)*basisTranspose*/,
		basis*direction/*direction*basisTranspose*/, -box.halfDims, box.halfDims);
}

template<typename T, typename V>
	requires std::floating_point<T> // #TODO Add vector (V) constraints
inline bool testLineNSphere(const V& origin, const V& direction, const V& center, T radius) noexcept
{
	V diff = origin - center;
	T a = dot(direction, direction);
	T b = T(2)*dot(direction, diff);
	T c = dot(diff, diff) - radius*radius;
	return !((b*b - T(4)*a*c) < T(0));
}

template<typename T, typename V>
	requires std::floating_point<T> // #TODO Add vector (V) constraints
inline bool testNormalizedLineNSphere(const V& origin, const V& direction, const V& center, T radius) noexcept
{
	V diff = origin - center;
	T halfB = dot(direction, diff);
	T c = dot(diff, diff) - radius*radius;
	return !((halfB*halfB - c) < T(0));
}

template<typename O, typename T, typename V>
	requires std::floating_point<T> // #TODO Add vector (V) constraints
inline O findLineNSphere(const V& origin, const V& direction, const V& center, T radius) noexcept
{
	V diff = origin - sphere.center;
	T a = dot(direction, direction);
	T b = T(2)*dot(direction, diff);
	T c = dot(diff, diff) - sphere.radius*sphere.radius;
	T delta = b*b - T(4)*a*c;

	if (delta < T(0))
	{
		return {};
	}
	else if (delta > T(0))
	{
		delta = std::sqrt(delta);
		a = T(0.5)/a;
		return detail::interval<O>((-b - delta)*a, (-b + delta)*a);
	}
	else // delta == 0
	{
		return detail::interval<O>(-b*T(0.5)/a);
	}
}

template<typename O, typename T, typename V>
	requires std::floating_point<T> // #TODO Add vector (V) constraints
inline O findNormalizedLineNSphere(const V& origin, const V& direction, const V& center, T radius) noexcept
{
	V diff = origin - sphere.center;
	T halfB = dot(direction, diff);
	T c = dot(diff, diff) - sphere.radius*sphere.radius;
	T delta = halfB*halfB - c;

	if (delta < T(0))
	{
		return {};
	}
	else if (delta > T(0))
	{
		delta = std::sqrt(delta);
		return detail::interval<O>(-halfB - delta, -halfB + delta);
	}
	else // delta == 0
	{
		return detail::interval<O>(-halfB);
	}
}

template<typename T>
	requires std::floating_point<T>
inline bool testLineEllipsoid(const Vector3<T>& origin, const Vector3<T>& direction, const Vector3<T>& center, 
	const Matrix3<T>& matrix) noexcept
{
	// http://www.geometrictools.com/

	Vector3<T> diff = origin - center;
	Vector3<T> matDir = direction*matrix;
	Vector3<T> matDiff = diff*matrix;
	T a2 = dot(direction, matDir);
	T a1 = dot(direction, matDiff);
	T a0 = dot(diff, matDiff) - T(1);

	T discr = (a1*a1 - a0*a2);
	return (discr >= T(0));
}

template<typename O, typename T>
	requires std::floating_point<T>
inline O findLineEllipsoid(const Vector3<T>& origin, const Vector3<T>& direction, const Vector3<T>& center, 
	const Matrix3<T>& matrix) noexcept
{
	// http://www.geometrictools.com/

	Vector3<T> diff = origin - center;
	Vector3<T> matDir = direction*matrix;
	Vector3<T> matDiff = diff*matrix;
	T a2 = dot(direction, matDir);
	T a1 = dot(direction, matDiff);
	T a0 = dot(diff, matDiff) - T(1);

	T discr = a1*a1 - a0*a2;
	if (discr < T(0))
	{
		return {};
	}
	else if (discr > Constants<T>::TOLERANCE)
	{
		T root = std::sqrt(discr);
		T inv = T(1)/a2;
		return detail::interval<O>((-a1 - root)*inv, (-a1 + root)*inv);
	}
	else
	{
		return detail::interval<O>(-a1/a2);
	}
}

template<typename T>
	requires std::floating_point<T>
inline bool testRayEllipsoid(const Vector3<T>& origin, const Vector3<T>& direction, const Vector3<T>& center, 
	const Matrix3<T>& matrix) noexcept
{
	// http://www.geometrictools.com/

	Vector3<T> diff = origin - center;
	Vector3<T> matDir = direction*matrix;
	Vector3<T> matDiff = diff*matrix;
	T a2 = dot(direction, matDir);
	T a1 = dot(direction, matDiff);
	T a0 = dot(diff, matDiff) - T(1);

	T discr = a1*a1 - a0*a2;
	if (discr < T(0))
		return false;
	if (a0 <= T(0))
		return true;

	return (a1 < T(0));
}

template<typename O, typename T>
	requires std::floating_point<T>
inline O findRayEllipsoid(const Vector3<T>& origin, const Vector3<T>& direction, const Vector3<T>& center, 
	const Matrix3<T>& matrix) noexcept
{
	// http://www.geometrictools.com/

	Vector3<T> diff = origin - center;
	Vector3<T> matDir = direction*matrix;
	Vector3<T> matDiff = diff*matrix;
	T a2 = dot(direction, matDir);
	T a1 = dot(direction, matDiff);
	T a0 = dot(diff, matDiff) - T(1);

	T discr = a1*a1 - a0*a2;
	if (discr < T(0))
	{
		return {};
	}
	else if (discr > T(0))
	{
		T root = std::sqrt(discr);
		T inv = T(1)/a2;
		T t0 = (-a1 - root)*inv;
		T t1 = (-a1 + root)*inv;

		if (t0 >= T(0))
			return detail::interval<O>(t0, t1);
		else if (t1 >= T(0))
			return detail::interval<O>(t1);
		else
			return {};
	}
	else
	{
		T t0 = -a1/a2;
		if (t0 >= T(0))
			return detail::interval<O>(t0);
		else
			return {};
	}
}

template<typename T>
	requires std::floating_point<T>
inline bool testSegmentEllipsoid(const Vector3<T>& start, const Vector3<T>& end, const Vector3<T>& center, 
	const Matrix3<T>& matrix) noexcept
{
	// http://www.geometrictools.com/

	Vector3<T> direction = end - start;
	Vector3<T> diff = (start + end)*T(0.5) - center;
	Vector3<T> matDir = direction*matrix;
	Vector3<T> matDiff = diff*matrix;
	T a2 = dot(direction, matDir);
	T a1 = dot(direction, matDiff);
	T a0 = dot(diff, matDiff) - T(1);

	T discr = a1*a1 - a0*a2;
	if (discr < T(0))
		return false;
	if (a0 <= T(0))
		return true;

	constexpr T E = T(0.5);
	if (a1 >= T(0))
	{
		if (((a0 + E*(T(-2)*a1 + a2*E)) <= T(0)) || ((a1 - a2*E) < T(0)))
			return true;
	}
	else
	{
		if (((a0 + E*(T(2)*a1 + a2*E)) <= T(0)) || ((a1 + a2*E) < T(0)))
			return true;
	}

	return false;
}

template<typename O, typename T>
	requires std::floating_point<T>
inline O findSegmentEllipsoid(const Vector3<T>& start, const Vector3<T>& end, const Vector3<T>& center, 
	const Matrix3<T>& matrix) noexcept
{
	// http://www.geometrictools.com/

	Vector3<T> direction = end - start;
	Vector3<T> diff = (start + end)*T(0.5) - center;
	Vector3<T> matDir = direction*matrix;
	Vector3<T> matDiff = diff*matrix;
	T a2 = dot(direction, matDir);
	T a1 = dot(direction, matDiff);
	T a0 = dot(diff, matDiff) - T(1);

	constexpr T E = T(0.5);
	T discr = a1*a1 - a0*a2;
	if (discr < T(0))
	{
		return {};
	}
	else if (discr > Constants<T>::TOLERANCE)
	{
		T root = std::sqrt(discr);
		T inv = T(1)/a2;
		T t0 = (-a1 - root)*inv;
		T t1 = (-a1 + root)*inv;

		auto intersection = Interval<T>(t0, t1).findIntersection(Interval<T>(-e, e));
		if (intersection.has_value())
			return detail::interval<O>(intersection.value().minimum + E, intersection.value().maximum + E);
		else
			return {};
	}
	else
	{
		T t0 = -a1/a2;
		if (std::fabs(t0) <= E)
			return detail::interval<O>(t0 + E);
		else
			return {};
	}
}

template<typename T>
	requires std::floating_point<T>
inline bool testAxisAlignedBoxHalfSpace(const Vector3<T>& center, const Vector3<T>& halfDims, const Vector3<T>& normal, T constant) noexcept
{
	T r = sum(abs(halfDims*normal));
	return ((dot(normal, center) + constant) <= r);
}

template<typename T>
	requires std::floating_point<T>
inline int classifyAxisAlignedBoxHalfSpace(const Vector3<T>& center, const Vector3<T>& halfDims, const Vector3<T>& normal, T constant) noexcept
{
	T r = sum(abs(halfDims*normal));
	T a = dot(normal, center) + constant;
	return (a <= r) ? ((a <= -r) ? -1 : 0) : 1;
}

template<typename T>
	requires std::floating_point<T>
inline bool testOrientedBoxHalfSpace(const Vector3<T>& center, const Matrix3<T>& basis, const Vector3<T>& halfDims, const Vector3<T>& normal, 
	T constant) noexcept
{
	//Matrix3 basisTranspose = transpose(basis);
	//T r = sum(abs(halfDims*(normal*basisTranspose)));
	T r = sum(abs(halfDims*(basis*normal)));
	return ((dot(normal, center) + constant) <= r);
}

template<typename T>
	requires std::floating_point<T>
inline int classifyOrientedBoxHalfSpace(const Vector3<T>& center, const Matrix3<T>& basis, const Vector3<T>& halfDims, const Vector3<T>& normal,
	T constant) noexcept
{
	//Matrix3 basisTranspose = transpose(basis);
	//T r = sum(abs(halfDims*(normal*basisTranspose)));
	T r = sum(abs(halfDims*(basis*normal)));
	T a = dot(normal, center) + constant;
	return (a <= r) ? ((a <= -r) ? -1 : 0) : 1;
}

template<typename T>
	requires std::floating_point<T>
inline bool testAxisAlignedBoxPlane(const Vector3<T>& center, const Vector3<T>& halfDims, const Vector3<T>& normal, T constant) noexcept
{
	T r = sum(abs(halfDims*normal));
	return (std::fabs(dot(normal, center) + constant) <= r);
}

template<typename T>
	requires std::floating_point<T>
inline bool testOrientedBoxPlane(const Vector3<T>& center, const Matrix3<T>& basis, const Vector3<T>& halfDims, const Vector3<T>& normal, T constant) noexcept
{
	//Matrix3 basisTranspose = transpose(basis);
	//T r = sum(abs(halfDims*(normal*basisTranspose)));
	T r = sum(abs(halfDims*(basis*normal)));
	return (std::fabs(dot(normal, center) + constant) <= r);
}

template<typename T>
	requires std::floating_point<T>
bool testAxisAlignedBoxTriangle(const Vector3<T>& center, const Vector3<T>& halfDims, const Vector3<T>& vertex0, 
	const Vector3<T>& vertex1, const Vector3<T>& vertex2) noexcept
{
	// http://fileadmin.cs.lth.se/cs/Personal/Tomas_Akenine-Moller/code/

	// Move everything so that the boxcenter is in (0,0,0)
	Vector3<T> v0 = vertex0 - center;
	Vector3<T> v1 = vertex1 - center;
	Vector3<T> v2 = vertex2 - center;

	// Compute triangle edges
	Vector3<T> e0 = v1 - v0;
	Vector3<T> e1 = v2 - v1;
	Vector3<T> e2 = v0 - v2;

	T fex = std::fabs(e0.x);
	T fey = std::fabs(e0.y);
	T fez = std::fabs(e0.z);
	if (!detail::testAxisX01(e0.z, e0.y, fez, fey, halfDims, v0, v1, v2))
		return false;
	if (!detail::testAxisY02(e0.z, e0.x, fez, fex, halfDims, v0, v1, v2))
		return false;
	if (!detail::testAxisZ12(e0.y, e0.x, fey, fex, halfDims, v0, v1, v2))
		return false;

	fex = std::fabs(e1.x);
	fey = std::fabs(e1.y);
	fez = std::fabs(e1.z);
	if (!detail::testAxisX01(e1.z, e1.y, fez, fey, halfDims, v0, v1, v2))
		return false;
	if (!detail::testAxisY02(e1.z, e1.x, fez, fex, halfDims, v0, v1, v2))
		return false;
	if (!detail::testAxisZ0(e1.y, e1.x, fey, fex, halfDims, v0, v1, v2))
		return false;

	fex = std::fabs(e2.x);
	fey = std::fabs(e2.y);
	fez = std::fabs(e2.z);
	if (!detail::testAxisX2(e2.z, e2.y, fez, fey, halfDims, v0, v1, v2))
		return false;
	if (!detail::testAxisY1(e2.z, e2.x, fez, fex, halfDims, v0, v1, v2))
		return false;
	if (!detail::testAxisZ12(e2.y, e2.x, fey, fex, halfDims, v0, v1, v2))
		return false;

	T minimum = std::min(std::min(v0.x, v1.x), v2.x);
	T maximum = std::max(std::max(v0.x, v1.x), v2.x);
	if ((minimum > halfDims.x) || (maximum < -halfDims.x))
		return false;

	minimum = std::min(std::min(v0.y, v1.y), v2.y);
	maximum = std::max(std::max(v0.y, v1.y), v2.y);
	if ((minimum > halfDims.y) || (maximum < -halfDims.y))
		return false;

	minimum = std::min(std::min(v0.z, v1.z), v2.z);
	maximum = std::max(std::max(v0.z, v1.z), v2.z);
	if ((minimum > halfDims.z) || (maximum < -halfDims.z))
		return false;

	Vector3<T> normal = cross(e0, e1);
	//return (std::fabs(dot(normal, v0)) <= std::fabs(halfDims.x*normal.x) + std::fabs(halfDims.y*normal.y) + std::fabs(halfDims.z*normal.z));

	Vector3<T> minV = -halfDims;
	Vector3<T> maxV = halfDims;

	if (normal.x <= T(0))
		std::swap(minV.x, maxV.x);
	if (normal.y <= T(0))
		std::swap(minV.y, maxV.y);
	if (normal.z <= T(0))
		std::swap(minV.z, maxV.z);

	minV -= v0;
	maxV -= v0;

	if (dot(normal, minV) > T(0))
		return false;
	if (dot(normal, maxV) >= T(0))
		return true;

	return false;
}

template<typename T>
	requires std::floating_point<T>
bool testAxisAlignedBoxTriangle(/*const Vector3<T>& center,*/ const Vector3<T>& halfDims, const Vector3<T>& v0, const Vector3<T>& v1,
	const Vector3<T>& v2) noexcept
{
	// http://fileadmin.cs.lth.se/cs/Personal/Tomas_Akenine-Moller/code/

	// Move everything so that the boxcenter is in (0,0,0)
	//Vector3<T> v0 = vertex0 - center;
	//Vector3<T> v1 = vertex1 - center;
	//Vector3<T> v2 = vertex2 - center;

	// Compute triangle edges
	Vector3<T> e0 = v1 - v0;
	Vector3<T> e1 = v2 - v1;
	Vector3<T> e2 = v0 - v2;

	T fex = std::fabs(e0.x);
	T fey = std::fabs(e0.y);
	T fez = std::fabs(e0.z);
	if (!detail::testAxisX01(e0.z, e0.y, fez, fey, halfDims, v0, v1, v2))
		return false;
	if (!detail::testAxisY02(e0.z, e0.x, fez, fex, halfDims, v0, v1, v2))
		return false;
	if (!detail::testAxisZ12(e0.y, e0.x, fey, fex, halfDims, v0, v1, v2))
		return false;

	fex = std::fabs(e1.x);
	fey = std::fabs(e1.y);
	fez = std::fabs(e1.z);
	if (!detail::testAxisX01(e1.z, e1.y, fez, fey, halfDims, v0, v1, v2))
		return false;
	if (!detail::testAxisY02(e1.z, e1.x, fez, fex, halfDims, v0, v1, v2))
		return false;
	if (!detail::testAxisZ0(e1.y, e1.x, fey, fex, halfDims, v0, v1, v2))
		return false;

	fex = std::fabs(e2.x);
	fey = std::fabs(e2.y);
	fez = std::fabs(e2.z);
	if (!detail::testAxisX2(e2.z, e2.y, fez, fey, halfDims, v0, v1, v2))
		return false;
	if (!detail::testAxisY1(e2.z, e2.x, fez, fex, halfDims, v0, v1, v2))
		return false;
	if (!detail::testAxisZ12(e2.y, e2.x, fey, fex, halfDims, v0, v1, v2))
		return false;

	T minimum = std::min(std::min(v0.x, v1.x), v2.x);
	T maximum = std::max(std::max(v0.x, v1.x), v2.x);
	if ((minimum > halfDims.x) || (maximum < -halfDims.x))
		return false;

	minimum = std::min(std::min(v0.y, v1.y), v2.y);
	maximum = std::max(std::max(v0.y, v1.y), v2.y);
	if ((minimum > halfDims.y) || (maximum < -halfDims.y))
		return false;

	minimum = std::min(std::min(v0.z, v1.z), v2.z);
	maximum = std::max(std::max(v0.z, v1.z), v2.z);
	if ((minimum > halfDims.z) || (maximum < -halfDims.z))
		return false;

	Vector3<T> normal = cross(e0, e1);
	//return (std::fabs(dot(normal, v0)) <= std::fabs(halfDims.x*normal.x) + std::fabs(halfDims.y*normal.y) + std::fabs(halfDims.z*normal.z));

	Vector3<T> minV = -halfDims;
	Vector3<T> maxV = halfDims;

	if (normal.x <= T(0))
		std::swap(minV.x, maxV.x);
	if (normal.y <= T(0))
		std::swap(minV.y, maxV.y);
	if (normal.z <= T(0))
		std::swap(minV.z, maxV.z);

	minV -= v0;
	maxV -= v0;

	if (dot(normal, minV) > T(0))
		return false;
	if (dot(normal, maxV) >= T(0))
		return true;

	return false;
}

template<typename T>
	requires std::floating_point<T>
inline bool testOrientedBoxTriangle(const Vector3<T>& center, const Matrix3<T>& basis, const Vector3<T>& halfDims, 
	const Vector3<T>& vertex0, const Vector3<T>& vertex1, const Vector3<T>& vertex2) noexcept
{
	//Matrix3 basisTranspose = transpose(basis);
	//return detail::testAxisAlignedBoxTriangle(/*Vector3<T>::ZERO,*/ halfDims, (vertex0 - center)*basisTranspose, 
	//	(vertex1 - center)*basisTranspose, (vertex2 - center)*basisTranspose);
	return detail::testAxisAlignedBoxTriangle(/*Vector3<T>::ZERO,*/ halfDims, basis*(vertex0 - center),
		basis*(vertex1 - center), basis*(vertex2 - center));
}

template<typename T>
	requires std::floating_point<T>
inline bool testAxisAlignedRectangleCircle(const Vector2<T>& minimum, const Vector2<T>& maximum, const Vector2<T>& center, 
	T radius) noexcept
{
	T d = T(0);
	if (center.x < minimum.x)
		d += sqr(center.x - minimum.x);
	else if (center.x > maximum.x)
		d += sqr(center.x - maximum.x);
	if (center.y < minimum.y)
		d += sqr(center.y - minimum.y);
	else if (center.y > maximum.y)
		d += sqr(center.y - maximum.y);
	return (d <= radius*radius);
}

template<typename T>
	requires std::floating_point<T>
inline bool testAxisAlignedBoxSphere(const Vector3<T>& minimum, const Vector3<T>& maximum, const Vector3<T>& center, T radius) noexcept
{
	T d = T(0);
	if (center.x < minimum.x)
		d += sqr(center.x - minimum.x);
	else if (center.x > maximum.x)
		d += sqr(center.x - maximum.x);
	if (center.y < minimum.y)
		d += sqr(center.y - minimum.y);
	else if (center.y > maximum.y)
		d += sqr(center.y - maximum.y);
	if (center.z < minimum.z)
		d += sqr(center.z - minimum.z);
	else if (center.z > maximum.z)
		d += sqr(center.z - maximum.z);
	return (d <= radius*radius);
}

template<typename T>
	requires std::floating_point<T>
bool testOrientedBoxOrientedBox(const Vector3<T>& centerA, const Matrix3<T>& basisA, const Vector3<T>& halfDimsA,
	const Vector3<T>& centerB, const Matrix3<T>& basisB, const Vector3<T>& halfDimsB) noexcept
{
	// https://www.geometrictools.com/

	constexpr T CUTOFF = T(1) - Constants<T>::TOLERANCE;
	bool existsParallelPair = false;
	Vector3<T> d = centerB - centerA;
	T c[3][3];
	T absC[3][3];
	T aD[3];
	T r0, r1, r, r01;

	for (int i = 0; i < 3; ++i)
	{
		c[0][i] = dot(basisA[0], basisB[i]);
		absC[0][i] = std::fabs(c[0][i]);
		if (absC[0][i] > CUTOFF)
			existsParallelPair = true;
	}

	aD[0] = dot(basisA[0], d);
	r = std::fabs(aD[0]);
	r1 = dot(halfDimsB, Vector3<T>(absC[0][0], absC[0][1], absC[0][2]));
	r01 = halfDimsA.x + r1;
	if (r > r01)
		return false;

	for (int i = 0; i < 3; ++i)
	{
		c[1][i] = dot(basisA[1], basisB[i]);
		absC[1][i] = std::fabs(c[1][i]);
		if (absC[1][i] > CUTOFF)
			existsParallelPair = true;
	}

	aD[1] = dot(basisA[1], d);
	r = std::fabs(aD[1]);
	r1 = dot(halfDimsB, Vector3<T>(absC[1][0], absC[1][1], absC[1][2]));
	r01 = halfDimsA.y + r1;
	if (r > r01)
		return false;

	for (int i = 0; i < 3; ++i)
	{
		c[2][i] = dot(basisA[2], basisB[i]);
		absC[2][i] = std::fabs(c[2][i]);
		if (absC[2][i] > CUTOFF)
			existsParallelPair = true;
	}

	aD[2] = dot(basisA[2], d);
	r = std::fabs(aD[2]);
	r1 = dot(halfDimsB, Vector3<T>(absC[2][0], absC[2][1], absC[2][2]));
	r01 = halfDimsA.z + r1;
	if (r > r01)
		return false;

	r = std::fabs(dot(basisB[0], d));
	r0 = dot(halfDimsA, Vector3<T>(absC[0][0], absC[1][0], absC[2][0]));
	r01 = r0 + halfDimsB.x;
	if (r > r01)
		return false;

	r = std::fabs(dot(basisB[1], d));
	r0 = dot(halfDimsA, Vector3<T>(absC[0][1], absC[1][1], absC[2][1]));
	r01 = r0 + halfDimsB.y;
	if (r > r01)
		return false;

	r = std::fabs(dot(basisB[2], d));
	r0 = dot(halfDimsA, Vector3<T>(absC[0][2], absC[1][2], absC[2][2]));
	r01 = r0 + halfDimsB.z;
	if (r > r01)
		return false;

	if (existsParallelPair)
		return true;

	r = std::fabs(aD[2]*c[1][0] - aD[1]*c[2][0]);
	r0 = halfDimsA.y*absC[2][0] + halfDimsA.z*absC[1][0];
	r1 = halfDimsB.y*absC[0][2] + halfDimsB.z*absC[0][1];
	r01 = r0 + r1;
	if (r > r01)
		return false;

	r = std::fabs(aD[2]*c[1][1] - aD[1]*c[2][1]);
	r0 = halfDimsA.y*absC[2][1] + halfDimsA.z*absC[1][1];
	r1 = halfDimsB.x*absC[0][2] + halfDimsB.z*absC[0][0];
	r01 = r0 + r1;
	if (r > r01)
		return false;

	r = std::fabs(aD[2]*c[1][2] - aD[1]*c[2][2]);
	r0 = halfDimsA.y*absC[2][2] + halfDimsA.z*absC[1][2];
	r1 = halfDimsB.x*absC[0][1] + halfDimsB.y*absC[0][0];
	r01 = r0 + r1;
	if (r > r01)
		return false;

	r = std::fabs(aD[0]*c[2][0] - aD[2]*c[0][0]);
	r0 = halfDimsA.x*absC[2][0] + halfDimsA.z*absC[0][0];
	r1 = halfDimsB.y*absC[1][2] + halfDimsB.z*absC[1][1];
	r01 = r0 + r1;
	if (r > r01)
		return false;

	r = std::fabs(aD[0]*c[2][1] - aD[2]*c[0][1]);
	r0 = halfDimsA.x*absC[2][1] + halfDimsA.z*absC[0][1];
	r1 = halfDimsB.x*absC[1][2] + halfDimsB.z*absC[1][0];
	r01 = r0 + r1;
	if (r > r01)
		return false;

	r = std::fabs(aD[0]*c[2][2] - aD[2]*c[0][2]);
	r0 = halfDimsA.x*absC[2][2] + halfDimsA.z*absC[0][2];
	r1 = halfDimsB.x*absC[1][1] + halfDimsB.y*absC[1][0];
	r01 = r0 + r1;
	if (r > r01)
		return false;

	r = std::fabs(aD[1]*c[0][0] - aD[0]*c[1][0]);
	r0 = halfDimsA.x*absC[1][0] + halfDimsA.y*absC[0][0];
	r1 = halfDimsB.y*absC[2][2] + halfDimsB.z*absC[2][1];
	r01 = r0 + r1;
	if (r > r01)
		return false;

	r = std::fabs(aD[1]*c[0][1] - aD[0]*c[1][1]);
	r0 = halfDimsA.x*absC[1][1] + halfDimsA.y*absC[0][1];
	r1 = halfDimsB.x*absC[2][2] + halfDimsB.z*absC[2][0];
	r01 = r0 + r1;
	if (r > r01)
		return false;

	r = std::fabs(aD[1]*c[0][2] - aD[0]*c[1][2]);
	r0 = halfDimsA.x*absC[1][2] + halfDimsA.y*absC[0][2];
	r1 = halfDimsB.x*absC[2][1] + halfDimsB.y*absC[2][0];
	r01 = r0 + r1;
	if (r > r01)
		return false;

	return true;
}

template<typename T>
	requires std::floating_point<T>
inline bool testOrientedBoxAxisAlignedBox(const Vector3<T>& centerA, const Matrix3<T>& basisA, const Vector3<T>& halfDimsA,
	const Vector3<T>& centerB, const Vector3<T>& halfDimsB) noexcept
{
	return testOrientedBoxOrientedBox(centerA, basisA, halfDimsA, centerB, Matrix3<T>::IDENTITY, halfDimsB);
}

template<typename T>
	requires std::floating_point<T>
inline bool testOrientedBoxSphere(const Vector3<T>& centerA, const Matrix3<T>& basisA, const Vector3<T>& halfDimsA, 
	const Vector3<T>& centerB, T radiusB) noexcept
{
	//Matrix3<T> boxBasisT(transpose(basisA));
	return testAxisAlignedBoxSphere(-halfDimsA, halfDimsA,
		basisA*(centerB - centerA)/*(centerB - centerA)*boxBasisT*/, radiusB);
}

template<typename T>
	requires std::floating_point<T>
bool testOrientedBoxSymmetricFrustum(const Vector3<T>& centerA, const Matrix3<T>& basisA, const Vector3<T>& halfDimsA,
	const Vector3<T>& originB, const Matrix3<T>& basisB, const Vector2<T>& halfDimsB, T depthMin, T depthMax) noexcept
{
	// http://www.geometrictools.com/

	Vector3<T> diff = centerA - originB;
	T rBound = halfDimsB.x;
	T uBound = halfDimsB.y;
	T a[3];
	T b[3];
	T c[3];
	T d[3];
	T nA[3];
	T nB[3];
	T nC[3];
	T nD[3];
	T rC[3];
	T rD[3];
	T uC[3];
	T uD[3];
	T nApRC[3];
	T nAmRC[3];
	T nBpUC[3];
	T nBmUC[3];
	T rBpUA[3];
	T rBmUA[3];
	T ddD, radius, p, fMin, fMax, mTwoUF, mTwoRF, tmp;

	d[2] = dot(diff, basisB[2]);
	for (int i = 0; i < 3; ++i)
		c[i] = dot(basisA[i], basisB[2]);
	radius = halfDimsA[0]*std::fabs(c[0]) + halfDimsA[1]*std::fabs(c[1]) + halfDimsA[2]*std::fabs(c[2]);
	if (d[2] + radius < depthMin || d[2] - radius > depthMax)
		return false;

	for (int i = 0; i < 3; ++i)
	{
		a[i] = dot(basisA[i], basisB[0]);
		rC[i] = rBound*c[i];
		nA[i] = depthMin*a[i];
		nAmRC[i] = nA[i] - rC[i];
	}

	d[0] = dot(diff, basisB[0]);
	radius = halfDimsA[0]*std::fabs(nAmRC[0]) + halfDimsA[1]*std::fabs(nAmRC[1]) + halfDimsA[2]*std::fabs(nAmRC[2]);
	nD[0] = depthMin*d[0];
	rD[2] = rBound*d[2];
	ddD = nD[0] - rD[2];
	mTwoRF = -2.0f*rBound*depthMax;
	if (ddD + radius < mTwoRF || ddD > radius)
		return false;

	for (int i = 0; i < 3; ++i)
		nApRC[i] = nA[i] + rC[i];
	radius = halfDimsA[0]*std::fabs(nApRC[0]) + halfDimsA[1]*std::fabs(nApRC[1]) + halfDimsA[2]*std::fabs(nApRC[2]);
	ddD = -(nD[0] + rD[2]);
	if (ddD + radius < mTwoRF || ddD > radius)
		return false;

	for (int i = 0; i < 3; ++i)
	{
		b[i] = dot(basisA[i], basisB[1]);
		uC[i] = uBound*c[i];
		nB[i] = depthMin*b[i];
		nBmUC[i] = nB[i] - uC[i];
	}

	d[1] = dot(diff, basisB[1]);
	radius = halfDimsA[0]*std::fabs(nBmUC[0]) + halfDimsA[1]*std::fabs(nBmUC[1]) + halfDimsA[2]*std::fabs(nBmUC[2]);
	nD[1] = depthMin*d[1];
	uD[2] = uBound*d[2];
	ddD = nD[1] - uD[2];
	mTwoUF = -2.0f*uBound*depthMax;
	if (ddD + radius < mTwoUF || ddD > radius)
		return false;

	for (int i = 0; i < 3; ++i)
		nBpUC[i] = nB[i] + uC[i];
	radius = halfDimsA[0]*std::fabs(nBpUC[0]) + halfDimsA[1]*std::fabs(nBpUC[1]) + halfDimsA[2]*std::fabs(nBpUC[2]);
	ddD = -(nD[1] + uD[2]);
	if (ddD + radius < mTwoUF || ddD > radius)
		return false;

	T dRatio = depthMax/depthMin;

	for (int i = 0; i < 3; ++i)
	{
		p = rBound*std::fabs(a[i]) + uBound*std::fabs(b[i]);
		nC[i] = depthMin*c[i];
		fMin = nC[i] - p;
		if (fMin < T(0))
			fMin *= dRatio;
		fMax = nC[i] + p;
		if (fMax > T(0))
			fMax *= dRatio;
		ddD = a[i]*d[0] + b[i]*d[1] + c[i]*d[2];
		if (ddD + halfDimsA[i] < fMin || ddD - halfDimsA[i] > fMax)
			return false;
	}

	for (int i = 0; i < 3; ++i)
	{
		p = uBound*std::fabs(c[i]);
		fMin = -nB[i] - p;
		if (fMin < T(0))
			fMin *= dRatio;
		fMax = -nB[i] + p;
		if (fMax > T(0))
			fMax *= dRatio;
		ddD = c[i]*d[1] - b[i]*d[2];
		radius = halfDimsA[0]*std::fabs(b[i]*c[0] - b[0]*c[i]) + halfDimsA[1]*std::fabs(b[i]*c[1] - b[1]*c[i]) + halfDimsA[2]*std::fabs(b[i]*c[2] - b[2]*c[i]);
		if (ddD + radius < fMin || ddD - radius > fMax)
			return false;
	}

	for (int i = 0; i < 3; ++i)
	{
		p = rBound*std::fabs(c[i]);
		fMin = nA[i] - p;
		if (fMin < T(0))
			fMin *= dRatio;
		fMax = nA[i] + p;
		if (fMax > T(0))
			fMax *= dRatio;
		ddD = -c[i]*d[0] + a[i]*d[2];
		radius = halfDimsA[0]*std::fabs(a[i]*c[0] - a[0]*c[i]) + halfDimsA[1]*std::fabs(a[i]*c[1] - a[1]*c[i]) + halfDimsA[2]*std::fabs(a[i]*c[2] - a[2]*c[i]);
		if (ddD + radius < fMin || ddD - radius > fMax)
			return false;
	}

	for (int i = 0; i < 3; ++i)
	{
		T fRB = rBound*b[i];
		T fUA = uBound*a[i];
		rBpUA[i] = fRB + fUA;
		rBmUA[i] = fRB - fUA;
	}

	for (int i = 0; i < 3; ++i)
	{
		p = rBound*std::fabs(nBmUC[i]) + uBound*std::fabs(nAmRC[i]);
		tmp = -depthMin*rBmUA[i];
		fMin = tmp - p;
		if (fMin < T(0))
			fMin *= dRatio;
		fMax = tmp + p;
		if (fMax > T(0))
			fMax *= dRatio;
		ddD = d[0]*nBmUC[i] - d[1]*nAmRC[i] - d[2]*rBmUA[i];
		radius = T(0);
		for (int j = 0; j < 3; j++)
			radius += halfDimsA[j]*std::fabs(a[j]*nBmUC[i] - b[j]*nAmRC[i] - c[j]*rBmUA[i]);
		if (ddD + radius < fMin || ddD - radius > fMax)
			return false;
	}

	for (int i = 0; i < 3; ++i)
	{
		p = rBound*std::fabs(nBpUC[i]) + uBound*std::fabs(nAmRC[i]);
		tmp = -depthMin*rBpUA[i];
		fMin = tmp - p;
		if (fMin < T(0))
			fMin *= dRatio;
		fMax = tmp + p;
		if (fMax > T(0))
			fMax *= dRatio;
		ddD = d[0]*nBpUC[i] - d[1]*nAmRC[i] - d[2]*rBpUA[i];
		radius = T(0);
		for (int j = 0; j < 3; ++j)
			radius += halfDimsA[j]*std::fabs(a[j]*nBpUC[i] - b[j]*nAmRC[i] - c[j]*rBpUA[i]);
		if (ddD + radius < fMin || ddD - radius > fMax)
			return false;
	}

	for (int i = 0; i < 3; ++i)
	{
		p = rBound*std::fabs(nBmUC[i]) + uBound*std::fabs(nApRC[i]);
		tmp = depthMin*rBpUA[i];
		fMin = tmp - p;
		if (fMin < T(0))
			fMin *= dRatio;
		fMax = tmp + p;
		if (fMax > T(0))
			fMax *= dRatio;
		ddD = d[0]*nBmUC[i] - d[1]*nApRC[i] + d[2]*rBpUA[i];
		radius = T(0);
		for (int j = 0; j < 3; ++j)
			radius += halfDimsA[j]*std::fabs(a[j]*nBmUC[i] - b[j]*nApRC[i] + c[j]*rBpUA[i]);
		if (ddD + radius < fMin || ddD - radius > fMax)
			return false;
	}

	for (int i = 0; i < 3; ++i)
	{
		p = rBound*std::fabs(nBpUC[i]) + uBound*std::fabs(nApRC[i]);
		tmp = depthMin*rBmUA[i];
		fMin = tmp - p;
		if (fMin < T(0))
			fMin *= dRatio;
		fMax = tmp + p;
		if (fMax > T(0))
			fMax *= dRatio;
		ddD = d[0]*nBpUC[i] - d[1]*nApRC[i] + d[2]*rBmUA[i];
		radius = T(0);
		for (int j = 0; j < 3; ++j)
			radius += halfDimsA[j]*std::fabs(a[j]*nBpUC[i] - b[j]*nApRC[i] + c[j]*rBmUA[i]);
		if (ddD + radius < fMin || ddD - radius > fMax)
			return false;
	}

	return true;
}

template<typename T>
	requires std::floating_point<T>
inline bool testAxisAlignedBoxSymmetricFrustum(const Vector3<T>& centerA, const Vector3<T>& halfDimsA,
	const Vector3<T>& originB, const Matrix3<T>& basisB, const Vector2<T>& halfDimsB, T depthMin, T depthMax) noexcept
{
	return testOrientedBoxSymmetricFrustum(centerA, Matrix3<T>::IDENTITY, halfDimsA, originB, basisB, halfDimsB, depthMin, depthMax);
}

template<typename T>
	requires std::floating_point<T>
inline bool testEllipsoidPlane(const Vector3<T>& center, const Matrix3<T>& inverseMatrix, const Vector3<T>& normal, T constant) noexcept
{
	return (distances::getPointPlaneSquared(center, normal, constant) <= std::fabs(dot(normal, normal*inverseMatrix)));
}

template<typename T>
	requires std::floating_point<T>
inline bool testConeSphere(const Vector3<T>& vertex, const Vector3<T>& axis, T height, T baseRadius, const Vector3<T>& center, 
	T radius) noexcept
{
	// http://www.geometrictools.com/

	T slantHeight = std::sqrt(radius*radius + height*height);
	T sinAngle = radius/slantHeight;
	T cosAngle = height/slantHeight;
	T invSin = slantHeight/radius; // T(1)/sinAngle;
	T cosSqr = cosAngle*cosAngle;
	Vector3<T> cmV = sphere.center - vertex;
	Vector3<T> d = cmV + (sphere.radius*invSin)*axis;
	T dSqrLen = d.getMagnitudeSquared();
	T e = dot(d, axis);

	if ((e > T(0)) && (e*e >= dSqrLen*cosSqr))
	{
		T sinSqr = sinAngle*sinAngle;
		dSqrLen = cmV.getMagnitudeSquared();
		e = -dot(cmV, axis);
		if ((e > T(0)) && (e*e >= dSqrLen*sinSqr) && (dSqrLen > sphere.radius*sphere.radius))
			return false;
		if (/*(height < std::numeric_limits<T>::max()) &&*/ !((dot(axis, center) - dot(axis, vertex + height*axis)) <= radius))
			return false;

		return true;
	}

	return false;
}

} // namespace mathematics::intersections
