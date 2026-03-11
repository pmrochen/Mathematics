/*
 *	Name: Distances
 *	Author: Pawel Mrochen
 */

#pragma once

#include <algorithm>
#include <cmath>
#include "Scalar.hpp"
#include "Vector3.hpp"

namespace mathematics::distances {

template<typename T>
    requires std::floating_point<T>
inline T getPointHalfSpace(const Vector3<T>& point, const Vector3<T>& normal, T constant) noexcept
{
    return std::max((dot(normal, point) + constant)/normal.getMagnitude(), T(0));
}

template<typename T>
    requires std::floating_point<T>
inline T getPointNormalizedHalfSpace(const Vector3<T>& point, const Vector3<T>& normal, T constant) noexcept
{
    return std::max(dot(normal, point) + constant, T(0));
}

template<typename T>
    requires std::floating_point<T>
inline T getPointPlane(const Vector3<T>& point, const Vector3<T>& normal, T constant) noexcept
{
    return std::fabs((dot(normal, point) + constant)/normal.getMagnitude());
}

template<typename T>
    requires std::floating_point<T>
inline T getPointNormalizedPlane(const Vector3<T>& point, const Vector3<T>& normal, T constant) noexcept
{
    return std::fabs(dot(normal, point) + constant);
}

template<typename T>
    requires std::floating_point<T>
inline T getPointPlaneSigned(const Vector3<T>& point, const Vector3<T>& normal, T constant) noexcept
{
    return (dot(normal, point) + constant)/normal.getMagnitude();
}

template<typename T>
    requires std::floating_point<T>
inline T getPointNormalizedPlaneSigned(const Vector3<T>& point, const Vector3<T>& normal, T constant) noexcept
{
    return dot(normal, point) + constant;
}

template<typename T>
    requires std::floating_point<T>
inline T getPointPlaneSquared(const Vector3<T>& point, const Vector3<T>& normal, T constant) noexcept
{
    return sqr((dot(normal, point) + constant)/normal.getMagnitude());
}

template<typename T>
    requires std::floating_point<T>
inline T getPointNormalizedPlaneSquared(const Vector3<T>& point, const Vector3<T>& normal, T constant) noexcept
{
    return sqr(dot(normal, point) + constant);
}

template<typename T>
    requires std::floating_point<T>
T getPointTriangleSquared(const Vector3<T>& point, const Vector3<T>& v0, const Vector3<T>& v1, const Vector3<T>& v2, 
    Vector3<T>* closestPoint = nullptr) noexcept
{
	// http://www.geometrictools.com/
	
    Vector3<T> diff = v0 - point;
    Vector3<T> edge0 = v1 - v0;
    Vector3<T> edge1 = v2 - v0;
    T a00 = lengthSquared(edge0);
    T a01 = dot(edge0, edge1);
    T a11 = lengthSquared(edge1);
    T b0 = dot(diff, edge0);
    T b1 = dot(diff, edge1);
    T c = lengthSquared(diff);
    T det = std::fabs(a00*a11 - a01*a01);
    T s = a01*b1 - a11*b0;
    T t = a01*b0 - a00*b1;
    T sqrDistance;

    if (s + t <= det)
    {
        if (s < T(0))
        {
            if (t < T(0))
            {
                if (b0 < T(0))
                {
                    t = T(0);
                    if (-b0 >= a00)
                    {
                        s = T(1);
                        sqrDistance = a00 + T(2)*b0 + c;
                    }
                    else
                    {
                        s = -b0/a00;
                        sqrDistance = b0*s + c;
                    }
                }
                else
                {
                    s = T(0);
                    if (b1 >= T(0))
                    {
                        t = T(0);
                        sqrDistance = c;
                    }
                    else if (-b1 >= a11)
                    {
                        t = T(1);
                        sqrDistance = a11 + T(2)*b1 + c;
                    }
                    else
                    {
                        t = -b1/a11;
                        sqrDistance = b1*t + c;
                    }
                }
            }
            else
            {
                s = T(0);
                if (b1 >= T(0))
                {
                    t = T(0);
                    sqrDistance = c;
                }
                else if (-b1 >= a11)
                {
                    t = T(1);
                    sqrDistance = a11 + T(2)*b1 + c;
                }
                else
                {
                    t = -b1/a11;
                    sqrDistance = b1*t + c;
                }
            }
        }
        else if (t < T(0))
        {
            t = T(0);
            if (b0 >= T(0))
            {
                s = T(0);
                sqrDistance = c;
            }
            else if (-b0 >= a00)
            {
                s = T(1);
                sqrDistance = a00 + T(2)*b0 + c;
            }
            else
            {
                s = -b0/a00;
                sqrDistance = b0*s + c;
            }
        }
        else
        {
            // Minimum at interior point
            T invDet = T(1)/det;
            s *= invDet;
            t *= invDet;
            sqrDistance = s*(a00*s + a01*t + T(2)*b0) + t*(a01*s + a11*t + T(2)*b1) + c;
        }
    }
    else
    {
        if (s < T(0))
        {
            T tmp0 = a01 + b0;
            T tmp1 = a11 + b1;
            if (tmp1 > tmp0)
            {
                T numer = tmp1 - tmp0;
                T denom = a00 - T(2)*a01 + a11;
                if (numer >= denom)
                {
                    s = T(1);
                    t = T(0);
                    sqrDistance = a00 + T(2)*b0 + c;
                }
                else
                {
                    s = numer/denom;
                    t = T(1) - s;
                    sqrDistance = s*(a00*s + a01*t + T(2)*b0) + t*(a01*s + a11*t + T(2)*b1) + c;
                }
            }
            else
            {
                s = T(0);
                if (tmp1 <= T(0))
                {
                    t = T(1);
                    sqrDistance = a11 + T(2)*b1 + c;
                }
                else if (b1 >= T(0))
                {
                    t = T(0);
                    sqrDistance = c;
                }
                else
                {
                    t = -b1/a11;
                    sqrDistance = b1*t + c;
                }
            }
        }
        else if (t < T(0))
        {
            T tmp0 = a01 + b1;
            T tmp1 = a00 + b0;
            if (tmp1 > tmp0)
            {
                T numer = tmp1 - tmp0;
                T denom = a00 - T(2)*a01 + a11;
                if (numer >= denom)
                {
                    t = T(1);
                    s = T(0);
                    sqrDistance = a11 + T(2)*b1 + c;
                }
                else
                {
                    t = numer/denom;
                    s = T(1) - t;
                    sqrDistance = s*(a00*s + a01*t + T(2)*b0) + t*(a01*s + a11*t + T(2)*b1) + c;
                }
            }
            else
            {
                t = T(0);
                if (tmp1 <= T(0))
                {
                    s = T(1);
                    sqrDistance = a00 + T(2)*b0 + c;
                }
                else if (b0 >= T(0))
                {
                    s = T(0);
                    sqrDistance = c;
                }
                else
                {
                    s = -b0/a00;
                    sqrDistance = b0*s + c;
                }
            }
        }
        else
        {
            T numer = a11 + b1 - a01 - b0;
            if (numer <= T(0))
            {
                s = T(0);
                t = T(1);
                sqrDistance = a11 + T(2)*b1 + c;
            }
            else
            {
                T denom = a00 - T(2)*a01 + a11;
                if (numer >= denom)
                {
                    s = T(1);
                    t = T(0);
                    sqrDistance = a00 + T(2)*b0 + c;
                }
                else
                {
                    s = numer/denom;
                    t = T(1) - s;
                    sqrDistance = s*(a00*s + a01*t + T(2)*b0) + t*(a01*s + a11*t + T(2)*b1) + c;
                }
            }
        }
    }

    // Account for numerical round-off error
    if (sqrDistance < T(0))
        sqrDistance = T(0);

	if (closestPoint)
		*closestPoint = v0 + s*edge0 + t*edge1;

    //triangleBary[1] = s;
    //triangleBary[2] = t;
    //triangleBary[0] = T(1) - s - t;
    
	return sqrDistance;
}

template<typename T>
    requires std::floating_point<T>
inline T getPointTriangle(const Vector3<T>& point, const Vector3<T>& v0, const Vector3<T>& v1, const Vector3<T>& v2,
    Vector3<T>* closestPoint = nullptr) noexcept
{
    return std::sqrt(getPointTriangleSquared(point, v0, v1, v2, closestPoint));
}

} // namespace mathematics::distances
