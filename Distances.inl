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
    return square((dot(normal, point) + constant)/normal.getMagnitude());
}

template<typename T>
    requires std::floating_point<T>
inline T getPointNormalizedPlaneSquared(const Vector3<T>& point, const Vector3<T>& normal, T constant) noexcept
{
    return square(dot(normal, point) + constant);
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

template<typename T>
    requires std::floating_point<T>
T getPointSymmetricFrustumSquared(const Vector3<T>& point, const Vector3<T>& origin, const Matrix3<T>& basis, const Vector2<T>& halfDims,
    T depthMin, T depthMax, Vector3<T>* closestPoint = nullptr) noexcept
{
    // http://www.geometrictools.com/

	Vector3<T> diff = point - origin;
	Vector3<T> test = basis*diff;
	//Vector3<T> test(dot(diff, basis[0]), dot(diff, basis[1]), dot(diff, basis[2]));

	bool rSignChange;
	if (test.x < T(0))
	{
		rSignChange = true;
		test.x = -test.x;
	}
	else
	{
		rSignChange = false;
	}

	bool uSignChange;
	if (test.y < T(0))
	{
		uSignChange = true;
		test.y = -test.y;
	}
	else
	{
		uSignChange = false;
	}

	T dRatio = depthMax/depthMin;
	T rMin = halfDims.x;
	T rMax = dRatio*rMin;
	T uMin = halfDims.y;
	T uMax = dRatio*uMin;
	T rMinSqr = rMin*rMin;
	T uMinSqr = uMin*uMin;
	T dminSqr = depthMin*depthMin;
	T minRDDot = rMinSqr + dminSqr;
	T minUDDot = uMinSqr + dminSqr;
	T minRUDDot = rMinSqr + minUDDot;
	T maxRDDot = dRatio*minRDDot;
	T maxUDDot = dRatio*minUDDot;
	T maxRUDDot = dRatio*minRUDDot;

	Vector3<T> closest(Uninitialized());
	T rDot, uDot, rdDot, udDot, rudDot, rEdgeDot, uEdgeDot, t;
	if (test.z >= depthMax)
	{
		if (test.x <= rMax)
		{
			if (test.y <= uMax)
				closest.set(test.x, test.y, depthMax);
			else
				closest.set(test.x, uMax, depthMax);
		}
		else
		{
			if (test.y <= uMax)
				closest.set(rMax, test.y, depthMax);
			else
				closest.set(rMax, uMax, depthMax);
		}
	}
	else if (test.z <= depthMin)
	{
		if (test.x <= rMin)
		{
			if (test.y <= uMin)
			{
				closest.set(test.x, test.y, depthMin);
			}
			else
			{
				udDot = uMin*test.y + depthMin*test.z;
				if (udDot >= maxUDDot)
				{
					closest.set(test.x, uMax, depthMax);
				}
				else if (udDot >= minUDDot)
				{
					uDot = depthMin*test.y - uMin*test.z;
					t = uDot/minUDDot;
					closest.set(test.x, test.y - t*depthMin, test.z + t*uMin);
				}
				else
				{
					closest.set(test.x, uMin, depthMin);
				}
			}
		}
		else
		{
			if (test.y <= uMin)
			{
				rdDot = rMin*test.x + depthMin*test.z;
				if (rdDot >= maxRDDot)
				{
					closest.set(rMax, test.y, depthMax);
				}
				else if (rdDot >= minRDDot)
				{
					rDot = depthMin*test.x - rMin*test.z;
					t = rDot/minRDDot;
					closest.set(test.x - t*depthMin, test.y, test.z + t*rMin);
				}
				else
				{
					closest.set(rMin, test.y, depthMin);
				}
			}
			else
			{
				rudDot = rMin*test.x + uMin*test.y + depthMin*test.z;
				rEdgeDot = uMin*rudDot - minRUDDot*test.y;
				if (rEdgeDot >= T(0))
				{
					rdDot = rMin*test.x + depthMin*test.z;
					if (rdDot >= maxRDDot)
					{
						closest.set(rMax, test.y, depthMax);
					}
					else if (rdDot >= minRDDot)
					{
						rDot = depthMin*test.x - rMin*test.z;
						t = rDot/minRDDot;
						closest.set(test.x - t*depthMin, test.y, test.z + t*rMin);
					}
					else
					{
						closest.set(rMin, test.y, depthMin);
					}
				}
				else
				{
					uEdgeDot = rMin*rudDot - minRUDDot*test.x;
					if (uEdgeDot >= T(0))
					{
						udDot = uMin*test.y + depthMin*test.z;
						if (udDot >= maxUDDot)
						{
							closest.set(test.x, uMax, depthMax);
						}
						else if (udDot >= minUDDot)
						{
							uDot = depthMin*test.y - uMin*test.z;
							t = uDot/minUDDot;
							closest.set(test.x, test.y - t*depthMin, test.z + t*uMin);
						}
						else
						{
							closest.set(test.x, uMin, depthMin);
						}
					}
					else
					{
						if (rudDot >= maxRUDDot)
						{
							closest.set(rMax, uMax, depthMax);
						}
						else if (rudDot >= minRUDDot)
						{
							t = rudDot/minRUDDot;
							closest.set(t*rMin, t*uMin, t*depthMin);
						}
						else
						{
							closest.set(rMin, uMin, depthMin);
						}
					}
				}
			}
		}
	}
	else
	{
		rDot = depthMin*test.x - rMin*test.z;
		uDot = depthMin*test.y - uMin*test.z;
		if (rDot <= T(0))
		{
			if (uDot <= T(0))
			{
				closest = test;
			}
			else
			{
				udDot = uMin*test.y + depthMin*test.z;
				if (udDot >= maxUDDot)
				{
					closest.set(test.x, uMax, depthMax);
				}
				else
				{
					t = uDot/minUDDot;
					closest.set(test.x, test.y - t*depthMin, test.z + t*uMin);
				}
			}
		}
		else
		{
			if (uDot <= T(0))
			{
				rdDot = rMin*test.x + depthMin*test.z;
				if (rdDot >= maxRDDot)
				{
					closest.set(rMax, test.y, depthMax);
				}
				else
				{
					t = rDot/minRDDot;
					closest.set(test.x - t*depthMin, test.y, test.z + t*rMin);
				}
			}
			else
			{
				rudDot = rMin*test.x + uMin*test.y + depthMin*test.z;
				rEdgeDot = uMin*rudDot - minRUDDot*test.y;
				if (rEdgeDot >= T(0))
				{
					rdDot = rMin*test.x + depthMin*test.z;
					if (rdDot >= maxRDDot)
					{
						closest.set(rMax, test.y, depthMax);
					}
					else
					{
						t = rDot/minRDDot;
						closest.set(test.x - t*depthMin, test.y, test.z + t*rMin);
					}
				}
				else
				{
					uEdgeDot = rMin*rudDot - minRUDDot*test.x;
					if (uEdgeDot >= T(0))
					{
						udDot = uMin*test.y + depthMin*test.z;
						if (udDot >= maxUDDot)
						{
							closest.set(test.x, uMax, depthMax);
						}
						else
						{
							t = uDot/minUDDot;
							closest.set(test.x, test.y - t*depthMin, test.z + t*uMin);
						}
					}
					else
					{
						if (rudDot >= maxRUDDot)
						{
							closest.set(rMax, uMax, depthMax);
						}
						else
						{
							t = rudDot/minRUDDot;
							closest.set(t*rMin, t*uMin, t*depthMin);
						}
					}
				}
			}
		}
	}

	diff = test - closest;

	if (rSignChange)
		closest.x = -closest.x;
	if (uSignChange)
		closest.y = -closest.y;

	if (closestPoint)
		*closestPoint = origin + closest*basis;

	return diff.getLengthSquared();
}

template<typename T>
    requires std::floating_point<T>
inline T getPointSymmetricFrustum(const Vector3<T>& point, const Vector3<T>& origin, const Matrix3<T>& basis, const Vector2<T>& halfDims,
    T depthMin, T depthMax, Vector3<T>* closestPoint = nullptr) noexcept
{
    return std::sqrt(getPointSymmetricFrustumSquared(point, origin, basis, halfDims, depthMin, depthMax, closestPoint));
}

} // namespace mathematics::distances
