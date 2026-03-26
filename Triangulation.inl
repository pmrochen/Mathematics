/*
 *	Name: Triangulation
 *	Author: Pawel Mrochen
 */

#pragma once

#include <limits>
#include <type_traits>
#include <concepts>
#include <iterator>
#include <cstddef>
#include <cmath>
#include <malloc.h>
#include "Vector2.hpp"
#include "Vector3.hpp"

namespace mathematics::triangulation {
namespace detail {

template<std::floating_point T, std::random_access_iterator<Vector2<T>> I>
inline T computeSignedPolygonArea(I firstVertex, I lastVertex) noexcept
{
	std::size_t nVertices = (std::size_t)std::distance(firstVertex, lastVertex);
	if (nVertices < 3u)
		return T(0);

	T area = firstVertex[nVertices - 1].x*firstVertex[0].y - firstVertex[0].x*firstVertex[nVertices - 1].y;
	for (std::size_t i = 0; i != (nVertices - 1); i++)
		area += firstVertex[i].x*firstVertex[i + 1].y - firstVertex[i + 1].x*firstVertex[i].y;
	return area*T(0.5);
}

template<std::floating_point T, std::random_access_iterator<Vector2<T>> I, std::integral U, std::random_access_iterator<U> J>
inline bool snip(I firstVertex, I lastVertex, J firstIndex, std::size_t u, std::size_t v, std::size_t w) noexcept
{
    const Vector2<T>& a = firstVertex[firstIndex[u]];
    const Vector2<T>& b = firstVertex[firstIndex[v]];
    const Vector2<T>& c = firstVertex[firstIndex[w]];
    if ((((b.x - a.x)*(c.y - a.y)) - ((b.y - a.y)*(c.x - a.x))) < std::numeric_limits<T>::epsilon())
        return false;

	std::size_t nVertices = (std::size_t)std::distance(firstVertex, lastVertex);
    for (std::size_t i = 0; i != nVertices; i++)
    {
        if ((i == u) || (i == v) || (i == w))
            continue;

        const Vector2<T>& p = firstVertex[firstIndex[i]];
		T axbp = (c.x - b.x)*(p.y - b.y) - (c.y - b.y)*(p.x - b.x);
		T cxap = (b.x - a.x)*(p.y - a.y) - (b.y - a.y)*(p.x - a.x);
		T bxcp = (a.x - c.x)*(p.y - c.y) - (a.y - c.y)*(p.x - c.x);
		if ((axbp >= T(0)) && (bxcp >= T(0)) && (cxap >= T(0)))
            return false;
    }

    return true;
}

} // namespace detail

template<std::floating_point T, std::random_access_iterator<Vector2<T>> I, std::integral U, std::output_iterator<U> O>
O triangulate2(I firstVertex, I lastVertex, O outIndex)
{
	//using IndexType = typename std::iterator_traits<O>::value_type;

	std::size_t nVertices = (std::size_t)std::distance(firstVertex, lastVertex);
	if (nVertices <= 3u)
	{
		for (std::size_t i = 0; i < nVertices; i++)
			*outIndex++ = U(i);
		
		return outIndex;
	}

	//constexpr std::size_t BUFFER_SIZE = 64;
	//U indexBuffer[BUFFER_SIZE];
	//U* vertexIndices = (nVertices > BUFFER_SIZE) ? new U[nVertices] : indexBuffer;
	U* vertexIndices = (U*)alloca(nVertices*sizeof(U));

	bool reverse = (computeSignedPolygonArea(firstVertex, lastVertex) < T(0));
	if (reverse)
	{
		for (std::size_t i = 0; i != nVertices; i++)
			vertexIndices[i] = U(nVertices - 1 - i);
	}
	else
	{
		for (std::size_t i = 0; i != nVertices; i++)
			vertexIndices[i] = U(i);
	}

	std::size_t nVerticesRemaining = nVertices;
	std::ptrdiff_t error = nVerticesRemaining*2;

	// Remove nVertices-2 vertices, creating 1 triangle every time
	for (std::size_t v = nVerticesRemaining - 1; nVerticesRemaining > 2; )
	{
		// If we loop, it is probably a non-simple polygon
		if (0 >= (error--))
		{
			//throw std::runtime_error("triangulate2() : bad polygon");

			for (std::size_t i = 0, n = nVertices - 2; i != n; i++)
			{
				*outIndex++ = U(i);
				*outIndex++ = U(i + 1);
				*outIndex++ = U(i + 2);
			}

			//if (nVertices > BUFFER_SIZE)
			//	delete[] vertexIndices;

			return outIndex;
		}

		// Three consecutive vertices in current polygon, <u,v,w>
		std::size_t u = v;
		if (u >= nVerticesRemaining)
			u = 0;
		v = u + 1;
		if (v >= nVerticesRemaining)
			v = 0;
		std::size_t w = v + 1;
		if (w >= nVerticesRemaining)
			w = 0;

		if (snip(firstVertex, firstVertex + nVerticesRemaining, vertexIndices, u, v, w))
		{
			// Output triangle
			if (reverse)
			{
				*outIndex++ = vertexIndices[w];
				*outIndex++ = vertexIndices[v];
				*outIndex++ = vertexIndices[u];
			}
			else
			{
				*outIndex++ = vertexIndices[u];
				*outIndex++ = vertexIndices[v];
				*outIndex++ = vertexIndices[w];
			}

			// Remove v from remaining polygon
			for (std::size_t i = v; (i + 1) < nVerticesRemaining; i++)
				vertexIndices[i] = vertexIndices[i + 1];
			nVerticesRemaining--;

			// Resest error detection counter
			error = nVerticesRemaining*2;
		}
	}

	//if (nVertices > BUFFER_SIZE)
	//	delete[] vertexIndices;

	return outIndex;
}

//template<std::floating_point T, std::random_access_iterator<Vector2<T>> I, std::integral U, std::random_access_iterator<U> J,
//	std::output_iterator<U> O>
//O triangulate2(I firstVertex, I lastVertex, J firstIndex, J lastIndex, O outIndex)
//{
//}
//
//template<std::floating_point T, std::random_access_iterator<Vector3<T>> I, std::integral U, std::output_iterator<U> O>
//O triangulate3(I firstVertex, I lastVertex, O outIndex)
//{
//}
//
//template<std::floating_point T, std::random_access_iterator<Vector3<T>> I, std::integral U, std::random_access_iterator<U> J, 
//	std::output_iterator<U> O>
//O triangulate3(I firstVertex, I lastVertex, J firstIndex, J lastIndex, O outIndex)
//{
//}

} // namespace mathematics::triangulation
