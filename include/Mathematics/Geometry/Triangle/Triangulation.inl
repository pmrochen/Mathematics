/*
 *	Name: Triangulation
 *	Author: Pawel Mrochen
 */

#pragma once

#include <limits>
#include <type_traits>
#include <concepts>
#include <iterator>
#include <utility>
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
	std::ptrdiff_t nVertices = std::distance(firstVertex, lastVertex);
	if (nVertices < 3)
		return T(0);

	T area = cross(firstVertex[nVertices - 1], firstVertex[0]);
	for (std::ptrdiff_t i = 0, n = nVertices - 1; i < n; i++)
		area += cross(firstVertex[i], firstVertex[i + 1]);
	
	return area*T(0.5);
}

template<std::floating_point T, std::random_access_iterator<Vector2<T>> I, std::integral U, std::random_access_iterator<U> J>
inline bool snip(I firstVertex, I lastVertex, J firstIndex, std::ptrdiff_t u, std::ptrdiff_t v, std::ptrdiff_t w) noexcept
{
    const Vector2<T>& a = firstVertex[firstIndex[u]];
    const Vector2<T>& b = firstVertex[firstIndex[v]];
    const Vector2<T>& c = firstVertex[firstIndex[w]];
    if (cross(b - a, c - a) < std::numeric_limits<T>::epsilon())
        return false;

	Vector2<T> bc = c - b;
	Vector2<T> ab = b - a;
	Vector2<T> ca = a - c;

	std::ptrdiff_t nVertices = std::distance(firstVertex, lastVertex);
    for (std::ptrdiff_t i = 0; i < nVertices; i++)
    {
        if ((i == u) || (i == v) || (i == w))
            continue;

        const Vector2<T>& p = firstVertex[firstIndex[i]];
		if ((cross(bc, p - b) >= T(0)) && (cross(ca, p - c) >= T(0)) && (cross(ab, p - a) >= T(0)))
            return false;
    }

    return true;
}

} // namespace detail

template<std::floating_point T, std::random_access_iterator<Vector2<T>> I, std::integral U, std::output_iterator<U> O>
std::pair<O, bool> triangulate2(I firstVertex, I lastVertex, O outIndex)
{
	//using IndexType = typename std::iterator_traits<O>::value_type;

	std::ptrdiff_t nVertices = std::distance(firstVertex, lastVertex);
	if (nVertices <= 3)
	{
		for (std::ptrdiff_t i = 0; i < nVertices; i++)
			*outIndex++ = U(i);
		
		return { outIndex, false };
	}

	//constexpr std::size_t BUFFER_SIZE = 64;
	//U indexBuffer[BUFFER_SIZE];
	//U* vertexIndices = (nVertices > BUFFER_SIZE) ? new U[nVertices] : indexBuffer;
	U* vertexIndices = (U*)alloca(nVertices*sizeof(U));

	bool reverse = (computeSignedPolygonArea(firstVertex, lastVertex) < T(0));
	if (reverse)
	{
		for (std::ptrdiff_t i = 0; i < nVertices; i++)
			vertexIndices[i] = U(nVertices - 1 - i);
	}
	else
	{
		for (std::ptrdiff_t i = 0; i < nVertices; i++)
			vertexIndices[i] = U(i);
	}

	std::ptrdiff_t nVerticesRemaining = nVertices;
	std::ptrdiff_t error = nVerticesRemaining*2;

	// Remove nVertices-2 vertices, creating 1 triangle every time
	for (std::ptrdiff_t v = nVerticesRemaining - 1; nVerticesRemaining > 2; )
	{
		// If we loop, it is probably a non-simple polygon
		if (std::ptrdiff_t(0) >= (error--))
		{
			//throw std::runtime_error("triangulate2() : bad polygon");

			for (std::ptrdiff_t i = 0, n = nVertices - 2; i < n; i++)
			{
				*outIndex++ = U(i);
				*outIndex++ = U(i + 1);
				*outIndex++ = U(i + 2);
			}

			//if (nVertices > BUFFER_SIZE)
			//	delete[] vertexIndices;

			return { outIndex, false };
		}

		// Three consecutive vertices in current polygon, <u,v,w>
		std::ptrdiff_t u = v;
		if (u >= nVerticesRemaining)
			u = 0;
		v = u + 1;
		if (v >= nVerticesRemaining)
			v = 0;
		std::ptrdiff_t w = v + 1;
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
			for (std::ptrdiff_t i = v; (i + 1) < nVerticesRemaining; i++)
				vertexIndices[i] = vertexIndices[i + 1];
			nVerticesRemaining--;

			// Resest error detection counter
			error = nVerticesRemaining*2;
		}
	}

	//if (nVertices > BUFFER_SIZE)
	//	delete[] vertexIndices;

	return { outIndex, true };
}

//template<std::floating_point T, std::random_access_iterator<Vector2<T>> I, std::integral U, std::random_access_iterator<U> J,
//	std::output_iterator<U> O>
//std::pair<O, bool> triangulate2(I firstVertex, I lastVertex, J firstIndex, J lastIndex, O outIndex)
//{
//}
//
//template<std::floating_point T, std::random_access_iterator<Vector3<T>> I, std::integral U, std::output_iterator<U> O>
//std::pair<O, bool> triangulate3(I firstVertex, I lastVertex, O outIndex)
//{
//}
//
//template<std::floating_point T, std::random_access_iterator<Vector3<T>> I, std::integral U, std::random_access_iterator<U> J, 
//	std::output_iterator<U> O>
//std::pair<O, bool> triangulate3(I firstVertex, I lastVertex, J firstIndex, J lastIndex, O outIndex)
//{
//}

} // namespace mathematics::triangulation
