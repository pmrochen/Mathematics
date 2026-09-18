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
#include <tuple>
#include <cstddef>
#include <cmath>
#include <malloc.h>

namespace mathematics::triangulation {
namespace detail {

template<std::random_access_iterator I>
inline typename std::tuple_element<0, typename std::iterator_traits<I>::value_type>::type computeSignedPolygonArea(I firstVertex, I lastVertex) noexcept
{
	using VectorType = typename std::iterator_traits<I>::value_type;
	using ScalarType = typename std::tuple_element<0, VectorType>::type;

	std::ptrdiff_t nVertices = std::distance(firstVertex, lastVertex);
	if (nVertices < 3)
		return ScalarType(0);

	ScalarType area = cross(firstVertex[nVertices - 1], firstVertex[0]);
	for (std::ptrdiff_t i = 0, n = nVertices - 1; i < n; i++)
		area += cross(firstVertex[i], firstVertex[i + 1]);
	
	return area*ScalarType(0.5);
}

template<std::random_access_iterator I, std::random_access_iterator J>
inline bool snip(I firstVertex, I lastVertex, J firstIndex, std::ptrdiff_t u, std::ptrdiff_t v, std::ptrdiff_t w) noexcept
{
	static_assert(std::is_integral_v<typename std::iterator_traits<J>::value_type>);
	using VectorType = typename std::iterator_traits<I>::value_type;
	using ScalarType = typename std::tuple_element<0, VectorType>::type;

    const VectorType& a = firstVertex[firstIndex[u]];
    const VectorType& b = firstVertex[firstIndex[v]];
    const VectorType& c = firstVertex[firstIndex[w]];
    if (cross(b - a, c - a) < std::numeric_limits<ScalarType>::epsilon())
        return false;

	VectorType bc = c - b;
	VectorType ab = b - a;
	VectorType ca = a - c;

	std::ptrdiff_t nVertices = std::distance(firstVertex, lastVertex);
    for (std::ptrdiff_t i = 0; i < nVertices; i++)
    {
        if ((i == u) || (i == v) || (i == w))
            continue;

        const VectorType& p = firstVertex[firstIndex[i]];
		if ((cross(bc, p - b) >= ScalarType(0)) && (cross(ca, p - c) >= ScalarType(0)) && (cross(ab, p - a) >= ScalarType(0)))
            return false;
    }

    return true;
}

} // namespace detail

template<std::random_access_iterator I, std::integral U, std::output_iterator<U> O>
std::pair<O, bool> triangulate(I firstVertex, I lastVertex, O outIndex)
{
	using VectorType = typename std::iterator_traits<I>::value_type;
	using ScalarType = typename std::tuple_element<0, VectorType>::type;
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

	bool reverse = (computeSignedPolygonArea(firstVertex, lastVertex) < ScalarType(0));
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

//template<std::random_access_iterator I, std::random_access_iterator J, std::integral U,
//	std::output_iterator<U> O>
//std::pair<O, bool> triangulate(I firstVertex, I lastVertex, J firstIndex, J lastIndex, O outIndex)
//{
//}

} // namespace mathematics::triangulation
