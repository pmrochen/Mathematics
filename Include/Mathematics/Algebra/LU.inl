/*
 *	Name: LU
 *	Author: Pawel Mrochen
 */

#pragma once

#include <limits>
#include <type_traits>
#include <concepts>
#include <utility>
#include <optional>
#include <array>
#include <algorithm>
#include <cmath>

namespace mathematics::lu {

template<int Order, std::floating_point T>
using Decomposition = std::pair<std::array<T, Order*Order>, std::array<int, Order>>;

// Returns no value if the matrix is singular, i.e. if a pivot relative to the largest element of its row (implicit
// scaling) does not exceed the tolerance. The default tolerance only rejects matrices that are singular up to rounding.
template<typename M, int Order = M::DIMENSION, std::floating_point T = typename M::Real>
std::optional<Decomposition<Order, T>> decompose(const M& matrix,
	std::type_identity_t<T> tolerance = T(Order)*std::numeric_limits<T>::epsilon()) noexcept
{
	// https://numerical.recipes/

	std::array<T, Order*Order> a;
	if constexpr (std::is_same_v<M, std::array<T, Order*Order>>)
	{
		a = matrix;
	}
	else
	{
		for (int i = 0; i < Order; i++)
		{
			for (int j = 0; j < Order; j++)
				a[i*Order + j] = matrix[i][j];
		}
	}

	T vv[Order];
	for (int i = 0; i < Order; i++)
	{
		T big = T(0);
		for (int j = 0; j < Order; j++)
		{
			T t = std::fabs(a[i*Order + j]);
			if (t > big)
				big = t;
		}

		if (big == T(0))
			return {}; // the matrix is singular

		vv[i] = T(1)/big;
	}

	std::array<int, Order> index;
	for (int j = 0; j < Order; j++)
	{
		int iMax = j;
		for (int i = 0; i < j; i++)
		{
			T sum = a[i*Order + j];
			for (int k = 0; k < i; k++)
				sum -= a[i*Order + k]*a[k*Order + j];
			a[i*Order + j] = sum;
		}

		T big = T(0);
		for (int i = j; i < Order; i++)
		{
			T sum = a[i*Order + j];
			for (int k = 0; k < j; k++)
				sum -= a[i*Order + k]*a[k*Order + j];
			a[i*Order + j] = sum;

			T dum = vv[i]*std::fabs(sum);
			if (dum >= big)
			{
				big = dum;
				iMax = i;
			}
		}

		if (!(big > tolerance)) // relative pivot too small (also rejects NaN)
			return {};

		if (j != iMax)
		{
			for (int k = 0; k < Order; k++)
				std::swap(a[iMax*Order + k], a[j*Order + k]);
			vv[iMax] = vv[j];
		}

		index[j] = iMax;

		if (j != (Order - 1))
		{
			T dum = T(1)/a[j*Order + j];
			for (int i = j + 1; i < Order; i++)
				a[i*Order + j] *= dum;
		}
	}

	return { { a, index } };
}

template<std::size_t Order, std::floating_point T> // std::size_t, so that Order is deducible from the decomposition
T determinant(const std::pair<std::array<T, Order*Order>, std::array<int, Order>>& decomposition) noexcept
{
	const auto& [a, index] = decomposition;
	T det = T(1);
	for (int i = 0; i < int(Order); i++)
	{
		det *= a[i*int(Order) + i];
		if (index[i] != i) // each recorded row swap flips the sign
			det = -det;
	}

	return det;
}

template<typename V, int Order = V::NUM_COMPONENTS, std::floating_point T = typename V::Real>
V solve(const std::type_identity_t<Decomposition<Order, T>>& decomposition, const V& vector) noexcept
{
	// https://numerical.recipes/

	const auto& [a, index] = decomposition;
	std::array<T, Order> b;
	for (int i = 0; i < Order; i++)
		b[i] = vector[i];

	int ii = -1;
	for (int i = 0; i < Order; i++)
	{
		int ip = index[i];
		T sum = b[ip];
		b[ip] = b[i];

		if (ii >= 0)
		{
			for (int j = ii; j < i; j++)
				sum -= a[i*Order + j]*b[j];
		}
		else if (sum != T(0))
		{
			ii = i;
		}

		b[i] = sum;
	}

	for (int i = Order - 1; i >= 0; i--)
	{
		T sum = b[i];
		for (int j = i + 1; j < Order; j++)
			sum -= a[i*Order + j]*b[j];
		b[i] = sum/a[i*Order + i];
	}

	if constexpr (std::is_same_v<V, std::array<T, Order>>)
		return b;
	else
		return V(&b[0]);
}

template<typename V, int Order = V::NUM_COMPONENTS, std::floating_point T = typename V::Real>
V solve(const std::type_identity_t<Decomposition<Order, T>>& decomposition) noexcept
{
	return solve<V, Order, T>(decomposition, V::ZERO);
}

template<typename M, int Order = M::DIMENSION, std::floating_point T = typename M::Real>
M inverse(const std::type_identity_t<Decomposition<Order, T>>& decomposition) noexcept
{
	std::array<T, Order*Order> a;
	std::array<T, Order> col;
	for (int j = 0; j < Order; j++)
	{
		for (int i = 0; i < Order; i++)
			col[i] = (i == j) ? T(1) : T(0);
		
		col = solve<std::array<T, Order>, Order, T>(decomposition, col);
		for (int i = 0; i < Order; i++)
			a[i*Order + j] = col[i];
	}

	if constexpr (std::is_same_v<M, std::array<T, Order*Order>>)
		return a;
	else
		return M(&a[0]);
}

} // namespace mathematics::lu
