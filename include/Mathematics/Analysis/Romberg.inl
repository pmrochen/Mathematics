/*
 *	Name: Romberg
 *	Author: Pawel Mrochen
 */

#pragma once

#include <concepts>

namespace mathematics::romberg {

template<int Order, std::floating_point T, typename F>
T estimate(T a, T b, F f) noexcept
{
	T h = b - a;

	T rom[2][Order];
	rom[0][0] = T(0.5)*h*(f(a) + f(b));
	for (int i0 = 2, ip0 = 1; i0 <= Order; i0++, ip0 *= 2, h *= T(0.5))
	{
		T sum = T(0.0);
		for (int i1 = 1; i1 <= ip0; i1++)
			sum += f(a + h*(i1 - T(0.5)));

		rom[1][0] = T(0.5)*(rom[0][0] + h*sum);
		for (int i2 = 1, ip2 = 4; i2 < i0; i2++, ip2 *= 4)
			rom[1][i2] = (ip2*rom[1][i2 - 1] - rom[0][i2 - 1])/(ip2 - 1);

		for (int i1 = 0; i1 < i0; i1++)
			rom[0][i1] = rom[1][i1];
	}

	return rom[0][Order - 1];
}

} // namespace mathematics::romberg
