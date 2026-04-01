/*
 *	Name: EulerOrder
 *	Author: Pawel Mrochen
 */

#pragma once

#include <tuple>
#include "Axis.hpp"

namespace mathematics {

enum class EulerOrder
{
    UNSPECIFIED = 0,
    XYZ	= ((((((0 << 1) | 0) << 1) | 0) << 1) | 1),
    XZY	= ((((((0 << 1) | 1) << 1) | 0) << 1) | 1),
    YZX	= ((((((1 << 1) | 0) << 1) | 0) << 1) | 1),
    YXZ	= ((((((1 << 1) | 1) << 1) | 0) << 1) | 1),
    ZXY	= ((((((2 << 1) | 0) << 1) | 0) << 1) | 1),	// YawPitchRoll
    ZYX	= ((((((2 << 1) | 1) << 1) | 0) << 1) | 1),
};

inline std::tuple<Axis, Axis, Axis> unpack(EulerOrder order) noexcept
{
	static const int safe[] = { 0, 1, 2, 0 };
	static const int next[] = { 1, 2, 0, 1 };
	unsigned int o = (unsigned int)order;
	int f = o & 1; o >>= 1;
	int s = o & 1; o >>= 1;
	int n = o & 1; o >>= 1;
	int i = safe[o & 3];
	int j = next[i + n];
	int k = next[i + 1 - n];

	return { (Axis)(f ? i : (s ? i : k)), (Axis)j, (Axis)(f ? (s ? i : k) : i) };
	//staticFrame = !f;
}

inline EulerOrder pack(Axis axis1, Axis axis2) noexcept
{
	static const EulerOrder orders[6] =
	{
		EulerOrder::XYZ,
		EulerOrder::XZY,
		EulerOrder::YZX,
		EulerOrder::YXZ,
		EulerOrder::ZXY,
		EulerOrder::ZYX
	};

	for (int i = 0; i < 6; i++)
	{
		auto [a1, a2, a3] = unpack(orders[i]);
		if ((a1 == axis1) && (a2 == axis2))
			return orders[i];
	}

	return EulerOrder::UNSPECIFIED;
}

inline EulerOrder pack(Axis axis1, Axis axis2, Axis axis3) noexcept
{
	static const EulerOrder orders[6] =
	{
		EulerOrder::XYZ,
		EulerOrder::XZY,
		EulerOrder::YZX,
		EulerOrder::YXZ,
		EulerOrder::ZXY,
		EulerOrder::ZYX
	};

	for (int i = 0; i < 6; i++)
	{
		auto [a1, a2, a3] = unpack(orders[i]);
		if ((a1 == axis1) && (a2 == axis2) && (a3 == axis3))
			return orders[i];
	}

	return EulerOrder::UNSPECIFIED;
}

} // namespace mathematics
