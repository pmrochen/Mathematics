/*
 *	Name: Perlin
 *	Author: Pawel Mrochen
 */

#pragma once

#include <concepts>
#include <cmath>
#include "Scalar.hpp"
#include "Vector2.hpp"
#include "Vector3.hpp"

namespace mathematics::perlin {
namespace detail {

inline int perm(int x) noexcept
{
	static const unsigned char lut[512] = 
	{ 
		151,160,137,91,90,15,
		131,13,201,95,96,53,194,233,7,225,140,36,103,30,69,142,8,99,37,240,21,10,23,
		190, 6,148,247,120,234,75,0,26,197,62,94,252,219,203,117,35,11,32,57,177,33,
		88,237,149,56,87,174,20,125,136,171,168, 68,175,74,165,71,134,139,48,27,166,
		77,146,158,231,83,111,229,122,60,211,133,230,220,105,92,41,55,46,245,40,244,
		102,143,54, 65,25,63,161, 1,216,80,73,209,76,132,187,208, 89,18,169,200,196,
		135,130,116,188,159,86,164,100,109,198,173,186, 3,64,52,217,226,250,124,123,
		5,202,38,147,118,126,255,82,85,212,207,206,59,227,47,16,58,17,182,189,28,42,
		223,183,170,213,119,248,152, 2,44,154,163, 70,221,153,101,155,167, 43,172,9,
		129,22,39,253, 19,98,108,110,79,113,224,232,178,185, 112,104,218,246,97,228,
		251,34,242,193,238,210,144,12,191,179,162,241, 81,51,145,235,249,14,239,107,
		49,192,214, 31,181,199,106,157,184, 84,204,176,115,121,50,45,127, 4,150,254,
		138,236,205,93,222,114,67,29,24,72,243,141,128,195,78,66,215,61,156,180,
		151,160,137,91,90,15,
		131,13,201,95,96,53,194,233,7,225,140,36,103,30,69,142,8,99,37,240,21,10,23,
		190, 6,148,247,120,234,75,0,26,197,62,94,252,219,203,117,35,11,32,57,177,33,
		88,237,149,56,87,174,20,125,136,171,168, 68,175,74,165,71,134,139,48,27,166,
		77,146,158,231,83,111,229,122,60,211,133,230,220,105,92,41,55,46,245,40,244,
		102,143,54, 65,25,63,161, 1,216,80,73,209,76,132,187,208, 89,18,169,200,196,
		135,130,116,188,159,86,164,100,109,198,173,186, 3,64,52,217,226,250,124,123,
		5,202,38,147,118,126,255,82,85,212,207,206,59,227,47,16,58,17,182,189,28,42,
		223,183,170,213,119,248,152, 2,44,154,163, 70,221,153,101,155,167, 43,172,9,
		129,22,39,253, 19,98,108,110,79,113,224,232,178,185, 112,104,218,246,97,228,
		251,34,242,193,238,210,144,12,191,179,162,241, 81,51,145,235,249,14,239,107,
		49,192,214, 31,181,199,106,157,184, 84,204,176,115,121,50,45,127, 4,150,254,
		138,236,205,93,222,114,67,29,24,72,243,141,128,195,78,66,215,61,156,180 
	};

	return (int)lut[x];
}

template<std::floating_point T>
inline T fade(T t) noexcept
{ 
	return t*t*t*(t*(t*T(6) - T(15)) + T(10)); 
}

template<std::floating_point T>
inline T grad(int hash, T x, T y, T z) noexcept
{
	int h = hash & 15;
	T u = (h < 8) ? x : y;
	T v = (h < 4) ? y : (h == 12) || (h == 14) ? x : z;
	return (((h & 1) == 0) ? u : -u) + (((h & 2) == 0) ? v : -v);
}

} // namespace detail

template<std::floating_point T>
inline T noise(T x) noexcept
{
	// http://mrl.nyu.edu/~perlin/noise/

	using namespace detail;

	T fx = std::floor(x);
	int ix = (int)fx & 255;
    x -= fx;
	T u = fade(x);
    
	int a = perm(ix);
	int aa = perm(a);
    int b = perm(ix + 1);
	int ba = perm(b);

    return lerp(grad(perm(aa), x, T(0), T(0)), grad(perm(ba), x - T(1), T(0), T(0)), u);
}

template<std::floating_point T>
inline T noise(T x, T y) noexcept
{
	// http://mrl.nyu.edu/~perlin/noise/

	using namespace detail;

	T fx = std::floor(x);
	T fy = std::floor(y);
	int ix = (int)fx & 255;
    int iy = (int)fy & 255;
    x -= fx;
    y -= fy;
	T u = fade(x);
    T v = fade(y);
    
	int a = perm(ix) + iy;
	int aa = perm(a);
	int ab = perm(a + 1);
    int b = perm(ix + 1) + iy;
	int ba = perm(b);
	int bb = perm(b + 1);

    return lerp(lerp(grad(perm(aa), x, y, T(0)), grad(perm(ba), x - T(1), y, T(0)), u),
                lerp(grad(perm(ab), x, y - T(1), T(0)), grad(perm(bb), x - T(1), y - T(1), T(0)), u), v);
}

template<std::floating_point T>
inline T noise(T x, T y, T z) noexcept
{
	// http://mrl.nyu.edu/~perlin/noise/

	using namespace detail;

	T fx = std::floor(x);
	T fy = std::floor(y);
	T fz = std::floor(z);
	int ix = (int)fx & 255;
    int iy = (int)fy & 255;
    int iz = (int)fz & 255;
    x -= fx;
    y -= fy;
    z -= fz;
	T u = fade(x);
    T v = fade(y);
    T w = fade(z);
    
	int a = perm(ix) + iy;
	int aa = perm(a) + iz;
	int ab = perm(a + 1) + iz;
    int b = perm(ix + 1) + iy;
	int ba = perm(b) + iz;
	int bb = perm(b + 1) + iz;

    return lerp(lerp(lerp(grad(perm(aa), x, y, z), grad(perm(ba), x - T(1), y, z), u),
                     lerp(grad(perm(ab), x, y - T(1), z), grad(perm(bb), x - T(1), y - T(1), z), u), v),
                lerp(lerp(grad(perm(aa + 1), x, y, z - T(1)), grad(perm(ba + 1), x - T(1), y, z - T(1)), u),
                     lerp(grad(perm(ab + 1), x, y - T(1), z - T(1)), grad(perm(bb + 1), x - T(1), y - T(1), z - T(1)), u), v), w);
}

template<std::floating_point T>
inline T noise(const templates::Vector2<T>& v) noexcept
{
	return noise(v.x, v.y);
}

template<std::floating_point T>
inline T noise(const templates::Vector3<T>& v) noexcept
{
	return noise(v.x, v.y, v.z);
}

} // namespace mathematics::perlin
