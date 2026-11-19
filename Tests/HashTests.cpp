#include "TestUtilities.hpp"
#include <functional>
#include <initializer_list>
#include <unordered_set>
#include <vector>

namespace mathematics::templates {

namespace {

template<typename V>
std::size_t hashOf(const V& v)
{
	return std::hash<V>()(v);
}

// Checks that pairwise distinct values produce pairwise distinct hashes, and that copies hash equally.
template<typename V>
void checkDistinct(std::initializer_list<V> values)
{
	std::unordered_set<std::size_t> hashes;
	for (const V& v : values)
	{
		const V copy(v);
		CHECK(hashOf(copy) == hashOf(v));
		hashes.insert(hashOf(v));
	}

	CHECK(hashes.size() == values.size());
}

} // namespace

// Regression tests: the SSE implementation of simd::hash used to XOR all the lanes together and then XOR the two
// (identical) halves of the result, returning 0 for every input.
TEST_SUITE("Hash")
{

TEST_CASE_TEMPLATE("Vector2", T, TEST_REAL_TYPES)
{
	using V = Vector2<T>;
	checkDistinct<V>({ V(T(0), T(0)), V(T(1), T(2)), V(T(2), T(1)), V(T(1), T(1)), V(T(2), T(2)), V(T(-1), T(2)), V(T(1), T(0)), 
		V(T(0), T(1)) });
	CHECK(V(T(0), T(-0.0)) == V(T(0), T(0)));
	CHECK(hashOf(V(T(0), T(-0.0))) == hashOf(V(T(0), T(0))));
}

TEST_CASE_TEMPLATE("Vector3", T, TEST_REAL_TYPES)
{
	using V = Vector3<T>;
	checkDistinct<V>({ V(T(0), T(0), T(0)), V(T(1), T(2), T(3)), V(T(4), T(5), T(6)), V(T(3), T(2), T(1)), V(T(2), T(1), T(3)), 
		V(T(1), T(1), T(1)), V(T(2), T(2), T(2)), V(T(1), T(2), T(-3)), V(T(1), T(0), T(0)), V(T(0), T(1), T(0)), 
		V(T(0), T(0), T(1)) });
	CHECK(V(T(-0.0), T(1), T(-0.0)) == V(T(0), T(1), T(0)));
	CHECK(hashOf(V(T(-0.0), T(1), T(-0.0))) == hashOf(V(T(0), T(1), T(0))));

	// Equal vectors built in different ways
	CHECK(hashOf(V(Vector2<T>(T(1), T(2)), T(3))) == hashOf(V(T(1), T(2), T(3))));
	CHECK(hashOf(V(Vector2<T>(T(1), T(2)))) == hashOf(V(T(1), T(2), T(0))));
	CHECK(hashOf(V(T(2))) == hashOf(V(T(2), T(2), T(2))));
	CHECK(hashOf(V(T(1), T(2), T(3)) + V(T(3), T(3), T(3))) == hashOf(V(T(4), T(5), T(6))));
}

TEST_CASE_TEMPLATE("Vector4", T, TEST_REAL_TYPES)
{
	using V = Vector4<T>;
	checkDistinct<V>({ V(T(0), T(0), T(0), T(0)), V(T(1), T(2), T(3), T(4)), V(T(4), T(3), T(2), T(1)), V(T(3), T(4), T(1), T(2)), 
		V(T(2), T(1), T(4), T(3)), V(T(1), T(1), T(1), T(1)), V(T(2), T(2), T(2), T(2)), V(T(1), T(1), T(2), T(2)), 
		V(T(2), T(2), T(1), T(1)), V(T(1), T(2), T(1), T(2)), V(T(0), T(0), T(0), T(1)), V(T(1), T(0), T(0), T(0)) });
	CHECK(hashOf(V(T(-0.0), T(1), T(2), T(-0.0))) == hashOf(V(T(0), T(1), T(2), T(0))));
}

TEST_CASE_TEMPLATE("Quaternion", T, TEST_REAL_TYPES)
{
	using Q = Quaternion<T>;
	checkDistinct<Q>({ Q(T(0), T(0), T(0), T(1)), Q(T(1), T(2), T(3), T(4)), Q(T(4), T(3), T(2), T(1)), Q(T(1), T(1), T(1), T(1)), 
		Q(T(1), T(0), T(0), T(0)), Q(T(0), T(1), T(0), T(0)) });
}

TEST_CASE_TEMPLATE("Plane and HalfSpace", T, TEST_REAL_TYPES)
{
	using P = Plane<T>;
	checkDistinct<P>({ P(T(1), T(0), T(0), T(0)), P(T(0), T(1), T(0), T(0)), P(T(0), T(0), T(1), T(0)), P(T(0), T(0), T(1), T(2)), 
		P(T(1), T(0), T(0), T(1)) });

	using H = HalfSpace<T>;
	checkDistinct<H>({ H(T(1), T(0), T(0), T(0)), H(T(0), T(1), T(0), T(0)), H(T(0), T(0), T(1), T(0)), H(T(0), T(0), T(1), T(2)), 
		H(T(1), T(0), T(0), T(1)) });
}

TEST_CASE_TEMPLATE("matrices", T, TEST_REAL_TYPES)
{
	using M2 = Matrix2<T>;
	checkDistinct<M2>({ M2::ZERO, M2::IDENTITY, M2(T(1), T(2), T(3), T(4)), M2(T(3), T(4), T(1), T(2)), M2(T(2), T(1), T(4), T(3)), 
		M2(T(1), T(1), T(1), T(1)), M2(T(1), T(2), T(1), T(2)) });
	CHECK(hashOf(M2(T(1), T(2), T(3), T(4))*M2::IDENTITY) == hashOf(M2(T(1), T(2), T(3), T(4))));

	using M3 = Matrix3<T>;
	const M3 m3(T(1), T(2), T(3), T(4), T(5), T(6), T(7), T(8), T(9));
	checkDistinct<M3>({ M3::ZERO, M3::IDENTITY, m3, transpose(m3), M3(T(4), T(5), T(6), T(1), T(2), T(3), T(7), T(8), T(9)), 
		M3(T(1), T(1), T(1), T(1), T(1), T(1), T(1), T(1), T(1)) });
	CHECK(hashOf(m3*M3::IDENTITY) == hashOf(m3));
	CHECK(hashOf(transpose(transpose(m3))) == hashOf(m3));

	using M4 = Matrix4<T>;
	const M4 m4(T(1), T(2), T(3), T(4), T(5), T(6), T(7), T(8), T(9), T(10), T(11), T(12), T(13), T(14), T(15), T(16));
	checkDistinct<M4>({ M4::ZERO, M4::IDENTITY, m4, transpose(m4) });
	CHECK(hashOf(transpose(transpose(m4))) == hashOf(m4));

	using A = AffineTransform<T>;
	const A a(Vector3<T>(T(1), T(2), T(3)), Vector3<T>(T(4), T(5), T(6)), Vector3<T>(T(7), T(8), T(9)), 
		Vector3<T>(T(10), T(11), T(12)));
	const A b(Vector3<T>(T(4), T(5), T(6)), Vector3<T>(T(1), T(2), T(3)), Vector3<T>(T(7), T(8), T(9)), 
		Vector3<T>(T(10), T(11), T(12)));
	checkDistinct<A>({ A::IDENTITY, a, b });
}

#if MATHEMATICS_HAS_SIMD && SIMD_HAS_FLOAT4

// The lanes of the SIMD register not occupied by a component may hold arbitrary data; such values still compare 
// equal, so they have to hash equally.
TEST_CASE("unused SIMD lanes are ignored")
{
	const Vector2<float> u0(1.f, 2.f);
	const Vector2<float> u1(simd::set4(1.f, 2.f, 7.f, -8.f));
	REQUIRE(u1 == u0);
	CHECK(hashOf(u1) == hashOf(u0));

	const Vector3<float> v0(1.f, 2.f, 3.f);
	const Vector3<float> v1(simd::set4(1.f, 2.f, 3.f, 9.f));
	REQUIRE(v1 == v0);
	CHECK(hashOf(v1) == hashOf(v0));

	const Matrix2<float> m0(1.f, 2.f, 3.f, 4.f);
	const Matrix2<float> m1(simd::set4(1.f, 2.f, 5.f, 6.f), simd::set4(3.f, 4.f, 7.f, 8.f));
	REQUIRE(m1 == m0);
	CHECK(hashOf(m1) == hashOf(m0));

	const Matrix3<float> n0(1.f, 2.f, 3.f, 4.f, 5.f, 6.f, 7.f, 8.f, 9.f);
	Matrix3<float> n1(n0);
	n1.row0 = simd::set4(1.f, 2.f, 3.f, -1.f);
	n1.row1 = simd::set4(4.f, 5.f, 6.f, -2.f);
	n1.row2 = simd::set4(7.f, 8.f, 9.f, -3.f);
	REQUIRE(n1 == n0);
	CHECK(hashOf(n1) == hashOf(n0));

	const AffineTransform<float> a0(Vector3<float>(1.f, 2.f, 3.f), Vector3<float>(4.f, 5.f, 6.f), Vector3<float>(7.f, 8.f, 9.f),
		Vector3<float>(10.f, 11.f, 12.f));
	const AffineTransform<float> a1(simd::set4(1.f, 2.f, 3.f, -1.f), simd::set4(4.f, 5.f, 6.f, -2.f), 
		simd::set4(7.f, 8.f, 9.f, -3.f), simd::set4(10.f, 11.f, 12.f, -4.f));
	REQUIRE(a1 == a0);
	CHECK(hashOf(a1) == hashOf(a0));
}

TEST_CASE("simd::float4")
{
	const std::hash<simd::float4> hasher;
	std::unordered_set<std::size_t> hashes;
	const float values[] = { 0.f, 1.f, 2.f, -1.f };
	for (float x : values)
		for (float y : values)
			for (float z : values)
				for (float w : values)
					hashes.insert(hasher(simd::set4(x, y, z, w)));

	CHECK(hashes.size() == 256);
	CHECK(hasher(simd::set4(-0.f, 1.f, -0.f, 2.f)) == hasher(simd::set4(0.f, 1.f, 0.f, 2.f)));
}

#endif

} // TEST_SUITE

} // namespace mathematics::templates
