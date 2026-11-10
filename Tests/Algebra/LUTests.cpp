#include "TestUtilities.hpp"
#include <Mathematics/Algebra/LU.inl>

namespace mathematics::templates {

TEST_SUITE("LU")
{

TEST_CASE_TEMPLATE("3x3 system", T, TEST_REAL_TYPES)
{
	// Leading zero forces partial pivoting.
	Matrix3<T> m(T(0), T(2), T(1), T(1), T(1), T(1), T(2), T(1), T(3));
	auto d = lu::decompose(m);
	REQUIRE(d.has_value());

	Vector3<T> b(T(1), T(2), T(3));
	Vector3<T> x = lu::solve(*d, b);
	CHECK_APPROX(m*x, b); // solves m*x = b for a column vector x
	CHECK(lu::solve<Vector3<T>>(*d) == Vector3<T>::ZERO);

	CHECK_APPROX(lu::inverse<Matrix3<T>>(*d), inverse(m));
}

TEST_CASE_TEMPLATE("4x4 system", T, TEST_REAL_TYPES)
{
	Matrix4<T> m(T(4), T(-2), T(1), T(0),
		T(3), T(6), T(-4), T(2),
		T(2), T(1), T(8), T(-5),
		T(1), T(0), T(2), T(7));
	auto d = lu::decompose(m);
	REQUIRE(d.has_value());

	Vector4<T> b(T(1), T(-2), T(3), T(0.5));
	CHECK_APPROX(m*lu::solve(*d, b), b);
	CHECK_APPROX(lu::inverse<Matrix4<T>>(*d), inverse(m));
}

TEST_CASE_TEMPLATE("raw arrays", T, TEST_REAL_TYPES)
{
	std::array<T, 4> m{ T(0), T(1), T(2), T(3) };
	auto d = lu::decompose<std::array<T, 4>, 2, T>(m);
	REQUIRE(d.has_value());
	std::array<T, 2> x = lu::solve<std::array<T, 2>, 2, T>(*d, std::array<T, 2>{ T(1), T(8) });
	CHECK_APPROX(x[0], T(2.5));
	CHECK_APPROX(x[1], T(1));
	std::array<T, 4> inv = lu::inverse<std::array<T, 4>, 2, T>(*d);
	CHECK_APPROX(inv[0], T(-1.5));
	CHECK_APPROX(inv[1], T(0.5));
	CHECK_APPROX(inv[2], T(1));
	CHECK_APPROX(inv[3], T(0));
}

TEST_CASE_TEMPLATE("matrix with a zero row is rejected", T, TEST_REAL_TYPES)
{
	Matrix3<T> m(T(1), T(2), T(3), T(0), T(0), T(0), T(4), T(5), T(6));
	CHECK_FALSE(lu::decompose(m).has_value());
}

TEST_CASE_TEMPLATE("singular matrices", T, TEST_REAL_TYPES)
{
	// Linearly dependent rows without a zero row
	CHECK_FALSE(lu::decompose(Matrix3<T>(T(1), T(2), T(3), T(2), T(4), T(6), T(1), T(0), T(1))).has_value());
	CHECK_FALSE(lu::decompose(Matrix3<T>(T(1), T(2), T(3), T(4), T(5), T(6), T(7), T(8), T(9))).has_value());
	CHECK_FALSE(lu::decompose(Matrix4<T>(T(1), T(0), T(0), T(0), T(0), T(1), T(0), T(0), T(1), T(1), T(0), T(0), T(0), T(0), T(0), T(1))).has_value());

	// Two parallel plane normals, as in ConvexPolyhedron::computeVertices()
	CHECK_FALSE(lu::decompose(Matrix3<T>(Vector3<T>::UNIT_X, Vector3<T>::MINUS_UNIT_X, Vector3<T>::UNIT_Y)).has_value());

	// Non-finite input
	CHECK_FALSE(lu::decompose(Matrix3<T>(T(1), T(0), T(0), T(0), std::numeric_limits<T>::quiet_NaN(), T(0), T(0), T(0), T(1))).has_value());
}

TEST_CASE_TEMPLATE("tolerance", T, TEST_REAL_TYPES)
{
	// The relative pivot of the last column equals delta.
	const T delta = std::numeric_limits<T>::epsilon();
	Matrix3<T> m(T(1), T(0), T(0), T(0), T(1), T(0), T(1), T(1), delta);
	CHECK_FALSE(lu::decompose(m).has_value()); // below the default tolerance of 3 epsilon

	auto d = lu::decompose(m, delta*T(0.5));
	REQUIRE(d.has_value());
	CHECK_APPROX(lu::solve(*d, Vector3<T>(T(1), T(0), T(1) + delta)), Vector3<T>(T(1), T(0), T(1))); // 1 + epsilon is exact

	// A looser tolerance rejects ill-conditioned matrices (e.g. nearly parallel planes).
	Matrix3<T> nearlyParallel(Vector3<T>::UNIT_X, normalize(Vector3<T>(T(1), T(1e-3), T(0))), Vector3<T>::UNIT_Z);
	CHECK(lu::decompose(nearlyParallel).has_value());
	CHECK_FALSE(lu::decompose(nearlyParallel, T(1e-2)).has_value());
}

TEST_CASE_TEMPLATE("row scaling does not affect singularity detection", T, TEST_REAL_TYPES)
{
	Matrix3<T> m(T(0), T(2), T(1), T(1), T(1), T(1), T(2), T(1), T(3));
	for (T scale : { T(1e-20), T(1e20) })
	{
		CAPTURE(scale);
		Matrix3<T> scaled = m*scale;
		auto d = lu::decompose(scaled);
		REQUIRE(d.has_value());
		CHECK_APPROX(scaled*lu::inverse<Matrix3<T>>(*d), Matrix3<T>::IDENTITY);
	}

	// Rows of very different magnitudes: the solution is still recovered accurately
	Matrix3<T> rows(m[0]*T(1e-15), m[1], m[2]*T(1e15));
	auto d = lu::decompose(rows);
	REQUIRE(d.has_value());
	Vector3<T> x(T(1), T(-2), T(3));
	CHECK_APPROX(lu::solve(*d, rows*x), x);

	// ...and a singular matrix stays singular (power of two scales keep the scaled matrix exact)
	Matrix3<T> singular(T(1), T(2), T(3), T(4), T(5), T(6), T(7), T(8), T(9));
	CHECK_FALSE(lu::decompose(Matrix3<T>(singular[0]*std::ldexp(T(1), -50), singular[1], singular[2]*std::ldexp(T(1), 50))).has_value());
}

TEST_CASE_TEMPLATE("determinant", T, TEST_REAL_TYPES)
{
	Matrix3<T> m3(T(0), T(2), T(1), T(1), T(1), T(1), T(2), T(1), T(3)); // requires a row swap
	CHECK_APPROX(lu::determinant(*lu::decompose(m3)), m3.getDeterminant());
	CHECK_APPROX(lu::determinant(*lu::decompose(Matrix3<T>::IDENTITY)), T(1));

	Matrix4<T> m4(T(4), T(-2), T(1), T(0),
		T(3), T(6), T(-4), T(2),
		T(2), T(1), T(8), T(-5),
		T(1), T(0), T(2), T(7));
	CHECK_APPROX(lu::determinant(*lu::decompose(m4)), m4.getDeterminant());
	Matrix4<T> p(T(0), T(1), T(0), T(0), T(0), T(0), T(1), T(0), T(1), T(0), T(0), T(0), T(0), T(0), T(0), T(1));
	CHECK_APPROX(lu::determinant(*lu::decompose(p)), p.getDeterminant()); // odd/even permutation sign
	CHECK_APPROX(lu::determinant(*lu::decompose(Matrix4<T>(T(0), T(1), T(0), T(0), T(1), T(0), T(0), T(0), T(0), T(0), T(1), T(0),
		T(0), T(0), T(0), T(1)))), T(-1));
}

} // TEST_SUITE

} // namespace mathematics::templates
