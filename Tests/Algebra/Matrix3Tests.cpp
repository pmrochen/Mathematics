#include "TestUtilities.hpp"
#include <unordered_set>

namespace mathematics::templates {

TEST_SUITE("Matrix3")
{

TEST_CASE_TEMPLATE("construction and element access", T, TEST_REAL_TYPES)
{
	using M = Matrix3<T>;
	using V = Vector3<T>;

	CHECK(M().isZero());
	CHECK(M(Identity()).isIdentity());

	M m(T(1), T(2), T(3), T(4), T(5), T(6), T(7), T(8), T(9));
	CHECK(m.m00 == T(1));
	CHECK(m.m12 == T(6));
	CHECK(m.m21 == T(8));
	CHECK(m[1] == V(T(4), T(5), T(6)));
	CHECK(m[2][0] == T(7));
	CHECK(m.getRow(0) == V(T(1), T(2), T(3)));
	CHECK(m.getColumn(0) == V(T(1), T(4), T(7)));
	CHECK(m.getColumn(2) == V(T(3), T(6), T(9)));

	CHECK(M(V(T(1), T(2), T(3)), V(T(4), T(5), T(6)), V(T(7), T(8), T(9))) == m);
	CHECK(M(std::tuple<V, V, V>(V(T(1), T(2), T(3)), V(T(4), T(5), T(6)), V(T(7), T(8), T(9)))) == m);
	const T a[] = { T(1), T(2), T(3), T(4), T(5), T(6), T(7), T(8), T(9) };
	CHECK(M(a) == m);
	CHECK(M(Matrix2<T>(T(1), T(2), T(3), T(4))) == M(T(1), T(2), T(0), T(3), T(4), T(0), T(0), T(0), T(1)));

	std::tuple<V, V, V> t = m;
	CHECK(std::get<2>(t) == V(T(7), T(8), T(9)));

	m[1][2] = T(0);
	CHECK(m.m12 == T(0));
}

TEST_CASE_TEMPLATE("arithmetic operators", T, TEST_REAL_TYPES)
{
	using M = Matrix3<T>;
	M a(T(1), T(2), T(3), T(4), T(5), T(6), T(7), T(8), T(10));
	M b(T(0), T(1), T(0), T(-1), T(0), T(0), T(0), T(0), T(2));

	CHECK(+a == a);
	CHECK(-a == a*T(-1));
	CHECK(a + b == M(T(1), T(3), T(3), T(3), T(5), T(6), T(7), T(8), T(12)));
	CHECK(a - b == M(T(1), T(1), T(3), T(5), T(5), T(6), T(7), T(8), T(8)));
	CHECK(T(2)*a == a*T(2));
	CHECK(a/T(2) == a*T(0.5));
	CHECK(a*b == M(T(-2), T(1), T(6), T(-5), T(4), T(12), T(-8), T(7), T(20)));
	CHECK(b*a == M(T(4), T(5), T(6), T(-1), T(-2), T(-3), T(14), T(16), T(20)));

	M c = a;
	c += b;
	c -= b;
	CHECK(c == a);
	c *= T(2);
	c /= T(2);
	CHECK(c == a);
	c *= b;
	CHECK(c == a*b);
	CHECK(a != b);
}

TEST_CASE_TEMPLATE("predicates", T, TEST_REAL_TYPES)
{
	using M = Matrix3<T>;
	using V = Vector3<T>;
	CHECK(M::ZERO.isZero());
	CHECK(M::IDENTITY.isIdentity());
	CHECK_FALSE(M::ZERO.isIdentity());
	CHECK((M::ZERO + M(Constants<T>::TOLERANCE*T(0.5), T(0), T(0), T(0), T(0), T(0), T(0), T(0), T(0))).isApproxZero());
	CHECK((M::IDENTITY + M(T(0), T(0), T(0), T(0), T(0), Constants<T>::TOLERANCE*T(0.5), T(0), T(0), T(0))).isApproxIdentity());
	CHECK_FALSE((M::IDENTITY*T(1.1)).isApproxIdentity());
	CHECK(M::IDENTITY.approxEquals(M::IDENTITY*T(1.05), T(0.1)));
	CHECK_FALSE(M::IDENTITY.approxEquals(M::IDENTITY*T(1.05)));

	CHECK(M::makeRotation(V(T(1), T(2), T(3)), T(0.4)).isApproxOrthogonal());
	CHECK_FALSE(M::makeScaling(T(2)).isApproxOrthogonal());
	CHECK(M::makeScaling(T(2)).hasApproxUniformScaling());
	CHECK((M::makeRotation(Axis::Y, T(0.4))*T(3)).hasApproxUniformScaling());
	CHECK_FALSE(M::makeScaling(V(T(1), T(2), T(1))).hasApproxUniformScaling());

	CHECK(M::IDENTITY.isFinite());
	CHECK_FALSE(M(T(1), T(0), T(0), T(0), T(1), T(0), T(0), T(0), std::numeric_limits<T>::quiet_NaN()).isFinite());
}

TEST_CASE_TEMPLATE("scalar properties", T, TEST_REAL_TYPES)
{
	using M = Matrix3<T>;
	using V = Vector3<T>;
	M m(T(1), T(2), T(3), T(4), T(5), T(6), T(7), T(8), T(10));
	CHECK(m.getTrace() == T(16));
	CHECK_APPROX(m.getDeterminant(), T(-3));
	CHECK_FALSE(m.isSingular());
	CHECK(M(T(1), T(2), T(3), T(4), T(5), T(6), T(7), T(8), T(9)).isSingular());
	CHECK_APPROX(M::makeRotation(Axis::Z, T(0.3)).getDeterminant(), T(1));
	CHECK_APPROX(M::makeScalingRotation(V(T(2), T(3), T(4)), V(T(1), T(1), T(0)), T(0.3)).getScaleComponent(), V(T(2), T(3), T(4)));
}

TEST_CASE_TEMPLATE("setters and factories", T, TEST_REAL_TYPES)
{
	using M = Matrix3<T>;
	using V = Vector3<T>;
	M m(T(1), T(2), T(3), T(4), T(5), T(6), T(7), T(8), T(9));
	CHECK(M(m).setZero().isZero());
	CHECK(M(m).setIdentity().isIdentity());
	CHECK(M(m).set(V(T(9), T(8), T(7)), V(T(6), T(5), T(4)), V(T(3), T(2), T(1))) == M(T(9), T(8), T(7), T(6), T(5), T(4), T(3), T(2), T(1)));
	CHECK(M(m).set(T(9), T(8), T(7), T(6), T(5), T(4), T(3), T(2), T(1)) == M(T(9), T(8), T(7), T(6), T(5), T(4), T(3), T(2), T(1)));

	CHECK(M::makeScaling(V(T(2), T(3), T(4))) == M(T(2), T(0), T(0), T(0), T(3), T(0), T(0), T(0), T(4)));
	CHECK(M::makeScaling(T(2)) == M::IDENTITY*T(2));
	CHECK(M::makeShearing(T(1), T(2), T(3), T(4), T(5), T(6)) == M(T(1), T(1), T(2), T(3), T(1), T(4), T(5), T(6), T(1)));
	CHECK(M::makeRotation(Axis::X, T(0)).isIdentity());
	CHECK(M::makeRotation(V::UNIT_X, T(0)).isIdentity());
	CHECK(M::makeRotation(V::ZERO, T(1)).isIdentity());
}

TEST_CASE_TEMPLATE("rotations", T, TEST_REAL_TYPES)
{
	using M = Matrix3<T>;
	using V = Vector3<T>;
	const T halfPi = Constants<T>::HALF_PI;

	// Active right-handed rotations of row vectors.
	CHECK_APPROX(V::UNIT_Y*M::makeRotation(Axis::X, halfPi), V::UNIT_Z);
	CHECK_APPROX(V::UNIT_Z*M::makeRotation(Axis::Y, halfPi), V::UNIT_X);
	CHECK_APPROX(V::UNIT_X*M::makeRotation(Axis::Z, halfPi), V::UNIT_Y);

	// Axis-angle rotation agrees with the principal axis rotations.
	for (Axis axis : { Axis::X, Axis::Y, Axis::Z })
	{
		CAPTURE(int(axis));
		CHECK_APPROX(M::makeRotation(V(axis)*T(3), T(0.8)), M::makeRotation(axis, T(0.8)));
	}

	// Rotation about an arbitrary axis leaves the axis unchanged.
	V axis = normalize(V(T(1), T(-2), T(0.5)));
	M r = M::makeRotation(axis, T(1.1));
	CHECK_APPROX(axis*r, axis);
	CHECK(r.isApproxOrthogonal());
	CHECK_APPROX(r.getDeterminant(), T(1));
	CHECK_APPROX(r*M::makeRotation(axis, T(-1.1)), M::IDENTITY);

	V s(T(2), T(3), T(4));
	CHECK_APPROX(M::makeScalingRotation(s, Axis::Z, T(0.3)), M::makeScaling(s)*M::makeRotation(Axis::Z, T(0.3)));
	CHECK_APPROX(M::makeScalingRotation(s, axis, T(0.3)), M::makeScaling(s)*M::makeRotation(axis, T(0.3)));
}

TEST_CASE_TEMPLATE("look rotations", T, TEST_REAL_TYPES)
{
	using M = Matrix3<T>;
	using V = Vector3<T>;

	SUBCASE("forward")
	{
		for (V forward : { V(T(0), T(0), T(2)), V(T(1), T(0), T(0)), V(T(1), T(2), T(3)), V(T(0), T(5), T(0)), V(T(0), T(-1), T(0)) })
		{
			CAPTURE(forward);
			M m(forward);
			CHECK(m.isApproxOrthogonal());
			CHECK_APPROX(m.getDeterminant(), T(1));
			CHECK_APPROX(m.getRow(2), normalize(forward)); // the third basis vector is the forward direction
		}

		CHECK(M(V::UNIT_Z).isApproxIdentity());
		CHECK(M(V::ZERO).isIdentity());
	}

	SUBCASE("up and forward")
	{
		M m(V::UNIT_Y, V(T(1), T(0), T(1)));
		CHECK(m.isApproxOrthogonal());
		CHECK_APPROX(m.getRow(1), V::UNIT_Y);
		CHECK_APPROX(m.getRow(2), normalize(V(T(1), T(0), T(1))));
		CHECK(M(V::UNIT_Y, V::UNIT_Z).isApproxIdentity());
		CHECK(M(V::ZERO, V::UNIT_Z).isIdentity());
	}
}

TEST_CASE_TEMPLATE("transpose, inverse and adjoint", T, TEST_REAL_TYPES)
{
	using M = Matrix3<T>;
	M m(T(2), T(0), T(1), T(1), T(3), T(2), T(1), T(1), T(2));
	M mt(T(2), T(1), T(1), T(0), T(3), T(1), T(1), T(2), T(2));

	CHECK(transpose(m) == mt);
	CHECK(M(m).transpose() == mt);
	CHECK(M().setTranspose(m) == mt);

	M inv = inverse(m);
	CHECK_APPROX(m*inv, M::IDENTITY);
	CHECK_APPROX(inv*m, M::IDENTITY);
	CHECK_APPROX(M(m).invert(), inv);
	CHECK_APPROX(M().setInverse(m), inv);
	CHECK_APPROX(inverseTranspose(m), transpose(inv));
	CHECK_APPROX(M().setInverseTranspose(m), transpose(inv));

	CHECK_APPROX(adjoint(m)*m, M::IDENTITY*m.getDeterminant());
	CHECK_APPROX(adjoint(m), inv*m.getDeterminant());

	// A permutation matrix has zeros on the diagonal.
	M p(T(0), T(1), T(0), T(0), T(0), T(1), T(1), T(0), T(0));
	CHECK_APPROX(inverse(p), transpose(p));
}

TEST_CASE_TEMPLATE("concatenation and scaling", T, TEST_REAL_TYPES)
{
	using M = Matrix3<T>;
	using V = Vector3<T>;
	M a = M::makeRotation(Axis::X, T(0.5));
	M b(T(1), T(2), T(0), T(-1), T(3), T(1), T(0), T(1), T(2));
	V v(T(0.5), T(-2), T(1));

	CHECK_APPROX(v*M(a).concatenate(b), (v*a)*b);
	CHECK_APPROX(v*M(a).preConcatenate(b), (v*b)*a);
	CHECK_APPROX(concatenate(a, b), a*b);

	V s(T(2), T(3), T(4));
	CHECK_APPROX(v*M(b).preScale(s), (v*s)*b);
	CHECK_APPROX(v*M(b).scale(s), (v*b)*s);
	CHECK_APPROX(M(b).preScale(T(2)), b*T(2));
	CHECK_APPROX(M(b).scale(T(2)), b*T(2));
	CHECK(M(b).negate() == -b);
}

TEST_CASE_TEMPLATE("orthonormalize", T, TEST_REAL_TYPES)
{
	using M = Matrix3<T>;
	M m(T(2), T(0.1), T(0), T(0.3), T(3), T(0.2), T(0.1), T(-0.2), T(1.5));
	m.orthonormalize();
	CHECK(m.isApproxOrthogonal());
	CHECK_APPROX(m.getRow(0), normalize(Vector3<T>(T(2), T(0.1), T(0))));
}

TEST_CASE_TEMPLATE("streaming and hashing", T, TEST_REAL_TYPES)
{
	using M = Matrix3<T>;
	M m(T(1), T(2), T(3), T(4), T(5), T(6), T(7), T(8), T(9.5));
	std::stringstream s;
	s << m;
	M n;
	s >> n;
	CHECK(n == m);

	std::unordered_set<M> set{ m, m, M::IDENTITY };
	CHECK(set.size() == 2);
}

} // TEST_SUITE

} // namespace mathematics::templates
