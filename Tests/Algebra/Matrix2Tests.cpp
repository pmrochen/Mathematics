#include "TestUtilities.hpp"
#include <unordered_set>

namespace mathematics::templates {

TEST_SUITE("Matrix2")
{

TEST_CASE_TEMPLATE("construction and element access", T, TEST_REAL_TYPES)
{
	using M = Matrix2<T>;
	using V = Vector2<T>;

	CHECK(M().isZero());
	CHECK(M(Identity()).isIdentity());

	M m(T(1), T(2), T(3), T(4));
	CHECK(m.m00 == T(1));
	CHECK(m.m01 == T(2));
	CHECK(m.m10 == T(3));
	CHECK(m.m11 == T(4));
	CHECK(m[0] == V(T(1), T(2)));
	CHECK(m[1][0] == T(3));
	CHECK(m.getRow(1) == V(T(3), T(4)));
	CHECK(m.getColumn(0) == V(T(1), T(3)));
	CHECK(m.getColumn(1) == V(T(2), T(4)));

	CHECK(M(V(T(1), T(2)), V(T(3), T(4))) == m);
	CHECK(M(std::pair<V, V>(V(T(1), T(2)), V(T(3), T(4)))) == m);
	CHECK(M(std::tuple<V, V>(V(T(1), T(2)), V(T(3), T(4)))) == m);
	const T a[] = { T(1), T(2), T(3), T(4) };
	CHECK(M(a) == m);

	std::tuple<V, V> t = m;
	CHECK(std::get<1>(t) == V(T(3), T(4)));

	m[0][1] = T(5);
	CHECK(m.m01 == T(5));
}

TEST_CASE_TEMPLATE("arithmetic operators", T, TEST_REAL_TYPES)
{
	using M = Matrix2<T>;
	M a(T(1), T(2), T(3), T(4));
	M b(T(0), T(1), T(-1), T(2));

	CHECK(+a == a);
	CHECK(-a == M(T(-1), T(-2), T(-3), T(-4)));
	CHECK(a + b == M(T(1), T(3), T(2), T(6)));
	CHECK(a - b == M(T(1), T(1), T(4), T(2)));
	CHECK(a*T(2) == M(T(2), T(4), T(6), T(8)));
	CHECK(T(2)*a == M(T(2), T(4), T(6), T(8)));
	CHECK(a/T(2) == M(T(0.5), T(1), T(1.5), T(2)));
	CHECK(a*b == M(T(-2), T(5), T(-4), T(11)));
	CHECK(b*a == M(T(3), T(4), T(5), T(6)));

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
	using M = Matrix2<T>;
	CHECK(M::ZERO.isZero());
	CHECK(M::IDENTITY.isIdentity());
	CHECK_FALSE(M::ZERO.isIdentity());
	CHECK(M(Constants<T>::TOLERANCE*T(0.5), T(0), T(0), T(0)).isApproxZero());
	CHECK(M(T(1), Constants<T>::TOLERANCE*T(0.5), T(0), T(1)).isApproxIdentity());
	CHECK_FALSE(M(T(1), T(0.1), T(0), T(1)).isApproxIdentity());
	CHECK(M(T(1), T(2), T(3), T(4)).approxEquals(M(T(1), T(2), T(3), T(4.05)), T(0.1)));
	CHECK_FALSE(M(T(1), T(2), T(3), T(4)).approxEquals(M(T(1), T(2), T(3), T(4.05))));

	CHECK(M::makeRotation(T(0.4)).isApproxOrthogonal());
	CHECK_FALSE(M::makeScaling(T(2)).isApproxOrthogonal());
	CHECK(M::makeScaling(T(2)).hasApproxUniformScaling());
	CHECK((M::makeRotation(T(0.4))*T(3)).hasApproxUniformScaling());
	CHECK_FALSE(M::makeScaling(Vector2<T>(T(1), T(2))).hasApproxUniformScaling());

	CHECK(M::IDENTITY.isFinite());
	CHECK_FALSE(M(T(1), std::numeric_limits<T>::infinity(), T(0), T(1)).isFinite());
}

TEST_CASE_TEMPLATE("scalar properties", T, TEST_REAL_TYPES)
{
	using M = Matrix2<T>;
	M m(T(1), T(2), T(3), T(4));
	CHECK(m.getTrace() == T(5));
	CHECK(m.getDeterminant() == T(-2));
	CHECK_FALSE(m.isSingular());
	CHECK(M(T(1), T(2), T(2), T(4)).isSingular());
	CHECK_APPROX(M::makeScalingRotation(Vector2<T>(T(2), T(3)), T(0.3)).getScaleComponent(), Vector2<T>(T(2), T(3)));
}

TEST_CASE_TEMPLATE("setters and factories", T, TEST_REAL_TYPES)
{
	using M = Matrix2<T>;
	using V = Vector2<T>;
	M m(T(1), T(2), T(3), T(4));
	CHECK(M(m).setZero().isZero());
	CHECK(M(m).setIdentity().isIdentity());
	CHECK(M(m).set(V(T(5), T(6)), V(T(7), T(8))) == M(T(5), T(6), T(7), T(8)));
	CHECK(M(m).set(T(5), T(6), T(7), T(8)) == M(T(5), T(6), T(7), T(8)));

	CHECK(M::makeScaling(V(T(2), T(3))) == M(T(2), T(0), T(0), T(3)));
	CHECK(M::makeScaling(T(2)) == M(T(2), T(0), T(0), T(2)));
	CHECK(M::makeShearing(T(2), T(3)) == M(T(1), T(2), T(3), T(1)));
	CHECK(M::makeRotation(T(0)).isIdentity());

	// Row vectors: UNIT_X*R(90 deg) = UNIT_Y
	CHECK_APPROX(V::UNIT_X*M::makeRotation(Constants<T>::HALF_PI), V::UNIT_Y);
	CHECK_APPROX(V(T(2), T(1))*M::makeScalingRotation(V(T(2), T(3)), Constants<T>::HALF_PI), V(T(-3), T(4)));
}

TEST_CASE_TEMPLATE("transpose, inverse and adjoint", T, TEST_REAL_TYPES)
{
	using M = Matrix2<T>;
	M m(T(1), T(2), T(3), T(4));
	M mt(T(1), T(3), T(2), T(4));
	M mi(T(-2), T(1), T(1.5), T(-0.5));

	CHECK(transpose(m) == mt);
	CHECK(M(m).transpose() == mt);
	CHECK(M().setTranspose(m) == mt);

	CHECK_APPROX(inverse(m), mi);
	CHECK_APPROX(M(m).invert(), mi);
	CHECK_APPROX(M().setInverse(m), mi);
	CHECK_APPROX(m*inverse(m), M::IDENTITY);
	CHECK_APPROX(inverseTranspose(m), transpose(mi));
	CHECK_APPROX(M().setInverseTranspose(m), transpose(mi));

	CHECK_APPROX(adjoint(m)*m, M::IDENTITY*m.getDeterminant());
	CHECK_APPROX(adjoint(m), inverse(m)*m.getDeterminant());
}

TEST_CASE_TEMPLATE("concatenation and scaling", T, TEST_REAL_TYPES)
{
	using M = Matrix2<T>;
	using V = Vector2<T>;
	M a = M::makeRotation(T(0.5));
	M b(T(1), T(2), T(-1), T(3));
	V v(T(0.5), T(-2));

	// concatenate appends a transformation, preConcatenate prepends one.
	CHECK_APPROX(v*M(a).concatenate(b), (v*a)*b);
	CHECK_APPROX(v*M(a).preConcatenate(b), (v*b)*a);
	CHECK_APPROX(concatenate(a, b), a*b);

	V s(T(2), T(3));
	CHECK_APPROX(v*M(b).preScale(s), (v*s)*b);
	CHECK_APPROX(v*M(b).scale(s), (v*b)*s);
	CHECK_APPROX(M(b).preScale(T(2)), b*T(2));
	CHECK_APPROX(M(b).scale(T(2)), b*T(2));
	CHECK(M(b).negate() == -b);
}

TEST_CASE_TEMPLATE("orthonormalize", T, TEST_REAL_TYPES)
{
	using M = Matrix2<T>;
	M m(T(2), T(0.1), T(0.3), T(3));
	m.orthonormalize();
	CHECK(m.isApproxOrthogonal());
	CHECK_APPROX(m.getRow(0), Vector2<T>(T(2), T(0.1))/std::sqrt(T(4.01)));
}

TEST_CASE_TEMPLATE("streaming and hashing", T, TEST_REAL_TYPES)
{
	using M = Matrix2<T>;
	M m(T(1), T(2), T(3), T(4.5));
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
