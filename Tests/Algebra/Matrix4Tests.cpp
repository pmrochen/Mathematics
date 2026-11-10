#include "TestUtilities.hpp"
#include <unordered_set>

namespace mathematics::templates {

namespace {

template<typename T>
Matrix4<T> makeSample()
{
	return Matrix4<T>(T(2), T(0), T(1), T(0),
		T(1), T(3), T(2), T(0),
		T(1), T(1), T(2), T(0),
		T(4), T(-2), T(3), T(1));
}

} // namespace

TEST_SUITE("Matrix4")
{

TEST_CASE_TEMPLATE("construction and element access", T, TEST_REAL_TYPES)
{
	using M = Matrix4<T>;
	using V = Vector4<T>;

	CHECK(M().isZero());
	CHECK(M(Identity()).isIdentity());

	M m(T(1), T(2), T(3), T(4), T(5), T(6), T(7), T(8), T(9), T(10), T(11), T(12), T(13), T(14), T(15), T(16));
	CHECK(m.m00 == T(1));
	CHECK(m.m13 == T(8));
	CHECK(m.m32 == T(15));
	CHECK(m[1] == V(T(5), T(6), T(7), T(8)));
	CHECK(m[3][0] == T(13));
	CHECK(m.getRow(2) == V(T(9), T(10), T(11), T(12)));
	CHECK(m.getColumn(1) == V(T(2), T(6), T(10), T(14)));
	CHECK(m.getColumn(3) == V(T(4), T(8), T(12), T(16)));

	CHECK(M(V(T(1), T(2), T(3), T(4)), V(T(5), T(6), T(7), T(8)), V(T(9), T(10), T(11), T(12)), V(T(13), T(14), T(15), T(16))) == m);
	CHECK(M(std::tuple<V, V, V, V>(m[0], m[1], m[2], m[3])) == m);
	const T a[] = { T(1), T(2), T(3), T(4), T(5), T(6), T(7), T(8), T(9), T(10), T(11), T(12), T(13), T(14), T(15), T(16) };
	CHECK(M(a) == m);

	std::tuple<V, V, V, V> t = m;
	CHECK(std::get<3>(t) == V(T(13), T(14), T(15), T(16)));

	CHECK(M(Matrix2<T>(T(1), T(2), T(3), T(4))) == M(T(1), T(2), T(0), T(0), T(3), T(4), T(0), T(0), T(0), T(0), T(1), T(0), T(0), T(0), T(0), T(1)));
	CHECK(M(Matrix3<T>(T(1), T(2), T(3), T(4), T(5), T(6), T(7), T(8), T(9))) ==
		M(T(1), T(2), T(3), T(0), T(4), T(5), T(6), T(0), T(7), T(8), T(9), T(0), T(0), T(0), T(0), T(1)));
	AffineTransform<T> affine(Matrix3<T>(T(1), T(2), T(3), T(4), T(5), T(6), T(7), T(8), T(9)), Vector3<T>(T(10), T(11), T(12)));
	CHECK(M(affine) == M(T(1), T(2), T(3), T(0), T(4), T(5), T(6), T(0), T(7), T(8), T(9), T(0), T(10), T(11), T(12), T(1)));
}

TEST_CASE_TEMPLATE("arithmetic operators", T, TEST_REAL_TYPES)
{
	using M = Matrix4<T>;
	M a = makeSample<T>();
	M b = M::makeTranslation(Vector3<T>(T(1), T(2), T(3)));

	CHECK(+a == a);
	CHECK(-a == a*T(-1));
	CHECK(a + b - b == a);
	CHECK((a + b).m30 == T(5));
	CHECK(T(2)*a == a*T(2));
	CHECK(a/T(2) == a*T(0.5));
	CHECK(M::IDENTITY*a == a);
	CHECK(a*M::IDENTITY == a);

	M ab = a*b;
	for (int i = 0; i < 4; i++)
	{
		for (int j = 0; j < 4; j++)
		{
			CAPTURE(i);
			CAPTURE(j);
			CHECK_APPROX(ab[i][j], dot(a.getRow(i), b.getColumn(j)));
		}
	}

	M c = a;
	c += b;
	c -= b;
	CHECK(c == a);
	c *= T(2);
	c /= T(2);
	CHECK(c == a);
	c *= b;
	CHECK(c == ab);
	CHECK(a != b);
}

TEST_CASE_TEMPLATE("predicates and scalar properties", T, TEST_REAL_TYPES)
{
	using M = Matrix4<T>;
	CHECK(M::ZERO.isZero());
	CHECK(M::IDENTITY.isIdentity());
	CHECK_FALSE(M::ZERO.isIdentity());
	CHECK((M::IDENTITY*T(1) + M::ZERO).isApproxIdentity());
	CHECK_FALSE((M::IDENTITY*T(1.1)).isApproxIdentity());
	CHECK((M::ZERO*T(2)).isApproxZero());
	CHECK(M::IDENTITY.approxEquals(M::IDENTITY*T(1.05), T(0.1)));
	CHECK_FALSE(M::IDENTITY.approxEquals(M::IDENTITY*T(1.05)));
	CHECK(M::IDENTITY.isFinite());
	CHECK_FALSE((M::IDENTITY*std::numeric_limits<T>::infinity()).isFinite());

	M m = makeSample<T>();
	CHECK(m.getTrace() == T(8));
	CHECK_APPROX(m.getDeterminant(), T(6)); // equals the determinant of the upper 3x3 block
	CHECK_FALSE(m.isSingular());
	CHECK(M(T(1), T(2), T(3), T(4), T(5), T(6), T(7), T(8), T(9), T(10), T(11), T(12), T(13), T(14), T(15), T(16)).isSingular());
	CHECK_APPROX(M::makeRotation(Axis::Y, T(0.7)).getDeterminant(), T(1));
	CHECK_APPROX(M::makeScaling(Vector3<T>(T(2), T(3), T(4))).getDeterminant(), T(24));
}

TEST_CASE_TEMPLATE("setters and factories", T, TEST_REAL_TYPES)
{
	using M = Matrix4<T>;
	using V3 = Vector3<T>;
	using V4 = Vector4<T>;
	M m = makeSample<T>();
	CHECK(M(m).setZero().isZero());
	CHECK(M(m).setIdentity().isIdentity());
	CHECK(M(m).set(m[3], m[2], m[1], m[0]).getRow(0) == m[3]);
	CHECK(M().set(T(1), T(2), T(3), T(4), T(5), T(6), T(7), T(8), T(9), T(10), T(11), T(12), T(13), T(14), T(15), T(16)).m23 == T(12));

	V4 p(T(1), T(2), T(3), T(1));
	V3 t(T(10), T(20), T(30));
	V3 s(T(2), T(3), T(4));

	CHECK_APPROX(p*M::makeTranslation(t), V4(T(11), T(22), T(33), T(1)));
	CHECK_APPROX(p*M::makeScaling(s), V4(T(2), T(6), T(12), T(1)));
	CHECK_APPROX(p*M::makeScalingTranslation(s, t), V4(T(12), T(26), T(42), T(1)));
	CHECK(M::makeShearing(T(1), T(2), T(3), T(4), T(5), T(6)) == M(Matrix3<T>::makeShearing(T(1), T(2), T(3), T(4), T(5), T(6))));

	for (Axis axis : { Axis::X, Axis::Y, Axis::Z })
	{
		CAPTURE(int(axis));
		Matrix3<T> r = Matrix3<T>::makeRotation(axis, T(0.6));
		CHECK_APPROX(M::makeRotation(axis, T(0.6)), M(r));
		CHECK_APPROX(M::makeRotation(V3(axis), T(0.6)), M(r));
		CHECK_APPROX(M::makeRotationTranslation(axis, T(0.6), t), M(AffineTransform<T>(r, t)));
		CHECK_APPROX(M::makeRotationTranslation(V3(axis), T(0.6), t), M(AffineTransform<T>(r, t)));
		CHECK_APPROX(M::makeScalingRotation(s, axis, T(0.6)), M(Matrix3<T>::makeScaling(s)*r));
		CHECK_APPROX(M::makeScalingRotation(s, V3(axis), T(0.6)), M(Matrix3<T>::makeScaling(s)*r));
		CHECK_APPROX(M::makeScalingRotationTranslation(s, axis, T(0.6), t), M(AffineTransform<T>(Matrix3<T>::makeScaling(s)*r, t)));
		CHECK_APPROX(M::makeScalingRotationTranslation(s, V3(axis), T(0.6), t), M(AffineTransform<T>(Matrix3<T>::makeScaling(s)*r, t)));
	}
}

TEST_CASE_TEMPLATE("reflection", T, TEST_REAL_TYPES)
{
	using M = Matrix4<T>;
	using V4 = Vector4<T>;
	Plane<T> plane(Vector3<T>::UNIT_Y, T(-2)); // y = 2
	M r = M::makeReflection(plane);
	CHECK_APPROX(V4(T(1), T(5), T(3), T(1))*r, V4(T(1), T(-1), T(3), T(1)));
	CHECK_APPROX(V4(T(1), T(2), T(3), T(1))*r, V4(T(1), T(2), T(3), T(1)));
	CHECK_APPROX(r*r, M::IDENTITY);
}

TEST_CASE_TEMPLATE("transpose, inverse and adjoint", T, TEST_REAL_TYPES)
{
	using M = Matrix4<T>;
	M m = makeSample<T>();
	m.m03 = T(1); // make it a general (projective) matrix

	M mt = transpose(m);
	CHECK(mt.getRow(0) == m.getColumn(0));
	CHECK(mt.getRow(3) == m.getColumn(3));
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

	SUBCASE("zero leading element")
	{
		M p(T(0), T(1), T(0), T(0),
			T(0), T(0), T(1), T(0),
			T(1), T(0), T(0), T(0),
			T(0), T(0), T(0), T(1));
		CHECK_APPROX(inverse(p), transpose(p));

		M r = M::makeRotationTranslation(Axis::Z, Constants<T>::HALF_PI, Vector3<T>(T(1), T(2), T(3)));
		CHECK_APPROX(r*inverse(r), M::IDENTITY);
	}
}

TEST_CASE_TEMPLATE("concatenation, translation and scaling", T, TEST_REAL_TYPES)
{
	using M = Matrix4<T>;
	using V3 = Vector3<T>;
	using V4 = Vector4<T>;
	M a = M::makeRotation(Axis::X, T(0.5));
	M b = makeSample<T>();
	V4 v(T(0.5), T(-2), T(1), T(1));

	CHECK_APPROX(v*M(a).concatenate(b), (v*a)*b);
	CHECK_APPROX(v*M(a).preConcatenate(b), (v*b)*a);
	CHECK_APPROX(concatenate(a, b), a*b);

	V3 t(T(1), T(2), T(3));
	CHECK_APPROX(v*M(b).translate(t), v*b + V4(t, T(0)));

	V3 s(T(2), T(3), T(4));
	CHECK_APPROX(v*M(b).preScale(s), (v*V4(s, T(1)))*b);
	CHECK_APPROX(v*M(b).scale(s), (v*b)*V4(s, T(1)));
	CHECK(M(b).negate() == -b);
}

TEST_CASE_TEMPLATE("streaming and hashing", T, TEST_REAL_TYPES)
{
	using M = Matrix4<T>;
	M m = makeSample<T>();
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
