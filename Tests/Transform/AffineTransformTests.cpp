#include "TestUtilities.hpp"
#include <unordered_set>

namespace mathematics::templates {

namespace {

template<typename T>
AffineTransform<T> makeGeneral()
{
	return AffineTransform<T>(Matrix3<T>(T(2), T(0), T(1), T(1), T(3), T(2), T(1), T(1), T(2)), Vector3<T>(T(4), T(-2), T(3)));
}

// The reference result: transform a point using the equivalent homogeneous 4x4 matrix.
template<typename T>
Vector3<T> transformByMatrix(const Vector3<T>& p, const AffineTransform<T>& a)
{
	return (Vector4<T>(p, T(1))*Matrix4<T>(a)).xyz();
}

} // namespace

TEST_SUITE("AffineTransform")
{

TEST_CASE_TEMPLATE("construction and element access", T, TEST_REAL_TYPES)
{
	using A = AffineTransform<T>;
	using V = Vector3<T>;
	using M = Matrix3<T>;

	CHECK(A() == A::ZERO);
	CHECK(A(Identity()).isIdentity());

	A a(T(1), T(2), T(3), T(4), T(5), T(6), T(7), T(8), T(9), T(10), T(11), T(12));
	CHECK(a.m00 == T(1));
	CHECK(a.m12 == T(6));
	CHECK(a.m31 == T(11));
	CHECK(a[3] == V(T(10), T(11), T(12)));
	CHECK(a.getRow(1) == V(T(4), T(5), T(6)));
	CHECK(a.getColumn(0) == Vector4<T>(T(1), T(4), T(7), T(10)));
	CHECK(a.getColumn(2) == Vector4<T>(T(3), T(6), T(9), T(12)));
	CHECK(a.getOrigin() == V(T(10), T(11), T(12)));
	CHECK(a.getBasis() == M(T(1), T(2), T(3), T(4), T(5), T(6), T(7), T(8), T(9)));

	CHECK(A(V(T(1), T(2), T(3)), V(T(4), T(5), T(6)), V(T(7), T(8), T(9)), V(T(10), T(11), T(12))) == a);
	CHECK(A(std::tuple<V, V, V, V>(a[0], a[1], a[2], a[3])) == a);
	CHECK(A(M(T(1), T(2), T(3), T(4), T(5), T(6), T(7), T(8), T(9)), V(T(10), T(11), T(12))) == a);
	CHECK(A(M(T(1), T(2), T(3), T(4), T(5), T(6), T(7), T(8), T(9))).getOrigin() == V::ZERO);
	CHECK(A(Matrix2<T>(T(1), T(2), T(3), T(4))) == A(M(Matrix2<T>(T(1), T(2), T(3), T(4))), V::ZERO));
	const T p[] = { T(1), T(2), T(3), T(4), T(5), T(6), T(7), T(8), T(9), T(10), T(11), T(12) };
	CHECK(A(p) == a);

	std::tuple<V, V, V, V> t = a;
	CHECK(std::get<3>(t) == a.getOrigin());

	a.setOrigin(V(T(-1)));
	CHECK(a.m30 == T(-1));
	a.setBasis(M::IDENTITY);
	CHECK(a.getBasis().isIdentity());
	CHECK(a[3] == V(T(-1)));
}

TEST_CASE_TEMPLATE("look transforms", T, TEST_REAL_TYPES)
{
	using A = AffineTransform<T>;
	using V = Vector3<T>;
	using M = Matrix3<T>;
	V origin(T(1), T(2), T(3));

	for (V forward : { V(T(1), T(2), T(3)), V(T(0), T(1), T(0)), V(T(0), T(0), T(-1)) })
	{
		CAPTURE(forward);
		A a(forward, origin);
		CHECK_APPROX(a.getBasis(), M(forward));
		CHECK(a.getOrigin() == origin);
	}

	A b(V::UNIT_Y, V(T(1), T(0), T(1)), origin);
	CHECK_APPROX(b.getBasis(), M(V::UNIT_Y, V(T(1), T(0), T(1))));
	CHECK(b.getOrigin() == origin);

	CHECK(A(V::ZERO, origin).getBasis().isIdentity());
	CHECK(A(V::ZERO, V::UNIT_Z, origin).getBasis().isIdentity());
}

TEST_CASE_TEMPLATE("predicates", T, TEST_REAL_TYPES)
{
	using A = AffineTransform<T>;
	using V = Vector3<T>;
	CHECK(A::IDENTITY.isIdentity());
	CHECK_FALSE(A::makeTranslation(V::UNIT_X).isIdentity());
	CHECK(A::makeTranslation(V(Constants<T>::TOLERANCE*T(0.5))).isApproxIdentity());
	CHECK_FALSE(A::makeTranslation(V(T(0.1))).isApproxIdentity());
	CHECK(A::makeRotationTranslation(V(T(1), T(1), T(0)), T(0.4), V(T(3))).isApproxOrthogonal());
	CHECK_FALSE(A::makeScaling(T(2)).isApproxOrthogonal());
	CHECK(A::makeScaling(T(2)).hasApproxUniformScaling());
	CHECK_FALSE(A::makeScaling(V(T(1), T(2), T(3))).hasApproxUniformScaling());
	CHECK(makeGeneral<T>().approxEquals(makeGeneral<T>()));
	CHECK(makeGeneral<T>().approxEquals(makeGeneral<T>().translate(V(T(0.05))), T(0.1)));
	CHECK_FALSE(makeGeneral<T>().approxEquals(makeGeneral<T>().translate(V(T(0.05)))));
	CHECK(makeGeneral<T>().isFinite());
	CHECK_FALSE(A::makeTranslation(V::INF).isFinite());
	CHECK(makeGeneral<T>() != A::IDENTITY);
	CHECK_APPROX(A::makeScalingRotation(V(T(2), T(3), T(4)), Axis::X, T(0.3)).getScaleComponent(), V(T(2), T(3), T(4)));
}

TEST_CASE_TEMPLATE("setters and factories", T, TEST_REAL_TYPES)
{
	using A = AffineTransform<T>;
	using V = Vector3<T>;
	using M = Matrix3<T>;
	using M4 = Matrix4<T>;
	V t(T(10), T(20), T(30));
	V s(T(2), T(3), T(4));

	CHECK(A(makeGeneral<T>()).setZero() == A::ZERO);
	CHECK(A(makeGeneral<T>()).setIdentity().isIdentity());
	CHECK(A().set(M::IDENTITY, t) == A::makeTranslation(t));
	CHECK(A().set(V::UNIT_X, V::UNIT_Y, V::UNIT_Z, t) == A::makeTranslation(t));
	CHECK(A().set(T(1), T(0), T(0), T(0), T(1), T(0), T(0), T(0), T(1), T(10), T(20), T(30)) == A::makeTranslation(t));

	CHECK(M4(A::makeTranslation(t)) == M4::makeTranslation(t));
	CHECK(M4(A::makeScaling(s)) == M4::makeScaling(s));
	CHECK(M4(A::makeScaling(T(2))) == M4::makeScaling(V(T(2))));
	CHECK(M4(A::makeScalingTranslation(s, t)) == M4::makeScalingTranslation(s, t));
	CHECK(M4(A::makeShearing(T(1), T(2), T(3), T(4), T(5), T(6))) == M4::makeShearing(T(1), T(2), T(3), T(4), T(5), T(6)));

	for (Axis axis : { Axis::X, Axis::Y, Axis::Z })
	{
		CAPTURE(int(axis));
		CHECK_APPROX(M4(A::makeRotation(axis, T(0.6))), M4::makeRotation(axis, T(0.6)));
		CHECK_APPROX(M4(A::makeRotation(V(axis), T(0.6))), M4::makeRotation(V(axis), T(0.6)));
		CHECK_APPROX(M4(A::makeRotationTranslation(axis, T(0.6), t)), M4::makeRotationTranslation(axis, T(0.6), t));
		CHECK_APPROX(M4(A::makeRotationTranslation(V(axis), T(0.6), t)), M4::makeRotationTranslation(V(axis), T(0.6), t));
		CHECK_APPROX(M4(A::makeScalingRotation(s, axis, T(0.6))), M4::makeScalingRotation(s, axis, T(0.6)));
		CHECK_APPROX(M4(A::makeScalingRotation(s, V(axis), T(0.6))), M4::makeScalingRotation(s, V(axis), T(0.6)));
		CHECK_APPROX(M4(A::makeScalingRotationTranslation(s, axis, T(0.6), t)), M4::makeScalingRotationTranslation(s, axis, T(0.6), t));
		CHECK_APPROX(M4(A::makeScalingRotationTranslation(s, V(axis), T(0.6), t)), M4::makeScalingRotationTranslation(s, V(axis), T(0.6), t));
	}

	Plane<T> plane(normalize(V(T(1), T(1), T(0))), T(-1));
	CHECK_APPROX(M4(A::makeReflection(plane)), M4::makeReflection(plane));
}

TEST_CASE_TEMPLATE("point transformation", T, TEST_REAL_TYPES)
{
	using V = Vector3<T>;
	AffineTransform<T> a = makeGeneral<T>();
	V p(T(1), T(-1), T(2));
	CHECK_APPROX(transform(p, a), transformByMatrix(p, a));
	CHECK_APPROX(V(p).transform(a), transformByMatrix(p, a));
	CHECK_APPROX(transform(p, a), p*a.getBasis() + a.getOrigin());
}

TEST_CASE_TEMPLATE("inversion", T, TEST_REAL_TYPES)
{
	using A = AffineTransform<T>;
	using V = Vector3<T>;
	A a = makeGeneral<T>();
	V p(T(1), T(-1), T(2));

	A inv = inverse(a);
	CHECK_APPROX(transform(transform(p, a), inv), p);
	CHECK_APPROX(concatenate(a, inv), A::IDENTITY);
	CHECK_APPROX(A().setInverse(a), inv);
	CHECK_APPROX(A(a).invert(), inv);
	CHECK_APPROX(inverse(A(a)), inv); // rvalue overload
	CHECK_APPROX(Matrix4<T>(inv), inverse(Matrix4<T>(a)));

	A r = A::makeRotationTranslation(V(T(1), T(2), T(3)), T(0.8), V(T(5), T(-1), T(2)));
	A rInv = inverse(r);
	CHECK_APPROX(inverseOrthogonal(r), rInv);
	CHECK_APPROX(inverseOrthogonal(A(r)), rInv); // rvalue overload
	CHECK_APPROX(A().setInverseOrthogonal(r), rInv);
	CHECK_APPROX(A(r).invertOrthogonal(), rInv);
}

TEST_CASE_TEMPLATE("concatenation, translation and scaling", T, TEST_REAL_TYPES)
{
	using A = AffineTransform<T>;
	using V = Vector3<T>;
	using M4 = Matrix4<T>;
	A a = makeGeneral<T>();
	A b = A::makeRotationTranslation(Axis::Y, T(0.7), V(T(1), T(2), T(3)));
	A c = A::makeScalingTranslation(V(T(2), T(1), T(0.5)), V(T(-1), T(0), T(1)));
	A d = A::makeRotation(Axis::Z, T(-0.4));
	V p(T(1), T(-1), T(2));

	// concatenate appends a transformation, preConcatenate prepends one.
	CHECK_APPROX(transform(p, concatenate(a, b)), transform(transform(p, a), b));
	CHECK_APPROX(transform(p, concatenate(a, b, c)), transform(transform(transform(p, a), b), c));
	CHECK_APPROX(transform(p, concatenate(a, b, c, d)), transform(transform(transform(transform(p, a), b), c), d));
	CHECK_APPROX(A(a).concatenate(b), concatenate(a, b));
	CHECK_APPROX(A(a).preConcatenate(b), concatenate(b, a));
	CHECK_APPROX(M4(concatenate(a, b)), M4(a)*M4(b));

	V t(T(1), T(2), T(3));
	CHECK_APPROX(transform(p, A(a).translate(t)), transform(p, a) + t);

	V s(T(2), T(3), T(4));
	CHECK_APPROX(transform(p, A(a).preScale(s)), transform(p*s, a));
	CHECK_APPROX(transform(p, A(a).preScale(T(2))), transform(p*T(2), a));
	CHECK_APPROX(transform(p, A(a).scale(s)), transform(p, a)*s);
	CHECK_APPROX(transform(p, A(a).scale(T(2))), transform(p, a)*T(2));
}

TEST_CASE_TEMPLATE("orthonormalize", T, TEST_REAL_TYPES)
{
	using A = AffineTransform<T>;
	A a(Matrix3<T>(T(2), T(0.1), T(0), T(0.3), T(3), T(0.2), T(0.1), T(-0.2), T(1.5)), Vector3<T>(T(1), T(2), T(3)));
	a.orthonormalize();
	CHECK(a.isApproxOrthogonal());
	CHECK(a.getOrigin() == Vector3<T>(T(1), T(2), T(3)));
}

TEST_CASE_TEMPLATE("streaming and hashing", T, TEST_REAL_TYPES)
{
	using A = AffineTransform<T>;
	A a = makeGeneral<T>();
	std::stringstream s;
	s << a;
	A b;
	s >> b;
	CHECK(b == a);

	std::unordered_set<A> set{ a, a, A::IDENTITY };
	CHECK(set.size() == 2);
}

} // TEST_SUITE

} // namespace mathematics::templates
