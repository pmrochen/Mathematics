#include "TestUtilities.hpp"
#include <unordered_set>

namespace mathematics::templates {

TEST_SUITE("Vector4")
{

TEST_CASE_TEMPLATE("construction", T, TEST_REAL_TYPES)
{
	using V = Vector4<T>;

	CHECK(V() == V(T(0), T(0), T(0), T(0)));
	CHECK(V(T(2)) == V(T(2), T(2), T(2), T(2)));

	V v(T(1), T(2), T(3), T(4));
	CHECK(v.x == T(1));
	CHECK(v.y == T(2));
	CHECK(v.z == T(3));
	CHECK(v.w == T(4));

	// Lower dimensional vectors are promoted to homogeneous points.
	CHECK(V(Vector2<T>(T(1), T(2))) == V(T(1), T(2), T(0), T(1)));
	CHECK(V(Vector3<T>(T(1), T(2), T(3))) == V(T(1), T(2), T(3), T(1)));
	CHECK(V(Vector2<T>(T(1), T(2)), T(3), T(4)) == v);
	CHECK(V(Vector2<T>(T(1), T(2)), Vector2<T>(T(3), T(4))) == v);
	CHECK(V(Vector3<T>(T(1), T(2), T(3)), T(4)) == v);
	CHECK(V(std::tuple<T, T, T, T>(T(1), T(2), T(3), T(4))) == v);

	const T a[] = { T(1), T(2), T(3), T(4) };
	CHECK(V(a) == v);

	CHECK(V(Axis::X) == V::UNIT_X);
	CHECK(V(Axis::Y) == V::UNIT_Y);
	CHECK(V(Axis::Z) == V::UNIT_Z);
	CHECK(V(Axis::W) == V::UNIT_W);

	CHECK(V(Vector4<int>(1, 2, 3, 4)) == v);
	CHECK(V(Vector4<double>(1.0, 2.0, 3.0, 4.0)) == v);
}

TEST_CASE_TEMPLATE("conversion and element access", T, TEST_REAL_TYPES)
{
	using V = Vector4<T>;
	V v(T(3), T(4), T(5), T(6));

	std::tuple<T, T, T, T> t = v;
	CHECK(t == std::tuple<T, T, T, T>(T(3), T(4), T(5), T(6)));

	CHECK(v[3] == T(6));
	v[3] = T(7);
	CHECK(v.w == T(7));
	CHECK(static_cast<const T*>(v)[3] == T(7));

	CHECK(v.xy() == Vector2<T>(T(3), T(4)));
	CHECK(v.xyz() == Vector3<T>(T(3), T(4), T(5)));

	CHECK(v.template get<3>() == T(7));
	CHECK(get<0>(v) == T(3));
	auto [x, y, z, w] = v;
	CHECK(x + y + z + w == T(19));
	static_assert(std::tuple_size_v<V> == 4);
}

TEST_CASE_TEMPLATE("arithmetic operators", T, TEST_REAL_TYPES)
{
	using V = Vector4<T>;
	V a(T(1), T(2), T(4), T(-1));
	V b(T(3), T(-4), T(2), T(2));

	CHECK(+a == a);
	CHECK(-a == V(T(-1), T(-2), T(-4), T(1)));
	CHECK(a + b == V(T(4), T(-2), T(6), T(1)));
	CHECK(a - b == V(T(-2), T(6), T(2), T(-3)));
	CHECK(a*b == V(T(3), T(-8), T(8), T(-2)));
	CHECK(a*T(2) == V(T(2), T(4), T(8), T(-2)));
	CHECK(T(2)*a == V(T(2), T(4), T(8), T(-2)));
	CHECK(b/a == V(T(3), T(-2), T(0.5), T(-2)));
	CHECK(a/T(2) == V(T(0.5), T(1), T(2), T(-0.5)));
	CHECK(T(4)/a == V(T(4), T(2), T(1), T(-4)));

	V c = a;
	c += b;
	c -= b;
	CHECK(c == a);
	c *= b;
	c /= b;
	CHECK(c == a);
	c *= T(4);
	c /= T(4);
	CHECK(c == a);
	CHECK(a != V(T(1), T(2), T(4), T(-2))); // differs only in w
}

TEST_CASE_TEMPLATE("predicates", T, TEST_REAL_TYPES)
{
	using V = Vector4<T>;
	CHECK(V::ZERO.isZero());
	CHECK_FALSE(V::UNIT_W.isZero());
	CHECK(V(Constants<T>::TOLERANCE*T(0.5)).isApproxZero());
	CHECK_FALSE(V(T(0), T(0), T(0), Constants<T>::TOLERANCE*T(2)).isApproxZero());
	CHECK(V(T(1)).approxEquals(V(T(1), T(1), T(1), T(1) + Constants<T>::TOLERANCE*T(0.5))));
	CHECK_FALSE(V(T(1)).approxEquals(V(T(1), T(1), T(1), T(1.1))));
	CHECK(V(T(1)).approxEquals(V(T(1), T(1), T(1), T(1.05)), T(0.1)));

	V a(T(1), T(2), T(3), T(4));
	CHECK(a.allLessThan(V(T(5))));
	CHECK_FALSE(a.allLessThan(V(T(4))));
	CHECK(a.allLessThanEqual(a));
	CHECK(a.allGreaterThan(V(T(0))));
	CHECK(a.allGreaterThanEqual(a));
	CHECK(a.anyLessThan(V(T(0), T(0), T(0), T(5))));
	CHECK_FALSE(a.anyLessThan(a));
	CHECK(a.anyLessThanEqual(V(T(0), T(0), T(0), T(4))));
	CHECK(a.anyGreaterThan(V(T(5), T(5), T(5), T(3))));
	CHECK_FALSE(a.anyGreaterThan(a));
	CHECK(a.anyGreaterThanEqual(V(T(5), T(5), T(5), T(4))));

	CHECK(a.isFinite());
	CHECK_FALSE(V(T(0), T(0), T(0), std::numeric_limits<T>::infinity()).isFinite());
}

TEST_CASE_TEMPLATE("magnitude and modifiers", T, TEST_REAL_TYPES)
{
	using V = Vector4<T>;
	V v(T(1), T(2), T(2), T(4));
	CHECK_APPROX(v.getMagnitude(), T(5));
	CHECK_APPROX(v.getLength(), T(5));
	CHECK(v.getMagnitudeSquared() == T(25));
	CHECK(v.getLengthSquared() == T(25));
	CHECK(v.getMinComponent() == T(1));
	CHECK(v.getMaxComponent() == T(4));

	v.setMagnitude(T(10));
	CHECK_APPROX(v, V(T(2), T(4), T(4), T(8)));
	v.setLength(T(5));
	CHECK_APPROX(v, V(T(1), T(2), T(2), T(4)));
	v.normalize();
	CHECK_APPROX(v, V(T(0.2), T(0.4), T(0.4), T(0.8)));

	CHECK(v.setZero().isZero());
	CHECK(v.set(T(1), T(2), T(3), T(4)) == V(T(1), T(2), T(3), T(4)));
	CHECK(v.negate() == V(T(-1), T(-2), T(-3), T(-4)));
	CHECK(V().setMinimum(V(T(1), T(5), T(0), T(3)), V(T(2), T(0), T(1), T(3))) == V(T(1), T(0), T(0), T(3)));
	CHECK(V().setMaximum(V(T(1), T(5), T(0), T(3)), V(T(2), T(0), T(1), T(3))) == V(T(2), T(5), T(1), T(3)));

	V zero;
	zero.normalize();
	CHECK(zero.isZero());
}

TEST_CASE_TEMPLATE("transformation", T, TEST_REAL_TYPES)
{
	using V = Vector4<T>;
	Matrix4<T> m = Matrix4<T>::makeTranslation(Vector3<T>(T(1), T(2), T(3)));
	V p(T(1), T(1), T(1), T(1));
	V d(T(1), T(1), T(1), T(0));

	CHECK_APPROX(p*m, V(T(2), T(3), T(4), T(1)));
	CHECK_APPROX(d*m, d); // directions are not translated
	CHECK_APPROX(transform(p, m), V(T(2), T(3), T(4), T(1)));
	CHECK_APPROX(V(p).transform(m), V(T(2), T(3), T(4), T(1)));
	V q = p;
	q *= m;
	CHECK_APPROX(q, V(T(2), T(3), T(4), T(1)));
	CHECK_APPROX(m*p, p*transpose(m));
}

TEST_CASE_TEMPLATE("free functions", T, TEST_REAL_TYPES)
{
	using V = Vector4<T>;
	V a(T(1), T(-2), T(2), T(0));
	V b(T(4), T(2), T(0), T(0));

	CHECK(abs(a) == V(T(1), T(2), T(2), T(0)));
	CHECK(sum(a) == T(1));
	CHECK(dot(a, b) == T(0));
	CHECK(dot(a, a) == T(9));
	CHECK_APPROX(distance(a, b), std::sqrt(T(29)));
	CHECK(distanceSquared(a, b) == T(29));
	CHECK_APPROX(length(a), T(3));
	CHECK(lengthSquared(a) == T(9));
	CHECK(min(a, b) == V(T(1), T(-2), T(0), T(0)));
	CHECK(max(a, b) == V(T(4), T(2), T(2), T(0)));
	CHECK(clamp(V(T(-5), T(5), T(0.5), T(2)), V(T(-1)), V(T(1))) == V(T(-1), T(1), T(0.5), T(1)));
	CHECK_APPROX(lerp(a, b, T(0.5)), V(T(2.5), T(0), T(1), T(0)));
	CHECK_APPROX(normalize(V(T(0), T(0), T(0), T(-3))), -V::UNIT_W);
	CHECK_APPROX(normalize(V(T(0), T(0), T(0), T(-3)) + V::ZERO), -V::UNIT_W); // rvalue overload

	CHECK_APPROX(slerp(V::UNIT_X, V::UNIT_W, T(0.5)), V(std::sqrt(T(0.5)), T(0), T(0), std::sqrt(T(0.5))));

	Matrix4<T> t = tensor(V(T(1), T(2), T(3), T(4)), V::UNIT_W);
	CHECK(t.getColumn(3) == V(T(1), T(2), T(3), T(4)));
	CHECK(t.getColumn(0) == V::ZERO);
}

TEST_CASE_TEMPLATE("four dimensional cross product", T, TEST_REAL_TYPES)
{
	using V = Vector4<T>;
	V c = cross(V::UNIT_X, V::UNIT_Y, V::UNIT_Z);
	CHECK_APPROX(std::fabs(c.w), T(1));
	CHECK_APPROX(c.xyz().getMagnitude(), T(0));

	V a(T(1), T(2), T(-1), T(0.5));
	V b(T(0), T(1), T(3), T(-2));
	V d(T(2), T(-1), T(1), T(1));
	V n = cross(a, b, d);
	CHECK_APPROX(dot(n, a), T(0));
	CHECK_APPROX(dot(n, b), T(0));
	CHECK_APPROX(dot(n, d), T(0));
	CHECK_FALSE(n.isApproxZero());
	CHECK(cross(a, a, d).isApproxZero()); // linearly dependent inputs
}

TEST_CASE_TEMPLATE("constants", T, TEST_REAL_TYPES)
{
	using V = Vector4<T>;
	CHECK(V::ZERO == V(T(0)));
	CHECK(V::UNIT_X == V(T(1), T(0), T(0), T(0)));
	CHECK(V::UNIT_Y == V(T(0), T(1), T(0), T(0)));
	CHECK(V::UNIT_Z == V(T(0), T(0), T(1), T(0)));
	CHECK(V::UNIT_W == V(T(0), T(0), T(0), T(1)));
	CHECK(V::ONE == V(T(1)));
	CHECK(V::INF == V(Constants<T>::INF));
}

TEST_CASE_TEMPLATE("streaming and hashing", T, TEST_REAL_TYPES)
{
	using V = Vector4<T>;
	std::stringstream s;
	s << V(T(1.5), T(-2), T(3), T(4));
	CHECK(s.str() == "1.5 -2 3 4");
	V v;
	s >> v;
	CHECK(v == V(T(1.5), T(-2), T(3), T(4)));

	std::unordered_set<V> set{ V(T(1)), V(T(1)), V(T(2)) };
	CHECK(set.size() == 2);
}

TEST_CASE("integer vector")
{
	using V = IntVector4;
	V a(1, -2, 3, -4);
	V b(4, 2, -1, 0);
	CHECK(V() == V::ZERO);
	CHECK(a + b == V(5, 0, 2, -4));
	CHECK(a - b == V(-3, -4, 4, -4));
	CHECK(a*b == V(4, -4, -3, 0));
	CHECK(a*2 == V(2, -4, 6, -8));
	CHECK(-a == V(-1, 2, -3, 4));
	CHECK(abs(a) == V(1, 2, 3, 4));
	CHECK(sum(a) == -2);
	CHECK(min(a, b) == V(1, -2, -1, -4));
	CHECK(max(a, b) == V(4, 2, 3, 0));
	CHECK(a.getMinComponent() == -4);
	CHECK(a.getMaxComponent() == 3);
	CHECK(Vector4f(a) == Vector4f(1.f, -2.f, 3.f, -4.f));
}

} // TEST_SUITE

} // namespace mathematics::templates
