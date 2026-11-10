#include "TestUtilities.hpp"
#include <unordered_set>

namespace mathematics::templates {

TEST_SUITE("Vector3")
{

TEST_CASE_TEMPLATE("construction", T, TEST_REAL_TYPES)
{
	using V = Vector3<T>;

	CHECK(V() == V(T(0), T(0), T(0)));
	CHECK(V(T(2)) == V(T(2), T(2), T(2)));

	V v(T(1), T(2), T(3));
	CHECK(v.x == T(1));
	CHECK(v.y == T(2));
	CHECK(v.z == T(3));

	CHECK(V(Vector2<T>(T(1), T(2))) == V(T(1), T(2), T(0)));
	CHECK(V(Vector2<T>(T(1), T(2)), T(3)) == v);
	CHECK(V(std::tuple<T, T, T>(T(1), T(2), T(3))) == v);
	CHECK(V(std::tuple<int, int, int>(1, 2, 3)) == v);

	const T a[] = { T(1), T(2), T(3) };
	CHECK(V(a) == v);

	CHECK(V(Axis::X) == V::UNIT_X);
	CHECK(V(Axis::Y) == V::UNIT_Y);
	CHECK(V(Axis::Z) == V::UNIT_Z);

	CHECK(V(Vector3<int>(1, 2, 3)) == v);
	CHECK(V(Vector3<double>(1.0, 2.0, 3.0)) == v);
	CHECK(V(Vector3<float>(1.f, 2.f, 3.f)) == v);
}

TEST_CASE_TEMPLATE("conversion and element access", T, TEST_REAL_TYPES)
{
	using V = Vector3<T>;
	V v(T(3), T(4), T(5));

	std::tuple<T, T, T> t = v;
	CHECK(t == std::tuple<T, T, T>(T(3), T(4), T(5)));
	CHECK(std::tuple<int, int, int>(v) == std::tuple<int, int, int>(3, 4, 5));

	CHECK(v[0] == T(3));
	CHECK(v[1] == T(4));
	CHECK(v[2] == T(5));
	v[2] = T(6);
	CHECK(v.z == T(6));

	const T* ptr = static_cast<const T*>(v);
	CHECK(ptr[2] == T(6));

	CHECK(v.xy() == Vector2<T>(T(3), T(4)));
	CHECK(v.xz() == Vector2<T>(T(3), T(6)));
	CHECK(v.zy() == Vector2<T>(T(6), T(4)));

	CHECK(v.template get<0>() == T(3));
	CHECK(get<2>(v) == T(6));
	auto [x, y, z] = v;
	CHECK(x == T(3));
	CHECK(y == T(4));
	CHECK(z == T(6));
	static_assert(std::tuple_size_v<V> == 3);
}

TEST_CASE_TEMPLATE("arithmetic operators", T, TEST_REAL_TYPES)
{
	using V = Vector3<T>;
	V a(T(1), T(2), T(4));
	V b(T(3), T(-4), T(2));

	CHECK(+a == a);
	CHECK(-a == V(T(-1), T(-2), T(-4)));
	CHECK(a + b == V(T(4), T(-2), T(6)));
	CHECK(a - b == V(T(-2), T(6), T(2)));
	CHECK(a*b == V(T(3), T(-8), T(8)));
	CHECK(a*T(2) == V(T(2), T(4), T(8)));
	CHECK(T(2)*a == V(T(2), T(4), T(8)));
	CHECK(b/a == V(T(3), T(-2), T(0.5)));
	CHECK(a/T(2) == V(T(0.5), T(1), T(2)));
	CHECK(T(4)/a == V(T(4), T(2), T(1)));

	V c = a;
	c += b;
	CHECK(c == V(T(4), T(-2), T(6)));
	c -= b;
	CHECK(c == a);
	c *= b;
	CHECK(c == V(T(3), T(-8), T(8)));
	c /= b;
	CHECK(c == a);
	c *= T(4);
	CHECK(c == V(T(4), T(8), T(16)));
	c /= T(4);
	CHECK(c == a);

	CHECK(a != b);
	CHECK(a != V(T(1), T(2), T(5))); // differs only in z
}

TEST_CASE_TEMPLATE("predicates", T, TEST_REAL_TYPES)
{
	using V = Vector3<T>;
	CHECK(V::ZERO.isZero());
	CHECK_FALSE(V::UNIT_Z.isZero());
	CHECK(V(Constants<T>::TOLERANCE*T(0.5)).isApproxZero());
	CHECK_FALSE(V(T(0), T(0), Constants<T>::TOLERANCE*T(2)).isApproxZero());
	CHECK(V(T(1), T(2), T(3)).approxEquals(V(T(1), T(2), T(3) + Constants<T>::TOLERANCE*T(0.5))));
	CHECK_FALSE(V(T(1), T(2), T(3)).approxEquals(V(T(1), T(2), T(3.1))));
	CHECK(V(T(1), T(2), T(3)).approxEquals(V(T(1), T(2), T(3.05)), T(0.1)));

	V a(T(1), T(2), T(3));
	CHECK(a.allLessThan(V(T(2), T(3), T(4))));
	CHECK_FALSE(a.allLessThan(V(T(2), T(3), T(3))));
	CHECK(a.allLessThanEqual(a));
	CHECK(a.allGreaterThan(V(T(0))));
	CHECK_FALSE(a.allGreaterThan(V(T(0), T(0), T(3))));
	CHECK(a.allGreaterThanEqual(a));
	CHECK(a.anyLessThan(V(T(0), T(0), T(4))));
	CHECK_FALSE(a.anyLessThan(a));
	CHECK(a.anyLessThanEqual(V(T(0), T(0), T(3))));
	CHECK(a.anyGreaterThan(V(T(5), T(5), T(2))));
	CHECK_FALSE(a.anyGreaterThan(a));
	CHECK(a.anyGreaterThanEqual(V(T(5), T(5), T(3))));

	CHECK(a.isFinite());
	CHECK_FALSE(V::MINUS_INF.isFinite());
	CHECK_FALSE(V(T(0), T(0), std::numeric_limits<T>::quiet_NaN()).isFinite());
}

TEST_CASE_TEMPLATE("magnitude and components", T, TEST_REAL_TYPES)
{
	using V = Vector3<T>;
	V v(T(2), T(3), T(6));
	CHECK_APPROX(v.getMagnitude(), T(7));
	CHECK_APPROX(v.getLength(), T(7));
	CHECK(v.getMagnitudeSquared() == T(49));
	CHECK(v.getLengthSquared() == T(49));

	v.setMagnitude(T(14));
	CHECK_APPROX(v, V(T(4), T(6), T(12)));
	v.setLength(T(7));
	CHECK_APPROX(v, V(T(2), T(3), T(6)));

	CHECK(V(T(1), T(-3), T(2)).getMajorAxis() == Axis::Y);
	CHECK(V(T(1), T(2), T(-3)).getMajorAxis() == Axis::Z);
	CHECK(V(T(-4), T(2), T(3)).getMajorAxis() == Axis::X);
	CHECK(V(T(1), T(2), T(-3)).getMinComponent() == T(-3));
	CHECK(V(T(1), T(2), T(-3)).getMaxComponent() == T(2));
}

TEST_CASE_TEMPLATE("modifiers", T, TEST_REAL_TYPES)
{
	using V = Vector3<T>;
	V v(T(1), T(2), T(3));
	CHECK(v.setZero().isZero());
	CHECK(v.set(T(3), T(4), T(5)) == V(T(3), T(4), T(5)));
	CHECK(v.setMinimum(V(T(1), T(5), T(0)), V(T(2), T(0), T(1))) == V(T(1), T(0), T(0)));
	CHECK(v.setMaximum(V(T(1), T(5), T(0)), V(T(2), T(0), T(1))) == V(T(2), T(5), T(1)));
	CHECK(v.negate() == V(T(-2), T(-5), T(-1)));

	v.set(T(2), T(3), T(6));
	v.normalize();
	CHECK_APPROX(v, V(T(2), T(3), T(6))/T(7));

	V zero;
	zero.normalize();
	CHECK(zero.isZero());
}

TEST_CASE_TEMPLATE("rotation", T, TEST_REAL_TYPES)
{
	using V = Vector3<T>;
	const T halfPi = Constants<T>::HALF_PI;

	// Right-handed rotations.
	CHECK_APPROX(V(V::UNIT_Y).rotate(Axis::X, halfPi), V::UNIT_Z);
	CHECK_APPROX(V(V::UNIT_Z).rotate(Axis::Y, halfPi), V::UNIT_X);
	CHECK_APPROX(V(V::UNIT_X).rotate(Axis::Z, halfPi), V::UNIT_Y);
	CHECK(V(T(1), T(2), T(3)).rotate(Axis::X, T(0)) == V(T(1), T(2), T(3)));

	V v(T(0.5), T(-1.5), T(2));
	for (Axis axis : { Axis::X, Axis::Y, Axis::Z })
	{
		CAPTURE(int(axis));
		V expected = v*Matrix3<T>::makeRotation(axis, T(0.7));
		CHECK_APPROX(V(v).rotate(axis, T(0.7)), expected);
		CHECK_APPROX(V(v).rotate(Quaternion<T>::fromAxisAngle(V(axis), T(0.7))), expected);
		CHECK_APPROX(rotate(v, Quaternion<T>::fromAxisAngle(V(axis), T(0.7))), expected);
	}
}

TEST_CASE_TEMPLATE("transformation", T, TEST_REAL_TYPES)
{
	using V = Vector3<T>;
	Matrix3<T> m(T(1), T(2), T(3), T(4), T(5), T(6), T(7), T(8), T(10));
	V v(T(1), T(-1), T(2));

	V expected(T(1) - T(4) + T(14), T(2) - T(5) + T(16), T(3) - T(6) + T(20)); // row vector times matrix
	CHECK_APPROX(v*m, expected);
	CHECK_APPROX(transform(v, m), expected);
	CHECK_APPROX(V(v).transform(m), expected);
	V w = v;
	w *= m;
	CHECK_APPROX(w, expected);
	CHECK_APPROX(m*v, v*transpose(m));

	AffineTransform<T> a(m, V(T(10), T(20), T(30)));
	CHECK_APPROX(transform(v, a), expected + V(T(10), T(20), T(30)));
	CHECK_APPROX(V(v).transform(a), expected + V(T(10), T(20), T(30)));
}

TEST_CASE_TEMPLATE("free functions", T, TEST_REAL_TYPES)
{
	using V = Vector3<T>;
	V a(T(1), T(-2), T(2));
	V b(T(4), T(2), T(0));

	CHECK(abs(a) == V(T(1), T(2), T(2)));
	CHECK(abs(V(T(-1), T(-2), T(-3))) == V(T(1), T(2), T(3)));
	CHECK(sum(a) == T(1));
	CHECK(dot(a, b) == T(0));
	CHECK(dot(a, a) == T(9));
	CHECK(cross(V::UNIT_X, V::UNIT_Y) == V::UNIT_Z);
	CHECK(cross(V::UNIT_Y, V::UNIT_Z) == V::UNIT_X);
	CHECK(cross(V::UNIT_Z, V::UNIT_X) == V::UNIT_Y);
	CHECK(cross(a, b) == V(T(-4), T(8), T(10)));
	CHECK_APPROX(distance(a, b), std::sqrt(T(29)));
	CHECK(distanceSquared(a, b) == T(29));
	CHECK_APPROX(length(a), T(3));
	CHECK(lengthSquared(a) == T(9));
	CHECK(min(a, b) == V(T(1), T(-2), T(0)));
	CHECK(max(a, b) == V(T(4), T(2), T(2)));
	CHECK(clamp(V(T(-5), T(5), T(0.5)), V(T(-1)), V(T(1))) == V(T(-1), T(1), T(0.5)));
	CHECK(lerp(a, b, T(0)) == a);
	CHECK(lerp(a, b, T(1)) == b);
	CHECK_APPROX(lerp(a, b, T(0.5)), V(T(2.5), T(0), T(1)));
	CHECK_APPROX(normalize(V(T(0), T(0), T(-3))), V::MINUS_UNIT_Z);
	CHECK_APPROX(normalize(V(T(0), T(0), T(-3)) + V::ZERO), V::MINUS_UNIT_Z); // rvalue overload
	CHECK_APPROX(project(V(T(2), T(3), T(4)), V(T(0), T(5), T(0))), V(T(0), T(3), T(0)));
	CHECK(project(a, V::ZERO) == V::ZERO);
	CHECK_APPROX(angle(V::UNIT_X, V::UNIT_Z), Constants<T>::HALF_PI);
	CHECK_APPROX(angle(V(T(2), T(0), T(0)), V(T(-3), T(0), T(0))), Constants<T>::PI);
	CHECK_APPROX_TOL(angle(a, a), T(0), T(1e-3)); // acos is ill-conditioned near 1

	Matrix3<T> t = tensor(V(T(1), T(2), T(3)), V(T(4), T(5), T(6)));
	CHECK(t == Matrix3<T>(T(4), T(5), T(6), T(8), T(10), T(12), T(12), T(15), T(18)));
}

TEST_CASE_TEMPLATE("perpendicular", T, TEST_REAL_TYPES)
{
	using V = Vector3<T>;
	for (V v : { V::UNIT_X, V::UNIT_Y, V::UNIT_Z, V(T(1), T(2), T(3)), V(T(-3), T(0.5), T(0.1)), V(T(0.2), T(-4), T(1)) })
	{
		CAPTURE(v);
		V p = perpendicular(v);
		CHECK_APPROX(dot(p, v), T(0));
		CHECK_APPROX(p.getMagnitude(), v.getMagnitude());
	}
}

TEST_CASE_TEMPLATE("reflection and refraction", T, TEST_REAL_TYPES)
{
	using V = Vector3<T>;
	V i = normalize(V(T(1), T(-1), T(0)));
	CHECK_APPROX(reflect(i, V::UNIT_Y), normalize(V(T(1), T(1), T(0))));
	CHECK_APPROX(reflect(V::MINUS_UNIT_Y, V::UNIT_Y), V::UNIT_Y);

	// No bending with equal indices.
	CHECK_APPROX(refract(i, V::UNIT_Y, T(1)), i);
	// Snell's law: sin(theta_t) = eta*sin(theta_i)
	V r = refract(i, V::UNIT_Y, T(0.5));
	CHECK_APPROX(r.getMagnitude(), T(1));
	CHECK_APPROX(r.x, T(0.5)*std::sqrt(T(0.5)));
	CHECK(r.y < T(0));
	// Total internal reflection yields zero.
	CHECK(refract(normalize(V(T(1), T(-0.1), T(0))), V::UNIT_Y, T(1.5)).isZero());
}

TEST_CASE_TEMPLATE("slerp", T, TEST_REAL_TYPES)
{
	using V = Vector3<T>;
	CHECK_APPROX(slerp(V::UNIT_X, V::UNIT_Z, T(0)), V::UNIT_X);
	CHECK_APPROX(slerp(V::UNIT_X, V::UNIT_Z, T(1)), V::UNIT_Z);
	CHECK_APPROX(slerp(V::UNIT_X, V::UNIT_Z, T(0.5)), V(std::sqrt(T(0.5)), T(0), std::sqrt(T(0.5))));
	V s = slerp(V::UNIT_Y, V::UNIT_Z, T(0.25));
	CHECK_APPROX(s.getMagnitude(), T(1));
	CHECK_APPROX(angle(V::UNIT_Y, s), Constants<T>::HALF_PI*T(0.25));
	CHECK_APPROX(slerp(V::UNIT_X, V::UNIT_X, T(0.5)), V::UNIT_X); // degenerate case
}

TEST_CASE_TEMPLATE("constants", T, TEST_REAL_TYPES)
{
	using V = Vector3<T>;
	CHECK(V::ZERO == V(T(0)));
	CHECK(V::UNIT_X == V(T(1), T(0), T(0)));
	CHECK(V::UNIT_Y == V(T(0), T(1), T(0)));
	CHECK(V::UNIT_Z == V(T(0), T(0), T(1)));
	CHECK(V::MINUS_UNIT_X == -V::UNIT_X);
	CHECK(V::MINUS_UNIT_Y == -V::UNIT_Y);
	CHECK(V::MINUS_UNIT_Z == -V::UNIT_Z);
	CHECK(V::ONE == V(T(1)));
	CHECK(V::TOLERANCE == V(Constants<T>::TOLERANCE));
	CHECK(V::INF == V(Constants<T>::INF));
	CHECK(V::MINUS_INF == V(Constants<T>::MINUS_INF));
}

TEST_CASE_TEMPLATE("streaming and hashing", T, TEST_REAL_TYPES)
{
	using V = Vector3<T>;
	std::stringstream s;
	s << V(T(1.5), T(-2), T(3));
	CHECK(s.str() == "1.5 -2 3");
	V v;
	s >> v;
	CHECK(v == V(T(1.5), T(-2), T(3)));

	std::unordered_set<V> set{ V(T(1), T(2), T(3)), V(T(1), T(2), T(3)), V(T(3), T(2), T(1)) };
	CHECK(set.size() == 2);
}

TEST_CASE("integer vector")
{
	using V = IntVector3;
	V a(1, -2, 3);
	V b(4, 2, -1);
	CHECK(V() == V::ZERO);
	CHECK(V(3) == V(3, 3, 3));
	CHECK(V(IntVector2(1, 2), 3) == V(1, 2, 3));
	CHECK(a + b == V(5, 0, 2));
	CHECK(a - b == V(-3, -4, 4));
	CHECK(a*b == V(4, -4, -3));
	CHECK(a*2 == V(2, -4, 6));
	CHECK(2*a == V(2, -4, 6));
	CHECK(V(4, 2, -6)/2 == V(2, 1, -3));
	CHECK(-a == V(-1, 2, -3));
	CHECK(abs(a) == V(1, 2, 3));
	CHECK(sum(a) == 2);
	CHECK(min(a, b) == V(1, -2, -1));
	CHECK(max(a, b) == V(4, 2, 3));
	CHECK(a.getMinComponent() == -2);
	CHECK(a.getMaxComponent() == 3);
	CHECK(V().setMinimum(a, b) == V(1, -2, -1));
	CHECK(V().setMaximum(a, b) == V(4, 2, 3));
	auto [x, y, z] = a;
	CHECK(x + y + z == 2);
	CHECK(Vector3f(a) == Vector3f(1.f, -2.f, 3.f));
}

} // TEST_SUITE

} // namespace mathematics::templates
