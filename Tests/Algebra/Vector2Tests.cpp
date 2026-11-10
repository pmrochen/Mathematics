#include "TestUtilities.hpp"
#include <unordered_set>

namespace mathematics::templates {

TEST_SUITE("Vector2")
{

TEST_CASE_TEMPLATE("construction", T, TEST_REAL_TYPES)
{
	using V = Vector2<T>;

	CHECK(V() == V(T(0), T(0)));
	CHECK(V(T(2)) == V(T(2), T(2)));

	V v(T(1), T(2));
	CHECK(v.x == T(1));
	CHECK(v.y == T(2));

	CHECK(V(std::pair<T, T>(T(1), T(2))) == v);
	CHECK(V(std::pair<int, int>(1, 2)) == v);
	CHECK(V(std::tuple<T, T>(T(1), T(2))) == v);
	CHECK(V(std::tuple<int, int>(1, 2)) == v);

	const T a[] = { T(1), T(2) };
	CHECK(V(a) == v);

	CHECK(V(Axis::X) == V::UNIT_X);
	CHECK(V(Axis::Y) == V::UNIT_Y);

	CHECK(V(Vector2<int>(1, 2)) == v);
	CHECK(V(Vector2<double>(1.0, 2.0)) == v);
	CHECK(V(Vector2<float>(1.f, 2.f)) == v);
}

TEST_CASE_TEMPLATE("conversion and element access", T, TEST_REAL_TYPES)
{
	using V = Vector2<T>;
	V v(T(3), T(4));

	std::pair<T, T> p = v;
	CHECK(p == std::pair<T, T>(T(3), T(4)));
	std::tuple<T, T> t = v;
	CHECK(t == std::tuple<T, T>(T(3), T(4)));
	CHECK(std::pair<int, int>(v) == std::pair<int, int>(3, 4));
	CHECK(std::tuple<int, int>(v) == std::tuple<int, int>(3, 4));

	CHECK(v[0] == T(3));
	CHECK(v[1] == T(4));
	v[1] = T(5);
	CHECK(v.y == T(5));

	const T* ptr = static_cast<const T*>(v);
	CHECK(ptr[0] == T(3));
	CHECK(ptr[1] == T(5));

	CHECK(v.template get<0>() == T(3));
	CHECK(get<1>(v) == T(5));
	auto [x, y] = v;
	CHECK(x == T(3));
	CHECK(y == T(5));
	static_assert(std::tuple_size_v<V> == 2);
	static_assert(std::is_same_v<std::tuple_element_t<0, V>, T>);
}

TEST_CASE_TEMPLATE("arithmetic operators", T, TEST_REAL_TYPES)
{
	using V = Vector2<T>;
	V a(T(1), T(2));
	V b(T(3), T(-4));

	CHECK(+a == a);
	CHECK(-a == V(T(-1), T(-2)));
	CHECK(a + b == V(T(4), T(-2)));
	CHECK(a - b == V(T(-2), T(6)));
	CHECK(a*b == V(T(3), T(-8)));
	CHECK(a*T(2) == V(T(2), T(4)));
	CHECK(T(2)*a == V(T(2), T(4)));
	CHECK(b/a == V(T(3), T(-2)));
	CHECK(a/T(2) == V(T(0.5), T(1)));
	CHECK(T(2)/a == V(T(2), T(1)));

	V c = a;
	c += b;
	CHECK(c == V(T(4), T(-2)));
	c -= b;
	CHECK(c == a);
	c *= b;
	CHECK(c == V(T(3), T(-8)));
	c /= b;
	CHECK(c == a);
	c *= T(4);
	CHECK(c == V(T(4), T(8)));
	c /= T(4);
	CHECK(c == a);

	CHECK(a != b);
	CHECK_FALSE(a == b);
}

TEST_CASE_TEMPLATE("predicates", T, TEST_REAL_TYPES)
{
	using V = Vector2<T>;
	CHECK(V::ZERO.isZero());
	CHECK_FALSE(V::UNIT_X.isZero());
	CHECK(V(Constants<T>::TOLERANCE*T(0.5)).isApproxZero());
	CHECK_FALSE(V(T(0), Constants<T>::TOLERANCE*T(2)).isApproxZero());
	CHECK(V(T(1), T(2)).approxEquals(V(T(1), T(2) + Constants<T>::TOLERANCE*T(0.5))));
	CHECK_FALSE(V(T(1), T(2)).approxEquals(V(T(1), T(2.1))));
	CHECK(V(T(1), T(2)).approxEquals(V(T(1.05), T(2)), T(0.1)));

	V a(T(1), T(2));
	CHECK(a.allLessThan(V(T(2), T(3))));
	CHECK_FALSE(a.allLessThan(V(T(2), T(2))));
	CHECK(a.allLessThanEqual(V(T(1), T(2))));
	CHECK(a.allGreaterThan(V(T(0), T(1))));
	CHECK_FALSE(a.allGreaterThan(V(T(1), T(1))));
	CHECK(a.allGreaterThanEqual(V(T(1), T(2))));
	CHECK(a.anyLessThan(V(T(2), T(0))));
	CHECK_FALSE(a.anyLessThan(V(T(1), T(2))));
	CHECK(a.anyLessThanEqual(V(T(1), T(0))));
	CHECK(a.anyGreaterThan(V(T(0), T(5))));
	CHECK_FALSE(a.anyGreaterThan(V(T(1), T(2))));
	CHECK(a.anyGreaterThanEqual(V(T(5), T(2))));

	CHECK(a.isFinite());
	CHECK_FALSE(V::INF.isFinite());
	CHECK_FALSE(V(T(0), std::numeric_limits<T>::quiet_NaN()).isFinite());
}

TEST_CASE_TEMPLATE("magnitude", T, TEST_REAL_TYPES)
{
	using V = Vector2<T>;
	V v(T(3), T(4));
	CHECK_APPROX(v.getMagnitude(), T(5));
	CHECK_APPROX(v.getLength(), T(5));
	CHECK(v.getMagnitudeSquared() == T(25));
	CHECK(v.getLengthSquared() == T(25));

	v.setMagnitude(T(10));
	CHECK_APPROX(v, V(T(6), T(8)));
	v.setLength(T(5));
	CHECK_APPROX(v, V(T(3), T(4)));

	V zero;
	zero.setMagnitude(T(3)); // must not produce NaNs
	CHECK(zero.isZero());
}

TEST_CASE_TEMPLATE("components", T, TEST_REAL_TYPES)
{
	using V = Vector2<T>;
	CHECK(V(T(3), T(-4)).getMajorAxis() == Axis::Y);
	CHECK(V(T(-5), T(4)).getMajorAxis() == Axis::X);
	CHECK(V(T(1), T(-3)).getMajorAxis() == Axis::Y);
	CHECK(V(T(2), T(1)).getMajorAxis() == Axis::X);
	CHECK(V(T(3), T(-4)).getMinComponent() == T(-4));
	CHECK(V(T(3), T(-4)).getMaxComponent() == T(3));
}

TEST_CASE_TEMPLATE("modifiers", T, TEST_REAL_TYPES)
{
	using V = Vector2<T>;
	V v(T(1), T(2));
	CHECK(v.setZero().isZero());
	CHECK(v.set(T(3), T(4)) == V(T(3), T(4)));
	CHECK(v.setMinimum(V(T(1), T(5)), V(T(2), T(0))) == V(T(1), T(0)));
	CHECK(v.setMaximum(V(T(1), T(5)), V(T(2), T(0))) == V(T(2), T(5)));
	CHECK(v.negate() == V(T(-2), T(-5)));

	v.set(T(3), T(4));
	v.normalize();
	CHECK_APPROX(v, V(T(0.6), T(0.8)));
	CHECK_APPROX(v.getMagnitude(), T(1));

	V zero;
	zero.normalize();
	CHECK(zero.isZero());
}

TEST_CASE_TEMPLATE("rotation and transformation", T, TEST_REAL_TYPES)
{
	using V = Vector2<T>;
	const T halfPi = Constants<T>::HALF_PI;

	CHECK_APPROX(V(V::UNIT_X).rotate(halfPi), V::UNIT_Y);
	CHECK_APPROX(V(V::UNIT_Y).rotate(halfPi), V::MINUS_UNIT_X);
	CHECK_APPROX(V(T(1), T(2)).rotate(Constants<T>::PI), V(T(-1), T(-2)));
	CHECK(V(T(1), T(2)).rotate(T(0)) == V(T(1), T(2)));

	// Row vector convention: v' = v*M. A rotation matrix must agree with rotate().
	Matrix2<T> r = Matrix2<T>::makeRotation(T(0.3));
	V v(T(1.5), T(-0.5));
	CHECK_APPROX(v*r, V(v).rotate(T(0.3)));
	CHECK_APPROX(transform(v, r), v*r);
	CHECK_APPROX(V(v).transform(r), v*r);
	V w = v;
	w *= r;
	CHECK_APPROX(w, v*r);

	// Column vector multiplication is the transposed operation.
	CHECK_APPROX(r*v, v*transpose(r));
}

TEST_CASE_TEMPLATE("free functions", T, TEST_REAL_TYPES)
{
	using V = Vector2<T>;
	V a(T(1), T(-2));
	V b(T(4), T(2));

	CHECK(abs(a) == V(T(1), T(2)));
	CHECK(abs(V(T(-3), T(4))) == V(T(3), T(4)));
	CHECK(sum(a) == T(-1));
	CHECK(dot(a, b) == T(0));
	CHECK(dot(a, a) == T(5));
	CHECK(cross(V::UNIT_X, V::UNIT_Y) == T(1));
	CHECK(cross(V::UNIT_Y, V::UNIT_X) == T(-1));
	CHECK(cross(a, b) == T(10));
	CHECK_APPROX(distance(a, b), T(5));
	CHECK(distanceSquared(a, b) == T(25));
	CHECK_APPROX(length(b), std::sqrt(T(20)));
	CHECK(lengthSquared(b) == T(20));
	CHECK(min(a, b) == V(T(1), T(-2)));
	CHECK(max(a, b) == V(T(4), T(2)));
	CHECK(clamp(V(T(-5), T(5)), V(T(-1)), V(T(1))) == V(T(-1), T(1)));
	CHECK(lerp(a, b, T(0)) == a);
	CHECK(lerp(a, b, T(1)) == b);
	CHECK_APPROX(lerp(a, b, T(0.5)), V(T(2.5), T(0)));
	CHECK_APPROX(normalize(V(T(0), T(-3))), V::MINUS_UNIT_Y);
	CHECK_APPROX(normalize(V(T(0), T(-3)) + V(T(0))), V::MINUS_UNIT_Y); // rvalue overload
	CHECK(perpendicular(V(T(1), T(2))) == V(T(-2), T(1)));
	CHECK(dot(perpendicular(a), a) == T(0));
	CHECK_APPROX(project(V(T(2), T(3)), V(T(5), T(0))), V(T(2), T(0)));
	CHECK(project(a, V::ZERO) == V::ZERO);
	CHECK_APPROX(angle(V::UNIT_X, V::UNIT_Y), Constants<T>::HALF_PI);
	CHECK_APPROX(angle(V(T(2), T(0)), V(T(-3), T(0))), Constants<T>::PI);
	CHECK_APPROX(angle(V(T(1), T(1)), V(T(1), T(0))), Constants<T>::PI/T(4));
}

TEST_CASE_TEMPLATE("slerp", T, TEST_REAL_TYPES)
{
	using V = Vector2<T>;
	CHECK_APPROX(slerp(V::UNIT_X, V::UNIT_Y, T(0)), V::UNIT_X);
	CHECK_APPROX(slerp(V::UNIT_X, V::UNIT_Y, T(1)), V::UNIT_Y);
	CHECK_APPROX(slerp(V::UNIT_X, V::UNIT_Y, T(0.5)), V(std::sqrt(T(0.5))));
	CHECK_APPROX(slerp(V::UNIT_X, V::UNIT_Y, T(1)/T(3)), V(std::cos(Constants<T>::PI/T(6)), std::sin(Constants<T>::PI/T(6))));
	CHECK_APPROX(slerp(V::UNIT_X, V::UNIT_X, T(0.5)), V::UNIT_X); // degenerate case
}

TEST_CASE_TEMPLATE("constants", T, TEST_REAL_TYPES)
{
	using V = Vector2<T>;
	CHECK(V::ZERO == V(T(0), T(0)));
	CHECK(V::UNIT_X == V(T(1), T(0)));
	CHECK(V::UNIT_Y == V(T(0), T(1)));
	CHECK(V::MINUS_UNIT_X == V(T(-1), T(0)));
	CHECK(V::MINUS_UNIT_Y == V(T(0), T(-1)));
	CHECK(V::ONE == V(T(1), T(1)));
	CHECK(V::TOLERANCE == V(Constants<T>::TOLERANCE));
	CHECK(V::INF == V(Constants<T>::INF));
	CHECK(V::MINUS_INF == V(Constants<T>::MINUS_INF));
}

TEST_CASE_TEMPLATE("streaming and hashing", T, TEST_REAL_TYPES)
{
	using V = Vector2<T>;
	std::stringstream s;
	s << V(T(1.5), T(-2));
	CHECK(s.str() == "1.5 -2");
	V v;
	s >> v;
	CHECK(v == V(T(1.5), T(-2)));

	std::unordered_set<V> set{ V(T(1), T(2)), V(T(1), T(2)), V(T(2), T(1)) };
	CHECK(set.size() == 2);
	CHECK(std::hash<V>()(V(T(1), T(2))) == std::hash<V>()(V(T(1), T(2))));
}

TEST_CASE("integer vector")
{
	using V = IntVector2;
	V a(1, -2);
	V b(4, 2);
	CHECK(V() == V::ZERO);
	CHECK(V(3) == V(3, 3));
	CHECK(a + b == V(5, 0));
	CHECK(a - b == V(-3, -4));
	CHECK(a*b == V(4, -4));
	CHECK(a*3 == V(3, -6));
	CHECK(3*a == V(3, -6));
	CHECK(b/2 == V(2, 1));
	CHECK(b/V(2, 1) == V(2, 2));
	CHECK(-a == V(-1, 2));
	CHECK(abs(a) == V(1, 2));
	CHECK(sum(b) == 6);
	CHECK(min(a, b) == V(1, -2));
	CHECK(max(a, b) == V(4, 2));
	CHECK(clamp(V(-5, 5), V(-1), V(1)) == V(-1, 1));
	CHECK(a.getMinComponent() == -2);
	CHECK(a.getMaxComponent() == 1);
	CHECK(a.allLessThan(V(2, 0)));
	CHECK(a.anyGreaterThan(V(0, 0)));
	CHECK(V(a).negate() == -a);
	CHECK(V().setMinimum(a, b) == V(1, -2));
	CHECK(V().setMaximum(a, b) == V(4, 2));

	V c = a;
	c += b;
	c -= b;
	c *= 2;
	c /= 2;
	CHECK(c == a);

	std::stringstream s;
	s << a;
	CHECK(s.str() == "1 -2");
	auto [x, y] = a;
	CHECK(x == 1);
	CHECK(y == -2);
	CHECK(Vector2f(a) == Vector2f(1.f, -2.f));
	CHECK(V(Vector2f(1.7f, -2.2f)) == V(1, -2));
}

} // TEST_SUITE

} // namespace mathematics::templates
