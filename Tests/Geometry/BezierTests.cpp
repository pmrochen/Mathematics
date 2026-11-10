#include "TestUtilities.hpp"
#include <unordered_set>
#include <vector>
#include <iterator>
#include <Mathematics/Geometry/Bezier.inl>
#include <Mathematics/Geometry/Hermite.inl>

namespace mathematics::templates {

TEST_SUITE("Bezier")
{

TEST_CASE_TEMPLATE("basis functions", T, TEST_REAL_TYPES)
{
	for (T t : { T(0), T(0.25), T(0.5), T(0.8), T(1) })
	{
		CAPTURE(t);
		Vector4<T> b = bezier::basis<Vector4<T>>(t);
		CHECK_APPROX(sum(b), T(1)); // partition of unity
		CHECK_APPROX(b, Vector4<T>(std::pow(T(1) - t, T(3)), T(3)*t*std::pow(T(1) - t, T(2)), T(3)*t*t*(T(1) - t), t*t*t));
		CHECK_APPROX(sum(bezier::derivativeBasis<Vector4<T>>(t)), T(0));

		Vector4<T> h = hermite::basis<Vector4<T>>(t);
		CHECK_APPROX(h.x + h.w, T(1)); // the position weights sum to one
		CHECK_APPROX(h, Vector4<T>(T(2)*t*t*t - T(3)*t*t + T(1), t*t*t - T(2)*t*t + t, t*t*t - t*t, T(-2)*t*t*t + T(3)*t*t));
	}

	CHECK_APPROX(bezier::evaluate(T(0), T(1), T(2), T(3), T(0.5)), T(1.5)); // evenly spaced control values are linear
	CHECK_APPROX(bezier::evaluate(T(0), T(1), T(1), T(0), T(0.5)), T(0.75));
	CHECK_APPROX(bezier::evaluate(T(2), T(5), T(-1), T(4), T(0)), T(2));
	CHECK_APPROX(bezier::evaluate(T(2), T(5), T(-1), T(4), T(1)), T(4));
}

} // TEST_SUITE

TEST_SUITE("BezierCurve2")
{

TEST_CASE_TEMPLATE("construction and control points", T, TEST_REAL_TYPES)
{
	using B = BezierCurve2<T>;
	using V = Vector2<T>;
	V p0(T(0), T(0)), p1(T(1), T(2)), p2(T(3), T(2)), p3(T(4), T(0));
	B b(p0, p1, p2, p3);
	CHECK(B().isZero());
	CHECK(B::ZERO.isZero());
	CHECK(B(std::tuple<V, V, V, V>(p0, p1, p2, p3)) == b);
	const V points[] = { p0, p1, p2, p3 };
	CHECK(B(points) == b);
	CHECK(b.getControlPoints() == std::tuple<V, V, V, V>(p0, p1, p2, p3));
	std::vector<V> copied;
	b.copyControlPoints(std::back_inserter(copied));
	CHECK(copied == std::vector<V>{ p0, p1, p2, p3 });
	CHECK(b.getControlPoint(2) == p2);
	CHECK(b.getControlPoint(4) == V::ZERO);

	B c = b;
	c.setControlPoint(3, V(T(5), T(0)));
	CHECK(c.points[3] == V(T(5), T(0)));
	CHECK_THROWS_AS(c.setControlPoint(-1, V::ZERO), std::out_of_range);
	CHECK(B(b).setZero().isZero());
	CHECK(B().set(p0, p1, p2, p3) == b);
	CHECK(B(V(Constants<T>::TOLERANCE*T(0.5)), V::ZERO, V::ZERO, V::ZERO).isApproxZero());
	CHECK(b.isFinite());
	CHECK_FALSE(B(p0, p1, p2, V::INF).isFinite());
	CHECK(b.approxEquals(B(p0, p1, p2, p3 + V(T(0.05))), T(0.1)));
	CHECK_FALSE(b.approxEquals(B(p0, p1, p2, p3 + V(T(0.05)))));
	CHECK(b.approxEquals(b));
	CHECK(b != c);
}

TEST_CASE_TEMPLATE("transformation", T, TEST_REAL_TYPES)
{
	using B = BezierCurve2<T>;
	using V = Vector2<T>;
	B b(V(T(0), T(0)), V(T(1), T(2)), V(T(3), T(2)), V(T(4), T(0)));
	CHECK(B(b).translate(V::ONE) == B(V(T(1), T(1)), V(T(2), T(3)), V(T(4), T(3)), V(T(5), T(1))));
	Matrix2<T> r = Matrix2<T>::makeRotation(T(0.5));
	B t = B(b).transform(r);
	CHECK_APPROX(t, B(b.points[0]*r, b.points[1]*r, b.points[2]*r, b.points[3]*r));
	CHECK_APPROX(t.evaluate(T(0.3)), b.evaluate(T(0.3))*r); // affine invariance
}

TEST_CASE_TEMPLATE("evaluation", T, TEST_REAL_TYPES)
{
	using B = BezierCurve2<T>;
	using V = Vector2<T>;
	B b(V(T(0), T(0)), V(T(1), T(2)), V(T(3), T(2)), V(T(4), T(0)));
	CHECK_APPROX(b.evaluate(T(0)), b.points[0]);
	CHECK_APPROX(b.evaluate(T(1)), b.points[3]);
	CHECK_APPROX(b.evaluate(T(0.5)), (b.points[0] + T(3)*b.points[1] + T(3)*b.points[2] + b.points[3])/T(8));
	CHECK_APPROX(b.template evaluate<0>(T(0.5)), b.evaluate(T(0.5)).x);
	CHECK_APPROX(b.template evaluate<1>(T(0.5)), b.evaluate(T(0.5)).y);

	CHECK_APPROX(b.calculateDerivative(T(0)), T(3)*(b.points[1] - b.points[0]));
	CHECK_APPROX(b.calculateDerivative(T(1)), T(3)*(b.points[3] - b.points[2]));
	CHECK_APPROX(b.calculateTangent(T(0)), normalize(V(T(1), T(2))));
	CHECK_APPROX(b.calculateSpeed(T(1)), T(3)*std::sqrt(T(5)));

	// Numerical derivative
	const T h = std::is_same_v<T, float> ? T(1e-2) : T(1e-5);
	CHECK_APPROX_TOL(b.calculateDerivative(T(0.4)), (b.evaluate(T(0.4) + h) - b.evaluate(T(0.4) - h))/(T(2)*h), T(1e-2));
}

TEST_CASE_TEMPLATE("arc length", T, TEST_REAL_TYPES)
{
	using B = BezierCurve2<T>;
	using V = Vector2<T>;

	// Evenly spaced collinear control points: uniform speed
	B line(V(T(0), T(0)), V(T(1), T(1)), V(T(2), T(2)), V(T(3), T(3)));
	T length = T(3)*std::sqrt(T(2));
	CHECK_APPROX(line.calculateLength(), length);
	CHECK_APPROX(line.calculateLength(T(0.5)), length*T(0.5));
	CHECK(line.calculateLength(T(-1)) == T(0));
	CHECK_APPROX(line.calculateLength(T(2)), length); // clamped
	REQUIRE(line.calculateTime(length*T(0.25)).has_value());
	CHECK_APPROX(*line.calculateTime(length*T(0.25)), T(0.25));
	CHECK(*line.calculateTime(T(-1)) == T(0));
	CHECK(*line.calculateTime(length*T(2)) == T(1));

	// A curved segment: compare to a dense polyline approximation
	B curve(V(T(0), T(0)), V(T(1), T(2)), V(T(3), T(2)), V(T(4), T(0)));
	T polyline = T(0);
	for (int i = 0; i < 2000; i++)
		polyline += distance(curve.evaluate(T(i)/T(2000)), curve.evaluate(T(i + 1)/T(2000)));
	CHECK_APPROX_TOL(curve.template calculateLength<10>(), polyline, T(1e-4));
	auto time = curve.template calculateTime<10, 64>(polyline*T(0.3), T(1e-5));
	REQUIRE(time.has_value());
	CHECK_APPROX_TOL(curve.template calculateLength<10>(*time), polyline*T(0.3), T(1e-4));
}

TEST_CASE_TEMPLATE("streaming and hashing", T, TEST_REAL_TYPES)
{
	using B = BezierCurve2<T>;
	using V = Vector2<T>;
	B b(V(T(0), T(0)), V(T(1), T(2)), V(T(3), T(2)), V(T(4), T(0.5)));
	std::stringstream s;
	s << b;
	CHECK(s.str() == "0 0 1 2 3 2 4 0.5");
	B read;
	s >> read;
	CHECK(read == b);
	std::unordered_set<B> set{ b, b, B() };
	CHECK(set.size() == 2);
}

} // TEST_SUITE

TEST_SUITE("BezierCurve3")
{

TEST_CASE_TEMPLATE("construction and control points", T, TEST_REAL_TYPES)
{
	using B = BezierCurve3<T>;
	using V = Vector3<T>;
	V p0(T(0), T(0), T(0)), p1(T(1), T(2), T(0)), p2(T(3), T(2), T(1)), p3(T(4), T(0), T(1));
	B b(p0, p1, p2, p3);
	CHECK(B().isZero());
	CHECK(B::ZERO.isZero());
	CHECK(B(std::tuple<V, V, V, V>(p0, p1, p2, p3)) == b);
	const V points[] = { p0, p1, p2, p3 };
	CHECK(B(points) == b);
	CHECK(b.getControlPoints() == std::tuple<V, V, V, V>(p0, p1, p2, p3));
	std::vector<V> copied;
	b.copyControlPoints(std::back_inserter(copied));
	CHECK(copied == std::vector<V>{ p0, p1, p2, p3 });
	CHECK(b.getControlPoint(1) == p1);
	CHECK(b.getControlPoint(-1) == V::ZERO);

	B c = b;
	c.setControlPoint(0, V::ONE);
	CHECK(c.points[0] == V::ONE);
	CHECK_THROWS_AS(c.setControlPoint(4, V::ZERO), std::out_of_range);
	CHECK(B(b).setZero().isZero());
	CHECK(B().set(p0, p1, p2, p3) == b);
	CHECK(B().isApproxZero());
	CHECK(b.isFinite());
	CHECK_FALSE(B(p0, p1, V::INF, p3).isFinite());
	CHECK(b.approxEquals(B(p0, p1, p2, p3 + V(T(0.05))), T(0.1)));
	CHECK_FALSE(b.approxEquals(B(p0, p1, p2, p3 + V(T(0.05)))));
	CHECK(b != c);
}

TEST_CASE_TEMPLATE("transformation and evaluation", T, TEST_REAL_TYPES)
{
	using B = BezierCurve3<T>;
	using V = Vector3<T>;
	B b(V(T(0), T(0), T(0)), V(T(1), T(2), T(0)), V(T(3), T(2), T(1)), V(T(4), T(0), T(1)));
	CHECK(B(b).translate(V::ONE).points[3] == V(T(5), T(1), T(2)));
	Matrix3<T> m = Matrix3<T>::makeRotation(Axis::X, T(0.5));
	CHECK_APPROX(B(b).transform(m).evaluate(T(0.3)), b.evaluate(T(0.3))*m);
	AffineTransform<T> a = AffineTransform<T>::makeRotationTranslation(Axis::Y, T(0.5), V::ONE);
	CHECK_APPROX(B(b).transform(a).evaluate(T(0.7)), transform(b.evaluate(T(0.7)), a));

	CHECK_APPROX(b.evaluate(T(0)), b.points[0]);
	CHECK_APPROX(b.evaluate(T(1)), b.points[3]);
	CHECK_APPROX(b.evaluate(T(0.5)), (b.points[0] + T(3)*b.points[1] + T(3)*b.points[2] + b.points[3])/T(8));
	CHECK_APPROX(b.template evaluate<2>(T(0.5)), b.evaluate(T(0.5)).z);
	CHECK_APPROX(b.calculateDerivative(T(0)), T(3)*(b.points[1] - b.points[0]));
	CHECK_APPROX(b.calculateDerivative(T(1)), T(3)*(b.points[3] - b.points[2]));
	CHECK_APPROX(b.calculateTangent(T(1)), normalize(b.points[3] - b.points[2]));
	CHECK_APPROX(b.calculateSpeed(T(0)), T(3)*std::sqrt(T(5)));

	B line(V(T(0)), V(T(1), T(2), T(2)), V(T(2), T(4), T(4)), V(T(3), T(6), T(6)));
	CHECK_APPROX(line.calculateLength(), T(9));
	REQUIRE(line.calculateTime(T(3)).has_value());
	CHECK_APPROX(*line.calculateTime(T(3)), T(1)/T(3));
}

TEST_CASE_TEMPLATE("streaming and hashing", T, TEST_REAL_TYPES)
{
	using B = BezierCurve3<T>;
	using V = Vector3<T>;
	B b(V(T(0)), V(T(1)), V(T(2)), V(T(3.5)));
	std::stringstream s;
	s << b;
	B read;
	s >> read;
	CHECK(read == b);
	std::unordered_set<B> set{ b, b, B() };
	CHECK(set.size() == 2);
}

} // TEST_SUITE

} // namespace mathematics::templates
