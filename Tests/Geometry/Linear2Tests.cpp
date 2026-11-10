#include "TestUtilities.hpp"
#include <unordered_set>
#include <vector>
#include <iterator>

namespace mathematics::templates {

TEST_SUITE("Line2")
{

TEST_CASE_TEMPLATE("construction and properties", T, TEST_REAL_TYPES)
{
	using L = Line2<T>;
	using V = Vector2<T>;

	L line(V(T(1), T(2)), V(T(3), T(4)));
	CHECK(line.origin == V(T(1), T(2)));
	CHECK(line.direction == V(T(3), T(4)));
	CHECK(line.getOrigin() == V(T(1), T(2)));
	CHECK(line.getDirection() == V(T(3), T(4)));
	CHECK(L(Ray2<T>(V(T(1), T(2)), V(T(3), T(4)))) == line);
	CHECK(line.asRay() == Ray2<T>(V(T(1), T(2)), V(T(3), T(4))));

	L angled(V(T(1), T(1)), Constants<T>::PI/T(4));
	CHECK_APPROX(angled.direction, V(std::sqrt(T(0.5))));
	CHECK_APPROX(angled.getInclinationAngle(), Constants<T>::PI/T(4));
	CHECK_APPROX(angled.getSlope(), T(1));
	CHECK_APPROX(line.getSlope(), T(4)/T(3));

	line.setOrigin(V(T(0)));
	line.setDirection(V::UNIT_Y);
	CHECK(line == L(V::ZERO, V::UNIT_Y));
	CHECK(L(line).set(V::ONE, V::UNIT_X) == L(V::ONE, V::UNIT_X));
	CHECK(L(line).setZero().isZero());
	CHECK(L(V(Constants<T>::TOLERANCE*T(0.5)), V::ZERO).isApproxZero());
	CHECK(line.approxEquals(L(V::ZERO, V(T(0), T(1) + Constants<T>::TOLERANCE*T(0.5)))));
	CHECK(line.approxEquals(L(V::ZERO, V(T(0), T(1.05))), T(0.1)));
	CHECK(line != L(V::ZERO, V::UNIT_X));
}

TEST_CASE_TEMPLATE("transformation and evaluation", T, TEST_REAL_TYPES)
{
	using L = Line2<T>;
	using V = Vector2<T>;
	L line(V(T(1), T(2)), V(T(3), T(4)));
	CHECK(line.evaluate(T(0)) == line.origin);
	CHECK(line.evaluate(T(2)) == V(T(7), T(10)));
	CHECK(L(line).translate(V(T(1), T(-1))) == L(V(T(2), T(1)), V(T(3), T(4))));
	CHECK_APPROX(L(line).normalize().direction, V(T(0.6), T(0.8)));
	CHECK_APPROX(normalize(line).direction, V(T(0.6), T(0.8)));
	CHECK_APPROX(normalize(L(line)).direction, V(T(0.6), T(0.8)));

	Matrix2<T> r = Matrix2<T>::makeRotation(Constants<T>::HALF_PI);
	L t = L(line).transform(r);
	CHECK_APPROX(t.origin, line.origin*r);
	CHECK_APPROX(t.direction, line.direction*r);
}

TEST_CASE_TEMPLATE("parallelism", T, TEST_REAL_TYPES)
{
	using L = Line2<T>;
	using V = Vector2<T>;
	L a(V(T(0), T(0)), V::UNIT_X);
	L b(V(T(0), T(1)), V(T(-2), T(0)));
	L c(V(T(5), T(0)), V(T(3), T(0)));
	L d(V(T(0), T(0)), V(T(1), T(1)));

	CHECK(a.isParallel(b));
	CHECK(a.isParallel(c));
	CHECK_FALSE(a.isParallel(d));
	CHECK_FALSE(a.isCoincident(b));
	CHECK(a.isCoincident(c));
	CHECK(a.isCoincident(a));
}

TEST_CASE_TEMPLATE("closest point and distance", T, TEST_REAL_TYPES)
{
	using L = Line2<T>;
	using V = Vector2<T>;
	L line(V(T(1), T(1)), V::UNIT_X);
	V p(T(4), T(3));
	CHECK_APPROX(line.getClosestPoint(p), V(T(4), T(1)));
	CHECK_APPROX(line.getDistanceTo(p), T(2));
	CHECK_APPROX(line.getSignedDistanceTo(p), T(2)); // positive on the left side
	CHECK_APPROX(line.getSignedDistanceTo(V(T(4), T(-1))), T(-2));

	L unnormalized(V(T(1), T(1)), V(T(5), T(0)));
	CHECK_APPROX(unnormalized.template getClosestPoint<Unnormalized>(p), V(T(4), T(1)));
	CHECK_APPROX(unnormalized.template getDistanceTo<Unnormalized>(p), T(2));
	CHECK_APPROX(line.template getClosestPoint<Normalized>(p), V(T(4), T(1)));
	CHECK_APPROX(line.template getDistanceTo<Normalized>(p), T(2));

	L parallel(V(T(-3), T(4)), V::UNIT_X);
	CHECK_APPROX(line.getDistanceTo(parallel), T(3));
	CHECK_APPROX(std::fabs(line.getSignedDistanceTo(parallel)), T(3));
	CHECK(line.getDistanceTo(L(V::ZERO, V::UNIT_Y)) == T(0)); // intersecting lines
}

TEST_CASE_TEMPLATE("intersections", T, TEST_REAL_TYPES)
{
	using L = Line2<T>;
	using V = Vector2<T>;
	L line(V(T(0), T(0)), V::UNIT_X);

	SUBCASE("line")
	{
		L other(V(T(3), T(-1)), V(T(0), T(2)));
		CHECK(line.intersects(other));
		auto t = line.findIntersection(other);
		REQUIRE(t.has_value());
		CHECK_APPROX(line.evaluate(*t), V(T(3), T(0))); // parameter along this line

		L parallel(V(T(0), T(1)), V::UNIT_X);
		CHECK_FALSE(line.intersects(parallel));
		CHECK_FALSE(line.findIntersection(parallel).has_value());

		L coincident(V(T(5), T(0)), V(T(-1), T(0)));
		CHECK(line.intersects(coincident));
		CHECK(line.findIntersection(coincident).has_value());
	}

	SUBCASE("ray")
	{
		Ray2<T> hit(V(T(2), T(-1)), V(T(0), T(3)));
		CHECK(line.intersects(hit));
		auto t = line.findIntersection(hit);
		REQUIRE(t.has_value());
		CHECK_APPROX(*t, T(2)); // parameter along this line
		CHECK_APPROX(line.evaluate(*t), V(T(2), T(0)));
		auto u = L(V(T(1), T(0)), V(T(2), T(0))).findIntersection(hit); // unnormalized direction
		REQUIRE(u.has_value());
		CHECK_APPROX(*u, T(0.5));

		Ray2<T> away(V(T(2), T(-1)), V(T(0), T(-1)));
		CHECK_FALSE(line.intersects(away));
		CHECK_FALSE(line.findIntersection(away).has_value());
		CHECK_FALSE(line.intersects(Ray2<T>(V(T(0), T(1)), V::UNIT_X))); // parallel

		// A collinear ray lies on the line; its origin is reported.
		auto c = line.findIntersection(Ray2<T>(V(T(5), T(0)), V::MINUS_UNIT_X));
		REQUIRE(c.has_value());
		CHECK_APPROX(*c, T(5));
		Ray2<T> awayFromOrigin(V(T(5), T(0)), V(T(2), T(0))); // does not contain the line origin
		auto a = line.findIntersection(awayFromOrigin);
		REQUIRE(a.has_value());
		CHECK_APPROX(*a, T(5));
		CHECK_APPROX(awayFromOrigin.getDistanceTo(line.evaluate(*a)), T(0));
		CHECK_APPROX(*L(V(T(1), T(0)), V(T(-2), T(0))).findIntersection(awayFromOrigin), T(-2)); // unnormalized, opposite direction
		CHECK(*awayFromOrigin.findIntersection(line) == T(0)); // the ray reports its own origin
	}

	SUBCASE("segment")
	{
		LineSegment2<T> crossing(V(T(2), T(-1)), V(T(2), T(3)));
		CHECK(line.intersects(crossing));
		auto t = line.findIntersection(crossing);
		REQUIRE(t.has_value());
		CHECK_APPROX(*t, T(2)); // parameter along this line
		CHECK_APPROX(line.evaluate(*t), V(T(2), T(0)));
		auto u = L(V(T(4), T(0)), V(T(-4), T(0))).findIntersection(crossing); // unnormalized, opposite direction
		REQUIRE(u.has_value());
		CHECK_APPROX(*u, T(0.5));

		CHECK_FALSE(line.intersects(LineSegment2<T>(V(T(2), T(1)), V(T(2), T(3)))));
		CHECK_FALSE(line.findIntersection(LineSegment2<T>(V(T(2), T(1)), V(T(2), T(3)))).has_value());

		// A collinear segment lies on the line; its start is reported.
		auto c = line.findIntersection(LineSegment2<T>(V(T(-3), T(0)), V(T(6), T(0))));
		REQUIRE(c.has_value());
		CHECK_APPROX(*c, T(-3));
		LineSegment2<T> awayFromOrigin(V(T(3), T(0)), V(T(6), T(0))); // does not contain the line origin
		auto a = line.findIntersection(awayFromOrigin);
		REQUIRE(a.has_value());
		CHECK_APPROX(*a, T(3));
		CHECK_APPROX(awayFromOrigin.getDistanceTo(line.evaluate(*a)), T(0));
		CHECK_APPROX(*L(V(T(1), T(0)), V(T(-2), T(0))).findIntersection(awayFromOrigin), T(-1)); // unnormalized, opposite direction
		CHECK(*awayFromOrigin.findIntersection(line) == T(0)); // the segment reports its own start
	}

	SUBCASE("rectangle")
	{
		AxisAlignedRectangle<T> rectangle(V(T(1), T(-1)), V(T(3), T(1)));
		CHECK(line.intersects(rectangle));
		auto i = line.findIntersection(rectangle);
		REQUIRE(i.has_value());
		CHECK_APPROX(*i, Interval<T>(T(1), T(3)));

		L diagonal(V(T(0), T(-2)), V(T(1), T(1)));
		auto j = diagonal.findIntersection(rectangle);
		REQUIRE(j.has_value());
		CHECK_APPROX(*j, Interval<T>(T(1), T(3)));

		CHECK_FALSE(L(V(T(0), T(5)), V::UNIT_X).intersects(rectangle));
	}

	SUBCASE("circle")
	{
		Circle2<T> circle(V(T(3), T(0)), T(1));
		CHECK(line.intersects(circle));
		CHECK(line.template intersects<Normalized>(circle));
		auto i = line.findIntersection(circle);
		REQUIRE(i.has_value());
		CHECK_APPROX(*i, Interval<T>(T(2), T(4)));
		auto n = line.template findIntersection<Normalized>(circle);
		REQUIRE(n.has_value());
		CHECK_APPROX(*n, Interval<T>(T(2), T(4)));

		L scaled(V::ZERO, V(T(2), T(0)));
		auto s = scaled.template findIntersection<Unnormalized>(circle);
		REQUIRE(s.has_value());
		CHECK_APPROX(*s, Interval<T>(T(1), T(2)));
		CHECK(scaled.template intersects<Unnormalized>(circle));

		auto tangent = L(V(T(0), T(1)), V::UNIT_X).findIntersection(circle);
		REQUIRE(tangent.has_value());
		CHECK_APPROX(tangent->getLength(), T(0));

		CHECK_FALSE(L(V(T(0), T(2)), V::UNIT_X).intersects(circle));
		CHECK_FALSE(L(V(T(0), T(2)), V::UNIT_X).findIntersection(circle).has_value());
	}
}

TEST_CASE_TEMPLATE("intersection parameters refer to the queried object", T, TEST_REAL_TYPES)
{
	using V = Vector2<T>;
	// Offset origins and unnormalized directions so that parameters of different objects never coincide
	Line2<T> line(V(T(-3), T(1)), V(T(2), T(0.5)));
	Ray2<T> ray(V(T(1), T(-4)), V(T(0.25), T(1.5)));
	LineSegment2<T> segment(V(T(-2), T(5)), V(T(4), T(-3)));
	LineSegment2<T> other(V(T(-1), T(-2)), V(T(3), T(6)));

	auto onLine = [&](const V& p) { return std::fabs(cross(p - line.origin, line.direction)); };
	auto onRay = [&](const V& p) { return std::fabs(cross(p - ray.origin, ray.direction)); };
	auto onSegment = [](const LineSegment2<T>& s, const V& p) { return s.getDistanceTo(p); };

	auto lr = line.findIntersection(ray);
	REQUIRE(lr.has_value());
	CHECK_APPROX(onRay(line.evaluate(*lr)), T(0));
	auto ls = line.findIntersection(segment);
	REQUIRE(ls.has_value());
	CHECK_APPROX(onSegment(segment, line.evaluate(*ls)), T(0));
	auto ll = line.findIntersection(Line2<T>(ray.origin, ray.direction));
	REQUIRE(ll.has_value());
	CHECK_APPROX(onRay(line.evaluate(*ll)), T(0));
	auto rl = ray.findIntersection(line);
	REQUIRE(rl.has_value());
	CHECK_APPROX(onLine(ray.evaluate(*rl)), T(0));
	auto sl = segment.findIntersection(line);
	REQUIRE(sl.has_value());
	CHECK_APPROX(onLine(segment.evaluate(*sl)), T(0));
	auto ss = segment.findIntersection(other);
	REQUIRE(ss.has_value());
	CHECK_APPROX(onSegment(other, segment.evaluate(*ss)), T(0));
	auto os = other.findIntersection(segment);
	REQUIRE(os.has_value());
	CHECK_APPROX(onSegment(segment, other.evaluate(*os)), T(0));
	CHECK_APPROX(segment.evaluate(*ss), other.evaluate(*os)); // both report the same point
}

TEST_CASE_TEMPLATE("streaming and hashing", T, TEST_REAL_TYPES)
{
	using L = Line2<T>;
	using V = Vector2<T>;
	L line(V(T(1), T(2)), V(T(3), T(4)));
	std::stringstream s;
	s << line;
	CHECK(s.str() == "1 2 3 4");
	L read;
	s >> read;
	CHECK(read == line);
	std::unordered_set<L> set{ line, line, L() };
	CHECK(set.size() == 2);
}

} // TEST_SUITE

TEST_SUITE("Ray2")
{

TEST_CASE_TEMPLATE("construction and properties", T, TEST_REAL_TYPES)
{
	using R = Ray2<T>;
	using V = Vector2<T>;
	R ray(V(T(1), T(2)), V(T(3), T(4)));
	CHECK(ray.getOrigin() == V(T(1), T(2)));
	CHECK(ray.getDirection() == V(T(3), T(4)));
	CHECK(R(Line2<T>(V(T(1), T(2)), V(T(3), T(4)))) == ray);
	CHECK(ray.asLine() == Line2<T>(V(T(1), T(2)), V(T(3), T(4))));
	CHECK_APPROX(ray.getSlope(), T(4)/T(3));
	CHECK_APPROX(ray.getInclinationAngle(), std::atan2(T(4), T(3)));

	ray.setOrigin(V::ZERO);
	ray.setDirection(V::UNIT_X);
	CHECK(ray == R(V::ZERO, V::UNIT_X));
	CHECK(R(ray).set(V::ONE, V::UNIT_Y) == R(V::ONE, V::UNIT_Y));
	CHECK(R(ray).setZero().isZero());
	CHECK(R().isApproxZero());
	CHECK(ray.approxEquals(R(V::ZERO, V(T(1) + Constants<T>::TOLERANCE*T(0.5), T(0)))));
	CHECK(ray.approxEquals(R(V::ZERO, V(T(1.05), T(0))), T(0.1)));
	CHECK(ray != R(V::ZERO, V::UNIT_Y));

	CHECK(ray.evaluate(T(2)) == V(T(2), T(0)));
	CHECK(R(ray).translate(V::ONE) == R(V::ONE, V::UNIT_X));
	CHECK_APPROX(R(V::ZERO, V(T(0), T(3))).normalize().direction, V::UNIT_Y);
	Matrix2<T> r = Matrix2<T>::makeRotation(Constants<T>::HALF_PI);
	CHECK_APPROX(R(V::UNIT_X, V::UNIT_X).transform(r), R(V::UNIT_Y, V::UNIT_Y));
}

TEST_CASE_TEMPLATE("closest point and distance", T, TEST_REAL_TYPES)
{
	using R = Ray2<T>;
	using V = Vector2<T>;
	R ray(V(T(1), T(1)), V::UNIT_X);
	CHECK_APPROX(ray.getClosestPoint(V(T(4), T(3))), V(T(4), T(1)));
	CHECK_APPROX(ray.getClosestPoint(V(T(-2), T(3))), V(T(1), T(1))); // behind the origin
	CHECK_APPROX(ray.getDistanceTo(V(T(4), T(3))), T(2));
	CHECK_APPROX(ray.getDistanceTo(V(T(-2), T(5))), T(5));

	R unnormalized(V(T(1), T(1)), V(T(3), T(0)));
	CHECK_APPROX(unnormalized.template getClosestPoint<Unnormalized>(V(T(4), T(3))), V(T(4), T(1)));
	CHECK_APPROX(unnormalized.template getDistanceTo<Unnormalized>(V(T(-2), T(5))), T(5));
	CHECK_APPROX(ray.template getClosestPoint<Normalized>(V(T(4), T(3))), V(T(4), T(1)));
	CHECK_APPROX(ray.template getDistanceTo<Normalized>(V(T(4), T(3))), T(2));
}

TEST_CASE_TEMPLATE("intersections", T, TEST_REAL_TYPES)
{
	using R = Ray2<T>;
	using V = Vector2<T>;
	R ray(V(T(0), T(0)), V::UNIT_X);

	SUBCASE("line")
	{
		Line2<T> line(V(T(3), T(-1)), V::UNIT_Y);
		CHECK(ray.intersects(line));
		auto t = ray.findIntersection(line);
		REQUIRE(t.has_value());
		CHECK_APPROX(*t, T(3));
		CHECK_FALSE(ray.intersects(Line2<T>(V(T(-3), T(-1)), V::UNIT_Y))); // behind the ray
	}

	SUBCASE("rectangle")
	{
		AxisAlignedRectangle<T> ahead(V(T(1), T(-1)), V(T(3), T(1)));
		auto a = ray.findIntersection(ahead);
		REQUIRE(a.has_value());
		CHECK_APPROX(*a, Interval<T>(T(1), T(3)));
		CHECK(ray.intersects(ahead));

		AxisAlignedRectangle<T> around(V(T(-1), T(-1)), V(T(3), T(1)));
		auto b = ray.findIntersection(around);
		REQUIRE(b.has_value());
		CHECK_APPROX(*b, Interval<T>(T(0), T(3))); // clipped to the origin

		CHECK_FALSE(ray.intersects(AxisAlignedRectangle<T>(V(T(-4), T(-1)), V(T(-2), T(1)))));
	}

	SUBCASE("circle")
	{
		Circle2<T> ahead(V(T(3), T(0)), T(1));
		auto a = ray.findIntersection(ahead);
		REQUIRE(a.has_value());
		CHECK_APPROX(*a, Interval<T>(T(2), T(4)));
		CHECK(ray.intersects(ahead));
		CHECK(ray.template intersects<Normalized>(ahead));
		auto n = ray.template findIntersection<Normalized>(ahead);
		REQUIRE(n.has_value());
		CHECK_APPROX(*n, Interval<T>(T(2), T(4)));

		Circle2<T> around(V(T(0), T(0)), T(2));
		auto b = ray.findIntersection(around);
		REQUIRE(b.has_value());
		CHECK_APPROX(*b, Interval<T>(T(0), T(2)));
		auto u = R(V::ZERO, V(T(2), T(0))).template findIntersection<Unnormalized>(around);
		REQUIRE(u.has_value());
		CHECK_APPROX(*u, Interval<T>(T(0), T(1)));

		CHECK_FALSE(ray.intersects(Circle2<T>(V(T(-3), T(0)), T(1))));
		CHECK_FALSE(ray.template intersects<Normalized>(Circle2<T>(V(T(-3), T(0)), T(1))));
	}
}

TEST_CASE_TEMPLATE("streaming and hashing", T, TEST_REAL_TYPES)
{
	using R = Ray2<T>;
	using V = Vector2<T>;
	R ray(V(T(1), T(2)), V(T(3), T(4)));
	std::stringstream s;
	s << ray;
	R read;
	s >> read;
	CHECK(read == ray);
	std::unordered_set<R> set{ ray, ray, R() };
	CHECK(set.size() == 2);
}

} // TEST_SUITE

TEST_SUITE("LineSegment2")
{

TEST_CASE_TEMPLATE("construction and properties", T, TEST_REAL_TYPES)
{
	using S = LineSegment2<T>;
	using V = Vector2<T>;
	S s(V(T(1), T(1)), V(T(4), T(5)));
	CHECK(s.getStart() == V(T(1), T(1)));
	CHECK(s.getEnd() == V(T(4), T(5)));
	CHECK(S(std::pair<V, V>(V(T(1), T(1)), V(T(4), T(5)))) == s);
	CHECK(S(std::tuple<V, V>(V(T(1), T(1)), V(T(4), T(5)))) == s);
	CHECK(S(Line2<T>(V(T(1), T(1)), V(T(3), T(4)))) == s);
	CHECK(S(Ray2<T>(V(T(1), T(1)), V(T(3), T(4)))) == s);
	CHECK_APPROX(S(Line2<T>(V(T(1), T(1)), V(T(0.6), T(0.8))), Interval<T>(T(0), T(5))), s);
	CHECK_APPROX(S(Ray2<T>(V(T(1), T(1)), V(T(0.6), T(0.8))), Interval<T>(T(0), T(5))), s);

	CHECK_APPROX(s.getDirection(), V(T(0.6), T(0.8)));
	CHECK_APPROX(s.getLength(), T(5));
	CHECK_APPROX(s.getCenter(), V(T(2.5), T(3)));
	CHECK_APPROX(s.getSlope(), T(4)/T(3));
	CHECK_APPROX(s.getInclinationAngle(), std::atan2(T(4), T(3)));
	CHECK(s.getEndpoints() == std::pair<V, V>(V(T(1), T(1)), V(T(4), T(5))));
	std::vector<V> endpoints;
	s.copyEndpoints(std::back_inserter(endpoints));
	CHECK(endpoints == std::vector<V>{ V(T(1), T(1)), V(T(4), T(5)) });

	CHECK(s.isFinite());
	CHECK_FALSE(S(V::ZERO, V::INF).isFinite());
	CHECK(S(s).setZero().isZero());
	CHECK(S().isApproxZero());
	CHECK(S().set(V::ZERO, V::ONE) == S(V::ZERO, V::ONE));
	S t = s;
	t.setStart(V::ZERO);
	t.setEnd(V::ONE);
	CHECK(t == S(V::ZERO, V::ONE));
	CHECK(s.approxEquals(S(V(T(1), T(1)), V(T(4), T(5) + Constants<T>::TOLERANCE*T(0.5)))));
	CHECK(s.approxEquals(S(V(T(1), T(1)), V(T(4), T(5.05))), T(0.1)));
	CHECK(s != t);
}

TEST_CASE_TEMPLATE("transformation, evaluation and closest point", T, TEST_REAL_TYPES)
{
	using S = LineSegment2<T>;
	using V = Vector2<T>;
	S s(V(T(0), T(0)), V(T(4), T(0)));
	CHECK(s.evaluate(T(0.25)) == V(T(1), T(0)));
	CHECK(S(s).translate(V::ONE) == S(V::ONE, V(T(5), T(1))));
	Matrix2<T> r = Matrix2<T>::makeRotation(Constants<T>::HALF_PI);
	CHECK_APPROX(S(s).transform(r), S(V::ZERO, V(T(0), T(4))));

	CHECK_APPROX(s.getClosestPoint(V(T(2), T(3))), V(T(2), T(0)));
	CHECK_APPROX(s.getClosestPoint(V(T(-2), T(3))), V(T(0), T(0)));
	CHECK_APPROX(s.getClosestPoint(V(T(7), T(-3))), V(T(4), T(0)));
	CHECK_APPROX(s.getDistanceTo(V(T(7), T(4))), T(5));
}

TEST_CASE_TEMPLATE("intersections", T, TEST_REAL_TYPES)
{
	using S = LineSegment2<T>;
	using V = Vector2<T>;
	S s(V(T(0), T(0)), V(T(4), T(0)));

	SUBCASE("line")
	{
		auto t = s.findIntersection(Line2<T>(V(T(1), T(-1)), V::UNIT_Y));
		REQUIRE(t.has_value());
		CHECK_APPROX(*t, T(0.25));
		CHECK(s.intersects(Line2<T>(V(T(1), T(-1)), V::UNIT_Y)));
		CHECK_FALSE(s.intersects(Line2<T>(V(T(5), T(-1)), V::UNIT_Y)));
	}

	SUBCASE("segment")
	{
		auto t = s.findIntersection(S(V(T(3), T(-1)), V(T(3), T(1))));
		REQUIRE(t.has_value());
		CHECK_APPROX(*t, T(0.75));
		CHECK(s.intersects(S(V(T(3), T(-1)), V(T(3), T(1)))));
		CHECK_FALSE(s.intersects(S(V(T(3), T(1)), V(T(3), T(2)))));
		CHECK_FALSE(s.intersects(S(V(T(0), T(1)), V(T(4), T(1))))); // parallel

		auto overlap = s.findIntersection(S(V(T(2), T(0)), V(T(6), T(0)))); // collinear overlap
		REQUIRE(overlap.has_value());
		CHECK_APPROX(*overlap, T(0.5));
		CHECK_FALSE(s.intersects(S(V(T(5), T(0)), V(T(6), T(0))))); // collinear, disjoint
	}

	SUBCASE("rectangle")
	{
		auto i = s.findIntersection(AxisAlignedRectangle<T>(V(T(1), T(-1)), V(T(2), T(1))));
		REQUIRE(i.has_value());
		CHECK_APPROX(*i, Interval<T>(T(0.25), T(0.5)));
		auto j = s.findIntersection(AxisAlignedRectangle<T>(V(T(-1), T(-1)), V(T(10), T(1))));
		REQUIRE(j.has_value());
		CHECK_APPROX(*j, Interval<T>(T(0), T(1))); // clipped to the segment
		CHECK(s.intersects(AxisAlignedRectangle<T>(V(T(1), T(-1)), V(T(2), T(1)))));
		CHECK_FALSE(s.intersects(AxisAlignedRectangle<T>(V(T(5), T(-1)), V(T(6), T(1)))));
	}

	SUBCASE("circle")
	{
		auto i = s.findIntersection(Circle2<T>(V(T(2), T(0)), T(1)));
		REQUIRE(i.has_value());
		CHECK_APPROX(*i, Interval<T>(T(0.25), T(0.75)));
		auto j = s.findIntersection(Circle2<T>(V(T(0), T(0)), T(2)));
		REQUIRE(j.has_value());
		CHECK_APPROX(*j, Interval<T>(T(0), T(0.5)));
		CHECK(s.intersects(Circle2<T>(V(T(2), T(0)), T(1))));
		CHECK_FALSE(s.intersects(Circle2<T>(V(T(7), T(0)), T(1))));
	}
}

TEST_CASE_TEMPLATE("streaming and hashing", T, TEST_REAL_TYPES)
{
	using S = LineSegment2<T>;
	using V = Vector2<T>;
	S seg(V(T(1), T(2)), V(T(3), T(4)));
	std::stringstream s;
	s << seg;
	S read;
	s >> read;
	CHECK(read == seg);
	std::unordered_set<S> set{ seg, seg, S() };
	CHECK(set.size() == 2);
}

} // TEST_SUITE

} // namespace mathematics::templates
