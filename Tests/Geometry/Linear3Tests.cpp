#include "TestUtilities.hpp"
#include <unordered_set>
#include <memory>
#include <vector>
#include <iterator>

namespace mathematics::templates {

namespace {

template<typename T>
struct Scene
{
	using V = Vector3<T>;
	Plane<T> plane{ V::UNIT_X, T(-3) };	// x = 3
	Triangle3<T> triangle{ V(T(3), T(-1), T(-1)), V(T(3), T(2), T(-1)), V(T(3), T(-1), T(2)) };
	AxisAlignedBox<T> box{ V(T(2), T(-1), T(-1)), V(T(4), T(1), T(1)) };
	OrientedBox<T> orientedBox{ V(T(3), T(0), T(0)), Matrix3<T>::makeRotation(Axis::X, T(0.4)), V(T(1)) };
	Sphere<T> sphere{ V(T(3), T(0), T(0)), T(1) };
	Ellipsoid<T> ellipsoid{ V(T(5), T(0), T(0)), Matrix3<T>::IDENTITY, V(T(2), T(1), T(1)) };
	std::unique_ptr<TriangleMesh<T>> mesh{ TriangleMesh<T>::makeBox(V(T(2))) }; // [-1, 1]^3
};

} // namespace

TEST_SUITE("Line3")
{

TEST_CASE_TEMPLATE("construction and properties", T, TEST_REAL_TYPES)
{
	using L = Line3<T>;
	using V = Vector3<T>;
	L line(V(T(1), T(2), T(3)), V(T(0), T(3), T(4)));
	CHECK(line.getOrigin() == V(T(1), T(2), T(3)));
	CHECK(line.getDirection() == V(T(0), T(3), T(4)));
	CHECK(L(Ray3<T>(line.origin, line.direction)) == line);
	CHECK(line.asRay() == Ray3<T>(line.origin, line.direction));

	L l = line;
	l.setOrigin(V::ZERO);
	l.setDirection(V::UNIT_X);
	CHECK(l == L(V::ZERO, V::UNIT_X));
	CHECK(L(line).set(V::ONE, V::UNIT_Z) == L(V::ONE, V::UNIT_Z));
	CHECK(L(line).setZero().isZero());
	CHECK(L().isApproxZero());
	CHECK(line.approxEquals(L(line.origin, line.direction + V(Constants<T>::TOLERANCE*T(0.5)))));
	CHECK(line.approxEquals(L(line.origin, line.direction + V(T(0.05))), T(0.1)));
	CHECK(line != l);

	CHECK(line.evaluate(T(2)) == V(T(1), T(8), T(11)));
	CHECK(L(line).translate(V::ONE) == L(V(T(2), T(3), T(4)), line.direction));
	CHECK_APPROX(L(line).normalize().direction, V(T(0), T(0.6), T(0.8)));
	CHECK_APPROX(normalize(line).direction, V(T(0), T(0.6), T(0.8)));
	CHECK_APPROX(normalize(L(line)).direction, V(T(0), T(0.6), T(0.8)));

	Matrix3<T> m = Matrix3<T>::makeRotation(Axis::Z, T(0.5));
	CHECK_APPROX(L(line).transform(m), L(line.origin*m, line.direction*m));
	AffineTransform<T> a = AffineTransform<T>::makeRotationTranslation(Axis::Z, T(0.5), V(T(1), T(2), T(3)));
	CHECK_APPROX(L(line).transform(a), L(transform(line.origin, a), line.direction*a.getBasis()));
}

TEST_CASE_TEMPLATE("closest point and distance", T, TEST_REAL_TYPES)
{
	using L = Line3<T>;
	using V = Vector3<T>;
	L line(V(T(1), T(1), T(1)), V::UNIT_X);
	V p(T(5), T(4), T(5));
	CHECK_APPROX(line.getClosestPoint(p), V(T(5), T(1), T(1)));
	CHECK_APPROX(line.getDistanceTo(p), T(5));
	L u(V(T(1), T(1), T(1)), V(T(3), T(0), T(0)));
	CHECK_APPROX(u.template getClosestPoint<Unnormalized>(p), V(T(5), T(1), T(1)));
	CHECK_APPROX(u.template getDistanceTo<Unnormalized>(p), T(5));
	CHECK_APPROX(line.template getClosestPoint<Normalized>(p), V(T(5), T(1), T(1)));
	CHECK_APPROX(line.template getDistanceTo<Normalized>(p), T(5));
}

TEST_CASE_TEMPLATE("intersections", T, TEST_REAL_TYPES)
{
	using L = Line3<T>;
	using V = Vector3<T>;
	Scene<T> s;
	L line(V(T(0), T(0), T(0)), V::UNIT_X);
	L miss(V(T(0), T(5), T(0)), V::UNIT_X);
	L parallel(V(T(0), T(0), T(0)), V::UNIT_Y);

	CHECK(line.intersects(s.plane));
	CHECK_FALSE(parallel.intersects(s.plane));
	REQUIRE(line.findIntersection(s.plane).has_value());
	CHECK_APPROX(*line.findIntersection(s.plane), T(3));
	CHECK_FALSE(parallel.findIntersection(s.plane).has_value());

	CHECK(line.intersects(s.triangle));
	REQUIRE(line.findIntersection(s.triangle).has_value());
	CHECK_APPROX(*line.findIntersection(s.triangle), T(3));
	CHECK(L(V(T(5), T(0), T(0)), V::MINUS_UNIT_X).intersects(s.triangle)); // back face
	CHECK_FALSE(L(V(T(0), T(1), T(1)), V::UNIT_X).intersects(s.triangle));

	CHECK(line.intersects(s.box));
	REQUIRE(line.findIntersection(s.box).has_value());
	CHECK_APPROX(*line.findIntersection(s.box), Interval<T>(T(2), T(4)));
	CHECK_FALSE(miss.intersects(s.box));

	CHECK(line.intersects(s.orientedBox));
	REQUIRE(line.findIntersection(s.orientedBox).has_value());
	CHECK_APPROX(*line.findIntersection(s.orientedBox), Interval<T>(T(2), T(4))); // rotation about the line direction
	CHECK_FALSE(miss.intersects(s.orientedBox));

	CHECK(line.intersects(s.sphere));
	CHECK(line.template intersects<Normalized>(s.sphere));
	CHECK(L(V::ZERO, V(T(2), T(0), T(0))).template intersects<Unnormalized>(s.sphere));
	REQUIRE(line.findIntersection(s.sphere).has_value());
	CHECK_APPROX(*line.findIntersection(s.sphere), Interval<T>(T(2), T(4)));
	CHECK_APPROX(*line.template findIntersection<Normalized>(s.sphere), Interval<T>(T(2), T(4)));
	CHECK_APPROX(*L(V::ZERO, V(T(2), T(0), T(0))).template findIntersection<Unnormalized>(s.sphere), Interval<T>(T(1), T(2)));
	CHECK_FALSE(miss.intersects(s.sphere));
	CHECK_FALSE(miss.findIntersection(s.sphere).has_value());

	CHECK(line.intersects(s.ellipsoid));
	REQUIRE(line.findIntersection(s.ellipsoid).has_value());
	CHECK_APPROX(*line.findIntersection(s.ellipsoid), Interval<T>(T(3), T(7)));
	CHECK_FALSE(miss.intersects(s.ellipsoid));
	CHECK_FALSE(miss.findIntersection(s.ellipsoid).has_value());

	L through(V(T(-5), T(0.3), T(0.6)), V::UNIT_X);
	CHECK(through.intersects(s.mesh.get()));
	std::vector<T> hits = through.findIntersections(s.mesh.get());
	std::sort(hits.begin(), hits.end());
	REQUIRE(hits.size() == 2);
	CHECK_APPROX(hits[0], T(4));
	CHECK_APPROX(hits[1], T(6));
	CHECK_FALSE(miss.intersects(s.mesh.get()));
	CHECK(miss.findIntersections(s.mesh.get()).empty());
	CHECK_FALSE(line.intersects(static_cast<const TriangleMesh<T>*>(nullptr)));
}

TEST_CASE_TEMPLATE("streaming and hashing", T, TEST_REAL_TYPES)
{
	using L = Line3<T>;
	using V = Vector3<T>;
	L line(V(T(1), T(2), T(3)), V(T(4), T(5), T(6)));
	std::stringstream s;
	s << line;
	CHECK(s.str() == "1 2 3 4 5 6");
	L read;
	s >> read;
	CHECK(read == line);
	std::unordered_set<L> set{ line, line, L() };
	CHECK(set.size() == 2);
}

} // TEST_SUITE

TEST_SUITE("Ray3")
{

TEST_CASE_TEMPLATE("construction and properties", T, TEST_REAL_TYPES)
{
	using R = Ray3<T>;
	using V = Vector3<T>;
	R ray(V(T(1), T(2), T(3)), V(T(0), T(3), T(4)));
	CHECK(R(Line3<T>(ray.origin, ray.direction)) == ray);
	CHECK(ray.asLine() == Line3<T>(ray.origin, ray.direction));
	R r = ray;
	r.setOrigin(V::ZERO);
	r.setDirection(V::UNIT_X);
	CHECK(r == R(V::ZERO, V::UNIT_X));
	CHECK(r.getOrigin() == V::ZERO);
	CHECK(r.getDirection() == V::UNIT_X);
	CHECK(R(ray).set(V::ONE, V::UNIT_Z) == R(V::ONE, V::UNIT_Z));
	CHECK(R(ray).setZero().isZero());
	CHECK(R().isApproxZero());
	CHECK(ray.approxEquals(R(ray.origin, ray.direction + V(T(0.05))), T(0.1)));
	CHECK_FALSE(ray.approxEquals(R(ray.origin, ray.direction + V(T(0.05)))));
	CHECK(ray != r);

	CHECK(ray.evaluate(T(2)) == V(T(1), T(8), T(11)));
	CHECK(R(ray).translate(V::ONE) == R(V(T(2), T(3), T(4)), ray.direction));
	CHECK_APPROX(R(ray).normalize().direction, V(T(0), T(0.6), T(0.8)));
	CHECK_APPROX(normalize(ray).direction, V(T(0), T(0.6), T(0.8)));
	CHECK_APPROX(normalize(R(ray)).direction, V(T(0), T(0.6), T(0.8)));
	Matrix3<T> m = Matrix3<T>::makeRotation(Axis::Z, T(0.5));
	CHECK_APPROX(R(ray).transform(m), R(ray.origin*m, ray.direction*m));
	AffineTransform<T> a = AffineTransform<T>::makeRotationTranslation(Axis::Z, T(0.5), V(T(1), T(2), T(3)));
	CHECK_APPROX(R(ray).transform(a), R(transform(ray.origin, a), ray.direction*a.getBasis()));
}

TEST_CASE_TEMPLATE("closest point and distance", T, TEST_REAL_TYPES)
{
	using R = Ray3<T>;
	using V = Vector3<T>;
	R ray(V(T(1), T(1), T(1)), V::UNIT_X);
	CHECK_APPROX(ray.getClosestPoint(V(T(5), T(4), T(5))), V(T(5), T(1), T(1)));
	CHECK_APPROX(ray.getClosestPoint(V(T(-5), T(4), T(5))), V(T(1), T(1), T(1)));
	CHECK_APPROX(ray.getDistanceTo(V(T(5), T(4), T(5))), T(5));
	CHECK_APPROX(ray.getDistanceTo(V(T(1), T(4), T(5))), T(5));
	R u(V(T(1), T(1), T(1)), V(T(3), T(0), T(0)));
	CHECK_APPROX(u.template getClosestPoint<Unnormalized>(V(T(5), T(4), T(5))), V(T(5), T(1), T(1)));
	CHECK_APPROX(u.template getDistanceTo<Unnormalized>(V(T(-2), T(5), T(1))), T(5));
	CHECK_APPROX(ray.template getClosestPoint<Normalized>(V(T(5), T(4), T(5))), V(T(5), T(1), T(1)));
	CHECK_APPROX(ray.template getDistanceTo<Normalized>(V(T(5), T(4), T(5))), T(5));
}

TEST_CASE_TEMPLATE("intersections", T, TEST_REAL_TYPES)
{
	using R = Ray3<T>;
	using V = Vector3<T>;
	Scene<T> s;
	R ray(V(T(0), T(0), T(0)), V::UNIT_X);
	R inside(V(T(3), T(0), T(0)), V::UNIT_X);
	R away(V(T(0), T(0), T(0)), V::MINUS_UNIT_X);

	CHECK(ray.intersects(s.plane));
	CHECK_APPROX(*ray.findIntersection(s.plane), T(3));
	CHECK_FALSE(away.intersects(s.plane));

	CHECK(ray.intersects(s.triangle));
	CHECK_APPROX(*ray.findIntersection(s.triangle), T(3));
	CHECK_FALSE(away.intersects(s.triangle));

	CHECK_APPROX(*ray.findIntersection(s.box), Interval<T>(T(2), T(4)));
	CHECK_APPROX(*inside.findIntersection(s.box), Interval<T>(T(0), T(1)));
	CHECK(ray.intersects(s.box));
	CHECK_FALSE(away.intersects(s.box));

	CHECK_APPROX(*ray.findIntersection(s.orientedBox), Interval<T>(T(2), T(4)));
	CHECK_APPROX(*inside.findIntersection(s.orientedBox), Interval<T>(T(0), T(1)));
	CHECK(ray.intersects(s.orientedBox));
	CHECK_FALSE(away.intersects(s.orientedBox));

	CHECK_APPROX(*ray.findIntersection(s.sphere), Interval<T>(T(2), T(4)));
	CHECK_APPROX(*inside.findIntersection(s.sphere), Interval<T>(T(0), T(1)));
	CHECK_APPROX(*ray.template findIntersection<Normalized>(s.sphere), Interval<T>(T(2), T(4)));
	CHECK_APPROX(*R(V::ZERO, V(T(2), T(0), T(0))).template findIntersection<Unnormalized>(s.sphere), Interval<T>(T(1), T(2)));
	CHECK(ray.intersects(s.sphere));
	CHECK(ray.template intersects<Normalized>(s.sphere));
	CHECK_FALSE(away.intersects(s.sphere));
	CHECK_FALSE(away.template intersects<Normalized>(s.sphere));

	CHECK_APPROX(*ray.findIntersection(s.ellipsoid), Interval<T>(T(3), T(7)));
	CHECK_APPROX(*R(V(T(5), T(0), T(0)), V::UNIT_X).findIntersection(s.ellipsoid), Interval<T>(T(0), T(2)));
	CHECK(ray.intersects(s.ellipsoid));
	CHECK_FALSE(away.intersects(s.ellipsoid));
	CHECK_FALSE(away.findIntersection(s.ellipsoid).has_value());

	R through(V(T(-5), T(0.3), T(0.6)), V::UNIT_X);
	CHECK(through.intersects(s.mesh.get()));
	CHECK(through.findIntersections(s.mesh.get()).size() == 2);
	R fromInside(V(T(0), T(0.3), T(0.6)), V::UNIT_X);
	std::vector<T> hits = fromInside.findIntersections(s.mesh.get());
	REQUIRE(hits.size() == 1);
	CHECK_APPROX(hits[0], T(1));
	CHECK_FALSE(R(V(T(5), T(0.3), T(0.6)), V::UNIT_X).intersects(s.mesh.get()));
	CHECK_FALSE(ray.intersects(static_cast<const TriangleMesh<T>*>(nullptr)));
}

TEST_CASE_TEMPLATE("streaming and hashing", T, TEST_REAL_TYPES)
{
	using R = Ray3<T>;
	using V = Vector3<T>;
	R ray(V(T(1), T(2), T(3)), V(T(4), T(5), T(6)));
	std::stringstream s;
	s << ray;
	R read;
	s >> read;
	CHECK(read == ray);
	std::unordered_set<R> set{ ray, ray, R() };
	CHECK(set.size() == 2);
}

} // TEST_SUITE

TEST_SUITE("LineSegment3")
{

TEST_CASE_TEMPLATE("construction and properties", T, TEST_REAL_TYPES)
{
	using S = LineSegment3<T>;
	using V = Vector3<T>;
	S s(V(T(1), T(1), T(1)), V(T(1), T(4), T(5)));
	CHECK(S(std::pair<V, V>(s.start, s.end)) == s);
	CHECK(S(std::tuple<V, V>(s.start, s.end)) == s);
	CHECK(S(Line3<T>(s.start, V(T(0), T(3), T(4)))) == s);
	CHECK(S(Ray3<T>(s.start, V(T(0), T(3), T(4)))) == s);
	CHECK_APPROX(S(Line3<T>(s.start, V(T(0), T(0.6), T(0.8))), Interval<T>(T(0), T(5))), s);
	CHECK_APPROX(S(Ray3<T>(s.start, V(T(0), T(0.6), T(0.8))), Interval<T>(T(0), T(5))), s);

	CHECK(s.getStart() == V(T(1), T(1), T(1)));
	CHECK(s.getEnd() == V(T(1), T(4), T(5)));
	CHECK_APPROX(s.getDirection(), V(T(0), T(0.6), T(0.8)));
	CHECK_APPROX(s.getLength(), T(5));
	CHECK_APPROX(s.getCenter(), V(T(1), T(2.5), T(3)));
	CHECK(s.getEndpoints() == std::pair<V, V>(s.start, s.end));
	std::vector<V> endpoints;
	s.copyEndpoints(std::back_inserter(endpoints));
	CHECK(endpoints == std::vector<V>{ s.start, s.end });

	S t = s;
	t.setStart(V::ZERO);
	t.setEnd(V::ONE);
	CHECK(t == S(V::ZERO, V::ONE));
	CHECK(S().set(V::ZERO, V::ONE) == t);
	CHECK(S(s).setZero().isZero());
	CHECK(S().isApproxZero());
	CHECK(s.isFinite());
	CHECK_FALSE(S(V::ZERO, V::INF).isFinite());
	CHECK(s.approxEquals(S(s.start, s.end + V(T(0.05))), T(0.1)));
	CHECK_FALSE(s.approxEquals(S(s.start, s.end + V(T(0.05)))));
	CHECK(s != t);
}

TEST_CASE_TEMPLATE("transformation, evaluation and closest point", T, TEST_REAL_TYPES)
{
	using S = LineSegment3<T>;
	using V = Vector3<T>;
	S s(V(T(0), T(0), T(0)), V(T(4), T(0), T(0)));
	CHECK(s.evaluate(T(0.25)) == V(T(1), T(0), T(0)));
	CHECK(S(s).translate(V::ONE) == S(V::ONE, V(T(5), T(1), T(1))));
	Matrix3<T> m = Matrix3<T>::makeRotation(Axis::Z, Constants<T>::HALF_PI);
	CHECK_APPROX(S(s).transform(m), S(V::ZERO, V(T(0), T(4), T(0))));
	AffineTransform<T> a = AffineTransform<T>::makeRotationTranslation(Axis::Z, Constants<T>::HALF_PI, V::ONE);
	CHECK_APPROX(S(s).transform(a), S(V::ONE, V(T(1), T(5), T(1))));

	CHECK_APPROX(s.getClosestPoint(V(T(2), T(3), T(4))), V(T(2), T(0), T(0)));
	CHECK_APPROX(s.getClosestPoint(V(T(-2), T(3), T(0))), V::ZERO);
	CHECK_APPROX(s.getClosestPoint(V(T(7), T(-3), T(0))), V(T(4), T(0), T(0)));
	CHECK_APPROX(s.getDistanceTo(V(T(7), T(4), T(0))), T(5));
}

TEST_CASE_TEMPLATE("intersections", T, TEST_REAL_TYPES)
{
	using S = LineSegment3<T>;
	using V = Vector3<T>;
	Scene<T> sc;
	S seg(V(T(0), T(0), T(0)), V(T(10), T(0), T(0)));
	S shortSeg(V(T(0), T(0), T(0)), V(T(1), T(0), T(0)));

	CHECK(seg.intersects(HalfSpace<T>(V::UNIT_X, T(-3)))); // start inside x <= 3
	CHECK(seg.intersects(HalfSpace<T>(V::MINUS_UNIT_X, T(3)))); // end inside x >= 3
	CHECK_FALSE(S(V(T(5), T(0), T(0)), V(T(6), T(0), T(0))).intersects(HalfSpace<T>(V::UNIT_X, T(-3))));

	CHECK(seg.intersects(sc.plane));
	CHECK_APPROX(*seg.findIntersection(sc.plane), T(0.3));
	CHECK_FALSE(shortSeg.intersects(sc.plane));

	CHECK(seg.intersects(sc.triangle));
	CHECK_APPROX(*seg.findIntersection(sc.triangle), T(0.3));
	CHECK_FALSE(shortSeg.intersects(sc.triangle));

	CHECK_APPROX(*seg.findIntersection(sc.box), Interval<T>(T(0.2), T(0.4)));
	CHECK_APPROX(*S(V(T(3), T(0), T(0)), V(T(13), T(0), T(0))).findIntersection(sc.box), Interval<T>(T(0), T(0.1)));
	CHECK(seg.intersects(sc.box));
	CHECK_FALSE(shortSeg.intersects(sc.box));

	CHECK_APPROX(*seg.findIntersection(sc.orientedBox), Interval<T>(T(0.2), T(0.4)));
	CHECK(seg.intersects(sc.orientedBox));
	CHECK_FALSE(shortSeg.intersects(sc.orientedBox));

	CHECK_APPROX(*seg.findIntersection(sc.sphere), Interval<T>(T(0.2), T(0.4)));
	CHECK_APPROX(*S(V(T(3), T(0), T(0)), V(T(13), T(0), T(0))).findIntersection(sc.sphere), Interval<T>(T(0), T(0.1)));
	CHECK(seg.intersects(sc.sphere));
	CHECK_FALSE(shortSeg.intersects(sc.sphere));

	CHECK_APPROX(*seg.findIntersection(sc.ellipsoid), Interval<T>(T(0.3), T(0.7)));
	CHECK_APPROX(*S(V(T(5), T(0), T(0)), V(T(15), T(0), T(0))).findIntersection(sc.ellipsoid), Interval<T>(T(0), T(0.2)));
	CHECK(seg.intersects(sc.ellipsoid));
	CHECK(S(V(T(5), T(0), T(0)), V(T(5.5), T(0), T(0))).intersects(sc.ellipsoid)); // fully inside
	CHECK_FALSE(shortSeg.intersects(sc.ellipsoid));
	CHECK_FALSE(shortSeg.findIntersection(sc.ellipsoid).has_value());
	// Approaching from either side (both branches of the test), reaching in or stopping short
	CHECK(S(V(T(0), T(0), T(0)), V(T(4), T(0), T(0))).intersects(sc.ellipsoid));
	CHECK(S(V(T(10), T(0), T(0)), V(T(6), T(0), T(0))).intersects(sc.ellipsoid));
	CHECK_FALSE(S(V(T(10), T(0), T(0)), V(T(9), T(0), T(0))).intersects(sc.ellipsoid));
	CHECK_FALSE(S(V(T(9), T(0), T(0)), V(T(10), T(0), T(0))).intersects(sc.ellipsoid));

	S through(V(T(-5), T(0.3), T(0.6)), V(T(5), T(0.3), T(0.6)));
	CHECK(through.intersects(sc.mesh.get()));
	std::vector<T> hits = through.findIntersections(sc.mesh.get());
	std::sort(hits.begin(), hits.end());
	REQUIRE(hits.size() == 2);
	CHECK_APPROX(hits[0], T(0.4));
	CHECK_APPROX(hits[1], T(0.6));
	CHECK_FALSE(S(V(T(-5), T(0.3), T(0.6)), V(T(-2), T(0.3), T(0.6))).intersects(sc.mesh.get()));
	CHECK_FALSE(seg.intersects(static_cast<const TriangleMesh<T>*>(nullptr)));
}

TEST_CASE_TEMPLATE("streaming and hashing", T, TEST_REAL_TYPES)
{
	using S = LineSegment3<T>;
	using V = Vector3<T>;
	S seg(V(T(1), T(2), T(3)), V(T(4), T(5), T(6)));
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
