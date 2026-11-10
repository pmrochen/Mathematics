#include "TestUtilities.hpp"
#include <unordered_set>
#include <vector>
#include <iterator>

namespace mathematics::templates {

TEST_SUITE("Circle2")
{

TEST_CASE_TEMPLATE("construction and properties", T, TEST_REAL_TYPES)
{
	using C = Circle2<T>;
	using V = Vector2<T>;
	CHECK(C() == C(V::ZERO, T(0)));
	C c(V(T(1), T(2)), T(3));
	CHECK(c.getCenter() == V(T(1), T(2)));
	CHECK(c.getRadius() == T(3));
	CHECK(c.getDiameter() == T(6));
	CHECK_APPROX(c.getCircumference(), T(6)*Constants<T>::PI);
	CHECK_APPROX(c.getArea(), T(9)*Constants<T>::PI);

	C d = c;
	d.setCenter(V::ZERO);
	d.setRadius(T(2));
	CHECK(d == C(V::ZERO, T(2)));
	d.setDiameter(T(10));
	CHECK(d.radius == T(5));
	CHECK(C().set(V::ONE, T(1)) == C(V::ONE, T(1)));
	CHECK(C(c).translate(V(T(1), T(-2))) == C(V(T(2), T(0)), T(3)));

	CHECK(c.isFinite());
	CHECK_FALSE(C(V::ZERO, std::numeric_limits<T>::infinity()).isFinite());
	CHECK(c.approxEquals(C(V(T(1), T(2)), T(3) + Constants<T>::TOLERANCE*T(0.5))));
	CHECK(c.approxEquals(C(V(T(1), T(2)), T(3.05)), T(0.1)));
	CHECK_FALSE(c.approxEquals(C(V(T(1), T(2)), T(3.05))));
	CHECK(c != d);

	CHECK(c.getCircumscribedRectangle() == AxisAlignedRectangle<T>(V(T(-2), T(-1)), V(T(4), T(5))));
}

TEST_CASE_TEMPLATE("distance, containment and intersection", T, TEST_REAL_TYPES)
{
	using C = Circle2<T>;
	using V = Vector2<T>;
	C c(V(T(0), T(0)), T(2));
	CHECK_APPROX(c.getDistanceTo(V(T(3), T(4))), T(3));
	CHECK(c.getDistanceTo(V(T(1), T(0))) == T(0));
	CHECK_APPROX(c.getSignedDistanceTo(V(T(1), T(0))), T(-1));
	CHECK_APPROX(c.getSignedDistanceTo(V(T(3), T(4))), T(3));

	CHECK(c.contains(V(T(1), T(1))));
	CHECK(c.contains(V(T(2), T(0))));
	CHECK_FALSE(c.contains(V(T(2), T(1))));

	CHECK(c.intersects(C(V(T(3), T(0)), T(1.5))));
	CHECK(c.intersects(C(V(T(4), T(0)), T(2)))); // touching
	CHECK_FALSE(c.intersects(C(V(T(5), T(0)), T(2))));

	CHECK(c.intersects(AxisAlignedRectangle<T>(V(T(1), T(1)), V(T(3), T(3)))));
	CHECK(c.intersects(AxisAlignedRectangle<T>(V(T(-5), T(-5)), V(T(5), T(5))))); // containing rectangle
	CHECK_FALSE(c.intersects(AxisAlignedRectangle<T>(V(T(1.5), T(1.5)), V(T(3), T(3)))));
}

TEST_CASE_TEMPLATE("streaming and hashing", T, TEST_REAL_TYPES)
{
	using C = Circle2<T>;
	C c(Vector2<T>(T(1), T(2)), T(3.5));
	std::stringstream s;
	s << c;
	CHECK(s.str() == "1 2 3.5");
	C read;
	s >> read;
	CHECK(read == c);
	std::unordered_set<C> set{ c, c, C() };
	CHECK(set.size() == 2);
}

} // TEST_SUITE

TEST_SUITE("AxisAlignedRectangle")
{

TEST_CASE_TEMPLATE("construction and properties", T, TEST_REAL_TYPES)
{
	using R = AxisAlignedRectangle<T>;
	using V = Vector2<T>;
	CHECK(R().isEmpty());
	CHECK(R() == R::EMPTY);
	R r(V(T(1), T(2)), V(T(5), T(4)));
	CHECK_FALSE(r.isEmpty());
	CHECK(R(V(T(4), T(2))) == R(V(T(-2), T(-1)), V(T(2), T(1))));
	CHECK(R(std::pair<V, V>(V(T(1), T(2)), V(T(5), T(4)))) == r);
	CHECK(R(std::tuple<V, V>(V(T(1), T(2)), V(T(5), T(4)))) == r);

	CHECK(r.getMinimum() == V(T(1), T(2)));
	CHECK(r.getMaximum() == V(T(5), T(4)));
	CHECK(r.getDimensions() == V(T(4), T(2)));
	CHECK(r.getHalfDimensions() == V(T(2), T(1)));
	CHECK(r.getCenter() == V(T(3), T(3)));
	CHECK_APPROX(r.getDiagonal(), std::sqrt(T(20)));
	CHECK(r.getPerimeter() == T(12));
	CHECK(r.getArea() == T(8));

	R s = r;
	s.setDimensions(V(T(2), T(6)));
	CHECK(s == R(V(T(2), T(0)), V(T(4), T(6))));
	s.setHalfDimensions(V(T(1), T(1)));
	CHECK(s == R(V(T(2), T(2)), V(T(4), T(4))));
	s.setCenter(V::ZERO);
	CHECK(s == R(V(T(-1)), V(T(1))));
	s.setMinimum(V(T(-2)));
	s.setMaximum(V(T(2)));
	CHECK(s == R(V(T(-2)), V(T(2))));
	CHECK(R(r).set(V::ZERO, V::ONE) == R(V::ZERO, V::ONE));
	CHECK(R(r).clear().isEmpty());

	CHECK(r.isFinite());
	CHECK_FALSE(R::EMPTY.isFinite());
	CHECK(r.approxEquals(R(V(T(1), T(2)), V(T(5), T(4) + Constants<T>::TOLERANCE*T(0.5)))));
	CHECK(r.approxEquals(R(V(T(1), T(2)), V(T(5), T(4.05))), T(0.1)));
	CHECK_FALSE(r.approxEquals(R(V(T(1), T(2)), V(T(5), T(4.05)))));
	CHECK(r != s);
}

TEST_CASE_TEMPLATE("vertices and circumscribed circle", T, TEST_REAL_TYPES)
{
	using R = AxisAlignedRectangle<T>;
	using V = Vector2<T>;
	R r(V(T(1), T(2)), V(T(5), T(4)));
	std::array<V, 4> expected{ V(T(1), T(2)), V(T(5), T(2)), V(T(1), T(4)), V(T(5), T(4)) };
	CHECK(r.getVertices() == expected);
	std::vector<V> copied(4);
	CHECK(r.copyVertices(copied.begin()) == copied.end());
	CHECK(std::equal(copied.begin(), copied.end(), expected.begin()));

	Circle2<T> c = r.getCircumscribedCircle();
	CHECK(c.center == V(T(3), T(3)));
	CHECK_APPROX(c.radius, std::sqrt(T(5)));
	for (const V& v : expected)
		CHECK_APPROX(distance(v, c.center), c.radius);
}

TEST_CASE_TEMPLATE("transformation, union and intersection", T, TEST_REAL_TYPES)
{
	using R = AxisAlignedRectangle<T>;
	using V = Vector2<T>;
	R r(V(T(0), T(0)), V(T(2), T(2)));
	CHECK(R(r).inflate(V(T(1), T(2))) == R(V(T(-1), T(-2)), V(T(3), T(4))));
	CHECK(R(r).translate(V(T(1), T(-1))) == R(V(T(1), T(-1)), V(T(3), T(1))));
	CHECK(R(r).scaleAroundCenter(T(2)) == R(V(T(-1)), V(T(3))));

	R s(V(T(1), T(1)), V(T(3), T(4)));
	CHECK(R::makeUnion(r, s) == R(V(T(0), T(0)), V(T(3), T(4))));
	CHECK(R::makeIntersection(r, s) == R(V(T(1), T(1)), V(T(2), T(2))));
	CHECK(R::makeIntersection(r, R(V(T(5)), V(T(6)))).isEmpty());

	R e;
	e.extendBy(V(T(1), T(1)));
	e.extendBy(V(T(-1), T(3)));
	CHECK(e == R(V(T(-1), T(1)), V(T(1), T(3))));
	e.extendBy(R(V(T(0), T(0)), V(T(2), T(2))));
	CHECK(e == R(V(T(-1), T(0)), V(T(2), T(3))));
	e.extendBy(Circle2<T>(V(T(5), T(5)), T(1)));
	CHECK(e == R(V(T(-1), T(0)), V(T(6), T(6))));
}

TEST_CASE_TEMPLATE("closest point, containment and intersection", T, TEST_REAL_TYPES)
{
	using R = AxisAlignedRectangle<T>;
	using V = Vector2<T>;
	R r(V(T(0), T(0)), V(T(2), T(2)));
	CHECK(r.getClosestPoint(V(T(1), T(1))) == V(T(1), T(1)));
	CHECK(r.getClosestPoint(V(T(5), T(-1))) == V(T(2), T(0)));
	CHECK_APPROX(r.getDistanceTo(V(T(5), T(6))), T(5));

	CHECK(r.contains(V(T(1), T(1))));
	CHECK(r.contains(V(T(2), T(0))));
	CHECK_FALSE(r.contains(V(T(2.1), T(0))));
	CHECK(r.contains(R(V(T(0.5)), V(T(1.5)))));
	CHECK_FALSE(r.contains(R(V(T(0.5)), V(T(2.5)))));
	CHECK(r.contains(Circle2<T>(V(T(1), T(1)), T(1))));
	CHECK_FALSE(r.contains(Circle2<T>(V(T(1), T(1)), T(1.5))));

	CHECK(r.intersects(R(V(T(1)), V(T(3)))));
	CHECK(r.intersects(R(V(T(2)), V(T(3))))); // touching corner
	CHECK_FALSE(r.intersects(R(V(T(2.5)), V(T(3)))));
	CHECK(r.intersects(Circle2<T>(V(T(3), T(1)), T(1.5))));
	CHECK_FALSE(r.intersects(Circle2<T>(V(T(3), T(3)), T(1))));
}

TEST_CASE_TEMPLATE("streaming and hashing", T, TEST_REAL_TYPES)
{
	using R = AxisAlignedRectangle<T>;
	R r(Vector2<T>(T(1), T(2)), Vector2<T>(T(5), T(4)));
	std::stringstream s;
	s << r;
	CHECK(s.str() == "1 2 5 4");
	R read;
	s >> read;
	CHECK(read == r);
	std::unordered_set<R> set{ r, r, R() };
	CHECK(set.size() == 2);
}

} // TEST_SUITE

} // namespace mathematics::templates
