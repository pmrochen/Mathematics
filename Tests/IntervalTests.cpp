#include "TestUtilities.hpp"
#include <unordered_set>

namespace mathematics::templates {

TEST_SUITE("Interval")
{

TEST_CASE_TEMPLATE("construction", T, TEST_REAL_TYPES)
{
	SUBCASE("default is empty")
	{
		Interval<T> i;
		CHECK(i.isEmpty());
		CHECK(i == Interval<T>::EMPTY);
	}

	SUBCASE("from bounds, value, pair and tuple")
	{
		Interval<T> a(T(1), T(3));
		CHECK(a.minimum == T(1));
		CHECK(a.maximum == T(3));
		CHECK_FALSE(a.isEmpty());

		Interval<T> b(T(2));
		CHECK(b.isPoint());
		CHECK(b.getLength() == T(0));

		CHECK(Interval<T>(std::pair<T, T>(T(1), T(3))) == a);
		CHECK(Interval<T>(std::tuple<T, T>(T(1), T(3))) == a);
	}

	SUBCASE("zero")
	{
		CHECK(Interval<T>::ZERO.isZero());
		CHECK(Interval<T>(Constants<T>::TOLERANCE*T(0.5), -Constants<T>::TOLERANCE*T(0.5)).isApproxZero());
		CHECK(Interval<T>(T(1), T(2)).setZero().isZero());
	}
}

TEST_CASE_TEMPLATE("properties", T, TEST_REAL_TYPES)
{
	Interval<T> i(T(2), T(6));
	CHECK(i.getLength() == T(4));
	CHECK(i.getCenter() == T(4));
	CHECK(i.getMinimum() == T(2));
	CHECK(i.getMaximum() == T(6));
	CHECK(i.isFinite());
	CHECK_FALSE(Interval<T>::EMPTY.isFinite());

	i.setLength(T(2));
	CHECK(i == Interval<T>(T(3), T(5)));
	i.setCenter(T(0));
	CHECK(i == Interval<T>(T(-1), T(1)));
	i.setMinimum(T(-2));
	i.setMaximum(T(7));
	CHECK(i == Interval<T>(T(-2), T(7)));
	CHECK(i.set(T(4)).isPoint());
	CHECK(i.set(T(1), T(2)) == Interval<T>(T(1), T(2)));
	CHECK(i.clear().isEmpty());
}

TEST_CASE_TEMPLATE("approximate comparison", T, TEST_REAL_TYPES)
{
	Interval<T> a(T(1), T(2));
	CHECK(a.approxEquals(Interval<T>(T(1), T(2) + Constants<T>::TOLERANCE*T(0.5))));
	CHECK_FALSE(a.approxEquals(Interval<T>(T(1), T(2.1))));
	CHECK(a.approxEquals(Interval<T>(T(1.05), T(2)), T(0.1)));
	CHECK(a != Interval<T>(T(1), T(2.1)));
}

TEST_CASE_TEMPLATE("interpolation and transformation", T, TEST_REAL_TYPES)
{
	Interval<T> i(T(2), T(6));
	CHECK(i.interpolate(T(0)) == T(2));
	CHECK(i.interpolate(T(1)) == T(6));
	CHECK(i.interpolate(T(0.5)) == T(4));

	CHECK(Interval<T>(i).translate(T(1)) == Interval<T>(T(3), T(7)));
	CHECK(Interval<T>(i).scale(T(2)) == Interval<T>(T(4), T(12)));
	CHECK(Interval<T>(i).scaleAroundCenter(T(0.5)) == Interval<T>(T(3), T(5)));
}

TEST_CASE_TEMPLATE("union and intersection", T, TEST_REAL_TYPES)
{
	Interval<T> a(T(0), T(2));
	Interval<T> b(T(1), T(3));
	Interval<T> c(T(4), T(5));

	CHECK(Interval<T>::makeUnion(a, b) == Interval<T>(T(0), T(3)));
	CHECK(Interval<T>::makeIntersection(a, b) == Interval<T>(T(1), T(2)));
	CHECK(Interval<T>::makeIntersection(a, c).isEmpty());

	Interval<T> e;
	e.extendBy(T(3));
	CHECK(e.isPoint());
	e.extendBy(T(-1));
	CHECK(e == Interval<T>(T(-1), T(3)));
	e.extendBy(c);
	CHECK(e == Interval<T>(T(-1), T(5)));
}

TEST_CASE_TEMPLATE("containment and intersection", T, TEST_REAL_TYPES)
{
	Interval<T> a(T(0), T(2));
	CHECK(a.contains(T(0)));
	CHECK(a.contains(T(2)));
	CHECK(a.contains(T(1)));
	CHECK_FALSE(a.contains(T(2.5)));
	CHECK(a.contains(Interval<T>(T(0.5), T(1.5))));
	CHECK_FALSE(a.contains(Interval<T>(T(0.5), T(2.5))));

	CHECK(a.intersects(Interval<T>(T(1), T(3))));
	CHECK(a.intersects(Interval<T>(T(2), T(3)))); // touching
	CHECK_FALSE(a.intersects(Interval<T>(T(2.5), T(3))));

	SUBCASE("findIntersection")
	{
		auto overlap = a.findIntersection(Interval<T>(T(1), T(3)));
		REQUIRE(overlap.has_value());
		CHECK(*overlap == Interval<T>(T(1), T(2)));

		auto inside = a.findIntersection(Interval<T>(T(0.5), T(1.5)));
		REQUIRE(inside.has_value());
		CHECK(*inside == Interval<T>(T(0.5), T(1.5)));

		auto touchingAbove = a.findIntersection(Interval<T>(T(2), T(3)));
		REQUIRE(touchingAbove.has_value());
		CHECK(*touchingAbove == Interval<T>(T(2)));

		auto touchingBelow = a.findIntersection(Interval<T>(T(-1), T(0)));
		REQUIRE(touchingBelow.has_value());
		CHECK(*touchingBelow == Interval<T>(T(0)));

		CHECK_FALSE(a.findIntersection(Interval<T>(T(3), T(4))).has_value());
	}
}

TEST_CASE_TEMPLATE("streaming and hashing", T, TEST_REAL_TYPES)
{
	Interval<T> a(T(1.5), T(-2));
	std::stringstream s;
	s << a;
	CHECK(s.str() == "1.5 -2");

	Interval<T> b(Uninitialized{});
	s >> b;
	CHECK(b == a);

	std::unordered_set<Interval<T>> set{ a, a, Interval<T>(T(0), T(1)) };
	CHECK(set.size() == 2);
}

} // TEST_SUITE

} // namespace mathematics::templates
