#include "TestUtilities.hpp"
#include <unordered_set>

namespace mathematics::templates {

TEST_SUITE("Plane")
{

TEST_CASE_TEMPLATE("construction and element access", T, TEST_REAL_TYPES)
{
	using P = Plane<T>;
	using V = Vector3<T>;
	CHECK(P().isZero());
	P p(T(0), T(1), T(0), T(-2));
	CHECK(p.a == T(0));
	CHECK(p.b == T(1));
	CHECK(p.c == T(0));
	CHECK(p.d == T(-2));
	CHECK(P(V::UNIT_Y, T(-2)) == p);
	CHECK(P(V::UNIT_Y, V(T(5), T(2), T(-3))) == p);
	CHECK(P(std::tuple<T, T, T, T>(T(0), T(1), T(0), T(-2))) == p);
	CHECK(P(HalfSpace<T>(T(0), T(1), T(0), T(-2))) == p);
	CHECK(p.asHalfSpace() == HalfSpace<T>(T(0), T(1), T(0), T(-2)));

	// Counter-clockwise vertices define the normal direction.
	P q(V(T(0), T(2), T(0)), V(T(0), T(2), T(1)), V(T(1), T(2), T(0)));
	CHECK_APPROX(q, p);

	std::tuple<T, T, T, T> t = p;
	CHECK(std::get<3>(t) == T(-2));
	CHECK(p.template get<1>() == T(1));
	CHECK(get<3>(p) == T(-2));
	auto [a, b, c, d] = p;
	CHECK(a + b + c + d == T(-1));

	CHECK(p.getNormal() == V::UNIT_Y);
	CHECK(p.getConstant() == T(-2));
	P r = p;
	r.setNormal(V::UNIT_X);
	r.setConstant(T(3));
	CHECK(r == P(V::UNIT_X, T(3)));
	CHECK(P(p).setZero().isZero());
	CHECK(P().set(V::UNIT_Z, T(1)) == P(T(0), T(0), T(1), T(1)));
	CHECK(P().set(T(1), T(2), T(3), T(4)) == P(T(1), T(2), T(3), T(4)));
	CHECK(P(p).flip() == P(T(0), T(-1), T(0), T(2)));

	CHECK(p.isFinite());
	CHECK_FALSE(P(T(0), T(1), T(0), std::numeric_limits<T>::infinity()).isFinite());
	CHECK(p.approxEquals(P(T(0), T(1), T(0), T(-2) + Constants<T>::TOLERANCE*T(0.5))));
	CHECK(p.approxEquals(P(T(0), T(1), T(0), T(-2.05)), T(0.1)));
	CHECK_FALSE(p.approxEquals(P(T(0), T(1), T(0), T(-2.05))));
	CHECK(p != r);
}

TEST_CASE_TEMPLATE("normalization and distance", T, TEST_REAL_TYPES)
{
	using P = Plane<T>;
	using V = Vector3<T>;
	P p(T(0), T(3), T(4), T(-10));
	P n = normalize(p);
	CHECK_APPROX(n, P(T(0), T(0.6), T(0.8), T(-2)));
	CHECK_APPROX(normalize(P(p)), n); // rvalue overload
	CHECK_APPROX(P(p).normalize(), n);
	CHECK(P().normalize().isZero());

	V point(T(1), T(5), T(5)); // 3*5 + 4*5 - 10 = 25
	CHECK_APPROX(n.getSignedDistanceTo(point), T(5));
	CHECK_APPROX(n.getDistanceTo(point), T(5));
	CHECK_APPROX(n.getSignedDistanceTo(-point), T(-9));
	CHECK_APPROX(n.getDistanceTo(-point), T(9));
	CHECK_APPROX(n.template getSignedDistanceTo<Normalized>(point), T(5));
	CHECK_APPROX(n.template getDistanceTo<Normalized>(-point), T(9));
	CHECK_APPROX(p.template getSignedDistanceTo<Unnormalized>(point), T(5));
	CHECK_APPROX(p.template getDistanceTo<Unnormalized>(-point), T(9));

	CHECK(n.contains(V(T(7), T(0), T(2.5))));
	CHECK_FALSE(n.contains(point));
}

TEST_CASE_TEMPLATE("transformation", T, TEST_REAL_TYPES)
{
	using P = Plane<T>;
	using V = Vector3<T>;
	P p(normalize(V(T(1), T(2), T(2))), T(-3));
	const V onPlane[] = { p.getNormal()*T(3), p.getNormal()*T(3) + V(T(2), T(-1), T(0)), p.getNormal()*T(3) + V(T(0), T(1), T(-1)) };
	for (const V& v : onPlane)
		REQUIRE(std::fabs(p.getSignedDistanceTo(v)) < T(10)*test::tolerance<T>());

	SUBCASE("translation")
	{
		P t = P(p).translate(V(T(1), T(0), T(0)));
		for (const V& v : onPlane)
			CHECK_APPROX(t.getSignedDistanceTo(v + V(T(1), T(0), T(0))), T(0));
		CHECK_APPROX(t.getNormal(), p.getNormal());
	}

	SUBCASE("orthogonal matrix")
	{
		Matrix3<T> m = Matrix3<T>::makeRotation(V(T(1), T(-1), T(0.5)), T(0.8));
		P t = P(p).transform(m, true);
		for (const V& v : onPlane)
			CHECK_APPROX(t.getSignedDistanceTo(v*m), T(0));
		CHECK_APPROX(t.getNormal(), p.getNormal()*m);
	}

	SUBCASE("general matrix")
	{
		Matrix3<T> m = Matrix3<T>::makeScalingRotation(V(T(2), T(0.5), T(3)), V(T(0), T(1), T(1)), T(0.6));
		P t = P(p).transform(m);
		CHECK_APPROX(t.getNormal().getMagnitude(), T(1));
		for (const V& v : onPlane)
			CHECK_APPROX(t.getSignedDistanceTo(v*m), T(0));
		CHECK(t.getSignedDistanceTo((p.getNormal()*T(4))*m) > T(0)); // the positive side is preserved
	}

	SUBCASE("affine transformation")
	{
		AffineTransform<T> a = AffineTransform<T>::makeScalingRotationTranslation(V(T(2), T(0.5), T(3)), Axis::Y, T(0.6),
			V(T(1), T(-2), T(4)));
		P t = P(p).transform(a);
		for (const V& v : onPlane)
			CHECK_APPROX(t.getSignedDistanceTo(transform(v, a)), T(0));

		AffineTransform<T> r = AffineTransform<T>::makeRotationTranslation(Axis::X, T(1.1), V(T(1), T(-2), T(4)));
		P u = P(p).transform(r, true);
		for (const V& v : onPlane)
			CHECK_APPROX(u.getSignedDistanceTo(transform(v, r)), T(0));
	}
}

TEST_CASE_TEMPLATE("intersections", T, TEST_REAL_TYPES)
{
	using P = Plane<T>;
	using V = Vector3<T>;
	P p(V::UNIT_Y, T(-1)); // y = 1

	CHECK(p.intersects(Triangle3<T>(V(T(0), T(0), T(0)), V(T(1), T(2), T(0)), V(T(0), T(0), T(1)))));
	CHECK_FALSE(p.intersects(Triangle3<T>(V(T(0), T(2), T(0)), V(T(1), T(3), T(0)), V(T(0), T(2), T(1)))));
	CHECK(p.intersects(AxisAlignedBox<T>(V(T(0)), V(T(2)))));
	CHECK_FALSE(p.intersects(AxisAlignedBox<T>(V(T(0), T(1.5), T(0)), V(T(2)))));
	CHECK(p.intersects(OrientedBox<T>(V(T(0), T(1.9), T(0)), Matrix3<T>::makeRotation(Axis::Z, Constants<T>::PI/T(4)), V(T(1)))));
	CHECK_FALSE(p.intersects(OrientedBox<T>(V(T(0), T(2.5), T(0)), Matrix3<T>::makeRotation(Axis::Z, Constants<T>::PI/T(4)), V(T(1)))));
	CHECK(p.intersects(Sphere<T>(V(T(0), T(1.5), T(0)), T(1))));
	CHECK_FALSE(p.intersects(Sphere<T>(V(T(0), T(3), T(0)), T(1))));
	CHECK(p.intersects(Ellipsoid<T>(V(T(0), T(3), T(0)), Matrix3<T>::IDENTITY, V(T(1), T(3), T(1))))); // spans 0 <= y <= 6
	CHECK_FALSE(p.intersects(Ellipsoid<T>(V(T(0), T(3), T(0)), Matrix3<T>::IDENTITY, V(T(3), T(1), T(1))))); // spans 2 <= y <= 4
	CHECK(p.intersects(Ellipsoid<T>(V(T(0), T(3), T(0)), Matrix3<T>::makeRotation(Axis::Z, Constants<T>::HALF_PI), V(T(3), T(1), T(1))))); // rotated
}

TEST_CASE_TEMPLATE("constants, streaming and hashing", T, TEST_REAL_TYPES)
{
	using P = Plane<T>;
	CHECK(P::EMPTY.isZero());
	P p(T(0), T(1), T(0), T(-2.5));
	std::stringstream s;
	s << p;
	CHECK(s.str() == "0 1 0 -2.5");
	P read;
	s >> read;
	CHECK(read == p);
	std::unordered_set<P> set{ p, p, P() };
	CHECK(set.size() == 2);
}

} // TEST_SUITE

TEST_SUITE("HalfSpace")
{

TEST_CASE_TEMPLATE("construction and properties", T, TEST_REAL_TYPES)
{
	using H = HalfSpace<T>;
	using V = Vector3<T>;
	CHECK(H().isZero());
	H h(T(0), T(1), T(0), T(-2));
	CHECK(H(V::UNIT_Y, T(-2)) == h);
	CHECK(H(V::UNIT_Y, V(T(5), T(2), T(-3))) == h);
	CHECK(H(std::tuple<T, T, T, T>(T(0), T(1), T(0), T(-2))) == h);
	CHECK(H(Plane<T>(T(0), T(1), T(0), T(-2))) == h);
	CHECK(h.asPlane() == Plane<T>(T(0), T(1), T(0), T(-2)));
	CHECK_APPROX(H(V(T(0), T(2), T(0)), V(T(0), T(2), T(1)), V(T(1), T(2), T(0))), h);

	CHECK(h.template get<1>() == T(1));
	CHECK(get<3>(h) == T(-2));
	CHECK(h.getNormal() == V::UNIT_Y);
	CHECK(h.getConstant() == T(-2));
	H r = h;
	r.setNormal(V::UNIT_X);
	r.setConstant(T(3));
	CHECK(r == H(V::UNIT_X, T(3)));
	CHECK(H(h).setZero().isZero());
	CHECK(H().set(V::UNIT_Z, T(1)) == H(T(0), T(0), T(1), T(1)));
	CHECK(H().set(T(1), T(2), T(3), T(4)) == H(T(1), T(2), T(3), T(4)));
	CHECK(H(h).flip() == H(T(0), T(-1), T(0), T(2)));
	CHECK(h.isFinite());
	CHECK(h.approxEquals(H(T(0), T(1), T(0), T(-2.05)), T(0.1)));
	CHECK_FALSE(h.approxEquals(H(T(0), T(1), T(0), T(-2.05))));
	CHECK(H::EMPTY.isZero());

	CHECK_APPROX(normalize(H(T(0), T(2), T(0), T(-4))), h);
	CHECK_APPROX(normalize(H(H(T(0), T(2), T(0), T(-4)))), h);
	CHECK_APPROX(H(T(0), T(2), T(0), T(-4)).normalize(), h);
}

TEST_CASE_TEMPLATE("distance and containment", T, TEST_REAL_TYPES)
{
	using H = HalfSpace<T>;
	using V = Vector3<T>;
	H h(V::UNIT_Y, T(-2)); // the inside is y <= 2
	CHECK(h.contains(V(T(5), T(1), T(-3))));
	CHECK(h.contains(V(T(5), T(2), T(-3))));
	CHECK_FALSE(h.contains(V(T(5), T(3), T(-3))));
	CHECK_APPROX(h.getSignedDistanceTo(V(T(0), T(5), T(0))), T(3));
	CHECK_APPROX(h.getSignedDistanceTo(V(T(0), T(-1), T(0))), T(-3));
	CHECK_APPROX(h.getDistanceTo(V(T(0), T(5), T(0))), T(3));
	CHECK(h.getDistanceTo(V(T(0), T(-1), T(0))) == T(0)); // inside

	H u(V(T(0), T(2), T(0)), T(-4));
	CHECK_APPROX(u.template getSignedDistanceTo<Unnormalized>(V(T(0), T(5), T(0))), T(3));
	CHECK_APPROX(u.template getDistanceTo<Unnormalized>(V(T(0), T(5), T(0))), T(3));
	CHECK(u.template getDistanceTo<Unnormalized>(V(T(0), T(-5), T(0))) == T(0));
	CHECK_APPROX(h.template getSignedDistanceTo<Normalized>(V(T(0), T(5), T(0))), T(3));
	CHECK_APPROX(h.template getDistanceTo<Normalized>(V(T(0), T(5), T(0))), T(3));
}

TEST_CASE_TEMPLATE("transformation", T, TEST_REAL_TYPES)
{
	using H = HalfSpace<T>;
	using V = Vector3<T>;
	H h(normalize(V(T(1), T(2), T(2))), T(-3));
	V inside = V::ZERO, outside = h.getNormal()*T(5);
	REQUIRE(h.contains(inside));
	REQUIRE_FALSE(h.contains(outside));

	CHECK_APPROX(H(h).translate(V(T(1), T(0), T(0))), H(h.getNormal(), h.d - h.getNormal().x));

	Matrix3<T> m = Matrix3<T>::makeScalingRotation(V(T(2), T(0.5), T(3)), V(T(0), T(1), T(1)), T(0.6));
	H t = H(h).transform(m);
	CHECK(t.contains(inside*m));
	CHECK_FALSE(t.contains(outside*m));
	CHECK_APPROX(t.getSignedDistanceTo((h.getNormal()*T(3))*m), T(0));

	Matrix3<T> r = Matrix3<T>::makeRotation(Axis::Z, T(0.6));
	CHECK_APPROX(H(h).transform(r, true).getSignedDistanceTo((h.getNormal()*T(3))*r), T(0));

	AffineTransform<T> a = AffineTransform<T>::makeScalingRotationTranslation(V(T(2), T(0.5), T(3)), Axis::Y, T(0.6),
		V(T(1), T(-2), T(4)));
	H ta = H(h).transform(a);
	CHECK(ta.contains(transform(inside, a)));
	CHECK_FALSE(ta.contains(transform(outside, a)));
	CHECK_APPROX(ta.getSignedDistanceTo(transform(h.getNormal()*T(3), a)), T(0));

	AffineTransform<T> ra = AffineTransform<T>::makeRotationTranslation(Axis::X, T(1.1), V(T(1), T(-2), T(4)));
	CHECK_APPROX(H(h).transform(ra, true).getSignedDistanceTo(transform(h.getNormal()*T(3), ra)), T(0));
}

TEST_CASE_TEMPLATE("intersections", T, TEST_REAL_TYPES)
{
	using H = HalfSpace<T>;
	using V = Vector3<T>;
	H h(V::UNIT_Y, T(-1)); // the inside is y <= 1

	CHECK(h.intersects(Triangle3<T>(V(T(0), T(0), T(0)), V(T(1), T(0), T(0)), V(T(0), T(0), T(1))))); // fully inside
	CHECK(h.intersects(Triangle3<T>(V(T(0), T(0), T(0)), V(T(1), T(2), T(0)), V(T(0), T(0), T(1))))); // straddling
	CHECK_FALSE(h.intersects(Triangle3<T>(V(T(0), T(2), T(0)), V(T(1), T(3), T(0)), V(T(0), T(2), T(1)))));
	CHECK(h.intersects(AxisAlignedBox<T>(V(T(0)), V(T(0.5)))));
	CHECK(h.intersects(AxisAlignedBox<T>(V(T(0)), V(T(2)))));
	CHECK_FALSE(h.intersects(AxisAlignedBox<T>(V(T(0), T(1.5), T(0)), V(T(2)))));
	CHECK(h.intersects(OrientedBox<T>(V(T(0), T(1.9), T(0)), Matrix3<T>::makeRotation(Axis::Z, Constants<T>::PI/T(4)), V(T(1)))));
	CHECK_FALSE(h.intersects(OrientedBox<T>(V(T(0), T(2.5), T(0)), Matrix3<T>::makeRotation(Axis::Z, Constants<T>::PI/T(4)), V(T(1)))));
	CHECK(h.intersects(Sphere<T>(V(T(0), T(-3), T(0)), T(1))));
	CHECK(h.intersects(Sphere<T>(V(T(0), T(1.5), T(0)), T(1))));
	CHECK_FALSE(h.intersects(Sphere<T>(V(T(0), T(3), T(0)), T(1))));
}

TEST_CASE_TEMPLATE("streaming and hashing", T, TEST_REAL_TYPES)
{
	using H = HalfSpace<T>;
	H h(T(0), T(1), T(0), T(-2.5));
	std::stringstream s;
	s << h;
	CHECK(s.str() == "0 1 0 -2.5");
	H read;
	s >> read;
	CHECK(read == h);
	std::unordered_set<H> set{ h, h, H() };
	CHECK(set.size() == 2);
}

} // TEST_SUITE

} // namespace mathematics::templates
