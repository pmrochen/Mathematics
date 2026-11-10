#include "TestUtilities.hpp"
#include <unordered_set>
#include <vector>
#include <iterator>
#include <cstdint>

namespace mathematics::templates {

namespace {

// Frustum looking along +z from the origin, 2x1 at depth 1 and 8x4 at depth 4
template<typename T>
SymmetricFrustum<T> makeFrustum()
{
	return SymmetricFrustum<T>(Vector3<T>::ZERO, Matrix3<T>::IDENTITY, Vector2<T>(T(1), T(0.5)), T(1), T(4));
}

template<typename T>
bool containsApprox(const Sphere<T>& s, const Vector3<T>& p)
{
	return distance(s.center, p) <= s.radius*(T(1) + test::tolerance<T>());
}

} // namespace

TEST_SUITE("Sphere")
{

TEST_CASE_TEMPLATE("construction and properties", T, TEST_REAL_TYPES)
{
	using S = Sphere<T>;
	using V = Vector3<T>;
	CHECK(S() == S(V::ZERO, T(0)));
	S s(V(T(1), T(2), T(3)), T(2));
	CHECK(s.getCenter() == V(T(1), T(2), T(3)));
	CHECK(s.getRadius() == T(2));
	CHECK(s.getDiameter() == T(4));
	CHECK_APPROX(s.getSurfaceArea(), T(16)*Constants<T>::PI);
	CHECK_APPROX(s.getVolume(), T(32)/T(3)*Constants<T>::PI);
	S t = s;
	t.setCenter(V::ZERO);
	t.setRadius(T(1));
	CHECK(t == S(V::ZERO, T(1)));
	t.setDiameter(T(6));
	CHECK(t.radius == T(3));
	CHECK(S().set(V::ONE, T(1)) == S(V::ONE, T(1)));
	CHECK(S(Ellipsoid<T>(V::ONE, Matrix3<T>::IDENTITY, V(T(1), T(3), T(2)))) == S(V::ONE, T(3)));
	CHECK(s.isFinite());
	CHECK_FALSE(S(V::ZERO, std::numeric_limits<T>::infinity()).isFinite());
	CHECK(s.approxEquals(S(s.center, T(2.05)), T(0.1)));
	CHECK_FALSE(s.approxEquals(S(s.center, T(2.05))));
	CHECK(s != t);
	CHECK(s.getCircumscribedBox() == AxisAlignedBox<T>(V(T(-1), T(0), T(1)), V(T(3), T(4), T(5))));
}

TEST_CASE_TEMPLATE("transformation", T, TEST_REAL_TYPES)
{
	using S = Sphere<T>;
	using V = Vector3<T>;
	S s(V(T(1), T(0), T(0)), T(2));
	CHECK(S(s).translate(V::ONE) == S(V(T(2), T(1), T(1)), T(2)));

	Matrix3<T> r = Matrix3<T>::makeRotation(Axis::Z, Constants<T>::HALF_PI);
	CHECK_APPROX(S(s).transform(r, true), S(V(T(0), T(1), T(0)), T(2)));
	CHECK_APPROX(S(s).transform(Matrix3<T>::makeScaling(V(T(1), T(3), T(1)))), S(V(T(1), T(0), T(0)), T(6))); // bounding sphere
	AffineTransform<T> a(r, V(T(0), T(0), T(5)));
	CHECK_APPROX(S(s).transform(a, true), S(V(T(0), T(1), T(5)), T(2)));
	CHECK_APPROX(S(s).transform(AffineTransform<T>::makeScalingTranslation(V(T(2)), V::ONE)), S(V(T(3), T(1), T(1)), T(4)));
}

TEST_CASE_TEMPLATE("distance, containment and intersections", T, TEST_REAL_TYPES)
{
	using S = Sphere<T>;
	using V = Vector3<T>;
	S s(V::ZERO, T(2));
	CHECK_APPROX(s.getDistanceTo(V(T(0), T(3), T(4))), T(3));
	CHECK(s.getDistanceTo(V::ONE) == T(0));
	CHECK_APPROX(s.getSignedDistanceTo(V(T(1), T(0), T(0))), T(-1));
	CHECK(s.contains(V::ONE));
	CHECK(s.contains(V(T(2), T(0), T(0))));
	CHECK_FALSE(s.contains(V(T(2), T(1), T(0))));

	CHECK(s.intersects(HalfSpace<T>(V::UNIT_X, T(-3))));	// x <= 3 contains the sphere
	CHECK(s.intersects(HalfSpace<T>(V::UNIT_X, T(1))));		// x <= -1 cuts the sphere
	CHECK_FALSE(s.intersects(HalfSpace<T>(V::UNIT_X, T(3)))); // x <= -3
	CHECK(s.intersects(Plane<T>(V::UNIT_X, T(-1.5))));
	CHECK_FALSE(s.intersects(Plane<T>(V::UNIT_X, T(-3))));
	CHECK(s.intersects(Triangle3<T>(V(T(1), T(-5), T(0)), V(T(1), T(5), T(0)), V(T(1), T(0), T(5)))));
	CHECK_FALSE(s.intersects(Triangle3<T>(V(T(3), T(-5), T(0)), V(T(3), T(5), T(0)), V(T(3), T(0), T(5)))));
	CHECK(s.intersects(AxisAlignedBox<T>(V(T(1)), V(T(3)))));
	CHECK_FALSE(s.intersects(AxisAlignedBox<T>(V(T(1.5)), V(T(3)))));
	CHECK(s.intersects(OrientedBox<T>(V(T(2.5), T(0), T(0)), Matrix3<T>::makeRotation(Axis::Z, T(0.7)), V(T(1)))));
	CHECK_FALSE(s.intersects(OrientedBox<T>(V(T(4), T(0), T(0)), Matrix3<T>::makeRotation(Axis::Z, T(0.7)), V(T(1)))));
	CHECK(s.intersects(S(V(T(3), T(0), T(0)), T(1.5))));
	CHECK(s.intersects(S(V(T(3), T(0), T(0)), T(1))));
	CHECK_FALSE(s.intersects(S(V(T(4), T(0), T(0)), T(1))));
	CHECK(s.intersects(Cone<T>(V(T(0), T(0), T(1.5)), V::UNIT_Z, T(2), T(1)))); // the apex is inside the sphere
	CHECK_FALSE(s.intersects(Cone<T>(V(T(0), T(0), T(2.5)), V::UNIT_Z, T(2), T(1)))); // pointing away from the sphere
	CHECK(s.intersects(Cone<T>(V(T(0), T(0), T(5)), V::MINUS_UNIT_Z, T(4), T(2)))); // the base reaches into the sphere
	CHECK(S(V(T(0), T(0), T(2)), T(0.5)).intersects(makeFrustum<T>()));
	CHECK_FALSE(S(V(T(0), T(0), T(-1)), T(0.5)).intersects(makeFrustum<T>()));
}

TEST_CASE_TEMPLATE("streaming and hashing", T, TEST_REAL_TYPES)
{
	using S = Sphere<T>;
	S s(Vector3<T>(T(1), T(2), T(3)), T(4.5));
	std::stringstream str;
	str << s;
	CHECK(str.str() == "1 2 3 4.5");
	S read;
	str >> read;
	CHECK(read == s);
	std::unordered_set<S> set{ s, s, S() };
	CHECK(set.size() == 2);
}

} // TEST_SUITE

TEST_SUITE("Ellipsoid")
{

TEST_CASE_TEMPLATE("construction and properties", T, TEST_REAL_TYPES)
{
	using E = Ellipsoid<T>;
	using V = Vector3<T>;
	using M = Matrix3<T>;
	CHECK(E().basis.isIdentity());
	M r = M::makeRotation(Axis::Z, Constants<T>::HALF_PI);
	E e(V(T(1), T(2), T(3)), r, V(T(3), T(2), T(1)));
	CHECK(e.getCenter() == V(T(1), T(2), T(3)));
	CHECK(e.getBasis() == r);
	CHECK(e.getRadii() == V(T(3), T(2), T(1)));
	CHECK(e.getDiameters() == V(T(6), T(4), T(2)));
	CHECK_APPROX(e.getVolume(), T(8)*Constants<T>::PI);
	CHECK(e.getCircumscribedSphere() == Sphere<T>(e.center, T(3)));

	CHECK(E(Sphere<T>(V::ONE, T(2))) == E(V::ONE, M::IDENTITY, V(T(2))));
	CHECK_APPROX(E(Sphere<T>(V::ONE, T(2)), M::makeScaling(V(T(1), T(2), T(3)))), E(V(T(1), T(2), T(3)), M::IDENTITY, V(T(2), T(4), T(6))));
	CHECK_APPROX(E(Sphere<T>(V::ONE, T(2)), r, true), E(V::ONE*r, r, V(T(2))));
	CHECK_APPROX(E(Sphere<T>(V::ONE, T(2)), AffineTransform<T>::makeScalingTranslation(V(T(2)), V(T(1), T(0), T(0)))),
		E(V(T(3), T(2), T(2)), M::IDENTITY, V(T(4))));
	CHECK_APPROX(E(Sphere<T>(V::ONE, T(2)), AffineTransform<T>(r, V::ZERO), true), E(V::ONE*r, r, V(T(2))));

	E f = e;
	f.setCenter(V::ZERO);
	f.setBasis(M::IDENTITY);
	f.setRadii(V::ONE);
	CHECK(f == E(V::ZERO, M::IDENTITY, V::ONE));
	CHECK(E().set(V::ZERO, M::IDENTITY, V::ONE) == f);
	CHECK(e.isFinite());
	CHECK_FALSE(E(V::INF, M::IDENTITY, V::ONE).isFinite());
	CHECK(e.approxEquals(E(e.center, e.basis, e.radii + V(T(0.05))), T(0.1)));
	CHECK_FALSE(e.approxEquals(E(e.center, e.basis, e.radii + V(T(0.05)))));
	CHECK(e != f);
}

TEST_CASE_TEMPLATE("evaluation and containment", T, TEST_REAL_TYPES)
{
	using E = Ellipsoid<T>;
	using V = Vector3<T>;
	using M = Matrix3<T>;
	// The first axis is rotated onto Y, so the ellipsoid spans 3 units along Y and 2 along X.
	E e(V(T(1), T(2), T(3)), M::makeRotation(Axis::Z, Constants<T>::HALF_PI), V(T(3), T(2), T(1)));
	CHECK_APPROX(e.evaluate(e.center), T(-1));
	CHECK_APPROX(e.evaluate(e.center + V(T(0), T(3), T(0))), T(0));
	CHECK_APPROX(e.evaluate(e.center + V(T(2), T(0), T(0))), T(0));
	CHECK_APPROX(e.evaluate(e.center + V(T(0), T(0), T(2))), T(3));
	CHECK(e.contains(e.center + V(T(0), T(2.9), T(0))));
	CHECK_FALSE(e.contains(e.center + V(T(2.1), T(0), T(0))));

	CHECK_APPROX(e.getMatrix()*e.getInverseMatrix(), M::IDENTITY);
	V d = V(T(0.3), T(-0.2), T(0.4));
	CHECK_APPROX(dot(d, d*e.getMatrix()) - T(1), e.evaluate(e.center + d));
}

TEST_CASE_TEMPLATE("transformation", T, TEST_REAL_TYPES)
{
	using E = Ellipsoid<T>;
	using V = Vector3<T>;
	using M = Matrix3<T>;
	E e(V(T(1), T(0), T(0)), M::IDENTITY, V(T(3), T(2), T(1)));
	CHECK(E(e).translate(V::ONE).center == V(T(2), T(1), T(1)));
	M r = M::makeRotation(Axis::Z, Constants<T>::HALF_PI);
	CHECK_APPROX(E(e).transform(r, true), E(V(T(0), T(1), T(0)), r, e.radii));
	CHECK_APPROX(E(e).transform(M::makeScaling(V(T(2), T(1), T(1)))), E(V(T(2), T(0), T(0)), M::IDENTITY, V(T(6), T(2), T(1))));
	AffineTransform<T> a(r, V(T(0), T(0), T(5)));
	CHECK_APPROX(E(e).transform(a, true), E(V(T(0), T(1), T(5)), r, e.radii));
	CHECK_APPROX(E(e).transform(AffineTransform<T>::makeScaling(T(2))), E(V(T(2), T(0), T(0)), M::IDENTITY, e.radii*T(2)));

	E n(V::ZERO, M::makeScaling(V(T(2), T(3), T(4))), V::ONE);
	CHECK_APPROX(n.orthonormalize(), E(V::ZERO, M::IDENTITY, V(T(2), T(3), T(4))));
}

TEST_CASE_TEMPLATE("streaming and hashing", T, TEST_REAL_TYPES)
{
	using E = Ellipsoid<T>;
	using V = Vector3<T>;
	E e(V(T(1), T(2), T(3)), Matrix3<T>::IDENTITY, V(T(4), T(5), T(6)));
	std::stringstream s;
	s << e;
	E read;
	s >> read;
	CHECK(read == e);
	std::unordered_set<E> set{ e, e, E() };
	CHECK(set.size() == 2);
}

} // TEST_SUITE

TEST_SUITE("Torus")
{

TEST_CASE_TEMPLATE("construction and properties", T, TEST_REAL_TYPES)
{
	using To = Torus<T>;
	using V = Vector3<T>;
	To t(V(T(1), T(2), T(3)), V::UNIT_Y, T(3), T(1));
	CHECK(To(V(T(1), T(2), T(3)), Axis::Y, T(3), T(1)) == t);
	CHECK(t.getCenter() == V(T(1), T(2), T(3)));
	CHECK(t.getAxis() == V::UNIT_Y);
	CHECK(t.getMajorRadius() == T(3));
	CHECK(t.getMinorRadius() == T(1));
	CHECK(t.getInnerRadius() == T(2));
	CHECK(t.getOuterRadius() == T(4));
	CHECK(t.getAspectRatio() == T(3));
	CHECK_APPROX(t.getSurfaceArea(), T(12)*Constants<T>::PI_SQR);
	CHECK_APPROX(t.getVolume(), T(6)*Constants<T>::PI_SQR);

	To u = t;
	u.setCenter(V::ZERO);
	u.setAxis(V::UNIT_Z);
	u.setMajorRadius(T(5));
	u.setMinorRadius(T(2));
	CHECK(u == To(V::ZERO, V::UNIT_Z, T(5), T(2)));
	CHECK(To().set(V::ZERO, V::UNIT_Z, T(5), T(2)) == u);
	CHECK(To().set(V::ZERO, Axis::Z, T(5), T(2)) == u);
	CHECK(To(t).translate(V::ONE).center == V(T(2), T(3), T(4)));
	CHECK(t.isFinite());
	CHECK_FALSE(To(V::ZERO, V::UNIT_Z, std::numeric_limits<T>::infinity(), T(1)).isFinite());
	CHECK(t.approxEquals(To(t.center, t.axis, T(3.05), T(1)), T(0.1)));
	CHECK_FALSE(t.approxEquals(To(t.center, t.axis, T(3.05), T(1))));
	CHECK(t != u);

	OrientedBox<T> box = t.getCircumscribedBox();
	CHECK(box.center == t.center);
	CHECK_APPROX(box.getVolume(), T(8)*T(16)*T(1));
	CHECK(box.contains(t.center + V(T(3.99), T(0), T(0))));
	CHECK(box.contains(t.center + V(T(0), T(0.99), T(-3.99))));
	CHECK_FALSE(box.contains(t.center + V(T(0), T(1.01), T(0))));

	std::stringstream s;
	s << t;
	To read;
	s >> read;
	CHECK(read == t);
	std::unordered_set<To> set{ t, t, u };
	CHECK(set.size() == 2);
}

} // TEST_SUITE

TEST_SUITE("Cylinder")
{

TEST_CASE_TEMPLATE("construction and properties", T, TEST_REAL_TYPES)
{
	using C = Cylinder<T>;
	using V = Vector3<T>;
	C c(V(T(1), T(2), T(3)), V::UNIT_Y, T(4), T(1));
	CHECK(C(V(T(1), T(2), T(3)), Axis::Y, T(4), T(1)) == c);
	CHECK_APPROX(C(V(T(1), T(0), T(3)), V(T(1), T(4), T(3)), T(1)), c);
	CHECK(c.getHeight() == T(4));
	CHECK(c.getHalfHeight() == T(2));
	CHECK(c.getRadius() == T(1));
	CHECK_APPROX(c.getSurfaceArea(), T(10)*Constants<T>::PI);
	CHECK_APPROX(c.getVolume(), T(4)*Constants<T>::PI);

	C d = c;
	d.setCenter(V::ZERO);
	d.setAxis(V::UNIT_Z);
	d.setHeight(T(2));
	d.setRadius(T(3));
	CHECK(d == C(V::ZERO, V::UNIT_Z, T(2), T(3)));
	CHECK(C().set(V::ZERO, V::UNIT_Z, T(2), T(3)) == d);
	CHECK(C().set(V::ZERO, Axis::Z, T(2), T(3)) == d);
	CHECK(C(c).translate(V::ONE).center == V(T(2), T(3), T(4)));
	CHECK_APPROX(C(V::ZERO, V(T(0), T(2), T(0)), T(2), T(1)).normalize(), C(V::ZERO, V::UNIT_Y, T(4), T(1)));
	CHECK(c.isFinite());
	CHECK_FALSE(C(V::ZERO, V::UNIT_Z, T(1), std::numeric_limits<T>::infinity()).isFinite());
	CHECK(c.approxEquals(C(c.center, c.axis, T(4.05), T(1)), T(0.1)));
	CHECK_FALSE(c.approxEquals(C(c.center, c.axis, T(4.05), T(1))));
	CHECK(c != d);

	std::stringstream s;
	s << c;
	C read;
	s >> read;
	CHECK(read == c);
	std::unordered_set<C> set{ c, c, d };
	CHECK(set.size() == 2);
}

TEST_CASE_TEMPLATE("containment and bounding box", T, TEST_REAL_TYPES)
{
	using C = Cylinder<T>;
	using V = Vector3<T>;
	C c(V(T(1), T(2), T(3)), normalize(V(T(1), T(1), T(0))), T(4), T(1));
	V up = c.axis;
	V side = normalize(cross(up, V::UNIT_Z));
	CHECK(c.contains(c.center));
	CHECK(c.contains(c.center + up*T(1.9) + side*T(0.9)));
	CHECK_FALSE(c.contains(c.center + up*T(2.1)));
	CHECK_FALSE(c.contains(c.center + side*T(1.1)));

	OrientedBox<T> box = c.getCircumscribedBox();
	CHECK_APPROX(box.getVolume(), T(16));
	CHECK(box.contains(c.center + up*T(1.99) + side*T(0.99)));
	CHECK(box.contains(c.center - up*T(1.99) + V::UNIT_Z*T(0.99)));
	CHECK_FALSE(box.contains(c.center + up*T(2.01)));
}

} // TEST_SUITE

TEST_SUITE("Capsule")
{

TEST_CASE_TEMPLATE("construction and properties", T, TEST_REAL_TYPES)
{
	using C = Capsule<T>;
	using V = Vector3<T>;
	C c(V(T(1), T(2), T(3)), V::UNIT_Y, T(4), T(1));
	CHECK(C(V(T(1), T(2), T(3)), Axis::Y, T(4), T(1)) == c);
	CHECK_APPROX(C(V(T(1), T(0), T(3)), V(T(1), T(4), T(3)), T(1)), c);
	CHECK(c.getHeight() == T(4));
	CHECK(c.getHalfHeight() == T(2));
	CHECK(c.getTotalHeight() == T(6));
	CHECK(c.getRadius() == T(1));
	CHECK_APPROX(c.getSurfaceArea(), T(12)*Constants<T>::PI);
	CHECK_APPROX(c.getVolume(), (T(4)/T(3) + T(4))*Constants<T>::PI);

	C d = c;
	d.setCenter(V::ZERO);
	d.setAxis(V::UNIT_Z);
	d.setHeight(T(2));
	d.setRadius(T(3));
	CHECK(d == C(V::ZERO, V::UNIT_Z, T(2), T(3)));
	d.setTotalHeight(T(10));
	CHECK(d.height == T(4));
	CHECK(C().set(V::ZERO, V::UNIT_Z, T(4), T(3)) == d);
	CHECK(C().set(V::ZERO, Axis::Z, T(4), T(3)) == d);
	CHECK(C(c).translate(V::ONE).center == V(T(2), T(3), T(4)));
	CHECK_APPROX(C(V::ZERO, V(T(0), T(2), T(0)), T(2), T(1)).normalize(), C(V::ZERO, V::UNIT_Y, T(4), T(1)));
	CHECK(c.isFinite());
	CHECK_FALSE(C(V::ZERO, V::UNIT_Z, T(1), std::numeric_limits<T>::infinity()).isFinite());
	CHECK(c.approxEquals(C(c.center, c.axis, T(4.05), T(1)), T(0.1)));
	CHECK_FALSE(c.approxEquals(C(c.center, c.axis, T(4.05), T(1))));
	CHECK(c != d);

	OrientedBox<T> box = c.getCircumscribedBox();
	CHECK_APPROX(box.getVolume(), T(2)*T(6)*T(2));
	CHECK(box.contains(c.center + V(T(0), T(2.99), T(0))));
	CHECK_FALSE(box.contains(c.center + V(T(0), T(3.01), T(0))));

	std::stringstream s;
	s << c;
	C read;
	s >> read;
	CHECK(read == c);
	std::unordered_set<C> set{ c, c, d };
	CHECK(set.size() == 2);
}

} // TEST_SUITE

TEST_SUITE("Cone")
{

TEST_CASE_TEMPLATE("construction and properties", T, TEST_REAL_TYPES)
{
	using C = Cone<T>;
	using V = Vector3<T>;
	C c(V::ZERO, V::UNIT_Z, T(2), T(1));
	CHECK_APPROX(C(V::ZERO, V(T(0), T(0), T(2)), T(1)), c);
	CHECK(c.getVertex() == V::ZERO);
	CHECK(c.getAxis() == V::UNIT_Z);
	CHECK(c.getHeight() == T(2));
	CHECK(c.getHalfHeight() == T(1));
	CHECK(c.getBaseCenter() == V(T(0), T(0), T(2)));
	CHECK(c.getBaseRadius() == T(1));
	CHECK_APPROX(c.getHalfAngle(), std::atan(T(0.5)));
	CHECK_APPROX(c.getAngle(), T(2)*std::atan(T(0.5)));
	CHECK_APPROX(c.getSolidAngle(), Constants<T>::TWO_PI*(T(1) - T(2)/std::sqrt(T(5))));
	CHECK_APPROX(c.getSlantHeight(), std::sqrt(T(5)));
	CHECK_APPROX(c.getSurfaceArea(), Constants<T>::PI*(T(1) + std::sqrt(T(5))));
	CHECK_APPROX(c.getVolume(), T(2)/T(3)*Constants<T>::PI);

	C d = c;
	d.setVertex(V::ONE);
	d.setAxis(V::UNIT_X);
	d.setHeight(T(3));
	d.setBaseRadius(T(2));
	CHECK(d == C(V::ONE, V::UNIT_X, T(3), T(2)));
	C e = c;
	e.setBaseCenter(V(T(0), T(4), T(0)));
	CHECK_APPROX(e, C(V::ZERO, V::UNIT_Y, T(4), T(1)));
	e.setHalfAngle(Constants<T>::PI/T(4));
	CHECK_APPROX(e.radius, T(4));
	e.setAngle(Constants<T>::HALF_PI);
	CHECK_APPROX(e.radius, T(4));
	CHECK(C().set(V::ONE, V::UNIT_X, T(3), T(2)) == d);
	CHECK(C(c).translate(V::ONE).vertex == V::ONE);
	CHECK_APPROX(C(V::ZERO, V(T(0), T(0), T(2)), T(1), T(1)).normalize(), C(V::ZERO, V::UNIT_Z, T(2), T(1)));
	CHECK(c.isFinite());
	CHECK_FALSE(C(V::ZERO, V::UNIT_Z, T(1), std::numeric_limits<T>::infinity()).isFinite());
	CHECK(c.approxEquals(C(c.vertex, c.axis, T(2.05), T(1)), T(0.1)));
	CHECK_FALSE(c.approxEquals(C(c.vertex, c.axis, T(2.05), T(1))));
	CHECK(c != d);

	std::stringstream s;
	s << c;
	C read;
	s >> read;
	CHECK(read == c);
	std::unordered_set<C> set{ c, c, d };
	CHECK(set.size() == 2);
}

TEST_CASE_TEMPLATE("bounding volumes", T, TEST_REAL_TYPES)
{
	using C = Cone<T>;
	using V = Vector3<T>;
	C c(V(T(1), T(1), T(1)), normalize(V(T(0), T(1), T(1))), T(2), T(1));
	V side = normalize(cross(c.axis, V::UNIT_X));
	std::vector<V> extremes{ c.vertex, c.getBaseCenter() + V::UNIT_X, c.getBaseCenter() - V::UNIT_X, c.getBaseCenter() + side,
		c.getBaseCenter() - side };

	OrientedBox<T> box = c.getCircumscribedBox();
	CHECK_APPROX(box.getVolume(), T(8));
	for (const V& v : extremes)
		CHECK(box.contains(lerp(box.center, v, T(0.999))));

	Sphere<T> sphere = c.getCircumscribedSphere();
	CHECK_APPROX(sphere.radius, T(1.25));
	for (const V& v : extremes)
		CHECK(containsApprox(sphere, v));

	C flat(V::ZERO, V::UNIT_Z, T(0.5), T(2));
	Sphere<T> flatSphere = flat.getCircumscribedSphere();
	CHECK_APPROX(flatSphere, Sphere<T>(V(T(0), T(0), T(0.5)), T(2)));
	CHECK(containsApprox(flatSphere, flat.vertex));
}

TEST_CASE_TEMPLATE("containment and intersection", T, TEST_REAL_TYPES)
{
	using C = Cone<T>;
	using V = Vector3<T>;
	C c(V::ZERO, V::UNIT_Z, T(2), T(1));
	CHECK(c.contains(V(T(0), T(0), T(1))));
	CHECK(c.contains(V(T(0.4), T(0), T(1))));
	CHECK_FALSE(c.contains(V(T(0.6), T(0), T(1))));
	CHECK_FALSE(c.contains(V(T(0), T(0), T(2.1))));
	CHECK_FALSE(c.contains(V(T(0), T(0), T(-0.1))));

	CHECK(c.intersects(Sphere<T>(V(T(0), T(0), T(1)), T(0.1))));
	CHECK(c.intersects(Sphere<T>(V(T(2), T(0), T(1)), T(1.5))));
	CHECK_FALSE(c.intersects(Sphere<T>(V(T(2), T(0), T(1)), T(0.5))));
	CHECK(c.intersects(Sphere<T>(V(T(0), T(0), T(-1)), T(1.1))));
	CHECK_FALSE(c.intersects(Sphere<T>(V(T(0), T(0), T(-1)), T(0.5))));
	CHECK(c.intersects(Sphere<T>(V(T(0), T(0), T(3)), T(1.2))));
	CHECK_FALSE(c.intersects(Sphere<T>(V(T(0), T(0), T(3)), T(0.5))));
}

} // TEST_SUITE

TEST_SUITE("SymmetricFrustum")
{

TEST_CASE_TEMPLATE("construction and properties", T, TEST_REAL_TYPES)
{
	using F = SymmetricFrustum<T>;
	using V = Vector3<T>;
	using V2 = Vector2<T>;
	F f = makeFrustum<T>();
	CHECK(F(V::ZERO, Matrix3<T>::IDENTITY, V2(T(1), T(0.5)), Interval<T>(T(1), T(4))) == f);
	CHECK(f.getOrigin() == V::ZERO);
	CHECK(f.getBasis().isIdentity());
	CHECK(f.getHalfDimensions() == V2(T(1), T(0.5)));
	CHECK(f.getDimensions() == V2(T(2), T(1)));
	CHECK(f.getBaseHalfDimensions() == V2(T(4), T(2)));
	CHECK(f.getBaseDimensions() == V2(T(8), T(4)));
	CHECK(f.getDepthRange() == Interval<T>(T(1), T(4)));
	CHECK(f.getMinDepth() == T(1));
	CHECK(f.getMaxDepth() == T(4));
	CHECK(f.getDepthRatio() == T(4));

	F g = f;
	g.setOrigin(V::ONE);
	g.setBasis(Matrix3<T>::makeRotation(Axis::Y, T(0.3)));
	g.setHalfDimensions(V2(T(2)));
	g.setDepthRange(Interval<T>(T(2), T(3)));
	CHECK(g == F(V::ONE, Matrix3<T>::makeRotation(Axis::Y, T(0.3)), V2(T(2)), T(2), T(3)));
	g.setBaseHalfDimensions(V2(T(6)));
	CHECK_APPROX(g.halfDims, V2(T(4)));
	CHECK(F().set(V::ZERO, Matrix3<T>::IDENTITY, V2(T(1), T(0.5)), Interval<T>(T(1), T(4))) == f);
	CHECK(F(f).translate(V::ONE).origin == V::ONE);
	CHECK(f.isFinite());
	CHECK(f.approxEquals(F(V::ZERO, Matrix3<T>::IDENTITY, V2(T(1.05), T(0.5)), T(1), T(4)), T(0.1)));
	CHECK_FALSE(f.approxEquals(F(V::ZERO, Matrix3<T>::IDENTITY, V2(T(1.05), T(0.5)), T(1), T(4))));
	CHECK(f != g);

	F n(V::ZERO, Matrix3<T>::makeScaling(V(T(2), T(3), T(4))), V2(T(1)), T(1), T(2));
	n.orthonormalize();
	CHECK_APPROX(n, F(V::ZERO, Matrix3<T>::IDENTITY, V2(T(2), T(3)), T(4), T(8)));

	std::stringstream s;
	s << f;
	F read;
	s >> read;
	CHECK(read == f);
	std::unordered_set<F> set{ f, f, g };
	CHECK(set.size() == 2);
}

TEST_CASE_TEMPLATE("vertices, primitives and half spaces", T, TEST_REAL_TYPES)
{
	using F = SymmetricFrustum<T>;
	using V = Vector3<T>;
	F f(V(T(1), T(2), T(3)), Matrix3<T>::makeRotation(V(T(1), T(2), T(0)), T(0.6)), Vector2<T>(T(1), T(0.5)), T(1), T(4));
	std::vector<V> vertices = f.getVertices();
	REQUIRE(vertices.size() == 8);
	for (std::size_t i = 0; i < 8; i++)
	{
		T depth = (i & 4) ? T(4) : T(1);
		V local((i & 1) ? depth : -depth, ((i & 2) ? depth : -depth)*T(0.5), depth);
		CHECK_APPROX(vertices[i], local*f.basis + f.origin);
	}
	std::vector<V> copied;
	f.copyVertices(std::back_inserter(copied));
	REQUIRE(copied.size() == 8);
	for (std::size_t i = 0; i < 8; i++)
		CHECK_APPROX(copied[i], vertices[i]);

	V center = f.origin + f.basis[2]*T(2.5);
	auto triangles = f.template getPrimitives<std::uint32_t>(3);
	REQUIRE(triangles.second - triangles.first == 36);
	for (const std::uint32_t* i = triangles.first; i != triangles.second; i += 3)
	{
		Triangle3<T> t(vertices[i[0]], vertices[i[1]], vertices[i[2]]);
		CHECK(dot(t.getNormal(), (t.vertices[0] + t.vertices[1] + t.vertices[2])/T(3) - center) > T(0));
	}
	CHECK(f.template getPrimitives<std::uint32_t>(2).second - f.template getPrimitives<std::uint32_t>(2).first == 24);
	CHECK(f.template getPrimitives<std::uint32_t>(4).second - f.template getPrimitives<std::uint32_t>(4).first == 24);
	CHECK(f.template getPrimitives<std::uint32_t>(5).first == nullptr);
	CHECK(f.getPrimitiveCount(2) == 12);
	CHECK(f.getPrimitiveCount(3) == 12);
	CHECK(f.getPrimitiveCount(4) == 6);
	CHECK(f.getPrimitiveCount(0) == 0);

	std::vector<HalfSpace<T>> halfSpaces = f.getHalfSpaces();
	REQUIRE(halfSpaces.size() == 6);
	std::vector<HalfSpace<T>> copiedHalfSpaces;
	f.copyHalfSpaces(std::back_inserter(copiedHalfSpaces));
	REQUIRE(copiedHalfSpaces.size() == 6);
	int enumerated = 0;
	CHECK(f.enumerateHalfSpaces([&](const HalfSpace<T>& h) { CHECK_APPROX(h, halfSpaces[enumerated++]); return true; }));
	CHECK(enumerated == 6);
	CHECK_FALSE(f.enumerateHalfSpaces([](const HalfSpace<T>&) { return false; }));
	for (std::size_t i = 0; i < 6; i++)
	{
		CAPTURE(i);
		CHECK_APPROX(copiedHalfSpaces[i], halfSpaces[i]);
		CHECK(halfSpaces[i].contains(center)); // normals point outwards
		for (const V& v : vertices)
			CHECK(halfSpaces[i].getSignedDistanceTo(v) <= T(10)*test::tolerance<T>());
	}

	// A left-handed basis must not flip the half spaces inside out.
	F mirrored(V::ZERO, Matrix3<T>::makeScaling(V(T(-1), T(1), T(1))), Vector2<T>(T(1), T(0.5)), T(1), T(4));
	for (const HalfSpace<T>& h : mirrored.getHalfSpaces())
		CHECK(h.contains(V(T(0), T(0), T(2.5))));

	// An infinite far plane omits the last half space.
	F infinite(V::ZERO, Matrix3<T>::IDENTITY, Vector2<T>(T(1)), T(1), std::numeric_limits<T>::max());
	CHECK(infinite.getHalfSpaces().size() == 5);
	std::vector<HalfSpace<T>> infiniteHalfSpaces;
	infinite.copyHalfSpaces(std::back_inserter(infiniteHalfSpaces));
	CHECK(infiniteHalfSpaces.size() == 5);
}

TEST_CASE_TEMPLATE("bounding volumes", T, TEST_REAL_TYPES)
{
	using F = SymmetricFrustum<T>;
	using V = Vector3<T>;
	for (F f : { makeFrustum<T>(), F(V(T(1), T(2), T(3)), Matrix3<T>::makeRotation(Axis::X, T(0.4)), Vector2<T>(T(2), T(3)), T(1), T(1.5)) })
	{
		std::vector<V> vertices = f.getVertices();
		V center = f.origin + f.basis[2]*f.depthRange.getCenter();

		OrientedBox<T> box = f.getCircumscribedBox();
		for (const V& v : vertices)
			CHECK(box.contains(lerp(center, v, T(0.999))));

		Sphere<T> sphere = f.getCircumscribedSphere();
		for (const V& v : vertices)
			CHECK(containsApprox(sphere, v));

		Cone<T> cone = f.getCircumscribedCone();
		for (const V& v : vertices)
			CHECK(cone.contains(lerp(center, v, T(0.999))));
	}
}

TEST_CASE_TEMPLATE("closest point, distance and containment", T, TEST_REAL_TYPES)
{
	using V = Vector3<T>;
	auto f = makeFrustum<T>();
	CHECK(f.contains(V(T(0), T(0), T(2))));
	CHECK(f.contains(V(T(3.9), T(1.9), T(3.99))));
	CHECK_FALSE(f.contains(V(T(0), T(0), T(0.5))));
	CHECK_FALSE(f.contains(V(T(0), T(0), T(5))));
	CHECK_FALSE(f.contains(V(T(2.5), T(0), T(2))));

	CHECK_APPROX(f.getClosestPoint(V(T(0), T(0), T(2))), V(T(0), T(0), T(2)));
	CHECK(f.getDistanceTo(V(T(0), T(0), T(2))) == T(0));
	CHECK_APPROX(f.getClosestPoint(V(T(0), T(0), T(0))), V(T(0), T(0), T(1)));
	CHECK_APPROX(f.getDistanceTo(V(T(0), T(0), T(0))), T(1));
	CHECK_APPROX(f.getClosestPoint(V(T(0), T(0), T(6))), V(T(0), T(0), T(4)));
	CHECK_APPROX(f.getClosestPoint(V(T(3), T(0), T(2))), V(T(2.5), T(0), T(2.5)));
	CHECK_APPROX(f.getDistanceTo(V(T(3), T(0), T(2))), std::sqrt(T(0.5)));
	CHECK_APPROX(f.getClosestPoint(V(T(10), T(10), T(10))), V(T(4), T(2), T(4))); // far corner
}

TEST_CASE_TEMPLATE("intersections", T, TEST_REAL_TYPES)
{
	using V = Vector3<T>;
	auto f = makeFrustum<T>();
	Matrix3<T> r = Matrix3<T>::makeRotation(V(T(1), T(1), T(1)), T(0.5));

	CHECK(f.intersects(Sphere<T>(V(T(3), T(0), T(2)), T(0.8))));
	CHECK_FALSE(f.intersects(Sphere<T>(V(T(3), T(0), T(2)), T(0.6))));
	CHECK(f.intersects(AxisAlignedBox<T>(V(T(-0.2), T(-0.2), T(1.8)), V(T(0.2), T(0.2), T(2.2)))));
	CHECK(f.intersects(AxisAlignedBox<T>(V(T(-10)), V(T(10))))); // containing box
	CHECK_FALSE(f.intersects(AxisAlignedBox<T>(V(T(2.8), T(-0.2), T(1.8)), V(T(3.2), T(0.2), T(2.2)))));
	CHECK_FALSE(f.intersects(AxisAlignedBox<T>(V(T(-1), T(-1), T(-2)), V(T(1), T(1), T(0.5)))));
	CHECK(f.intersects(OrientedBox<T>(V(T(0), T(0), T(2)), r, V(T(0.2)))));
	CHECK(f.intersects(OrientedBox<T>(V(T(2.4), T(0), T(2)), r, V(T(0.2)))));
	CHECK_FALSE(f.intersects(OrientedBox<T>(V(T(3.2), T(0), T(2)), r, V(T(0.2)))));
	CHECK_FALSE(f.intersects(OrientedBox<T>(V(T(0), T(0), T(5)), r, V(T(0.2)))));
}

} // TEST_SUITE

} // namespace mathematics::templates
