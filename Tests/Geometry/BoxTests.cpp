#include "TestUtilities.hpp"
#include <unordered_set>
#include <vector>
#include <iterator>
#include <cstdint>

namespace mathematics::templates {

namespace {

// Checks that the primitive index tables describe the surface of a box with the given vertices.
template<typename T>
void checkPrimitives(const std::vector<Vector3<T>>& vertices, const Vector3<T>& center, std::pair<const std::uint32_t*, const std::uint32_t*> triangles,
	std::pair<const std::uint32_t*, const std::uint32_t*> quads, std::pair<const std::uint32_t*, const std::uint32_t*> edges, T totalEdgeLength)
{
	REQUIRE(triangles.second - triangles.first == 36);
	for (const std::uint32_t* i = triangles.first; i != triangles.second; i += 3)
	{
		Triangle3<T> t(vertices[i[0]], vertices[i[1]], vertices[i[2]]);
		CHECK(dot(t.getNormal(), (t.vertices[0] + t.vertices[1] + t.vertices[2])/T(3) - center) > T(0)); // counter-clockwise, facing outwards
	}

	REQUIRE(quads.second - quads.first == 24);
	for (const std::uint32_t* i = quads.first; i != quads.second; i += 4)
	{
		Quadrilateral3<T> q(vertices[i[0]], vertices[i[1]], vertices[i[2]], vertices[i[3]]);
		CHECK(dot(q.getNormal(), (q.vertices[0] + q.vertices[2])*T(0.5) - center) > T(0));
	}

	REQUIRE(edges.second - edges.first == 24);
	T length = T(0);
	for (const std::uint32_t* i = edges.first; i != edges.second; i += 2)
		length += distance(vertices[i[0]], vertices[i[1]]);
	CHECK_APPROX(length, totalEdgeLength);
}

} // namespace

TEST_SUITE("AxisAlignedBox")
{

TEST_CASE_TEMPLATE("construction and properties", T, TEST_REAL_TYPES)
{
	using B = AxisAlignedBox<T>;
	using V = Vector3<T>;
	CHECK(B().isEmpty());
	CHECK(B() == B::EMPTY);
	B b(V(T(1), T(2), T(3)), V(T(3), T(6), T(4)));
	CHECK_FALSE(b.isEmpty());
	CHECK(B(V(T(2), T(4), T(6))) == B(V(T(-1), T(-2), T(-3)), V(T(1), T(2), T(3))));
	CHECK(B(std::pair<V, V>(b.minimum, b.maximum)) == b);
	CHECK(B(std::tuple<V, V>(b.minimum, b.maximum)) == b);

	CHECK(b.getMinimum() == V(T(1), T(2), T(3)));
	CHECK(b.getMaximum() == V(T(3), T(6), T(4)));
	CHECK(b.getDimensions() == V(T(2), T(4), T(1)));
	CHECK(b.getHalfDimensions() == V(T(1), T(2), T(0.5)));
	CHECK(b.getCenter() == V(T(2), T(4), T(3.5)));
	CHECK_APPROX(b.getDiagonal(), std::sqrt(T(21)));
	CHECK(b.getSurfaceArea() == T(28));
	CHECK(b.getVolume() == T(8));

	B c = b;
	c.setDimensions(V(T(2)));
	CHECK(c == B(V(T(1), T(3), T(2.5)), V(T(3), T(5), T(4.5))));
	c.setHalfDimensions(V(T(1)));
	CHECK(c == B(V(T(1), T(3), T(2.5)), V(T(3), T(5), T(4.5))));
	c.setCenter(V::ZERO);
	CHECK(c == B(V(T(-1)), V(T(1))));
	c.setMinimum(V(T(-2)));
	c.setMaximum(V(T(2)));
	CHECK(c == B(V(T(-2)), V(T(2))));
	CHECK(B(b).set(V::ZERO, V::ONE) == B(V::ZERO, V::ONE));
	CHECK(B(b).clear().isEmpty());

	CHECK(b.isFinite());
	CHECK_FALSE(B::EMPTY.isFinite());
	CHECK(b.approxEquals(B(b.minimum, b.maximum + V(T(0.05))), T(0.1)));
	CHECK_FALSE(b.approxEquals(B(b.minimum, b.maximum + V(T(0.05)))));
	CHECK(b != c);
}

TEST_CASE_TEMPLATE("vertices, primitives and half spaces", T, TEST_REAL_TYPES)
{
	using B = AxisAlignedBox<T>;
	using V = Vector3<T>;
	B b(V(T(1), T(2), T(3)), V(T(3), T(6), T(4)));
	std::vector<V> vertices = b.getVertices();
	REQUIRE(vertices.size() == 8);
	for (std::size_t i = 0; i < 8; i++) // bit i of the index selects the maximum along axis i
		CHECK(vertices[i] == V((i & 1) ? b.maximum.x : b.minimum.x, (i & 2) ? b.maximum.y : b.minimum.y, (i & 4) ? b.maximum.z : b.minimum.z));
	std::vector<V> copied(8);
	CHECK(b.copyVertices(copied.begin()) == copied.end());
	CHECK(copied == vertices);

	checkPrimitives(vertices, b.getCenter(), b.template getPrimitives<std::uint32_t>(3), b.template getPrimitives<std::uint32_t>(4),
		b.template getPrimitives<std::uint32_t>(2), T(4)*(T(2) + T(4) + T(1)));
	CHECK(b.template getPrimitives<int>(5).first == nullptr);
	CHECK(b.getPrimitiveCount(2) == 12);
	CHECK(b.getPrimitiveCount(3) == 12);
	CHECK(b.getPrimitiveCount(4) == 6);
	CHECK(b.getPrimitiveCount(1) == 0);

	std::vector<HalfSpace<T>> halfSpaces = b.getHalfSpaces();
	REQUIRE(halfSpaces.size() == 6);
	std::vector<HalfSpace<T>> copiedHalfSpaces;
	b.copyHalfSpaces(std::back_inserter(copiedHalfSpaces));
	CHECK(copiedHalfSpaces == halfSpaces);
	for (const HalfSpace<T>& h : halfSpaces)
	{
		CHECK(h.contains(b.getCenter()));
		CHECK_FALSE(h.contains(b.getCenter() + h.getNormal()*T(10)));
		CHECK_APPROX(h.getSignedDistanceTo(b.getCenter()), -std::fabs(dot(h.getNormal(), b.getHalfDimensions())));
	}

	Sphere<T> s = b.getCircumscribedSphere();
	CHECK(s.center == b.getCenter());
	for (const V& v : vertices)
		CHECK_APPROX(distance(v, s.center), s.radius);
}

TEST_CASE_TEMPLATE("transformation, union and intersection", T, TEST_REAL_TYPES)
{
	using B = AxisAlignedBox<T>;
	using V = Vector3<T>;
	B b(V(T(0)), V(T(2)));
	CHECK(B(b).inflate(V(T(1), T(2), T(3))) == B(V(T(-1), T(-2), T(-3)), V(T(3), T(4), T(5))));
	CHECK(B(b).translate(V(T(1), T(-1), T(0))) == B(V(T(1), T(-1), T(0)), V(T(3), T(1), T(2))));
	CHECK(B(b).scaleAroundCenter(T(2)) == B(V(T(-1)), V(T(3))));

	// The transformed box encloses all transformed vertices
	Matrix3<T> m = Matrix3<T>::makeScalingRotation(V(T(1), T(2), T(1)), Axis::Z, Constants<T>::PI/T(4));
	B t = B(b).transform(m);
	for (const V& v : b.getVertices())
		CHECK(t.approxEquals(B::makeUnion(t, B(v*m, v*m)), test::tolerance<T>()));
	CHECK_APPROX(B(b).transform(Matrix3<T>::makeRotation(Axis::Z, Constants<T>::HALF_PI), true), B(V(T(-2), T(0), T(0)), V(T(0), T(2), T(2))));
	AffineTransform<T> a = AffineTransform<T>::makeRotationTranslation(Axis::Z, Constants<T>::HALF_PI, V(T(10), T(0), T(0)));
	CHECK_APPROX(B(b).transform(a), B(V(T(8), T(0), T(0)), V(T(10), T(2), T(2))));
	CHECK_APPROX(B(b).transform(a, true), B(V(T(8), T(0), T(0)), V(T(10), T(2), T(2))));

	B c(V(T(1)), V(T(3)));
	CHECK(B::makeUnion(b, c) == B(V(T(0)), V(T(3))));
	CHECK(B::makeIntersection(b, c) == B(V(T(1)), V(T(2))));
	CHECK(B::makeIntersection(b, B(V(T(5)), V(T(6)))).isEmpty());

	B e;
	e.extendBy(V(T(1), T(1), T(1)));
	e.extendBy(V(T(-1), T(3), T(0)));
	CHECK(e == B(V(T(-1), T(1), T(0)), V(T(1), T(3), T(1))));
	e.extendBy(B(V(T(0)), V(T(2))));
	CHECK(e == B(V(T(-1), T(0), T(0)), V(T(2), T(3), T(2))));
	e.extendBy(Sphere<T>(V(T(5)), T(1)));
	CHECK(e == B(V(T(-1), T(0), T(0)), V(T(6))));

	OrientedBox<T> o(V(T(1), T(1), T(1)), Matrix3<T>::makeRotation(Axis::Z, Constants<T>::PI/T(4)), V(T(1), T(1), T(1)));
	CHECK_APPROX(B(o), B(V(T(1) - std::sqrt(T(2)), T(1) - std::sqrt(T(2)), T(0)), V(T(1) + std::sqrt(T(2)), T(1) + std::sqrt(T(2)), T(2))));
}

TEST_CASE_TEMPLATE("closest point, containment and classification", T, TEST_REAL_TYPES)
{
	using B = AxisAlignedBox<T>;
	using V = Vector3<T>;
	B b(V(T(0)), V(T(2)));
	CHECK(b.getClosestPoint(V(T(1))) == V(T(1)));
	CHECK(b.getClosestPoint(V(T(5), T(-1), T(1))) == V(T(2), T(0), T(1)));
	CHECK_APPROX(b.getDistanceTo(V(T(5), T(6), T(1))), T(5));

	CHECK(b.contains(V(T(1))));
	CHECK(b.contains(V(T(2), T(0), T(1))));
	CHECK_FALSE(b.contains(V(T(2.1), T(0), T(1))));
	CHECK(b.contains(B(V(T(0.5)), V(T(1.5)))));
	CHECK_FALSE(b.contains(B(V(T(0.5)), V(T(2.5)))));
	CHECK(b.contains(Sphere<T>(V(T(1)), T(1))));
	CHECK_FALSE(b.contains(Sphere<T>(V(T(1)), T(1.5))));

	CHECK(b.classify(HalfSpace<T>(V::UNIT_X, T(-5))) == -1); // x <= 5: inside
	CHECK(b.classify(HalfSpace<T>(V::UNIT_X, T(5))) == 1);	// x <= -5: outside
	CHECK(b.classify(HalfSpace<T>(V::UNIT_X, T(-1))) == 0);	// x <= 1: partial
}

TEST_CASE_TEMPLATE("intersections", T, TEST_REAL_TYPES)
{
	using B = AxisAlignedBox<T>;
	using V = Vector3<T>;
	B b(V(T(0)), V(T(2)));
	Matrix3<T> r = Matrix3<T>::makeRotation(V(T(1), T(1), T(0)), T(0.7));

	CHECK(b.intersects(HalfSpace<T>(V::UNIT_X, T(-1))));
	CHECK(b.intersects(HalfSpace<T>(V::UNIT_X, T(-5))));
	CHECK_FALSE(b.intersects(HalfSpace<T>(V::UNIT_X, T(5))));
	CHECK(b.intersects(Plane<T>(V::UNIT_X, T(-1))));
	CHECK_FALSE(b.intersects(Plane<T>(V::UNIT_X, T(-5))));
	CHECK(b.intersects(Triangle3<T>(V(T(1), T(1), T(-5)), V(T(1), T(1), T(5)), V(T(1), T(5), T(0)))));
	CHECK_FALSE(b.intersects(Triangle3<T>(V(T(3), T(3), T(-5)), V(T(3), T(3), T(5)), V(T(3), T(5), T(0)))));
	CHECK(b.intersects(B(V(T(1)), V(T(3)))));
	CHECK(b.intersects(B(V(T(2)), V(T(3)))));
	CHECK_FALSE(b.intersects(B(V(T(2.5)), V(T(3)))));
	CHECK(b.intersects(OrientedBox<T>(V(T(2.5), T(2.5), T(1)), r, V(T(1)))));
	CHECK_FALSE(b.intersects(OrientedBox<T>(V(T(5)), r, V(T(1)))));
	CHECK(b.intersects(Sphere<T>(V(T(3), T(1), T(1)), T(1.5))));
	CHECK_FALSE(b.intersects(Sphere<T>(V(T(3)), T(1))));
}

TEST_CASE_TEMPLATE("streaming and hashing", T, TEST_REAL_TYPES)
{
	using B = AxisAlignedBox<T>;
	using V = Vector3<T>;
	B b(V(T(1), T(2), T(3)), V(T(4), T(5), T(6.5)));
	std::stringstream s;
	s << b;
	CHECK(s.str() == "1 2 3 4 5 6.5");
	B read;
	s >> read;
	CHECK(read == b);
	std::unordered_set<B> set{ b, b, B() };
	CHECK(set.size() == 2);
}

} // TEST_SUITE

TEST_SUITE("OrientedBox")
{

TEST_CASE_TEMPLATE("construction and properties", T, TEST_REAL_TYPES)
{
	using O = OrientedBox<T>;
	using V = Vector3<T>;
	using M = Matrix3<T>;
	CHECK(O().basis.isIdentity());
	M r = M::makeRotation(Axis::Z, Constants<T>::HALF_PI);
	O o(V(T(1), T(2), T(3)), r, V(T(1), T(2), T(3)));
	CHECK(o.getCenter() == V(T(1), T(2), T(3)));
	CHECK(o.getBasis() == r);
	CHECK(o.getHalfDimensions() == V(T(1), T(2), T(3)));
	CHECK(o.getDimensions() == V(T(2), T(4), T(6)));
	CHECK_APPROX(o.getDiagonal(), T(2)*std::sqrt(T(14)));
	CHECK(o.getSurfaceArea() == T(88));
	CHECK(o.getVolume() == T(48));

	AxisAlignedBox<T> b(V(T(0)), V(T(2), T(4), T(6)));
	O fromBox(b);
	CHECK(fromBox == O(V(T(1), T(2), T(3)), M::IDENTITY, V(T(1), T(2), T(3))));
	O rotated(b, r);
	CHECK_APPROX(rotated, O(V(T(1), T(2), T(3))*r, r, V(T(1), T(2), T(3))));
	O scaled(b, M::makeScaling(V(T(2), T(1), T(1))));
	CHECK_APPROX(scaled, O(V(T(2), T(2), T(3)), M::IDENTITY, V(T(2), T(2), T(3)))); // scale moves into the half dimensions
	O affine(b, AffineTransform<T>(r, V(T(10), T(0), T(0))), true);
	CHECK_APPROX(affine, O(V(T(1), T(2), T(3))*r + V(T(10), T(0), T(0)), r, V(T(1), T(2), T(3))));

	O p = o;
	p.setCenter(V::ZERO);
	p.setBasis(M::IDENTITY);
	p.setHalfDimensions(V::ONE);
	CHECK(p == O(V::ZERO, M::IDENTITY, V::ONE));
	CHECK(O().set(V::ZERO, M::IDENTITY, V::ONE) == p);
	CHECK(o.isFinite());
	CHECK_FALSE(O(V::INF, M::IDENTITY, V::ONE).isFinite());
	CHECK(o.approxEquals(O(o.center, o.basis, o.halfDims + V(T(0.05))), T(0.1)));
	CHECK_FALSE(o.approxEquals(O(o.center, o.basis, o.halfDims + V(T(0.05)))));
	CHECK(o != p);
}

TEST_CASE_TEMPLATE("vertices, primitives and half spaces", T, TEST_REAL_TYPES)
{
	using O = OrientedBox<T>;
	using V = Vector3<T>;
	O o(V(T(1), T(2), T(3)), Matrix3<T>::makeRotation(V(T(1), T(2), T(3)), T(0.8)), V(T(1), T(2), T(0.5)));
	std::vector<V> vertices = o.getVertices();
	REQUIRE(vertices.size() == 8);
	for (std::size_t i = 0; i < 8; i++)
	{
		V local((i & 1) ? o.halfDims.x : -o.halfDims.x, (i & 2) ? o.halfDims.y : -o.halfDims.y, (i & 4) ? o.halfDims.z : -o.halfDims.z);
		CHECK_APPROX(vertices[i], local*o.basis + o.center);
	}
	std::vector<V> copied(8);
	CHECK(o.copyVertices(copied.begin()) == copied.end());
	for (std::size_t i = 0; i < 8; i++)
		CHECK_APPROX(copied[i], vertices[i]);

	checkPrimitives(vertices, o.center, o.template getPrimitives<std::uint32_t>(3), o.template getPrimitives<std::uint32_t>(4),
		o.template getPrimitives<std::uint32_t>(2), T(8)*(T(1) + T(2) + T(0.5)));
	CHECK(o.template getPrimitives<int>(1).first == nullptr);
	CHECK(o.getPrimitiveCount(2) == 12);
	CHECK(o.getPrimitiveCount(3) == 12);
	CHECK(o.getPrimitiveCount(4) == 6);
	CHECK(o.getPrimitiveCount(5) == 0);

	std::vector<HalfSpace<T>> halfSpaces = o.getHalfSpaces();
	REQUIRE(halfSpaces.size() == 6);
	std::vector<HalfSpace<T>> copiedHalfSpaces;
	o.copyHalfSpaces(std::back_inserter(copiedHalfSpaces));
	CHECK(copiedHalfSpaces.size() == 6);
	for (const HalfSpace<T>& h : halfSpaces)
	{
		CHECK(h.contains(o.center));
		for (const V& v : vertices)
			CHECK(h.getSignedDistanceTo(v) <= T(10)*test::tolerance<T>());
	}

	AxisAlignedBox<T> box = o.getCircumscribedBox();
	for (const V& v : vertices)
		CHECK(box.approxEquals(AxisAlignedBox<T>::makeUnion(box, AxisAlignedBox<T>(v, v)), test::tolerance<T>()));
	CHECK_APPROX(box, AxisAlignedBox<T>(o));

	Sphere<T> sphere = o.getCircumscribedSphere();
	for (const V& v : vertices)
		CHECK_APPROX(distance(v, sphere.center), sphere.radius);
}

TEST_CASE_TEMPLATE("transformation", T, TEST_REAL_TYPES)
{
	using O = OrientedBox<T>;
	using V = Vector3<T>;
	using M = Matrix3<T>;
	O o(V(T(1), T(2), T(3)), M::IDENTITY, V(T(1), T(2), T(3)));
	CHECK(O(o).translate(V::ONE).center == V(T(2), T(3), T(4)));

	M m = M::makeScalingRotation(V(T(2), T(1), T(1)), Axis::Z, T(0.5));
	O t = O(o).transform(m);
	CHECK(t.basis.isApproxOrthogonal());
	for (const V& v : o.getVertices())
	{
		std::vector<V> transformed = t.getVertices();
		bool found = false;
		for (const V& u : transformed)
			found = found || u.approxEquals(v*m, T(10)*test::tolerance<T>());
		CHECK(found);
	}

	M r = M::makeRotation(Axis::X, T(0.5));
	CHECK_APPROX(O(o).transform(r, true), O(o.center*r, r, o.halfDims));
	AffineTransform<T> a(r, V(T(1), T(0), T(0)));
	CHECK_APPROX(O(o).transform(a), O(transform(o.center, a), r, o.halfDims));
	CHECK_APPROX(O(o).transform(a, true), O(transform(o.center, a), r, o.halfDims));

	O n(V::ZERO, M::makeScaling(V(T(2), T(3), T(4))), V::ONE);
	n.orthonormalize();
	CHECK_APPROX(n, O(V::ZERO, M::IDENTITY, V(T(2), T(3), T(4))));
}

TEST_CASE_TEMPLATE("closest point, containment and classification", T, TEST_REAL_TYPES)
{
	using O = OrientedBox<T>;
	using V = Vector3<T>;
	O o(V(T(1), T(1), T(0)), Matrix3<T>::makeRotation(Axis::Z, Constants<T>::PI/T(4)), V(T(1), T(1), T(1)));
	T s = std::sqrt(T(2));

	CHECK(o.contains(V(T(1), T(1), T(0))));
	CHECK(o.contains(V(T(1) + s*T(0.99), T(1), T(0)))); // inside a corner
	CHECK_FALSE(o.contains(V(T(2), T(2), T(0)))); // outside along the diagonal
	CHECK_APPROX(o.getClosestPoint(V(T(1), T(1), T(0.5))), V(T(1), T(1), T(0.5)));
	CHECK_APPROX(o.getClosestPoint(V(T(3), T(3), T(0))), V(T(1) + s/T(2), T(1) + s/T(2), T(0)));
	CHECK_APPROX(o.getDistanceTo(V(T(3), T(3), T(0))), T(2)*s - T(1));

	CHECK(o.classify(HalfSpace<T>(V::UNIT_X, T(-5))) == -1);
	CHECK(o.classify(HalfSpace<T>(V::UNIT_X, T(5))) == 1);
	CHECK(o.classify(HalfSpace<T>(V::UNIT_X, T(-1))) == 0);
}

TEST_CASE_TEMPLATE("intersections", T, TEST_REAL_TYPES)
{
	using O = OrientedBox<T>;
	using V = Vector3<T>;
	O o(V(T(0)), Matrix3<T>::makeRotation(Axis::Z, Constants<T>::PI/T(4)), V(T(1)));
	T s = std::sqrt(T(2));

	CHECK(o.intersects(HalfSpace<T>(V::UNIT_X, T(-1.3))));		// x <= 1.3
	CHECK(o.intersects(HalfSpace<T>(V::UNIT_X, s - T(0.1))));	// x <= -(sqrt(2) - 0.1)
	CHECK_FALSE(o.intersects(HalfSpace<T>(V::UNIT_X, s + T(0.1))));
	CHECK(o.intersects(Plane<T>(V::UNIT_X, T(-1.3))));
	CHECK_FALSE(o.intersects(Plane<T>(V::UNIT_X, T(-1.5))));
	CHECK(o.intersects(Triangle3<T>(V(T(1.3), T(0), T(-5)), V(T(1.3), T(0), T(5)), V(T(5), T(0), T(0)))));
	CHECK_FALSE(o.intersects(Triangle3<T>(V(T(1.5), T(0), T(-5)), V(T(1.5), T(0), T(5)), V(T(5), T(0), T(0)))));
	CHECK(o.intersects(AxisAlignedBox<T>(V(T(1.3), T(-0.1), T(-0.1)), V(T(2), T(0.1), T(0.1)))));
	CHECK_FALSE(o.intersects(AxisAlignedBox<T>(V(T(1), T(1), T(-1)), V(T(2), T(2), T(1))))); // the corner region is outside
	CHECK(o.intersects(O(V(T(2.5), T(0), T(0)), Matrix3<T>::makeRotation(Axis::Z, Constants<T>::PI/T(4)), V(T(1)))));
	CHECK_FALSE(o.intersects(O(V(T(3), T(0), T(0)), Matrix3<T>::makeRotation(Axis::Z, Constants<T>::PI/T(4)), V(T(1)))));
	CHECK(o.intersects(Sphere<T>(V(T(1.8), T(0), T(0)), T(0.5))));
	CHECK_FALSE(o.intersects(Sphere<T>(V(T(1.2), T(1.2), T(0)), T(0.3))));
}

TEST_CASE_TEMPLATE("streaming and hashing", T, TEST_REAL_TYPES)
{
	using O = OrientedBox<T>;
	using V = Vector3<T>;
	O o(V(T(1), T(2), T(3)), Matrix3<T>::IDENTITY, V(T(4), T(5), T(6)));
	std::stringstream s;
	s << o;
	O read;
	s >> read;
	CHECK(read == o);
	std::unordered_set<O> set{ o, o, O() };
	CHECK(set.size() == 2);
}

} // TEST_SUITE

} // namespace mathematics::templates
