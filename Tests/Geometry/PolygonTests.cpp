#include "TestUtilities.hpp"
#include <unordered_set>
#include <vector>
#include <iterator>
#include <cstdint>

namespace mathematics::templates {

namespace {

// Sum of triangle areas produced by an index list
template<typename T, typename I>
T triangulatedArea(const std::vector<Vector3<T>>& vertices, const std::vector<I>& indices)
{
	T area = T(0);
	for (std::size_t i = 0; i + 2 < indices.size(); i += 3)
		area += Triangle3<T>(vertices[indices[i]], vertices[indices[i + 1]], vertices[indices[i + 2]]).getArea();
	return area;
}

} // namespace

TEST_SUITE("Triangle3")
{

TEST_CASE_TEMPLATE("construction and properties", T, TEST_REAL_TYPES)
{
	using Tr = Triangle3<T>;
	using V = Vector3<T>;
	V a(T(0), T(0), T(0)), b(T(4), T(0), T(0)), c(T(0), T(3), T(0));
	Tr t(a, b, c);
	CHECK(t.vertices[1] == b);
	CHECK(Tr(std::tuple<V, V, V>(a, b, c)) == t);
	const V v[] = { a, b, c };
	CHECK(Tr(v) == t);
	CHECK(t.getVertices() == std::tuple<V, V, V>(a, b, c));
	std::vector<V> copied;
	t.copyVertices(std::back_inserter(copied));
	CHECK(copied == std::vector<V>{ a, b, c });

	CHECK(t.getVertex(2) == c);
	CHECK(t.getVertex(3) == V::ZERO); // out of range
	Tr u = t;
	u.setVertex(2, V(T(0), T(0), T(3)));
	CHECK(u.vertices[2] == V(T(0), T(0), T(3)));
	CHECK_THROWS_AS(u.setVertex(3, V::ZERO), std::out_of_range);

	CHECK(Tr(t).setZero().isZero());
	CHECK(Tr().set(a, b, c) == t);
	CHECK(Tr(V(Constants<T>::TOLERANCE*T(0.5)), V::ZERO, V::ZERO).isApproxZero());
	CHECK(t.isFinite());
	CHECK_FALSE(Tr(a, b, V::INF).isFinite());
	CHECK(t.approxEquals(Tr(a, b, c + V(T(0.05))), T(0.1)));
	CHECK_FALSE(t.approxEquals(Tr(a, b, c + V(T(0.05)))));
	CHECK(t != u);

	CHECK_APPROX(t.getPerimeter(), T(12));
	CHECK_APPROX(t.getArea(), T(6));
	CHECK_APPROX(Tr(a + V(T(5)), b + V(T(5)), c + V(T(5))).getArea(), T(6)); // translation invariant
	CHECK_APPROX(t.getNormal(), V::UNIT_Z);
	CHECK_APPROX(Tr::computeNormal(a, c, b), V::MINUS_UNIT_Z);
	CHECK(Tr(t).flip() == Tr(c, b, a));
}

TEST_CASE_TEMPLATE("tangent space", T, TEST_REAL_TYPES)
{
	using Tr = Triangle3<T>;
	using V = Vector3<T>;
	using V2 = Vector2<T>;
	Tr t(V(T(0), T(0), T(0)), V(T(2), T(0), T(0)), V(T(0), T(2), T(0)));
	Matrix3<T> basis = t.getTangentSpaceBasis(V2(T(0), T(0)), V2(T(1), T(0)), V2(T(0), T(1)), false);
	CHECK_APPROX(basis.getRow(0), V::UNIT_X); // tangent follows u
	CHECK_APPROX(basis.getRow(1), V::UNIT_Y); // bitangent follows v
	CHECK_APPROX(basis.getRow(2), V::UNIT_Z);

	Matrix3<T> weighted = Tr::computeTangentSpaceBasis(t.vertices[0], t.vertices[1], t.vertices[2],
		V2(T(0), T(0)), V2(T(1), T(0)), V2(T(0), T(1)), true);
	CHECK_APPROX(weighted.getRow(0), V(T(0.5), T(0), T(0)));

	Matrix3<T> degenerate = t.getTangentSpaceBasis(V2(T(0)), V2(T(0)), V2(T(0)), false);
	CHECK(degenerate.getRow(0).isZero());
	CHECK(degenerate.getRow(1).isZero());
}

TEST_CASE_TEMPLATE("inscribed and circumscribed shapes", T, TEST_REAL_TYPES)
{
	using Tr = Triangle3<T>;
	using V = Vector3<T>;
	Tr t(V(T(0), T(0), T(0)), V(T(4), T(0), T(0)), V(T(0), T(3), T(0)));

	auto [incenter, inradius] = t.getInscribedCircle();
	CHECK_APPROX(inradius, T(1)); // area/semi-perimeter = 6/6
	CHECK_APPROX(incenter, V(T(1), T(1), T(0)));
	auto [degenerateCenter, degenerateRadius] = Tr(V::ZERO, V::UNIT_X, V::UNIT_X*T(2)).getInscribedCircle();
	CHECK(degenerateRadius == T(0));

	auto [circumcenter, circumradius] = t.getCircumscribedCircle();
	CHECK_APPROX(circumcenter, V(T(2), T(1.5), T(0))); // midpoint of the hypotenuse
	CHECK_APPROX(circumradius, T(2.5));

	Sphere<T> sphere = t.getCircumscribedSphere();
	CHECK_APPROX(sphere.center, V(T(2), T(1.5), T(0)));
	CHECK_APPROX(sphere.radius, T(2.5));

	CHECK(t.getCircumscribedBox() == AxisAlignedBox<T>(V(T(0), T(0), T(0)), V(T(4), T(3), T(0))));
}

TEST_CASE_TEMPLATE("transformation", T, TEST_REAL_TYPES)
{
	using Tr = Triangle3<T>;
	using V = Vector3<T>;
	Tr t(V(T(0), T(0), T(0)), V(T(4), T(0), T(0)), V(T(0), T(3), T(0)));
	CHECK(Tr(t).translate(V::ONE) == Tr(V::ONE, V(T(5), T(1), T(1)), V(T(1), T(4), T(1))));
	Matrix3<T> m = Matrix3<T>::makeRotation(Axis::X, T(0.4));
	CHECK_APPROX(Tr(t).transform(m), Tr(t.vertices[0]*m, t.vertices[1]*m, t.vertices[2]*m));
	AffineTransform<T> a = AffineTransform<T>::makeRotationTranslation(Axis::X, T(0.4), V::ONE);
	CHECK_APPROX(Tr(t).transform(a), Tr(transform(t.vertices[0], a), transform(t.vertices[1], a), transform(t.vertices[2], a)));
}

TEST_CASE_TEMPLATE("evaluation, barycentric coordinates and distance", T, TEST_REAL_TYPES)
{
	using Tr = Triangle3<T>;
	using V = Vector3<T>;
	Tr t(V(T(0), T(0), T(0)), V(T(4), T(0), T(0)), V(T(0), T(4), T(0)));
	CHECK(t.evaluate(T(0.25), T(0.5)) == V(T(1), T(2), T(0)));
	CHECK(t.evaluate(T(0.25), T(0.25), T(0.5)) == V(T(1), T(2), T(0)));

	for (V p : { V(T(1), T(2), T(0)), V(T(0), T(0), T(0)), V(T(5), T(1), T(0)), V(T(-1), T(-2), T(0)) })
	{
		CAPTURE(p);
		V bary = t.getBarycentricCoords(p);
		CHECK_APPROX(sum(bary), T(1));
		CHECK_APPROX(t.evaluate(bary.x, bary.y, bary.z), p);
	}
	CHECK_APPROX(t.getBarycentricCoords(V(T(1), T(2), T(7))), V(T(0.25), T(0.25), T(0.5))); // projected onto the plane

	CHECK_APPROX(t.getClosestPoint(V(T(1), T(1), T(3))), V(T(1), T(1), T(0)));
	CHECK_APPROX(t.getClosestPoint(V(T(-2), T(-2), T(0))), V(T(0), T(0), T(0)));
	CHECK_APPROX(t.getClosestPoint(V(T(4), T(4), T(0))), V(T(2), T(2), T(0)));
	CHECK_APPROX(t.getDistanceTo(V(T(1), T(1), T(3))), T(3));
	CHECK_APPROX(t.getDistanceTo(V(T(4), T(4), T(0))), std::sqrt(T(8)));
}

TEST_CASE_TEMPLATE("intersections", T, TEST_REAL_TYPES)
{
	using Tr = Triangle3<T>;
	using V = Vector3<T>;
	Tr t(V(T(0), T(0), T(0)), V(T(4), T(0), T(0)), V(T(0), T(4), T(0)));

	CHECK(t.intersects(HalfSpace<T>(V::UNIT_X, T(-1))));		// partially inside x <= 1
	CHECK_FALSE(t.intersects(HalfSpace<T>(V::UNIT_X, T(1))));	// x <= -1
	CHECK(t.intersects(Plane<T>(V::UNIT_X, T(-1))));
	CHECK(t.intersects(Plane<T>(V::UNIT_Z, T(0))));				// coplanar
	CHECK_FALSE(t.intersects(Plane<T>(V::UNIT_X, T(-5))));
	CHECK_FALSE(t.intersects(Plane<T>(V::UNIT_Z, T(-1))));

	CHECK(t.intersects(AxisAlignedBox<T>(V(T(1), T(1), T(-1)), V(T(2), T(2), T(1)))));
	CHECK(t.intersects(AxisAlignedBox<T>(V(T(-1)), V(T(10))))); // containing box
	CHECK_FALSE(t.intersects(AxisAlignedBox<T>(V(T(3), T(3), T(-1)), V(T(4), T(4), T(1)))));
	CHECK_FALSE(t.intersects(AxisAlignedBox<T>(V(T(1), T(1), T(1)), V(T(2), T(2), T(2)))));

	Matrix3<T> r = Matrix3<T>::makeRotation(Axis::Z, Constants<T>::PI/T(4));
	CHECK(t.intersects(OrientedBox<T>(V(T(1), T(1), T(0)), r, V(T(0.5)))));
	CHECK_FALSE(t.intersects(OrientedBox<T>(V(T(4), T(4), T(0)), r, V(T(0.5)))));

	CHECK(t.intersects(Sphere<T>(V(T(1), T(1), T(0.5)), T(1))));
	CHECK(t.intersects(Sphere<T>(V(T(3), T(3), T(0)), T(1.5))));
	CHECK_FALSE(t.intersects(Sphere<T>(V(T(3), T(3), T(0)), T(1))));
}

TEST_CASE_TEMPLATE("triangulation", T, TEST_REAL_TYPES)
{
	using V = Vector3<T>;

	SUBCASE("triangle")
	{
		std::vector<V> vertices{ V(T(0), T(0), T(0)), V(T(1), T(0), T(0)), V(T(0), T(1), T(0)) };
		std::vector<int> indices;
		auto [out, result] = Triangle3<T>::triangulate(vertices.begin(), vertices.end(), std::back_inserter(indices));
		CHECK_FALSE(result);
		CHECK(indices == std::vector<int>{ 0, 1, 2 });
	}

	SUBCASE("quadrilateral")
	{
		std::vector<V> vertices{ V(T(0), T(0), T(0)), V(T(3), T(0), T(0)), V(T(3), T(1), T(0)), V(T(0), T(1), T(0)) };
		std::vector<std::uint32_t> indices;
		auto [out, result] = Triangle3<T>::template triangulate<std::uint32_t>(vertices.begin(), vertices.end(), std::back_inserter(indices));
		CHECK(result);
		REQUIRE(indices.size() == 6);
		CHECK_APPROX(triangulatedArea(vertices, indices), T(3));
	}

	SUBCASE("convex polygon")
	{
		std::vector<V> vertices;
		for (int i = 0; i < 7; i++)
		{
			T angle = Constants<T>::TWO_PI*T(i)/T(7);
			vertices.emplace_back(std::cos(angle), std::sin(angle), T(2));
		}
		std::vector<std::size_t> indices;
		auto [out, result] = Triangle3<T>::triangulate(vertices.begin(), vertices.end(), std::back_inserter(indices));
		CHECK(result);
		REQUIRE(indices.size() == 15);
		CHECK_APPROX(triangulatedArea(vertices, indices), T(3.5)*std::sin(Constants<T>::TWO_PI/T(7))); // (n/2)*sin(2*pi/n)
		for (std::size_t i = 0; i < indices.size(); i += 3) // consistent counter-clockwise winding
			CHECK(dot(Triangle3<T>(vertices[indices[i]], vertices[indices[i + 1]], vertices[indices[i + 2]]).getNormal(), V::UNIT_Z) > T(0));
	}

	SUBCASE("large polygon")
	{
		// More vertices than the triangulation's stack buffers hold
		const int n = 100;
		std::vector<V> vertices;
		for (int i = 0; i < n; i++)
		{
			T angle = Constants<T>::TWO_PI*T(i)/T(n);
			vertices.emplace_back(T(2)*std::cos(angle), T(0), -T(2)*std::sin(angle)); // counter-clockwise around +y
		}
		std::vector<int> indices;
		auto [out, result] = Triangle3<T>::triangulate(vertices.begin(), vertices.end(), std::back_inserter(indices));
		CHECK(result);
		REQUIRE(indices.size() == std::size_t(3*(n - 2)));
		CHECK_APPROX_TOL(triangulatedArea(vertices, indices), T(2*n)*std::sin(Constants<T>::TWO_PI/T(n)), T(1e-3)); // (n/2)*r^2*sin(2*pi/n)

		std::vector<int> polygon(n);
		for (int i = 0; i < n; i++)
			polygon[i] = n - 1 - i; // indexed in reverse order
		std::vector<int> indexed(3*(n - 2));
		auto [indexedOut, indexedResult] = Triangle3<T>::template triangulate<int>(vertices.begin(), vertices.end(), polygon.begin(),
			polygon.end(), indexed.begin());
		CHECK(indexedResult);
		CHECK(indexedOut == indexed.end());
		CHECK_APPROX_TOL(triangulatedArea(vertices, indexed), T(2*n)*std::sin(Constants<T>::TWO_PI/T(n)), T(1e-3));
	}

	SUBCASE("concave polygon in a tilted plane")
	{
		// An L-shaped hexagon with area 3
		std::vector<V> flat{ V(T(0), T(0), T(0)), V(T(2), T(0), T(0)), V(T(2), T(1), T(0)), V(T(1), T(1), T(0)), V(T(1), T(2), T(0)),
			V(T(0), T(2), T(0)) };
		Matrix3<T> tilt = Matrix3<T>::makeRotation(V(T(1), T(1), T(0)), T(0.7));
		std::vector<V> vertices;
		for (const V& v : flat)
			vertices.push_back(v*tilt);
		std::vector<int> indices;
		auto [out, result] = Triangle3<T>::triangulate(vertices.begin(), vertices.end(), std::back_inserter(indices));
		CHECK(result);
		REQUIRE(indices.size() == 12);
		CHECK_APPROX(triangulatedArea(vertices, indices), T(3));
	}

	SUBCASE("indexed polygon")
	{
		std::vector<V> vertices{ V(T(9), T(9), T(9)), V(T(0), T(0), T(0)), V(T(2), T(0), T(0)), V(T(2), T(1), T(0)), V(T(1), T(1), T(0)),
			V(T(1), T(2), T(0)), V(T(0), T(2), T(0)) };
		const int polygon[] = { 1, 2, 3, 4, 5, 6 };
		int indices[12];
		auto [out, result] = Triangle3<T>::template triangulate<int>(vertices.begin(), vertices.end(), std::begin(polygon), std::end(polygon),
			&indices[0]);
		CHECK(result);
		CHECK(out == std::end(indices));
		std::vector<int> list(std::begin(indices), std::end(indices));
		CHECK(std::find(list.begin(), list.end(), 0) == list.end()); // the unused vertex is not referenced
		CHECK_APPROX(triangulatedArea(vertices, list), T(3));

		const int quad[] = { 1, 2, 3, 6 };
		int quadIndices[6];
		auto [quadOut, quadResult] = Triangle3<T>::template triangulate<int>(vertices.begin(), vertices.end(), std::begin(quad), std::end(quad),
			&quadIndices[0]);
		CHECK(quadResult);
		CHECK_APPROX(triangulatedArea(vertices, std::vector<int>(std::begin(quadIndices), std::end(quadIndices))), T(3));
	}
}

TEST_CASE_TEMPLATE("streaming and hashing", T, TEST_REAL_TYPES)
{
	using Tr = Triangle3<T>;
	using V = Vector3<T>;
	Tr t(V(T(1)), V(T(2)), V(T(3.5)));
	std::stringstream s;
	s << t;
	Tr read;
	s >> read;
	CHECK(read == t);
	std::unordered_set<Tr> set{ t, t, Tr() };
	CHECK(set.size() == 2);
}

} // TEST_SUITE

TEST_SUITE("Quadrilateral3")
{

TEST_CASE_TEMPLATE("construction and properties", T, TEST_REAL_TYPES)
{
	using Q = Quadrilateral3<T>;
	using V = Vector3<T>;
	V a(T(0), T(0), T(0)), b(T(3), T(0), T(0)), c(T(3), T(2), T(0)), d(T(0), T(2), T(0));
	Q q(a, b, c, d);
	CHECK(q.vertices[2] == c);
	CHECK(Q(std::tuple<V, V, V, V>(a, b, c, d)) == q);
	const V v[] = { a, b, c, d };
	CHECK(Q(v) == q);
	CHECK(q.getVertices() == std::tuple<V, V, V, V>(a, b, c, d));
	std::vector<V> copied;
	q.copyVertices(std::back_inserter(copied));
	CHECK(copied == std::vector<V>{ a, b, c, d });
	CHECK(q.getVertex(3) == d);
	CHECK(q.getVertex(4) == V::ZERO);
	Q r = q;
	r.setVertex(3, V(T(0), T(3), T(0)));
	CHECK(r.vertices[3] == V(T(0), T(3), T(0)));
	CHECK_THROWS_AS(r.setVertex(4, V::ZERO), std::out_of_range);

	CHECK(Q(q).setZero().isZero());
	CHECK(Q().set(a, b, c, d) == q);
	CHECK(Q().isApproxZero());
	CHECK(q.isFinite());
	CHECK_FALSE(Q(a, b, c, V::INF).isFinite());
	CHECK(q.approxEquals(Q(a, b, c, d + V(T(0.05))), T(0.1)));
	CHECK_FALSE(q.approxEquals(Q(a, b, c, d + V(T(0.05)))));
	CHECK(q != r);

	CHECK_APPROX(q.getPerimeter(), T(10));
	CHECK_APPROX(q.getArea(), T(6));
	CHECK_APPROX(q.getNormal(), V::UNIT_Z);
	CHECK_APPROX(Q::computeNormal(d, c, b, a), V::MINUS_UNIT_Z);
	CHECK(Q(q).flip() == Q(d, c, b, a));
	CHECK(q.getCircumscribedBox() == AxisAlignedBox<T>(a, c));
}

TEST_CASE_TEMPLATE("transformation, streaming and hashing", T, TEST_REAL_TYPES)
{
	using Q = Quadrilateral3<T>;
	using V = Vector3<T>;
	Q q(V(T(0)), V(T(3), T(0), T(0)), V(T(3), T(2), T(0)), V(T(0), T(2), T(0)));
	CHECK(Q(q).translate(V::ONE).vertices[2] == V(T(4), T(3), T(1)));
	Matrix3<T> m = Matrix3<T>::makeRotation(Axis::X, T(0.4));
	CHECK_APPROX(Q(q).transform(m), Q(q.vertices[0]*m, q.vertices[1]*m, q.vertices[2]*m, q.vertices[3]*m));
	AffineTransform<T> a = AffineTransform<T>::makeRotationTranslation(Axis::X, T(0.4), V::ONE);
	CHECK_APPROX(Q(q).transform(a), Q(transform(q.vertices[0], a), transform(q.vertices[1], a), transform(q.vertices[2], a),
		transform(q.vertices[3], a)));

	std::stringstream s;
	s << q;
	Q read;
	s >> read;
	CHECK(read == q);
	std::unordered_set<Q> set{ q, q, Q() };
	CHECK(set.size() == 2);
}

} // TEST_SUITE

TEST_SUITE("Tetrahedron")
{

TEST_CASE_TEMPLATE("construction and properties", T, TEST_REAL_TYPES)
{
	using Te = Tetrahedron<T>;
	using V = Vector3<T>;
	V a(T(0), T(0), T(0)), b(T(2), T(0), T(0)), c(T(0), T(2), T(0)), d(T(0), T(0), T(2));
	Te t(a, b, c, d);
	CHECK(Te(std::tuple<V, V, V, V>(a, b, c, d)) == t);
	const V v[] = { a, b, c, d };
	CHECK(Te(v) == t);
	CHECK(t.getVertices() == std::tuple<V, V, V, V>(a, b, c, d));
	std::vector<V> copied;
	t.copyVertices(std::back_inserter(copied));
	CHECK(copied == std::vector<V>{ a, b, c, d });
	CHECK(t.getVertex(3) == d);
	CHECK(t.getVertex(4) == V::ZERO);
	Te u = t;
	u.setVertex(3, V(T(0), T(0), T(3)));
	CHECK(u.vertices[3] == V(T(0), T(0), T(3)));
	CHECK_THROWS_AS(u.setVertex(4, V::ZERO), std::out_of_range);

	CHECK(Te(t).setZero().isZero());
	CHECK(Te().set(a, b, c, d) == t);
	CHECK(Te().isApproxZero());
	CHECK(t.isFinite());
	CHECK_FALSE(Te(a, b, c, V::INF).isFinite());
	CHECK(t.approxEquals(Te(a, b, c, d + V(T(0.05))), T(0.1)));
	CHECK_FALSE(t.approxEquals(Te(a, b, c, d + V(T(0.05)))));
	CHECK(t != u);

	CHECK_APPROX(t.getCentroid(), V(T(0.5)));
	CHECK_APPROX(t.getVolume(), T(8)/T(6));
	CHECK_APPROX(Te(t).flip().getVolume(), T(8)/T(6));
	CHECK(Te(t).flip() == Te(b, a, c, d));
	CHECK(t.getCircumscribedBox() == AxisAlignedBox<T>(V::ZERO, V(T(2))));
}

TEST_CASE_TEMPLATE("transformation, streaming and hashing", T, TEST_REAL_TYPES)
{
	using Te = Tetrahedron<T>;
	using V = Vector3<T>;
	Te t(V(T(0)), V(T(2), T(0), T(0)), V(T(0), T(2), T(0)), V(T(0), T(0), T(2)));
	CHECK(Te(t).translate(V::ONE).vertices[3] == V(T(1), T(1), T(3)));
	Matrix3<T> m = Matrix3<T>::makeRotation(Axis::X, T(0.4));
	CHECK_APPROX(Te(t).transform(m), Te(t.vertices[0]*m, t.vertices[1]*m, t.vertices[2]*m, t.vertices[3]*m));
	CHECK_APPROX(Te(t).transform(m).getVolume(), t.getVolume());
	AffineTransform<T> a = AffineTransform<T>::makeScalingTranslation(V(T(2)), V::ONE);
	CHECK_APPROX(Te(t).transform(a).getVolume(), t.getVolume()*T(8));

	std::stringstream s;
	s << t;
	Te read;
	s >> read;
	CHECK(read == t);
	std::unordered_set<Te> set{ t, t, Te() };
	CHECK(set.size() == 2);
}

} // TEST_SUITE

} // namespace mathematics::templates
