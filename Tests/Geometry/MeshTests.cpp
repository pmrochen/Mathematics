#include "TestUtilities.hpp"
#include <memory>
#include <vector>
#include <iterator>
#include <cstdint>

namespace mathematics::templates {

namespace {

// Every triangle of a closed convex mesh must face away from the centroid.
template<typename T, typename U>
void checkOutwardWinding(const TriangleMesh<T, U>& mesh)
{
	Vector3<T> centroid = mesh.computeCentroid();
	for (const Triangle3<T>& t : mesh.getTriangles())
		CHECK(dot(t.getNormal(), (t.vertices[0] + t.vertices[1] + t.vertices[2])/T(3) - centroid) > T(0));
}

// A 2x2x2 box whose top consists of y <= 1 (for x <= 0) and a plane tilted by the given angle (for x >= 0). The two top
// planes meet at the ridge x = 0, y = 1, so the ridge vertices come from a nearly parallel plane pair.
template<typename T>
ConvexPolyhedron<T> makeWedge(T angle)
{
	using V = Vector3<T>;
	V tilted(std::sin(angle), std::cos(angle), T(0));
	return { HalfSpace<T>(V::MINUS_UNIT_X, T(-1)), HalfSpace<T>(V::UNIT_X, T(-1)), HalfSpace<T>(V::MINUS_UNIT_Y, T(-1)),
		HalfSpace<T>(V::UNIT_Y, T(-1)), HalfSpace<T>(tilted, V::UNIT_Y), HalfSpace<T>(V::MINUS_UNIT_Z, T(-1)),
		HalfSpace<T>(V::UNIT_Z, T(-1)) };
}

} // namespace

TEST_SUITE("ConvexPolyhedron")
{

TEST_CASE_TEMPLATE("construction and half spaces", T, TEST_REAL_TYPES)
{
	using P = ConvexPolyhedron<T>;
	using H = HalfSpace<T>;
	using V = Vector3<T>;
	AxisAlignedBox<T> box(V(T(-1)), V(T(1)));
	std::vector<H> halfSpaces = box.getHalfSpaces();

	CHECK(P().isEmpty());
	CHECK(P::getEmpty()->isEmpty());
	CHECK(P(halfSpaces).getHalfSpaceCount() == 6);
	CHECK(P(std::vector<H>(halfSpaces)).getHalfSpaceCount() == 6);
	CHECK(P(halfSpaces.data(), halfSpaces.size()).getHalfSpaceCount() == 6);
	CHECK(P(halfSpaces.begin(), halfSpaces.end()).getHalfSpaceCount() == 6);
	CHECK(P({ halfSpaces[0], halfSpaces[1] }).getHalfSpaceCount() == 2);
	std::vector<Plane<T>> planes(halfSpaces.begin(), halfSpaces.end());
	CHECK(P(planes).getHalfSpaceCount() == 6);
	CHECK(P(planes.data(), planes.size()).getHalfSpaceCount() == 6);

	std::unique_ptr<P> p(P::from(box));
	CHECK(p->halfSpaces == halfSpaces);
	std::unique_ptr<P> q(p->clone());
	CHECK(q->halfSpaces == halfSpaces);
	P copy(*p);
	P moved(std::move(copy));
	CHECK(moved.halfSpaces == halfSpaces);
	P assigned;
	assigned = moved;
	CHECK(assigned.halfSpaces == halfSpaces);
	assigned = std::move(moved);
	CHECK(assigned.halfSpaces == halfSpaces);

	CHECK(p->getHalfSpaceMask() == 0x3Full);
	CHECK(p->template getHalfSpaceMask<std::uint8_t>() == 0x3F);
	CHECK(p->getHalfSpace(1) == halfSpaces[1]);
	CHECK(p->getHalfSpace(6) == H::EMPTY);
	CHECK(p->getHalfSpaces().size() == 6);
	CHECK(p->getPlanes() == planes);
	std::vector<Plane<T>> copiedPlanes;
	p->copyPlanes(std::back_inserter(copiedPlanes));
	CHECK(copiedPlanes == planes);

	P r;
	r.assign(p.get());
	CHECK(r.halfSpaces == halfSpaces);
	r.append(p.get());
	CHECK(r.getHalfSpaceCount() == 12);
	r.assign(nullptr);
	CHECK(r.isEmpty());
	r.addHalfSpace(halfSpaces[0]);
	r.addHalfSpaces(halfSpaces.data() + 1, 2);
	r.addHalfSpaces(std::vector<H>{ halfSpaces[3] });
	r.addHalfSpaces(halfSpaces.begin() + 4, halfSpaces.end());
	CHECK(r.halfSpaces == halfSpaces);
	r.setHalfSpace(0, halfSpaces[5]);
	CHECK(r.halfSpaces[0] == halfSpaces[5]);
	CHECK_THROWS_AS(r.setHalfSpace(6, halfSpaces[0]), std::out_of_range);
	r.insertHalfSpace(6, halfSpaces[0]);
	CHECK(r.getHalfSpaceCount() == 7);
	CHECK_THROWS_AS(r.insertHalfSpace(8, halfSpaces[0]), std::out_of_range);
	r.removeHalfSpaceAt(0);
	r.removeHalfSpaceAt(100); // ignored
	CHECK(r.getHalfSpaceCount() == 6);
	r.setHalfSpaceCount(3);
	CHECK(r.getHalfSpaceCount() == 3);
	r.removeAllHalfSpaces();
	CHECK(r.isEmpty());
	r.addHalfSpace(halfSpaces[0]);
	r.clear();
	CHECK(r.isEmpty());

	std::unique_ptr<P> fromOrientedBox(P::from(OrientedBox<T>(box)));
	CHECK(fromOrientedBox->getHalfSpaceCount() == 6);
	std::unique_ptr<P> fromFrustum(P::from(SymmetricFrustum<T>(V::ZERO, Matrix3<T>::IDENTITY, Vector2<T>(T(1)), T(1), T(2))));
	CHECK(fromFrustum->getHalfSpaceCount() == 6);
}

TEST_CASE_TEMPLATE("vertices and bounding box", T, TEST_REAL_TYPES)
{
	using P = ConvexPolyhedron<T>;
	using V = Vector3<T>;
	AxisAlignedBox<T> box(V(T(-1), T(-2), T(-3)), V(T(1), T(2), T(3)));
	std::unique_ptr<P> p(P::from(box));
	std::vector<V> vertices = p->computeVertices();
	CHECK(vertices.size() == 8);
	for (const V& v : box.getVertices())
		CHECK(std::any_of(vertices.begin(), vertices.end(), [&v](const V& u) { return u.approxEquals(v, test::tolerance<T>()); }));
	CHECK_APPROX(p->computeAxisAlignedBoundingBox(), box);

	// Cutting off a corner adds three vertices and removes one
	p->addHalfSpace(HalfSpace<T>(normalize(V(T(1), T(1), T(1))), -T(5)/std::sqrt(T(3)))); // x + y + z <= 5
	CHECK(p->computeVertices().size() == 10);

	SUBCASE("nearly parallel planes")
	{
		const T angle = T(0.01);
		ConvexPolyhedron<T> wedge = makeWedge(angle);
		std::vector<V> vertices = wedge.computeVertices();
		CHECK(vertices.size() == 10); // the default pivot tolerance keeps the ridge vertices of a 0.57 degree wedge
		for (T z : { T(-1), T(1) })
		{
			CAPTURE(z);
			CHECK(std::any_of(vertices.begin(), vertices.end(), [z](const V& v) { return v.approxEquals(V(T(0), T(1), z), T(1e-4)); }));
			CHECK(std::any_of(vertices.begin(), vertices.end(), [&](const V& v) { return v.approxEquals(V(T(1), T(1) - std::tan(angle), z), T(1e-4)); }));
		}
		CHECK_APPROX(wedge.computeAxisAlignedBoundingBox(), AxisAlignedBox<T>(V(T(-1)), V(T(1))));

		// A pivot tolerance above the angle between the planes drops the ridge vertices.
		CHECK(wedge.computeVertices(Constants<T>::TOLERANCE*T(1000), T(0.1)).size() == 8);
		CHECK(wedge.computeAxisAlignedBoundingBox(Constants<T>::TOLERANCE*T(1000), T(0.1)).isFinite());

		// Default pivot tolerance of sqrt(epsilon): exactly parallel planes never generate vertices.
		ConvexPolyhedron<T> box = makeWedge(T(0));
		CHECK(box.computeVertices().size() == 8);
	}
}

TEST_CASE_TEMPLATE("transformation", T, TEST_REAL_TYPES)
{
	using P = ConvexPolyhedron<T>;
	using V = Vector3<T>;
	AxisAlignedBox<T> box(V(T(-1)), V(T(1)));

	std::unique_ptr<P> t(P::from(box));
	t->translate(V(T(5), T(0), T(0)));
	CHECK(t->contains(V(T(5.5), T(0), T(0))));
	CHECK_FALSE(t->contains(V::ZERO));
	t->translate(V::ZERO);

	Matrix3<T> m = Matrix3<T>::makeScalingRotation(V(T(2), T(1), T(1)), Axis::Z, Constants<T>::PI/T(6));
	std::unique_ptr<P> s(P::from(box));
	s->transform(m);
	for (const V& v : box.getVertices())
		CHECK(s->contains(v*m*T(0.999)));
	CHECK_FALSE(s->contains(V(T(1.5), T(0), T(0))*m*T(1.1)));
	s->transform(Matrix3<T>::IDENTITY); // no-op

	Matrix3<T> r = Matrix3<T>::makeRotation(Axis::X, T(0.4));
	std::unique_ptr<P> o(P::from(box));
	o->transform(r, true);
	for (const V& v : box.getVertices())
		CHECK(o->contains(v*r*T(0.999)));

	AffineTransform<T> a = AffineTransform<T>::makeScalingRotationTranslation(V(T(1), T(3), T(1)), Axis::Y, T(0.3), V(T(1), T(2), T(3)));
	std::unique_ptr<P> pa(P::from(box));
	pa->transform(a);
	for (const V& v : box.getVertices())
		CHECK(pa->contains(transform(v*T(0.999), a)));
	CHECK_FALSE(pa->contains(transform(V(T(1.1), T(0), T(0)), a)));

	std::unique_ptr<P> pt(P::from(box));
	pt->transform(AffineTransform<T>::makeTranslation(V(T(5), T(0), T(0))));
	CHECK(pt->contains(V(T(5.5), T(0), T(0))));
	pt->transform(AffineTransform<T>::makeRotationTranslation(Axis::Z, Constants<T>::HALF_PI, V::ZERO), true);
	CHECK(pt->contains(V(T(0), T(5.5), T(0))));
	pt->transform(AffineTransform<T>::IDENTITY); // no-op
	CHECK(pt->contains(V(T(0), T(5.5), T(0))));
}

TEST_CASE_TEMPLATE("containment and intersections", T, TEST_REAL_TYPES)
{
	using P = ConvexPolyhedron<T>;
	using V = Vector3<T>;
	std::unique_ptr<P> p(P::from(AxisAlignedBox<T>(V(T(-1)), V(T(1)))));
	CHECK(p->contains(V::ZERO));
	CHECK(p->contains(V(T(1))));
	CHECK_FALSE(p->contains(V(T(1.01), T(0), T(0))));
	CHECK(p->contains(V(T(1.01), T(0), T(0)), T(0.1)));
	CHECK_FALSE(p->contains(V(T(1.2), T(0), T(0)), T(0.1)));

	CHECK(p->intersects(AxisAlignedBox<T>(V(T(0.5)), V(T(2)))));
	CHECK_FALSE(p->intersects(AxisAlignedBox<T>(V(T(1.5)), V(T(2)))));
	CHECK(p->intersects(OrientedBox<T>(V(T(1.5), T(0), T(0)), Matrix3<T>::makeRotation(Axis::Z, T(0.5)), V(T(1)))));
	CHECK_FALSE(p->intersects(OrientedBox<T>(V(T(3), T(0), T(0)), Matrix3<T>::makeRotation(Axis::Z, T(0.5)), V(T(1)))));
	CHECK(p->intersects(Sphere<T>(V(T(1.5), T(0), T(0)), T(1))));
	CHECK_FALSE(p->intersects(Sphere<T>(V(T(3), T(0), T(0)), T(1))));

	CHECK(Sphere<T>(V(T(1.5), T(0), T(0)), T(1)).intersects(p.get()));
	CHECK_FALSE(Sphere<T>(V(T(1.5), T(0), T(0)), T(1)).intersects(static_cast<const P*>(nullptr)));
	CHECK(AxisAlignedBox<T>(V(T(0.5)), V(T(2))).intersects(p.get()));
	CHECK(OrientedBox<T>(V(T(1.5), T(0), T(0)), Matrix3<T>::IDENTITY, V(T(1))).intersects(p.get()));

	SUBCASE("hierarchical culling with a plane mask")
	{
		unsigned int mask = p->template getHalfSpaceMask<unsigned int>();
		CHECK(p->intersects(AxisAlignedBox<T>(V(T(-0.5)), V(T(0.5))), mask));
		CHECK(mask == 0u); // fully inside: no planes need to be tested for children

		mask = p->template getHalfSpaceMask<unsigned int>();
		CHECK(p->intersects(AxisAlignedBox<T>(V(T(0.5), T(-0.5), T(-0.5)), V(T(1.5), T(0.5), T(0.5))), mask));
		CHECK(mask == 0x2u); // only the +x plane straddles the box

		mask = p->template getHalfSpaceMask<unsigned int>();
		CHECK_FALSE(p->intersects(AxisAlignedBox<T>(V(T(2)), V(T(3))), mask));

		mask = p->template getHalfSpaceMask<unsigned int>();
		OrientedBox<T> inside(V::ZERO, Matrix3<T>::makeRotation(Axis::Z, T(0.5)), V(T(0.2)));
		CHECK(p->intersects(inside, mask));
		CHECK(mask == 0u);
		mask = p->template getHalfSpaceMask<unsigned int>();
		CHECK_FALSE(p->intersects(OrientedBox<T>(V(T(3)), Matrix3<T>::IDENTITY, V(T(0.5))), mask));
	}
}

} // TEST_SUITE

TEST_SUITE("TriangleMesh")
{

TEST_CASE_TEMPLATE("construction", T, TEST_REAL_TYPES)
{
	using M = TriangleMesh<T>;
	using V = Vector3<T>;
	std::vector<V> vertices{ V(T(0)), V(T(1), T(0), T(0)), V(T(0), T(1), T(0)), V(T(0), T(0), T(1)) };
	std::vector<std::uint32_t> indices{ 0, 2, 1, 0, 1, 3, 0, 3, 2, 1, 2, 3 };
	std::vector<int> intIndices(indices.begin(), indices.end());

	CHECK(M().isEmpty());
	CHECK(M(vertices.begin(), vertices.end()).vertices == vertices);
	CHECK(M(vertices.data(), vertices.size()).vertices == vertices);
	CHECK(M(vertices).vertices == vertices);
	CHECK(M(std::vector<V>(vertices)).vertices == vertices);
	CHECK(M({ V(T(0)), V(T(1)) }).getVertexCount() == 2);

	M a(vertices.data(), vertices.size(), indices.data(), indices.size());
	CHECK(a.vertices == vertices);
	CHECK(a.indices == indices);
	CHECK(M(vertices.data(), vertices.size(), intIndices.data(), intIndices.size()).indices == indices);
	CHECK(M(vertices, indices).indices == indices);
	CHECK(M(vertices, intIndices).indices == indices);
	CHECK(M(std::vector<V>(vertices), std::vector<std::uint32_t>(indices)).indices == indices);

	std::vector<Triangle3<T>> triangles{ Triangle3<T>(vertices[0], vertices[1], vertices[2]), Triangle3<T>(vertices[1], vertices[2], vertices[3]) };
	M b(triangles);
	CHECK(b.getVertexCount() == 6);
	CHECK(b.getTriangleCount() == 2);
	CHECK(b.getTriangle(1) == triangles[1]);
	CHECK(M(triangles.data(), triangles.size()).getVertexCount() == 6);

	M copy(a);
	CHECK(copy.indices == indices);
	M moved(std::move(copy));
	CHECK(moved.indices == indices);
	M assigned;
	assigned = moved;
	CHECK(assigned.vertices == vertices);
	assigned = std::move(moved);
	CHECK(assigned.vertices == vertices);
	std::unique_ptr<M> clone(a.clone());
	CHECK(clone->indices == indices);
}

TEST_CASE_TEMPLATE("vertices and indices", T, TEST_REAL_TYPES)
{
	using M = TriangleMesh<T>;
	using V = Vector3<T>;
	M m;
	m.addVertex(V(T(0)));
	m.addVertices(std::vector<V>{ V(T(1)), V(T(2)) });
	const V more[] = { V(T(3)), V(T(4)) };
	m.addVertices(more, 2);
	m.addVertices(std::begin(more), std::end(more));
	CHECK(m.getVertexCount() == 7);
	CHECK(m.getVertex(2) == V(T(2)));
	CHECK(m.getVertex(100) == V::ZERO);
	m.setVertex(0, V(T(-1)));
	CHECK(m.vertices[0] == V(T(-1)));
	CHECK_THROWS_AS(m.setVertex(7, V::ZERO), std::out_of_range);
	m.insertVertex(0, V(T(9)));
	CHECK(m.vertices[0] == V(T(9)));
	CHECK_THROWS_AS(m.insertVertex(100, V::ZERO), std::out_of_range);
	m.removeVertexAt(0);
	m.removeVertexAt(100); // ignored
	CHECK(m.getVertexCount() == 7);
	m.setVertexCount(6);
	CHECK(m.getVertexCount() == 6);

	// Without indices, every three consecutive vertices form a triangle.
	CHECK(m.getTriangleCount() == 2);
	CHECK(m.getVertexIndex(1, 2) == 5u);
	CHECK(m.computeVertexIndices());
	CHECK(m.indices == std::vector<std::uint32_t>{ 0, 1, 2, 3, 4, 5 });
	CHECK_FALSE(m.computeVertexIndices()); // already indexed

	m.addVertexIndex(0);
	m.addVertexIndices(1u, 2u, 3u);
	const std::uint32_t moreIndices[] = { 4, 5 };
	m.addVertexIndices(moreIndices, 2);
	m.addVertexIndices(std::vector<std::uint32_t>{ 0, 1 });
	m.addVertexIndices(std::vector<int>{ 2 });
	m.addVertexIndices(std::begin(moreIndices), std::end(moreIndices));
	CHECK(m.getVertexIndexCount() == 17);
	CHECK(m.getVertexIndex(6) == 0u);
	CHECK(m.getVertexIndex(100) == 0u);
	m.setVertexIndex(6, 5);
	CHECK(m.indices[6] == 5u);
	CHECK_THROWS_AS(m.setVertexIndex(17, 0), std::out_of_range);
	m.insertVertexIndex(0, 3);
	CHECK(m.indices[0] == 3u);
	CHECK_THROWS_AS(m.insertVertexIndex(100, 0), std::out_of_range);
	m.insertVertexIndices(1, 4, 4, 4);
	CHECK(m.indices[3] == 4u);
	CHECK_THROWS_AS(m.insertVertexIndices(100, 0, 0, 0), std::out_of_range);
	m.removeVertexIndexAt(0);
	m.removeVertexIndexAt(1000); // ignored
	CHECK(m.getVertexIndexCount() == 20);
	CHECK(m.getTriangleCount() == 6);
	CHECK(m.hasValidVertexIndices());
	m.addVertexIndex(6);
	CHECK_FALSE(m.hasValidVertexIndices());
	m.setVertexIndexCount(3);
	CHECK(m.getVertexIndexCount() == 3);
	m.removeAllVertexIndices();
	CHECK(m.indices.empty());
	m.removeAllVertices();
	CHECK(m.isEmpty());

	M n({ V(T(0)), V(T(1)), V(T(2)) });
	n.addVertexIndices(0u, 1u, 2u);
	n.clear();
	CHECK(n.isEmpty());
}

TEST_CASE_TEMPLATE("assign and append", T, TEST_REAL_TYPES)
{
	using M = TriangleMesh<T>;
	std::unique_ptr<M> box(M::makeBox(Vector3<T>(T(2))));
	std::unique_ptr<M> tetrahedron(M::makeTetrahedron(T(1)));

	M m;
	m.assign(box.get());
	CHECK(m.vertices == box->vertices);
	CHECK(m.indices == box->indices);
	m.append(tetrahedron.get());
	CHECK(m.getVertexCount() == 12);
	CHECK(m.getTriangleCount() == 16);
	CHECK(m.hasValidVertexIndices());
	CHECK(m.getTriangle(12) == tetrahedron->getTriangle(0)); // appended indices are offset
	CHECK_APPROX(m.computeVolume(), box->computeVolume() + tetrahedron->computeVolume());
	m.append(nullptr);
	CHECK(m.getVertexCount() == 12);
	m.assign(nullptr);
	CHECK(m.isEmpty());
}

TEST_CASE_TEMPLATE("factories and measurements", T, TEST_REAL_TYPES)
{
	using M = TriangleMesh<T>;
	using V = Vector3<T>;

	SUBCASE("box")
	{
		std::unique_ptr<M> box(M::makeBox(V(T(2), T(4), T(6))));
		REQUIRE(box);
		CHECK(box->getVertexCount() == 8);
		CHECK(box->getTriangleCount() == 12);
		CHECK(box->hasValidVertexIndices());
		CHECK(box->isConvex());
		CHECK_APPROX(box->computeVolume(), T(48));
		CHECK_APPROX(box->computeSurfaceArea(), T(88));
		CHECK_APPROX(box->computeCentroid(), V::ZERO);
		CHECK_APPROX(box->computeAxisAlignedBoundingBox(), AxisAlignedBox<T>(V(T(-1), T(-2), T(-3)), V(T(1), T(2), T(3))));
		checkOutwardWinding(*box);
		CHECK(M::makeBox(V(T(1), T(0), T(1))) == nullptr);

		std::unique_ptr<M> fromBox(M::from(AxisAlignedBox<T>(V(T(-1)), V(T(1)))));
		CHECK_APPROX(fromBox->computeVolume(), T(8));
		std::unique_ptr<M> fromOrientedBox(M::from(OrientedBox<T>(V::ONE, Matrix3<T>::makeRotation(Axis::Y, T(0.3)), V(T(1)))));
		CHECK_APPROX(fromOrientedBox->computeVolume(), T(8));
		checkOutwardWinding(*fromOrientedBox);
	}

	SUBCASE("tetrahedron")
	{
		std::unique_ptr<M> t(M::makeTetrahedron(T(2)));
		REQUIRE(t);
		CHECK(t->getTriangleCount() == 4);
		CHECK(t->isConvex());
		for (const Triangle3<T>& tri : t->getTriangles())
		{
			CHECK_APPROX(distance(tri.vertices[0], tri.vertices[1]), T(2));
			CHECK_APPROX(distance(tri.vertices[1], tri.vertices[2]), T(2));
		}
		CHECK_APPROX(t->computeVolume(), T(8)/(T(6)*std::sqrt(T(2)))); // a^3/(6*sqrt(2))
		CHECK_APPROX(t->computeSurfaceArea(), std::sqrt(T(3))*T(4));
		checkOutwardWinding(*t);
		CHECK(M::makeTetrahedron(T(0)) == nullptr);
	}

	SUBCASE("frustum")
	{
		SymmetricFrustum<T> f(V::ZERO, Matrix3<T>::IDENTITY, Vector2<T>(T(1)), T(1), T(2));
		std::unique_ptr<M> m(M::from(f));
		CHECK(m->getTriangleCount() == 12);
		CHECK_APPROX(m->computeVolume(), T(28)/T(3)); // (h/3)*(A1 + A2 + sqrt(A1*A2)) = (1/3)*(4 + 16 + 8)
		checkOutwardWinding(*m);
	}

	SUBCASE("convex polyhedron")
	{
		std::unique_ptr<ConvexPolyhedron<T>> p(ConvexPolyhedron<T>::from(AxisAlignedBox<T>(V(T(-1), T(-2), T(-3)), V(T(1), T(2), T(3)))));
		std::unique_ptr<M> m(M::from(p.get()));
		REQUIRE(m);
		CHECK(m->getVertexCount() == 8);
		CHECK(m->getTriangleCount() == 12);
		CHECK(m->hasValidVertexIndices());
		CHECK_APPROX(m->computeVolume(), T(48));
		CHECK_APPROX(m->computeSurfaceArea(), T(88));
		checkOutwardWinding(*m);

		// Frustum: faces with four vertices each, non-parallel side planes
		std::unique_ptr<ConvexPolyhedron<T>> fp(ConvexPolyhedron<T>::from(SymmetricFrustum<T>(V::ZERO, Matrix3<T>::IDENTITY,
			Vector2<T>(T(1)), T(1), T(2))));
		std::unique_ptr<M> fm(M::from(fp.get()));
		CHECK(fm->getVertexCount() == 8);
		CHECK_APPROX(fm->computeVolume(), T(28)/T(3));
		checkOutwardWinding(*fm);

		// A thin wedge (two nearly parallel top planes)
		const T angle = T(0.01);
		ConvexPolyhedron<T> wedge = makeWedge(angle);
		std::unique_ptr<M> wm(M::from(&wedge));
		CHECK(wm->getVertexCount() == 10);
		CHECK(wm->hasValidVertexIndices());
		CHECK_APPROX_TOL(wm->computeVolume(), T(8) - std::tan(angle), T(1e-4));
		checkOutwardWinding(*wm);
		std::unique_ptr<M> coarse(M::from(&wedge, Constants<T>::TOLERANCE*T(1000), T(0.1)));
		CHECK(coarse->getVertexCount() == 8); // the ridge vertices are dropped

		CHECK(M::from(static_cast<const ConvexPolyhedron<T>*>(nullptr)) == nullptr);
		ConvexPolyhedron<T> open{ HalfSpace<T>(V::UNIT_X, T(0)), HalfSpace<T>(V::UNIT_Y, T(0)) };
		std::unique_ptr<M> empty(M::from(&open));
		CHECK(empty->isEmpty());
	}
}

TEST_CASE_TEMPLATE("triangles, half spaces and planes", T, TEST_REAL_TYPES)
{
	using M = TriangleMesh<T>;
	using V = Vector3<T>;
	std::unique_ptr<M> box(M::makeBox(V(T(2))));
	std::vector<Triangle3<T>> triangles = box->getTriangles();
	REQUIRE(triangles.size() == 12);
	std::vector<Triangle3<T>> copied;
	box->copyTriangles(std::back_inserter(copied));
	CHECK(copied == triangles);
	CHECK(box->getTriangle(3) == triangles[3]);

	std::vector<HalfSpace<T>> halfSpaces = box->getHalfSpaces();
	REQUIRE(halfSpaces.size() == 12);
	std::vector<HalfSpace<T>> copiedHalfSpaces;
	box->copyHalfSpaces(std::back_inserter(copiedHalfSpaces));
	CHECK(copiedHalfSpaces == halfSpaces);
	for (const HalfSpace<T>& h : halfSpaces)
		CHECK(h.contains(V::ZERO));

	std::vector<Plane<T>> planes = box->getPlanes();
	REQUIRE(planes.size() == 12);
	std::vector<Plane<T>> copiedPlanes;
	box->copyPlanes(std::back_inserter(copiedPlanes));
	CHECK(copiedPlanes == planes);
	for (std::size_t i = 0; i < 12; i++)
		CHECK(planes[i].contains(triangles[i].vertices[0]));
}

TEST_CASE_TEMPLATE("bounding sphere", T, TEST_REAL_TYPES)
{
	using M = TriangleMesh<T>;
	using V = Vector3<T>;
	CHECK(M().computeBoundingSphere() == Sphere<T>());
	CHECK(M({ V(T(1), T(2), T(3)) }).computeBoundingSphere() == Sphere<T>(V(T(1), T(2), T(3)), T(0)));
	CHECK(M().computeCentroid() == V::ZERO);
	CHECK(M().computeAxisAlignedBoundingBox().isEmpty());

	std::vector<V> points;
	for (int i = 0; i < 50; i++)
		points.emplace_back(std::sin(T(i)*T(1.3))*T(3), std::cos(T(i)*T(0.7))*T(2), std::sin(T(i)*T(0.31) + T(1))*T(1.5));
	M cloud(points);
	Sphere<T> s = cloud.computeBoundingSphere();
	for (const V& p : points)
		CHECK(distance(p, s.center) <= s.radius*(T(1) + test::tolerance<T>()));
	CHECK(s.radius <= cloud.computeAxisAlignedBoundingBox().getCircumscribedSphere().radius*T(1.05)); // reasonably tight
}

TEST_CASE_TEMPLATE("transformation", T, TEST_REAL_TYPES)
{
	using M = TriangleMesh<T>;
	using V = Vector3<T>;
	std::unique_ptr<M> box(M::makeBox(V(T(2))));
	M m(*box);
	m.translate(V(T(1), T(2), T(3)));
	CHECK_APPROX(m.computeCentroid(), V(T(1), T(2), T(3)));
	m.translate(V::ZERO);
	CHECK_APPROX(m.computeCentroid(), V(T(1), T(2), T(3)));

	M s(*box);
	s.transform(Matrix3<T>::makeScaling(V(T(2), T(1), T(3))));
	CHECK_APPROX(s.computeVolume(), T(48));
	s.transform(Matrix3<T>::IDENTITY);
	CHECK_APPROX(s.computeVolume(), T(48));

	M a(*box);
	a.transform(AffineTransform<T>::makeScalingRotationTranslation(V(T(2)), Axis::X, T(0.4), V(T(5), T(0), T(0))));
	CHECK_APPROX(a.computeVolume(), T(64));
	CHECK_APPROX(a.computeCentroid(), V(T(5), T(0), T(0)));
	a.transform(AffineTransform<T>::makeTranslation(V(T(-5), T(0), T(0))));
	CHECK_APPROX(a.computeCentroid(), V::ZERO);
	a.transform(AffineTransform<T>::IDENTITY);
	CHECK_APPROX(a.computeCentroid(), V::ZERO);
}

TEST_CASE_TEMPLATE("welding and containment", T, TEST_REAL_TYPES)
{
	using M = TriangleMesh<T>;
	using V = Vector3<T>;
	std::unique_ptr<M> box(M::makeBox(V(T(2))));
	M soup(box->getTriangles()); // 36 unshared vertices
	CHECK(soup.getVertexCount() == 36);
	CHECK(soup.weldVertices());
	CHECK(soup.getVertexCount() == 8);
	CHECK(soup.getTriangleCount() == 12);
	CHECK(soup.hasValidVertexIndices());
	CHECK_APPROX(soup.computeVolume(), T(8));
	CHECK_FALSE(soup.weldVertices()); // nothing left to weld
	CHECK_FALSE(M({ V(T(0)) }).weldVertices());

	CHECK(soup.contains(V::ZERO));
	CHECK(soup.contains(V(T(0.99))));
	CHECK_FALSE(soup.contains(V(T(1.1), T(0), T(0))));
	CHECK(soup.isConvex());

	M concave({ V(T(0), T(0), T(0)), V(T(2), T(0), T(0)), V(T(1), T(1), T(0)), V(T(1), T(0.5), T(-1)) });
	concave.indices = { 0, 1, 2, 0, 3, 1, 1, 3, 2, 2, 3, 0 };
	CHECK(concave.isConvex());
	concave.vertices[3] = V(T(1), T(0.5), T(1)); // flip the apex through the base
	CHECK_FALSE(concave.isConvex());
}

} // TEST_SUITE

TEST_SUITE("TetrahedralMesh")
{

TEST_CASE_TEMPLATE("construction, editing and measurements", T, TEST_REAL_TYPES)
{
	using M = TetrahedralMesh<T>;
	using V = Vector3<T>;
	// A unit cube split into five tetrahedra
	std::vector<V> vertices;
	for (int i = 0; i < 8; i++)
		vertices.emplace_back(T(i & 1), T((i >> 1) & 1), T((i >> 2) & 1));
	std::vector<std::uint32_t> indices{ 0, 1, 2, 4, 1, 2, 3, 7, 1, 4, 5, 7, 2, 4, 6, 7, 1, 2, 4, 7 };

	M m(vertices, indices);
	CHECK(m.getVertexCount() == 8);
	CHECK(m.getTetrahedronCount() == 5);
	CHECK(m.hasValidVertexIndices());
	CHECK_APPROX(m.computeVolume(), T(1));
	CHECK_APPROX(m.computeCentroid(), V(T(0.5)));
	CHECK(m.computeAxisAlignedBoundingBox() == AxisAlignedBox<T>(V::ZERO, V::ONE));
	Sphere<T> s = m.computeBoundingSphere();
	for (const V& v : vertices)
		CHECK(distance(v, s.center) <= s.radius*(T(1) + test::tolerance<T>()));
	CHECK(m.getTetrahedron(4) == Tetrahedron<T>(vertices[1], vertices[2], vertices[4], vertices[7]));

	std::vector<Tetrahedron<T>> tetrahedra = m.getTetrahedrons();
	REQUIRE(tetrahedra.size() == 5);
	std::vector<Tetrahedron<T>> copied;
	m.copyTetrahedrons(std::back_inserter(copied));
	CHECK(copied == tetrahedra);
	T volume = T(0);
	for (const Tetrahedron<T>& t : tetrahedra)
		volume += t.getVolume();
	CHECK_APPROX(volume, T(1));

	CHECK(M(vertices.data(), vertices.size(), indices.data(), indices.size()).indices == indices);
	CHECK(M(std::vector<V>(vertices), std::vector<std::uint32_t>(indices)).indices == indices);
	CHECK(M(vertices, std::vector<int>(indices.begin(), indices.end())).indices == indices);
	CHECK(M(tetrahedra).getVertexCount() == 20);
	CHECK_APPROX(M(tetrahedra.data(), tetrahedra.size()).computeVolume(), T(1));
	CHECK(M(vertices.begin(), vertices.end()).getVertexCount() == 8);
	CHECK(M(vertices.data(), vertices.size()).getTetrahedronCount() == 2);

	M soup(tetrahedra);
	CHECK(soup.weldVertices());
	CHECK(soup.getVertexCount() == 8);
	CHECK_APPROX(soup.computeVolume(), T(1));

	M e;
	e.addVertex(V::ZERO);
	const V axes[] = { V::UNIT_X, V::UNIT_Y, V::UNIT_Z };
	e.addVertices(axes, 3);
	e.addVertexIndices(0u, 1u, 2u, 3u);
	CHECK(e.getTetrahedronCount() == 1);
	CHECK_APPROX(e.computeVolume(), T(1)/T(6));
	e.insertVertexIndices(0, 3u, 2u, 1u, 0u);
	CHECK(e.getTetrahedronCount() == 2);
	CHECK_THROWS_AS(e.insertVertexIndices(100, 0u, 0u, 0u, 0u), std::out_of_range);
	e.clear();
	CHECK(e.isEmpty());

	M t(m);
	t.translate(V::ONE);
	CHECK_APPROX(t.computeCentroid(), V(T(1.5)));
	t.transform(Matrix3<T>::makeScaling(T(2)));
	CHECK_APPROX(t.computeVolume(), T(8));
	t.transform(AffineTransform<T>::makeScalingTranslation(V(T(0.5)), V::ZERO));
	CHECK_APPROX(t.computeVolume(), T(1));

	M a;
	a.assign(&m);
	a.append(&m);
	CHECK(a.getTetrahedronCount() == 10);
	CHECK(a.hasValidVertexIndices());
	CHECK_APPROX(a.computeVolume(), T(2));
}

} // TEST_SUITE

} // namespace mathematics::templates
