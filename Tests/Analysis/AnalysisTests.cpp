#include "TestUtilities.hpp"
#include <vector>
#include <random>
#include <cstdint>
#include <Mathematics/Analysis/Perlin.inl>
#include <Mathematics/Analysis/Romberg.inl>
#include <Mathematics/Analysis/Distances.inl>

namespace mathematics::templates {

namespace {

// Reference: distance from a point to a convex planar polygon (counter-clockwise vertices)
template<typename T>
T distanceToConvexPolygon(const Vector3<T>& p, const std::vector<Vector3<T>>& polygon)
{
	Vector3<T> normal = normalize(cross(polygon[1] - polygon[0], polygon[2] - polygon[0]));
	Vector3<T> projected = p - dot(p - polygon[0], normal)*normal;
	bool inside = true;
	for (std::size_t i = 0; i < polygon.size(); i++)
		inside = inside && (dot(cross(polygon[(i + 1) % polygon.size()] - polygon[i], projected - polygon[i]), normal) >= T(0));
	if (inside)
		return distance(p, projected);

	T best = std::numeric_limits<T>::max();
	for (std::size_t i = 0; i < polygon.size(); i++)
		best = std::min(best, LineSegment3<T>(polygon[i], polygon[(i + 1) % polygon.size()]).getDistanceTo(p));
	return best;
}

template<typename T>
std::vector<Vector3<T>> samplePoints(std::uint32_t seed, int count, T extent)
{
	std::mt19937 g(seed);
	std::uniform_real_distribution<double> d(-double(extent), double(extent));
	std::vector<Vector3<T>> points;
	for (int i = 0; i < count; i++)
		points.emplace_back(T(d(g)), T(d(g)), T(d(g)));
	return points;
}

} // namespace

TEST_SUITE("Perlin")
{

TEST_CASE_TEMPLATE("lattice points", T, TEST_REAL_TYPES)
{
	// Gradient noise vanishes at integer coordinates
	for (int i = -3; i < 4; i++)
	{
		CHECK(perlin::noise(T(i)) == T(0));
		CHECK(perlin::noise(T(i), T(i*2)) == T(0));
		CHECK(perlin::noise(T(i), T(-i), T(i + 7)) == T(0));
	}
}

TEST_CASE_TEMPLATE("range, continuity and periodicity", T, TEST_REAL_TYPES)
{
	std::mt19937 g(42);
	std::uniform_real_distribution<double> d(-50.0, 50.0);
	bool nonZero = false;
	for (int i = 0; i < 500; i++)
	{
		T x = T(d(g)), y = T(d(g)), z = T(d(g));
		T n1 = perlin::noise(x), n2 = perlin::noise(x, y), n3 = perlin::noise(x, y, z);
		CHECK(std::fabs(n1) <= T(1.01));
		CHECK(std::fabs(n2) <= T(1.01));
		CHECK(std::fabs(n3) <= T(1.01));
		nonZero = nonZero || (std::fabs(n3) > T(0.01));

		const T h = T(1e-3);
		CHECK(std::fabs(perlin::noise(x + h, y, z) - n3) < T(0.02));
		CHECK(std::fabs(perlin::noise(x, y + h) - n2) < T(0.02));

		// The permutation table repeats every 256 units (exact shifts keep the fractional part unchanged)
		T px = std::floor(x) + T(0.375), py = std::floor(y) + T(0.625), pz = std::floor(z) + T(0.125);
		CHECK(perlin::noise(px + T(256), py, pz) == perlin::noise(px, py, pz));
		CHECK(perlin::noise(px, py - T(256)) == perlin::noise(px, py));
		CHECK(perlin::noise(px + T(512)) == perlin::noise(px));
	}
	CHECK(nonZero);

	CHECK(perlin::noise(Vector2<T>(T(1.5), T(2.25))) == perlin::noise(T(1.5), T(2.25)));
	CHECK(perlin::noise(Vector3<T>(T(1.5), T(2.25), T(-0.75))) == perlin::noise(T(1.5), T(2.25), T(-0.75)));
}

} // TEST_SUITE

TEST_SUITE("Romberg")
{

TEST_CASE_TEMPLATE("integration", T, TEST_REAL_TYPES)
{
	CHECK_APPROX(romberg::estimate<8>(T(0), T(2), [](T x) { return x*x; }), T(8)/T(3));
	CHECK_APPROX(romberg::estimate<6>(T(-1), T(3), [](T x) { return T(3)*x*x*x - x + T(2); }), T(64));
	CHECK_APPROX(romberg::estimate<10>(T(0), Constants<T>::PI, [](T x) { return std::sin(x); }), T(2));
	CHECK_APPROX(romberg::estimate<10>(T(1), T(2), [](T x) { return T(1)/x; }), std::log(T(2)));
	CHECK_APPROX(romberg::estimate<1>(T(0), T(1), [](T x) { return x; }), T(0.5)); // plain trapezoid rule
	CHECK(romberg::estimate<4>(T(1), T(1), [](T x) { return x; }) == T(0));
}

} // TEST_SUITE

TEST_SUITE("Distances")
{

TEST_CASE_TEMPLATE("planes and half spaces", T, TEST_REAL_TYPES)
{
	using V = Vector3<T>;
	V n(T(0), T(3), T(4)); // |n| = 5
	V p(T(1), T(2), T(2)); // dot(n, p) = 14
	CHECK_APPROX(distances::getPointPlane(p, n, T(-4)), T(2));
	CHECK_APPROX(distances::getPointPlane(p, n, T(-24)), T(2));
	CHECK_APPROX(distances::getPointPlaneSigned(p, n, T(-24)), T(-2));
	CHECK_APPROX(distances::getPointPlaneSquared(p, n, T(-24)), T(4));
	CHECK_APPROX(distances::getPointHalfSpace(p, n, T(-4)), T(2));
	CHECK(distances::getPointHalfSpace(p, n, T(-24)) == T(0)); // inside

	V u = n/T(5);
	CHECK_APPROX(distances::getPointNormalizedPlane(p, u, T(-0.8)), T(2));
	CHECK_APPROX(distances::getPointNormalizedPlaneSigned(p, u, T(-4.8)), T(-2));
	CHECK_APPROX(distances::getPointNormalizedPlaneSquared(p, u, T(-4.8)), T(4));
	CHECK_APPROX(distances::getPointNormalizedHalfSpace(p, u, T(-0.8)), T(2));
	CHECK(distances::getPointNormalizedHalfSpace(p, u, T(-4.8)) == T(0));
}

TEST_CASE_TEMPLATE("point and triangle", T, TEST_REAL_TYPES)
{
	using V = Vector3<T>;
	const std::vector<std::vector<V>> triangles{
		{ V(T(0), T(0), T(0)), V(T(4), T(0), T(0)), V(T(0), T(3), T(0)) },			// right angle
		{ V(T(-1), T(0), T(0)), V(T(5), T(0.5), T(0)), V(T(4), T(1), T(0.5)) },		// obtuse
		{ V(T(1), T(2), T(3)), V(T(-2), T(1), T(0)), V(T(0), T(-3), T(2)) } };		// arbitrary orientation
	for (const auto& tri : triangles)
	{
		for (const V& p : samplePoints<T>(7, 300, T(8)))
		{
			CAPTURE(p);
			V closest;
			T d2 = distances::getPointTriangleSquared(p, tri[0], tri[1], tri[2], &closest);
			T reference = distanceToConvexPolygon(p, tri);
			CHECK_APPROX_TOL(std::sqrt(d2), reference, T(1e-3));
			CHECK_APPROX_TOL(distance(p, closest), reference, T(1e-3));
			CHECK_APPROX_TOL(distances::getPointTriangle(p, tri[0], tri[1], tri[2]), reference, T(1e-3));
		}
	}

	// A point on the triangle has zero distance
	V on(T(1), T(1), T(0));
	CHECK(distances::getPointTriangleSquared(on, triangles[0][0], triangles[0][1], triangles[0][2]) == T(0));
}

TEST_CASE_TEMPLATE("point and symmetric frustum", T, TEST_REAL_TYPES)
{
	using V = Vector3<T>;
	const std::vector<SymmetricFrustum<T>> frustums{
		SymmetricFrustum<T>(V::ZERO, Matrix3<T>::IDENTITY, Vector2<T>(T(1), T(0.5)), T(1), T(4)),
		SymmetricFrustum<T>(V(T(1), T(-2), T(0.5)), Matrix3<T>::makeRotation(V(T(1), T(2), T(-1)), T(0.9)), Vector2<T>(T(0.3), T(0.8)),
			T(0.5), T(2)) };
	for (const SymmetricFrustum<T>& f : frustums)
	{
		std::vector<V> vertices = f.getVertices();
		auto quads = f.template getPrimitives<std::uint32_t>(4);
		for (const V& p : samplePoints<T>(11, 300, T(6)))
		{
			CAPTURE(p);
			T reference = T(0);
			if (!f.contains(p))
			{
				reference = std::numeric_limits<T>::max();
				for (const std::uint32_t* q = quads.first; q != quads.second; q += 4)
					reference = std::min(reference, distanceToConvexPolygon(p, { vertices[q[0]], vertices[q[1]], vertices[q[2]], vertices[q[3]] }));
			}

			V closest;
			T d2 = distances::getPointSymmetricFrustumSquared(p, f.origin, f.basis, f.halfDims, f.depthRange.minimum, f.depthRange.maximum, &closest);
			CHECK_APPROX_TOL(std::sqrt(d2), reference, T(1e-3));
			CHECK_APPROX_TOL(distance(p, closest), reference, T(1e-3));
			CHECK_APPROX_TOL(f.getDistanceTo(p), reference, T(1e-3));
		}
	}
}

} // TEST_SUITE

TEST_SUITE("Intersections")
{

TEST_CASE_TEMPLATE("helper predicates", T, TEST_REAL_TYPES)
{
	using V2 = Vector2<T>;
	using V = Vector3<T>;
	CHECK(intersections::testLineLine(V2::ZERO, V2::UNIT_X, V2::UNIT_Y, V2::UNIT_Y));
	CHECK_FALSE(intersections::testLineLine(V2::ZERO, V2::UNIT_X, V2::UNIT_Y, V2::UNIT_X));
	CHECK(intersections::testLineLine(V2::ZERO, V2::UNIT_X, V2(T(5), T(0)), V2::MINUS_UNIT_X)); // coincident
	CHECK(intersections::testLinePlane(V::UNIT_X, V::UNIT_X));
	CHECK_FALSE(intersections::testLinePlane(V::UNIT_Y, V::UNIT_X));

	// Non-optional output types report a miss with infinity.
	CHECK(std::isinf(intersections::findLinePlane<T>(V::ZERO, V::UNIT_Y, V::UNIT_X, T(-1))));
	CHECK_APPROX(intersections::findLinePlane<T>(V::ZERO, V::UNIT_X, V::UNIT_X, T(-1)), T(1));
	CHECK(std::isinf(intersections::findLineLine<T>(V2::ZERO, V2::UNIT_X, V2::UNIT_Y, V2::UNIT_X)));
	CHECK_APPROX(intersections::findLineLine<T>(V2::ZERO, V2::UNIT_X, V2(T(2), T(1)), V2::UNIT_Y), T(2));
	Interval<T> i = intersections::findLineNSphere<Interval<T>>(V::ZERO, V::UNIT_X, V(T(3), T(0), T(0)), T(1));
	CHECK_APPROX(i, Interval<T>(T(2), T(4)));
	CHECK_APPROX(intersections::findNormalizedLineNSphere<Interval<T>>(V::ZERO, V::UNIT_X, V(T(3), T(1), T(0)), T(1)), Interval<T>(T(3)));
}

} // TEST_SUITE

} // namespace mathematics::templates
