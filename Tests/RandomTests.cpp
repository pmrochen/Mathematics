#include "TestUtilities.hpp"
#include <cstdint>
#include <set>

namespace mathematics {

namespace {

// A uniform random bit generator that always yields one of its bounds; used to check the interval end points.
template<bool Max>
struct ExtremeGenerator
{
	using result_type = std::uint32_t;
	static constexpr result_type min() { return 0; }
	static constexpr result_type max() { return 0xFFFFFFFFu; }
	result_type operator()() { return Max ? max() : min(); }
};

constexpr int SAMPLES = 2000;

} // namespace

TEST_SUITE("Random")
{

TEST_CASE_TEMPLATE("scalars are within their ranges", T, TEST_REAL_TYPES)
{
	using R = Random<T>;
	std::mt19937 g(1234);
	T sum = T(0);
	for (int i = 0; i < SAMPLES; i++)
	{
		T a = R::get(T(-2), T(3), g);
		CHECK((a >= T(-2) && a <= T(3)));
		T b = R::getRightOpen(T(-2), T(3), g);
		CHECK((b >= T(-2) && b < T(3)));
		T c = R::get01(g);
		CHECK((c >= T(0) && c <= T(1)));
		T d = R::get01RightOpen(g);
		CHECK((d >= T(0) && d < T(1)));
		sum += c;
	}
	CHECK_APPROX_TOL(sum/T(SAMPLES), T(0.5), T(0.05)); // roughly uniform
}

TEST_CASE_TEMPLATE("interval end points", T, TEST_REAL_TYPES)
{
	using Low = Random<T, ExtremeGenerator<false>>;
	using High = Random<T, ExtremeGenerator<true>>;
	ExtremeGenerator<false> low;
	ExtremeGenerator<true> high;

	CHECK(Low::get(T(-2), T(3), low) == T(-2));
	CHECK(High::get(T(-2), T(3), high) == T(3));
	CHECK(Low::getRightOpen(T(-2), T(3), low) == T(-2));
	CHECK(High::getRightOpen(T(-2), T(3), high) < T(3));
	CHECK(Low::get01(low) == T(0));
	CHECK(High::get01(high) == T(1));
	CHECK(Low::get01RightOpen(low) == T(0));
	CHECK(High::get01RightOpen(high) < T(1));
}

TEST_CASE_TEMPLATE("integers", T, int, long long, unsigned int)
{
	using R = Random<T>;
	std::mt19937 g(99);
	std::set<T> seen;
	for (int i = 0; i < SAMPLES; i++)
	{
		T v = R::get(T(1), T(6), g);
		CHECK((v >= T(1) && v <= T(6)));
		seen.insert(v);
	}
	CHECK(seen.size() == 6); // every value occurs

	ExtremeGenerator<false> low;
	ExtremeGenerator<true> high;
	CHECK(Random<T, ExtremeGenerator<false>>::get(T(1), T(6), low) == T(1));
	CHECK(Random<T, ExtremeGenerator<true>>::get(T(1), T(6), high) == T(6));
}

TEST_CASE_TEMPLATE("vectors", T, TEST_REAL_TYPES)
{
	using R = Random<T>;
	using V3 = templates::Vector3<T>;
	std::mt19937 g(4321);
	for (int i = 0; i < SAMPLES/4; i++)
	{
		CHECK_APPROX(R::getVectorOnUnitCircle(g).getMagnitude(), T(1));
		CHECK(R::getVectorOnUnitDisk(g).getMagnitude() <= T(1));
		CHECK(R::getVectorInUnitSphere(g).getMagnitude() <= T(1));
		CHECK_APPROX(R::getVectorOnUnitSphere(g).getMagnitude(), T(1));

		V3 c = R::getVectorOnUnitSphericalCone(T(0.3), g);
		CHECK_APPROX(c.getMagnitude(), T(1));
		CHECK(angle(c, V3::UNIT_Z) <= T(0.3) + test::tolerance<T>());

		V3 direction = normalize(V3(T(1), T(-2), T(0.5)));
		V3 d = R::getVectorOnUnitSphericalCone(direction, T(0.2), g);
		CHECK_APPROX(d.getMagnitude(), T(1));
		CHECK(angle(d, direction) <= T(0.2) + T(10)*test::tolerance<T>());

		V3 e = R::getVectorOnUnitSphericalCone(V3::UNIT_Z, T(0.2), g); // direction parallel to the helper axis
		CHECK(angle(e, V3::UNIT_Z) <= T(0.2) + T(10)*test::tolerance<T>());
	}
}

TEST_CASE_TEMPLATE("rotations", T, TEST_REAL_TYPES)
{
	using R = Random<T>;
	std::mt19937 g(777);
	for (int i = 0; i < SAMPLES/10; i++)
	{
		templates::Matrix3<T> m = R::getRotationMatrix(g);
		CHECK((m*transpose(m)).approxEquals(templates::Matrix3<T>::IDENTITY, T(10)*test::tolerance<T>()));
		CHECK_APPROX(m.getDeterminant(), T(1));
		CHECK_APPROX(R::getQuaternion(g).getMagnitude(), T(1));
	}
}

TEST_CASE("default generator")
{
	using R = Random<float>;
	R::getDefaultGenerator().seed(5);
	float a = R::get(0.f, 1.f);
	R::getDefaultGenerator().seed(5);
	CHECK(R::get(0.f, 1.f) == a); // reproducible
	CHECK(R::getRightOpen(0.f, 1.f) < 1.f);
	CHECK(R::get01() <= 1.f);
	CHECK(R::get01RightOpen() < 1.f);
	CHECK_APPROX(R::getVectorOnUnitCircle().getMagnitude(), 1.f);
	CHECK(R::getVectorOnUnitDisk().getMagnitude() <= 1.f);
	CHECK(R::getVectorInUnitSphere().getMagnitude() <= 1.f);
	CHECK_APPROX(R::getVectorOnUnitSphere().getMagnitude(), 1.f);
	CHECK_APPROX(R::getVectorOnUnitSphericalCone(0.5f).getMagnitude(), 1.f);
	CHECK_APPROX(R::getVectorOnUnitSphericalCone(templates::Vector3<float>::UNIT_X, 0.5f).getMagnitude(), 1.f);
	CHECK_APPROX(R::getRotationMatrix().getDeterminant(), 1.f);
	CHECK_APPROX(R::getQuaternion().getMagnitude(), 1.f);

	int v = random(1, 3);
	CHECK((v >= 1 && v <= 3));
	double d = random(-1.0, 1.0);
	CHECK((d >= -1.0 && d <= 1.0));
}

} // TEST_SUITE

} // namespace mathematics
