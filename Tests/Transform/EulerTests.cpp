#include "TestUtilities.hpp"
#include <unordered_set>

namespace mathematics::templates {

namespace {

const EulerOrder ALL_ORDERS[] = { EulerOrder::XYZ, EulerOrder::XZY, EulerOrder::YZX, EulerOrder::YXZ, EulerOrder::ZXY, EulerOrder::ZYX };

// Reference rotation matrix: the three principal rotations applied in the order given by the name of the Euler order
// (row vectors, so the first rotation is the leftmost factor).
template<typename T>
Matrix3<T> composeRotations(const Euler<T>& e)
{
	auto [a1, a2, a3] = unpack(e.order);
	return Matrix3<T>::makeRotation(a1, e[int(a1)])*Matrix3<T>::makeRotation(a2, e[int(a2)])*
		Matrix3<T>::makeRotation(a3, e[int(a3)]);
}

template<typename T>
bool sameRotation(const Quaternion<T>& a, const Quaternion<T>& b)
{
	return a.approxEquals(b, test::tolerance<T>()) || a.approxEquals(-b, test::tolerance<T>());
}

} // namespace

TEST_SUITE("EulerOrder")
{

TEST_CASE("unpack")
{
	CHECK(unpack(EulerOrder::XYZ) == std::tuple(Axis::X, Axis::Y, Axis::Z));
	CHECK(unpack(EulerOrder::XZY) == std::tuple(Axis::X, Axis::Z, Axis::Y));
	CHECK(unpack(EulerOrder::YZX) == std::tuple(Axis::Y, Axis::Z, Axis::X));
	CHECK(unpack(EulerOrder::YXZ) == std::tuple(Axis::Y, Axis::X, Axis::Z));
	CHECK(unpack(EulerOrder::ZXY) == std::tuple(Axis::Z, Axis::X, Axis::Y));
	CHECK(unpack(EulerOrder::ZYX) == std::tuple(Axis::Z, Axis::Y, Axis::X));
}

TEST_CASE("pack")
{
	for (EulerOrder order : ALL_ORDERS)
	{
		auto [a1, a2, a3] = unpack(order);
		CHECK(pack(a1, a2) == order);
		CHECK(pack(a1, a2, a3) == order);
	}

	CHECK(pack(Axis::X, Axis::X) == EulerOrder::UNSPECIFIED);
	CHECK(pack(Axis::X, Axis::Y, Axis::X) == EulerOrder::UNSPECIFIED);
}

} // TEST_SUITE

TEST_SUITE("Euler")
{

TEST_CASE_TEMPLATE("construction and element access", T, TEST_REAL_TYPES)
{
	using E = Euler<T>;
	CHECK(E().isZero());
	CHECK(E().order == EulerOrder::UNSPECIFIED);
	CHECK(E(EulerOrder::XZY) == E::ZERO_XZY);

	E e(T(1), T(2), T(3), EulerOrder::YXZ);
	CHECK(e.x == T(1));
	CHECK(e.y == T(2));
	CHECK(e.z == T(3));
	CHECK(e.order == EulerOrder::YXZ);
	CHECK(E(e, EulerOrder::XYZ) == E(T(1), T(2), T(3), EulerOrder::XYZ));
	CHECK(E(std::tuple<T, T, T, EulerOrder>(T(1), T(2), T(3), EulerOrder::YXZ)) == e);
	CHECK(E(std::tuple<T, T, T>(T(1), T(2), T(3)), EulerOrder::YXZ) == e);
	const T a[] = { T(1), T(2), T(3) };
	CHECK(E(a, EulerOrder::YXZ) == e);

	E r = YawPitchRoll<T>(T(1), T(2), T(3));
	CHECK(r == E(T(2), T(1), T(3), EulerOrder::ZXY));

	CHECK(e[1] == T(2));
	e[1] = T(5);
	CHECK(e.y == T(5));
	CHECK(static_cast<const T*>(e)[2] == T(3));
	CHECK(e.template get<0>() == T(1));
	CHECK(e.template get<3>() == EulerOrder::YXZ);
	CHECK(get<2>(e) == T(3));
	auto [x, y, z, order] = e;
	CHECK(x + y + z == T(9));
	CHECK(order == EulerOrder::YXZ);
}

TEST_CASE_TEMPLATE("arithmetic", T, TEST_REAL_TYPES)
{
	using E = Euler<T>;
	E a(T(1), T(2), T(3), EulerOrder::XYZ);
	E b(T(0.5), T(-1), T(2), EulerOrder::XYZ);

	CHECK(+a == a);
	CHECK(-a == E(T(-1), T(-2), T(-3), EulerOrder::XYZ));
	CHECK(a + b == E(T(1.5), T(1), T(5), EulerOrder::XYZ));
	CHECK(a - b == E(T(0.5), T(3), T(1), EulerOrder::XYZ));
	CHECK(a*T(2) == E(T(2), T(4), T(6), EulerOrder::XYZ));
	CHECK(T(2)*a == E(T(2), T(4), T(6), EulerOrder::XYZ));
	CHECK(a/T(2) == E(T(0.5), T(1), T(1.5), EulerOrder::XYZ));
	CHECK(T(6)/a == E(T(6), T(3), T(2), EulerOrder::XYZ));

	E c = a;
	c += b;
	c -= b;
	CHECK(c == a);
	c *= T(2);
	c /= T(2);
	CHECK(c == a);

	// Mixing orders is an error.
	E d(T(1), T(2), T(3), EulerOrder::ZYX);
	CHECK_THROWS_AS(c += d, std::invalid_argument);
	CHECK_THROWS_AS(c -= d, std::invalid_argument);
	CHECK(a != d); // same angles, different order
}

TEST_CASE_TEMPLATE("predicates and modifiers", T, TEST_REAL_TYPES)
{
	using E = Euler<T>;
	E e(T(1), T(2), T(3), EulerOrder::XYZ);
	CHECK(E(Constants<T>::TOLERANCE*T(0.5), T(0), T(0), EulerOrder::XYZ).isApproxZero());
	CHECK(e.approxEquals(E(T(1), T(2), T(3) + Constants<T>::TOLERANCE*T(0.5), EulerOrder::XYZ)));
	CHECK(e.approxEquals(E(T(1), T(2), T(3.05), EulerOrder::XYZ), T(0.1)));
	CHECK_FALSE(e.approxEquals(E(T(1), T(2), T(3.05), EulerOrder::XYZ)));
	CHECK(e.isFinite());
	CHECK_FALSE(E(T(0), std::numeric_limits<T>::infinity(), T(0), EulerOrder::XYZ).isFinite());

	CHECK(E(e).setZero() == E::ZERO_XYZ);
	CHECK(E(e).setZero(EulerOrder::ZYX) == E::ZERO_ZYX);
	CHECK(E(e).set(T(4), T(5), T(6)) == E(T(4), T(5), T(6), EulerOrder::XYZ));
	CHECK(E(e).set(T(4), T(5), T(6), EulerOrder::YZX) == E(T(4), T(5), T(6), EulerOrder::YZX));
	CHECK(E(e).negate() == -e);
}

TEST_CASE_TEMPLATE("rotation matrices", T, TEST_REAL_TYPES)
{
	using E = Euler<T>;
	using M = Matrix3<T>;

	for (EulerOrder order : ALL_ORDERS)
	{
		CAPTURE(int(order));
		E e(T(0.3), T(-0.7), T(1.1), order);
		M expected = composeRotations(e);
		CHECK_APPROX(M(e), expected);
		CHECK_APPROX(M::makeRotation(e), expected);
		CHECK_APPROX(Matrix4<T>::makeRotation(e), Matrix4<T>(expected));
		CHECK_APPROX(AffineTransform<T>::makeRotation(e).getBasis(), expected);

		// Round trip through the matrix representation
		CHECK_APPROX(E(expected, order), e);
	}

	CHECK(M(E::ZERO_XYZ).isIdentity());
	CHECK(M(E(T(1), T(2), T(3), EulerOrder::UNSPECIFIED)).isIdentity());
	CHECK(E(M::IDENTITY, EulerOrder::UNSPECIFIED) == E::ZERO);
}

TEST_CASE_TEMPLATE("gimbal lock", T, TEST_REAL_TYPES)
{
	using E = Euler<T>;
	using M = Matrix3<T>;
	for (EulerOrder order : ALL_ORDERS)
	{
		CAPTURE(int(order));
		E e(T(0), T(0), T(0), order);
		e[int(std::get<1>(unpack(order)))] = Constants<T>::HALF_PI; // middle rotation of 90 degrees
		e[int(std::get<0>(unpack(order)))] = T(0.4);
		M m(e);
		CHECK_APPROX(M(E(m, order)), m); // the angles are not unique but must represent the same rotation
	}
}

TEST_CASE_TEMPLATE("quaternions", T, TEST_REAL_TYPES)
{
	using E = Euler<T>;
	using Q = Quaternion<T>;
	using M = Matrix3<T>;

	for (EulerOrder order : ALL_ORDERS)
	{
		CAPTURE(int(order));
		E e(T(-0.5), T(0.2), T(0.9), order);
		Q q(e);
		CHECK(sameRotation(q, Q(M(e))));
		CHECK_APPROX(E(q, order), e);
		CHECK_APPROX(Matrix4<T>::makeRotation(q), Matrix4<T>(M(e)));
	}

	CHECK(Q(E::ZERO_XYZ).isIdentity());
	CHECK(E(Q::IDENTITY, EulerOrder::YZX) == E::ZERO_YZX);
}

TEST_CASE_TEMPLATE("streaming and hashing", T, TEST_REAL_TYPES)
{
	using E = Euler<T>;
	for (EulerOrder order : ALL_ORDERS)
	{
		E e(T(1.5), T(-2), T(3), order);
		std::stringstream s;
		s << e;
		E f;
		s >> f;
		CHECK(f == e);
	}

	std::stringstream s;
	s << E(T(1), T(2), T(3), EulerOrder::UNSPECIFIED);
	CHECK(s.str() == "1 2 3 Unspecified");
	E u(T(0), T(0), T(0), EulerOrder::XYZ);
	s >> u;
	CHECK(u.order == EulerOrder::UNSPECIFIED);

	std::stringstream n("1 2 3 17");
	n >> u;
	CHECK(u.order == EulerOrder::ZXY); // numeric order values are accepted too

	std::unordered_set<E> set{ E(T(1), T(2), T(3), EulerOrder::XYZ), E(T(1), T(2), T(3), EulerOrder::XYZ),
		E(T(1), T(2), T(3), EulerOrder::ZYX) };
	CHECK(set.size() == 2);
}

} // TEST_SUITE

TEST_SUITE("YawPitchRoll")
{

TEST_CASE_TEMPLATE("construction, element access and arithmetic", T, TEST_REAL_TYPES)
{
	using R = YawPitchRoll<T>;
	CHECK(R().isZero());
	R r(T(1), T(2), T(3));
	CHECK(r.yaw == T(1));
	CHECK(r.pitch == T(2));
	CHECK(r.roll == T(3));
	CHECK(R(std::tuple<T, T, T>(T(1), T(2), T(3))) == r);
	CHECK(r.template get<0>() == T(1));
	CHECK(get<2>(r) == T(3));
	auto [yaw, pitch, roll] = r;
	CHECK(yaw + pitch + roll == T(6));

	R s(T(0.5), T(-1), T(2));
	CHECK(+r == r);
	CHECK(-r == R(T(-1), T(-2), T(-3)));
	CHECK(r + s == R(T(1.5), T(1), T(5)));
	CHECK(r - s == R(T(0.5), T(3), T(1)));
	CHECK(r*T(2) == R(T(2), T(4), T(6)));
	CHECK(T(2)*r == R(T(2), T(4), T(6)));
	CHECK(r/T(2) == R(T(0.5), T(1), T(1.5)));
	R c = r;
	c += s;
	c -= s;
	c *= T(2);
	c /= T(2);
	CHECK(c == r);
	CHECK(r != s);

	CHECK(R(Constants<T>::TOLERANCE*T(0.5), T(0), T(0)).isApproxZero());
	CHECK(r.approxEquals(R(T(1), T(2), T(3) + Constants<T>::TOLERANCE*T(0.5))));
	CHECK(r.approxEquals(R(T(1), T(2), T(3.05)), T(0.1)));
	CHECK_FALSE(r.approxEquals(R(T(1), T(2), T(3.05))));
	CHECK(r.isFinite());
	CHECK_FALSE(R(T(0), std::numeric_limits<T>::quiet_NaN(), T(0)).isFinite());
	CHECK(R(r).setZero().isZero());
	CHECK(R().set(T(1), T(2), T(3)) == r);
	CHECK(R(r).negate() == -r);
}

TEST_CASE_TEMPLATE("rotation conventions", T, TEST_REAL_TYPES)
{
	using R = YawPitchRoll<T>;
	using M = Matrix3<T>;
	using V = Vector3<T>;
	R r(T(0.4), T(-0.3), T(1.2));

	// Roll about Z is applied first, then pitch about X and yaw about Y.
	M expected = M::makeRotation(Axis::Z, r.roll)*M::makeRotation(Axis::X, r.pitch)*M::makeRotation(Axis::Y, r.yaw);
	CHECK_APPROX(M(r), expected);
	CHECK_APPROX(M::makeRotation(r), expected);
	CHECK_APPROX(M(Euler<T>(r)), expected); // Euler ZXY is the same rotation
	CHECK_APPROX(Matrix4<T>::makeRotation(r), Matrix4<T>(expected));
	CHECK_APPROX(Matrix4<T>::makeRotationTranslation(r, V(T(1), T(2), T(3))), Matrix4<T>(AffineTransform<T>(expected, V(T(1), T(2), T(3)))));
	CHECK_APPROX(AffineTransform<T>::makeRotation(r).getBasis(), expected);
	CHECK(sameRotation(Quaternion<T>(r), Quaternion<T>(expected)));
	CHECK(M(R()).isIdentity());
	CHECK(Quaternion<T>(R()).isIdentity());

	// The forward direction is the transformed Z axis.
	CHECK_APPROX(r.toForward(), V::UNIT_Z*expected);
}

TEST_CASE_TEMPLATE("conversions", T, TEST_REAL_TYPES)
{
	using R = YawPitchRoll<T>;
	using M = Matrix3<T>;
	R r(T(0.4), T(-0.3), T(1.2));

	CHECK_APPROX(R(M(r)), r);
	CHECK_APPROX(R(Quaternion<T>(r)), r);
	CHECK_APPROX(R(Euler<T>(r)), r);
	CHECK(R(Quaternion<T>::IDENTITY).isZero());
	CHECK(R(Euler<T>()).isZero());

	// From another Euler order representing the same rotation
	Euler<T> e(M(r), EulerOrder::XYZ);
	CHECK_APPROX(R(e), r);
}

TEST_CASE_TEMPLATE("forward direction", T, TEST_REAL_TYPES)
{
	using R = YawPitchRoll<T>;
	using V = Vector3<T>;
	CHECK_APPROX(R().toForward(), V::UNIT_Z);
	CHECK_APPROX(R(Constants<T>::HALF_PI, T(0), T(0)).toForward(), V::UNIT_X);
	CHECK_APPROX(R(T(0), Constants<T>::HALF_PI, T(0)).toForward(), V::MINUS_UNIT_Y);

	for (V forward : { V(T(1), T(2), T(3)), V(T(-1), T(0.5), T(-2)), V(T(-3), T(-1), T(1)) })
	{
		CAPTURE(forward);
		R r = R::fromForward(forward, T(0.25));
		CHECK_APPROX(r.toForward(), normalize(forward));
		CHECK(r.roll == T(0.25));
	}

	CHECK_APPROX(R::fromForward(V::UNIT_Y).pitch, -Constants<T>::HALF_PI);
	CHECK(R::fromForward(V::ZERO).isZero());
}

TEST_CASE_TEMPLATE("streaming and hashing", T, TEST_REAL_TYPES)
{
	using R = YawPitchRoll<T>;
	std::stringstream s;
	s << R(T(1.5), T(-2), T(3));
	CHECK(s.str() == "1.5 -2 3");
	R r;
	s >> r;
	CHECK(r == R(T(1.5), T(-2), T(3)));

	std::unordered_set<R> set{ r, r, R() };
	CHECK(set.size() == 2);
}

} // TEST_SUITE

} // namespace mathematics::templates
