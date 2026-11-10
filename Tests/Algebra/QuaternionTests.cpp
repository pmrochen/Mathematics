#include "TestUtilities.hpp"
#include <unordered_set>

namespace mathematics::templates {

namespace {

// q and -q represent the same rotation
template<typename T>
bool sameRotation(const Quaternion<T>& a, const Quaternion<T>& b)
{
	return a.approxEquals(b, test::tolerance<T>()) || a.approxEquals(-b, test::tolerance<T>());
}

} // namespace

TEST_SUITE("Quaternion")
{

TEST_CASE_TEMPLATE("construction and element access", T, TEST_REAL_TYPES)
{
	using Q = Quaternion<T>;
	using V = Vector3<T>;

	CHECK(Q().isZero());
	CHECK(Q(Identity()).isIdentity());
	CHECK(Q(T(2)) == Q(T(0), T(0), T(0), T(2)));

	Q q(T(1), T(2), T(3), T(4));
	CHECK(q.x == T(1));
	CHECK(q.y == T(2));
	CHECK(q.z == T(3));
	CHECK(q.w == T(4));
	CHECK(Q(V(T(1), T(2), T(3))) == Q(T(1), T(2), T(3), T(0)));
	CHECK(Q(V(T(1), T(2), T(3)), T(4)) == q);
	CHECK(Q(std::tuple<T, T, T, T>(T(1), T(2), T(3), T(4))) == q);
	CHECK(Q(std::tuple<int, int, int, int>(1, 2, 3, 4)) == q);
	const T a[] = { T(1), T(2), T(3), T(4) };
	CHECK(Q(a) == q);

	std::tuple<T, T, T, T> t = q;
	CHECK(t == std::tuple<T, T, T, T>(T(1), T(2), T(3), T(4)));
	CHECK(q[3] == T(4));
	CHECK(static_cast<const T*>(q)[2] == T(3));
	CHECK(q.template get<1>() == T(2));
	CHECK(get<3>(q) == T(4));
	auto [x, y, z, w] = q;
	CHECK(x + y + z + w == T(10));

	CHECK(q.getVector() == V(T(1), T(2), T(3)));
	CHECK(q.getScalar() == T(4));
	q.setVector(V(T(5), T(6), T(7)));
	q.setScalar(T(8));
	CHECK(q == Q(T(5), T(6), T(7), T(8)));
}

TEST_CASE_TEMPLATE("arithmetic", T, TEST_REAL_TYPES)
{
	using Q = Quaternion<T>;
	using V = Vector3<T>;
	Q a(T(1), T(2), T(3), T(4));
	Q b(T(-1), T(0.5), T(2), T(1));

	CHECK(+a == a);
	CHECK(-a == Q(T(-1), T(-2), T(-3), T(-4)));
	CHECK(a + b == Q(T(0), T(2.5), T(5), T(5)));
	CHECK(a - b == Q(T(2), T(1.5), T(1), T(3)));
	CHECK(a*T(2) == Q(T(2), T(4), T(6), T(8)));
	CHECK(T(2)*a == Q(T(2), T(4), T(6), T(8)));
	CHECK(a/T(2) == Q(T(0.5), T(1), T(1.5), T(2)));

	// Hamilton product: i*j = k, j*k = i, k*i = j, i*i = -1
	Q i(T(1), T(0), T(0), T(0)), j(T(0), T(1), T(0), T(0)), k(T(0), T(0), T(1), T(0));
	CHECK(i*j == k);
	CHECK(j*k == i);
	CHECK(k*i == j);
	CHECK(j*i == -k);
	CHECK(i*i == Q(T(0), T(0), T(0), T(-1)));
	CHECK_APPROX(a*b, Q(T(-0.5), T(-1), T(13.5), T(-2)));

	// Products with pure vector quaternions
	V v(T(1), T(-2), T(0.5));
	CHECK_APPROX(a*v, a*Q(v));
	CHECK_APPROX(v*a, Q(v)*a);

	CHECK_APPROX(T(2)/a, inverse(a)*T(2));

	Q c = a;
	c += b;
	c -= b;
	CHECK(c == a);
	c *= T(2);
	c /= T(2);
	CHECK(c == a);
	c *= b;
	CHECK_APPROX(c, a*b);
	c = a;
	c *= v;
	CHECK_APPROX(c, a*v);
	CHECK(a != b);
}

TEST_CASE_TEMPLATE("predicates and magnitude", T, TEST_REAL_TYPES)
{
	using Q = Quaternion<T>;
	CHECK(Q::ZERO.isZero());
	CHECK(Q::IDENTITY.isIdentity());
	CHECK(Q(Constants<T>::TOLERANCE*T(0.5), T(0), T(0), T(0)).isApproxZero());
	CHECK(Q(T(0), T(0), Constants<T>::TOLERANCE*T(0.5), T(1)).isApproxIdentity());
	CHECK_FALSE(Q(T(0), T(0), T(0.1), T(1)).isApproxIdentity());
	CHECK(Q(T(1), T(2), T(3), T(4)).approxEquals(Q(T(1), T(2), T(3), T(4.05)), T(0.1)));
	CHECK_FALSE(Q(T(1), T(2), T(3), T(4)).approxEquals(Q(T(1), T(2), T(3), T(4.05))));
	CHECK(Q::IDENTITY.isFinite());
	CHECK_FALSE(Q(T(0), T(0), T(0), std::numeric_limits<T>::infinity()).isFinite());

	Q q(T(1), T(2), T(2), T(4));
	CHECK_APPROX(q.getMagnitude(), T(5));
	CHECK_APPROX(q.getAbsoluteValue(), T(5));
	CHECK(q.getNorm() == T(25));
	CHECK(q.getMagnitudeSquared() == T(25));
	CHECK_APPROX(normalize(q), q/T(5));
	CHECK_APPROX(normalize(Q(q) + Q::ZERO), q/T(5)); // rvalue overload
	CHECK_APPROX(Q(q).normalize().getMagnitude(), T(1));
	CHECK(Q().normalize().isZero());
	CHECK(dot(q, Q::IDENTITY) == T(4));
}

TEST_CASE_TEMPLATE("modifiers", T, TEST_REAL_TYPES)
{
	using Q = Quaternion<T>;
	using V = Vector3<T>;
	Q q(T(1), T(2), T(3), T(4));
	CHECK(Q(q).setZero().isZero());
	CHECK(Q(q).setIdentity().isIdentity());
	CHECK(Q().set(V(T(1), T(2), T(3)), T(4)) == q);
	CHECK(Q().set(T(1), T(2), T(3), T(4)) == q);
	CHECK(Q().setConjugate(q) == Q(T(-1), T(-2), T(-3), T(4)));
	CHECK(Q(q).conjugate() == Q(T(-1), T(-2), T(-3), T(4)));
	CHECK(conjugate(q) == Q(T(-1), T(-2), T(-3), T(4)));
	CHECK(conjugate(Q(q) + Q::ZERO) == Q(T(-1), T(-2), T(-3), T(4))); // rvalue overload
	CHECK(Q(q).negate() == -q);

	CHECK_APPROX(inverse(q)*q, Q::IDENTITY);
	CHECK_APPROX(q*inverse(q), Q::IDENTITY);
	CHECK_APPROX(Q(q).invert(), inverse(q));
	CHECK_APPROX(Q().setInverse(q), inverse(q));
}

TEST_CASE_TEMPLATE("axis-angle", T, TEST_REAL_TYPES)
{
	using Q = Quaternion<T>;
	using V = Vector3<T>;
	V axis = normalize(V(T(1), T(2), T(-2)));
	Q q = Q::fromAxisAngle(axis*T(3), T(1.2)); // the axis does not need to be normalized
	CHECK_APPROX(q.getMagnitude(), T(1));
	CHECK_APPROX(q.getAxis(), axis);
	CHECK_APPROX(q.getAngle(), T(1.2));
	CHECK_APPROX(q, Q(axis*std::sin(T(0.6)), std::cos(T(0.6))));

	CHECK(Q::fromAxisAngle(V::ZERO, T(1)).isIdentity());
	CHECK(Q::fromAxisAngle(V::UNIT_X, T(0)).isIdentity());
	CHECK(Q::IDENTITY.getAxis().isZero());
	CHECK(Q::IDENTITY.getAngle() == T(0));
}

TEST_CASE_TEMPLATE("rotation", T, TEST_REAL_TYPES)
{
	using Q = Quaternion<T>;
	using V = Vector3<T>;
	using M = Matrix3<T>;
	const T halfPi = Constants<T>::HALF_PI;

	CHECK_APPROX(rotate(V::UNIT_X, Q::fromAxisAngle(V::UNIT_Z, halfPi)), V::UNIT_Y);
	CHECK_APPROX(rotate(V::UNIT_Y, Q::fromAxisAngle(V::UNIT_X, halfPi)), V::UNIT_Z);
	CHECK_APPROX(rotate(V::UNIT_Z, Q::fromAxisAngle(V::UNIT_Y, halfPi)), V::UNIT_X);

	V axis = normalize(V(T(0.3), T(-1), T(0.6)));
	Q q = Q::fromAxisAngle(axis, T(0.9));
	V v(T(1.5), T(0.5), T(-2));
	V expected = v*M::makeRotation(axis, T(0.9));
	CHECK_APPROX(rotate(v, q), expected);
	CHECK_APPROX(V(v).rotate(q), expected);
	CHECK_APPROX((q*v*conjugate(q)).getVector(), expected);

	// Rotating a quaternion: q*p*q^-1
	Q p(v, T(0));
	CHECK_APPROX(Q(p).rotate(q).getVector(), expected);
}

TEST_CASE_TEMPLATE("conversion to and from matrices", T, TEST_REAL_TYPES)
{
	using Q = Quaternion<T>;
	using V = Vector3<T>;
	using M = Matrix3<T>;

	// Exercise all branches of the matrix-to-quaternion conversion (positive trace and each dominant diagonal element).
	const V axes[] = { V(T(1), T(2), T(3)), V::UNIT_X, V::UNIT_Y, V::UNIT_Z, V(T(-1), T(0.2), T(0.1)) };
	const T angles[] = { T(0.3), T(2.5), T(3), Constants<T>::PI };
	for (const V& axis : axes)
	{
		for (T angle : angles)
		{
			CAPTURE(axis);
			CAPTURE(angle);
			M m = M::makeRotation(axis, angle);
			Q q = Q::fromAxisAngle(axis, angle);
			CHECK(sameRotation(Q(m), q));
			CHECK_APPROX(M(q), m);
			CHECK_APPROX(M::makeRotation(q), m);
			CHECK_APPROX(M().setRotation(q), m);
		}
	}

	CHECK(Q(M::IDENTITY).isApproxIdentity());
	CHECK(M(Q::IDENTITY).isApproxIdentity());
}

TEST_CASE_TEMPLATE("concatenation", T, TEST_REAL_TYPES)
{
	using Q = Quaternion<T>;
	using V = Vector3<T>;
	using M = Matrix3<T>;
	Q a = Q::fromAxisAngle(V(T(1), T(0), T(1)), T(0.5));
	Q b = Q::fromAxisAngle(V(T(0), T(1), T(-1)), T(1.3));
	Q c = Q::fromAxisAngle(V(T(1), T(1), T(1)), T(-0.7));
	Q d = Q::fromAxisAngle(V::UNIT_Z, T(2.2));
	V v(T(1), T(2), T(3));

	// concatenate(a, b) applies a first, then b (consistent with a row vector matrix product).
	CHECK_APPROX(rotate(v, concatenate(a, b)), rotate(rotate(v, a), b));
	CHECK_APPROX(M(concatenate(a, b)), M(a)*M(b));
	CHECK_APPROX(rotate(v, concatenate(a, b, c)), rotate(rotate(rotate(v, a), b), c));
	CHECK_APPROX(rotate(v, concatenate(a, b, c, d)), rotate(rotate(rotate(rotate(v, a), b), c), d));
	CHECK_APPROX(Q(a).concatenate(b), concatenate(a, b));
	CHECK_APPROX(Q(a).preConcatenate(b), concatenate(b, a));
}

TEST_CASE_TEMPLATE("interpolation", T, TEST_REAL_TYPES)
{
	using Q = Quaternion<T>;
	using V = Vector3<T>;
	Q a = Q::fromAxisAngle(V::UNIT_Y, T(0.2));
	Q b = Q::fromAxisAngle(V::UNIT_Y, T(1.4));

	CHECK_APPROX(lerp(a, b, T(0)), a);
	CHECK_APPROX(lerp(a, b, T(1)), b);
	CHECK_APPROX(lerp(a, b, T(0.5)), (a + b)*T(0.5));

	CHECK(sameRotation(slerp(a, b, T(0)), a));
	CHECK(sameRotation(slerp(a, b, T(1)), b));
	CHECK(sameRotation(slerp(a, b, T(0.5)), Q::fromAxisAngle(V::UNIT_Y, T(0.8))));
	CHECK(sameRotation(slerp(a, b, T(0.25)), Q::fromAxisAngle(V::UNIT_Y, T(0.5))));
	CHECK_APPROX(slerp(a, b, T(0.3)).getMagnitude(), T(1));

	// Takes the shortest path even if the inputs lie in opposite hemispheres.
	CHECK(sameRotation(slerp(a, -b, T(0.5)), Q::fromAxisAngle(V::UNIT_Y, T(0.8))));
	// Nearly identical rotations
	CHECK(sameRotation(slerp(a, a, T(0.5)), a));
}

TEST_CASE_TEMPLATE("arc", T, TEST_REAL_TYPES)
{
	using Q = Quaternion<T>;
	using V = Vector3<T>;
	V u(T(1), T(2), T(-1));
	V v(T(-3), T(0.5), T(2));
	Q q = arc(u, v);
	CHECK_APPROX(q.getMagnitude(), T(1));
	CHECK_APPROX(rotate(normalize(u), q), normalize(v));
	CHECK(arc(u, u*T(2)).isApproxIdentity());
	CHECK_APPROX(rotate(V::UNIT_X, arc(V::UNIT_X, V::UNIT_Y)), V::UNIT_Y);
}

TEST_CASE_TEMPLATE("constants, streaming and hashing", T, TEST_REAL_TYPES)
{
	using Q = Quaternion<T>;
	CHECK(Q::ZERO == Q(T(0), T(0), T(0), T(0)));
	CHECK(Q::IDENTITY == Q(T(0), T(0), T(0), T(1)));

	std::stringstream s;
	s << Q(T(1.5), T(-2), T(3), T(4));
	CHECK(s.str() == "1.5 -2 3 4");
	Q q;
	s >> q;
	CHECK(q == Q(T(1.5), T(-2), T(3), T(4)));

	std::unordered_set<Q> set{ q, q, Q::IDENTITY };
	CHECK(set.size() == 2);
}

} // TEST_SUITE

} // namespace mathematics::templates
