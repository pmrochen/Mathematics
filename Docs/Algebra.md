# Algebra

[← Back to the README](../README.md)

- [Scalar functions](#scalar-functions)
- [Constants](#constants)
- [Vectors](#vectors)
- [Matrices](#matrices)
- [Quaternion](#quaternion)
- [LU decomposition](#lu-decomposition)

Read [Conventions](Conventions.md) first for the row vector convention, the order of composition and the direction of
rotations. Transformation-specific types (`AffineTransform`, `Matrix4`, Euler angles) are covered in
[Transformations](Transformations.md).

## Scalar functions

Header: `<Mathematics/Scalar.hpp>`, namespace `mathematics`.

| Function | Description |
| --- | --- |
| `abs(x)` | Absolute value (floating point and integers) |
| `sign(x)` | `-1`, `0` or `1` |
| `square(x)` | `x*x` |
| `frac(x)` | Fractional part, `x - floor(x)` |
| `approxEquals(a, b)`, `approxEquals(a, b, tolerance)` | `|b - a| < tolerance`; the default is `Constants<T>::TOLERANCE` |
| `log(x, b)` | Logarithm of `x` to the base `b` |
| `log2(x)` | Base-2 logarithm; for integers the position of the highest set bit (`log2(1024) == 10`) |
| `radians(x)`, `degrees(x)` | Degree and radian conversion |
| `lerp(a, b, t)` | Linear interpolation, `a + t*(b - a)` |
| `step(a, t)` | `1` when `t >= a`, otherwise `0` |
| `pulse(a, b, t)` | `1` when `a <= t < b`, otherwise `0` |
| `boxStep(a, b, t)` | Linear ramp from 0 at `a` to 1 at `b`, clamped |
| `smoothStep(a, b, t)` | Cubic Hermite ramp `3t² - 2t³`, clamped |
| `smootherStep(a, b, t)` | Quintic ramp `6t⁵ - 15t⁴ + 10t³`, clamped |

```cpp
float angle = radians(45.f);
float fade = smoothStep(0.f, 1.f, 0.25f);       // 0.15625
```

## Constants

Header: `<Mathematics/Constants.hpp>`. `Constants<T>` is specialized for `float` and `double`.

| Constant | Value |
| --- | --- |
| `TOLERANCE` | Default comparison tolerance: `1e-6f` / `1e-15` |
| `EPSILON` | `std::numeric_limits<T>::epsilon()` |
| `INF`, `MINUS_INF` | Positive and negative infinity |
| `E`, `LN2` | Euler's number, natural logarithm of 2 |
| `PI`, `HALF_PI`, `TWO_PI`, `INV_PI`, `PI_SQR` | π, π/2, 2π, 1/π, π² |
| `SQRT_2`, `SQRT_3` | Square roots of 2 and 3 |
| `DEG_TO_RAD`, `RAD_TO_DEG` | Conversion factors |

```cpp
float circumference = Constants<float>::TWO_PI*radius;
```

## Vectors

Headers: `<Mathematics/Algebra/Vector2.hpp>`, `Vector3.hpp`, `Vector4.hpp`.

| Type | Components | Aliases |
| --- | --- | --- |
| `templates::Vector2<T>` | `x, y` | `Vector2`, `Vector2f`, `Vector2d`, `IntVector2` |
| `templates::Vector3<T>` | `x, y, z` | `Vector3`, `Vector3f`, `Vector3d`, `IntVector3` |
| `templates::Vector4<T>` | `x, y, z, w` | `Vector4`, `Vector4f`, `Vector4d`, `IntVector4` |

`T` may be a floating point or an integral type. Integer vectors offer construction, arithmetic, comparisons, `min`,
`max`, `abs` and component access; everything that involves lengths, angles or matrices requires a floating point
type.

### Construction

```cpp
Vector3 zero;                           // (0, 0, 0)
Vector3 ones(1.f);                      // (1, 1, 1); explicit
Vector3 v(1.f, 2.f, 3.f);
Vector3 fromXy(Vector2(1.f, 2.f), 3.f); // extend a smaller vector
Vector4 point(v, 1.f);
Vector3 up(Axis::Y);                    // unit vector along an axis; explicit
Vector3d precise(v);                    // convert the component type; explicit
Vector3 scratch(UNINITIALIZED);         // uninitialized
```

Vectors also construct from a `const T*`, from a `std::tuple` and, for `Vector2`, from a `std::pair`.

### Constants

`ZERO`, `ONE`, `UNIT_X`, `UNIT_Y`, `UNIT_Z`, `UNIT_W`, `MINUS_UNIT_X`, `MINUS_UNIT_Y`, `MINUS_UNIT_Z`, `MINUS_UNIT_W`
(as far as the dimension allows), `TOLERANCE`, `INF` and `MINUS_INF`. Integer vectors provide `ZERO`.

### Operators

| Operator | Meaning |
| --- | --- |
| `a + b`, `a - b`, `-a` | Component-wise sum, difference, negation |
| `a*b`, `a/b` | **Component-wise** product and quotient |
| `a*s`, `s*a`, `a/s`, `s/a` | Scalar product and quotient |
| `v*m` | Transforms the row vector by a matrix of the same dimension |
| `m*v` | Transforms `v` as a column vector, i.e. by the transposed matrix |
| `a == b`, `a != b` | Exact comparison |
| `v[i]` | Component access without bounds checking |

All arithmetic operators have compound assignment forms.

### Member functions

| Member | Description |
| --- | --- |
| `set(x, y, ...)`, `setZero()` | Assign components |
| `getMagnitude()`, `getMagnitudeSquared()` | Length; `getLength()` and `getLengthSquared()` are synonyms |
| `setMagnitude(m)`, `setLength(m)` | Rescale to the given length |
| `normalize()` | Scale to unit length; a zero vector is left unchanged |
| `negate()` | Negate in place |
| `getMinComponent()`, `getMaxComponent()` | Smallest and largest component |
| `getMajorAxis()` | Axis of the component with the largest magnitude (`Vector2`, `Vector3`) |
| `setMinimum(a, b)`, `setMaximum(a, b)` | Component-wise minimum or maximum of two vectors |
| `isZero()`, `isApproxZero()`, `isFinite()` | Predicates |
| `approxEquals(v)`, `approxEquals(v, tolerance)` | Approximate comparison |
| `allLessThan(v)`, `allLessThanEqual(v)`, `allGreaterThan(v)`, `allGreaterThanEqual(v)` | True when the relation holds for every component |
| `anyLessThan(v)`, `anyLessThanEqual(v)`, `anyGreaterThan(v)`, `anyGreaterThanEqual(v)` | True when the relation holds for at least one component |
| `xy()`, `xz()`, `zy()` (`Vector3`); `xy()`, `xyz()` (`Vector4`) | Sub-vectors |
| `transform(m)` | Same as `*this *= m` |
| `rotate(angle)` (`Vector2`) | Rotate counter-clockwise (from +X toward +Y) |
| `rotate(axis, angle)`, `rotate(q)` (`Vector3`) | Rotate about a principal axis or by a quaternion |

### Free functions

| Function | `Vector2` | `Vector3` | `Vector4` | Description |
| --- | :-: | :-: | :-: | --- |
| `dot(a, b)` | ✓ | ✓ | ✓ | Dot product |
| `cross(a, b)` | ✓ | ✓ | | 3D: cross product. 2D: the scalar `a.x*b.y - a.y*b.x` |
| `cross(a, b, c)` | | | ✓ | 4D generalization of the cross product |
| `length(v)`, `lengthSquared(v)` | ✓ | ✓ | ✓ | Length |
| `distance(a, b)`, `distanceSquared(a, b)` | ✓ | ✓ | ✓ | Distance between two points |
| `normalize(v)` | ✓ | ✓ | ✓ | Unit-length copy |
| `abs(v)`, `sum(v)` | ✓ | ✓ | ✓ | Component-wise absolute value; sum of components |
| `min(a, b)`, `max(a, b)`, `clamp(v, low, high)` | ✓ | ✓ | ✓ | Component-wise |
| `lerp(a, b, t)` | ✓ | ✓ | ✓ | Linear interpolation |
| `slerp(a, b, t)` | ✓ | ✓ | ✓ | Spherical interpolation between unit vectors |
| `project(a, b)` | ✓ | ✓ | ✓ | Projection of `a` onto `b` |
| `angle(a, b)` | ✓ | ✓ | | Angle between two vectors, in radians |
| `perpendicular(v)` | ✓ | ✓ | | A vector perpendicular to `v`; in 2D, `v` rotated by +90° |
| `reflect(i, n)` | | ✓ | | Reflection of the incident vector `i` about the unit normal `n` |
| `refract(i, n, etaRatio)` | | ✓ | | Refraction; zero vector on total internal reflection |
| `rotate(v, q)` | | ✓ | | Rotation by a quaternion |
| `transform(v, m)` | ✓ | ✓ | ✓ | Same as `v*m`; `Vector3` also accepts an `AffineTransform` |
| `tensor(a, b)` | | ✓ | ✓ | Outer product, a matrix with rows `a[i]*b` |

```cpp
Vector3 a(1.f, 2.f, 3.f);

dot(a, a);                                              // 14
normalize(a);                                           // (0.267261, 0.534522, 0.801784)
project(a, Vector3::UNIT_X);                            // (1, 0, 0)
reflect(Vector3(1.f, -1.f, 0.f), Vector3::UNIT_Y);      // (1, 1, 0)
angle(Vector3::UNIT_X, Vector3::UNIT_Y);                // 1.5708
clamp(Vector3(5.f, -5.f, 0.5f), Vector3::ZERO, Vector3::ONE);   // (1, 0, 0.5)
```

The concept `templates::AnyVector<V, T>` is satisfied by `Vector2<T>`, `Vector3<T>` and `Vector4<T>`.

## Matrices

Headers: `<Mathematics/Algebra/Matrix2.hpp>`, `Matrix3.hpp`, `Matrix4.hpp`.

| Type | Size | Typical use | Aliases |
| --- | --- | --- | --- |
| `templates::Matrix2<T>` | 2×2 | 2D rotation and scale | `Matrix2`, `Matrix2f`, `Matrix2d` |
| `templates::Matrix3<T>` | 3×3 | 3D rotation and scale, orientation bases | `Matrix3`, `Matrix3f`, `Matrix3d` |
| `templates::Matrix4<T>` | 4×4 | General homogeneous transformations | `Matrix4`, `Matrix4f`, `Matrix4d` |

Matrices require a floating point component type. The elements are public members named `mRC` (`m00`, `m01`, ...);
`DIMENSION` gives the order.

### Construction

```cpp
Matrix3 zero;                                               // all zeros
Matrix3 identity{ Identity() };                             // or Matrix3::IDENTITY
Matrix3 elements(1.f, 0.f, 0.f,
                 0.f, 1.f, 0.f,
                 0.f, 0.f, 1.f);                            // row by row
Matrix3 rows(Vector3::UNIT_X, Vector3::UNIT_Y, Vector3::UNIT_Z);
Matrix3 fromQuaternion(Quaternion::fromAxisAngle(Vector3::UNIT_Y, radians(90.f)));
Matrix3 embedded(Matrix2::makeRotation(radians(30.f)));     // upper-left block, rest identity
```

`Matrix3` can also be built from an orientation:

| Constructor | Result |
| --- | --- |
| `Matrix3(forward)` | A basis whose Z axis is `forward`, with the X axis kept horizontal |
| `Matrix3(up, forward)` | A basis whose Z axis is `forward` and whose Y axis is as close to `up` as possible |
| `Matrix3(yawPitchRoll)`, `Matrix3(euler)`, `Matrix3(quaternion)` | Conversion from another rotation representation |

Neither argument needs to be normalized. Degenerate input (a zero vector) gives the identity.

### Factory functions

Each `makeX` has a matching in-place `setX`.

| Function | `Matrix2` | `Matrix3` | Description |
| --- | :-: | :-: | --- |
| `makeScaling(v)`, `makeScaling(f)` | ✓ | ✓ | Non-uniform or uniform scale |
| `makeRotation(angle)` | ✓ | | Rotation in the plane |
| `makeRotation(Axis, angle)` | | ✓ | Rotation about a principal axis |
| `makeRotation(axis, angle)` | | ✓ | Rotation about an arbitrary axis (normalized internally) |
| `makeRotation(yawPitchRoll)`, `makeRotation(euler)`, `makeRotation(quaternion)` | | ✓ | Conversion |
| `makeScalingRotation(s, ...)` | ✓ | ✓ | Scale followed by any of the rotations above |
| `makeShearing(xy, yx)` | ✓ | | Shear; the factor `ab` adds `ab` times the `a` coordinate to the `b` coordinate |
| `makeShearing(xy, xz, yx, yz, zx, zy)` | | ✓ | Shear |

`Matrix4` and `AffineTransform` add translation and reflection; see [Transformations](Transformations.md).

### Member functions

| Member | Description |
| --- | --- |
| `m[i]`, `getRow(i)` | Row `i` as a vector |
| `getColumn(i)` | Column `i` as a vector (by value) |
| `set(...)`, `setZero()`, `setIdentity()` | Assign |
| `getTrace()`, `getDeterminant()` | Scalar properties |
| `isSingular()` | True when the determinant is zero |
| `getScaleComponent()` | Lengths of the rows, i.e. the scale along each transformed axis (`Matrix2`, `Matrix3`) |
| `isZero()`, `isIdentity()`, `isFinite()` | Exact predicates |
| `isApproxZero()`, `isApproxIdentity()`, `approxEquals(m)`, `approxEquals(m, tolerance)` | Approximate predicates |
| `isApproxOrthogonal()` | True when the rows are approximately orthonormal (`Matrix2`, `Matrix3`) |
| `hasApproxUniformScaling()` | True when all rows have approximately the same length (`Matrix2`, `Matrix3`) |
| `transpose()`, `invert()`, `negate()` | In-place operations |
| `setTranspose(m)`, `setInverse(m)`, `setInverseTranspose(m)` | Assign from an operation on another matrix |
| `orthonormalize()` | Make the rows orthonormal (`Matrix2`, `Matrix3`) |
| `concatenate(m)`, `preConcatenate(m)` | Append or prepend a transformation |
| `scale(v)`, `preScale(v)` | Append or prepend a scale; `Matrix2` and `Matrix3` also accept a uniform factor |

### Operators and free functions

| Function | Description |
| --- | --- |
| `a + b`, `a - b`, `-a`, `a*s`, `s*a`, `a/s` | Element-wise arithmetic |
| `a*b` | Matrix product; applies `a`, then `b` |
| `concatenate(a, b)`, `concatenate(a, b, c)`, `concatenate(a, b, c, d)` | Product of two to four matrices, applied left to right |
| `transpose(m)`, `inverse(m)`, `inverseTranspose(m)` | Copies |
| `adjoint(m)` | Adjugate matrix; `inverse(m) == adjoint(m)/determinant` |

```cpp
Matrix3 rotation = Matrix3::makeRotation(Axis::Z, radians(30.f));

rotation[0];                                        // (0.866025, 0.5, 0): where the X axis ends up
rotation.isApproxOrthogonal();                      // true
inverse(rotation).approxEquals(transpose(rotation));// true for any pure rotation

Matrix3 combined = concatenate(Matrix3::makeScaling(2.f), rotation);    // scale, then rotate
Vector3 v = Vector3::UNIT_X*combined;               // (1.73205, 1, 0)
```

The concept `templates::AnyMatrix<M, T>` is satisfied by `Matrix2<T>`, `Matrix3<T>` and `Matrix4<T>`.

## Quaternion

Header: `<Mathematics/Algebra/Quaternion.hpp>`. Alias: `Quaternion`.

A quaternion has the vector part `x, y, z` and the scalar part `w`. Unit quaternions represent rotations.

### Construction

```cpp
Quaternion zero;                                                    // (0, 0, 0, 0)
Quaternion identity{ Identity() };                                  // (0, 0, 0, 1); or Quaternion::IDENTITY
Quaternion q = Quaternion::fromAxisAngle(Vector3::UNIT_Z, radians(90.f));
Quaternion fromMatrix(Matrix3::makeRotation(Axis::Y, 1.f));         // from a rotation matrix
Quaternion fromAngles(YawPitchRoll(radians(90.f), 0.f, 0.f));       // also from Euler
Quaternion raw(0.f, 0.f, 0.707107f, 0.707107f);                     // x, y, z, w
Quaternion parts(Vector3(0.f, 0.f, 0.707107f), 0.707107f);          // vector, scalar
```

The axis passed to `fromAxisAngle` does not need to be normalized.

### Member functions

| Member | Description |
| --- | --- |
| `getVector()`, `setVector(v)`, `getScalar()`, `setScalar(s)` | The two parts |
| `getAxis()`, `getAngle()` | Rotation axis and angle of a **unit** quaternion |
| `getMagnitude()` (synonym `getAbsoluteValue()`), `getMagnitudeSquared()` (synonym `getNorm()`) | Length and squared length |
| `set(x, y, z, w)`, `set(vector, scalar)`, `setZero()`, `setIdentity()` | Assign |
| `normalize()`, `negate()` | In place |
| `conjugate()`, `setConjugate(q)` | Conjugate; the inverse of a unit quaternion |
| `invert()`, `setInverse(q)` | General inverse |
| `concatenate(q)`, `preConcatenate(q)` | Append or prepend a rotation |
| `rotate(q)` | Conjugate this quaternion by `q`, i.e. rotate the rotation |
| `isZero()`, `isIdentity()`, `isFinite()`, `isApproxZero()`, `isApproxIdentity()`, `approxEquals(...)` | Predicates |

### Operators and free functions

| Function | Description |
| --- | --- |
| `q1 + q2`, `q1 - q2`, `-q`, `q*s`, `s*q` | Component-wise arithmetic |
| `q1*q2` | Hamilton product; **rotates by `q2` first, then by `q1`** |
| `concatenate(q1, q2)`, up to four arguments | Rotates by `q1`, then `q2`, ... |
| `rotate(v, q)` | Rotates a `Vector3` |
| `rotate(q1, q2)` | `q1` conjugated by `q2` |
| `conjugate(q)`, `inverse(q)`, `normalize(q)` | Copies |
| `dot(q1, q2)` | Four-dimensional dot product |
| `lerp(q1, q2, t)` | Component-wise interpolation (not normalized) |
| `slerp(q1, q2, t)` | Spherical linear interpolation |
| `arc(v1, v2)` | The shortest rotation that takes direction `v1` to direction `v2` |

```cpp
Quaternion aboutZ = Quaternion::fromAxisAngle(Vector3::UNIT_Z, radians(90.f));
Quaternion aboutX = Quaternion::fromAxisAngle(Vector3::UNIT_X, radians(90.f));

rotate(Vector3::UNIT_X, aboutZ);                            // (0, 1, 0)
rotate(Vector3::UNIT_X, concatenate(aboutZ, aboutX));       // (0, 0, 1): about Z, then about X
rotate(Vector3::UNIT_X, aboutX*aboutZ);                     // (0, 0, 1): the same rotation

aboutZ.getAxis();                                           // (0, 0, 1)
degrees(aboutZ.getAngle());                                 // 90

slerp(Quaternion::IDENTITY, aboutZ, 0.5f);                  // 45 degrees about Z
rotate(Vector3::UNIT_X, arc(Vector3::UNIT_X, Vector3::UNIT_Y));     // (0, 1, 0)
```

## LU decomposition

Header: `<Mathematics/Algebra/LU.inl>` (included by the mesh headers and therefore by the umbrella header). Namespace
`mathematics::lu`.

LU decomposition with partial pivoting for square matrices of any fixed order. It accepts the library's matrix types
or a flat `std::array<T, Order*Order>` in row-major order.

| Function | Description |
| --- | --- |
| `lu::decompose(matrix, tolerance)` | Returns `std::optional<Decomposition<Order, T>>`; empty when the matrix is singular. The default tolerance only rejects matrices that are singular up to rounding |
| `lu::determinant(decomposition)` | Determinant |
| `lu::solve<V>(decomposition, b)` | Solves `M*x = b` with `x` and `b` as column vectors |
| `lu::inverse<M>(decomposition)` | Inverse matrix |

```cpp
Matrix3 m(1.f, 2.f, 0.f,
          0.f, 1.f, 0.f,
          0.f, 0.f, 1.f);

if (auto decomposition = lu::decompose(m))
{
    float determinant = lu::determinant(*decomposition);                                // 1
    Vector3 x = lu::solve<Vector3>(*decomposition, Vector3(3.f, 1.f, 0.f));             // (1, 1, 0)
    Matrix3 inverted = lu::inverse<Matrix3>(*decomposition);
}
```

For 2×2 to 4×4 matrices the member `invert()` and the free `inverse()` are the simpler choice; `lu` is useful for
solving systems, for testing singularity with a tolerance and for larger orders.
