# Utilities

[← Back to the README](../README.md)

- [Interval](#interval)
- [Random](#random)
- [Perlin noise](#perlin-noise)
- [Romberg integration](#romberg-integration)
- [Spline bases](#spline-bases)
- [Axis and concepts](#axis-and-concepts)

## Interval

Header: `<Mathematics/Interval.hpp>`. Alias: `Interval`.

A closed range `[minimum, maximum]` of scalars. It is the result type of the intersection queries that return an
entry and an exit parameter, and the depth range of a `SymmetricFrustum`.

```cpp
Interval range(0.f, 2.f);

range.contains(1.f);                                        // true
range.getLength();                                          // 2
range.interpolate(0.25f);                                   // 0.5
range.intersects(Interval(2.f, 3.f));                       // true: the end points touch

Interval::makeIntersection(range, Interval(1.f, 5.f));      // [1, 2]
Interval::makeUnion(range, Interval(5.f, 6.f));             // [0, 6]: the hull of both
```

A default constructed interval is **empty** (`minimum > maximum`), so it can accumulate values:

```cpp
Interval extent;
extent.extendBy(3.f).extendBy(-1.f);                        // [-1, 3]
```

| Member | Description |
| --- | --- |
| `Interval()`, `Interval(minimum, maximum)`, `Interval(value)` | Empty, range, single point |
| `isEmpty()`, `clear()` | Empty state; `EMPTY` is the constant |
| `isPoint()` | `minimum == maximum` |
| `isZero()`, `isApproxZero()`, `setZero()`, `approxEquals(...)`, `isFinite()` | Predicates; `ZERO` is `[0, 0]` |
| `getMinimum()`, `getMaximum()`, `getLength()`, `getCenter()` and setters | The length and center setters keep the other property fixed |
| `set(minimum, maximum)`, `set(value)` | Assign |
| `interpolate(t)` | `minimum + t*(maximum - minimum)` |
| `translate(offset)` | Shifts both bounds |
| `scale(factor)` | Multiplies both bounds |
| `scaleAroundCenter(factor)` | Scales the length, keeping the center |
| `extendBy(value)`, `extendBy(interval)` | Grows to include the argument |
| `setUnion(a, b)`, `setIntersection(a, b)`, `static makeUnion(a, b)`, `static makeIntersection(a, b)` | Hull and overlap; an empty overlap yields an empty interval |
| `contains(value)`, `contains(interval)`, `intersects(interval)` | Tests |
| `findIntersection(interval)` | The overlap as `std::optional<Interval>` |

## Random

Header: `<Mathematics/Random.hpp>`.

`Random<T, Generator>` is a collection of static functions that draw uniformly distributed values of type `T` using a
standard library random number engine. `Generator` defaults to `std::mt19937`.

Every function has two forms: one that uses a default generator shared by all calls with the same template arguments,
and one that takes the generator explicitly as the last argument.

```cpp
float f = Random<float>::get(1.f, 2.f);             // uniform in [1, 2]
int die = Random<int>::get(1, 6);                   // uniform integer in [1, 6]
float u = Random<float>::get01();                   // uniform in [0, 1]

Vector3 direction = Random<float>::getVectorOnUnitSphere();
Quaternion orientation = Random<float>::getQuaternion();

std::mt19937 generator(42);                         // reproducible sequence
float g = Random<float>::get01(generator);
```

| Function | Result |
| --- | --- |
| `get(minimum, maximum)` | Scalar in `[minimum, maximum]` (floating point or integer) |
| `getRightOpen(minimum, maximum)` | Scalar in `[minimum, maximum)` |
| `get01()`, `get01RightOpen()` | Scalar in `[0, 1]` and `[0, 1)` |
| `getVectorOnUnitCircle()` | `Vector2` on the unit circle |
| `getVectorOnUnitDisk()` | `Vector2` inside the unit disk |
| `getVectorOnUnitSphere()` | `Vector3` on the surface of the unit sphere |
| `getVectorInUnitSphere()` | `Vector3` inside the unit sphere |
| `getVectorOnUnitSphericalCone(halfAngle)` | Unit `Vector3` within `halfAngle` of the +Z axis |
| `getVectorOnUnitSphericalCone(direction, halfAngle)` | Unit `Vector3` within `halfAngle` of the unit vector `direction` |
| `getRotationMatrix()` | Random rotation as a `Matrix3` |
| `getQuaternion()` | Random unit `Quaternion` |
| `getDefaultGenerator()` | Reference to the default generator, for example to seed it |

The vector, matrix and quaternion functions require a floating point `T`, and return types of that precision.

The default generator is constructed with the engine's default seed, so a program produces the same sequence on every
run unless it is seeded:

```cpp
Random<float>::getDefaultGenerator().seed(std::random_device()());
```

The default generators are static objects without synchronization. In multi-threaded code give each thread its own
generator and use the overloads that take one.

The free function `random(minimum, maximum)` is shorthand for `Random<T>::get(minimum, maximum)`.

## Perlin noise

Header: `<Mathematics/Analysis/Perlin.inl>` (not included by `Mathematics.hpp`). Namespace `mathematics::perlin`.

Ken Perlin's improved gradient noise in one, two and three dimensions.

```cpp
#include <Mathematics/Analysis/Perlin.inl>

float a = perlin::noise(x);
float b = perlin::noise(x, y);
float c = perlin::noise(x, y, z);
float d = perlin::noise(Vector2(x, y));
float e = perlin::noise(Vector3(x, y, z));
```

The result lies approximately in `[-1, 1]` and is zero at integer coordinates. The pattern repeats with a period of
256 along each axis. Sum several octaves at increasing frequency and decreasing amplitude for fractal noise:

```cpp
float fractal(Vector2 p, int octaves)
{
    float sum = 0.f, amplitude = 0.5f;
    for (int i = 0; i < octaves; i++, p *= 2.f, amplitude *= 0.5f)
        sum += amplitude*perlin::noise(p);
    return sum;
}
```

## Romberg integration

Header: `<Mathematics/Analysis/Romberg.inl>` (included by the Bezier curve headers). Namespace
`mathematics::romberg`.

`romberg::estimate<Order>(a, b, f)` computes the definite integral of `f` over `[a, b]` by Romberg's method, i.e.
the trapezoidal rule with Richardson extrapolation. `Order` is the number of refinement levels; the function is
evaluated `2^(Order - 1) + 1` times.

```cpp
float area = romberg::estimate<8>(0.f, Constants<float>::PI, [](float x) { return std::sin(x); });     // 2
```

## Spline bases

Headers: `<Mathematics/Geometry/Bezier.inl>` (included by the Bezier curve headers) and
`<Mathematics/Geometry/Hermite.inl>` (not included by `Mathematics.hpp`).

These functions return the four weights of a cubic segment as any type that can be brace-initialized from four
scalars, typically a `Vector4`. Combining the weights with control points of any type evaluates the curve.

| Function | Returns |
| --- | --- |
| `bezier::evaluate(x0, x1, x2, x3, t)` | A scalar cubic Bezier at `t` |
| `bezier::basis<V>(t)` | Bernstein weights of `p0, p1, p2, p3` |
| `bezier::derivativeBasis<V>(t)` | Weights of `p0, p1, p2, p3` for the first derivative |
| `hermite::basis<V>(t)` | Hermite weights, in the order start point, start tangent, end tangent, end point |

```cpp
#include <Mathematics/Geometry/Hermite.inl>

Vector4 w = bezier::basis<Vector4>(0.5f);                   // (0.125, 0.375, 0.375, 0.125)
Vector3 onBezier = p0*w.x + p1*w.y + p2*w.z + p3*w.w;

Vector4 h = hermite::basis<Vector4>(0.25f);
Vector3 onHermite = start*h.x + startTangent*h.y + endTangent*h.z + end*h.w;
```

For complete cubic Bezier curves with arc length support, use
[`BezierCurve2` and `BezierCurve3`](Geometry.md#bezier-curves).

## Axis and concepts

`<Mathematics/Geometry/Axis.hpp>` defines the enumeration used wherever a principal axis is selected:

```cpp
enum class Axis { X, Y, Z, W };
```

`<Mathematics/Concepts.hpp>` defines the tag types and concepts used across the library:

| Name | Kind | Purpose |
| --- | --- | --- |
| `Uninitialized`, `UNINITIALIZED` | Tag type and constant | Construct without initialization |
| `Identity` | Tag type | Construct an identity matrix, transformation or quaternion |
| `Normalized`, `Unnormalized` | Tag types | Select the variant of a query; see [Conventions](Conventions.md#normalized-and-unnormalized-queries) |
| `Normalization` | Concept | Satisfied by `Normalized` and `Unnormalized` |
| `Arithmetic` | Concept | Satisfied by floating point and integral types |
