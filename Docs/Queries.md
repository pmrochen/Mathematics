# Queries

[← Back to the README](../README.md)

- [Kinds of queries](#kinds-of-queries)
- [Containment](#containment)
- [Distance and closest point](#distance-and-closest-point)
- [Intersection tests](#intersection-tests)
- [Finding intersections](#finding-intersections)
- [Classification against half-spaces](#classification-against-half-spaces)
- [Low-level functions](#low-level-functions)

Queries are member functions of the primitives described in [Geometry](Geometry.md) and [Meshes](Meshes.md).

## Kinds of queries

| Member | Returns | Question answered |
| --- | --- | --- |
| `a.contains(b)` | `bool` | Is `b` entirely inside `a`? |
| `a.getDistanceTo(point)` | scalar | How far is the point from `a`? `0` when the point is inside a solid |
| `a.getSignedDistanceTo(point)` | scalar | The same, negative inside or behind |
| `a.getClosestPoint(point)` | point | Which point of `a` is nearest? |
| `a.intersects(b)` | `bool` | Do `a` and `b` have at least one point in common? |
| `a.findIntersection(b)` | `std::optional` | Where does the line, ray or segment `a` meet `b`? |
| `a.findIntersections(mesh)` | `std::vector` | All hits of a line, ray or segment on a triangle mesh |
| `a.classify(halfSpace)` | `int` | Is the box inside, outside or straddling the half-space? |

Solid shapes are treated as filled volumes: a sphere inside a box intersects it, and a point inside a sphere is at
distance zero.

> [!NOTE]
> Queries on lines, rays, planes and half-spaces assume **normalized** directions and normals unless the
> `<Unnormalized>` overload is used. See [Conventions](Conventions.md#normalized-and-unnormalized-queries).

## Containment

| Type | `contains(...)` accepts |
| --- | --- |
| `Interval` | value, `Interval` |
| `AxisAlignedRectangle` | point, `AxisAlignedRectangle`, `Circle2` |
| `Circle2` | point |
| `Plane` | point (lies on the plane) |
| `HalfSpace` | point |
| `AxisAlignedBox` | point, `AxisAlignedBox`, `Sphere` |
| `OrientedBox` | point |
| `Sphere` | point |
| `Ellipsoid` | point |
| `Cylinder` | point |
| `Cone` | point |
| `SymmetricFrustum` | point |
| `ConvexPolyhedron` | point, optionally with a tolerance |
| `TriangleMesh` | point (the mesh must be convex with outward-facing triangles) |

```cpp
AxisAlignedBox room(Vector3(0.f, 0.f, 0.f), Vector3(10.f, 3.f, 8.f));

room.contains(Vector3(5.f, 1.f, 4.f));                              // true
room.contains(Sphere(Vector3(5.f, 1.f, 4.f), 2.f));                 // false: the sphere pokes through the floor
room.contains(AxisAlignedBox(Vector3(1.f), Vector3(2.f)));          // true
```

## Distance and closest point

| Type | `getDistanceTo(point)` | `getSignedDistanceTo(point)` | `getClosestPoint(point)` |
| --- | :-: | :-: | :-: |
| `Line2` | ✓ | ✓ | ✓ |
| `Ray2`, `LineSegment2` | ✓ | | ✓ |
| `Line3`, `Ray3`, `LineSegment3` | ✓ | | ✓ |
| `Plane`, `HalfSpace` | ✓ | ✓ | |
| `Triangle3` | ✓ | | ✓ |
| `AxisAlignedRectangle`, `AxisAlignedBox`, `OrientedBox` | ✓ | | ✓ |
| `Circle2`, `Sphere` | ✓ | ✓ | |
| `SymmetricFrustum` | ✓ | | ✓ |

`Line2` additionally offers `getDistanceTo(line)` and `getSignedDistanceTo(line)` for two lines.

```cpp
Triangle3 triangle(Vector3(-1.f, -1.f, 0.f), Vector3(1.f, -1.f, 0.f), Vector3(0.f, 1.f, 0.f));

triangle.getClosestPoint(Vector3(0.f, 5.f, 0.f));       // (0, 1, 0): the top vertex
triangle.getDistanceTo(Vector3(0.f, 5.f, 0.f));         // 4

Sphere sphere(Vector3::ZERO, 1.f);
sphere.getDistanceTo(Vector3(0.5f, 0.f, 0.f));          // 0: inside
sphere.getSignedDistanceTo(Vector3(0.5f, 0.f, 0.f));    // -0.5
```

## Intersection tests

`a.intersects(b)` is available for the pairs marked below. Read each table as `row.intersects(column)`.

### 2D

| | `Line2` | `Ray2` | `LineSegment2` | `AxisAlignedRectangle` | `Circle2` |
| --- | :-: | :-: | :-: | :-: | :-: |
| **`Line2`** | ✓ | ✓ | ✓ | ✓ | ✓ |
| **`Ray2`** | ✓ | | | ✓ | ✓ |
| **`LineSegment2`** | ✓ | | ✓ | ✓ | ✓ |
| **`AxisAlignedRectangle`** | | | | ✓ | ✓ |
| **`Circle2`** | | | | ✓ | ✓ |

The linear types are the left-hand side: write `line.intersects(circle)`, not `circle.intersects(line)`.

### 3D: lines, rays and segments

| | `Plane` | `HalfSpace` | `Triangle3` | `AxisAlignedBox` | `OrientedBox` | `Sphere` | `Ellipsoid` | `TriangleMesh` |
| --- | :-: | :-: | :-: | :-: | :-: | :-: | :-: | :-: |
| **`Line3`** | ✓ | | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ |
| **`Ray3`** | ✓ | | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ |
| **`LineSegment3`** | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ |

Again the linear type is the left-hand side. Meshes are passed by pointer: `ray.intersects(mesh.get())`.

### 3D: surfaces and volumes

| | `Plane` | `HalfSpace` | `Triangle3` | `AxisAlignedBox` | `OrientedBox` | `Sphere` | `Ellipsoid` | `Cone` | `SymmetricFrustum` | `ConvexPolyhedron` |
| --- | :-: | :-: | :-: | :-: | :-: | :-: | :-: | :-: | :-: | :-: |
| **`Plane`** | | | ✓ | ✓ | ✓ | ✓ | ✓ | | | |
| **`HalfSpace`** | | | ✓ | ✓ | ✓ | ✓ | | | | |
| **`Triangle3`** | ✓ | ✓ | | ✓ | ✓ | ✓ | | | | |
| **`AxisAlignedBox`** | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | | | ✓ | ✓ |
| **`OrientedBox`** | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | | | ✓ | ✓ |
| **`Sphere`** | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | | ✓ | ✓ | ✓ |
| **`Ellipsoid`** | ✓ | | | | | | | | | |
| **`Cone`** | | | | | | ✓ | | | | |
| **`SymmetricFrustum`** | | | | ✓ | ✓ | ✓ | | | | |
| **`ConvexPolyhedron`** | | | | ✓ | ✓ | ✓ | | | | |

`ConvexPolyhedron` arguments are passed by pointer.

```cpp
AxisAlignedBox box(Vector3(-1.f), Vector3(1.f));
Sphere sphere(Vector3(1.5f, 0.f, 0.f), 1.f);
SymmetricFrustum view(Vector3::ZERO, Matrix3::IDENTITY, Vector2(1.f, 0.75f), 1.f, 100.f);

box.intersects(sphere);                                         // true
box.intersects(Plane(Vector3::UNIT_Y, 0.f));                    // true: the plane y = 0 cuts the box
view.intersects(Sphere(Vector3(0.f, 0.f, 200.f), 1.f));         // false: beyond the far plane
```

The `ConvexPolyhedron` tests check the object against each half-space in turn, as is usual for visibility culling.
They never reject an object that touches the volume, but may report an intersection for an object that is just
outside near an edge or corner.

## Finding intersections

Lines, rays and segments can report *where* they meet another primitive.

| Target | Result of `findIntersection` | Meaning |
| --- | --- | --- |
| Another line, ray or segment (2D) | `std::optional<T>` | Parameter `t` on the calling object |
| `Plane`, `Triangle3` | `std::optional<T>` | Parameter `t` of the single hit |
| `AxisAlignedRectangle`, `Circle2`, `AxisAlignedBox`, `OrientedBox`, `Sphere`, `Ellipsoid` | `std::optional<Interval<T>>` | The range of `t` that lies inside the shape: entry and exit |
| `TriangleMesh` (`findIntersections`) | `std::vector<T>` | One `t` per intersected triangle, in triangle order |

The parameters refer to the calling object's own parameterization, so `evaluate(t)` gives the point:

| Caller | Parameter range of the results |
| --- | --- |
| Line | unbounded, may be negative |
| Ray | `t >= 0`; a ray starting inside a solid reports an interval beginning at `0` |
| Line segment | `0 <= t <= 1` |

```cpp
Ray3 ray(Vector3(0.f, 0.f, -5.f), Vector3::UNIT_Z);
Sphere sphere(Vector3::ZERO, 1.f);

if (auto hit = ray.findIntersection(sphere))
{
    Vector3 entry = ray.evaluate(hit->getMinimum());        // (0, 0, -1), t = 4
    Vector3 exit = ray.evaluate(hit->getMaximum());         // (0, 0, 1), t = 6
}

Plane ground(Vector3::UNIT_Y, 0.f);
if (auto t = Ray3(Vector3(0.f, 5.f, 0.f), -Vector3::UNIT_Y).findIntersection(ground))
    float distance = *t;                                    // 5

LineSegment3 segment(Vector3(-5.f, 0.f, 0.f), Vector3(5.f, 0.f, 0.f));
if (auto span = segment.findIntersection(AxisAlignedBox(Vector3(-1.f), Vector3(1.f))))
    Vector3 first = segment.evaluate(span->getMinimum());   // (-1, 0, 0), t = 0.4
```

An empty optional means there is no intersection: the objects miss each other, a line is parallel to a plane, or the
hit lies outside the valid parameter range of a ray or segment.

### 2D example

```cpp
Line2 a(Vector2(0.f, 0.f), Vector2(1.f, 1.f));
Line2 b(Vector2(0.f, 2.f), Vector2(1.f, -1.f));

if (auto t = a.findIntersection(b))
    Vector2 crossing = a.evaluate(*t);                      // (1, 1)
```

### Picking in a triangle mesh

```cpp
std::unique_ptr<TriangleMesh> mesh(TriangleMesh::makeBox(Vector3(2.f, 2.f, 2.f)));
Ray3 ray(Vector3(0.25f, 0.5f, -5.f), Vector3::UNIT_Z);

std::vector<float> hits = ray.findIntersections(mesh.get());
if (!hits.empty())
{
    float nearest = *std::min_element(hits.begin(), hits.end());   // 4: the front face
    Vector3 point = ray.evaluate(nearest);
}
```

The results are not sorted, and a ray that passes exactly through an edge shared by two triangles reports both.

## Classification against half-spaces

`AxisAlignedBox::classify(halfSpace)` and `OrientedBox::classify(halfSpace)` return

| Value | Meaning |
| --- | --- |
| `-1` | The box is completely inside the half-space |
| `0` | The box straddles the boundary plane |
| `1` | The box is completely outside |

This is the building block for hierarchical culling. `ConvexPolyhedron::intersects(box, mask)` uses it with a bit
mask of the half-spaces that still need testing: on return, the bits of the half-spaces that contain the box entirely
are cleared, so they can be skipped for everything nested inside that box.

```cpp
std::unique_ptr<ConvexPolyhedron> volume(ConvexPolyhedron::from(view));

unsigned int mask = volume->getHalfSpaceMask<unsigned int>();      // all half-spaces active
if (volume->intersects(nodeBounds, mask))
{
    // visit the children with the reduced mask;
    // mask == 0 means the node is entirely inside and its children need no further tests
}
```

## Low-level functions

The member queries are thin wrappers over free functions that take raw parameters (points, directions, radii, ...)
instead of primitive objects. They are useful when the data is already stored in another form.

| Namespace | Header | Naming | Examples |
| --- | --- | --- | --- |
| `mathematics::containment` | `Analysis/Containment.inl` | `testAPoint` returning `bool` | `testCylinderPoint`, `testConePoint` |
| `mathematics::distances` | `Analysis/Distances.inl` | `getAB`, with `Squared`, `Signed` and `Normalized` variants | `getPointPlane`, `getPointNormalizedPlaneSigned`, `getPointTriangle`, `getPointSymmetricFrustum` |
| `mathematics::intersections` | `Analysis/Intersections.inl` | `testAB` returning `bool`; `findAB<Result>` returning the result type given as the first template argument | `testLinePlane`, `testOrientedBoxOrientedBox`, `findLineTriangle`, `findLineAxisAlignedBox` |

```cpp
float d = distances::getPointTriangle(point, v0, v1, v2);

auto t = intersections::findLineTriangle<std::optional<float>>(origin, direction, v0, v1, v2);
```

`getPointTriangle` and `getPointSymmetricFrustum` accept an optional trailing `Vector3*` that receives the closest
point. These headers are included automatically by the primitives that use them.
