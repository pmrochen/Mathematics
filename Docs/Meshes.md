# Meshes

[← Back to the README](../README.md)

- [Ownership](#ownership)
- [ConvexPolyhedron](#convexpolyhedron)
- [TriangleMesh](#trianglemesh)
- [TetrahedralMesh](#tetrahedralmesh)
- [Polygon triangulation](#polygon-triangulation)

Unlike the fixed-size primitives in [Geometry](Geometry.md), the types on this page own dynamically sized data.

| Type | Header | Stores |
| --- | --- | --- |
| `ConvexPolyhedron` | `<Mathematics/Geometry/ConvexPolyhedron.hpp>` | A list of half-spaces whose intersection is the volume |
| `TriangleMesh` | `<Mathematics/Geometry/TriangleMesh.hpp>` | Vertices and, optionally, triangle indices |
| `TetrahedralMesh` | `<Mathematics/Geometry/TetrahedralMesh.hpp>` | Vertices and, optionally, tetrahedron indices |

The aliases `ConvexPolyhedron`, `TriangleMesh` and `TetrahedralMesh` use the default precision and, for the meshes,
`std::uint32_t` indices. Other combinations are available through the templates, for example
`templates::TriangleMesh<double, std::uint16_t>`.

## Ownership

The three types are regular copyable and movable classes and can be used as values. They also derive from
`std::enable_shared_from_this`, and their factory functions (`from`, `makeBox`, `makeTetrahedron`, `clone`) return a
**raw pointer to a new object that the caller owns**. Wrap the result in a smart pointer immediately:

```cpp
std::unique_ptr<TriangleMesh> box(TriangleMesh::makeBox(Vector3(2.f, 2.f, 2.f)));
std::shared_ptr<ConvexPolyhedron> volume(ConvexPolyhedron::from(frustum));
```

Queries that take one of these types accept a pointer: `ray.findIntersections(box.get())`,
`sphere.intersects(volume.get())`.

## ConvexPolyhedron

A convex volume described as the intersection of half-spaces. The normals of the half-spaces point outwards. The
volume does not have to be bounded: two parallel half-spaces describe an infinite slab.

```cpp
ConvexPolyhedron slab{ HalfSpace(Vector3::UNIT_Y, -1.f), HalfSpace(-Vector3::UNIT_Y, -1.f) };   // -1 <= y <= 1

slab.contains(Vector3(100.f, 0.5f, 0.f));               // true
slab.intersects(Sphere(Vector3(0.f, 3.f, 0.f), 1.f));   // false
```

### Construction

| Source | Description |
| --- | --- |
| `ConvexPolyhedron{ h0, h1, ... }` | From an initializer list of half-spaces |
| `ConvexPolyhedron(first, last)` | From an iterator range of half-spaces |
| `ConvexPolyhedron(halfSpaces)`, `ConvexPolyhedron(planes)` | From a `std::vector` of half-spaces or planes, or a pointer and a count |
| `ConvexPolyhedron::from(box)` | The six faces of an `AxisAlignedBox` or `OrientedBox` |
| `ConvexPolyhedron::from(frustum)` | The six faces of a `SymmetricFrustum`, in the order left, right, bottom, top, near, far |
| `clone()` | A copy on the heap |

### Members

| Member | Description |
| --- | --- |
| `isEmpty()`, `clear()` | No half-spaces |
| `assign(polyhedron)`, `append(polyhedron)` | Replace with, or add, the half-spaces of another polyhedron; appending intersects the two volumes |
| `getHalfSpaces()`, `getHalfSpaceCount()`, `getHalfSpace(i)` | Access; the public member `halfSpaces` is the underlying `std::vector` |
| `setHalfSpace(i, h)`, `insertHalfSpace(i, h)` | Throw `std::out_of_range` for an invalid index |
| `addHalfSpace(h)`, `addHalfSpaces(...)`, `removeHalfSpaceAt(i)`, `removeAllHalfSpaces()`, `setHalfSpaceCount(n)` | Editing |
| `getPlanes()`, `copyPlanes(out)` | The boundary planes |
| `getHalfSpaceMask<U>()` | A bit mask with one bit set per half-space, for [masked culling](Queries.md#classification-against-half-spaces) |
| `computeVertices(tolerance, pivotTolerance)` | The corner points of a bounded polyhedron, found by intersecting triples of planes |
| `computeAxisAlignedBoundingBox(tolerance, pivotTolerance)` | Bounds of those vertices |
| `translate(offset)`, `transform(m, orthogonal = false)` | Moves every half-space |
| `contains(point)`, `contains(point, tolerance)` | Point containment |
| `intersects(box)`, `intersects(box, mask)`, `intersects(sphere)` | Volume tests against `AxisAlignedBox`, `OrientedBox` and `Sphere` |

```cpp
std::unique_ptr<ConvexPolyhedron> cube(ConvexPolyhedron::from(AxisAlignedBox(Vector3(-1.f), Vector3(1.f))));

cube->getHalfSpaceCount();              // 6
cube->computeVertices().size();         // 8
cube->computeAxisAlignedBoundingBox();  // (-1, -1, -1) to (1, 1, 1)
```

## TriangleMesh

A list of vertices with an optional index list.

- **Without indices** the mesh is a *triangle soup*: every three consecutive vertices form a triangle.
- **With indices** every three consecutive indices form a triangle, and vertices can be shared.

`getTriangleCount()`, `getTriangle(i)` and all computations work in both modes.

```cpp
std::vector<Vector3> vertices{ Vector3(0.f, 0.f, 0.f), Vector3(1.f, 0.f, 0.f), Vector3(1.f, 1.f, 0.f), Vector3(0.f, 1.f, 0.f) };
std::vector<unsigned int> indices{ 0, 1, 2,  0, 2, 3 };

TriangleMesh quad(vertices, indices);

quad.getTriangleCount();                    // 2
quad.computeSurfaceArea();                  // 1
quad.computeAxisAlignedBoundingBox();       // (0, 0, 0) to (1, 1, 0)
```

### Construction

| Source | Description |
| --- | --- |
| `TriangleMesh(vertices)` | A triangle soup from a `std::vector`, an initializer list, an iterator range, or a pointer and a count |
| `TriangleMesh(vertices, indices)` | An indexed mesh; index vectors of another integer type are converted |
| `TriangleMesh(triangles)` | From a `std::vector<Triangle3>` or a pointer and a count |
| `TriangleMesh::makeBox(dimensions)` | A box centered at the origin: 8 vertices, 12 triangles |
| `TriangleMesh::makeTetrahedron(edgeLength)` | A regular tetrahedron centered at the origin |
| `TriangleMesh::from(box)`, `from(frustum)` | The surface of an `AxisAlignedBox`, `OrientedBox` or `SymmetricFrustum` |
| `TriangleMesh::from(convexPolyhedron, ...)` | The surface of a bounded convex polyhedron |
| `clone()` | A copy on the heap |

The vector overloads also accept rvalues and then take over the storage without copying.

### Members

| Group | Members |
| --- | --- |
| State | `isEmpty()`, `clear()`, `assign(mesh)`, `append(mesh)`, `hasValidVertexIndices()` |
| Vertices | `getVertexCount()`, `getVertex(i)`, `setVertex(i, v)`, `addVertex(v)`, `addVertices(...)`, `insertVertex(i, v)`, `removeVertexAt(i)`, `removeAllVertices()`, `setVertexCount(n)` |
| Indices | `getVertexIndexCount()`, `getVertexIndex(i)`, `getVertexIndex(triangle, corner)`, `setVertexIndex(i, value)`, `addVertexIndex(value)`, `addVertexIndices(v0, v1, v2)`, `addVertexIndices(...)`, `insertVertexIndex(i, value)`, `insertVertexIndices(i, v0, v1, v2)`, `removeVertexIndexAt(i)`, `removeAllVertexIndices()`, `setVertexIndexCount(n)` |
| Triangles | `getTriangleCount()` (synonym `getPrimitiveCount()`), `getTriangle(i)`, `getTriangles()`, `copyTriangles(out)` |
| Planes | `getHalfSpaces()`, `copyHalfSpaces(out)`, `getPlanes()`, `copyPlanes(out)` - one per triangle |
| Measures | `computeCentroid()`, `computeSurfaceArea()`, `computeVolume()`, `computeAxisAlignedBoundingBox()`, `computeBoundingSphere()` |
| Topology | `isConvex()`, `computeVertexIndices()`, `weldVertices(tolerance)` |
| Transformation | `translate(offset)`, `transform(matrix)`, `transform(transformation)` |
| Queries | `contains(point)`; `intersects(mesh)` and `findIntersections(mesh)` on `Line3`, `Ray3` and `LineSegment3` |

The public members `vertices` and `indices` are the underlying `std::vector`s. Setters and inserters with an index
throw `std::out_of_range` when it is invalid.

Notes:

- `computeVolume()` and `computeCentroid()` are meaningful for closed meshes with consistent winding.
- `contains(point)` tests the point against the plane of every triangle and is therefore only valid for **convex**
  meshes whose triangles face outwards (counter-clockwise when seen from outside).
- `computeVertexIndices()` turns a triangle soup into a trivially indexed mesh (`0, 1, 2, ...`) and returns `false` if
  the mesh already has indices.
- `weldVertices(tolerance)` merges vertices that are approximately equal, creating or updating the indices, and
  returns `true` if anything was merged. It compares all pairs of vertices, so its cost grows quadratically with the
  vertex count.

```cpp
// Two triangles given as a soup of six vertices
TriangleMesh soup{ Vector3(0.f, 0.f, 0.f), Vector3(1.f, 0.f, 0.f), Vector3(1.f, 1.f, 0.f),
                   Vector3(0.f, 0.f, 0.f), Vector3(1.f, 1.f, 0.f), Vector3(0.f, 1.f, 0.f) };

soup.weldVertices();            // true
soup.getVertexCount();          // 4
soup.getVertexIndexCount();     // 6
```

```cpp
std::unique_ptr<TriangleMesh> box(TriangleMesh::makeBox(Vector3(2.f, 2.f, 2.f)));

box->computeVolume();           // 8
box->computeSurfaceArea();      // 24
box->isConvex();                // true
box->contains(Vector3(0.5f));   // true

box->transform(AffineTransform::makeTranslation(Vector3(10.f, 0.f, 0.f)));
box->computeCentroid();         // (10, 0, 0)
```

Ray casting against a mesh is shown in [Queries](Queries.md#picking-in-a-triangle-mesh).

## TetrahedralMesh

A volumetric mesh with the same structure as `TriangleMesh`, with four vertices or indices per element.

```cpp
std::vector<Tetrahedron> cells{ Tetrahedron(Vector3(0.f, 0.f, 0.f), Vector3(1.f, 0.f, 0.f),
                                            Vector3(0.f, 1.f, 0.f), Vector3(0.f, 0.f, 1.f)) };
TetrahedralMesh mesh(cells);

mesh.getTetrahedronCount();     // 1
mesh.computeVolume();           // 0.166667
mesh.computeCentroid();         // (0.25, 0.25, 0.25)
```

It offers the same construction options (vertices, vertices with indices, a list of `Tetrahedron`), the same vertex
and index editing members (`addVertexIndices` and `insertVertexIndices` take four indices), and:

| Member | Description |
| --- | --- |
| `getTetrahedronCount()` (synonym `getPrimitiveCount()`), `getTetrahedron(i)`, `getTetrahedrons()`, `copyTetrahedrons(out)` | Element access |
| `computeCentroid()`, `computeVolume()` | Measures |
| `computeAxisAlignedBoundingBox()`, `computeBoundingSphere()` | Bounds |
| `computeVertexIndices()`, `weldVertices(tolerance)` | As for `TriangleMesh` |
| `translate(offset)`, `transform(matrix)`, `transform(transformation)` | Transformation |

## Polygon triangulation

Simple polygons (convex or concave, without holes) are triangulated by ear clipping. The result is a list of
`(n - 2)*3` indices into the input vertices, and a flag that is `false` when the polygon could not be triangulated,
for example because it is self-intersecting.

### Planar polygons in 3D

`Triangle3::triangulate` projects the polygon onto its own plane first:

```cpp
std::vector<Vector3> polygon{ Vector3(0.f, 0.f, 0.f), Vector3(2.f, 0.f, 0.f), Vector3(2.f, 2.f, 0.f),
                              Vector3(1.f, 1.f, 0.f), Vector3(0.f, 2.f, 0.f) };   // concave pentagon
std::vector<std::size_t> indices;

auto [end, succeeded] = Triangle3::triangulate(polygon.begin(), polygon.end(), std::back_inserter(indices));
// succeeded == true; indices == { 1, 2, 3,  3, 4, 0,  0, 1, 3 }
```

The index type defaults to `std::size_t`; pass another as the first template argument:
`Triangle3::triangulate<unsigned int>(...)`.

A second overload triangulates a polygon whose outline is given by an index loop into a shared vertex array. The
output indices then refer to the vertex array as well. This overload needs a writable range of `(n - 2)*3` elements
rather than an insert iterator:

```cpp
std::vector<unsigned int> outline{ 4, 3, 2, 1, 0 };             // the polygon, as indices into 'polygon'
std::vector<unsigned int> triangles((outline.size() - 2)*3);

auto [last, ok] = Triangle3::triangulate<unsigned int>(polygon.begin(), polygon.end(),
                                                       outline.begin(), outline.end(), triangles.begin());
```

### Polygons in 2D

The underlying algorithm works on `Vector2` vertices and is available directly in `mathematics::triangulation`
(`<Mathematics/Geometry/Triangulation.inl>`):

```cpp
std::vector<Vector2> outline2{ Vector2(0.f, 0.f), Vector2(2.f, 0.f), Vector2(2.f, 2.f), Vector2(1.f, 1.f), Vector2(0.f, 2.f) };
std::vector<int> triangles2;

auto [last2, ok2] = triangulation::triangulate(outline2.begin(), outline2.end(), std::back_inserter(triangles2));
```
