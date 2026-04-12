/*
 * Project: Mathematics Library
 * Author: Pawel Mrochen (owner)
 * Refactored & Modularized by: hypernova-developer
 * * Description: 
 * This is a header-only C++20 library for linear algebra and geometry.
 * The directory structure has been reorganized into logical modules 
 * to improve scalability and ease of integration for modern C++ engines.
 */

#pragma once

// ============================================================================
// Core & Mathematical Foundations
#include "Core/Constants.hpp"
#include "Core/Scalar.hpp"
#include "Core/Interval.hpp"
#include "Core/SymmetricFrustum.hpp"

// ============================================================================
// Linear Algebra (Vectors, Matrices, Quaternions)
// ============================================================================
#include "Algebra/Vector2.hpp"
#include "Algebra/Vector3.hpp"
#include "Algebra/Vector4.hpp"
#include "Algebra/Matrix2.hpp"
#include "Algebra/Matrix3.hpp"
#include "Algebra/Matrix4.hpp"
#include "Algebra/Quaternion.hpp"

// ============================================================================
// Transformations & Coordinate Systems
// ============================================================================
#include "Transform/Axis.hpp"
#include "Transform/EulerOrder.hpp"
#include "Transform/AffineTransform.hpp"
#include "Transform/YawPitchRoll.hpp"
#include "Transform/Euler.hpp"

// ============================================================================
// Geometry - Primitive Elements
// ============================================================================
#include "Geometry/Primitives/Line2.hpp"
#include "Geometry/Primitives/Ray2.hpp"
#include "Geometry/Primitives/Segment2.hpp"
#include "Geometry/Primitives/Line3.hpp"
#include "Geometry/Primitives/Ray3.hpp"
#include "Geometry/Primitives/Segment3.hpp"
#include "Geometry/Primitives/HalfSpace.hpp"
#include "Geometry/Primitives/Plane.hpp"

// ============================================================================
// Geometry - Volumetric Shapes & Boxes
// ============================================================================
#include "Geometry/Boxes/AxisAlignedRectangle.hpp"
#include "Geometry/Boxes/AxisAlignedBox.hpp"
#include "Geometry/Boxes/OrientedBox.hpp"
#include "Geometry/Shapes/Circle2.hpp"
#include "Geometry/Shapes/Sphere.hpp"
#include "Geometry/Shapes/Ellipsoid.hpp"
#include "Geometry/Shapes/Torus.hpp"
#include "Geometry/Shapes/Cylinder.hpp"
#include "Geometry/Shapes/Capsule.hpp"
#include "Geometry/Shapes/Cone.hpp"
#include "Geometry/Shapes/Tetrahedron.hpp"
#include "Geometry/Shapes/Quadrilateral3.hpp"

// ============================================================================
// Geometry - Specialized Components (Triangle & Triangulation)
// ============================================================================
// Note: Kept in a separate module for algorithmic extensibility and 
// to prevent circular dependencies during complex intersection tests.
#include "Geometry/Triangle/Triangle3.hpp"

// ============================================================================
// Numerical Analysis & Algorithms (Inlined Implementations)
// ============================================================================
#include "Analysis/Containment.inl"
#include "Analysis/Distances.inl"
#include "Analysis/Intersections.inl"
#include "Analysis/Romberg.inl"
#include "Analysis/Perlin.inl"