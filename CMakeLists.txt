# Mathematics

cmake_minimum_required(VERSION 3.25)

if(NOT DEFINED CACHE{MATHEMATICS_LIB_VERSION})
	string(TIMESTAMP VER_YEAR "%y")
	string(TIMESTAMP VER_MONTH "%m")
	math(EXPR VER_MONTH "${VER_MONTH} + 0")
    set(MATHEMATICS_LIB_VERSION "${VER_YEAR}.${VER_MONTH}" CACHE STRING "Project version (auto-set at first configure)")
endif()

project(Mathematics
	VERSION ${MATHEMATICS_LIB_VERSION}
	LANGUAGES CXX
)

set(MATHEMATICS_LIB_NAME Mathematics)

option(MATHEMATICS_ENABLE_SIMD "Enable SIMD support" ON)
option(MATHEMATICS_ENABLE_QUICKHULL "Enable the Quickhull library" ON)

#set(MATHEMATICS_HAS_SIMD OFF)
set(MATHEMATICS_HAS_SIMD "${MATHEMATICS_ENABLE_SIMD}")
set(MATHEMATICS_HAS_QUICKHULL OFF)

set(SIMD_TARGET "Simd")
set(SIMD_DIR "Simd")
#set(SIMD_DIR "${CMAKE_CURRENT_SOURCE_DIR}/Simd")
set(QUICKHULL_DIR "quickhull")

#if(MATHEMATICS_ENABLE_SIMD)
#    if(EXISTS "${SIMD_DIR}/CMakeLists.txt")
#        if(NOT TARGET "${SIMD_TARGET}")
#            add_subdirectory("${SIMD_DIR}")
#        endif()
#        set(MATHEMATICS_HAS_SIMD ON)
#        message(STATUS "+++ SIMD support: ENABLED")
#    else()
#        message(STATUS "--- SIMD support: DISABLED (Source files not found)")
#    endif()
#else()
#    message(STATUS "--- SIMD support: DISABLED")
#endif()

if(MATHEMATICS_HAS_SIMD)
	if(NOT TARGET "${SIMD_TARGET}")
		add_subdirectory("${SIMD_DIR}")
	endif()
endif()

add_library(${MATHEMATICS_LIB_NAME} INTERFACE)

set(MATHEMATICS_HEADERS
	Include/Mathematics/Transform/AffineTransform.hpp
	Include/Mathematics/Geometry/Axis.hpp
	Include/Mathematics/Geometry/AxisAlignedBox.hpp
	Include/Mathematics/Geometry/AxisAlignedRectangle.hpp
	Include/Mathematics/Geometry/Bezier.inl
	Include/Mathematics/Geometry/BezierCurve2.hpp
	Include/Mathematics/Geometry/BezierCurve3.hpp
	Include/Mathematics/Geometry/Capsule.hpp
	Include/Mathematics/Geometry/Circle2.hpp
	Include/Mathematics/Geometry/Cone.hpp
	Include/Mathematics/Constants.hpp
	Include/Mathematics/Analysis/Containment.inl
	Include/Mathematics/Geometry/Cylinder.hpp
	Include/Mathematics/Analysis/Distances.inl
	Include/Mathematics/Geometry/Ellipsoid.hpp
	Include/Mathematics/Transform/Euler.hpp
	Include/Mathematics/Transform/EulerOrder.hpp
	Include/Mathematics/Geometry/HalfSpace.hpp
	Include/Mathematics/Geometry/Hermite.inl
	Include/Mathematics/Analysis/Intersections.inl
	Include/Mathematics/Interval.hpp
	Include/Mathematics/LineSegment2.hpp
	Include/Mathematics/LineSegment3.hpp
	Include/Mathematics/Line2.hpp
	Include/Mathematics/Line3.hpp
	Include/Mathematics/Algebra/LU.inl
	Include/Mathematics/Mathematics.hpp
	Include/Mathematics/Algebra/Matrix2.hpp
	Include/Mathematics/Algebra/Matrix3.hpp
	Include/Mathematics/Algebra/Matrix4.hpp
	Include/Mathematics/Geometry/OrientedBox.hpp
	Include/Mathematics/Analysis/Perlin.inl
	Include/Mathematics/Geometry/Plane.hpp
	Include/Mathematics/Geometry/Quadrilateral3.hpp
	Include/Mathematics/Algebra/Quaternion.hpp
	Include/Mathematics/Random.hpp
	Include/Mathematics/Geometry/Ray2.hpp
	Include/Mathematics/Geometry/Ray3.hpp
	Include/Mathematics/Analysis/Romberg.inl
	Include/Mathematics/Scalar.hpp
	Include/Mathematics/Geometry/Sphere.hpp
	Include/Mathematics/Geometry/SymmetricFrustum.hpp
	Include/Mathematics/Geometry/Tetrahedron.hpp
	Include/Mathematics/Geometry/Torus.hpp
	Include/Mathematics/Geometry/Triangle3.hpp
	Include/Mathematics/Geometry/Triangulation.inl
	Include/Mathematics/Algebra/Vector2.hpp
	Include/Mathematics/Algebra/Vector3.hpp
	Include/Mathematics/Algebra/Vector4.hpp
	Include/Mathematics/Transform/YawPitchRoll.hpp
)

target_sources(${MATHEMATICS_LIB_NAME}
	INTERFACE
    	FILE_SET HEADERS
		BASE_DIRS Include
    	FILES
			${MATHEMATICS_HEADERS}
)

#source_group(TREE "${CMAKE_CURRENT_SOURCE_DIR}/Include" FILES ${MATHEMATICS_HEADERS})

target_compile_features(${MATHEMATICS_LIB_NAME} 
    INTERFACE 
        cxx_std_20
)

target_compile_definitions(${MATHEMATICS_LIB_NAME}
    INTERFACE 
        MATHEMATICS_HAS_SIMD=$<BOOL:${MATHEMATICS_HAS_SIMD}>
        MATHEMATICS_HAS_QUICKHULL=$<BOOL:${MATHEMATICS_HAS_QUICKHULL}>
)
