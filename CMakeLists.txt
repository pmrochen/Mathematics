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

target_include_directories(${MATHEMATICS_LIB_NAME}		# to remove after moving header files to Include folder
	INTERFACE
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
        $<INSTALL_INTERFACE:include>
)

set(MATHEMATICS_HEADERS
	AffineTransform.hpp
	Axis.hpp
	AxisAlignedBox.hpp
	AxisAlignedRectangle.hpp
	Bezier.inl
	BezierCurve2.hpp
	BezierCurve3.hpp
	Capsule.hpp
	Circle2.hpp
	Cone.hpp
	Constants.hpp
	Containment.inl
	Cylinder.hpp
	Distances.inl
	Ellipsoid.hpp
	Euler.hpp
	EulerOrder.hpp
	HalfSpace.hpp
	Intersections.inl
	Interval.hpp
	LineSegment2.hpp
	LineSegment3.hpp
	Line2.hpp
	Line3.hpp
	LU.inl
	Mathematics.hpp
	Matrix2.hpp
	Matrix3.hpp
	Matrix4.hpp
	OrientedBox.hpp
	Perlin.inl
	Plane.hpp
	Quadrilateral3.hpp
	Quaternion.hpp
	Ray2.hpp
	Ray3.hpp
	Romberg.inl
	Scalar.hpp
	Sphere.hpp
	SymmetricFrustum.hpp
	Tetrahedron.hpp
	Torus.hpp
	Triangle3.hpp
	Triangulation.inl
	Vector2.hpp
	Vector3.hpp
	Vector4.hpp
	YawPitchRoll.hpp
)

target_sources(${MATHEMATICS_LIB_NAME}
	INTERFACE
    	FILE_SET HEADERS
		#BASE_DIRS Include
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
