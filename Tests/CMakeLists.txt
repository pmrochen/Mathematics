# Mathematics/Tests

cmake_minimum_required(VERSION 3.25)

include(FetchContent)

FetchContent_Declare(doctest
	URL https://raw.githubusercontent.com/doctest/doctest/v2.5.3/doctest/doctest.h
	DOWNLOAD_NO_EXTRACT TRUE
	SOURCE_DIR downloaded/doctest
)

FetchContent_MakeAvailable(doctest)

set(MATHEMATICS_TESTS_NAME MathematicsTests)

add_executable(${MATHEMATICS_TESTS_NAME} EXCLUDE_FROM_ALL)

set_target_properties(${MATHEMATICS_TESTS_NAME} PROPERTIES
	CXX_STANDARD 20
	CXX_STANDARD_REQUIRED ON
	CXX_EXTENSIONS OFF
)

# Linking the interface target brings in the Mathematics headers, the Simd dependency and the
# MATHEMATICS_* compile definitions, so the tests exercise the same configuration as consumers do.
target_link_libraries(${MATHEMATICS_TESTS_NAME} PRIVATE
	${MATHEMATICS_LIB_NAME}
)

target_include_directories(${MATHEMATICS_TESTS_NAME} PRIVATE
	${CMAKE_CURRENT_SOURCE_DIR} # TestUtilities.hpp
	${CMAKE_CURRENT_BINARY_DIR}/downloaded # doctest.h
)

target_sources(${MATHEMATICS_TESTS_NAME} PRIVATE
	Main.cpp
	TestUtilities.hpp
	ScalarTests.cpp
	IntervalTests.cpp
	RandomTests.cpp
	Algebra/Vector2Tests.cpp
	Algebra/Vector3Tests.cpp
	Algebra/Vector4Tests.cpp
	Algebra/Matrix2Tests.cpp
	Algebra/Matrix3Tests.cpp
	Algebra/Matrix4Tests.cpp
	Algebra/QuaternionTests.cpp
	Algebra/LUTests.cpp
	Transform/AffineTransformTests.cpp
	Transform/EulerTests.cpp
	Geometry/Linear2Tests.cpp
	Geometry/Linear3Tests.cpp
	Geometry/Shapes2Tests.cpp
	Geometry/PlaneTests.cpp
	Geometry/PolygonTests.cpp
	Geometry/BoxTests.cpp
	Geometry/SolidTests.cpp
	Geometry/BezierTests.cpp
	Geometry/MeshTests.cpp
	Analysis/AnalysisTests.cpp
)

add_test(
	NAME "${MATHEMATICS_TESTS_NAME}"
	COMMAND ${MATHEMATICS_TESTS_NAME} -ni -fc
)
