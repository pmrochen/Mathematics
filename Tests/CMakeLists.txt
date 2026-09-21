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

target_include_directories(${MATHEMATICS_TESTS_NAME} PRIVATE
	#${CMAKE_CURRENT_SOURCE_DIR}/../Include # Mathematics headers
	${CMAKE_CURRENT_BINARY_DIR}/downloaded # doctest.h
)

target_sources(${MATHEMATICS_TESTS_NAME} PRIVATE
	Main.cpp
)

add_test(
	NAME "${MATHEMATICS_TESTS_NAME}"
	COMMAND ${MATHEMATICS_TESTS_NAME} -ni -fc
)	
