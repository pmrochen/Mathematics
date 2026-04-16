cmake_minimum_required(VERSION 3.25)

project(Mathematics LANGUAGES CXX)

if(NOT TARGET "Simd")
    add_subdirectory(Simd)
endif()

set(MATHEMATICS_LIB_NAME "Mathematics")
add_library(${MATHEMATICS_LIB_NAME} INTERFACE)

target_include_directories(${MATHEMATICS_LIB_NAME}
	INTERFACE
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
        $<INSTALL_INTERFACE:.>
)

target_compile_features(${MATHEMATICS_LIB_NAME} 
    INTERFACE 
        cxx_std_20
)
