cmake_minimum_required(VERSION 3.25)

# Project metadata
project(Mathematics 
    VERSION 1.0.0 
    LANGUAGES CXX 
    DESCRIPTION "Modularized Linear Algebra and Geometry library"
)

# Check for SIMD submodule presence
if(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/Simd/Simd.hpp)
    add_subdirectory(Simd)
endif()

set(MATHEMATICS_LIB_NAME "Mathematics")

# Define header-only library using INTERFACE
add_library(${MATHEMATICS_LIB_NAME} INTERFACE)

# Enforce C++20 standard
target_compile_features(${MATHEMATICS_LIB_NAME} 
    INTERFACE 
        cxx_std_20
)

# Configure include directories for builds and installations
target_include_directories(${MATHEMATICS_LIB_NAME}
    INTERFACE
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
        $<INSTALL_INTERFACE:include>
)

# Organization for IDEs (Group files into virtual folders)
file(GLOB_RECURSE ALL_HEADERS "include/*.hpp" "include/*.inl")
source_group(TREE "${CMAKE_CURRENT_SOURCE_DIR}/include" FILES ${ALL_HEADERS})

# Metadata for tracking
set_target_properties(${MATHEMATICS_LIB_NAME} PROPERTIES
    INTERFACE_MAINTAINER "hypernova-developer"
)
