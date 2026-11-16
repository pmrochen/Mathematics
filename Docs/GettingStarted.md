# Getting Started

[← Back to the README](../README.md)

- [Getting the sources](#getting-the-sources)
- [Using the library from CMake](#using-the-library-from-cmake)
- [Using the library without CMake](#using-the-library-without-cmake)
- [Build options](#build-options)
- [Preprocessor definitions](#preprocessor-definitions)
- [Headers](#headers)
- [A first program](#a-first-program)
- [Building and running the tests](#building-and-running-the-tests)
- [Repository layout](#repository-layout)

## Getting the sources

The library depends on the [Simd](https://github.com/pmrochen/Simd) library, which is included as a Git submodule.
Clone both in one step:

```bash
git clone --recurse-submodules https://github.com/pmrochen/Mathematics.git
```

If the repository was already cloned without submodules:

```bash
git submodule update --init --recursive
```

## Using the library from CMake

Mathematics is an `INTERFACE` (header-only) target named `Mathematics`, with the alias `Core::Mathematics`. Linking it
gives the consumer:

- the `Include` directory (and the Simd include directory when SIMD is enabled),
- the `cxx_std_20` compile feature,
- the `MATHEMATICS_*` compile definitions that match the chosen [build options](#build-options),
- `-fno-strict-aliasing` on GCC, Clang, AppleClang and IntelLLVM (see
  [Conventions](Conventions.md#strict-aliasing)).

### As a subdirectory or submodule

```cmake
add_subdirectory(External/Mathematics)

add_executable(MyApplication Main.cpp)
target_link_libraries(MyApplication PRIVATE Core::Mathematics)
```

When Mathematics is not the top-level project, the tests are disabled by default.

### With FetchContent

```cmake
include(FetchContent)

FetchContent_Declare(Mathematics
    GIT_REPOSITORY https://github.com/pmrochen/Mathematics.git
    GIT_TAG main
)
FetchContent_MakeAvailable(Mathematics)

target_link_libraries(MyApplication PRIVATE Core::Mathematics)
```

`FetchContent` checks out submodules by default, so Simd is fetched as well.

### Setting options from the parent project

Options are ordinary cache variables and can be set before the library is added:

```cmake
set(MATHEMATICS_DOUBLE_PRECISION ON CACHE BOOL "" FORCE)
add_subdirectory(External/Mathematics)
```

If the parent project already defines a `Core::Simd` target, Mathematics uses it instead of adding its bundled copy.

### Choosing the instruction set

The library does not force an instruction set on its consumers. With SIMD enabled on x86-64, compile the code that
includes the library with at least SSE4.1:

```cmake
if(MSVC)
    target_compile_options(MyApplication PRIVATE /arch:AVX2)
else()
    target_compile_options(MyApplication PRIVATE -msse4.2 -mavx2)
endif()
```

On other architectures (for example arm64) the SIMD specializations are not active yet and the portable scalar code is
used automatically.

## Using the library without CMake

Because the library is header-only, any build system works:

1. Add `Mathematics/Include` to the include path.
2. For SIMD support also add `Mathematics/Simd/Include` and define `MATHEMATICS_HAS_SIMD=1`.
3. Compile as C++20.
4. On GCC and Clang add `-fno-strict-aliasing`.

```bash
g++ -std=c++20 -O2 -fno-strict-aliasing -msse4.2 -DMATHEMATICS_HAS_SIMD=1 -IMathematics/Include -IMathematics/Simd/Include Main.cpp
```

## Build options

| CMake option | Default | Description |
| --- | --- | --- |
| `MATHEMATICS_WITH_TESTS` | `ON` when Mathematics is the top-level project, otherwise `OFF` | Adds the `Tests` directory and the `MathematicsTests` executable. |
| `MATHEMATICS_DOUBLE_PRECISION` | `OFF` | Selects `double` for the unsuffixed type aliases. Explicit aliases such as `Vector3f` and `Vector3d` are always available. |
| `MATHEMATICS_ENABLE_SIMD` | `ON` | Links `Core::Simd` and defines `MATHEMATICS_HAS_SIMD=1`. |
| `MATHEMATICS_ENABLE_SSE4` | `OFF` | Adds `-msse4.2` (`/arch:SSE4` on MSVC) to the test build. |
| `MATHEMATICS_ENABLE_AVX` | `OFF` | Adds `-msse4.2 -mavx` (`/arch:AVX`) to the test build. |
| `MATHEMATICS_ENABLE_AVX2` | `OFF` | Adds `-msse4.2 -mavx2` (`/arch:AVX2`) to the test build. |
| `MATHEMATICS_ENABLE_QUICKHULL` | `ON` | Reserved for convex hull generation. The integration is not available yet and the option has no effect. |

The three instruction set options are mutually exclusive in effect: AVX2 takes precedence over AVX, which takes
precedence over SSE4. They apply only when both `MATHEMATICS_WITH_TESTS` and `MATHEMATICS_ENABLE_SIMD` are on.

## Preprocessor definitions

The CMake target defines these for every consumer. When building without CMake they may be defined by hand; all of them
default to `0` when left undefined, except where noted.

| Definition | Meaning |
| --- | --- |
| `MATHEMATICS_DOUBLE_PRECISION` | `1` makes `Vector3`, `Matrix3`, `Sphere`, ... aliases of the `double` instantiations. |
| `MATHEMATICS_HAS_SIMD` | `1` includes `<Simd/Simd.hpp>` and enables the SIMD specializations where the target supports them. |
| `MATHEMATICS_HAS_QUICKHULL` | Reserved. Always `0`. |
| `MATHEMATICS_FAST_NORMALIZE` | Defaults to `1`. The SIMD `normalize()` implementations use an approximate reciprocal square root. Define as `0` before including any header to use an exact square root and division instead. |

## Headers

`<Mathematics/Mathematics.hpp>` includes all the commonly used headers. Individual headers can be included instead to
reduce compile times; each one includes what it needs.

| Directory | Headers |
| --- | --- |
| `Mathematics/` | `Mathematics.hpp`, `Concepts.hpp`, `Constants.hpp`, `Scalar.hpp`, `Interval.hpp`, `Random.hpp` |
| `Mathematics/Algebra/` | `Vector2.hpp`, `Vector3.hpp`, `Vector4.hpp`, `Matrix2.hpp`, `Matrix3.hpp`, `Matrix4.hpp`, `Quaternion.hpp`, `LU.inl` |
| `Mathematics/Transform/` | `AffineTransform.hpp`, `Euler.hpp`, `EulerOrder.hpp`, `YawPitchRoll.hpp` |
| `Mathematics/Geometry/` | `Axis.hpp`, the primitive and mesh headers, `Bezier.inl`, `Hermite.inl`, `Triangulation.inl` |
| `Mathematics/Analysis/` | `Containment.inl`, `Distances.inl`, `Intersections.inl`, `Perlin.inl`, `Romberg.inl` |

Two headers are **not** pulled in by the umbrella header and must be included explicitly when needed:

```cpp
#include <Mathematics/Analysis/Perlin.inl>      // mathematics::perlin::noise
#include <Mathematics/Geometry/Hermite.inl>     // mathematics::hermite::basis
```

## A first program

```cpp
#include <Mathematics/Mathematics.hpp>
#include <iostream>

using namespace mathematics;

int main()
{
    Vector3 v(1.f, 2.f, 3.f);
    std::cout << length(v) << "\n";                         // 3.74166
    std::cout << normalize(v) << "\n";                      // 0.267261 0.534522 0.801784

    Matrix3 rotation = Matrix3::makeRotation(Axis::Z, radians(90.f));
    std::cout << Vector3::UNIT_X*rotation << "\n";          // ~0 1 0

    AxisAlignedBox box(Vector3(-1.f), Vector3(1.f));
    Sphere sphere(Vector3::ZERO, 1.f);
    std::cout << box.intersects(sphere) << "\n";            // 1
    std::cout << box.getVolume() << "\n";                   // 8
}
```

## Building and running the tests

The tests use [doctest](https://github.com/doctest/doctest), which CMake downloads at configure time, so the first
configuration needs network access.

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DMATHEMATICS_ENABLE_AVX2=ON
```

```bash
cmake --build build --target MathematicsTests --parallel
```

```bash
ctest --test-dir build --output-on-failure
```

`MathematicsTests` is excluded from the default `all` target, so it has to be named explicitly when building.

Useful variations:

| Goal | Configure arguments |
| --- | --- |
| SSE4.2 only | `-DMATHEMATICS_ENABLE_SSE4=ON` |
| Portable scalar code | `-DMATHEMATICS_ENABLE_SIMD=OFF` |
| Double precision | `-DMATHEMATICS_DOUBLE_PRECISION=ON` |
| arm64 and other non-x86 hosts | no instruction set options |

### Configure presets

`CMakePresets.json` defines Ninja-based configure presets that write to `out/build/<preset>`:

| Platform | Presets |
| --- | --- |
| Windows, MSVC | `windows-msvc-x64-debug`, `windows-msvc-x64-release`, `windows-msvc-arm64-debug`, `windows-msvc-arm64-release` |
| Windows, clang-cl | `windows-clang-x64-debug`, `windows-clang-x64-release` |
| Linux, GCC | `linux-gcc-x64-debug`, `linux-gcc-x64-release`, `linux-gcc-arm64-debug`, `linux-gcc-arm64-release` |
| Linux, Clang | `linux-clang-x64-debug`, `linux-clang-x64-release`, `linux-clang-arm64-debug`, `linux-clang-arm64-release` |
| WebAssembly | `emscripten-debug`, `emscripten-release` (requires the `EMSDK` environment variable) |

The x64 presets enable AVX2, and the Emscripten presets enable SSE4.2 emulation.

```bash
cmake --preset linux-gcc-x64-release
```

```bash
cmake --build out/build/linux-gcc-x64-release --target MathematicsTests
```

```bash
ctest --preset linux-gcc-x64-release
```

### Continuous integration

The [Tests workflow](../.github/workflows/tests.yml) builds and runs the test suite on every push and pull request
with the following configurations:

| Platform | Compiler | Configuration |
| --- | --- | --- |
| Ubuntu 24.04 | GCC 14 | AVX2, Release |
| Ubuntu 24.04 | Clang 18 | AVX2, Debug |
| Ubuntu 24.04 | Clang 18 | SSE4.2, Release |
| Ubuntu 24.04 | GCC 14 | No SIMD, double precision, Release |
| macOS (arm64) | AppleClang | Portable code, Release |
| Windows | MSVC | AVX2, Release |

## Repository layout

```
Mathematics/
├── CMakeLists.txt          Interface target, options
├── CMakePresets.json       Configure and test presets
├── Include/Mathematics/    Public headers
│   ├── Algebra/            Vectors, matrices, quaternions, LU decomposition
│   ├── Transform/          Affine transformation, Euler angles
│   ├── Geometry/           Primitives, curves, meshes
│   └── Analysis/           Low-level queries, noise, integration
├── Simd/                   Simd library (Git submodule)
├── Tests/                  doctest unit tests
└── Docs/                   This documentation
```
