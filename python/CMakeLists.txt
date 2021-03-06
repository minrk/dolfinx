cmake_minimum_required(VERSION 3.5.0)

PROJECT(dolfin_pybind11)

find_package(pybind11 REQUIRED CONFIG HINTS ${PYBIND11_DIR} ${PYBIND11_ROOT}
  $ENV{PYBIND11_DIR} $ENV{PYBIND11_ROOT})

find_package(DOLFIN REQUIRED)
include(${DOLFIN_USE_FILE})

# Strict compiler flags
#include(CheckCXXCompilerFlag)
#CHECK_CXX_COMPILER_FLAG("-Wall -Werror -pedantic" HAVE_PEDANTIC)
#if (HAVE_PEDANTIC)
#   set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Werror -pedantic")
#endif()

# Create the binding library
pybind11_add_module(cpp SHARED
  src/dolfin.cpp
  src/parameter.cpp
  src/common.cpp
  src/fem.cpp
  src/function.cpp
  src/generation.cpp
  src/geometry.cpp
  src/graph.cpp
  src/log.cpp
  src/math.cpp
  src/mesh.cpp
  src/io.cpp
  src/la.cpp
  src/nls.cpp
  src/refinement.cpp
  src/MPICommWrapper.cpp)

# Add DOLFIN libraries and other config
target_link_libraries(cpp PRIVATE pybind11::module dolfin)

# Add to CMake search path
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH}
  ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

# Check for petsc4py
find_package(PETSc4py REQUIRED)
if (PETSC4PY_FOUND)
  target_include_directories(cpp PRIVATE ${PETSC4PY_INCLUDE_DIRS})
  target_compile_definitions(cpp PRIVATE HAS_PYBIND11_PETSC4PY)
endif()

# Check for mpi4py
find_package(MPI4PY REQUIRED)
if (MPI4PY_FOUND)
  target_include_directories(cpp PRIVATE ${MPI4PY_INCLUDE_DIR})
  target_compile_definitions(cpp PRIVATE HAS_PYBIND11_MPI4PY)
endif()
