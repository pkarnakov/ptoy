cmake_minimum_required(VERSION 3.3.0)

project("ptoy" CXX)

option(USE_STATIC "Enable static build" OFF)
option(USE_AVX "Enable AVX vectorization" ON)
option(USE_WARNINGS "Enable compiler warnings" ON)

# C++14.
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

if (USE_AVX)
  add_compile_definitions(-DUSE_AVX)
  add_compile_options(-march=native)
endif()

set(T "geometry")
add_library(${T} src/geometry.cpp)

set(T "particles")
add_library(${T} src/particles.cpp)
target_link_libraries(${T} geometry)

# Executable.
set(EXE "ptoy")
add_executable(${EXE} src/main.cpp)
target_link_libraries(${EXE} geometry particles)

# Minimal executable for development.
set(EXEM "ptoymin")
add_executable(${EXEM} src/mainmin.cpp)
target_link_libraries(${EXEM} geometry particles view_gl)

# Executable.
set(T "sandbox")
add_executable(sandbox src/sandbox.cpp)
target_link_libraries(${T} geometry)

# Default build type.
set(BuildTypeValues None Debug Release RelWithDebInfo MinSizeRel)
if (NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release CACHE STRING
      "Choose the type of build, options are: ${BuildTypeValues}." FORCE)
endif ()
set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS ${BuildTypeValues})

# OpenMP.
if (NOT DEFINED OPENMP)
  set(OPENMP ON)
endif()
if (OPENMP)
  find_package(OpenMP)
  if (OPENMP_FOUND)
      set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
      set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
      set(CMAKE_EXE_LINKER_FLAGS 
          "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
  endif()
endif()

# GCC warnings and optimizations.
if (USE_WARNINGS)
  add_compile_options(-Wall -Wextra -pedantic -Wshadow)
endif()

if (USE_STATIC)
  SET(BUILD_SHARED_LIBRARIES OFF)
  set_property(TARGET ${EXE} PROPERTY LINK_SEARCH_START_STATIC 1)
  set_property(TARGET ${EXE} PROPERTY LINK_SEARCH_END_STATIC 1)
  if (UNIX)
    SET(CMAKE_FIND_LIBRARY_SUFFIXES .a)
    SET(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -static")
  else()
    SET(CMAKE_FIND_LIBRARY_SUFFIXES .lib)
  endif()
endif()

set(T "view_gl")
add_library(${T} src/view_gl.cpp)
target_link_libraries(${T} geometry particles)

# OpenGL and GLX.
set(OpenGL_GL_PREFERENCE "GLVND")
find_package(OpenGL REQUIRED)
include_directories(${OPENGL_INCLUDE_DIRS})
target_link_libraries(${EXE} ${OPENGL_LIBRARIES})
target_link_libraries(${T} ${OPENGL_LIBRARIES})

find_package(SDL2 REQUIRED)
include_directories(${SDL2_INCLUDE_DIRS})
target_link_libraries(${EXE} SDL2::SDL2)
target_link_libraries(${T} SDL2::SDL2)

find_package(GLEW REQUIRED)
target_link_libraries(${EXE} GLEW::GLEW)
target_link_libraries(${T} GLEW::GLEW)
