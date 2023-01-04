cmake_minimum_required(VERSION 3.3.0)

project("ptoy" CXX)

# Add target for the executable.
set(T "ptoy")
add_executable(${T} source/main.cpp source/geometry.cpp source/particles_system.cpp)

# C++14.
set_property(TARGET ${T} PROPERTY CXX_STANDARD 14)
set_property(TARGET ${T} PROPERTY CXX_STANDARD_REQUIRED ON)
set_property(TARGET ${T} PROPERTY CXX_EXTENSIONS OFF)

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
if (UNIX)
  target_compile_options(${T} PRIVATE -Wall -Wextra -pedantic -Wshadow)
  target_compile_options(${T} PRIVATE -march=native)
endif()

option(USE_STATIC "Enable static build" OFF)
if (USE_STATIC)
  SET(BUILD_SHARED_LIBRARIES OFF)
  set_property(TARGET ${T} PROPERTY LINK_SEARCH_START_STATIC 1)
  set_property(TARGET ${T} PROPERTY LINK_SEARCH_END_STATIC 1)
  if (UNIX)
    SET(CMAKE_FIND_LIBRARY_SUFFIXES .a)
    SET(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -static")
  else()
    SET(CMAKE_FIND_LIBRARY_SUFFIXES .lib)
  endif()
endif()

# OpenGL and GLX.
set(OpenGL_GL_PREFERENCE "GLVND")
find_package(OpenGL REQUIRED)
include_directories(${OPENGL_INCLUDE_DIRS})
target_link_libraries(${T} ${OPENGL_LIBRARIES})

find_package(SDL2 REQUIRED)
include_directories(${SDL2_INCLUDE_DIRS})
target_link_libraries(${T} SDL2::SDL2)

find_package(GLEW REQUIRED)
target_link_libraries(${T} GLEW::GLEW)

option(USE_AVX "Enable AVX vectorization" ON)
if (USE_AVX)
  add_definitions(-DUSE_AVX)
endif()
