cmake_minimum_required(VERSION 3.3.0 FATAL_ERROR)

project("ptoy" CXX)

set(CMAKE_VERBOSE_MAKEFILE on)

# Source folders
aux_source_directory(source SRC)

# Add target
set(T "ptoy")
add_executable(${T} ${SRC})

# C++11
set_property(TARGET ${T} PROPERTY CXX_STANDARD 11)
set_property(TARGET ${T} PROPERTY CXX_STANDARD_REQUIRED ON)

# OpenMP 
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

# GCC warnings and optimizations
if (UNIX)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -pedantic")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -falign-functions=64 -march=native -frename-registers")
endif()

# Debug build
if (DEBUG) 
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g -O0")
endif()

# Static build
if (DEFINED STATIC)
  set(STATIC_BUILD ${STATIC})   
else()
  set(STATIC_BUILD OFF)
endif()
if (STATIC_BUILD)
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

# OpenGL and GLX 
set(OpenGL_GL_PREFERENCE "GLVND")
find_package(OpenGL REQUIRED)
include_directories(${OPENGL_INCLUDE_DIRS})
target_link_libraries(${T} ${OPENGL_LIBRARIES})

find_package(SDL2 REQUIRED)
include_directories(${SDL2_INCLUDE_DIRS})
target_link_libraries(${T} SDL2::SDL2)

option(USE_AVX "Enable AVX vectorization" OFF)
if (USE_AVX)
  add_definitions(-DUSE_AVX)
endif()

# Multisampling flag
if (DEFINED MULTISAMPLE)
  set(MULTISAMPLE ${MULTISAMPLE})
else()
  set(MULTISAMPLE OFF)
endif()
if (MULTISAMPLE)
  add_definitions(-D_MULTISAMPLE_)
endif()
