cmake_minimum_required(VERSION 3.10.0)

project("ptoy" CXX)

option(USE_AVX "Enable AVX vectorization" ON)
option(USE_OPENMP "Enable OpenMP" ON)
option(USE_EXE "Build executable " ON)
option(USE_SANDBOX "Build sandbox executable for development" OFF)
option(USE_WASM "Build WebAssembly" OFF)
option(USE_BACKEND_SDL "Enable SDL2 backend" ON)
option(USE_BACKEND_TEXT "Enable text backend" OFF)
option(USE_WARNINGS "Enable compiler warnings" ON)

# Default build type.
set(BuildTypeValues None Debug Release RelWithDebInfo MinSizeRel)
if (NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release CACHE STRING
      "Choose the type of build, options are: ${BuildTypeValues}." FORCE)
endif ()
set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS ${BuildTypeValues})

# C++14.
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

set(T "geometry")
add_library(${T} src/geometry.cpp)

set(T "particles")
add_library(${T} src/particles.cpp)
target_link_libraries(${T} geometry)
if (USE_AVX)
  target_compile_definitions(${T} PRIVATE USE_AVX=1)
  target_compile_options(${T} PRIVATE -march=native)
else()
  target_compile_definitions(${T} PRIVATE USE_AVX=0)
endif()

set(T "control")
add_library(${T} src/control.cpp)

# Executable.
if (USE_EXE)
  set(EXE "ptoy")
  add_executable(${EXE} src/main.cpp)
  target_link_libraries(${EXE} geometry particles)
  if (USE_BACKEND_TEXT)
    target_compile_definitions(${EXE} PRIVATE USE_BACKEND_TEXT=1)
  else()
    target_compile_definitions(${EXE} PRIVATE USE_BACKEND_TEXT=0)
  endif()
  if (USE_BACKEND_SDL)
    target_link_libraries(${EXE} view_gl control)
    target_compile_definitions(${EXE} PRIVATE USE_BACKEND_SDL=1)
  else()
    target_compile_definitions(${EXE} PRIVATE USE_BACKEND_SDL=0)
  endif()
endif()

# Sandbox executable.
if (USE_SANDBOX)
  set(T "sandbox")
  add_executable(sandbox src/sandbox.cpp)
  target_link_libraries(${T} geometry)
endif()

# OpenMP.
if (USE_OPENMP)
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


if (USE_BACKEND_SDL)
  set(T "view_gl")
  add_library(${T} src/view_gl.cpp)
  target_link_libraries(${T} geometry particles)

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
endif()

if (USE_WASM)
  set(T ptoy)
  add_link_options(-sALLOW_MEMORY_GROWTH=1)
  add_executable(${T} src/ptoy_wasm.cpp)
  target_link_libraries(${T} geometry particles control)
  target_link_options(${T} PRIVATE "-sEXPORTED_FUNCTIONS=[\
    '_main', '_malloc' \
    , '_GetParticles', '_GetPortals', '_GetBonds', '_GetFrozen' \
    , '_SendKeyDown', '_SendMouseMotion', '_SendMouseDown', '_SendMouseUp' \
    , '_SetControlDebug', '_Init', '_SetPause', '_GetMouseMode' \
    , '_GetGravity', '_SetGravity', '_SetGravityVect' \
  ]")
  target_link_options(${T} PRIVATE "-sEXPORTED_RUNTIME_METHODS=['cwrap', 'ccall']")
  target_compile_options(${T} PRIVATE -fexceptions)
  target_link_options(${T} PRIVATE -fexceptions)
  set_target_properties(${T} PROPERTIES OUTPUT_NAME ptoy)
  configure_file(src/ptoy.html ptoy.html COPYONLY)
  configure_file(src/ptoy.css ptoy.css COPYONLY)
  configure_file(src/ptoy_inc.js ptoy_inc.js COPYONLY)
  configure_file(src/libs/lz-string.js libs/lz-string.js COPYONLY)
  configure_file(assets/favicon.png favicon.png COPYONLY)
  configure_file(src/libs/github_buttons.js libs/github_buttons.js COPYONLY)
endif()

