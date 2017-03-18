cmake_minimum_required(VERSION 3.1.0 FATAL_ERROR)

# Static build flag
if (DEFINED STATIC)
  set(STATIC_BUILD ${STATIC})   
else()
  set(STATIC_BUILD OFF)
endif()

project("ptoy")

set(BINARY_NAME "ptoy")

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)  

# OpenMP support                              
find_package(OpenMP)
if (OPENMP_FOUND)
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set(CMAKE_EXE_LINKER_FLAGS 
        "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()                              

# GCC warnings and optimizations
if (UNIX)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -pedantic")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O2")
  set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -s")
endif()

# Source folders
aux_source_directory(source SOURCES_ROOT)

set(SOURCE_FILES
  ${SOURCES_ROOT}
)

# Add target
add_executable(${BINARY_NAME} ${SOURCE_FILES})
     
# C++11
set_property(TARGET ${BINARY_NAME} PROPERTY CXX_STANDARD 11)
set_property(TARGET ${BINARY_NAME} PROPERTY CXX_STANDARD_REQUIRED ON)

# Static build
if (STATIC_BUILD)
  SET(BUILD_SHARED_LIBRARIES OFF)
  set_property(TARGET ${BINARY_NAME} PROPERTY LINK_SEARCH_START_STATIC 1)
  set_property(TARGET ${BINARY_NAME} PROPERTY LINK_SEARCH_END_STATIC 1)
  if (UNIX)
    SET(CMAKE_FIND_LIBRARY_SUFFIXES .a)
    SET(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -static")
  else()
    SET(CMAKE_FIND_LIBRARY_SUFFIXES .lib)
  endif()
endif()

# Threads
find_package(Threads REQUIRED)
set(LIBRARIES ${LIBRARIES} ${CMAKE_THREAD_LIBS_INIT})

# OpenGL and Glut
find_package(OpenGL REQUIRED)
find_package(GLUT REQUIRED)
include_directories( ${OPENGL_INCLUDE_DIRS}  ${GLUT_INCLUDE_DIRS} )
set(LIBRARIES ${LIBRARIES} ${OPENGL_LIBRARIES} ${GLUT_LIBRARY})

target_link_libraries(${BINARY_NAME} ${LIBRARIES})
