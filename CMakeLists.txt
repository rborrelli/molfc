cmake_minimum_required(VERSION 2.6)

#Prohibit in-source build

IF("${CMAKE_SOURCE_DIR}" STREQUAL "${CMAKE_BINARY_DIR}")
  MESSAGE(FATAL_ERROR "In-source build prohibited.")
ENDIF("${CMAKE_SOURCE_DIR}" STREQUAL "${CMAKE_BINARY_DIR}")

SET(CMAKE_Fortran_COMPILER_INIT ifort)

enable_language(Fortran)

SET(PROJECT_NAME MOLFC3)

SET(EXECUTABLE_OUTPUT_PATH ${PROJECT_BINARY_DIR}/bin)

SET(LIBRARY_OUTPUT_PATH ${PROJECT_BINARY_DIR}/lib)

SET(CMAKE_Fortran_MODULE_DIRECTORY ${PROJECT_BINARY_DIR}/mod)

#add_definitions(-DINTELF95)
#SET(CMAKE_Fortran_FLAGS "-O3 -xSSSE3 -heap-arrays 10 -openmp")
SET(CMAKE_Fortran_FLAGS "-ftrapuv -traceback -debug -O0 -g -save-temps -fpic -Wl,-no_pie")
#SET(CMAKE_Fortran_FLAGS " -debug -traceback -O0")

SET(CMAKE_Fortran_MODULE_DIRECTORY ${PROJECT_BINARY_DIR}/mod)

include_directories(src/include)

#####  SET THIS VALUE ACCORDIND TO YOUR SYSTEM ########
SET(INCLUDE_INTEL_DIRS "/opt/intel/mkl/include/intel64/lp64")
#include_directories(/opt/intel/mkl/include/intel64/lp64)
include_directories(${INCLUDE_INTEL_DIRS})

#####  SET THIS VALUE TO ON OR OFF ACCORDING TO YOUR SYSTEM ########
option(Enable_LAPACK "Compiling LAPACK and BLAS libraries" ON)
option(Enable_FFTW3 "Compiling FFTW libraries" ON)

####  Compile FFTW  #####
include(ExternalProject)

ExternalProject_Add(
  libfftw3
  SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/lib/fftw-3.3.5
  CONFIGURE_COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/lib/fftw-3.3.5/configure --prefix=${CMAKE_CURRENT_SOURCE_DIR}/build/
  PREFIX ${CMAKE_CURRENT_SOURCE_DIR}/build/lib/fftw-3.3.5
  BUILD_COMMAND ${MAKE}
  BUILD_IN_SOURCE 0
)

IF (Enable_FFTW3)
  add_library(fftw3 STATIC IMPORTED)
  set_property(TARGET fftw3 PROPERTY IMPORTED_LOCATION ${CMAKE_CURRENT_SOURCE_DIR}/build/lib/libfftw3.a)
  set_property(TARGET fftw3 PROPERTY INCLUDE_DIRECTORIES ${CMAKE_CURRENT_SOURCE_DIR}/build/lib/fftw-3.3.5/src/fftw3-build/api)
  include_directories(${CMAKE_CURRENT_SOURCE_DIR}/build/lib/fftw-3.3.5/src/fftw3-build/api)
ENDIF(Enable_FFTW3)

add_subdirectory(lib)
add_subdirectory(src)
add_subdirectory(utils)

