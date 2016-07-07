cmake_minimum_required(VERSION 2.8)
project(simul@atrophy)

# add cmake modules path to the default one:
MESSAGE("The directory is: ${CMAKE_SOURCE_DIR}")
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/Modules")

# Use cxx11 or not; If yes run cmake as: cmake -DUSE_CXX11:bool=true .. 
set(USE_CXX11 false CACHE BOOL "use cxx11: true or false")

if(USE_CXX11)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
  # This part taken from: http://stackoverflow.com/questions/10851247/how-to-activate-c-11-in-cmake
  include(CheckCXXCompilerFlag)
  CHECK_CXX_COMPILER_FLAG("-std=c++11" COMPILER_SUPPORTS_CXX11)
  CHECK_CXX_COMPILER_FLAG("-std=c++0x" COMPILER_SUPPORTS_CXX0X)
  if(COMPILER_SUPPORTS_CXX11)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
  elseif(COMPILER_SUPPORTS_CXX0X)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++0x")
  else()
    message(STATUS "The compiler ${CMAKE_CXX_COMPILER} has no C++11 support. Please use a different C++ compiler.")
  endif()
endif(USE_CXX11)


# Find the external library for PETSc
# PETSc must have been configured with --with-clanguage=cxx option!
find_package(PETSc COMPONENTS CXX REQUIRED)
set(CMAKE_CXX_COMPILER ${PETSC_COMPILER})

# Generate compile_commands.json file used by irony-mode in emacs code completion.
set(CMAKE_EXPORT_COMPILE_COMMANDS on)


# Find ITK
FIND_PACKAGE(ITK)
IF(ITK_FOUND)
  INCLUDE(${ITK_USE_FILE})
ELSE(ITK_FOUND)
  MESSAGE(FATAL_ERROR
    "ITK not found. Please set ITK_DIR.")
ENDIF(ITK_FOUND)

# Find boost
find_package(Boost COMPONENTS program_options REQUIRED)

add_subdirectory(src)


