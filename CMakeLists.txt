cmake_minimum_required(VERSION 2.8)
project(ADDeformationModeling)

set(PETSC_DIR /home/bkhanal/Documents/softwares/petsc-3.4.3)
set(PETSC_ARCH arch-linux2-cxx-debug)


#add cmake modules path to the default one:
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/Modules/")

#Find the external library for Petsc
find_package(PETSc COMPONENTS CXX REQUIRED)
set(CMAKE_CXX_COMPILER ${PETSC_COMPILER})

#ITK
set(ITK_DIR /home/bkhanal/Documents/softwares/InsightToolkit-4.4.1/build)
FIND_PACKAGE(ITK)
IF(ITK_FOUND)
INCLUDE(${ITK_USE_FILE})
ELSE(ITK_FOUND)
MESSAGE(FATAL_ERROR
"ITK not found. Please set ITK_DIR.")
ENDIF(ITK_FOUND)

#boost
find_package(Boost COMPONENTS program_options REQUIRED)
add_subdirectory(src)


