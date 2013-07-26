cmake_minimum_required(VERSION 2.8)
project(ADDeformationModeling)

#set(PETSC_DIR /home/bkhanal/Documents/softwares/petsc-3.3-p7)
set(PETSC_DIR /home/bkhanal/Documents/softwares/petsc-3.4.2)
set(PETSC_ARCH arch-linux2-cxx-debug)


#add cmake modules path to the default one:
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/Modules/")

#Find the external library for Petsc
find_package(PETSc COMPONENTS CXX REQUIRED)
set(CMAKE_CXX_COMPILER ${PETSC_COMPILER})


add_subdirectory(src)


