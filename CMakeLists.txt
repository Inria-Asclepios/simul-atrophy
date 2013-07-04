cmake_minimum_required(VERSION 2.8)
project(ADDeformationModeling)

#set(PETSC_DIR /home/bkhanal/Documents/softwares/petsc-3.3-p7)
set(PETSC_DIR /home/bkhanal/Documents/softwares/petsc-3.4.1)
set(PETSC_ARCH arch-linux2-cxx-debug)
#set(LD_LIBRARY_PATH "${PETSC_DIR}/${PETSC_ARCH}/lib:${LD_LIBRARY_PATH}")


#add cmake modules path to the default one:
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/Modules/")

#Find the external library for Petsc
find_package(PETSc COMPONENTS CXX REQUIRED)
set(CMAKE_CXX_COMPILER ${PETSC_COMPILER})
#find_package(PETSc)
#find_package(PackageMultipass)

#INCLUDE_DIRECTORIES(${PETSC_DIR}/include)
#INCLUDE_DIRECTORIES(${PETSC_DIR}/${PETSC_ARCH}/include)


add_subdirectory(src)


