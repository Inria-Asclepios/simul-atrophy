cmake_minimum_required(VERSION 2.8)
project(ADDeformationModeling)

# For cluster use argument: -DIS_CLUSTER:bool=true
set(IS_CLUSTER false CACHE BOOL "cmake run in the cluster: true or false")

#Home directory
if(IS_CLUSTER)
  set(MY_HOME /epi/asclepios2/bkhanal)
  set(ITK_DIR ${MY_HOME}/ITK-build)
  message("setting up for cluster")
else(IS_CLUSTER)
  set(MY_HOME /home/bkhanal)
  set(ITK_DIR ${MY_HOME}/Documents/softwares/ITK-build)
  message("setting up for delorme")
endif(IS_CLUSTER)

#add cmake modules path to the default one:
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${MY_HOME}/cmakeModules/petsc")

#Find the external library for Petsc
find_package(PETSc COMPONENTS CXX REQUIRED)
set(CMAKE_CXX_COMPILER ${PETSC_COMPILER})

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


