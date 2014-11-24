# cmake_minimum_required(VERSION 2.8)
# project(ADLemModelling)

# #Home directory
# set(MY_HOME /epi/asclepios2/bkhanal)
# #Petsc


# #add cmake modules path to the default one:
# set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${MY_HOME}/cmakeModules/petsc")

# #Find the external library for Petsc
# find_package(PETSc COMPONENTS CXX REQUIRED)
# set(CMAKE_CXX_COMPILER ${PETSC_COMPILER})


# #ITK
# set(ITK_DIR ${MY_HOME}/InsightToolkit-4.4.1/build)
# FIND_PACKAGE(ITK)
# IF(ITK_FOUND)
# INCLUDE(${ITK_USE_FILE})
# ELSE(ITK_FOUND)
# MESSAGE(FATAL_ERROR
# "ITK not found. Please set ITK_DIR.")
# ENDIF(ITK_FOUND)


# #boost
# find_package(Boost COMPONENTS program_options REQUIRED)

# add_subdirectory(src)

cmake_minimum_required(VERSION 2.8)
project(ADDeformationModeling)

set(MY_HOME /home/bkhanal)

#add cmake modules path to the default one:
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${MY_HOME}/cmakeModules/petsc")

set(PETSC_DIR ${MY_HOME}/Documents/softwares/petsc)
set(PETSC_ARCH arch-linux2-cxx-debug)

#Find the external library for Petsc
find_package(PETSc COMPONENTS CXX REQUIRED)
set(CMAKE_CXX_COMPILER ${PETSC_COMPILER})

#ITK
set(ITK_DIR ${MY_HOME}/Documents/softwares/ITK-build)
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


