cmake_minimum_required(VERSION 2.8)
project(ADDeformationModeling)

# For cluster use argument: -DIS_CLUSTER:bool=true
set(IS_CLUSTER false CACHE BOOL "cmake run in the cluster: true or false")

# For new cluster also use argument: -DIS_NEW_CLUSTER:bool=true
set(IS_NEW_CLUSTER false CACHE BOOL "cmake run in new cluster: true or false")

MESSAGE("The directory is: ${CMAKE_SOURCE_DIR}")
#add cmake modules path to the default one:
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/Modules")

#Find the external library for Petsc
find_package(PETSc COMPONENTS CXX REQUIRED)
#find_package(PETSc COMPONENTS)
set(CMAKE_CXX_COMPILER ${PETSC_COMPILER})

#Home directory
if(IS_CLUSTER)
  if(IS_NEW_CLUSTER)
    set(MY_HOME /data/asclepios/user/bkhanal/softwares)
    set(ITK_DIR ${MY_HOME}/ITK-build)
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
    message("setting up for new cluster")

  else(IS_NEW_CLUSTER)
    set(MY_HOME /epi/asclepios2/bkhanal)
    set(ITK_DIR ${MY_HOME}/ITK-build)
    message("setting up for legacy cluster")
  endif(IS_NEW_CLUSTER)

else(IS_CLUSTER)
  set(MY_HOME /home/bkhanal)
  set(ITK_DIR ${MY_HOME}/Documents/softwares/ITK-build)
  set(CMAKE_EXPORT_COMPILE_COMMANDS on) #Will generate compile_commands.json file used by irony-mode in emacs code completion.
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
  message("setting up for delorme")
endif(IS_CLUSTER)

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


