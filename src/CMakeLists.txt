#include_directories(${PETSC_TEST_SOURCE_DIR}/examples ${PETSC_INCLUDES})
include_directories(${PETSC_INCLUDES})


#link_directories(${PETSC_TEST_BINARY_DIR}/examples ${PETSC_LIBRARIES})
link_directories(${PETSC_LIBRARIES})

add_executable(PetscStokes2D PetscStokes2D.cxx)
target_link_libraries(PetscStokes2D ${PETSC_LIBRARIES})

