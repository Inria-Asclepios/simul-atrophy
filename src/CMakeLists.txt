#include_directories(${PETSC_TEST_SOURCE_DIR}/examples ${PETSC_INCLUDES})
include_directories(${PETSC_INCLUDES} ${PROJECT_SOURCE_DIR}/src/includes/)

#link_directories(${PETSC_TEST_BINARY_DIR}/examples ${PETSC_LIBRARIES})
link_directories(${PETSC_LIBRARIES})

#add_executable(PetscStokes2D PetscStokes2D.cxx)
#add_executable(PetscStokes2D PetscStokes2DTest.cxx)
#target_link_libraries(PetscStokes2D ${PETSC_LIBRARIES})

add_executable(AdLem2Dmain
    PetscAdLem2DMain.cxx
    AdLem2D.cxx
    PetscAdLem2D.cxx
    PetscAdLemTaras2D.cxx
)
target_link_libraries(AdLem2Dmain ${PETSC_LIBRARIES})
