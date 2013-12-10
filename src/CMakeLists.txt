#include_directories(${PETSC_TEST_SOURCE_DIR}/examples ${PETSC_INCLUDES})
include_directories(${PETSC_INCLUDES} ${PROJECT_SOURCE_DIR}/src/includes/)

#link_directories(${PETSC_TEST_BINARY_DIR}/examples ${PETSC_LIBRARIES})
link_directories(${PETSC_LIBRARIES})

#add_executable(PetscStokes2D PetscStokes2D.cxx)
#add_executable(PetscStokes2D PetscStokes2DTest.cxx)
#target_link_libraries(PetscStokes2D ${PETSC_LIBRARIES})

add_executable(AdLemMain
    PetscAdLemMain.cxx
    #AdLem2D.cxx
    #PetscAdLem2D.cxx
    #PetscAdLemTaras2D.cxx
    AdLem3D.cxx
    PetscAdLem3D.cxx
    PetscAdLemTaras3D.cxx
    #PetscAdLemTaras3D_computeMatrixRhsDrivenCavityTest.cxx
    #PetscAdLemTaras3D_computeMatrixRhsDrivenCavityTest1.cxx
)


target_link_libraries(AdLemMain ${PETSC_LIBRARIES} ${ITK_LIBRARIES})
