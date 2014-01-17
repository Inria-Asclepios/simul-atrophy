#include_directories(${PETSC_TEST_SOURCE_DIR}/examples ${PETSC_INCLUDES})
include_directories(${PETSC_INCLUDES} ${PROJECT_SOURCE_DIR}/src/includes/)

#link_directories(${PETSC_TEST_BINARY_DIR}/examples ${PETSC_LIBRARIES})
link_directories(${PETSC_LIBRARIES})

add_executable(AdLemMain
    #PetscAdLemMain.cxx
    #AdLem3D.cxx
    #PetscAdLem3D.cxx
    #PetscAdLemTaras3D.cxx
    HelmHoltzDecomposer.cxx
    HhDecomposeMain.cxx
    PoissonProblem.cxx
    PoissonSolver.cxx

)


target_link_libraries(AdLemMain ${PETSC_LIBRARIES} ${ITK_LIBRARIES})
