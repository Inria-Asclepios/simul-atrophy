#include_directories(${PETSC_TEST_SOURCE_DIR}/examples ${PETSC_INCLUDES})
include_directories(${PETSC_INCLUDES} ${PROJECT_SOURCE_DIR}/src/includes/)

#link_directories(${PETSC_TEST_BINARY_DIR}/examples ${PETSC_LIBRARIES})
link_directories(${PETSC_LIBRARIES})

add_executable(AdLemMain
    AdLem3D.cxx
    PetscAdLem3D.cxx
    PetscAdLemTaras3D.cxx
    PetscAdLemMain.cxx
)
target_link_libraries(AdLemMain ${PETSC_LIBRARIES} ${ITK_LIBRARIES})


add_executable(HhDecompose
    PoissonSolver.cxx
    PoissonProblem.cxx
    HelmHoltzDecomposer.cxx
    HhDecomposeMain.cxx
)
target_link_libraries(HhDecompose ${PETSC_LIBRARIES} ${ITK_LIBRARIES})


add_executable(UpdateMaskAndAtrophy
    UpdateMaskAndAtrophyMain.cxx
)
target_link_libraries(UpdateMaskAndAtrophy ${ITK_LIBRARIES})


include_directories( ${Boost_INCLUDE_DIR} )
add_executable(ModifyAtrophy
    testsMain.cxx
)
target_link_libraries(ModifyAtrophy ${ITK_LIBRARIES}
    ${Boost_SYSTEM_LIBRARY}
    ${Boost_PROGRAM_OPTIONS_LIBRARY}
)


include_directories( ${Boost_INCLUDE_DIR} )
add_executable(WarpImage
    WarpImage.cxx
)
target_link_libraries(WarpImage ${ITK_LIBRARIES}
    ${Boost_SYSTEM_LIBRARY}
    ${Boost_PROGRAM_OPTIONS_LIBRARY}
)

