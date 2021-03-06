#include_directories(${PETSC_TEST_SOURCE_DIR}/examples ${PETSC_INCLUDES})
include_directories(${PETSC_INCLUDES} ${PROJECT_SOURCE_DIR}/src/includes/)

#link_directories(${PETSC_TEST_BINARY_DIR}/examples ${PETSC_LIBRARIES})
link_directories(${PETSC_LIBRARIES})

# add_executable(AdLemMain
#     AdLem3D.cxx
#     PetscAdLem3D.cxx
#     PetscAdLemTaras3D.cxx
#     PetscAdLemMain.cxx
#     )

add_executable(simul_atrophy
  PetscAdLemMain.cxx
  )
target_link_libraries(simul_atrophy ${PETSC_LIBRARIES} ${ITK_LIBRARIES})


include_directories( ${Boost_INCLUDE_DIR} )
add_executable(WarpImage
  WarpImage.cxx
  )
target_link_libraries(WarpImage ${ITK_LIBRARIES}
  ${Boost_SYSTEM_LIBRARY}
  ${Boost_PROGRAM_OPTIONS_LIBRARY}
  )

include_directories( ${Boost_INCLUDE_DIR} )
add_executable(ExtractRegionOfDtiByMask
  ExtractRegionOfDtiByMask.cxx
  )
target_link_libraries(ExtractRegionOfDtiByMask ${ITK_LIBRARIES}
  ${Boost_SYSTEM_LIBRARY}
  ${Boost_PROGRAM_OPTIONS_LIBRARY}
  )

include_directories( ${Boost_INCLUDE_DIR} )
add_executable(BinarizeThreshold
  binarizeThreshold.cxx
  )
target_link_libraries(BinarizeThreshold ${ITK_LIBRARIES}
  ${Boost_SYSTEM_LIBRARY}
  ${Boost_PROGRAM_OPTIONS_LIBRARY}
  )

add_executable(CreateImageFromRegion
  CreateImageFromRegion.cxx
  )
target_link_libraries(CreateImageFromRegion ${ITK_LIBRARIES}
  ${Boost_SYSTEM_LIBRARY}
  ${Boost_PROGRAM_OPTIONS_LIBRARY}
  )

add_executable(composeDisplacementFields
  composeDisplacementFields.cxx
  )
target_link_libraries(composeDisplacementFields ${ITK_LIBRARIES}
  ${Boost_SYSTEM_LIBRARY}
  ${Boost_PROGRAM_OPTIONS_LIBRARY}
  )

add_executable(createDiffusionTensorImage
  createDiffusionTensorImage.cxx
  )
target_link_libraries(createDiffusionTensorImage ${ITK_LIBRARIES}
  ${Boost_SYSTEM_LIBRARY}
  ${Boost_PROGRAM_OPTIONS_LIBRARY}
  )

add_executable(blurWithGaussian
  blurWithGaussian.cxx
  )
target_link_libraries(blurWithGaussian ${ITK_LIBRARIES}
  ${Boost_SYSTEM_LIBRARY}
  ${Boost_PROGRAM_OPTIONS_LIBRARY}
  )

add_executable(diffuseImage
  diffuseImage.cxx
  )
target_link_libraries(diffuseImage ${ITK_LIBRARIES}
  ${Boost_SYSTEM_LIBRARY}
  ${Boost_PROGRAM_OPTIONS_LIBRARY}
  )

add_executable(SwitchUpperLowerTriangular
  SwitchUpperLowerTriangular.cxx
  )
target_link_libraries(SwitchUpperLowerTriangular ${ITK_LIBRARIES}
  )

add_executable(createImageFromLabelImage
  createImageFromLabelImage.cxx
  )
target_link_libraries(createImageFromLabelImage ${ITK_LIBRARIES}
  ${Boost_SYSTEM_LIBRARY}
  ${Boost_PROGRAM_OPTIONS_LIBRARY}
  )

add_executable(pasteImageToBiggerImage
  pasteImageToBiggerImage.cxx
  )
target_link_libraries(pasteImageToBiggerImage ${ITK_LIBRARIES}
  )

add_executable(signedDanielssonDistance
  signedDanielssonDistance.cxx
  )
target_link_libraries(signedDanielssonDistance ${ITK_LIBRARIES}
  )

add_executable(uniformizeWithRegionalMeans
  uniformizeWithRegionalMeans.cxx
  )
target_link_libraries(uniformizeWithRegionalMeans ${ITK_LIBRARIES}
  ${Boost_SYSTEM_LIBRARY}
  ${Boost_PROGRAM_OPTIONS_LIBRARY}
  )


