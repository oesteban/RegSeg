
Project(RegSeg-Examples)


#FILE( GLOB REG_SEG_HEADERS "${SOURCE_DIR}/*.h")
#FILE( GLOB TEST_IMPLEMENTATION "${EXAMPLES_DIR}/*.cxx")
#SOURCE_GROUP( "Tests" FILES ${TEST_HEADERS} ${TEST_IMPLEMENTATION} )

INCLUDE_DIRECTORIES( ${Boost_INCLUDE_DIRS} ${SOURCE_DIR} )

SET( resamplerEx "resampler")
ADD_EXECUTABLE(${resamplerEx} "${EXAMPLES_DIR}/DeformationFieldResampler.cxx" ) 
TARGET_LINK_LIBRARIES(${resamplerEx} ${ITK_LIBRARIES} ${Boost_PROGRAM_OPTIONS_LIBRARY})

SET( dfieldEx "dfield")
ADD_EXECUTABLE(${dfieldEx} "${EXAMPLES_DIR}/SparseDeformationField.cxx" ) 
TARGET_LINK_LIBRARIES(${dfieldEx} ${ITK_LIBRARIES} ${Boost_PROGRAM_OPTIONS_LIBRARY})

SET( gradient "gradient")
ADD_EXECUTABLE(${gradient} "${EXAMPLES_DIR}/GaussianModelGradientComputation.cxx" ) 
TARGET_LINK_LIBRARIES(${gradient} ${ITK_LIBRARIES} ${Boost_PROGRAM_OPTIONS_LIBRARY})