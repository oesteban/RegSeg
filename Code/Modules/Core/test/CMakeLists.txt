set(RSTKCoreTests
  ContourDisplacementFieldTest.cxx
  IDWMultivariateInterpolatorTest.cxx
)

ADD_EXECUTABLE(itkTriangleMeshToBinaryImageFilterTest itkTriangleMeshToBinaryImageFilterTest.cxx ) 
TARGET_LINK_LIBRARIES(itkTriangleMeshToBinaryImageFilterTest ${ITK_LIBRARIES} )

ADD_EXECUTABLE(ContourDisplacementFieldTest ContourDisplacementFieldTest.cxx ) 
TARGET_LINK_LIBRARIES(ContourDisplacementFieldTest ${ITK_LIBRARIES} )

ADD_EXECUTABLE(IDWMultivariateInterpolatorTest IDWMultivariateInterpolatorTest.cxx ) 
TARGET_LINK_LIBRARIES(IDWMultivariateInterpolatorTest ${ITK_LIBRARIES} )

#add_library(RSTKOptimizers ${RSTKOptimizers_SRC})
#target_link_libraries(RSTKOptimizers
#  ${RSTKEnergy_LIBRARIES}
#  )