set(RSTKCoreTests
  ContourDisplacementFieldTest.cxx
)

ADD_EXECUTABLE(ContourDisplacementFieldTest ContourDisplacementFieldTest.cxx ) 
TARGET_LINK_LIBRARIES(ContourDisplacementFieldTest ${ITK_LIBRARIES} )

#add_library(RSTKOptimizers ${RSTKOptimizers_SRC})
#target_link_libraries(RSTKOptimizers
#  ${RSTKEnergy_LIBRARIES}
#  )