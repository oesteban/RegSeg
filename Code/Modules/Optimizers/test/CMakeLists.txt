set(RSTKCoreTests
  GradientDescentLevelSetsOptimizerTest.cxx
)

ADD_EXECUTABLE(GradientDescentLevelSetsOptimizerTest GradientDescentLevelSetsOptimizerTest.cxx ) 
TARGET_LINK_LIBRARIES(GradientDescentLevelSetsOptimizerTest ${ITK_LIBRARIES} )

#add_library(RSTKOptimizers ${RSTKOptimizers_SRC})
#target_link_libraries(RSTKOptimizers
#  ${RSTKEnergy_LIBRARIES}
#  )