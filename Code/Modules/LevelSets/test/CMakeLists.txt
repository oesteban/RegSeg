set(RSTKLevelSetsTests
  MahalanobisLevelSetsTest.cxx
  LevelSetsGenerateTestObjects.cxx
)

ADD_EXECUTABLE(MahalanobisLevelSetsTest MahalanobisLevelSetsTest.cxx ) 
TARGET_LINK_LIBRARIES(MahalanobisLevelSetsTest ${ITK_LIBRARIES} )

ADD_EXECUTABLE(LevelSetsGenerateTestObjects LevelSetsGenerateTestObjects.cxx ) 
TARGET_LINK_LIBRARIES(LevelSetsGenerateTestObjects ${ITK_LIBRARIES} )

#add_library(RSTKOptimizers ${RSTKOptimizers_SRC})
#target_link_libraries(RSTKOptimizers
#  ${RSTKEnergy_LIBRARIES}
#  )
