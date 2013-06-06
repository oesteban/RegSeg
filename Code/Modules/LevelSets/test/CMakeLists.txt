#set(RSTKLevelSetsTests
#  MahalanobisLevelSetsTest.cxx
#  LevelSetsGenerateTestObjects.cxx
#  LevelSetsBaseTest.cxx
#)

#ADD_EXECUTABLE(LevelSetsBaseTest LevelSetsBaseTest.cxx )
#TARGET_LINK_LIBRARIES(  LevelSetsBaseTest gtest ${ITK_LIBRARIES} )
#ADD_TEST( NAME LevelSetsBaseTest COMMAND LevelSetsBaseTest WORKING_DIRECTORY "/home/oesteban/workspace/ACWE-Reg/Data/Ellipse")
#
#ADD_EXECUTABLE(MahalanobisLevelSetsTest MahalanobisLevelSetsTest.cxx ) 
#TARGET_LINK_LIBRARIES(MahalanobisLevelSetsTest ${ITK_LIBRARIES} )
#
#ADD_EXECUTABLE(LevelSetsGenerateTestObjects LevelSetsGenerateTestObjects.cxx ) 
#TARGET_LINK_LIBRARIES(LevelSetsGenerateTestObjects ${ITK_LIBRARIES} )

#add_library(RSTKOptimizers ${RSTKOptimizers_SRC})
#target_link_libraries(RSTKOptimizers
#  ${RSTKEnergy_LIBRARIES}
#  )
