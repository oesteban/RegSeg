#set(RSTKFunctionalTests
#  MahalanobisFunctionalTest.cxx
#  FunctionalGenerateTestObjects.cxx
#  FunctionalBaseTest.cxx
#)

#ADD_EXECUTABLE(FunctionalBaseTest FunctionalBaseTest.cxx )
#TARGET_LINK_LIBRARIES(  FunctionalBaseTest gtest ${ITK_LIBRARIES} )
#ADD_TEST( NAME FunctionalBaseTest COMMAND FunctionalBaseTest WORKING_DIRECTORY "/home/oesteban/workspace/ACWE-Reg/Data/Ellipse")
#
#ADD_EXECUTABLE(MahalanobisFunctionalTest MahalanobisFunctionalTest.cxx ) 
#TARGET_LINK_LIBRARIES(MahalanobisFunctionalTest ${ITK_LIBRARIES} )
#
#ADD_EXECUTABLE(FunctionalGenerateTestObjects FunctionalGenerateTestObjects.cxx ) 
#TARGET_LINK_LIBRARIES(FunctionalGenerateTestObjects ${ITK_LIBRARIES} )

#add_library(RSTKOptimizers ${RSTKOptimizers_SRC})
#target_link_libraries(RSTKOptimizers
#  ${RSTKEnergy_LIBRARIES}
#  )
