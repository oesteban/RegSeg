set(RSTKLevelSetsTests
  MahalanobisLevelSetsTest.cxx
)

ADD_EXECUTABLE(MahalanobisLevelSetsTest MahalanobisLevelSetsTest.cxx ) 
TARGET_LINK_LIBRARIES(MahalanobisLevelSetsTest ${ITK_LIBRARIES} )

#add_library(RSTKOptimizers ${RSTKOptimizers_SRC})
#target_link_libraries(RSTKOptimizers
#  ${RSTKEnergy_LIBRARIES}
#  )
