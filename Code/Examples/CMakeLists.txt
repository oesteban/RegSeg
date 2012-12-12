
set(RSTKExamples
  Model1Registration.cxx
  Model1Registration2Contour.cxx
)

ADD_EXECUTABLE(model1 Model1Registration.cxx ) 
TARGET_LINK_LIBRARIES(model1 ${ITK_LIBRARIES} )


ADD_EXECUTABLE(model1_2cont Model1Registration2Contour.cxx ) 
TARGET_LINK_LIBRARIES(model1_2cont ${ITK_LIBRARIES} )