
set(RSTKExamples
  Model1Registration.cxx
)

ADD_EXECUTABLE(model1 Model1Registration.cxx ) 
TARGET_LINK_LIBRARIES(model1 ${ITK_LIBRARIES} )