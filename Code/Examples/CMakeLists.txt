
set(RSTKExamples
  Model1Registration.cxx
  Model1Registration2Contour.cxx
  Model2Registration.cxx
)

ADD_EXECUTABLE(fstest FreeSurferTest.cxx ) 
TARGET_LINK_LIBRARIES(fstest ${ITK_LIBRARIES} ${Boost_LIBRARIES} )

ADD_EXECUTABLE(phantom1 phantom1.cxx ) 
TARGET_LINK_LIBRARIES(phantom1 ${ITK_LIBRARIES}  ${Boost_LIBRARIES}  )


#ADD_EXECUTABLE(model1 Model1Registration.cxx ) 
#TARGET_LINK_LIBRARIES(model1 ${ITK_LIBRARIES} )
 
# ADD_EXECUTABLE(model1_2cont Model1Registration2Contour.cxx ) 
# TARGET_LINK_LIBRARIES(model1_2cont ${ITK_LIBRARIES} )
# 
# ADD_EXECUTABLE(model2 Model2Registration.cxx ) 
# TARGET_LINK_LIBRARIES(model2 ${ITK_LIBRARIES} )