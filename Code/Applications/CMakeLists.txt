
set(RSTKApplication
  regseg.cxx
)


ADD_EXECUTABLE(regseg regseg.cxx ) 
TARGET_LINK_LIBRARIES(regseg ${ITK_LIBRARIES} ${Boost_LIBRARIES} ${JsonCpp_LIBRARY})

ADD_EXECUTABLE(bspline_field bspline_field.cxx ) 
TARGET_LINK_LIBRARIES(bspline_field ${ITK_LIBRARIES} ${Boost_LIBRARIES} ${JsonCpp_LIBRARY})