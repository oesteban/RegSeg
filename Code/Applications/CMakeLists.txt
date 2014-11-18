
set(RSTKApplication
  regseg.cxx
)


ADD_EXECUTABLE(regseg regseg.cxx ) 
TARGET_LINK_LIBRARIES(regseg ${ITK_LIBRARIES} ${Boost_LIBRARIES} ${JsonCpp_LIBRARY})

ADD_EXECUTABLE(bspline_field bspline_field.cxx ) 
TARGET_LINK_LIBRARIES(bspline_field ${ITK_LIBRARIES} ${Boost_LIBRARIES} ${JsonCpp_LIBRARY})

ADD_EXECUTABLE(bspline_invert bspline_invert.cxx ) 
TARGET_LINK_LIBRARIES(bspline_invert ${ITK_LIBRARIES} ${Boost_LIBRARIES} ${JsonCpp_LIBRARY})

ADD_EXECUTABLE(warp_image warp_image.cxx ) 
TARGET_LINK_LIBRARIES(warp_image ${ITK_LIBRARIES} ${Boost_LIBRARIES})

ADD_EXECUTABLE(regridavg regridavg.cxx ) 
TARGET_LINK_LIBRARIES(regridavg ${ITK_LIBRARIES} ${Boost_LIBRARIES})

INSTALL(TARGETS bspline_field bspline_invert regseg warp_image regridavg
  RUNTIME DESTINATION bin
)


ADD_EXECUTABLE(slice_contours slice_contours.cxx ) 


if(VTK_LIBRARIES)
	TARGET_LINK_LIBRARIES(slice_contours ${ITK_LIBRARIES} ${Boost_LIBRARIES} ${VTK_LIBRARIES})
else()
	TARGET_LINK_LIBRARIES(slice_contours ${ITK_LIBRARIES} ${Boost_LIBRARIES} vtkHybrid vtkWidgets)
endif()