
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

ADD_EXECUTABLE(binarize_surfaces binarize_surfaces.cxx )
TARGET_LINK_LIBRARIES(binarize_surfaces ${ITK_LIBRARIES} ${Boost_LIBRARIES})

ADD_EXECUTABLE(regridavg regridavg.cxx ) 
TARGET_LINK_LIBRARIES(regridavg ${ITK_LIBRARIES} ${Boost_LIBRARIES})

find_package(VTK COMPONENTS vtkRenderingOpenGL vtkRenderingFreeType vtkRenderingImage vtkIOLegacy NO_MODULE)

ADD_EXECUTABLE(slice_contours slice_contours.cxx )
set_property(TARGET slice_contours APPEND
  PROPERTY COMPILE_DEFINITIONS "${VTK_DEFINITIONS}")
TARGET_LINK_LIBRARIES(slice_contours ${ITK_LIBRARIES} ${Boost_LIBRARIES} ${VTK_LIBRARIES})
set_property(TARGET slice_contours APPEND
  PROPERTY INCLUDE_DIRECTORIES ${VTK_INCLUDE_DIRS})

INSTALL(TARGETS bspline_field bspline_invert regseg warp_image binarize_surfaces regridavg slice_contours
  RUNTIME DESTINATION bin
)
