
set(RSTKApplication
  regseg.cxx
)

FIND_PACKAGE( Boost 1.49 COMPONENTS program_options REQUIRED )
INCLUDE_DIRECTORIES( ${Boost_INCLUDE_DIR} )

ADD_EXECUTABLE(regseg regseg.cxx ) 
TARGET_LINK_LIBRARIES(regseg ${ITK_LIBRARIES} ${Boost_LIBRARIES})