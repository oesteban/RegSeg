
set(RSTKApplication
  regseg.cxx
)

FIND_PACKAGE( Boost 1.49 COMPONENTS program_options filesystem system REQUIRED )
INCLUDE_DIRECTORIES( ${Boost_INCLUDE_DIR} )

FIND_PATH(JsonCpp_INCLUDE_DIR     
  NAMES 
        jsoncpp/json/json.h   
  PATHS
        /usr/local/include
        /usr/include
)
INCLUDE_DIRECTORIES(${JsonCpp_INCLUDE_DIR})

find_library(JsonCpp_LIBRARY 
    NAMES 
        jsoncpp
        libjsoncpp
    PATHS
       /usr/lib 
       /usr/local/lib
)

ADD_EXECUTABLE(regseg regseg.cxx ) 
TARGET_LINK_LIBRARIES(regseg ${ITK_LIBRARIES} ${Boost_LIBRARIES} ${JsonCpp_LIBRARY})