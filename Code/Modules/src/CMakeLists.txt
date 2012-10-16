cmake_minimum_required(VERSION 2.8)

PROJECT( RegSeg )
OPTION( BUILD_SHARED_LIBS OFF )

OPTION(BUILD_EXAMPLES "Build the Examples directory." ON )

#-----------------------------------------------------------------------------
# Set a default build type if none was specified
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  message(STATUS "Setting build type to 'Release' as none was specified.")
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Choose the type of build." FORCE)
  # Set the possible values of build type for cmake-gui
  set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug" "Release" "MinSizeRel" "RelWithDebInfo")
endif()

# Setup build locations.
if(NOT CMAKE_RUNTIME_OUTPUT_DIRECTORY)
  set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/bin)
endif()
if(NOT CMAKE_LIBRARY_OUTPUT_DIRECTORY)
  set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/lib)
endif()
if(NOT CMAKE_ARCHIVE_OUTPUT_DIRECTORY)
  set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/lib)
endif()


IF (WIN32)
	SET( Boost_NO_SYSTEM_PATHS ON )
	SET ( BOOST_ROOT "C:/Boost" )
	SET ( BOOST_INCLUDEDIR "C:/Boost/include" )
	SET ( Boost_ADDITIONAL_VERSIONS "1.46" "1.46.1")
	#MESSAGE( "BOOST_ROOT=${BOOST_ROOT}")
ENDIF(WIN32)


SET( Boost_USE_STATIC_LIBS ON )
FIND_PACKAGE ( Boost COMPONENTS system filesystem program_options REQUIRED)

# Find ITK.
SET( ITK_USE_STATIC_LIBS ON )
FIND_PACKAGE(ITK REQUIRED)
INCLUDE(${ITK_USE_FILE})

SET( EXAMPLES_DIR "${CMAKE_CURRENT_SOURCE_DIR}/../Examples" )

INCLUDE_DIRECTORIES( ${Boost_INCLUDE_DIRS} )
LINK_DIRECTORIES   ( ${Boost_LIBRARY_DIRS} )


SET(SOURCE_DIR "${CMAKE_CURRENT_SOURCE_DIR}" )

SET( regseg "regseg")
ADD_EXECUTABLE(${regseg} "${SOURCE_DIR}/regseg.cxx" ) 
TARGET_LINK_LIBRARIES(${regseg} ${ITK_LIBRARIES} ${Boost_PROGRAM_OPTIONS_LIBRARY})


SET( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -D_CRT_SECURE_NO_DEPRECATE -D_SCL_SECURE_NO_WARNINGS")
SET( CPACK_DEBIAN_PACKAGE_DEPENDS "libc6 (>= 2.13-ubuntu13), libgcc1 (>= 1:4.5.2-8ubuntu4), libboost-system1.42.0 (>= 1.42.0-4ubuntu2), libboost-filesystem1.42.0 (>= 1.42.0-4ubuntu2), libboost-program-options1.42.0 (>= 1.42.0-4ubuntu2)")


if(BUILD_EXAMPLES)
  add_subdirectory(${EXAMPLES_DIR} "${CMAKE_CURRENT_BINARY_DIR}/examples")
endif()

#IF(CMAKE_COMPILER_IS_GNUCC)
#  SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fmessage-length=0")
#ENDIF(CMAKE_COMPILER_IS_GNUCC)
#IF(CMAKE_COMPILER_IS_GNUCXX)
#  SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fmessage-length=0")
#ENDIF(CMAKE_COMPILER_IS_GNUCXX)
#SET(CMAKE_VERBOSE_MAKEFILE ON)
