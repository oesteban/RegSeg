cmake_minimum_required(VERSION 2.8 FATAL_ERROR)
#######################################

FIND_PACKAGE(ITK)
IF(ITK_FOUND)
  INCLUDE(${ITK_USE_FILE})
ELSE(ITK_FOUND)
  MESSAGE(FATAL_ERROR
          "Please Set ITK_DIR")
ENDIF(ITK_FOUND)


SET( CORE_SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/Modules/Core )

SET( LEVELSETS_SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/Modules/LevelSets )
SET( LEVELSETS_LIBRARY_DIR ${CMAKE_CURRENT_BINARY_DIR}/Modules/LevelSets )

SET( OPTIMIZERS_SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/Modules/Optimizers )
SET( OPTIMIZERS_LIBRARY_DIR ${CMAKE_CURRENT_BINARY_DIR}/Modules/Optimizers )


INCLUDE_DIRECTORIES( ${OPTIMIZERS_SOURCE_DIR}/include
                     ${LEVELSETS_SOURCE_DIR}/include
                     ${CORE_SOURCE_DIR}/include
                     
)

ADD_SUBDIRECTORY( Modules/Core )

ADD_SUBDIRECTORY( Modules/LevelSets )
SET( LEVELSETS_LIB "RSTKLevelSets" )

ADD_SUBDIRECTORY( Modules/Optimizers )
SET( OPTIMIZERS_LIB "RSTKOptimizers" )




#######################################
# if(WIN32 OR NOT EXISTS /dev/urandom)
#   cmake_minimum_required(VERSION 2.8.5)
# endif()
# 
# project(ACWERegistration)
# 
# if( WIN32 )
#   string( LENGTH "${CMAKE_CURRENT_SOURCE_DIR}" n )
#   if( n GREATER 50 )
#     message(
#       FATAL_ERROR
#       "ACWERegistration source code directory path length is too long (${n} > 50)."
#       "Please move the ACWERegistration source code directory to a directory with a shorter path."
#       )
#   endif()
# 
#   string( LENGTH "${CMAKE_CURRENT_BINARY_DIR}" n )
#   if( n GREATER 50 )
#     message(
#       FATAL_ERROR
#       "ACWERegistration build directory path length is too long (${n} > 50)."
#       "Please set the ACWERegistration build directory to a directory with a shorter path."
#       )
#   endif()
# endif()
# 
# 
# set( CMAKE_MODULE_PATH /home/oesteban/workspace/ITK/src/CMake ${CMAKE_MODULE_PATH} )
# include(ITKModuleMacros)
# 
# SET (RSTK_VERSION_MAJOR 0)
# SET (RSTK_VERSION_MINOR 1)
# SET (RSTK_RELEASE "0")
# 
# if(NOT RSTK_INSTALL_RUNTIME_DIR)
#   set(RSTK_INSTALL_RUNTIME_DIR bin)
# endif()
# if(NOT RSTK_INSTALL_LIBRARY_DIR)
#   set(RSTK_INSTALL_LIBRARY_DIR lib)
# endif()
# if(NOT RSTK_INSTALL_ARCHIVE_DIR)
#   set(RSTK_INSTALL_ARCHIVE_DIR lib)
# endif()
# if(NOT RSTK_INSTALL_INCLUDE_DIR)
#   set(RSTK_INSTALL_INCLUDE_DIR include/RSTK-${RSTK_VERSION_MAJOR}.${RSTK_VERSION_MINOR})
# endif()
# if(NOT RSTK_INSTALL_DATA_DIR)
#   set(RSTK_INSTALL_DATA_DIR share/RSTK-${RSTK_VERSION_MAJOR}.${RSTK_VERSION_MINOR})
# endif()
# if(NOT RSTK_INSTALL_DOC_DIR)
#   set(RSTK_INSTALL_DOC_DIR share/doc/RSTK-${RSTK_VERSION_MAJOR}.${RSTK_VERSION_MINOR})
# endif()
# if(NOT RSTK_INSTALL_PACKAGE_DIR)
#   set(RSTK_INSTALL_PACKAGE_DIR "lib/cmake/RSTK-${RSTK_VERSION_MAJOR}.${RSTK_VERSION_MINOR}")
# endif()
# 
# configure_file (
#   "${PROJECT_SOURCE_DIR}/RSTKConfig.h.in"
#   "${PROJECT_BINARY_DIR}/RSTKConfig.h"
# )
# 
# FIND_PACKAGE(ITK)
# IF(ITK_FOUND)
#   INCLUDE(${ITK_USE_FILE})
# ELSE(ITK_FOUND)
#   MESSAGE(FATAL_ERROR
#           "Please Set ITK_DIR")
# ENDIF(ITK_FOUND)
# 
# #-----------------------------------------------------------------------------
# # Set a default build type if none was specified
# if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
#   message(STATUS "Setting build type to 'Release' as none was specified.")
#   set(CMAKE_BUILD_TYPE Release CACHE STRING "Choose the type of build." FORCE)
#   # Set the possible values of build type for cmake-gui
#   set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug" "Release" "MinSizeRel" "RelWithDebInfo")
# endif()
# 
# #-----------------------------------------------------------------------------
# # RSTK build configuration options.
# option(BUILD_SHARED_LIBS "Build RSTK with shared libraries." OFF)
# set(RSTK_BUILD_SHARED_LIBS ${BUILD_SHARED_LIBS})
# 
# 
# message(INFO " RSTK_BINARY_DIR=${RSTK_BINARY_DIR}" )
# # Setup build locations.
# if(NOT CMAKE_RUNTIME_OUTPUT_DIRECTORY)
#   set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${RSTK_BINARY_DIR}/bin)
# endif()
# if(NOT CMAKE_LIBRARY_OUTPUT_DIRECTORY)
#   set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${RSTK_BINARY_DIR}/lib)
# endif()
# if(NOT CMAKE_ARCHIVE_OUTPUT_DIRECTORY)
#   set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${RSTK_BINARY_DIR}/lib)
# endif()
# set(RSTK_MODULES_DIR "${RSTK_BINARY_DIR}/${RSTK_INSTALL_PACKAGE_DIR}/Modules")
# 
# #----------------------------------------------------------------------
# # Load the module DAG.
# set(RSTK_MODULES_ALL)
# file(GLOB meta RELATIVE "${PROJECT_SOURCE_DIR}"
#    "${PROJECT_SOURCE_DIR}/*/*/itk-module.cmake" # grouped modules
#   )
# 
# 
# foreach(f ${meta})
#   include(${PROJECT_SOURCE_DIR}/${f})
#   list(APPEND RSTK_MODULES_ALL ${itk-module})
#   get_filename_component(${itk-module}_BASE ${f} PATH)
#   set(${itk-module}_SOURCE_DIR ${PROJECT_SOURCE_DIR}/${${itk-module}_BASE})
#   message(INFO "Module detected: ${${itk-module}_SOURCE_DIR}" )
#   set(${itk-module}_BINARY_DIR ${RSTK_BINARY_DIR}/${${itk-module}_BASE})
#   if(BUILD_TESTING AND EXISTS ${${itk-module}_SOURCE_DIR}/test)
#     list(APPEND RSTK_MODULES_ALL ${itk-module-test})
#     set(${itk-module-test}_SOURCE_DIR ${${itk-module}_SOURCE_DIR}/test)
#     set(${itk-module-test}_BINARY_DIR ${${itk-module}_BINARY_DIR}/test)
#     set(${itk-module-test}_IS_TEST 1)
#     set(${itk-module}_TESTED_BY ${itk-module-test})
#     set(${itk-module-test}_TESTS_FOR ${itk-module})
#   endif()
# 
#   # Reject bad dependencies.
#   string(REGEX MATCHALL ";(RSTKDeprecated|RSTKReview|RSTKIntegratedTest);"
#     _bad_deps ";${RSTK_MODULE_${itk-module}_DEPENDS};${RSTK_MODULE_${itk-module-test}_DEPENDS};")
#   foreach(dep ${_bad_deps})
#     if(NOT "${itk-module}" MATCHES "^(${dep}|RSTKIntegratedTest)$")
#       message(FATAL_ERROR
#         "Module \"${itk-module}\" loaded from\n"
#         "  ${${itk-module}_BASE}/itk-module.cmake\n"
#         "may not depend on module \"${dep}\".")
#     endif()
#   endforeach()
# endforeach()
# # Clear variables set later in each module.
# unset(itk-module)
# unset(itk-module-test)
# 
#  
# SET( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -D_CRT_SECURE_NO_DEPRECATE -D_SCL_SECURE_NO_WARNINGS")
