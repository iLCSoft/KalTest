####################################
# cmake file for building KalFit libraries
# @author Jan Engels, DESY
# @author F. Gaede, DESY
####################################


########################################################
# CMake compatibility issues: don't modify this, please!
CMAKE_MINIMUM_REQUIRED( VERSION 2.4.6 )
#SET( CMAKE_BACKWARDS_COMPATIBILITY 2.4.6 )
MARK_AS_ADVANCED(CMAKE_BACKWARDS_COMPATIBILITY)
# allow more human readable "if then else" constructs
SET( CMAKE_ALLOW_LOOSE_LOOP_CONSTRUCTS TRUE )
########################################################

# project name
PROJECT( KalTest )

# project options
OPTION( BUILD_SHARED_LIBS "Set to OFF to build static libraries" ON )

# project version
SET( KALTEST_MAJOR_VERSION 0 )
SET( KALTEST_MINOR_VERSION 6 )
SET( KALTEST_PATCH_LEVEL 0 )

# set default install prefix to project root directory
IF( CMAKE_INSTALL_PREFIX STREQUAL "/usr/local" )
    SET( CMAKE_INSTALL_PREFIX "${PROJECT_SOURCE_DIR}" )
ENDIF()

##---------------- 32/64 bit issues ---------------------------------------
#OPTION( BUILD_32BIT_COMPATIBLE "Set to ON to build in 32 bit compatibility mode" OFF )
##INCLUDE( CheckTypeSize )
##CHECK_TYPE_SIZE( "void*" SIZE_OF_VOID_P )
#
#IF( CMAKE_SIZEOF_VOID_P EQUAL 4 )
#    MESSAGE( STATUS "32 bit architecture detected" )
#ENDIF()
#
#IF( CMAKE_SIZEOF_VOID_P EQUAL 8 )
#    MESSAGE( STATUS "64 bit architecture detected" )
#
#    IF( BUILD_32BIT_COMPATIBLE )
#        IF( COMMAND SET_PROPERTY )
#            SET_PROPERTY(GLOBAL PROPERTY FIND_LIBRARY_USE_LIB64_PATHS 0)
#            SET( CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -m32" )
#        ELSE()
#            MESSAGE( FATAL_ERROR "to build in 32 bit compatibility mode you need cmake >= 2.6" )
#        ENDIF()
#        MESSAGE( STATUS "Build in 32 bit compatibility mode" )
#    ELSE()
#        MESSAGE( STATUS "Build in native 64 bit mode" )
#    ENDIF()
#ENDIF()
##-------------------------------------------------------------------------

# library *nix style versioning
SET( KALTEST_SOVERSION "${KALTEST_MAJOR_VERSION}.${KALTEST_MINOR_VERSION}" )
SET( KALTEST_VERSION "${KALTEST_SOVERSION}.${KALTEST_PATCH_LEVEL}" )

# cmake output directories
SET( EXECUTABLE_OUTPUT_PATH "${PROJECT_BINARY_DIR}/bin" CACHE PATH
    "EXECUTABLE_OUTPUT_PATH" FORCE )
SET( LIBRARY_OUTPUT_PATH "${PROJECT_BINARY_DIR}/lib" CACHE PATH
    "LIBRARY_OUTPUT_PATH" FORCE )
MARK_AS_ADVANCED( EXECUTABLE_OUTPUT_PATH LIBRARY_OUTPUT_PATH )

###################################################################################

# require proper c++
ADD_DEFINITIONS( "-Wall -ansi -pedantic" )
#----- need long long for int64 for now ------
ADD_DEFINITIONS( "-Wno-long-long -fno-strict-aliasing" )

#-----------------------------------------------------

# macros for generating root dict sources with rootcint
SET( ROOT_DICT_CINT_DEFINITIONS "-DHANDLE_DICT_EXCEPTIONS=IGNORED_FOR_CINT" )
INCLUDE( cmake/MacroRootDict.cmake )

###################################################################################

# KalTest library
ADD_SUBDIRECTORY( src )

# build with 'make examples'
ADD_SUBDIRECTORY( examples/kaltest )
###################################################################################


# create uninstall configuration file 
CONFIGURE_FILE( "${PROJECT_SOURCE_DIR}/cmake/cmake_uninstall.cmake.in"
                "${PROJECT_BINARY_DIR}/cmake_uninstall.cmake"
                IMMEDIATE @ONLY )

# create uninstall target
ADD_CUSTOM_TARGET( uninstall
    "${CMAKE_COMMAND}" -P "${PROJECT_BINARY_DIR}/cmake_uninstall.cmake"
)

# create configuration file from .in file
CONFIGURE_FILE( "${PROJECT_SOURCE_DIR}/cmake/KalTestConfig.cmake.in"
                "${PROJECT_BINARY_DIR}/KalTestConfig.cmake" @ONLY
)

# install configuration file
INSTALL( FILES "${PROJECT_BINARY_DIR}/KalTestConfig.cmake" DESTINATION . )

# display status message for important variables
MESSAGE( STATUS )
MESSAGE( STATUS "-------------------------------------------------------------------------------" )
MESSAGE( STATUS "BUILD_SHARED_LIBS = ${BUILD_SHARED_LIBS}" )
MESSAGE( STATUS "CMAKE_INSTALL_PREFIX = ${CMAKE_INSTALL_PREFIX}" )
#IF( CMAKE_SIZEOF_VOID_P EQUAL 8 )
#    MESSAGE( STATUS "BUILD_32BIT_COMPATIBLE = ${BUILD_32BIT_COMPATIBLE}" )
#ENDIF()
MESSAGE( STATUS "Change a value with: cmake -D<Variable>=<Value>" )
MESSAGE( STATUS "-------------------------------------------------------------------------------" )
MESSAGE( STATUS )

# force some variables that could be defined in the command line
# to be written to cache
#SET( BUILD_SHARED_LIBS "${BUILD_SHARED_LIBS}" CACHE BOOL
#    "Set to OFF to build static libraries" FORCE )
#SET( CMAKE_INSTALL_PREFIX "${CMAKE_INSTALL_PREFIX}" CACHE PATH
#    "Where to install KALTEST" FORCE )

# export build settings
INCLUDE( CMakeExportBuildSettings )
CMAKE_EXPORT_BUILD_SETTINGS( "KALTESTBuildSettings.cmake" )
INSTALL( FILES "${PROJECT_BINARY_DIR}/KALTESTBuildSettings.cmake" DESTINATION lib/cmake )

# export library dependencies (keep this as the last line in the file)
EXPORT_LIBRARY_DEPENDENCIES( "KALTESTLibDeps.cmake" )
INSTALL( FILES "${PROJECT_BINARY_DIR}/KALTESTLibDeps.cmake" DESTINATION lib/cmake )

