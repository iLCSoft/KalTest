########################################################
# cmake file for building KalFit libraries
# @author Jan Engels, DESY
# @author F. Gaede, DESY
CMAKE_MINIMUM_REQUIRED( VERSION 2.6 FATAL_ERROR )
########################################################


# project name
PROJECT( KalTest )


# project version
SET( KALTEST_VERSION_MAJOR 1 )
SET( KALTEST_VERSION_MINOR 0 )
SET( KALTEST_VERSION_PATCH 0 )



### SETTINGS ################################################################

# ilcsoft_default_settings.cmake needs to be in CMAKE_MODULE_PATH
INCLUDE( ilcsoft_default_settings )
#INCLUDE( build_32bit_compatible )



### DEPENDENCIES ############################################################

# fix for backwards compatibility
IF( NOT ROOT_DIR AND DEFINED ROOT_HOME )
    SET( ROOT_DIR ${ROOT_HOME} )
ENDIF()

FIND_PACKAGE( ROOT REQUIRED )

# macros for generating root dict sources with rootcint
SET( ROOT_DICT_CINT_DEFINITIONS "-DHANDLE_DICT_EXCEPTIONS=IGNORED_FOR_CINT" )
INCLUDE( ${ROOT_DICT_MACROS_FILE} )



### LIBRARY AND TOOLS #######################################################


# require proper c++
ADD_DEFINITIONS( "-Wall -ansi -pedantic" )
#----- need long long for int64 for now ------
ADD_DEFINITIONS( "-Wno-long-long -fno-strict-aliasing" )

#-----------------------------------------------------


# KalTest library
ADD_SUBDIRECTORY( src )

# examples (build with 'make examples')
ADD_SUBDIRECTORY( examples/kaltest )


#############################################################################

# display some variables and write them to cache
DISPLAY_STD_VARIABLES()

# export library dependencies (keep this as the last line in the file)
EXPORT_LIBRARY_DEPENDENCIES( "KalTestLibDeps.cmake" )
INSTALL( FILES "${PROJECT_BINARY_DIR}/KalTestLibDeps.cmake" DESTINATION lib/cmake )

