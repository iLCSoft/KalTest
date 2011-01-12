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



### DEPENDENCIES ############################################################

# FindILCSOFT_CMAKE_MODULES.cmake needs to be in CMAKE_MODULE_PATH
FIND_PACKAGE( ILCSOFT_CMAKE_MODULES 0.11 REQUIRED )

# load default settings from ILCSOFT_CMAKE_MODULES
INCLUDE( ilcsoft_default_settings )
#INCLUDE( build_32bit_compatible )


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
ADD_SUBDIRECTORY( ./src )

# if this option is set to ON examples will be built with default target 'make'
# if set to OFF examples will only be built with target 'make examples'
OPTION( BUILD_EXAMPLES "Set to ON to build examples" OFF )

# examples (build with 'make examples')
ADD_SUBDIRECTORY( ./examples/kaltest )


#############################################################################

# display some variables and write them to cache
DISPLAY_STD_VARIABLES()

MESSAGE( STATUS "BUILD_EXAMPLES = ${BUILD_EXAMPLES}" )
MESSAGE( STATUS )

# generate and install following configuration files
GENERATE_PACKAGE_CONFIGURATION_FILES( KalTestConfig.cmake KalTestConfigVersion.cmake KalTestLibDeps.cmake )

