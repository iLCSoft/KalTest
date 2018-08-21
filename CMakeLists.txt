########################################################
# cmake file for building KalFit libraries
# @author Jan Engels, DESY
# @author F. Gaede, DESY
CMAKE_MINIMUM_REQUIRED( VERSION 2.6 FATAL_ERROR )
########################################################


# project name
PROJECT( KalTest )


# project version
SET( KalTest_VERSION_MAJOR 2 )
SET( KalTest_VERSION_MINOR 5 )
SET( KalTest_VERSION_PATCH 0 )

cmake_policy(SET CMP0008 NEW)  

### DEPENDENCIES ############################################################

FIND_PACKAGE( ILCUTIL COMPONENTS ILCSOFT_CMAKE_MODULES REQUIRED )
FIND_PACKAGE( ROOT REQUIRED COMPONENTS Eve)

LINK_LIBRARIES( ${ROOT_COMPONENT_LIBRARIES} )

# load default settings from ILCSOFT_CMAKE_MODULES
INCLUDE( ilcsoft_default_settings )


# macros for generating root dict sources with rootcint
SET( ROOT_DICT_CINT_DEFINITIONS "-DHANDLE_DICT_EXCEPTIONS=IGNORED_FOR_CINT" )

#  ROOTConfig.cmake uses a different name for rootcint than we used to in FINDROOT.cmake 
#SET( ROOT_CINT_EXECUTABLE ${ROOT_rootcint_CMD} )

INCLUDE( MacroRootDict )



### LIBRARY AND TOOLS #######################################################


#---- use 5 dimensional track state:
OPTION( BUILD_WITH_T0_FIT "Set to ON to build with t0 fit (kSdim=6)" OFF )
IF( NOT BUILD_WITH_T0_FIT)
 ADD_DEFINITIONS( -D __NOT0__ )
ENDIF()

ADD_SUBDIRECTORY( ./src )


# --------- kaltest examples -------------------------------------------------
# if this option is set to ON examples will be built with default target 'make'
# if set to OFF examples will only be built with target 'make examples'
OPTION( BUILD_EXAMPLES "Set to ON to build examples" OFF )
MESSAGE( STATUS "BUILD_EXAMPLES = ${BUILD_EXAMPLES}" )
ADD_SUBDIRECTORY( ./examples/kaltest )
# ----------------------------------------------------------------------------



# display some variables and write them to cache
DISPLAY_STD_VARIABLES()

# generate and install following configuration files
GENERATE_PACKAGE_CONFIGURATION_FILES( KalTestConfig.cmake KalTestConfigVersion.cmake KalTestLibDeps.cmake )

