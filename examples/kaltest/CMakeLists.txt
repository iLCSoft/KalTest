

SET( ROOT_DICT_INCLUDE_DIRS 
    ${PROJECT_SOURCE_DIR}/src/utils
    ${PROJECT_SOURCE_DIR}/src/kallib
    ${PROJECT_SOURCE_DIR}/src/geomlib
    ${PROJECT_SOURCE_DIR}/src/kaltracklib
)


# ====== macro for adding new examples within KalTest =======================
ADD_CUSTOM_TARGET( examples )
MACRO( ADD_KALTEST_EXAMPLE TESTNAME INPUT_DIR )

PREPARE_ROOT_DICT_HEADERS( ${INPUT_DIR} )

GEN_ROOT_DICT_SOURCES( ${TESTNAME}Dict.cxx )

FILE( GLOB sources "${INPUT_DIR}/*.cxx" )

INCLUDE_DIRECTORIES( ${ROOT_INCLUDE_DIRS} )
INCLUDE_DIRECTORIES( ${ROOT_DICT_INCLUDE_DIRS} )

ADD_EXECUTABLE( ${TESTNAME} EXCLUDE_FROM_ALL
    ${ROOT_DICT_OUTPUT_SOURCES}
    ${sources}
)

ADD_DEPENDENCIES( examples ${TESTNAME} )

TARGET_LINK_LIBRARIES( ${TESTNAME} KalTest )

ENDMACRO( ADD_KALTEST_EXAMPLE )
# ===========================================================================


ADD_KALTEST_EXAMPLE( cdc cdc )
ADD_KALTEST_EXAMPLE( ct ct )
ADD_KALTEST_EXAMPLE( simple simple )


# hybrid
ADD_SUBDIRECTORY( hybrid )

