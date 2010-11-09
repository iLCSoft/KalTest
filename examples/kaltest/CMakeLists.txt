

SET( ROOT_DICT_INCLUDE_DIRS 
    ${PROJECT_SOURCE_DIR}/src/utils
    ${PROJECT_SOURCE_DIR}/src/kallib
    ${PROJECT_SOURCE_DIR}/src/geomlib
    ${PROJECT_SOURCE_DIR}/src/kaltracklib
)

INCLUDE_DIRECTORIES( ${ROOT_INCLUDE_DIRS} )
INCLUDE_DIRECTORIES( BEFORE ${ROOT_DICT_INCLUDE_DIRS} )


ADD_CUSTOM_TARGET( examples )


# ====== macro for adding new examples within KalTest =======================
MACRO( ADD_KALTEST_EXAMPLE EXAMPLE )

    IF( BUILD_EXAMPLES )
        ADD_EXECUTABLE( ${EXAMPLE} ${ARGN} )
    ELSE()
        ADD_EXECUTABLE( ${EXAMPLE} EXCLUDE_FROM_ALL ${ARGN} )
    ENDIF()

    ADD_DEPENDENCIES( examples ${EXAMPLE} )

    TARGET_LINK_LIBRARIES( ${EXAMPLE} KalTest )
ENDMACRO( ADD_KALTEST_EXAMPLE )
# ===========================================================================



# ====== macro for adding new examples within KalTest =======================
MACRO( ADD_KALTEST_ROOT_DICT_EXAMPLE EXAMPLE INPUT_DIR )

    PREPARE_ROOT_DICT_HEADERS( ${INPUT_DIR} )

    GEN_ROOT_DICT_SOURCES( ${EXAMPLE}Dict.cxx )

    #SET( sources )
    #AUX_SOURCE_DIRECTORY( ${INPUT_DIR} sources )
    FILE( GLOB sources "${INPUT_DIR}/*.cxx" )

    ADD_KALTEST_EXAMPLE( ${EXAMPLE} ${ROOT_DICT_OUTPUT_SOURCES} ${sources} )

ENDMACRO( ADD_KALTEST_ROOT_DICT_EXAMPLE )
# ===========================================================================


ADD_KALTEST_ROOT_DICT_EXAMPLE( cdc cdc )
ADD_KALTEST_ROOT_DICT_EXAMPLE( ct ct )
ADD_KALTEST_ROOT_DICT_EXAMPLE( simple simple )


# hybrid
ADD_SUBDIRECTORY( hybrid )

