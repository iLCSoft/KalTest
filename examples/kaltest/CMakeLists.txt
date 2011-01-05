

SET( ROOT_DICT_INCLUDE_DIRS 
    ${PROJECT_SOURCE_DIR}/src/utils
    ${PROJECT_SOURCE_DIR}/src/kallib
    ${PROJECT_SOURCE_DIR}/src/geomlib
    ${PROJECT_SOURCE_DIR}/src/kaltracklib
)

INCLUDE_DIRECTORIES( ${ROOT_INCLUDE_DIRS} )
INCLUDE_DIRECTORIES( BEFORE ${ROOT_DICT_INCLUDE_DIRS} )



# ====== macro for adding new examples within KalTest =======================
ADD_CUSTOM_TARGET( examples )
MACRO( ADD_KALTEST_EXAMPLE _example )

    IF( BUILD_EXAMPLES )
        ADD_EXECUTABLE( ${_example} ${ARGN} )
    ELSE()
        ADD_EXECUTABLE( ${_example} EXCLUDE_FROM_ALL ${ARGN} )
    ENDIF()

    ADD_DEPENDENCIES( examples ${_example} )

    TARGET_LINK_LIBRARIES( ${_example} KalTest )
ENDMACRO( ADD_KALTEST_EXAMPLE )
# ===========================================================================



# ====== macro for adding new examples within KalTest =======================
MACRO( ADD_KALTEST_ROOT_DICT_EXAMPLE _example _input_dir )

    PREPARE_ROOT_DICT_HEADERS( ${_input_dir} )

    GEN_ROOT_DICT_SOURCES( ${_example}Dict.cxx )

    #SET( sources )
    #AUX_SOURCE_DIRECTORY( ${_input_dir} sources )
    FILE( GLOB sources "${_input_dir}/*.cxx" )

    ADD_KALTEST_EXAMPLE( ${_example} ${ROOT_DICT_OUTPUT_SOURCES} ${sources} )

ENDMACRO( ADD_KALTEST_ROOT_DICT_EXAMPLE )
# ===========================================================================


ADD_KALTEST_ROOT_DICT_EXAMPLE( cdc cdc )
ADD_KALTEST_ROOT_DICT_EXAMPLE( ct ct )
ADD_KALTEST_ROOT_DICT_EXAMPLE( simple simple )


# hybrid
ADD_SUBDIRECTORY( ./hybrid )

