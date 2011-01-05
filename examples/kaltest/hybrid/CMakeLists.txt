

SET( hybrid_input_dirs gen it kern tpc vtx )

FOREACH( _input_dir ${hybrid_input_dirs} )
    LIST( APPEND ROOT_DICT_INCLUDE_DIRS ${CMAKE_CURRENT_SOURCE_DIR}/${_input_dir} )
ENDFOREACH()


FOREACH( _input_dir ${hybrid_input_dirs} )

    AUX_SOURCE_DIRECTORY( ${_input_dir} hybrid_sources )

    PREPARE_ROOT_DICT_HEADERS( ${_input_dir} )

    GEN_ROOT_DICT_SOURCES( ${_input_dir}Dict.cxx )

    LIST( APPEND hybrid_sources ${ROOT_DICT_OUTPUT_SOURCES} )

ENDFOREACH()

INCLUDE_DIRECTORIES( ./main )
INCLUDE_DIRECTORIES( ${ROOT_DICT_INCLUDE_DIRS} )
INCLUDE_DIRECTORIES( BEFORE ${ROOT_INCLUDE_DIRS} )

AUX_SOURCE_DIRECTORY( ./main hybrid_sources )

ADD_KALTEST_EXAMPLE( hybrid ${hybrid_sources} )

