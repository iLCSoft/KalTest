

SET( hybrid_input_dirs gen it kern tpc vtx )

FOREACH( input_dir ${hybrid_input_dirs} )
    LIST( APPEND ROOT_DICT_INCLUDE_DIRS ${CMAKE_CURRENT_SOURCE_DIR}/${input_dir} )
ENDFOREACH()


FOREACH( input_dir ${hybrid_input_dirs} )

    AUX_SOURCE_DIRECTORY( ${input_dir} hybrid_sources )

    PREPARE_ROOT_DICT_HEADERS( ${input_dir} )
    GEN_ROOT_DICT_SOURCES( ${input_dir}Dict.cxx )

    LIST( APPEND hybrid_sources ${ROOT_DICT_OUTPUT_SOURCES} )

ENDFOREACH()


AUX_SOURCE_DIRECTORY( main hybrid_sources )

INCLUDE_DIRECTORIES( main )
INCLUDE_DIRECTORIES( ${ROOT_INCLUDE_DIRS} )
INCLUDE_DIRECTORIES( ${ROOT_DICT_INCLUDE_DIRS} )

ADD_EXECUTABLE( hybrid EXCLUDE_FROM_ALL ${hybrid_sources} )

ADD_DEPENDENCIES( examples hybrid )

TARGET_LINK_LIBRARIES( hybrid KalTest )

