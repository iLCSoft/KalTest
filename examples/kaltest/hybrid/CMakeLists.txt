

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

INCLUDE_DIRECTORIES( main )
INCLUDE_DIRECTORIES( ${ROOT_DICT_INCLUDE_DIRS} )
INCLUDE_DIRECTORIES( BEFORE ${ROOT_INCLUDE_DIRS} )

AUX_SOURCE_DIRECTORY( main hybrid_sources )

# for backwards compatibility with cmake 2.4
IF( "${CMAKE_MAJOR_VERSION}.${CMAKE_MINOR_VERSION}" LESS 2.6 )
    ADD_CUSTOM_TARGET( examples )
ENDIF()

ADD_KALTEST_EXAMPLE( hybrid ${hybrid_sources} )

