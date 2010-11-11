# build KalTest library

SET( lib_input_dirs geomlib kallib kaltracklib utils )

FOREACH( lib_input_dir ${lib_input_dirs} )
    LIST( APPEND ROOT_DICT_INCLUDE_DIRS ${CMAKE_CURRENT_SOURCE_DIR}/${lib_input_dir} )
ENDFOREACH()

#MESSAGE( STATUS "ROOT_DICT_INCLUDE_DIRS: ${ROOT_DICT_INCLUDE_DIRS}" )


FOREACH( lib_input_dir ${lib_input_dirs} )

    AUX_SOURCE_DIRECTORY( ${lib_input_dir} lib_sources )

    PREPARE_ROOT_DICT_HEADERS( ${lib_input_dir} )

    INSTALL( DIRECTORY ${lib_input_dir}/ DESTINATION "include/kaltest"
        FILES_MATCHING PATTERN "*.h" PATTERN "LinkDef.h" EXCLUDE )

    GEN_ROOT_DICT_SOURCES( ${lib_input_dir}Dict.cxx )
    
    LIST( APPEND lib_sources ${ROOT_DICT_OUTPUT_SOURCES} )

ENDFOREACH()

INCLUDE_DIRECTORIES( ${ROOT_DICT_INCLUDE_DIRS} )
INCLUDE_DIRECTORIES( ${ROOT_INCLUDE_DIRS} )

#MESSAGE( STATUS "KalTest lib sources: ${lib_sources}" )

ADD_LIBRARY( KalTest ${lib_sources} )

ADD_CUSTOM_TARGET( lib DEPENDS KalTest )

TARGET_LINK_LIBRARIES( KalTest ${ROOT_LIBRARIES} )

SET_TARGET_PROPERTIES( KalTest PROPERTIES
    # create *nix style library versions + symbolic links
    VERSION ${KALTEST_VERSION}
    SOVERSION ${KALTEST_SOVERSION}
    # allow creating static and shared libs without conflicts
    CLEAN_DIRECT_OUTPUT 1
)

IF( APPLE ) #---- need special linker flags for ROOT dictionary on MacOS
    SET_TARGET_PROPERTIES( KalTest PROPERTIES
        LINK_FLAGS "-single_module -undefined dynamic_lookup -bind_at_load"
    )
ENDIF( APPLE )

# install library
INSTALL( TARGETS KalTest DESTINATION lib PERMISSIONS
    OWNER_READ OWNER_WRITE OWNER_EXECUTE
    GROUP_READ GROUP_EXECUTE
    WORLD_READ WORLD_EXECUTE
)


