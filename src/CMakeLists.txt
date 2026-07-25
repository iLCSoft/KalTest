# build KalTest library

SET( lib_input_dirs geomlib kallib kaltracklib utils bfield )

SET( KalTest_PUBLIC_INCLUDE_DIR "${PROJECT_SOURCE_DIR}/include" )
SET( KalTest_PUBLIC_HEADER_DIR "${KalTest_PUBLIC_INCLUDE_DIR}/kaltest" )

# One dictionary and one PCM serve the whole KalTest library.  Every header is
# passed by its stable public spelling, kaltest/Header.h.
FILE( GLOB KALTEST_PUBLIC_HEADERS "${KalTest_PUBLIC_HEADER_DIR}/*.h" )
LIST( SORT KALTEST_PUBLIC_HEADERS )

SET( KALTEST_DICT_HEADERS )
SET( KALTEST_DICT_DEPENDS )
FOREACH( _header ${KALTEST_PUBLIC_HEADERS} )
    GET_FILENAME_COMPONENT( _header_name "${_header}" NAME )
    LIST( APPEND KALTEST_DICT_HEADERS "kaltest/${_header_name}" )
    LIST( APPEND KALTEST_DICT_DEPENDS "${_header}" )
ENDFOREACH()

SET( KALTEST_DICT_LINKDEF "${CMAKE_CURRENT_SOURCE_DIR}/KalTestLinkDef.h" )
SET( KALTEST_DICT_SOURCE "${CMAKE_BINARY_DIR}/rootdict/KalTestDict.cxx" )
SET( KALTEST_DICT_PCM "${CMAKE_BINARY_DIR}/rootdict/KalTestDict_rdict.pcm" )
ADD_CUSTOM_COMMAND(
    OUTPUT "${KALTEST_DICT_SOURCE}" "${KALTEST_DICT_PCM}"
    COMMAND "${CMAKE_COMMAND}" -E make_directory "${CMAKE_BINARY_DIR}/rootdict"
    COMMAND ${ROOT_CINT_WRAPPER} -f "${KALTEST_DICT_SOURCE}" -c
        ${ROOT_DICT_CINT_DEFINITIONS} -I${KalTest_PUBLIC_INCLUDE_DIR}
        ${KALTEST_DICT_HEADERS} "${KALTEST_DICT_LINKDEF}"
    WORKING_DIRECTORY "${KalTest_PUBLIC_INCLUDE_DIR}"
    DEPENDS ${KALTEST_DICT_DEPENDS} "${KALTEST_DICT_LINKDEF}"
    COMMENT "generating: ${KALTEST_DICT_SOURCE} ${KALTEST_DICT_PCM}"
)

SET( lib_sources )
FOREACH( lib_input_dir ${lib_input_dirs} )
    AUX_SOURCE_DIRECTORY( ${lib_input_dir} _component_sources )
    LIST( APPEND lib_sources ${_component_sources} )
ENDFOREACH()
LIST( APPEND lib_sources "${KALTEST_DICT_SOURCE}" )

INCLUDE_DIRECTORIES( SYSTEM ${ROOT_INCLUDE_DIRS} )

#MESSAGE( STATUS "KalTest lib sources: ${lib_sources}" )

ADD_SHARED_LIBRARY( KalTest ${lib_sources} )
INSTALL_SHARED_LIBRARY( KalTest DESTINATION lib )
TARGET_LINK_LIBRARIES( KalTest ${ROOT_LIBRARIES} )



IF( APPLE ) #---- need special linker flags for ROOT dictionary on MacOS
    SET_TARGET_PROPERTIES( KalTest PROPERTIES
        LINK_FLAGS "-single_module -undefined dynamic_lookup -bind_at_load"
    )
ENDIF( APPLE )


install( FILES "${KALTEST_DICT_PCM}" DESTINATION lib )
INSTALL( DIRECTORY "${KalTest_PUBLIC_HEADER_DIR}/" DESTINATION "include/kaltest"
    FILES_MATCHING PATTERN "*.h"
)
