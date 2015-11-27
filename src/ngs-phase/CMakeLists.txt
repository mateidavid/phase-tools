message(STATUS "Processing: ${CMAKE_CURRENT_SOURCE_DIR}")

add_executable(ngs-phase ngs-phase.cpp)
target_link_libraries(ngs-phase
    common
    thirdparty
    version
    ${HTSLib_LIBRARIES}
    ${ZLIB_LIBRARIES}
    )

install(TARGETS
    ngs-phase
    RUNTIME DESTINATION bin
    )
