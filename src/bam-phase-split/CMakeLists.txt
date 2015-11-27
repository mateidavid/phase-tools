message(STATUS "Processing: ${CMAKE_CURRENT_SOURCE_DIR}")

add_executable(bam-phase-split bam-phase-split.cpp)
target_link_libraries(bam-phase-split
    common
    thirdparty
    version
    ${HTSLib_LIBRARIES}
    ${ZLIB_LIBRARIES}
    )

install(TARGETS
    bam-phase-split
    RUNTIME DESTINATION bin
    )
