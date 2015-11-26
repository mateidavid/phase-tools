message(STATUS "Processing: ${CMAKE_CURRENT_SOURCE_DIR}")

file(COPY
    add-parent-phasing extend-phase-set random-flip-phase-set collapse-phase-sets
    DESTINATION .
    )

install(PROGRAMS
    add-parent-phasing extend-phase-set random-flip-phase-set collapse-phase-sets
    DESTINATION bin
    )
