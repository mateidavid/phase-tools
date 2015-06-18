find_package(Git)

if(GIT_FOUND)
    execute_process(COMMAND ${GIT_EXECUTABLE} describe --tags --always
        WORKING_DIRECTORY ${SRC_ROOT}
        OUTPUT_VARIABLE VERSION OUTPUT_STRIP_TRAILING_WHITESPACE
    )
    set(VERSION_SOURCE " (git)")
elseif(EXISTS ${SRC_ROOT}/GIT_VERSION)
    execute_process(COMMAND cat ${SRC_ROOT}/GIT_VERSION
        OUTPUT_VARIABLE VERSION OUTPUT_STRIP_TRAILING_WHITESPACE
    )
    set(VERSION_SOURCE " (GIT_VERSION)")
else()
    set(VERSION "unknown")
endif()

configure_file(${SRC_ROOT}/version.h.in ${BIN_ROOT}/version.h @ONLY)
message(STATUS "version: ${VERSION}${VERSION_SOURCE}")
