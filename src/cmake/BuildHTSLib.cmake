#
# Build HTSLib as an External Project
#
# Arguments:
# - PREFIX: External project prefix. (required)
# - URL: Tarball URL (optinal). Use github by default.
#
# Set variables:
# - HTSLib_INCLUDE_DIRS: HTSLib headers directory.
# - HTSLib_LIBRARIES: HTSLib library.
#

include(ExternalProject)

function(build_htslib)
    # parse arguments
    set(one_value_args PREFIX URL)
    cmake_parse_arguments(build_htslib "" "${one_value_args}" "" ${ARGN})

    # check arguments
    if (NOT build_htslib_PREFIX)
        status(FATAL_ERROR "PREFIX is required.")
    endif()
    if (build_htslib_URL)
        set(htslib_source URL ${build_htslib_URL})
        # tarballs have a configure script that sets prefix
        ExternalProject_Add(HTSLib
            PREFIX ${build_htslib_PREFIX}
            ${htslib_source}
            CONFIGURE_COMMAND ./configure --prefix=${build_htslib_PREFIX}
            BUILD_IN_SOURCE 1
            BUILD_COMMAND make lib-static
            INSTALL_COMMAND make install
            )
    else()
        set(htslib_source GIT_REPOSITORY https://github.com/samtools/htslib.git GIT_TAG 1.2.1)
        # git repo has no configure script, and prefix set in Makefile
        ExternalProject_Add(HTSLib
            PREFIX ${build_htslib_PREFIX}
            ${htslib_source}
            CONFIGURE_COMMAND ""
            BUILD_IN_SOURCE 1
            BUILD_COMMAND make lib-static
            INSTALL_COMMAND make install prefix=${build_htslib_PREFIX}
            )
    endif()
    message(STATUS "Building HTSLib (${htslib_source}) in: ${build_htslib_PREFIX}")
    ExternalProject_Get_Property(HTSLib INSTALL_DIR)
    set(HTSLib_INCLUDE_DIR ${INSTALL_DIR}/include
        CACHE INTERNAL "HTSLib headers")
    set(HTSLib_LIBRARY ${INSTALL_DIR}/lib/${CMAKE_STATIC_LIBRARY_PREFIX}hts${CMAKE_STATIC_LIBRARY_SUFFIX}
        CACHE INTERNAL "HTSLib library")
    mark_as_advanced(HTSLib_INCLUDE_DIR HTSLib_LIBRARY)
    set(HTSLib_INCLUDE_DIRS ${HTSLib_INCLUDE_DIR} PARENT_SCOPE)
    set(HTSLib_LIBRARIES ${HTSLib_LIBRARY} PARENT_SCOPE)
endfunction()
