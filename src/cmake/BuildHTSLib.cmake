#
# Build HTSLib as an External Project
#
# Arguments:
# - PREFIX: External project prefix. (required)
# - URL: Tarball URL (optinal). Use github by default.
#
# Set variables:
# - HTSLib_INCLUDE_DIR: HTSLib headers directory.
# - HTSLib_LIBRARY: HTSLib library.
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
    else()
        #set(htslib_source GIT_REPOSITORY https://github.com/samtools/htslib.git GIT_TAG 1.2.1)
        set(htslib_source URL http://github.com/samtools/htslib/releases/download/1.2.1/htslib-1.2.1.tar.bz2)
    endif()
    message(STATUS "Building HTSLib (${htslib_source}) in: ${build_htslib_PREFIX}")

    ExternalProject_Add(HTSLib
        PREFIX ${build_htslib_PREFIX}
        ${htslib_source}
        CONFIGURE_COMMAND ./configure --prefix=${build_htslib_PREFIX}
        BUILD_IN_SOURCE 1
        BUILD_COMMAND make lib-static
        INSTALL_COMMAND make install
        )
    ExternalProject_Get_Property(HTSLib INSTALL_DIR)
    set(HTSLib_INCLUDE_DIR ${INSTALL_DIR}/include PARENT_SCOPE)
    set(HTSLib_LIBRARY ${INSTALL_DIR}/lib/${CMAKE_STATIC_LIBRARY_PREFIX}hts${CMAKE_STATIC_LIBRARY_SUFFIX} PARENT_SCOPE)
endfunction()
