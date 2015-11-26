#
# Build Boost as an External Project
#
# Arguments:
# - PREFIX: External project prefix. (required)
# - VERSION: Boost version to download and build.
# - URL: Boost tarball URL. Either VERSION or URL must be given. URL
#     takes precedence over VERSION.
# - LIBRARIES: List of Boost libraries to build. If empty, only header-only
#     libraries are built.
# - B2_ARGUMENTS: Arguments to be passed to b2.
#     By default: "--build-dir=build variant=release link=static threading=multi"
#
# Set variables:
# - Boost_INCLUDE_DIRS: Boost headers directory.
# - Boost_<UPPER_CASE_COMPONENT>_LIBRARY: Boost library.
#

include(ExternalProject)

function(build_boost)
    # parse arguments
    set(one_value_args VERSION PREFIX URL)
    set(multi_value_args LIBRARIES B2_ARGUMENTS)
    cmake_parse_arguments(build_boost "" "${one_value_args}" "${multi_value_args}" ${ARGN})

    # check arguments
    if (NOT build_boost_PREFIX)
        status(FATAL_ERROR "PREFIX is required.")
    endif()
    if (NOT build_boost_VERSION AND NOT build_boost_URL)
        status(FATAL_ERROR "Either VERSION of URL is required.")
    endif()

    if(NOT build_boost_URL)
        string(REPLACE "." "_" BOOST_VERSION_UNDERSCORE ${build_boost_VERSION})
        set(build_boost_URL "http://sourceforge.net/projects/boost/files/boost/${build_boost_VERSION}/boost_${BOOST_VERSION_UNDERSCORE}.tar.bz2/download")
        message(STATUS "Building Boost (${build_boost_VERSION}) headers and libraries [${build_boost_LIBRARIES}] in: ${build_boost_PREFIX}")
    else()
        message(STATUS "Building Boost (${build_boost_URL}) headers and libraries [${build_boost_LIBRARIES}] in: ${build_boost_PREFIX}")
    endif()

    if(NOT build_boost_LIBRARIES)
        if(UNIX)
            set(HEADERS_INSTALL_COMMAND mkdir -p ${build_boost_PREFIX}/include && cp -r boost ${build_boost_PREFIX}/include)
        else()
            set(HEADERS_INSTALL_COMMAND "")
        endif()
        ExternalProject_Add(Boost
            PREFIX ${build_boost_PREFIX}
            URL ${build_boost_URL}
            CONFIGURE_COMMAND ""
            BUILD_IN_SOURCE 1
            BUILD_COMMAND ""
            INSTALL_COMMAND ${HEADERS_INSTALL_COMMAND}
            )
    else()
        if(UNIX)
            set(BOOST_BOOTSTRAP "./bootstrap.sh")
            set(BOOST_B2 "./b2")
        else()
            set(BOOST_BOOTSTRAP "bootstrap.bat")
            set(BOOST_B2 "b2")
        endif()
        if(NOT build_boost_B2_ARGUMENTS)
            set(build_boost_B2_ARGUMENTS --build-dir=build variant=release link=static threading=multi)
        endif()
        message(STATUS "Using Boost b2 arguments: ${build_boost_B2_ARGUMENTS}")
        string(REPLACE ";" "," BOOST_LIBRARIES_COMMA "${build_boost_LIBRARIES}")
        set(BOOST_BOOTSTRAP_COMMAND ${BOOST_BOOTSTRAP} --prefix=${build_boost_PREFIX} --with-libraries=${BOOST_LIBRARIES_COMMA})
        set(BOOST_B2_COMMAND ${BOOST_B2} ${build_boost_B2_ARGUMENTS} install)
        ExternalProject_Add(Boost
            PREFIX ${build_boost_PREFIX}
            URL ${build_boost_URL}
            CONFIGURE_COMMAND ${BOOST_BOOTSTRAP_COMMAND}
            BUILD_IN_SOURCE 1
            BUILD_COMMAND ${BOOST_B2_COMMAND}
            INSTALL_COMMAND ""
            )
        ExternalProject_Get_Property(Boost INSTALL_DIR)
        foreach(l ${build_boost_LIBRARIES})
            string(TOUPPER l l_upper)
            string(TOLOWER l l_lower)
            set(Boost_${l_upper}_LIBRARY
                ${build_boost_PREFIX}/lib/${CMAKE_STATIC_LIBRARY_PREFIX}boost_${l_lower}{CMAKE_STATIC_LIBRARY_SUFFIX}
                PARENT_SCOPE
                )
        endforeach()
    endif()
    ExternalProject_Get_Property(Boost INSTALL_DIR)
    set(Boost_INCLUDE_DIRS ${INSTALL_DIR}/include PARENT_SCOPE)
endfunction()
