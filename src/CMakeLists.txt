message(STATUS "Processing: ${CMAKE_CURRENT_SOURCE_DIR}")

# set build type
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE "Release" CACHE STRING
       "Choose the type of build, options are: Debug Test Release GProf GProfRel."
       FORCE)
endif()
message(STATUS "Build type: ${CMAKE_BUILD_TYPE}")
if(NOT ${CMAKE_BUILD_TYPE} AND NOT ${CMAKE_BUILD_TYPE} STREQUAL "Release")
    set(CMAKE_VERBOSE_MAKEFILE TRUE CACHE BOOL "Enable verbose make output.")
    message(STATUS "Enabling verbose make output.")
endif()

cmake_minimum_required(VERSION 2.8.12)
project(PHASE_TOOLS C CXX)
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake")
include(CheckLibraryExists)
include(CheckIncludeFileCXX)
include(ExternalProject)

# check main file exists; if not, assume the source dir is wrong
if(NOT EXISTS "${PROJECT_SOURCE_DIR}/bam-phase-split/bam-phase-split.cpp")
    message(FATAL_ERROR "${PROJECT_SOURCE_DIR}: source files not found")
endif()

# set project-related variables that appear in config.h.in
#set(PACKAGE_VERSION "0.10.10")
set(PACKAGE_BUGREPORT "mdavid@oicr.on.ca")
set(PACKAGE_URL "")

set(PACKAGE "${PROJECT_NAME}")
set(PACKAGE_NAME "${PROJECT_NAME}")
set(PACKAGE_TARNAME "${PROJECT_NAME}")
set(PACKAGE_STRING "${PACKAGE_NAME} ${PACKAGE_VERSION}")

set(OPTIONAL_SUBMODULE_PREFIX ${PROJECT_SOURCE_DIR}/optional
    CACHE INTERNAL "Directory for installing optional submodules.")
set(HEADER_SUBDIRS
    common thirdparty version
    CACHE INTERNAL "Subdirectories containing header files")
set(SUBDIRS
    bam-phase-split common ngs-phase scripts thirdparty version
    CACHE INTERNAL "Subdirectories to descend into")

if(NOT ${CMAKE_INSTALL_PREFIX})
    set(CMAKE_INSTALL_PREFIX ${PROJECT_BINARY_DIR})
endif()
message(STATUS "Using install prefix: ${CMAKE_INSTALL_PREFIX}")

### Resolve external dependencies
#
# prefer static libraries
#set(CMAKE_FIND_LIBRARY_SUFFIXES .a .so)

# bake-in RPATH to prevent library search problems
#set(CMAKE_INSTALL_RPATH_USE_LINK_PATH true)
#set(CMAKE_BUILD_WITH_INSTALL_RPATH true)

# find zlib
# => ZLIB_INCLUDE_DIRS, ZLIB_LIBRARIES
find_package(ZLIB REQUIRED)

# check for HTSLib
# => HTSLib_INCLUDE_DIR, HTSLib_LIBRARY
if(BUILD_HTSLIB)
    include(BuildHTSLib)
    build_htslib(PREFIX ${OPTIONAL_SUBMODULE_PREFIX} URL ${HTSLIB_URL})
else()
    find_package(HTSLib REQUIRED)
endif()

# check for Boost headers
# => Boost_INCLUDE_DIRS, Boost_<COMPONENT>_LIBRARY
if(BUILD_BOOST)
    include(BuildBoost)
    set(BOOST_VERSION 1.59.0 CACHE STRING "Boost version to build.")
    build_boost(PREFIX ${OPTIONAL_SUBMODULE_PREFIX}
        VERSION ${BOOST_VERSION}
        URL ${BOOST_URL}
        )
else()
    # Use static libraries to avoid LD_RUN_PATH issues
    #set(Boost_USE_STATIC_LIBS true)
    #find_package(Boost 1.57.0 EXACT REQUIRED COMPONENTS iostreams program_options)
    find_package(Boost 1.57.0 REQUIRED)
    if(NOT Boost_FOUND)
        message(FATAL_ERROR "Boost not found")
    endif()
endif()
message(STATUS "Using Boost headers in: ${Boost_INCLUDE_DIRS}")

### Prepare default compilation flags
#
# headers
include_directories(SYSTEM
    ${HTSLib_INCLUDE_DIRS}
    ${Boost_INCLUDE_DIRS}
    ${ZLIB_INCLUDE_DIRS}
    tclap/include
    hpptools/include)
include_directories(
    ${HEADER_SUBDIRS}
    ${PROJECT_BINARY_DIR})
get_directory_property(include_directories INCLUDE_DIRECTORIES)
message(STATUS "INCLUDE_DIRECTORIES='${include_directories}'")

### general compile flags
#
set(EXTRA_FLAGS "-std=c++11 -pthread -Wall -Wextra -pedantic")

# compiler-specific flags
if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
    set(EXTRA_FLAGS "${EXTRA_FLAGS} -fmax-errors=1")
elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
    set(EXTRA_FLAGS "${EXTRA_FLAGS} -ferror-limit=1")
endif()

# consolidate flags
set(CMAKE_CXX_FLAGS "$ENV{CXXFLAGS} ${EXTRA_FLAGS}")
message(STATUS "CMAKE_CXX_FLAGS='${CMAKE_CXX_FLAGS}'")

### build-specific flags
#
set(CMAKE_CXX_FLAGS_DEBUG "-O0 -g3 -ggdb -fno-inline -fno-eliminate-unused-debug-types")
set(CMAKE_CXX_FLAGS_TEST "-O3 -g3")
set(CMAKE_CXX_FLAGS_TEST_O2 "-O2 -g3")
set(CMAKE_CXX_FLAGS_TEST_O1 "-O1 -g3")
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -DNDEBUG -DDISABLE_ASSERTS")
set(CMAKE_CXX_FLAGS_GPROF "-O3 -g3 -pg")
set(CMAKE_CXX_FLAGS_GPROFREL "-O3 -DNDEBUG -DDISABLE_ASSERTS -pg")

# link flags
#set(CMAKE_EXE_LINKER_FLAGS "-Wl,-rpath=$ENV{LD_RUN_PATH} -Wl,--as-needed")
set(CMAKE_EXE_LINKER_FLAGS "-Wl,--as-needed")
message(STATUS "CMAKE_EXE_LINKER_FLAGS='${CMAKE_EXE_LINKER_FLAGS}'")

### Add target that updates version
#
add_custom_target(package_version
    ${PROJECT_SOURCE_DIR}/get-dir-version -v -d ${PROJECT_SOURCE_DIR} -r ${PROJECT_SOURCE_DIR}/.. -o ${PROJECT_BINARY_DIR}/package_version.h -t PACKAGE_VERSION
    )

### Enable testing targets
#
#enable_testing()

### Descend into subdirectories
#
foreach(dir ${SUBDIRS})
    add_subdirectory(${dir})
endforeach()
