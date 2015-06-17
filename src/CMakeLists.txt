message(STATUS "Processing: main")

cmake_minimum_required(VERSION 2.8.12)
set(CMAKE_VERBOSE_MAKEFILE TRUE CACHE BOOL "Enable verbose make output.") #TODO: remove
project(BAM_PHASE_SPLIT C CXX)
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/Modules")
include(CheckLibraryExists)
include(CheckIncludeFileCXX)
include(ExternalProject)

# check main file exists; if not, assume the source dir is wrong
if(NOT EXISTS "${PROJECT_SOURCE_DIR}/bam-phase-split.cpp")
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
set(OPTIONAL_SUBMODULE_DIR ${PROJECT_SOURCE_DIR}/optional CACHE INTERNAL "Directory for downloading and building optional submodules")
set(OPTIONAL_SUBMODULE_PREFIX ${PROJECT_SOURCE_DIR}/optional-prefix CACHE INTERNAL "Directory for installing optional submodules")

set(HEADER_SUBDIRS
    ${PROJECT_SOURCE_DIR}
    CACHE INTERNAL "Subdirectories containing header files")
#set(SUBDIRS
#    CACHE INTERNAL "Subdirectories to descend into")

### Resolve external dependencies
#
# prefer static libraries
#set(CMAKE_FIND_LIBRARY_SUFFIXES .a .so)

# bake-in RPATH to prevent library search problems
set(CMAKE_INSTALL_RPATH_USE_LINK_PATH true)
set(CMAKE_BUILD_WITH_INSTALL_RPATH true)

# find zlib
# => ZLIB_INCLUDE_DIRS, ZLIB_LIBRARIES
find_package(ZLIB REQUIRED)

# check for HTSlib
# => HTSlib_INCLUDE_DIR, HTSlib_LIBRARY
if(BUILD_HTSLIB)
    message(STATUS "Building HTSlib in: ${PROJECT_SOURCE_DIR}/htslib")
    ExternalProject_Add(HTSlib
        SOURCE_DIR ${PROJECT_SOURCE_DIR}/htslib
        CONFIGURE_COMMAND ""
        BUILD_IN_SOURCE 1
        BUILD_COMMAND make lib-static
        INSTALL_COMMAND ""
        )
    set(HTSlib_INCLUDE_DIR ${PROJECT_SOURCE_DIR}/htslib)
    set(HTSlib_LIBRARY ${PROJECT_SOURCE_DIR}/htslib/${CMAKE_STATIC_LIBRARY_PREFIX}hts${CMAKE_STATIC_LIBRARY_SUFFIX})
else()
    find_package(HTSlib REQUIRED)
endif()

# check for Boost 1.57
# => Boost_INCLUDE_DIRS, Boost_<COMPONENT>_LIBRARY
if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang" AND "${Boost_COMPILER}" STREQUAL "")
    message(WARNING
        "The cmake FindBoost package does not properly detect the compiler name when using CLang. Boost libraries built with 'versioned' layout (compiler name in library name, e.g. 'libboost_program_options-clang35-mt-1_57.so') might be incorrectly detected. To find such libraries, give the compiler name on the cmake command line, e.g. 'cmake -DBoost_COMPILER=\"-clang35\" [...]'. To check which libraries have been detected, use: 'grep libboost_program_options CMakeCache.txt'")
endif()
set(BOOST_ROOT "$ENV{BOOST_ROOT}" CACHE PATH "Path to Boost")
# Use static libraries to avoid LD_RUN_PATH issues
#set(Boost_USE_STATIC_LIBS true)
#find_package(Boost 1.57.0 EXACT REQUIRED COMPONENTS iostreams program_options)
find_package(Boost 1.57.0 EXACT REQUIRED)
if(NOT Boost_FOUND)
    message(FATAL_ERROR "Boost 1.57.0 not found")
endif()

### Build configuration file
#
# check for clock_gettime()
# => HAVE_CLOCK_GETTIME
check_library_exists(rt clock_gettime "" HAVE_CLOCK_GETTIME)

# generate configuration file
configure_file(config.h.cmake.in config.h)

### Prepare default compilation flags
#
# headers
include_directories(SYSTEM
    ${ZLIB_INCLUDE_DIRS}
    ${Boost_INCLUDE_DIRS}
    ${HTSlib_INCLUDE_DIR}
    tclap/include
    zstr/src)
include_directories(
    ${PROJECT_SOURCE_DIR}
    ${PROJECT_BINARY_DIR})
get_directory_property(include_directories INCLUDE_DIRECTORIES)
message(STATUS "INCLUDE_DIRECTORIES='${include_directories}'")

### general compile flags
#
set(EXTRA_FLAGS "-std=c++11 -Wall -Wextra -pedantic -pthread")

# compiler-specific flags
if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
    set(EXTRA_FLAGS "${EXTRA_FLAGS} -fmax-errors=1")
elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
    set(EXTRA_FLAGS "${EXTRA_FLAGS} -ferror-limit=1")
endif()

# consolidate flags
set(CMAKE_C_FLAGS "$ENV{CFLAGS} ${EXTRA_FLAGS}")
set(CMAKE_CXX_FLAGS "$ENV{CXXFLAGS} ${EXTRA_FLAGS}")
message(STATUS "CMAKE_CXX_FLAGS='${CMAKE_CXX_FLAGS}'")

### build-specific flags
#
# compile flags for Debug build
set(CMAKE_CXX_FLAGS_DEBUG "-O0 -g3 -ggdb -fno-inline -fno-eliminate-unused-debug-types")
set(CMAKE_C_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG}")
# compile flags for Test build
set(CMAKE_CXX_FLAGS_TEST "-O3 -g3")
set(CMAKE_C_FLAGS_TEST "${CMAKE_CXX_FLAGS_TEST}")
set(CMAKE_CXX_FLAGS_TEST_O2 "-O2 -g3")
set(CMAKE_C_FLAGS_TEST_O2 "${CMAKE_CXX_FLAGS_TEST_O2}")
set(CMAKE_CXX_FLAGS_TEST_O1 "-O1 -g3")
set(CMAKE_C_FLAGS_TEST_O1 "${CMAKE_CXX_FLAGS_TEST_O1}")
# compile flags for GProf build
set(CMAKE_CXX_FLAGS_GPROF "-O3 -g3 -pg")
set(CMAKE_C_FLAGS_GPROF "${CMAKE_CXX_FLAGS_GPROF}")
# compile flags for Release build
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -DNDEBUG -DDISABLE_ASSERTS")
set(CMAKE_C_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE}")
# compile flags for GProfRel build
set(CMAKE_CXX_FLAGS_GPROFREL "-O3 -DNDEBUG -DDISABLE_ASSERTS -pg")
set(CMAKE_C_FLAGS_GPROFREL "${CMAKE_CXX_FLAGS_GPROF}")

# link flags
#set(CMAKE_EXE_LINKER_FLAGS "-Wl,-rpath=$ENV{LD_RUN_PATH} -Wl,--as-needed -static-libgcc -static-libstdc++")
set(CMAKE_EXE_LINKER_FLAGS "-Wl,--as-needed")
message(STATUS "CMAKE_EXE_LINKER_FLAGS='${CMAKE_EXE_LINKER_FLAGS}'")

### Add target that updates version
#
add_custom_target(version ALL
    ${CMAKE_COMMAND} -DSRC_ROOT=${PROJECT_SOURCE_DIR} -DBIN_ROOT=${PROJECT_BINARY_DIR} -P ${PROJECT_SOURCE_DIR}/cmake/update-version.cmake)

### Enable testing targets
#
#enable_testing()

### Descend into subdirectories
#
#foreach(dir ${SUBDIRS})
#    add_subdirectory(${dir})
#endforeach()

add_executable(bam-phase-split bam-phase-split.cpp)
target_link_libraries(bam-phase-split
    ${HTSlib_LIBRARY}
    ${ZLIB_LIBRARIES}
)
add_dependencies(bam-phase-split version)

file(COPY add-parent-phasing DESTINATION .)
