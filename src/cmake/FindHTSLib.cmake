# Find HTSLib
# Uses hint:
#   HTSLIB_ROOT
# Sets:
#   HTSLIB_FOUND
#   HTSLib_INCLUDE_DIR
#   HTSLib_LIBRARY
# Saves:
#   HTSLIB_ROOT
#   HTSLib_INCLUDE_DIR_CACHED
#   HTSLib_LIBRARY_CACHED

if(NOT "${HTSLIB_ROOT}" STREQUAL "${OLD_HTSLIB_ROOT}")
    message(STATUS "Detecting HTSLib: redetecing with new HTSLIB_ROOT=${HTSLIB_ROOT} (OLD_HTSLIB_ROOT=${OLD_HTSLIB_ROOT}).")
    unset(HTSLib_INCLUDE_DIRS_CACHED CACHE)
    unset(HTSLib_LIBRARIES_CACHED CACHE)
else()
    message(STATUS "Detecting HTSLib: HTSLIB_ROOT=${HTSLIB_ROOT} is not new; using cached paths.")
    message(STATUS "HTSLib_INCLUDE_DIRS_CACHED=${HTSLib_INCLUDE_DIRS_CACHED}")
    message(STATUS "HTSLib_LIBRARIES_CACHED=${HTSLib_LIBRARIES_CACHED}")
endif()
set(OLD_HTSLIB_ROOT ${HTSLIB_ROOT} CACHE INTERNAL "Last used value of HTSLIB_ROOT")

find_path(HTSLib_INCLUDE_DIR_CACHED htslib/vcf.h PATHS ${HTSLIB_ROOT} ${HTSLIB_ROOT}/include)
find_library(HTSLib_LIBRARY_CACHED hts PATHS ${HTSLIB_ROOT} ${HTSLIB_ROOT}/lib)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(HTSLib
    REQUIRED HTSLib_INCLUDE_DIR_CACHED HTSLib_LIBRARY_CACHED
    #"HTSLib library (https://github.com/samtools/htslib.git) not found. Specify location with -DHTSLIB_ROOT=<path>"
    )
mark_as_advanced(HTSLib_INCLUDE_DIR_CACHED HTSLib_LIBRARY_CACHED)

if(HTSLIB_FOUND)
    set(HTSLib_INCLUDE_DIR ${HTSLib_INCLUDE_DIR_CACHED})
    set(HTSLib_LIBRARY ${HTSLib_LIBRARY_CACHED})
endif()
