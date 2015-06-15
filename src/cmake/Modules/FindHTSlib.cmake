# Find HTSlib
# Uses hint:
#   HTSLIB_ROOT
# Sets:
#   HTSLIB_FOUND
#   HTSlib_INCLUDE_DIR
#   HTSlib_LIBRARY
# Saves:
#   HTSLIB_ROOT
#   HTSlib_INCLUDE_DIR_CACHED
#   HTSlib_LIBRARY_CACHED


if(NOT HTSlib_INCLUDE_DIR_CACHED OR NOT HTSlib_LIBRARY_CACHED)
    set(HTSLIB_ROOT "$ENV{HTSLIB_ROOT}" CACHE PATH "Path to HTSlib")

    find_path(HTSlib_INCLUDE_DIR_CACHED htslib/vcf.h PATHS ${HTSLIB_ROOT} ${HTSLIB_ROOT}/include)
    find_library(HTSlib_LIBRARY_CACHED hts PATHS ${HTSLIB_ROOT} ${HTSLIB_ROOT}/lib)

    include(FindPackageHandleStandardArgs)
    find_package_handle_standard_args(HTSlib
        "HTSlib library (https://github.com/samtools/htslib.git) not found. Specify location with -DHTSLIB_ROOT=<path>"
        HTSlib_LIBRARY_CACHED HTSlib_INCLUDE_DIR_CACHED)
    mark_as_advanced(HTSlib_INCLUDE_DIR_CACHED HTSlib_LIBRARY_CACHED)
else()
    message(STATUS "Using HTSlib: ${HTSlib_LIBRARY_CACHED}")
    set(HTSLIB_FOUND TRUE)
endif()

if(HTSLIB_FOUND)
    set(HTSlib_INCLUDE_DIR ${HTSlib_INCLUDE_DIR_CACHED})
    set(HTSlib_LIBRARY ${HTSlib_LIBRARY_CACHED})
endif()
