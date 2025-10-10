# Try to locate the CFITSIO library
#
# Once done, this module defines
#  CFITSIO_FOUND - system has CFITSIO
#  CFITSIO_INCLUDE_DIRS - the CFITSIO include directory
#  CFITSIO_LIBRARIES - libraries needed to link CFITSIO
#  CFITSIO::CFITSIO - imported target for convenience

find_path(CFITSIO_INCLUDE_DIR
  NAMES fitsio.h
  PATH_SUFFIXES include
)

find_library(CFITSIO_LIBRARY
  NAMES cfitsio
)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(CFITSIO
  REQUIRED_VARS CFITSIO_LIBRARY CFITSIO_INCLUDE_DIR
)

if(CFITSIO_FOUND)
  set(CFITSIO_INCLUDE_DIRS "${CFITSIO_INCLUDE_DIR}")
  set(CFITSIO_LIBRARIES "${CFITSIO_LIBRARY}")
  if(NOT TARGET CFITSIO::CFITSIO)
    add_library(CFITSIO::CFITSIO UNKNOWN IMPORTED)
    set_target_properties(CFITSIO::CFITSIO PROPERTIES
      IMPORTED_LOCATION "${CFITSIO_LIBRARY}"
      INTERFACE_INCLUDE_DIRECTORIES "${CFITSIO_INCLUDE_DIR}"
    )
  endif()
endif()

mark_as_advanced(CFITSIO_INCLUDE_DIR CFITSIO_LIBRARY)