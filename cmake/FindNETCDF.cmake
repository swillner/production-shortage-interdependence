if(NETCDF_INCLUDE_DIR AND NETCDF_LIBRARY)
  set(NETCDF_FIND_QUIETLY TRUE)
endif()

find_path(NETCDF_INCLUDE_DIR NAMES netcdf.h)
find_library(NETCDF_LIBRARY NAMES netcdf)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(NETCDF DEFAULT_MSG NETCDF_INCLUDE_DIR NETCDF_LIBRARY)

mark_as_advanced(NETCDF_INCLUDE_DIR NETCDF_LIBRARY)

set(NETCDF_LIBRARIES ${NETCDF_LIBRARY} )
set(NETCDF_INCLUDE_DIRS ${NETCDF_INCLUDE_DIR} )

if(NETCDF_FOUND AND NOT TARGET netcdf)
  add_library(netcdf UNKNOWN IMPORTED)
  set_target_properties(netcdf PROPERTIES INTERFACE_INCLUDE_DIRECTORIES "${NETCDF_INCLUDE_DIRS}")
  set_target_properties(netcdf PROPERTIES IMPORTED_LINK_INTERFACE_LANGUAGES "C")
  set_target_properties(netcdf PROPERTIES IMPORTED_LOCATION "${NETCDF_LIBRARIES}")
endif()
