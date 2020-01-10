cmake_minimum_required(VERSION 3.0 FATAL_ERROR)

project(diffusion_1D LANGUAGES Fortran)


if (CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -fdefault-real-8 -ffpe-summary=none -ffree-line-length-none -Ofast -fopenmp")
elseif (CMAKE_Fortran_COMPILER_ID STREQUAL "PGI")
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -r8 -fast -fPIC")
endif ()

if (DEFINED ENV{NETCDF} AND (NOT DEFINED ENV{NETCDF_ROOT}))
  set(ENV{NETCDF_ROOT} $ENV{NETCDF})
endif ()
if (DEFINED ENV{NETCDF_ROOT})
  include_directories("$ENV{NETCDF_ROOT}/include")
  link_directories("$ENV{NETCDF_ROOT}/lib")
else ()
  find_package(PkgConfig)
  if (PKG_CONFIG_FOUND)
    pkg_search_module(NETCDF REQUIRED netcdf)
    if (NETCDF_FOUND)
      include_directories(${NETCDF_INCLUDEDIR})
      link_directories(${NETCDF_LIBDIR})
    else ()
      message(FATAL_ERROR "Pkg-config could not find netcdf library!")
    endif ()
  else ()
    message(FATAL_ERROR "Unable to find pkg-config library!")
  endif ()
endif ()

# add_subdirectory(lib/datetime)

# add_subdirectory(lib/container)
# include_directories(${CMAKE_BINARY_DIR}/fortran_container)

# add_subdirectory(lib/fft)

set(source_files
  src/main.F90
  src/namelist_mod.F90
  src/mesh_mod.F90
  src/state_mod.F90
  src/scheme_mod.F90
  src/history_mod.F90
)

set(CMAKE_Fortran_MODULE_DIRECTORY "${CMAKE_BINARY_DIR}/.mods")

add_executable(diffusion_1d.exe ${source_files})
target_link_libraries(diffusion_1d.exe ${MPI_Fortran_LIBRARIES} netcdf netcdff)

