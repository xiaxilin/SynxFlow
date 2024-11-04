cmake_minimum_required(VERSION 3.11)

include(ExternalProject)

# Specify install directory for the external projects
set(EXTERNAL_INSTALL_DIR ${CMAKE_BINARY_DIR}/external)

get_filename_component(NETCDF_PARENT_DIR ${netCDF_INCLUDE_DIR} DIRECTORY)

# External project for netCDF-CXX4
ExternalProject_Add(netcdf-cxx4
    GIT_REPOSITORY https://github.com/Unidata/netcdf-cxx4.git
    GIT_TAG 9ecd295a8f9a5e167daa313022ced453b4731f04
    CMAKE_ARGS -DBUILD_SHARED_LIBS=OFF
               -DCMAKE_POSITION_INDEPENDENT_CODE=ON
               -DHDF5_ROOT=${NETCDF_PARENT_DIR}
               -DHDF5_INCLUDE_DIRS=${netCDF_INCLUDE_DIR}
               -DHDF5_LIBRARIES=${netCDF_LIB_DIR}/libhdf5.so
               -DCMAKE_C_FLAGS=-I${netCDF_INCLUDE_DIR}
    BUILD_COMMAND ${CMAKE_COMMAND} --build . --config Release --verbose
    INSTALL_COMMAND ${CMAKE_COMMAND} --install . --prefix ${EXTERNAL_INSTALL_DIR}
    STEP_TARGETS build install
)

# Include directories for the headers
file(GLOB_RECURSE gc_headers "${CMAKE_CURRENT_SOURCE_DIR}/include/*.h")
set(GC_INCLUDE_DIRS "")
foreach (_header_file ${gc_headers})
    get_filename_component(_dir ${_header_file} DIRECTORY)
    list(APPEND GC_INCLUDE_DIRS ${_dir})
endforeach()
list(REMOVE_DUPLICATES GC_INCLUDE_DIRS)

# Gather source files
file(GLOB_RECURSE gc_srcs "${CMAKE_CURRENT_SOURCE_DIR}/src/*.cc")
file(GLOB_RECURSE gc_cuda_srcs "${CMAKE_CURRENT_SOURCE_DIR}/src/*.cu")

if (TARGET lib)
    message("lib is already defined")
else()
    include_directories(${GC_INCLUDE_DIRS})
    include_directories(${EXTERNAL_INSTALL_DIR}/include)

    link_directories(${EXTERNAL_INSTALL_DIR}/lib)
    link_directories(${netCDF_LIB_DIR})

    # Add the main static library
    add_library(lib STATIC ${gc_srcs})
    add_dependencies(lib netcdf-cxx4)  # Ensure netCDF-CXX4 is built before linking
    target_link_libraries(lib ${EXTERNAL_INSTALL_DIR}/lib/libnetcdf-cxx4.a netcdf)

    # Add CUDA static library
    cuda_add_library(cuda_lib STATIC ${gc_cuda_srcs})
    add_dependencies(cuda_lib netcdf-cxx4)  # Ensure netCDF-CXX4 is built before linking
    target_link_libraries(cuda_lib ${EXTERNAL_INSTALL_DIR}/lib/libnetcdf-cxx4.a netcdf)

endif()
