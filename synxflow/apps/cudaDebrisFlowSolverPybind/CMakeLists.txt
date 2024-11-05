cmake_minimum_required(VERSION 3.11)

set(EXE debris)
set(SRCS debris_solver_pybind.cc)
set(LIB lib_debris)
set(LIB_SRCS cuda_debris_flow_solver.cu)

#set (CUDA_NVCC_FLAGS "-arch=sm_20 -Xptxas=v")
#set (CUDA_SEPARABLE_COMPILATION ON)

file(GLOB_RECURSE gc_headers "${PROJECT_SOURCE_DIR}/synxflow/lib/include/*.h")

set (GC_INCLUDE_DIRS "")
foreach (_header_file ${gc_headers})
	get_filename_component(_dir ${_header_file} DIRECTORY)
	list (APPEND GC_INCLUDE_DIRS ${_dir})
endforeach()
list(REMOVE_DUPLICATES GC_INCLUDE_DIRS)

include_directories(${GC_INCLUDE_DIRS})
include_directories(${netCDF_INCLUDE_DIR})
include_directories("${pybind11_sources_SOURCE_DIR}/include")
include_directories(${PYTHON_INCLUDE_DIRS})
include_directories(${PROJECT_SOURCE_DIR}/synxflow/apps/cudaDebrisFlowSolverPybind)
link_directories(${netCDF_LIB_DIR})

cuda_add_library(${LIB} STATIC ${LIB_SRCS})
add_dependencies(${LIB} lib cuda_lib)

#set (CUDA_NVCC_FLAGS "-arch=sm_35")
set(EXTERNAL_INSTALL_DIR ${CMAKE_BINARY_DIR}/external)
include_directories(${EXTERNAL_INSTALL_DIR}/include)

if(EXISTS "${EXTERNAL_INSTALL_DIR}/lib64")
    set(LIB_DIR "${EXTERNAL_INSTALL_DIR}/lib64")
elseif(EXISTS "${EXTERNAL_INSTALL_DIR}/lib")
    set(LIB_DIR "${EXTERNAL_INSTALL_DIR}/lib")
else()
    message(FATAL_ERROR "Neither lib nor lib64 directory found in ${EXTERNAL_INSTALL_DIR}")
endif()


link_directories(${LIB_DIR})
link_directories(${netCDF_LIB_DIR})

message(STATUS "External install directory: ${EXTERNAL_INSTALL_DIR}")

target_link_libraries(${LIB} ${LIB_DIR}/libnetcdf-cxx4.a netcdf lib cuda_lib)

add_library(${EXE} MODULE ${SRCS})

target_link_libraries(${EXE} ${LIB})
target_link_libraries(${EXE} pybind11::module)

set_target_properties(${EXE} PROPERTIES PREFIX "${PYTHON_MODULE_PREFIX}"
                                        SUFFIX "${PYTHON_MODULE_EXTENSION}")
#set(EXECUTABLE_OUTPUT_PATH ${PROJECT_SOURCE_DIR}/release/bin)
#set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/release/lib)

install(
  TARGETS
    ${debris}
  LIBRARY
    DESTINATION ${CMAKE_LIBRARY_OUTPUT_DIRECTORY}/synxflow
  )
