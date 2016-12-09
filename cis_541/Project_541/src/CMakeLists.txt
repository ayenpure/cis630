cmake_minimum_required(VERSION 2.8.12.1)

PROJECT(project541)
SET(VTK_DIR /path/to/VTK6.0.0)

find_package(VTK REQUIRED)
include(${VTK_USE_FILE})

add_executable(project541 project541.cxx Triangle.h VectorOperations.h VectorOperations.cxx MatrixOperations.h)

if(VTK_LIBRARIES)
  target_link_libraries(project541 ${VTK_LIBRARIES})
else()
  target_link_libraries( vtkHybrid)
endif()
