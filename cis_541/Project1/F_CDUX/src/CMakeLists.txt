cmake_minimum_required(VERSION 2.8.12.1)

PROJECT(project1F)
SET(VTK_DIR /path/to/VTK6.0.0)

find_package(VTK REQUIRED)
include(${VTK_USE_FILE})

add_executable(qualitywork quality_work.cxx CameraPositions.h CameraPositions.cxx)
add_executable(qualitycamera quality_camera.cxx CameraPositions.h CameraPositions.cxx)
add_executable(qualitytriangles quality_triangle.cxx CameraPositions.h CameraPositions.cxx)

if(VTK_LIBRARIES)
  target_link_libraries(qualitywork ${VTK_LIBRARIES})
  target_link_libraries(qualitycamera ${VTK_LIBRARIES})
  target_link_libraries(qualitytriangles ${VTK_LIBRARIES})
else()
  target_link_libraries( vtkHybrid)
endif()
