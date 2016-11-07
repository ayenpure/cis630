cmake_minimum_required(VERSION 2.8.12.1)

PROJECT(bezier)
find_package(OpenGL REQUIRED)
find_package(GLUT REQUIRED)
include_directories( ${OPENGL_INCLUDE_DIRS}  ${GLUT_INCLUDE_DIRS} )
add_executable(bezier bezier)
target_link_libraries(bezier ${OPENGL_LIBRARIES} ${GLUT_LIBRARY} )
