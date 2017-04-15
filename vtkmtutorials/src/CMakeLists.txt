cmake_minimum_required(VERSION 3.3)
project(pointelev)
project(vecex)
find_package(OpenGL REQUIRED)
find_package(GLUT REQUIRED)
find_package(VTKm REQUIRED
COMPONENTS Rendering
OPTIONAL_COMPONENTS TBB   )
add_executable(pointelev elevationfilter.cxx)
add_executable(vecex vecex.cxx)
target_include_directories(pointelev PRIVATE ${VTKm_INCLUDE_DIRS})
target_link_libraries(pointelev ${VTKm_LIBRARIES})
target_link_libraries(pointelev ${OPENGL_LIBRARIES} ${GLUT_LIBRARY} )
target_compile_options(pointelev PRIVATE ${VTKm_COMPILE_OPTIONS})

target_include_directories(vecex PRIVATE ${VTKm_INCLUDE_DIRS})
target_link_libraries(vecex ${VTKm_LIBRARIES})
target_link_libraries(vecex ${OPENGL_LIBRARIES} ${GLUT_LIBRARY} )
target_compile_options(vecex PRIVATE ${VTKm_COMPILE_OPTIONS})

