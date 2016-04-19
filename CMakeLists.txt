cmake_minimum_required(VERSION 3.5)
project(PBF_15462)

find_package(PkgConfig REQUIRED)

pkg_search_module(GLFW REQUIRED glfw3)

find_package(GLEW REQUIRED)

if (GLEW_FOUND)
    include_directories($(GLEW_INCLUDE_DIRS))
    link_libraries(${GLEW_LIBRARY})
endif()

IF(APPLE)
    INCLUDE_DIRECTORIES ( /System/Library/Frameworks )
    FIND_LIBRARY(COCOA_LIBRARY Cocoa)
    FIND_LIBRARY(GLUT_LIBRARY GLUT )
    FIND_LIBRARY(OpenGL_LIBRARY OpenGL )
    MARK_AS_ADVANCED (COCOA_LIBRARY
            GLUT_LIBRARY
            OpenGL_LIBRARY)
    SET(EXTRA_LIBS ${COCOA_LIBRARY} ${GLUT_LIBRARY} ${OpenGL_LIBRARY})
ENDIF (APPLE)

include_directories(${GLFW_INCLUDE_DIRS})

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")

set(SOURCE_FILES main.cpp ParticleSystem.cpp ParticleSystem.h Particle.cpp Particle.h shader.cpp shader.hpp sphere.h
        controls.hpp controls.cpp)

add_executable(PBF_15462 ${SOURCE_FILES})

target_link_libraries(PBF_15462 ${GLFW_LIBRARIES})

target_link_libraries(PBF_15462 ${EXTRA_LIBS})


