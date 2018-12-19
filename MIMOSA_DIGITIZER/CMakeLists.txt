#----------------------------------------------------------------------------
# Project Name: MDIGIT
# Description: Classes for Amoresim, Cupsim framework, which are used to store
#              various information.
#----------------------------------------------------------------------------
# Setup the project and require CMake version above 2.8
#----------------------------------------------------------------------------
cmake_minimum_required(VERSION 2.8 FATAL_ERROR)
project(MDIGIT)

#----------------------------------------------------------------------------
# Find ROOT and setup it
#----------------------------------------------------------------------------
list(APPEND CMAKE_PREFIX_PATH $ENV{ROOTSYS})
find_package(ROOT REQUIRED)
if(NOT DEFINED ROOT_GENERATE_DICTIONARY)
  include(${ROOT_USE_FILE})
endif()

#----------------------------------------------------------------------------
# Locate sources and headers for this project
#----------------------------------------------------------------------------
file(GLOB MDIGIT_SOURCES ${PROJECT_SOURCE_DIR}/src/*.cc)
file(GLOB MDIGIT_HEADERS ${PROJECT_SOURCE_DIR}/include/*.hh)

#----------------------------------------------------------------------------
# Add the libraries and link it to the ROOT framework and MDIGIT library
#----------------------------------------------------------------------------
add_library(MDIGIT SHARED ${MDIGIT_SOURCES})
target_link_libraries(MDIGIT ${ROOT_LIBRARIES})

#----------------------------------------------------------------------------
# Setup global and local include directories for this project
#----------------------------------------------------------------------------
include_directories(${PROJECT_SOURCE_DIR})
include_directories(${ROOT_INCLUDE_DIRS})
target_include_directories(MDIGIT PUBLIC ${PROJECT_SOURCE_DIR})
target_include_directories(MDIGIT PUBLIC ${ROOT_INCLUDE_DIRS})

#----------------------------------------------------------------------------
# Expose this public includes to other subprojects through cache variable.
#----------------------------------------------------------------------------
set(${PROJECT_NAME}_INCLUDE_DIR ${PROJECT_SOURCE_DIR}/include
  CACHE INTERNAL "${PROJECT_NAME}: Include Directories" FORCE)
set(${PROJECT_NAME}_LIBRARY ${PROJECT_BINARY_DIR}/lib${PROJECT_NAME}.so
  CACHE INTERNAL "${PROJECT_NAME}: Path to library" FORCE)
set(${PROJECT_NAME} "Exist!"
  CACHE INTERNAL "${PROJECT_NAME}: Variable for checking dependencies" FORCE)

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

