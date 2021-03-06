#----------------------------------------------------------------------------
# Project Name: TestBeam
# Description: Simulation library & executable for TestBeam_BField package
#----------------------------------------------------------------------------
# Setup the project and require CMake version above 2.8
#----------------------------------------------------------------------------
cmake_minimum_required(VERSION 2.8 FATAL_ERROR)
project(TestBeam_BField)

#----------------------------------------------------------------------------
# Add variables for the git revision
#----------------------------------------------------------------------------
# Get the current working branch
execute_process(
    COMMAND git rev-parse --abbrev-ref HEAD
    WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
    OUTPUT_VARIABLE TestBeam_BField_GIT_BRANCH
    OUTPUT_STRIP_TRAILING_WHITESPACE
    )

# Get the latest abbreviated commit hash of the working branch
execute_process(
    COMMAND git log -1 --format=%h
    WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
    OUTPUT_VARIABLE TestBeam_BField_GIT_COMMIT_HASH
    OUTPUT_STRIP_TRAILING_WHITESPACE
    )

#----------------------------------------------------------------------------
# Find ROOT and setup it
#----------------------------------------------------------------------------
list(APPEND CMAKE_PREFIX_PATH $ENV{ROOTSYS})
find_package(ROOT REQUIRED)
include(${ROOT_USE_FILE})

#----------------------------------------------------------------------------
# Find Geant4 and setup it
#----------------------------------------------------------------------------
option(WITH_GEANT4_UIVIS "Build example with Geant4 UI and Vis drivers" ON)
if(WITH_GEANT4_UIVIS)
    find_package(Geant4 REQUIRED ui_all vis_all)
else()
    find_package(Geant4 REQUIRED)
endif()
include(${Geant4_USE_FILE})
#For backward compatibility
set(CMAKE_CXX_STANDARD 11)

#----------------------------------------------------------------------------
# Copy the various files to build directory
#----------------------------------------------------------------------------
file(COPY ${PROJECT_SOURCE_DIR}/config ${PROJECT_SOURCE_DIR}/Authentic DESTINATION ${PROJECT_BINARY_DIR})

#----------------------------------------------------------------------------
# Locate sources and headers for this project
#----------------------------------------------------------------------------
file(GLOB MIMOSA_SIMU_EXEC_SOURCE ${PROJECT_SOURCE_DIR}/MimosaSimu.cc)
file(GLOB MIMOSA_SIMU_SOURCES ${PROJECT_SOURCE_DIR}/src/*.cc)
file(GLOB MIMOSA_SIMU_HEADERS ${PROJECT_SOURCE_DIR}/include/*.hh)

#----------------------------------------------------------------------------
# Add definitions for git commit hash
#----------------------------------------------------------------------------
set_source_files_properties(${MIMOSA_SIMU_EXEC_SOURCE} PROPERTIES
    COMPILE_DEFINITIONS "MIMOSA_SIMU_GIT_BRANCH=${MIMOSA_SIMU_GIT_BRANCH};MIMOSA_SIMU_GIT_COMMIT_HASH=${MIMOSA_SIMU_GIT_COMMIT_HASH}" )

#----------------------------------------------------------------------------
# Add the executables and libraries and link it to the Geant4 framework and
# MCObjs library
#----------------------------------------------------------------------------
add_executable(MimosaSimu ${MIMOSA_SIMU_EXEC_SOURCE} ${MIMOSA_SIMU_SOURCES} ${MIMOSA_SIMU_HEADERS})
target_link_libraries(MimosaSimu ${Geant4_LIBRARIES} ${ROOT_LIBRARIES})

#----------------------------------------------------------------------------
# Check dependencies for this project and set include directories and libraries
#----------------------------------------------------------------------------
if(NOT DEFINED MDIGIT)
    message(SEND_ERROR "This project cannot be built without MDIGIT. Please check about this.")
else()
    include_directories(${MDIGIT_INCLUDE_DIR})
    target_include_directories(MimosaSimu PUBLIC ${MDIGIT_INCLUDE_DIR})
    target_link_libraries(MimosaSimu MDIGIT)
endif()

#----------------------------------------------------------------------------
# Setup global and local include directories for this project
#----------------------------------------------------------------------------
include_directories(${PROJECT_SOURCE_DIR}/include ${Geant4_INCLUDE_DIRS} ${ROOT_INCLUDE_DIRS})
target_include_directories(MimosaSimu PUBLIC ${PROJECT_SOURCE_DIR}/include ${Geant4_INCLUDE_DIRS} ${ROOT_INCLUDE_DIRS})
target_compile_options(MimosaSimu PUBLIC -fdiagnostics-color=always)

#----------------------------------------------------------------------------
# Expose this public includes and library to other subprojects through cache
# variable.
#----------------------------------------------------------------------------
set(${PROJECT_NAME}_INCLUDE_DIR ${PROJECT_SOURCE_DIR}/include
    CACHE INTERNAL "${PROJECT_NAME}: Include Directories" FORCE)
set(${PROJECT_NAME}_LIBRARY ${PROJECT_BINARY_DIR}/lib${PROJECT_NAME}.so
    CACHE INTERNAL "${PROJECT_NAME}: Path to library" FORCE)
set(${PROJECT_NAME} "Exist!"
    CACHE INTERNAL "${PROJECT_NAME}: Variable for checking dependencies" FORCE)
