#----------------------------------------------------------------------------
# Project Name: Mimosa(?) Framework
# Description: A framework for MIMOSA simulation and digitizer
# Author: BaseHardware(basehw); 2018/12/19
#----------------------------------------------------------------------------
# Setup the project and require CMake version above 2.8
#----------------------------------------------------------------------------
cmake_minimum_required(VERSION 2.8 FATAL_ERROR)
project(MIMOSA_FRWK)

#----------------------------------------------------------------------------
# Add variables for the git revision
#----------------------------------------------------------------------------
# Get the current working branch
execute_process(
  COMMAND git rev-parse --abbrev-ref HEAD
  WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
  OUTPUT_VARIABLE Master_GIT_BRANCH
  OUTPUT_STRIP_TRAILING_WHITESPACE
)

# Get the latest abbreviated commit hash of the working branch
execute_process(
  COMMAND git log -1 --format=%h
  WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
  OUTPUT_VARIABLE Master_GIT_COMMIT_HASH
  OUTPUT_STRIP_TRAILING_WHITESPACE
)

#----------------------------------------------------------------------------
# Add subdirectories
#----------------------------------------------------------------------------
add_subdirectory(MIMOSA_DIGITIZER)
#add_subdirectory(TestBeam_Geant4Simu_MagField)
