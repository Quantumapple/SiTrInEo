#!/bin/bash

### Please use this shell script when you have an account on KNU Tier3 machine ###

## Set base environments
echo "Install CMSSW environments and setup"
source /cvmfs/cms.cern.ch/cmsset_default.sh
cmsrel CMSSW_9_3_7
cd ./CMSSW_9_3_7/src
eval `scramv1 runtime -sh`
cd ~ 
echo "Your CMSSW base is $CMSSW_BASE and version is $CMSSW_VERSION !"
echo ""

## CMake setup ## 
echo "CMake setup"
export CXX=`which g++`
export CC=`which gcc`
echo "Your gcc version is `gcc --version` !"
source /cvmfs/geant4.cern.ch/geant4/10.5/x86_64-slc6-gcc63-opt/CMake-setup.sh
source /cvmfs/sft.cern.ch/lcg/contrib/gcc/6.3/x86_64-slc6-gcc63-opt/setup.sh


## Geant4 setup ##
source /cvmfs/geant4.cern.ch/geant4/10.5/x86_64-slc6-gcc63-opt/bin/geant4.sh

