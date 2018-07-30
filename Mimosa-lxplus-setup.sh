#!/bin/bash

export GEANT4DIR="/afs/cern.ch/sw/lcg/external/geant4/10.0.p02/x86_64-slc6-gcc47-opt/"
###export GEANT4DIR="/afs/cern.ch/sw/lcg/external/geant4/10.2/x86_64-slc6-gcc48-opt/"

source /afs/cern.ch/sw/lcg/contrib/gcc/4.7.2/x86_64-slc6/setup.sh
###source /afs/cern.ch/sw/lcg/contrib/gcc/4.8.0/x86_64-slc6/setup.sh

source $GEANT4DIR/CMake-setup.sh
source $GEANT4DIR/../setup_g4datasets.sh

source $GEANT4DIR/share/Geant4-10.0.2/geant4make/geant4make.sh
###source $GEANT4DIR/share/Geant4-10.2.0/geant4make/geant4make.sh

export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/afs/cern.ch/sw/lcg/external/XercesC/3.1.1p2/x86_64-slc6-gcc47-opt/lib:/afs/cern.ch/sw/lcg/external/clhep/2.1.4.1/x86_64-slc6-gcc47-opt/lib:
export LIBRARY_PATH=$LIBRARY_PATH:/afs/cern.ch/sw/lcg/external/XercesC/3.1.1p2/x86_64-slc6-gcc47-opt/lib

source /afs/cern.ch/sw/lcg/app/releases/ROOT/5.34.24/x86_64-slc6-gcc47-opt/root/bin/thisroot.sh
###source /afs/cern.ch/sw/lcg/app/releases/ROOT/6.06.08/x86_64-slc6-gcc48-opt/root/bin/thisroot.sh

export CXX=`which g++`
export CC=`which gcc`
