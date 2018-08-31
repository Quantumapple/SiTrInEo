#!/bin/bash

### Basic environment setup ####
export CXX=`which g++`
export CC=`which gcc`
export GEANT4DIR="/usr/local/Geant4"

#### Geant4 environment setup ####
source $GEANT4DIR/share/Geant4-10.4.2/geant4make/geant4make.sh
##source $GEANT4DIR/../setup_g4datasets.sh

#### Library path setup ####
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/usr/lib64
export LD_LIBRARY_PATH=/usr/local/lib:$LD_LIBRARY_PATH
export LIBRARY_PATH=$LIBRARY_PATH:/usr/lib64

#### Root 5.34 version setup ####
source /home/jongho/Software/root_version5/root/bin/thisroot.sh 

#### Qt environments setup ####
export QTDIR="/usr/lib64/qt-3.3"
export QTINC="/usr/lib64/qt-3.3/include"
export QTLIB="/usr/lib64/qt-3.3/lib"

#### CLHEP environment setup ####
export CLHEP_BASE_DIR="/usr/local"
export CLHEP_INCLUDE_DIR="/usr/local/include"
export CLHEP_LIB_DIR="/usr/local/lib"
