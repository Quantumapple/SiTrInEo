#! /bin/sh -f

clear

export DIGI_DIR="/afs/cern.ch/work/j/jongho/Strasbourg/MIMOSA_DIGITIZER"
export GEANT4_BUILD_DIR="/afs/cern.ch/sw/lcg/external/geant4/10.0.p02/x86_64-slc6-gcc47-opt/share/Geant4-10.0.2/geant4make"

echo
echo
[ -f geant4make.sh ] && echo "File geant4make.sh already exist." || echo "File geant4make.sh doesn't exist, creating symbolic link to $GEANT4_BUILD_DIR/geant4make.sh"
[ -f geant4make.sh ] && echo "" || ln -s $GEANT4_BUILD_DIR/geant4make.sh  geant4make.sh
ls -tlrh geant4make.sh
source geant4make.sh
echo "Setting DIGI_DIR = $DIGI_DIR"
echo
echo

