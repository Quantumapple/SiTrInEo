#! /bin/sh -f

clear

#export DIGI_DIR="/Users/jeromeb/rawcmos/SVNMAC/CMOS/MIMOSA_DIGITIZER/trunk/"
export DIGI_DIR="/home/jongho/Analysis/Strasbourg/MIMOSA_DIGITIZER/trunk/"
#export GEANT4_BUILD_DIR="/home/aperez/Geant4/geant4.10.01.p02-build"
export GEANT4_BUILD_DIR="/home/jongho/Software/geant4/build"
echo "Path to build directory $GEANT4_BUILD_DIR "

#export GEANT4_BUILD_DIR="/usr/local/share/Geant4-10.4.2"

echo
echo
[ -f geant4make.sh ] && echo "File geant4make.sh already exist." || echo "File geant4make.sh doesn't exist, creating symbolic link to $GEANT4_BUILD_DIR/geant4make.sh"
[ -f geant4make.sh ] && echo "" || ln -s $GEANT4_BUILD_DIR/geant4make.sh  geant4make.sh
ls -tlrh geant4make.sh
source geant4make.sh
echo "Setting DIGI_DIR = $DIGI_DIR"
echo
echo

