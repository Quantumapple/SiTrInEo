MimosaSimu framework
====================
MimosaSimu is a simulation framework for the SiTrInEo project which is built based on the Geant4. 
Unfortunately for now, using cmake tool to install MimosaSimu is only work on private machine.
But we will develop this to work on lxplus or other servers.

Pre-requisite
- Xerces-C, qt, qt-devel, qt3, libXi, libGL : for Geant4 graphic tools
- cmake > 3.9.0 (Recommended), ROOT: the latest version
- Geant4: the latest version

### Install pre-requisite libraries
The commands are based on Cern Centos7 (CC7).

```
yum -y update
yum -y groupinstall "Development Tools"
yum -y install expat.x86_64 expat-devel.x86_64
yum -y install qt qt-devel qt-config
yum -y install libXmu.x86_64 libXmu-devel.x86_64
yum -y install xerces-c
```

### Install CMake and set ccmake configuration
```
wget https://cmake.org/files/v3.9/cmake-3.9.5-Linux-x86_64.tar.gz
tar -zxvf cmake-3.9.5-Linux-x86_64.tar.gz
cd cmake-3.9.5-Linux-x86_64/bin
export ccmake=${PWD}/ccmake
```

## How to install Geant4 

Download the latest Geant4 version, [Geant4 download link](https://geant4.web.cern.ch/support/download)

```
export QT_QMAKE_EXECUTABLE=/usr/bin/qmake-qt4
```
![Exmample](./InstallExample/Geant4Configure.png)

### Install ROOT using binary version

Example with ROOT 6.16.00 (2019.01.28), [ROOT download link](https://root.cern.ch/content/release-61600)

After go the website via link, find a tab "Binary distribution".
Check your platform and compare with the uploaded files.
Here I recommend you use CC7. 
If you're platform is not listed, you should install using source distribution.
I assumed that you are using CC7.

```
wget https://root.cern/download/root_v6.16.00.Linux-centos7-x86_64-gcc4.8.tar.gz
tar -zxvf root_v6.16.00.Linux-centos7-x86_64-gcc4.8.tar.gz
cd root/bin
source thisroot.sh # -csh for csh
```

After setting ROOT environments, move on to install SiTrInEo.

### Install MimosaSimu

Clone "Develop" branch in git repository
```
git clone -b Develop git@github.com:Quantumapple/SiTrInEo.git
mkdir build
cd build
ccmake ../SiTrInEo
```
After ccmake command, you can see the new screen.
Just press "c" to cofigure and press agian.
If you can see "g" to generate at below, press g.
Then
```
make -j 2(4) # depend on your machine
```

## How to run MimosaSimu framework
I recommend you to setup your new command.
For example
```
export MimosaSimu=/path/to/the/build/directory/TestBeam_Geant4Simu_MagField/MimosaSimu
```
**Caution: You must check that before running MimosaSimu, you should set environment for ROOT and Geant4 together!**

Then
```
cd ../../SiTrInEo/TestBeam_Geant4Simu_MagField
source /path/to/Geant4/build/directory/geant4make.sh
MimosaSimu config/runXXX.cfg.Combine
```





