MimosaSimu framework
====================
MimosaSimu is a simulation framework for the SiTrInEo project which is built based on the Geant4. 
Unfortunately for now, using cmake tool to install MimosaSimu is only work on private machine.
But we will develop this to work on lxplus or other servers.

Pre-requisite
- Xerces-C, qt, qt-devel, qt4, libXi, libGL : for Geant4 graphic tools
- cmake > 3.9.0 (Necessary), ROOT: the latest version
- Geant4: 10.05 (If the version is > 10.05, you must use qt5 instead of qt4)

### Install pre-requisite libraries
The commands are based on Cern Centos7 (CC7).

```
yum -y update
yum -y groupinstall "Development Tools"
yum -y install expat.x86_64 expat-devel.x86_64
yum -y install qt qt-devel
yum -y install libXmu.x86_64 libXmu-devel.x86_64
yum -y install xerces-c-devel.x86_64
```

Commnads for Ubuntu (18.04, 20.04 LTS) machine:

```
apt-get update
apt-get upgrade
apt install build-essential
apt install libqt4-dev
apt install libxmu-dev
apt install libxerces-c-dev
apt install libcanberra-gtk-module libcanberra-gtk3-module
apt install cmake
```

### Install CMake and set ccmake configuration (CC7)
Download cmake 3.9.5 version.
```
wget https://cmake.org/files/v3.9/cmake-3.9.5-Linux-x86_64.tar.gz
tar -zxvf cmake-3.9.5-Linux-x86_64.tar.gz
```
Clone git repository.
```
git clone git@github.com:<github username>/SiTrInEo.git
cd SetupScripts
vi cmake_setup.sh # please use your prefer editor.
```
Change the path to your Cmake area.
```
source cmake_setup.sh
```

## How to install Geant4 

Download the latest Geant4 version, [Geant4 download link](https://github.com/Geant4/geant4/archive/v10.5.1.tar.gz)

```
export QT_QMAKE_EXECUTABLE=/usr/bin/qmake-qt4
```

Let's install Geant4 in user area
```
cd ~
mkdir Software && cd Software
mkdir Geant4 && cd Geant4
mkdir build src install
mv download.Geant4.source.tar.gz src // download.Geant4.source is just arbitrary name
cd src
tar -zxvf download.Geant4.source.tar.gz
cd ../build
ccmake ../src/geant4.10.05
```

Follow the configuration as below example:
![geant4configure](https://user-images.githubusercontent.com/35092541/53545505-07368b80-3b6d-11e9-9397-58262f1c127c.png)

CMAKE\_INSTALL\_PREFIX: path to where you want to install Geant4, recommended: /path/to/Software/Geant4/install  
Press 'c' until you can press 'g' to generate

```
make -j 4 && make install
source /path/to/Software/Geant4/install/share/Geant4-10.5.0/geant4make/geant4make.sh
```

### Install ROOT

Example with ROOT 6.16.00 (2019.01.28), [ROOT download link](https://root.cern/install/all_releases/)

After go the website via link, find a tab "Binary distribution".
Check your platform and compare with the uploaded files.
Here I recommend you use CC7. 
If you're platform is not listed, you should install using source distribution.
I assumed that you are using CC7.

```
wget https://root.cern/download/root_v6.16.00.Linux-centos7-x86_64-gcc4.8.tar.gz
tar -zxvf root_v6.16.00.Linux-centos7-x86_64-gcc4.8.tar.gz
cd root/bin
source thisroot.sh # tcsh for csh
```

Install root from source distribution:  
First please find this link to install pre-requisite libraries for ROOT framework.  
[ROOT pre-requisties](https://root.cern/install/dependencies/)

```
wget https://root.cern/download/root_v6.16.00.source.tar.gz
tar -zxvf root_v6.16.00.source.tar.gz
mkdir build src
mv root-6.16.00 src
cd build
ccmake ../src/root-6.16.00
```

Proceed the same process when you install Geant4.

After setting ROOT environments, move on to install SiTrInEo.

### Install MimosaSimu

Clone git repository
```
git clone git@github.com:<github username>/SiTrInEo.git
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
alias MimosaSimu=/path/to/the/build/directory/TestBeam_Geant4Simu_MagField/MimosaSimu
```
**Caution: You must check that before running MimosaSimu, you should set environment for ROOT and Geant4 together!**  

Then
```
cd ../../SiTrInEo/TestBeam_Geant4Simu_MagField
source /path/to/Software/Geant4/install/share/Geant4-10.5.0/geant4make/geant4make.sh  
source /path/to/root/bin/thisroot.sh  
MimosaSimu config/runXXX.cfg.Combine
```
