MimosaSimu framework
====================
MimosaSimu is a simulation framework for the SiTrInEo project which is built based on the Geant4. To run this, we assume that you already install Geant4 properly in your private machine. 
Unfortunately for now, using cmake tool to install MimosaSimu is only work on private machine.
But we will develop this to work on lxplus or other servers.

## How to install MimosaSimu 

Pre-requisite
- Xerces-C, qt, qt-devel, qt3, libXi, libGL : for Geant4 graphic tools
- cmake > 3.9.0 (Recommended), ROOT > 6.12.xx
- Geant4 > 4.10.2 (Recommended)

### Install ROOT using binary version

Example with ROOT 6.16.00 (2019.01.28), [ROOT download link](https://root.cern.ch/content/release-61600)

After go the website via link, find a tab "Binary distribution".
Check your platform and compare with the uploaded files.
Here I recommend you use Cern Centos7 (CC7). 
If you're platform is not listed, you should install using source distribution.
I assumed that you are using CC7.

```
wget https://root.cern/download/root_v6.16.00.Linux-centos7-x86_64-gcc4.8.tar.gz
tar -zxvf [downloaded_filename].tar
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
**You must check that before running MimosaSimu, you should set environment for ROOT and Geant4 together!**

Then
```
cd ../../SiTrInEo/TestBeam_Geant4Simu_MagField
source /path/to/Geant4/build/directory/geant4make.sh
MimosaSimu config/runXXX.cfg.Combine
```





