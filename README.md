How to install MimosaSimu 
=================================

There are two options to use MimosaSimu. <br>

Install
-------

### Lxplus
Run in Lxplus machine. To set up environment, type this: 
<blockquote>
<p> source Mimosa-lxplus-setup.sh</p>
</blockquote>

### Private machine
To run in your own machine. You have to install pre-requisite libraries and programs. 
<li> Xerces-C(libxerces-c), qt, qt-devel, qt3, libXi, libGL </li>
<li> CLHEP > 2.1.x.x, cmake > 3.9.0(Necessary!), ROOT 5.34(Not available with ROOT 6) </li>
<li> Install Geant4 (4.9.x or 4.10.x) </li> 
<li> Use 'ccmake' command and set install options as picture </li>

![Option](./image/screenshot.png)

<li> Open Mimosa-private-setup.sh and change path to fit your own machine</li>    
<li> Then, type command:</li>
<blockquote>
<p> source Mimosa-private-setup.sh</p>
</blockquote>

### KNU Tier3 server
Pre-requisite library: Xerces-C <br>
<blockquote>
<p> wget https://archive.apache.org/dist/xerces/c/3/sources/xerces-c-3.1.1.tar.gz </p>
<p> tar -zxvf xerces-c-3.1.1.tar.gz </p>
</blockquote>
Move to xerces-c directory, then
<blockquote>
<p> ./configure --prefix=path/to/build/directory </p>
<p> make -j3 </p>
<p> make install </p>
</blockquote>

Open ~/.bashrc and set environments
```export SCRAM_ARCH=slc6_amd64_gcc480''' 



### Common
<li> Below commands are the same to both Lxplus and private machine </li>
<blockquote>
<p> cd MIMOSA\_DIGITIZER/trunk</p>
<p> source Complie.sh</p>
<p> cd ../../TestBeam\_Geant4Simu/trunk</p>
<p> source Settings.sh</p>
<p> source Complie.sh</p>
</blockquote>

<li> Follow README in TestBeam\_Geant4Simu/trunk</li>


