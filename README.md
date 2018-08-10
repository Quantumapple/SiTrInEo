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

<li> **RECOMMAND**: option GEANT4\_BUILD\_MULTITHREADED set "ON" </li>    
<li> Open Mimosa-private-setup.sh and change path to fit your own machine</li>    
<li> Then, type command:</li>
<blockquote>
<p> source Mimosa-private-setup.sh</p>
</blockquote>

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


