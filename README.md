How to install and run MimosaSimu 
=================================

There are two options to use MimosaSimu. <br>

Install
-------

## Lxplus
Run in Lxplus machine. To set up environment, type this: 
<blockquote>
<p> source Mimosa-lxplus-setup.sh</p>
</blockquote>

## Private machine
To run in your own machine. You have to install pre-requisite libraries and programs. 
<li> Xerces-C(libxerces-c), qt, qt-devel, qt3, libXi, libGL </li>
<li> CLHEP > 2.1.x.x, cmake > 3.9.0(Necessary!), ROOT 5.34(Not available with ROOT 6) </li>
<li> Install Geant4 (4.9.x or 4.10.x) </li> 
<li> Use 'ccmake' command and set install options as picture </li>

![Option](./image/screenshot.png)

<li> Open Mimosa-private-setup.sh and change path to fit your own machine</li>    
<li> Then,type command:</li>
<blockquote>
<p> source Mimosa-private-setup.sh</p>
</blockquote>

Run
---

