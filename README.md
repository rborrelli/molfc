# MolFC

MolFC is a software package which can calculate Franck-Condon integrals
for non-linear molecular systems and which allows to study the
qunatum dynamical aspects on non-stationary molecular states
by using the methods described in the papers:  

1) R. Borrelli and A. Peluso J. Chem. Phys. 129, 064116 (2008).
2) R. Borrelli and A. Peluso J. Chem. Phys. 125, 194308 (2006).
3) R. Borrelli and A. Peluso J. Chem. Phys. 119, 16 (2003).

##Downloading the code
This installation works with git submodules, so you should be sure, that you got them all right.
The best way is to work with this repository is to use git with a version >= 1.6.5.
Then, to clone the repository, you can simply run
```
git clone --recursive git://github.com/rborrelli/molfc.git

```
To update to the latest version, run
```
git pull
```
##Prerequisites

**It is highly recommended** that you use the Intel Fortran Compiler releases 15.0 or laters.
**Gfortran compatibility is not guaranteed at the moment.

##Installing the package
The installation of the package is done via cmake.
Follow the instructions in the INSTALL file.

