# MolFC

MolFC is a software package that can calculate Franck-Condon integrals
for non-linear molecular systems and simulate electronic spectra in the harmonic approximation including 
normal mode displacements and Duschinsky effects.

1. R. Borrelli, A. Capobianco and A. Peluso J. Phys. Chem. A 116 9934 (2012).
2. R. Borrelli and A. Peluso J. Chem. Phys. 129, 064116 (2008).
3. R. Borrelli and A. Peluso J. Chem. Phys. 125, 194308 (2006).
4. R. Borrelli and A. Peluso J. Chem. Phys. 119, 16 (2003).

##Downloading the code
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

