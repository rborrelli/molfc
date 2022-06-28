# MolFC

MolFC is a software package that can calculate Franck-Condon integrals
for non-linear molecular systems and simulate electronic spectra in the harmonic approximation including 
normal mode displacements and Duschinsky effects.

1. R. Borrelli, A. Capobianco and A. Peluso J. Phys. Chem. A 116 9934 (2012).
3. R. Borrelli and A. Peluso J. Chem. Phys. 125, 194308 (2006).
2. R. Borrelli and A. Peluso J. Chem. Phys. 129, 064116 (2008).
4. R. Borrelli and A. Peluso J. Chem. Phys. 119, 16 (2003).
5. R. Borrelli, A. Capobianco, A. Peluso, Can. J. Chem., 91, 495 (2013).

##Downloading the code

To clone the repository, you can simply run
```
git clone --recursive git://github.com/rborrelli/molfc.git

```
To update to the latest version, run
```
git pull
```
##Prerequisites

**It is highly recommended** that you use the Intel Fortran Compiler releases 15.0 or laters.
**Gfortran compatibility is not guaranteed at the moment. The code won't compile with gfortran version <= 4.9. Minimum gfortran version is 5.0.

##Installing the package
The installation of the package is done via cmake.
Follow the instructions in the INSTALL file.

