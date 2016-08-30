#
# Fortran macros for PGI f90
#
.SUFFIXES: .SUFFIXES .f .f90 .f95 .F .F90 .F95
#
LIB_STD=$(FCQD_ROOT)/lib/
MOD_STD=$(FCQD_ROOT)/modules/
INC_STD=$(FCQD_ROOT)/include/
BIN_STD=$(FCQD_ROOT)/bin/
#
FC=pgf90
CFLAGS= 
FFLAGS= -O   $(CFLAGS)
FFLAGS_DEBUG= -g 
LDFLAGS=      $(CFLAGS)
OPT=
#
INC_PREFIX=-I
MOD_PREFIX=-I
LIB_PREFIX=-L
#
MOD_EXT=.mod
MOD_SEARCH_STD= $(MOD_PREFIX)$(MOD_STD) 
MOD_SEARCH= $(MOD_SEARCH_STD) $(MOD_SEARCH_OTHER)
INC_SEARCH= $(INC_PREFIX)$(INC_STD)
#
#
AR=ar
RANLIB=ranlib
#
CPP=/bin/cpp -P
COCO=$(BIN_STD)coco -I $(INC_STD)
DEFS= -DMAXMEM=$(MAXMEM)
#
# Experimental : the following deactivates an implicit rule
# which breaks havoc with the operation of this makefile
# It works at least with GNU make
%.o : %.mod
#
.F.o:
	$(FC) -c $(MOD_SEARCH) $(INC_SEARCH) $(FFLAGS)  $(DEFS) $<
.f.o:
	$(FC) -c $(MOD_SEARCH) $(INC_SEARCH) $(FFLAGS)   $<
.F90.o:
	$(CPP) $(DEFS) $< > $*.f90
	$(FC) -c $(MOD_SEARCH) $(INC_SEARCH) $(FFLAGS) $*.f90
	@rm -f $*.f90
.f90.o:
	$(FC) -c $(MOD_SEARCH) $(INC_SEARCH) $(FFLAGS)   $<
#

