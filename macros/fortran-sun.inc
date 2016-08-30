#
# Fortran macros for Sun Fortran compiler
#
.SUFFIXES: .SUFFIXES .f .f90 .f95 .F .F90 .F95 .fpp
#
OBJ_DIR=$(FCQD_ROOT)/obj/
LIB_STD=$(FCQD_ROOT)/lib/
MOD_STD=$(FCQD_ROOT)/obj/
INC_STD=$(FCQD_ROOT)/src/include/
BIN_STD=$(FCQD_ROOT)/bin/
#
######### Fortran Configuration ##############
# Laptop
SUNDIR=/opt
FC=$(SUNDIR)/sun/sunstudioceres/bin/sunf95
###########
#University
#SUNDIR=/usr/local
#FC=$(SUNDIR)/sun/sunstudio12/bin/sunf95
#########
FFLAGS_DEBUG= 
# type make FFLAGS_DEBUG= "-g" for debugging.
FFLAGS_DEBUG= -g
#FFLAGS= -cm -w95 -O3 -xW -tpp7 -static $(FFLAGS_DEBUG)
FFLAGS2= -cm -w95 -O3 -xW -tpp7 -module $(MOD_STD)
FFLAGS1= -cm -w95 -O3 -xW -tpp7 -module $(MOD_STD) $(FFLAGS_DEBUG)
#FFLAGS0= -cm -O0 -module $(MOD_STD) $(FFLAGS_DEBUG)
FFLAGS0= -O0 -moddir=$(MOD_STD) -M$(MOD_STD) -g
FFLAGS= $(FFLAGS0)
#LDFLAGS=  -static $(CFLAGS)
#LDFLAGS=  $(CFLAGS)
#LDFLAGS= -cm -O0 
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
RANLIB=echo
#
DEFS= -D_CPU_TIME_  -DMAXMEM=$(MAXMEM) -DSUNF95
#
CP=cp
MV=mv
#
# Experimental : the following deactivates an implicit rule
# which breaks havoc with the operation of this makefile
# It works at least with GNU make
%.o : %.mod
#
.F.o:
	$(FC) -c $(MOD_SEARCH) $(INC_SEARCH) $(FFLAGS)  $(OPT) $(DEFS) $<
$(OBJ_DIR)%.o: %.f
	$(FC) -c $(MOD_SEARCH) $(INC_SEARCH) $(FFLAGS)  $(OPT) $< -o $@  $(DEVNULL) 
#.f.o:
#	$(FC) -c $(MOD_SEARCH) $(INC_SEARCH) $(FFLAGS)  $(OPT) $< -o $@
#.F90.o:
$(OBJ_DIR)%.o: %.F90
	$(FC) -c $(MOD_SEARCH) $(INC_SEARCH) $(FFLAGS) $(OPT) $(DEFS) $< -o $@ $(DEVNULL)
#	$(FC) -c $(MOD_SEARCH) $(INC_SEARCH) $(FFLAGS) $(OPT) $(DEFS) $< -o $(OBJ_DIR)$@
#.f90.o:
#	$(FC) -c $(MOD_SEARCH) $(INC_SEARCH) $(FFLAGS) $(OPT)  $< 
##	$(FC) -c $(MOD_SEARCH) $(INC_SEARCH) $(FFLAGS) $(OPT)  $< -o $(OBJ_DIR)$@
$(OBJ_DIR)%.o: %.f90
	$(FC) -c $(MOD_SEARCH) $(INC_SEARCH) $(FFLAGS) $(OPT)  $< -o $@ $(DEVNULL)   
#	$(FC) -c $(MOD_SEARCH) $(INC_SEARCH) $(FFLAGS) $(OPT)  $<
#	$(FC) -c $(MOD_SEARCH) $(INC_SEARCH) $(FFLAGS) $(OPT)  $< -o $(OBJ_DIR)$@
.F95.o:
	$(FC) -c $(MOD_SEARCH) $(INC_SEARCH) $(FFLAGS) $(OPT) $(DEFS) $<
.f95.o:
	$(FC) -c $(MOD_SEARCH) $(INC_SEARCH) $(FFLAGS) $(OPT)  $<
#
.fpp.f90:
	$(FC) -c $(MOD_SEARCH) $(INC_SEARCH) $(FFLAGS) $(OPT) $(DEFS) $<
#
