#
# Fortran macros for Intel Fortran compiler
#
.SUFFIXES: .SUFFIXES .f .f90 .f95 .F .F90 .F95 .fpp
#
OBJ_DIR=$(FCQD_ROOT)/obj/
LIB_STD=$(FCQD_ROOT)/lib/
MOD_STD=$(FCQD_ROOT)/obj/
INC_STD=$(FCQD_ROOT)/src/include/
BIN_STD=$(FCQD_ROOT)/bin/
#Computer personale
#FC=/home/shared/intel/fc/9.1.036/bin/ifort
#Cluster
#FC=/home/shared/intel/fc/9.1.037/bin/ifort
#Portatile
#FC=/opt/intel/fc/9.1.037/bin/ifort
FC=/opt/intel/fc/10.0.023/bin/ifort
FFLAGS_DEBUG= 
# type make FFLAGS_DEBUG= "-g" for debugging with IDB
#FFLAGS_DEBUG= -g -debug all
FFLAGS_DEBUG= -g -debug all -traceback -CB
#FFLAGS_DEBUG= -g
###################################
#processor type
# T = core 2 duo
# W = pentium 4
PROCESSOR=T
###################################
# Flags for Optimization with Intel Fortran 9.1
#FFLAGS2_9_1= -cm -w95 -O3 -x$(PROCESSOR) -tpp7 -module $(MOD_STD)
#FFLAGS1_9_1= -cm -w95 -O3 -x$(PROCESSOR) -tpp7 -module $(MOD_STD) $(FFLAGS_DEBUG)
###################################
# Flags for Optimization with Intel Fortran 10.0
FFLAGS3_10PAR= -cm -w95 -O3 -x$(PROCESSOR) -openmp -module $(MOD_STD)
FFLAGS3_10= -cm -w95 -O3 -x$(PROCESSOR) -module $(MOD_STD)
FFLAGS2_10= -cm -w95 -O2 -x$(PROCESSOR) -module $(MOD_STD)
FFLAGS1_10= -cm -w95 -x$(PROCESSOR) -module $(MOD_STD) $(FFLAGS_DEBUG)
##################################
# Flags for debugging with Intel Fortran
FFLAGS0= -cm -O0 -module $(MOD_STD) $(FFLAGS_DEBUG)
FFLAGS0_10= $(FFLAGS0)
##################################
# Actually used compiler options
FFLAGS= $(FFLAGS0_10)
#LDFLAGS=  -static $(CFLAGS)
#LDFLAGS=  $(CFLAGS)
#LDFLAGS= -cm -O0 
OPT=
#
INC_PREFIX=-I
MOD_PREFIX=-I
LIB_PREFIX=-L
#
MOD_EXT=mod
MOD_SEARCH_STD= $(MOD_PREFIX)$(MOD_STD) 
MOD_SEARCH= $(MOD_SEARCH_STD) $(MOD_SEARCH_OTHER)
INC_SEARCH= $(INC_PREFIX)$(INC_STD)
#
#
AR=ar
RANLIB=echo
#
DEFS= -D_CPU_TIME_  -DMAXMEM=$(MAXMEM) -DINTELF95
#
CP=cp
MV=mv
# redirec to null device
#DEVNULL= > /dev/null 2>&1
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



