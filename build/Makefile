#
# makefile
#

BINROOT = ../execute
FC=undefined
FFLAGS=undefined
LDFLAGS=undefined
TYPE=undefined
EX=main.x
EXE = $(BINROOT)/$(EX)
SHELL=/bin/sh


OBJ_MOD = serial_version_mod.o parallel_version_mod.o parse_mod.o error_mod.o \
	utility_mod.o nanojet_mod.o fit_mod.o support_functions_mod.o \
	dynamic_refinement_mod.o viscoelastic_force_mod.o coulomb_force_mod.o \
	eom_mod.o driver_eom_mod.o integrator_mod.o statistic_mod.o io_mod.o

OBJ_SRC = main.o



all:
	@echo "Error - you should specify a target machine!"
	@echo "Possible choices are:              "
	@echo "                                   "
	@echo "gfortran                           "
	@echo "gfortran-mpi                       "
	@echo "cygwin                             "
	@echo "cygwin-mpi                         "
	@echo "intel                              "
	@echo "intel-mpi                          "
	@echo "intel-openmpi                      "
	@echo "                                   "
	@echo "Possible choices for debugging are:"
	@echo "                                   "
	@echo "gfortran-debugger                  "
	@echo "gfortran-mpidebugger               "
	@echo "intel-debugger                     "
	@echo "intel-mpidebugger                  "
	@echo "                                   "
	@echo "Please examine Makefile for further details "

help:
	@echo "help: List of possible target machine"
	@echo "Possible choices are:              "
	@echo "                                   "
	@echo "gfortran                           "
	@echo "gfortran-mpi                       "
	@echo "cygwin                             "
	@echo "cygwin-mpi                         "
	@echo "intel                              "
	@echo "intel-mpi                          "
	@echo "intel-openmpi                      "
	@echo "                                   "
	@echo "Possible choices for debugging are:"
	@echo "                                   "
	@echo "gfortran-debugger                  "
	@echo "gfortran-mpidebugger               "
	@echo "intel-debugger                     "
	@echo "intel-mpidebugger                  "
	@echo "                                   "
	@echo "Please examine Makefile for further details "

cygwin:
	$(MAKE) FC=gfortran \
	FFLAGS="-O2 -c" \
	LDFLAGS="-O2 -o" \
	TYPE="seq" \
	EX=$(EX) BINROOT=$(BINROOT) seq
	
cygwin-mpi:
	$(MAKE) FC=mpif90 \
	FFLAGS="-O2 -c" \
	LDFLAGS="-O2 -o" \
	TYPE=mpi \
	EX=$(EX) BINROOT=$(BINROOT) mpi
	
gfortran:
	$(MAKE) FC=gfortran \
	FFLAGS="-O2 -c" \
	LDFLAGS="-O2 -o" \
	TYPE="seq" \
	EX=$(EX) BINROOT=$(BINROOT) seq
	
gfortran-mpi:
	$(MAKE) FC=mpif90 \
	FFLAGS="-O2 -c" \
	LDFLAGS="-O2 -o" \
	TYPE=mpi \
	EX=$(EX) BINROOT=$(BINROOT) mpi

gfortran-mpidebugger:
	$(MAKE) FC=mpif90 \
	FFLAGS="-O0 -g -fcheck=all -c" \
	LDFLAGS="-O0 -g -fcheck=all -o" \
	TYPE=mpi \
	EX=$(EX) BINROOT=$(BINROOT) mpi

gfortran-debugger:
	$(MAKE) FC=gfortran \
	FFLAGS="-O0 -g -fcheck=all -c" \
	LDFLAGS="-O0 -g -fcheck=all -o" \
	TYPE=seq \
	EX=$(EX) BINROOT=$(BINROOT) seq

intel:
	$(MAKE) FC=ifort \
	FFLAGS="-O2 -c" \
	LDFLAGS="-O2 -o" \
	TYPE=seq \
	EX=$(EX) BINROOT=$(BINROOT) seq
	
intel-mpi:
	$(MAKE) FC=mpiifort \
	FFLAGS="-O2 -c" \
	LDFLAGS="-O2 -o" \
	TYPE=mpi \
	EX=$(EX) BINROOT=$(BINROOT) mpi

intel-openmpi:
	$(MAKE) FC=mpif90 \
	FFLAGS="-O2 -c" \
	LDFLAGS="-O2  -o" \
	TYPE=mpi \
	EX=$(EX) BINROOT=$(BINROOT) mpi

intel-debugger:
	$(MAKE) FC=ifort \
	FFLAGS="-O0 -g -check all -fpe0 -traceback -debug extended -c" \
	LDFLAGS="-O0 -g -check all -fpe0 -traceback -debug extended -o" \
	TYPE=seq \
	EX=$(EX) BINROOT=$(BINROOT) seq

intel-mpidebugger:
	$(MAKE) FC=mpiifort \
	FFLAGS="-O0 -g -check all -fpe0 -traceback -debug extended -c" \
	LDFLAGS="-O0 -g -check all -fpe0 -traceback -debug extended -o" \
	TYPE=mpi \
	EX=$(EX) BINROOT=$(BINROOT) mpi

seq:serial_version_mod.o parse_mod.o error_mod.o utility_mod.o nanojet_mod.o \
	fit_mod.o support_functions_mod.o dynamic_refinement_mod.o \
	viscoelastic_force_mod.o coulomb_force_mod.o eom_mod.o driver_eom_mod.o \
	integrator_mod.o statistic_mod.o io_mod.o main.o
	$(FC) $(LDFLAGS) $(EX) serial_version_mod.o parse_mod.o error_mod.o utility_mod.o \
	nanojet_mod.o fit_mod.o support_functions_mod.o dynamic_refinement_mod.o \
	viscoelastic_force_mod.o coulomb_force_mod.o eom_mod.o driver_eom_mod.o \
	integrator_mod.o statistic_mod.o io_mod.o main.o
	mv $(EX) $(EXE)
	
mpi:parallel_version_mod.o parse_mod.o error_mod.o utility_mod.o nanojet_mod.o \
	fit_mod.o support_functions_mod.o dynamic_refinement_mod.o \
	viscoelastic_force_mod.o coulomb_force_mod.o eom_mod.o driver_eom_mod.o \
	integrator_mod.o statistic_mod.o io_mod.o main.o
	$(FC) $(LDFLAGS) $(EX) parallel_version_mod.o parse_mod.o error_mod.o utility_mod.o \
	nanojet_mod.o fit_mod.o support_functions_mod.o dynamic_refinement_mod.o \
	viscoelastic_force_mod.o coulomb_force_mod.o eom_mod.o driver_eom_mod.o \
	integrator_mod.o statistic_mod.o io_mod.o main.o
	mv $(EX) $(EXE)

serial_version_mod.o:serial_version_mod.f90
	$(FC) $(FFLAGS) serial_version_mod.f90
	
parallel_version_mod.o:parallel_version_mod.f90
	$(FC) $(FFLAGS) parallel_version_mod.f90
	
parse_mod.o:parse_mod.f90
	$(FC) $(FFLAGS) parse_mod.f90

error_mod.o:error_mod.f90
	$(FC) $(FFLAGS) error_mod.f90

utility_mod.o:utility_mod.f90
	$(FC) $(FFLAGS) utility_mod.f90
	
nanojet_mod.o:nanojet_mod.f90
	$(FC) $(FFLAGS) nanojet_mod.f90
	
fit_mod.o:fit_mod.f90
	$(FC) $(FFLAGS) fit_mod.f90
	
support_functions_mod.o:support_functions_mod.f90
	$(FC) $(FFLAGS) support_functions_mod.f90
	
dynamic_refinement_mod.o:dynamic_refinement_mod.f90
	$(FC) $(FFLAGS) dynamic_refinement_mod.f90
	
viscoelastic_force_mod.o: viscoelastic_force_mod.f90
	$(FC) $(FFLAGS) viscoelastic_force_mod.f90
	
coulomb_force_mod.o: coulomb_force_mod.f90
	$(FC) $(FFLAGS) coulomb_force_mod.f90
	
eom_mod.o: eom_mod.f90
	$(FC) $(FFLAGS) eom_mod.f90
	
driver_eom_mod.o: driver_eom_mod.f90
	$(FC) $(FFLAGS) driver_eom_mod.f90

integrator_mod.o:integrator_mod.f90
	$(FC) $(FFLAGS) integrator_mod.f90

statistic_mod.o:statistic_mod.f90
	$(FC) $(FFLAGS) statistic_mod.f90

io_mod.o:io_mod.f90
	$(FC) $(FFLAGS) io_mod.f90
	
main.o:main.f90
	$(FC) $(FFLAGS) main.f90

clean:
	rm -rf $(OBJ_MOD) $(OBJ_SRC) *.mod






