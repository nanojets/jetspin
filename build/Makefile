#
# makefile
#

BINROOT = ../execute
FC=undefined
FFLAGS=undefined
FPPFLAGS ?= -cpp
FPPFLAGS_EXTRA ?=
LDFLAGS=undefined
TYPE=undefined
EX=main.x
EXE = $(BINROOT)/$(EX)
SHELL=/bin/sh
MPI_COMPAT_FLAG ?= -fallow-argument-mismatch
MPIFC ?= mpif90
NVFORTRAN ?= nvfortran
GPUCC ?= 80
CUDA_VERSION ?= 12.3
OPENACC_GPU_FLAGS = -gpu=cc$(GPUCC),cuda$(CUDA_VERSION),nofma
ACCOBJ ?= host_accelerator_mod.o
MAKEFILE_PATH := $(lastword $(MAKEFILE_LIST))


OBJ_MOD = serial_version_mod.o parallel_version_mod.o host_accelerator_mod.o \
	openacc_accelerator_mod.o parse_mod.o error_mod.o \
	profiling_mod.o utility_mod.o nanojet_mod.o fit_mod.o support_functions_mod.o \
	breaking_mod.o dynamic_refinement_mod.o viscoelastic_force_mod.o \
	coulomb_force_mod.o electric_field_mod.o eom_mod.o eom_ev_mod.o driver_eom_mod.o \
	integrator_mod.o integrator_kv_ev_mod.o statistic_mod.o io_mod.o

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
	@echo "nvfortran                          "
	@echo "nvfortran-mpi                      "
	@echo "nvfortran-openacc                  "
	@echo "nvfortran-openacc-host             "
	@echo "nvfortran-openacc-force-oracle     "
	@echo "nvfortran-openacc-coulomb-oracle   "
	@echo "nvfortran-openacc-host-akima       "
	@echo "nvfortran-openacc-compare-akima    "
	@echo "nvfortran-openacc-compare-refinement"
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
	@echo "nvfortran                          "
	@echo "nvfortran-mpi                      "
	@echo "nvfortran-openacc                  "
	@echo "nvfortran-openacc-host             "
	@echo "nvfortran-openacc-force-oracle     "
	@echo "nvfortran-openacc-coulomb-oracle   "
	@echo "nvfortran-openacc-host-akima       "
	@echo "nvfortran-openacc-compare-akima    "
	@echo "nvfortran-openacc-compare-refinement"
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
	$(MAKE) -f $(MAKEFILE_PATH) FC=gfortran \
	FFLAGS="-O2 -c" \
	LDFLAGS="-O2 -o" \
	TYPE="seq" \
	EX=$(EX) BINROOT=$(BINROOT) seq
	
cygwin-mpi:
	$(MAKE) -f $(MAKEFILE_PATH) FC=mpif90 \
	FFLAGS="-O2 -c" \
	LDFLAGS="-O2 -o" \
	TYPE=mpi \
	EX=$(EX) BINROOT=$(BINROOT) mpi
	
gfortran:
	$(MAKE) -f $(MAKEFILE_PATH) FC=gfortran \
	FFLAGS="-O2 -c" \
	LDFLAGS="-O2 -o" \
	TYPE="seq" \
	EX=$(EX) BINROOT=$(BINROOT) seq
	
gfortran-mpi:
	$(MAKE) -f $(MAKEFILE_PATH) FC=$(MPIFC) \
	FFLAGS="-O2 $(MPI_COMPAT_FLAG) -c" \
	LDFLAGS="-O2 -o" \
	TYPE=mpi \
	EX=$(EX) BINROOT=$(BINROOT) mpi

nvfortran:
	$(MAKE) -f $(MAKEFILE_PATH) FC=$(NVFORTRAN) \
	FFLAGS="-O3 -c" \
	FPPFLAGS="-Mpreprocess $(FPPFLAGS_EXTRA)" \
	LDFLAGS="-O3 -o" \
	ACCOBJ=host_accelerator_mod.o TYPE=seq \
	EX=$(EX) BINROOT=$(BINROOT) seq

nvfortran-mpi:
	$(MAKE) -f $(MAKEFILE_PATH) FC=$(MPIFC) \
	FFLAGS="-O3 -c" \
	FPPFLAGS="-Mpreprocess $(FPPFLAGS_EXTRA)" \
	LDFLAGS="-O3 -o" \
	ACCOBJ=host_accelerator_mod.o TYPE=mpi \
	EX=$(EX) BINROOT=$(BINROOT) mpi

nvfortran-openacc:
	$(MAKE) -f $(MAKEFILE_PATH) FC=$(NVFORTRAN) \
	FFLAGS="-O3 -acc=gpu $(OPENACC_GPU_FLAGS) -Minfo=accel -c" \
	FPPFLAGS="-Mpreprocess $(FPPFLAGS_EXTRA)" \
	LDFLAGS="-O3 -acc=gpu $(OPENACC_GPU_FLAGS) -o" \
	ACCOBJ=openacc_accelerator_mod.o TYPE=seq \
	EX=$(EX) BINROOT=$(BINROOT) seq

# Development-only OpenACC oracle: evaluate every supported force stage
# through the trusted host equations, with or without evaporation, upload the
# derivatives, and leave integration updates and topology on the device.
nvfortran-openacc-force-oracle:
	$(MAKE) -f $(MAKEFILE_PATH) FC=$(NVFORTRAN) \
	FFLAGS="-O3 -acc=gpu $(OPENACC_GPU_FLAGS) -Minfo=accel -c" \
	FPPFLAGS="-Mpreprocess -DJETSPIN_DEV_HOST_FORCE_ORACLE $(FPPFLAGS_EXTRA)" \
	LDFLAGS="-O3 -acc=gpu $(OPENACC_GPU_FLAGS) -o" \
	ACCOBJ=openacc_accelerator_mod.o TYPE=seq \
	EX=$(EX) BINROOT=$(BINROOT) seq

# Narrower development diagnostic: compute only the direct Coulomb sum on the
# host, with or without evaporation, then upload the resulting force array.
# All other force terms, integration updates, and topology remain on device.
nvfortran-openacc-coulomb-oracle:
	$(MAKE) -f $(MAKEFILE_PATH) FC=$(NVFORTRAN) \
	FFLAGS="-O3 -acc=gpu $(OPENACC_GPU_FLAGS) -Minfo=accel -c" \
	FPPFLAGS="-Mpreprocess -DJETSPIN_DEV_HOST_COULOMB_ORACLE $(FPPFLAGS_EXTRA)" \
	LDFLAGS="-O3 -acc=gpu $(OPENACC_GPU_FLAGS) -o" \
	ACCOBJ=openacc_accelerator_mod.o TYPE=seq \
	EX=$(EX) BINROOT=$(BINROOT) seq

# Development-only Akima oracle: retain the historical host coefficient and
# interpolation path while the surrounding refinement lifecycle stays active.
nvfortran-openacc-host-akima:
	$(MAKE) -f $(MAKEFILE_PATH) FC=$(NVFORTRAN) \
	FFLAGS="-O3 -acc=gpu $(OPENACC_GPU_FLAGS) -Minfo=accel -c" \
	FPPFLAGS="-Mpreprocess -DJETSPIN_DEV_HOST_AKIMA $(FPPFLAGS_EXTRA)" \
	LDFLAGS="-O3 -acc=gpu $(OPENACC_GPU_FLAGS) -o" \
	ACCOBJ=openacc_accelerator_mod.o TYPE=seq \
	EX=$(EX) BINROOT=$(BINROOT) seq

# Development-only A/B build: compute host reference coefficients and values,
# then run the device Akima path and report field-wise maximum differences.
nvfortran-openacc-compare-akima:
	$(MAKE) -f $(MAKEFILE_PATH) FC=$(NVFORTRAN) \
	FFLAGS="-O3 -acc=gpu $(OPENACC_GPU_FLAGS) -Minfo=accel -c" \
	FPPFLAGS="-Mpreprocess -DJETSPIN_COMPARE_AKIMA $(FPPFLAGS_EXTRA)" \
	LDFLAGS="-O3 -acc=gpu $(OPENACC_GPU_FLAGS) -o" \
	ACCOBJ=openacc_accelerator_mod.o TYPE=seq \
	EX=$(EX) BINROOT=$(BINROOT) seq

# Development-only A/B build: compute host reference volume/evaporation-
# volume reconstruction, conservation rescale, and density conversion, then
# run the device path and report field-wise maximum differences.
nvfortran-openacc-compare-refinement:
	$(MAKE) -f $(MAKEFILE_PATH) FC=$(NVFORTRAN) \
	FFLAGS="-O3 -acc=gpu $(OPENACC_GPU_FLAGS) -Minfo=accel -c" \
	FPPFLAGS="-Mpreprocess -DJETSPIN_COMPARE_REFINEMENT_ASSEMBLY $(FPPFLAGS_EXTRA)" \
	LDFLAGS="-O3 -acc=gpu $(OPENACC_GPU_FLAGS) -o" \
	ACCOBJ=openacc_accelerator_mod.o TYPE=seq \
	EX=$(EX) BINROOT=$(BINROOT) seq

nvfortran-openacc-host:
	$(MAKE) -f $(MAKEFILE_PATH) FC=$(NVFORTRAN) \
	FFLAGS="-O2 -acc=host -Minfo=accel -c" \
	FPPFLAGS="-Mpreprocess $(FPPFLAGS_EXTRA)" \
	LDFLAGS="-O2 -acc=host -o" \
	ACCOBJ=openacc_accelerator_mod.o TYPE=seq \
	EX=$(EX) BINROOT=$(BINROOT) seq

gfortran-mpidebugger:
	$(MAKE) -f $(MAKEFILE_PATH) FC=$(MPIFC) \
	FFLAGS="-O0 -g -fcheck=all $(MPI_COMPAT_FLAG) -c" \
	LDFLAGS="-O0 -g -fcheck=all -o" \
	TYPE=mpi \
	EX=$(EX) BINROOT=$(BINROOT) mpi

gfortran-debugger:
	$(MAKE) -f $(MAKEFILE_PATH) FC=gfortran \
	FFLAGS="-O0 -g -fcheck=all -c" \
	LDFLAGS="-O0 -g -fcheck=all -o" \
	TYPE=seq \
	EX=$(EX) BINROOT=$(BINROOT) seq

intel:
	$(MAKE) -f $(MAKEFILE_PATH) FC=ifort \
	FFLAGS="-O2 -c" \
	FPPFLAGS="-fpp" \
	LDFLAGS="-O2 -o" \
	TYPE=seq \
	EX=$(EX) BINROOT=$(BINROOT) seq
	
intel-mpi:
	$(MAKE) -f $(MAKEFILE_PATH) FC=mpiifort \
	FFLAGS="-O2 -c" \
	FPPFLAGS="-fpp" \
	LDFLAGS="-O2 -o" \
	TYPE=mpi \
	EX=$(EX) BINROOT=$(BINROOT) mpi

intel-openmpi:
	$(MAKE) -f $(MAKEFILE_PATH) FC=mpif90 \
	FFLAGS="-O2 -c" \
	FPPFLAGS="-fpp" \
	LDFLAGS="-O2  -o" \
	TYPE=mpi \
	EX=$(EX) BINROOT=$(BINROOT) mpi

intel-debugger:
	$(MAKE) -f $(MAKEFILE_PATH) FC=ifort \
	FFLAGS="-O0 -g -check all -fpe0 -traceback -debug extended -c" \
	FPPFLAGS="-fpp" \
	LDFLAGS="-O0 -g -check all -fpe0 -traceback -debug extended -o" \
	TYPE=seq \
	EX=$(EX) BINROOT=$(BINROOT) seq

intel-mpidebugger:
	$(MAKE) -f $(MAKEFILE_PATH) FC=mpiifort \
	FFLAGS="-O0 -g -check all -fpe0 -traceback -debug extended -c" \
	FPPFLAGS="-fpp" \
	LDFLAGS="-O0 -g -check all -fpe0 -traceback -debug extended -o" \
	TYPE=mpi \
	EX=$(EX) BINROOT=$(BINROOT) mpi

seq:serial_version_mod.o $(ACCOBJ) profiling_mod.o parse_mod.o error_mod.o utility_mod.o \
	nanojet_mod.o fit_mod.o support_functions_mod.o breaking_mod.o \
	dynamic_refinement_mod.o viscoelastic_force_mod.o coulomb_force_mod.o \
	electric_field_mod.o eom_mod.o eom_ev_mod.o driver_eom_mod.o integrator_mod.o \
	integrator_kv_ev_mod.o statistic_mod.o io_mod.o main.o
	$(FC) $(LDFLAGS) $(EX) serial_version_mod.o $(ACCOBJ) profiling_mod.o parse_mod.o error_mod.o \
	utility_mod.o nanojet_mod.o fit_mod.o support_functions_mod.o \
	breaking_mod.o dynamic_refinement_mod.o viscoelastic_force_mod.o \
	coulomb_force_mod.o electric_field_mod.o eom_mod.o eom_ev_mod.o driver_eom_mod.o \
	integrator_mod.o integrator_kv_ev_mod.o statistic_mod.o io_mod.o main.o
	mv $(EX) $(EXE)
	
mpi:parallel_version_mod.o $(ACCOBJ) profiling_mod.o parse_mod.o error_mod.o utility_mod.o \
	nanojet_mod.o fit_mod.o support_functions_mod.o breaking_mod.o \
	dynamic_refinement_mod.o viscoelastic_force_mod.o coulomb_force_mod.o \
	electric_field_mod.o eom_mod.o eom_ev_mod.o driver_eom_mod.o integrator_mod.o \
	integrator_kv_ev_mod.o statistic_mod.o io_mod.o main.o
	$(FC) $(LDFLAGS) $(EX) parallel_version_mod.o $(ACCOBJ) profiling_mod.o parse_mod.o error_mod.o \
	utility_mod.o nanojet_mod.o fit_mod.o support_functions_mod.o \
	breaking_mod.o dynamic_refinement_mod.o viscoelastic_force_mod.o \
	coulomb_force_mod.o electric_field_mod.o eom_mod.o eom_ev_mod.o driver_eom_mod.o \
	integrator_mod.o integrator_kv_ev_mod.o statistic_mod.o io_mod.o main.o
	mv $(EX) $(EXE)

serial_version_mod.o:serial_version_mod.f90
	$(FC) $(FFLAGS) serial_version_mod.f90
	
parallel_version_mod.o:parallel_version_mod.f90
	$(FC) $(FFLAGS) parallel_version_mod.f90

host_accelerator_mod.o:host_accelerator_mod.f90
	$(FC) $(FFLAGS) host_accelerator_mod.f90

openacc_accelerator_mod.o:openacc_accelerator_mod.f90
	$(FC) $(FFLAGS) $(FPPFLAGS) openacc_accelerator_mod.f90

profiling_mod.o:profiling_mod.f90
	$(FC) $(FFLAGS) profiling_mod.f90
	
parse_mod.o:parse_mod.f90
	$(FC) $(FFLAGS) parse_mod.f90

error_mod.o:error_mod.f90
	$(FC) $(FFLAGS) error_mod.f90

utility_mod.o:utility_mod.f90
	$(FC) $(FFLAGS) $(FPPFLAGS) utility_mod.f90
	
nanojet_mod.o:nanojet_mod.f90
	$(FC) $(FFLAGS) nanojet_mod.f90
	
fit_mod.o:fit_mod.f90
	$(FC) $(FFLAGS) $(FPPFLAGS) fit_mod.f90
	
support_functions_mod.o:support_functions_mod.f90
	$(FC) $(FFLAGS) support_functions_mod.f90
	
breaking_mod.o:breaking_mod.f90
	$(FC) $(FFLAGS) breaking_mod.f90
	
dynamic_refinement_mod.o:dynamic_refinement_mod.f90
	$(FC) $(FFLAGS) $(FPPFLAGS) dynamic_refinement_mod.f90
	
viscoelastic_force_mod.o: viscoelastic_force_mod.f90
	$(FC) $(FFLAGS) viscoelastic_force_mod.f90
	
coulomb_force_mod.o: coulomb_force_mod.f90
	$(FC) $(FFLAGS) $(FPPFLAGS) coulomb_force_mod.f90
	
electric_field_mod.o: electric_field_mod.f90
	$(FC) $(FFLAGS) electric_field_mod.f90
	
eom_mod.o: eom_mod.f90
	$(FC) $(FFLAGS) eom_mod.f90
	
eom_ev_mod.o: eom_ev_mod.f90
	$(FC) $(FFLAGS) $(FPPFLAGS) eom_ev_mod.f90
	
driver_eom_mod.o: driver_eom_mod.f90
	$(FC) $(FFLAGS) driver_eom_mod.f90

integrator_mod.o:integrator_mod.f90 statistic_mod.o
	$(FC) $(FFLAGS) $(FPPFLAGS) integrator_mod.f90

integrator_kv_ev_mod.o:integrator_kv_ev_mod.f90 statistic_mod.o
	$(FC) $(FFLAGS) $(FPPFLAGS) integrator_kv_ev_mod.f90

statistic_mod.o:statistic_mod.f90
	$(FC) $(FFLAGS) $(FPPFLAGS) statistic_mod.f90

io_mod.o:io_mod.f90
	$(FC) $(FFLAGS) io_mod.f90
	
main.o:main.f90
	$(FC) $(FFLAGS) $(FPPFLAGS) main.f90

clean:
	rm -rf $(OBJ_MOD) $(OBJ_SRC) *.mod
