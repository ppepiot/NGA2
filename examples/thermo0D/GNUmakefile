# 0D constant-volume NASG heat-addition reactor
# Run: mpirun -np 1 ./thermo0D.dp.gnu.opt.mpi.exe -i input

NGA_HOME ?= ../..
PRECISION = DOUBLE
USE_MPI   = TRUE
USE_CVODE = FALSE
COMP      = gnu
DEBUG     = FALSE
EXEBASE   = thermo0D

Udirs   := src
Upack   += $(foreach dir, $(Udirs), $(wildcard $(dir)/Make.package))
Ulocs   += $(foreach dir, $(Udirs), $(wildcard $(dir)))
include $(Upack)
INCLUDE_LOCATIONS += $(Ulocs)
VPATH_LOCATIONS   += $(Ulocs)

include $(NGA_HOME)/tools/GNUMake/Make.defs

# Libraries + thermo (NASG hierarchy). Filter to the minimal set we need.
Bdirs   := thermo libraries
Bpack   += $(foreach dir, $(Bdirs), $(NGA_HOME)/src/$(dir)/Make.package)
include $(Bpack)

f90EXE_sources := $(filter \
	precision_$(PREC).f90 string.f90 parallel.f90 messager.f90 param.f90 material_class.f90 thermorelax_class.f90 \
	ideal_gas_class.f90 stiffened_gas_class.f90 nasg_class.f90 relax_ig_nasg_class.f90, \
	$(f90EXE_sources))
f90EXE_sources += main.f90

all: $(executable)
	@echo "Build complete. Run: mpirun -np 1 ./$(executable) -i input"

include $(NGA_HOME)/tools/GNUMake/Make.rules
