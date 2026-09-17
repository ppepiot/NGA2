# Two-dimensional laminar slot burner (Selle, Poinsot & Ferret, Combust. Flame 158 (2011) 146-154, Figs. 6 and 11).
# The kinetics module src/chem_data_fc.f90 is generated from a Cantera YAML mechanism; each input file
# names the one it expects in its header. From this directory:
#   python3 ../../tools/scripts/chemistry/yaml2nga.py <mechanism>.yaml src/chem_data_fc.f90
# Build: make AMREX_DIR=... HYPRE_DIR=... HDF5_DIR=...   Run: mpiexec -n 4 ./nga2.dp.gnu.opt.mpi.exe -i input_ch4_phi1.0

# NGA location
NGA_HOME ?= ../..

# Compilation parameters
PRECISION = DOUBLE
USE_MPI   = TRUE
USE_FFTW  = FALSE
USE_HYPRE = TRUE
USE_LAPACK= FALSE
USE_AMREX = TRUE
USE_HDF5  = TRUE
PROFILE   = FALSE
DEBUG     = FALSE
COMP      = gnu
EXEBASE   = nga2

# Directories that contain user-defined code
Udirs   := src

# Include user-defined sources
Upack   += $(foreach dir, $(Udirs), $(wildcard $(dir)/Make.package))
Ulocs   += $(foreach dir, $(Udirs), $(wildcard $(dir)))
include $(Upack)
INCLUDE_LOCATIONS += $(Ulocs)
VPATH_LOCATIONS   += $(Ulocs)

# Kinetics module. src/chem_data_fc.f90 is generated for one specific mechanism and is NOT kept in the
# repository (see .gitignore) -- produce it with, from this directory:
#   python3 ../../tools/scripts/chemistry/yaml2nga.py <mechanism>.yaml src/chem_data_fc.f90
# The line below is searched after src/, so a generated module always wins. With none present the build
# falls back to the minimal placeholder module in tools/scripts/chemistry/skeleton: it compiles and links,
# but it carries only inert air species, so a real case stops at start-up instead of running silently
# against the wrong chemistry.
VPATH_LOCATIONS   += $(NGA_HOME)/tools/scripts/chemistry/skeleton

# NGA compilation definitions
include $(NGA_HOME)/tools/GNUMake/Make.defs

# Include NGA base code
Bdirs   := amrbase amrsolvers thermo chemistry libraries core
Bpack   += $(foreach dir, $(Bdirs), $(NGA_HOME)/src/$(dir)/Make.package)
include $(Bpack)

# Inform user of Make.packages used
ifdef Ulocs
   $(info Taking user code from: $(Ulocs))
endif
$(info Taking base code from: $(Bdirs))

# Target definition
all: $(executable)
	@echo COMPILATION SUCCESSFUL

# NGA compilation rules
include $(NGA_HOME)/tools/GNUMake/Make.rules
