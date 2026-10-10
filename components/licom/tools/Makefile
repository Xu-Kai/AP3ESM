# General settings
FC_SW = mpif90

# Flags for SW Architecture
FFLAGS_SW := -O2 -mieee -msimd
netcdf_path=/usr/sw/yyzlib/local-saca
INCL_SW := -I. -I${netcdf_path}/include -I/usr/sw/mpi/mpi_current/include
LDFLAGS_SW := -L. -L${netcdf_path}/lib
LIBS_SW := -lnetcdff -lnetcdf -lhdf5_hl -lhdf5 -ldl -lz

# The build/run direction
OUTOBJDIR = ./run
# Make sure the run direction is exist
$(shell mkdir -p $(OUTOBJDIR))
$(shell cp submit_sw $(OUTOBJDIR)/)

# Source Files
F_SOURCES := preprocess.F90
OBJECTS := $(F_SOURCES:.F90=.o)

# Executables
EXECUTABLE := $(OUTOBJDIR)/invertnc

# Targets
.PHONY: all clean sw

all: sw

sw: $(EXECUTABLE)

$(EXECUTABLE): $(OBJECTS)
	$(FC_SW) $(OBJECTS) $(LDFLAGS_SW) -o $@ $(LIBS_SW)

# Object files compilation SW
%.o: %.F90
	$(FC_SW) $(FFLAGS_SW) $(INCL_SW) -c $< -o $@

# Cleaning up
clean:
	rm -f *.o  *.mod $(EXECUTABLE)
