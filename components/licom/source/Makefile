#INC  = -I.  -I/usr/sw/swfile/bak/include -I/usr/sw/mpi/mpi_current/include/include
netcdf_path=/usr/sw/yyzlib/local-saca
SWLU_PATH := /home/export/online1/mdt00/shisuan/swiap/xukai/swlu

#INC  = -I.  -I${netcdf_path}/include -I/usr/sw/mpi/mpi_20210609_SEA/include
INC  = -I.  -I${netcdf_path}/include -I/usr/sw/mpi/mpi_current/include -I$(SWLU_PATH)/include -Ibuild
#LIB  = -L. -L/usr/sw/swfile/bak/lib -lnetcdf -lnetcdff
LIB  = -L. -L${netcdf_path}/lib  -lnetcdff -lnetcdf  -lhdf5_hl -lhdf5 -ldl -lz
#FFLAGS = -O2 -g -traceback -mcmodel=large -convert big_endian -assume byterecl
#FFLAGS  = $(CPPDEFS) -O2 -r8 -i4 -free -fp-model precise -convert big_endian -assume byterecl -no-vec -ftz
#FFLAGS  = $(CPPDEFS) -O0 -g -r8 -i4 -free -fp-model precise -convert big_endian -assume byterecl -no-vec -ftz -traceback
#FFLAGS = $(CPPDEFS) -O2 -ffree-form -ffree-line-length-none -fconvert=big-endian  -mcmodel=large 
FFLAGS = $(CPPDEFS) -O2  -fdefault-double-8 -fdefault-real-8  -fconvert='big-endian' -ffree-line-length-none  -mieee-with-inexact -mftz -faddress_align=64
#FFLAGS = -O2 -mcmodel=large -shared-intel -convert big_endian -assume byterecl
#FFLAGS = -O2 -mcmodel=large -shared-intel -convert big_endian -assume byterecl -fp-model precise -fp-speculation=safe -no-vec 
#FC  = /usr/sw/mpi/mpi_20210609_SEA/bin/mpif90 $(FFLAGS) 
FC  = mpif90 
#FC  = mpiifort -axAVX -m64 -mcmodel=medium $(FFLAGS)
#OBJS = mk_data_lhl1.o 
#LIB = -L/work1/jjr/soft/netcdf/lib -lnetcdf  -L/soft/mpi/openmpi/1.6.5/intel/lib -lmpi  
#FC  = $(MPIFC)

CC_HOST_FLAGS := -mhost  -mftz -mieee -g  -I/usr/sw/swgcc/swgcc710-tools-SEA/shared_include -I/usr/sw/mpi/mpi_current/include -faddress_align=64 -I/usr/sw/penv/include
CC_SLAVE_FLAGS := -mslave -mieee -msimd -g -O3 -faddress_align=64 -lswperf -I/usr/sw/penv/include -I/usr/sw/swgcc/swgcc710-tools-SEA/shared_include
CC_HYBRID_FLAGS := -mhybrid 
#-static 
#-lswperf
LDFLAGS = -L$(SWLU_PATH)/lib -lswlu_mpi -lm_slave -L/usr/sw/swgcc/swgcc710-tools-SEA-1307/usr/sw_64sw6a-sunway-linux-gnu/lib -L/usr/sw/penv/lib -lswperf 


#/usr/sw/penv/lib/libswperf.a

CC  = mpicc
PRG = ../run/licom3.exe
#PRG = ../run-10kmcpu/licom3.exe
#include ./Macros
OBJ =time.o shr_kind_mod.o shr_log_mod.o shr_mpi_mod.o shr_const_mod.o shr_sys_mod.o precision_mod.o param_mod.o msg_mod.o LICOM_Error_mod.o POP_BlocksMod.o \
      POP_CommMod.o POP_SpaceCurveMod.o POP_DistributionMod.o  POP_FieldMod.o  POP_GridHorzMod.o  POP_ReductionsMod.o POP_HaloMod.o \
			control_mod.o broadcast.o blocks.o distribution.o domain.o pconst_mod.o pmix_mod.o constant_mod.o cdf_mod.o gather_scatter.o global_reductions.o var_check.o \
      grid.o tracer_mod.o forc_mod.o  dyn_mod.o work_mod.o boundary.o operators.o output_mod.o grids.o buf_mod.o hmix_del2.o \
			hmix_del4.o diag_mod.o dailyavgini.o inirun.o isopyc_mod.o accumm_daily_avg.o accumm_daily.o smuvh.o barotr.o bclinc.o accumm.o advection.o tracer.o energy.o \
			ssave-cdf-month.o \
			ssave-cdf-instant.o  ssave-cdf-lonlat.o  readyt.o icesnow.o convadj.o addps.o nextstep.o vinteg.o utils.o isopyi.o isopyc.o rdriver.o setidealvalue.o \
			upwell.o  canuto_mod_2002.o  density.o dens.o thermalexpansion.o yy00.o invtri.o mm00.o k1_3.o k2_3.o k3_123.o isoadv.o isoflux.o jra_daily.o core_daily_ideal.o canuto_2002.o readyc.o canuto_ini_2002.o \
			licom_drv.o

SLAVEOBJ = readyt_cu_1_slave.o dens_cu_slave.o readyt_cu_2_slave.o readyt_cu_2_2_slave.o thermal_slave.o \
	   readyt_3_slave.o readyt_3_2_slave.o readyt_4_slave.o readyt_5_slave.o convadj_cu_slave.o readyc_cu_slave.o \
	   readyc_upwell_slave.o readyc_turb_cu_slave.o barotr1_cu_slave.o barotr1_cu1_slave.o barotr2_cu1_slave.o \
	   bhdel4_div_slave.o bhdel4_grad_slave.o bhdel4_zcurl_slave.o barotr2_udel4_0_slave.o barotr2_udel4_1_slave.o \
	   barotr2_tide_slave.o barotr3_cu_slave.o barotr3_tdel4_1_slave.o barotr3_cu1_slave.o barotr4_cu_slave.o bclinc11_cu_slave.o \
	   bclinc12_cu_slave.o bclinc2_cu_slave.o bclinc3_cu_slave.o bclinc4_cu_slave.o bclinc5_cu_slave.o \
	   bclinc71_cu_slave.o bclinc72_cu_slave.o  bclinc8_cu_slave.o hdel4_div_slave.o hdel4_grad_slave.o \
	   hdel4_zcurl_slave.o hdiffu_del4_0_slave.o hdiffu_del4_1_slave.o icesnowsw.o tracersw.o upwellsw.o \
	   invtrit_slave.o tracer5_cu_slave.o tracer6_cu_slave.o tracer7_cu_slave.o tracer_9_slave.o invtritsw.o\
	   tracer_del4_1_slave.o tracer_del4_2_slave.o tracer_advection_slave.o tracer_upwell_slave.o bclinc_vinteg_slave.o

HOST_C_OBJ := clock0f.o swlu_wrapper.o read_file.o athread_core.o

BASESOURCE:=$(basename $(OBJ))
OBJECTS:=$(addsuffix .F90.o, $(BASESOURCE))
OBJECTS:=$(addprefix build/, $(OBJECTS))

BASESOURCE:=$(basename $(HOST_C_OBJ))
HOST_C_OBJECTS:=$(addsuffix .host.o, $(BASESOURCE))
HOST_C_OBJECTS:=$(addprefix build/, $(HOST_C_OBJECTS))

BASESOURCE:=$(basename $(SLAVEOBJ))
SLAVE_C_OBJECTS:=$(addsuffix .slave.o, $(BASESOURCE))
SLAVE_C_OBJECTS:=$(addprefix build/, $(SLAVE_C_OBJECTS))


RM:=rm -f
MKDIR:=mkdir -p
TEST:=test

all: bld_dir $(PRG)

bld_dir:
	$(shell $(TEST) -d build || $(MKDIR) -p build)

$(PRG):$(OBJECTS) $(SLAVE_C_OBJECTS) $(HOST_C_OBJECTS)
	$(FC) $(INC)  -o $@ $(OBJECTS) $(SLAVE_C_OBJECTS) $(HOST_C_OBJECTS) $(LIB) $(LDFLAGS) 
		      
.SUFFIXES: .F90 .F .f90 .f .c .cpp .o

build/%.host.o:%.c  
	$(CC) $(CC_HOST_FLAGS) $(INC) -c $< -o $@  

build/%.slave.o : ./sunway/%.c 
	$(CC) $(CC_SLAVE_FLAGS) $(INC) -o $@ -c $<

build/%.F90.o: %.F90				  
	$(FC) $(FFLAGS) $(CC_HOST_FLAGS) $(INC) -c $< -fsyntax-only -J build
	$(FC) $(FFLAGS) $(CC_HOST_FLAGS) $(INC) -c $< -J build -o $@
clean:  
	rm -f $(OBJECTS) $(SLAVE_C_OBJECTS) $(HOST_C_OBJECTS)  $(PRG) build/*.mod
	-$(TEST) -d build && $(RM) -r build
