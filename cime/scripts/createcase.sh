#!/bin/bash

##########################################################################
# Set up case
##########################################################################
COMPSET=GIAF_LICOM  ## or CIAF_LICOM,AO1850NCLICOM
RES=T62_lt0.1v1       ## or T62_lt0.1v1,ne30_lt1v1
PES=80              ## or lt0.1v1 grid is 960  
DEBUG=FALSE   ## TRUE/FALSE
#STOPOPT=STOP_OPTION=nmonths,STOP_N=1,RESUBMIT=0
#WALLCLOCK=00:30:00  ### 12:00:00

MACHINE=sunway
COMPILER=gnu
# MACHINE=swx86
# COMPILER=swx86
#CESMDIR=/home/export/online1/mdt00/shisuan/swgbcm/xukai/GB25/CESM2.2/cesm2.2
#CASEROOT=${CESMDIR}/cime/cases

##########################################################################
#standard case name
CASE=${COMPSET}.${RES}.${PES}.01


##########################################################################
# create CASE

#cd $CASEROOT

./create_newcase --compset ${COMPSET} --res ${RES} --machine ${MACHINE} --compiler ${COMPILER} --case ${CASE} --run-unsupported --pecount ${PES} || { echo 'create_newcase failed' ; exit 1; }

echo ${CASE} || cd ${CASE}

###########################################################################
## modify case environment variables with xmlchange and 
## Specific mods for high res case
#./xmlchange ATM_NCPL=48 
#./xmlchange OCN_NCPL=24 
## ./xmlchange MAX_TASKS_PER_NODE=16,MAX_MPITASKS_PER_NODE=16
## ./xmlchange PIO_BUFFER_SIZE_LIMIT=512000000
## ./xmlchange PIO_NETCDF_FORMAT=64bit_data
#
###########################################################################
## Generic changes for case specific run (start time, length, etc)
#./xmlchange DEBUG=$DEBUG
#./xmlchange ${STOPOPT}
#./xmlchange REST_OPTION=ndays,REST_N=2
#./xmlchange DOUT_S=FALSE
#./xmlchange RUN_STARTDATE=0001-01-01
#
###########################################################################
## Case setup
#./*.setup
#
###########################################################################
## create user_nl namelists
#cat >> user_nl_cam << EOFcamse
#! Users should add all user specific namelist changes below in the form of 
#! namelist_var = new_namelist_value 
#
#ncdata = '/home/export/online3/cesm06/jflfy/CESM2.2/components/cam/tools/interpic_back/cami-ne480-210207-netcdf4-remap.nc'
#bnd_topo = '/home/export/online3/cesm06/CESM2-inputdata/inputdata/atm/cam/topo/se/ne480pg3_nc3000_Co120_Fi001_MulG_PF_Nsw085.nc'
#
#cpl_cdf64 = .true.
#se_statefreq      = 144                                                
#se_hypervis_subcycle_sponge = 3
#se_nu_top=1.25E5
#EOFcamse
#
#cat >> user_nl_clm << EOFcamse
#! Users should add all user specific namelist changes below in the form of
#! namelist_var = new_namelist_value
#!
#! EXCEPTIONS:
#! Set use_cndv           by the compset you use and the CLM_BLDNML_OPTS -dynamic_vegetation setting
#! Set use_vichydro       by the compset you use and the CLM_BLDNML_OPTS -vichydro           setting
#! Set use_cn             by the compset you use and CLM_BLDNML_OPTS -bgc  setting
#! Set use_crop           by the compset you use and CLM_BLDNML_OPTS -crop setting
#! Set spinup_state       by the CLM_BLDNML_OPTS -bgc_spinup      setting
#! Set co2_ppmv           with CCSM_CO2_PPMV                      option
#! Set fatmlndfrc         with LND_DOMAIN_PATH/LND_DOMAIN_FILE    options
#! Set finidat            with RUN_REFCASE/RUN_REFDATE/RUN_REFTOD options for hybrid or branch cases
#!                        (includes $inst_string for multi-ensemble cases)
#!                        or with CLM_FORCE_COLDSTART to do a cold start
#!                        or set it with an explicit filename here.
#! Set maxpatch_glcmec    with GLC_NEC                            option
#! Set glc_do_dynglacier  with GLC_TWO_WAY_COUPLING               env variable
#!----------------------------------------------------------------------------------
#
#fsurdat = '/home/export/online3/cesm06/CESM2-inputdata/inputdata/lnd/clm2/surfdata_map/release-clm5.0.30/surfdata_ne480np4_hist_78pfts_CMIP6_simyr1850_cdf5_c210202.nc'
#EOFcamse
###########################################################################
## build and submit case
#
## ./*.build
## ./*.submit
#

###########################################################################
## rebuild of cpe code and submit case 
#
## cp ../../../components/licom/SourceMods_cpe/licom_cpe/* SourceMods/src.licom/ && cp ../../../components/licom/SourceMods_cpe/src.cpe/common_type_mod.F90 SourceMods/src.licom/
## cp -r ../../../components/licom/SourceMods_cpe/src.cpe SourceMods/
#
### vim Macros.make, change as followed :
#  MPIFC:= /home/export/online1/mdt00/shisuan/swgbcm/wshuang/tools/xfort/xfort.py
### add in Line 49 as followed :
#  FFLAGS := $(FFLAGS)  -DLICOM3_CPE
### change in Line 74 as followed:
#  ifeq ($(MODEL),mct)
#    LDFLAGS := $(LDFLAGS)  -L/home/export/online1/mdt00/shisuan/swgbcm/xukai/GB25//CESM2.2/sw9gcc-libs/lib -llapack -lnetcdff -lnetcdf -lhdf5_hl -lhdf5 -lz  -Wl,--no-relax -lm_slave
#  else
#    LDFLAGS := $(LDFLAGS)  -L/home/export/online1/mdt00/shisuan/swgbcm/xukai/GB25//CESM2.2/sw9gcc-libs/lib -lcpe -llapack -lnetcdff -lnetcdf -lhdf5_hl -lhdf5 -lz  -Wl,--no-relax -lm_slave
#  endif
#
### vim Macros.cmake, change as followed:
#  set(MPIFC "/home/export/online1/mdt00/shisuan/swgbcm/wshuang/tools/xfort/xfort.py")
### add in Line 49 as followd :
#  set(FFLAGS "${FFLAGS}  -DLICOM3_CPE")
### change in Line 74, add -lcpe as followed : 
#  set(LDFLAGS "${LDFLAGS}  -L/home/export/online1/mdt00/shisuan/swgbcm/xukai/GB25//CESM2.2/sw9gcc-libs/lib -lcpe -llapack -lnetcdff -lnetcdf -lhdf5_hl -lhdf5 -lz  -Wl,--no-relax -lm_slave")
#
## ./*.build -clean
## ./*.build 
## ./*.submit
#

