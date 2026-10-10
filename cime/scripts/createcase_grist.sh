#!/bin/bash

##########################################################################
# Set up case
##########################################################################
COMPSET=AO1850NCLICOMGRIST 
RES=grist8_lt0.1v2  ## or grist6_lt1v1, grist9_lt0.05v2, grist10_lt0.03v2, grist11_lt0.02v2
PES=1560              ## or grist6_lt1v1 is 80, grist9_lt0.05v2 is 4620, grist10_lt0.03v2 is 10440, grist11_lt0.02v2 is 28800(or 19591)
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
## config of grist case
#
## cd ${CASE} 
## vim env_build.xml  ! change Line 352 to <entry id="CICE_CONFIG_OPTS" value=" -phys cice5 -ntr_aero 0 ">
## vim env_run.xml    ! change Line 113 to <entry id="RUN_STARTDATE" value="2023-07-18">
## cd ${CASE}/run
## mkdir result       ! this subdirectory is used for grist output
## ./prepare_ocn_lt*_link.csh
## ./copy_case.sh G*
## cd ..
## ./*.build -clean
## ./*.build 
## ./*.submit
#

