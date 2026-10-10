# Required version
cmake_minimum_required(VERSION 3.12)

######################################################################
#
# Project Version
#
######################################################################
project(LiCOM LANGUAGES C CXX Fortran)
set(VERSION_MAJOR "3")
set(VERSION_MINOR "1")
set(VERSION_BugFix "0")
set(PROJECT_VERSION ${VERSION_MAJOR}.${VERSION_MINOR}.${VERSION_BugFix})

######################################################################
#
# set up options
#
######################################################################
option(ENABLE_SINGLE    "Enable mix precision" OFF)
option(ENABLE_DOUBLE    "Enable double precision library" ON)
option(ARCH_SUNWAY_MPE  "Enable SUNWAY architecture compile" OFF)
option(ARCH_SUNWAY_CPE  "Enable SUNWAY architecture compile" ON)
option(ARCH_X86CPU      "Enable x86 cpu architecture compile" OFF)
option(ARCH_CUDA        "Enable NVIDIA GPU architecture compile" OFF)

######################################################################
#
# set up compiler and flags
#
######################################################################

set(LICOM_SOURCES)
set_property(GLOBAL PROPERTY SUNWAY_SOURCES)
set_property(GLOBAL PROPERTY SUNWAY_HOST)
set_property(GLOBAL PROPERTY FORTRAN_FILE)
enable_language(Fortran)
include_directories(SYSTEM ${CMAKE_CURRENT_SOURCE_DIR})

# GPTL
include_directories(/home/export/online1/mdt00/shisuan/swiap/local/sw/gptl/include)
link_directories(/home/export/online1/mdt00/shisuan/swiap/local/sw/gptl/lib)

set(NETCDF_PATH "/usr/sw/yyzlib/local-saca")
set(SWLU_PATH "/home/export/online1/mdt00/shisuan/swiap/xukai/swlu")
set(INC "-I.  -I${NETCDF_PATH}/include -I/usr/sw/mpi/mpi_current/include -I${SWLU_PATH}/include")
set(LIB "-L. -L${NETCDF_PATH}/lib  -lnetcdff -lnetcdf  -lhdf5_hl -lhdf5 -ldl -lz")

set(HEAD_FILE
    def-undef.h
)

set_property(GLOBAL APPEND PROPERTY FORTRAN_FILE
        time.F90
        shr_kind_mod.F90
        shr_log_mod.F90
        shr_mpi_mod.F90
        shr_const_mod.F90
        shr_sys_mod.F90
        precision_mod.F90
        param_mod.F90
        sunway_mod.F90
        sunway_init.F90
        msg_mod.F90
        LICOM_Error_mod.F90
        POP_BlocksMod.F90
        POP_CommMod.F90
        POP_SpaceCurveMod.F90
        POP_DistributionMod.F90
        POP_FieldMod.F90
        POP_GridHorzMod.F90
        POP_ReductionsMod.F90
        POP_HaloMod.F90
        control_mod.F90
        broadcast.F90
        blocks.F90
        distribution.F90
        domain.F90
        pconst_mod.F90
        pmix_mod.F90
        constant_mod.F90
        cdf_mod.F90
        module_add.F90
        MPIIO-module.F90
        gather_scatter.F90
        global_reductions.F90
        grid.F90
        tracer_mod.F90
        forc_mod.F90
        dyn_mod.F90
        work_mod.F90
        licom_drv.F90
        boundary.F90
        operators.F90
        output_mod.F90
        grids.F90
        buf_mod.F90
        hmix_del2.F90
        hmix_del4.F90
        diag_mod.F90
        dailyavgini.F90
        inirun.F90
        isopyc_mod.F90
        accumm_daily_avg.F90
        accumm_3hour.F90
        accumm_3hour_avg.F90
        accumm_daily.F90
        avg3hourini.F90
        smuvh.F90
        barotr.F90
        bclinc.F90
        accumm.F90
        advection.F90
        tracer.F90
        energy.F90
        ssave-cdf-month.F90
        ssave-cdf-instant.F90
        ssave-cdf-lonlat.F90
	      ssave-cdf-instantstep.F90
        ssave-cdf-budgetdailyavg.F90
        ssave-cdf-dailyavg.F90
        ssave-cdf-monthavg.F90
        ssave-cdf-avg3hour.F90
        readyt.F90
        icesnow.F90
        convadj.F90
        addps.F90
        nextstep.F90
        vinteg.F90
        utils.F90
        isopyi.F90
        isopyc.F90
        rdriver.F90
        setidealvalue.F90
        upwell.F90
        canuto_mod_2002.F90
        density.F90
        dens.F90
        thermalexpansion.F90
        yy00.F90
        invtri.F90
        mm00.F90
        k1_3.F90
        k2_3.F90
        k3_123.F90
        isoadv.F90
        isoflux.F90
        jra_daily.F90
        gfs_forcing.F90
        core_daily_ideal.F90
        canuto_2002.F90
        readyc.F90
	    sw_absor.F90
	    sw_mod.F90
	    canuto_ini_2002.F90
	    canuto_2010_mod.F90
        var_check.F90
        tidal_force.F90)

if(ARCH_SUNWAY_CPE)
    # add sub directory
    add_subdirectory(sunway)
    add_subdirectory(host-sunway)

    set (CMAKE_C_COMPILER "/usr/sw/mpi/mpi_current/bin/mpicc")
    set (CMAKE_CXX_COMPILER "/usr/sw/mpi/mpi_current/bin/mpic++")
    set(CMAKE_Fortran_COMPILER "/usr/sw/mpi/mpi_current/bin/mpif90")
    set(SUNWEY_HOST_FLAG "-O3 -mhost  -mftz -mieee -fPIC -funroll-loops --param bbss_data_min_size=1 -g  -I/usr/sw/swgcc/swgcc710-tools-SEA/shared_include -I/usr/sw/mpi/mpi_current/include -I/usr/sw/penv/include  -faddress_align=64 ${INC}")
	set(SUNWEY_SLAVE_FLAG "-mslave -msimd -g -O3 -fPIC -faddress_align=64 -lswperf  -I/usr/sw/swgcc/swgcc710-tools-SEA/shared_include -I/usr/sw/penv/include ${INC}")
    set(SUNWEY_HYBRID_FLAG "-mhybrid -static")
    set(LDFLAGS "-lm_slave -L/usr/sw/swgcc/swgcc710-tools-SEA-1307/usr/sw_64sw6a-sunway-linux-gnu/lib -L/usr/sw/penv/lib -L${SWLU_PATH}/lib -lswlu_mpi -lswperf ")

    set(Fortran_HOST_FLAGS "-fdefault-double-8 -fdefault-real-8  -fconvert='big-endian' -ffree-line-length-none   ${SUNWEY_HOST_FLAG}")
 #   set(CMAKE_C_FLAGS "${SUNWEY_HOST_FLAG}")
    get_property(SUNWAY_SOURCES GLOBAL PROPERTY SUNWAY_SOURCES)
    set_source_files_properties(${SUNWAY_SOURCES} PROPERTIES COMPILE_FLAGS "${SUNWEY_SLAVE_FLAG}")
    set_property(GLOBAL APPEND PROPERTY SUNWAY_HOST 
	    ${CMAKE_CURRENT_SOURCE_DIR}/clock0f.c
	    ${CMAKE_CURRENT_SOURCE_DIR}/write_file.c
        ${CMAKE_CURRENT_SOURCE_DIR}/swlu_wrapper.c
        ${CMAKE_CURRENT_SOURCE_DIR}/gptl_timer.c
        ${CMAKE_CURRENT_SOURCE_DIR}/read_file.c    )
    get_property(SUNWAY_HOST GLOBAL PROPERTY SUNWAY_HOST)
    message("-- SUNWAY_HOST='${SUNWAY_HOST}'")
    set_source_files_properties(${SUNWAY_HOST} PROPERTIES COMPILE_FLAGS "${SUNWEY_HOST_FLAG}")
    get_property(FORTRAN_FILE GLOBAL PROPERTY FORTRAN_FILE)
    set_source_files_properties(${FORTRAN_FILE} PROPERTIES COMPILE_FLAGS "${Fortran_HOST_FLAGS}")
    set(CMAKE_EXE_LINKER_FLAGS " ${SUNWEY_HYBRID_FLAG} ")
    SET(MPI_Fortran_LINK_FLAGS "${CMAKE_EXE_LINKER_FLAGS}")
    message("-- MPI_Fortran_LINK_FLAGS='${MPI_Fortran_LINK_FLAGS}'")
    list(APPEND LICOM_SOURCES ${HEAD_FILE} ${FORTRAN_FILE} ${SUNWAY_HOST} ${SUNWAY_SOURCES})
endif()

add_executable(licom3 ${LICOM_SOURCES})
target_link_libraries(licom3 PRIVATE "${LDFLAGS} ${LIB}" gptl)
#target_link_options(Licom PRIVATE "${LDFLAGS}")
#target_link_options(Licom PRIVATE "${LIB}")
#set_target_properties(Licom PROPERTIES LINK_FLAGS "${SUNWEY_HYBRID_FLAG} ${LDFLAGS} ${LIB}" )

