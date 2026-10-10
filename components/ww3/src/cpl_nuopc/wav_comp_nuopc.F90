module wav_comp_nuopc

  !/ ------------------------------------------------------------------- /
  !/
  !/                  +-----------------------------------+
  !/                  | WAVEWATCH III           NOAA/NCEP |
  !/                  |           H. L. Tolman            |
  !/                  |                        FORTRAN 90 |
  !/                  | Last update :         29-May-2009 |
  !/                  +-----------------------------------+
  !/
  !/    Copyright 2009 National Weather Service (NWS),
  !/       National Oceanic and Atmospheric Administration.  All rights
  !/       reserved.  WAVEWATCH III is a trademark of the NWS.
  !/       No unauthorized use without permission.
  !/
  !  1. Purpose :
  !
  !     A generic nuopc interface for WAVEWATCH III
  !     using input fields from CMEPS.
  !
  !  2. Method :
  !
  !     NUOPC component for the actual wave model (W3WAVE).
  !
  !  3. Parameters :
  !
  !     Local parameters.
  !     ----------------------------------------------------------------
  !       NHMAX   I.P.  Maximum number of homogeneous fields.
  !
  !       NDSO    Int.  General output unit number (shell only).
  !       NDSE    Int.  Error output unit number (shell only).
  !       NDST    Int.  Test output unit number (shell only).
  !       FLH     L.A.  Flags for homogeneous fields.
  !       NH      I.A.  Number of times for homogeneous fields.
  !       THO     I.A.  Times of homogeneous fields.
  !       TIME0   I.A.  Starting time.
  !       TIMEN   I.A.  Ending time.
  !     ----------------------------------------------------------------
  !
  !       NDS, NTRACE, ..., see W3WAVE
  !
  !  4. Subroutines used :
  !
  !      Name      Type  Module   Description
  !     ----------------------------------------------------------------
  !      W3NMOD    Subr. W3GDATMD Set nummber of data structures
  !      W3SETG    Subr.   Id.    Point to data structure.
  !      W3NDAT    Subr. W3WDATMD Set nummber of data structures
  !      W3SETW    Subr.   Id.    Point to data structure.
  !      W3NMOD    Subr. W3ADATMD Set nummber of data structures
  !      W3NAUX    Subr.   Id.    Point to data structure.
  !      W3NOUT    Subr. W3ODATMD Set nummber of data structures
  !      W3SETO    Subr.   Id.    Point to data structure.
  !      W3NINP    Subr. W3IDATMD Set nummber of data structures
  !      W3SETI    Subr.   Id.    Point to data structure.
  !
  !      STME21    Subr. W3TIMEMD Print date and time readable.
  !      DSEC21    Func.   Id.    Difference between times.
  !      TICK21    Subr.   Id.    Increment time.
  !
  !      W3INIT    Subr. W3INITMD Wave model initialization.
  !      W3WAVE    Subr. W3WAVEMD Wave model.
  !     ----------------------------------------------------------------
  !
  !  5. Called by :
  !
  !     None, stand-alone program.
  !
  !  6. Error messages :
  !
  !     - Checks on I-O.
  !     - Check on time interval.
  !
  !  7. Remarks :
  !
  !     - A rigourous input check is made in W3INIT.
  !     - See W3WDAS for documentation on the set-up of the data
  !       assimilation.
  !
  !  8. Structure :
  !
  !     ----------------------------------------------------------------
  !
  !     wav_comp_init
  !
  !        0.   Set up data structures.                ( W3NMOD, etc. )
  !        1.   I-O setup.
  !          a  For shell.
  !          b  For WAVEWATCH III.
  !          c  Local parameters.
  !        2.   Define input fields
  !        3.   Set time frame.
  !        4.   Define output
  !          a  Loop over types, do
  !        +--------------------------------------------------------+
  !        | b    Process standard line                             |
  !        | c    If type 1: fields of mean wave parameters         |
  !        | d    If type 2: point output                           |
  !        | e    If type 3: track output                           |
  !        | f    If type 4: restart files                          |
  !        | g    If type 5: boundary output                        |
  !        | h    If type 6: separated wave fields                  |
  !        +--------------------------------------------------------+
  !        5.   Initialzations
  !          a  Wave model.                              ( W3INIT )
  !          d  Set field times.
  !
  !     wav_comp_run
  !
  !        7.   Run model with input
  !             Do until end time is reached
  !        +--------------------------------------------------------+
  !        | a  Determine next time interval and input fields.      |
  !        |   1  Preparation                                       |
  !        |      Loop over input fields                            |
  !        | +------------------------------------------------------|
  !        | | 2  Check if update is needed                         |
  !        | | 4  Update next ending time                           |
  !        | +------------------------------------------------------|
  !        | b  Run wave model.                          ( W3WAVE ) |
  !        | c  If requested, data assimilation.         ( W3WDAS ) |
  !        | d  Final output if needed.                  ( W3WAVE ) |
  !        | e  Check time                                          |
  !        +--------------------------------------------------------+
  !
  !     wav_comp_fin
  !
  !     ----------------------------------------------------------------
  !
  !  9. Switches :
  !
  !       !/SHRD  Switch for shared / distributed memory architecture.
  !       !/DIST  Id.
  !       !/MPI   Id.
  !
  !       !/LLG   Spherical grid.
  !       !/XYG   Cartesian grid.
  !       !/MGW   Moving grid wind correction.
  !       !/MGP   Moving grid propagation correction.
  !
  !       !/T     Enable test output.
  !       !/O7    Echo input homogeneous fields.
  !
  !       !/NCO   NCEP NCO modifications for operational implementation.
  !
  !       !/F90   Timer function included for F90.
  !
  ! 10. Source code :
  !
  !/ ------------------------------------------------------------------- /
  
  use w3gdatmd
  use w3wdatmd              , only: time, w3ndat, w3dimw, w3setw
  use w3adatmd
  use w3idatmd              , only: flags, w3seti, w3ninp
  use w3idatmd              , only: TC0, CX0, CY0, TCN, CXN, CYN
  use w3idatmd              , only: TW0, WX0, WY0, DT0, TWN, WXN, WYN, DTN
  use w3idatmd              , only: TIN, ICEI, TLN, WLEV, HML
  use w3odatmd              , only: w3nout, w3seto, naproc, iaproc, napout, naperr
  use w3odatmd              , only: nogrd, idout, fnmpre, iostyp
  use w3initmd
  use w3wavemd
  use w3iopomd
  use w3timemd
  use w3cesmmd              , only : casename, initfile, rstwr, runtype, histwr, outfreq
  use w3cesmmd              , only : inst_index, inst_name, inst_suffix
  use ESMF
  use NUOPC                 , only : NUOPC_CompDerive, NUOPC_CompSetEntryPoint, NUOPC_CompSpecialize
  use NUOPC                 , only : NUOPC_CompFilterPhaseMap, NUOPC_IsUpdated, NUOPC_IsAtTime
  use NUOPC                 , only : NUOPC_CompAttributeGet, NUOPC_Advertise
  use NUOPC                 , only : NUOPC_SetAttribute, NUOPC_CompAttributeGet, NUOPC_CompAttributeSet
  use NUOPC_Model           , only : model_routine_SS           => SetServices
  use NUOPC_Model           , only : model_label_Advance        => label_Advance
  use NUOPC_Model           , only : model_label_DataInitialize => label_DataInitialize
  use NUOPC_Model           , only : model_label_SetRunClock    => label_SetRunClock
  use NUOPC_Model           , only : model_label_Finalize       => label_Finalize
  use NUOPC_Model           , only : NUOPC_ModelGet, SetVM
  use shr_sys_mod           , only : shr_sys_flush, shr_sys_abort
  use shr_nl_mod            , only : shr_nl_find_group_name
  use shr_mpi_mod           , only : shr_mpi_bcast
  use shr_kind_mod          , only : r8=>shr_kind_r8, i8=>shr_kind_i8, cl=>shr_kind_cl, cs=>shr_kind_cs
  use shr_file_mod          , only : shr_file_getLogUnit, shr_file_setLogUnit, shr_file_getunit
  use shr_cal_mod           , only : shr_cal_ymd2date
  use wav_import_export     , only : advertise_fields, realize_fields, import_fields, export_fields
  use wav_import_export     , only : state_getfldptr
  use wav_shr_methods       , only : chkerr, state_setscalar, state_getscalar, state_diagnose, alarmInit
  use wav_shr_methods       , only : set_component_logging, get_component_instance, log_clock_advance

  implicit none
  private ! except

  public  :: SetServices
  private :: InitializeP0
  private :: InitializeAdvertise
  private :: InitializeRealize
  private :: ModelSetRunClock
  private :: ModelAdvance
  private :: ModelFinalize

  integer :: stdout
  include "mpif.h"

  !--------------------------------------------------------------------------
  ! Private module data
  !--------------------------------------------------------------------------

  character(len=CL)       :: flds_scalar_name = ''
  integer                 :: flds_scalar_num = 0
  integer                 :: flds_scalar_index_nx = 0
  integer                 :: flds_scalar_index_ny = 0
  integer                 :: flds_scalar_index_precip_factor = 0._r8

  logical                 :: masterproc 
  integer     , parameter :: debug = 1
  character(*), parameter :: modName =  "(wav_comp_nuopc)"
  character(*), parameter :: u_FILE_u = &
       __FILE__

!===============================================================================
contains
!===============================================================================

  subroutine SetServices(gcomp, rc)
    type(ESMF_GridComp)  :: gcomp
    integer, intent(out) :: rc

    character(len=*),parameter  :: subname=trim(modName)//':(SetServices) '

    rc = ESMF_SUCCESS
    call ESMF_LogWrite(subname//' called', ESMF_LOGMSG_INFO)

    ! the NUOPC gcomp component will register the generic methods
    call NUOPC_CompDerive(gcomp, model_routine_SS, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! switching to IPD versions
    call ESMF_GridCompSetEntryPoint(gcomp, ESMF_METHOD_INITIALIZE, &
         userRoutine=InitializeP0, phase=0, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! set entry point for methods that require specific implementation
    call NUOPC_CompSetEntryPoint(gcomp, ESMF_METHOD_INITIALIZE, &
         phaseLabelList=(/"IPDv01p1"/), userRoutine=InitializeAdvertise, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call NUOPC_CompSetEntryPoint(gcomp, ESMF_METHOD_INITIALIZE, &
         phaseLabelList=(/"IPDv01p3"/), userRoutine=InitializeRealize, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! attach specializing method(s)

    call NUOPC_CompSpecialize(gcomp, specLabel=model_label_DataInitialize, &
         specRoutine=DataInitialize, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call NUOPC_CompSpecialize(gcomp, specLabel=model_label_Advance, &
         specRoutine=ModelAdvance, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call ESMF_MethodRemove(gcomp, label=model_label_SetRunClock, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call NUOPC_CompSpecialize(gcomp, specLabel=model_label_SetRunClock, &
         specRoutine=ModelSetRunClock, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call NUOPC_CompSpecialize(gcomp, specLabel=model_label_Finalize, &
         specRoutine=ModelFinalize, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call ESMF_LogWrite(subname//' done', ESMF_LOGMSG_INFO)

  end subroutine SetServices

  !===============================================================================

  subroutine InitializeP0(gcomp, importState, exportState, clock, rc)
    type(ESMF_GridComp)   :: gcomp
    type(ESMF_State)      :: importState, exportState
    type(ESMF_Clock)      :: clock
    integer, intent(out)  :: rc
    !-------------------------------------------------------------------------------

    rc = ESMF_SUCCESS

    ! Switch to IPDv01 by filtering all other phaseMap entries

    call NUOPC_CompFilterPhaseMap(gcomp, ESMF_METHOD_INITIALIZE, acceptStringList=(/"IPDv01p"/), rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

  end subroutine InitializeP0

  !===============================================================================

  subroutine InitializeAdvertise(gcomp, importState, exportState, clock, rc)

    ! input/output arguments
    type(ESMF_GridComp)  :: gcomp
    type(ESMF_State)     :: importState, exportState
    type(ESMF_Clock)     :: clock
    integer, intent(out) :: rc

    ! local variables
    character(len=CL) :: logmsg
    logical           :: isPresent, isSet
    character(len=CL) :: cvalue
    character(len=*), parameter :: subname=trim(modName)//':(InitializeAdvertise) '
    !-------------------------------------------------------------------------------

    rc = ESMF_SUCCESS
    call ESMF_LogWrite(subname//' called', ESMF_LOGMSG_INFO)

    !----------------------------------------------------------------------------
    ! advertise fields
    !----------------------------------------------------------------------------

    call NUOPC_CompAttributeGet(gcomp, name="ScalarFieldName", value=cvalue, isPresent=isPresent, isSet=isSet, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    if (isPresent .and. isSet) then
       flds_scalar_name = trim(cvalue)
       call ESMF_LogWrite(trim(subname)//' flds_scalar_name = '//trim(flds_scalar_name), ESMF_LOGMSG_INFO)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
    else
       call shr_sys_abort(subname//'Need to set attribute ScalarFieldName')
    endif

    call NUOPC_CompAttributeGet(gcomp, name="ScalarFieldCount", value=cvalue, isPresent=isPresent, isSet=isSet, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    if (isPresent .and. isSet) then
       read(cvalue, *) flds_scalar_num
       write(logmsg,*) flds_scalar_num
       call ESMF_LogWrite(trim(subname)//' flds_scalar_num = '//trim(logmsg), ESMF_LOGMSG_INFO)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
    else
       call shr_sys_abort(subname//'Need to set attribute ScalarFieldCount')
    endif

    call NUOPC_CompAttributeGet(gcomp, name="ScalarFieldIdxGridNX", value=cvalue, isPresent=isPresent, isSet=isSet, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    if (isPresent .and. isSet) then
       read(cvalue,*) flds_scalar_index_nx
       write(logmsg,*) flds_scalar_index_nx
       call ESMF_LogWrite(trim(subname)//' : flds_scalar_index_nx = '//trim(logmsg), ESMF_LOGMSG_INFO)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
    else
       call shr_sys_abort(subname//'Need to set attribute ScalarFieldIdxGridNX')
    endif

    call NUOPC_CompAttributeGet(gcomp, name="ScalarFieldIdxGridNY", value=cvalue, isPresent=isPresent, isSet=isSet, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    if (isPresent .and. isSet) then
       read(cvalue,*) flds_scalar_index_ny
       write(logmsg,*) flds_scalar_index_ny
       call ESMF_LogWrite(trim(subname)//' : flds_scalar_index_ny = '//trim(logmsg), ESMF_LOGMSG_INFO)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
    else
       call shr_sys_abort(subname//'Need to set attribute ScalarFieldIdxGridNY')
    endif

    call advertise_fields(importState, exportState, flds_scalar_name, rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call ESMF_LogWrite(subname//' done', ESMF_LOGMSG_INFO)

  end subroutine InitializeAdvertise

  !========================================================================

  subroutine InitializeRealize(gcomp, importState, exportState, clock, rc)

    ! input/output variables
    type(ESMF_GridComp)  :: gcomp
    type(ESMF_State)     :: importState
    type(ESMF_State)     :: exportState
    type(ESMF_Clock)     :: clock
    integer, intent(out) :: rc

    ! local variables
    integer, parameter             :: nhmax = 200
    real                           :: a(nhmax,4)
    type(ESMF_DistGrid)            :: distGrid
    type(ESMF_Mesh)                :: Emesh, EmeshTemp
    type(ESMF_VM)                  :: vm
    type(ESMF_Time)                :: ETime
    type(ESMF_TimeInterval)        :: TimeStep
    character(CL)                  :: cvalue
    integer                        :: shrlogunit
    integer                        :: yy,mm,dd,hh,ss
    integer                        :: dtime_sync        ! integer timestep size
    integer                        :: start_ymd         ! start date (yyyymmdd)
    integer                        :: start_tod         ! start time of day (sec)
    integer                        :: stop_ymd          ! stop date (yyyymmdd)
    integer                        :: stop_tod          ! stop time of day (sec)
    integer                        :: ix, iy
    character(CL)                  :: starttype
    integer                        :: unitn            ! namelist unit number
    integer                        :: ndso, ndse, nds(13), ntrace(2), time0(2)
    integer                        :: timen(2), odat(30), nh(4), iprt(6)
    integer                        :: i,j,npts
    integer                        :: ierr
    real, allocatable              :: x(:), y(:)
    integer                        :: n, jsea,isea, ncnt
    integer                        :: ntotal, nlnd
    integer                        :: nlnd_global, nlnd_local
    integer                        :: my_lnd_start, my_lnd_end
    integer, allocatable, target   :: mask_global(:)
    integer, allocatable, target   :: mask_local(:)
    integer, allocatable           :: gindex_lnd(:)
    integer, allocatable           :: gindex_sea(:)
    integer, allocatable           :: gindex(:) 
    logical                        :: flgrd(nogrd), prtfrm, flt
    character(len=23)              :: dtme21
    integer                        :: iam, mpi_comm
    character(len=10), allocatable :: pnames(:)
    character(len=*),parameter :: subname = '(wav_comp_nuopc:InitializeRealize)'
    ! -------------------------------------------------------------------

    namelist /ww3_inparm/ initfile, outfreq

    !--------------------------------------------------------------------
    ! Set up data structures
    !--------------------------------------------------------------------

    call w3nmod ( 1, 6, 6 )
    call w3ndat (    6, 6 )
    call w3naux (    6, 6 )
    call w3nout (    6, 6 )
    call w3ninp (    6, 6 )

    call w3setg ( 1, 6, 6 )
    call w3setw ( 1, 6, 6 )
    call w3seta ( 1, 6, 6 )
    call w3seto ( 1, 6, 6 )
    call w3seti ( 1, 6, 6 )

    !----------------------------------------------------------------------------
    ! Generate local mpi comm
    !----------------------------------------------------------------------------

    call ESMF_GridCompGet(gcomp, vm=vm, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call ESMF_VMGet(vm, mpiCommunicator=mpi_comm, peCount=naproc, localPet=iam, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    iaproc = iam + 1

    !--------------------------------------------------------------------
    ! IO set-up
    !--------------------------------------------------------------------

    ! 1.b For WAVEWATCH III (See W3INIT) ??? ask adrean if i am missing something
    !
    ! The following units are referenced in module w3initmd
    ! NDS(1) ! OUTPUT LOG: General output unit number ("log file") (NDS0)
    ! NDS(2) ! OUTPUT LOG: Error output unit number (NDSE)
    ! NDS(3) ! OUTPUT LOG: Test output unit number (NDST)
    ! NDS(4) ! OUTPUT LOG: Unit for 'direct' output (SCREEN)
    !
    ! NDS(5) ! INPUT: mod_def.ww3 file (model definition) unit number
    ! NDS(9) ! INPUT: unit for read in boundary conditions (based on FLBPI)
    !
    ! The following units are referenced in module w3wavemd for output
    ! NDS( 6) ! OUTPUT DATA: restart(N).ww3 file (model restart) unit number
    ! NDS( 7) ! unit for output for FLOUT(1) flag
    ! NDS( 8) ! unit for output for FLOUT(2) flag
    ! NDS(11) ! unit for output for FLOUT(3) flag
    ! NDS(12) ! unit for output for FLOUT(3) flag
    ! NDS(13) ! unit for output for FLOUT(6) flag
    ! NDS(10) ! unit for output for FLOUT(5) flag

    !----------------------------------------------------------------------------
    ! determine instance information
    !----------------------------------------------------------------------------

    call get_component_instance(gcomp, inst_suffix, inst_index, rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    inst_name = "WAV"//trim(inst_suffix)

    !----------------------------------------------------------------------------
    ! reset shr logging to my log file
    !----------------------------------------------------------------------------

    napout = 1
    naperr = 1

    if (iaproc == napout) then
       masterproc = .true.
    else
       masterproc = .false.
    end if

    call set_component_logging(gcomp, masterproc, stdout, shrlogunit, rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    nds( 1) = stdout
    nds( 2) = stdout
    nds( 3) = stdout
    nds( 4) = stdout
    nds( 5) = shr_file_getunit()
    nds( 6) = shr_file_getunit()
    nds( 7) = shr_file_getunit()
    nds( 8) = shr_file_getunit()
    nds( 9) = shr_file_getunit()
    nds(10) = shr_file_getunit()
    nds(11) = shr_file_getunit()
    nds(12) = shr_file_getunit()
    nds(13) = shr_file_getunit()

    ndso      =  stdout
    ndse      =  stdout
    ntrace(1) =  nds(3)
    ntrace(2) =  10

    ! Redirect share output to wav log
    call shr_file_getLogUnit (shrlogunit)
    call shr_file_setLogUnit (ndso)

    if ( iaproc == napout ) write (ndso,900)
    call shr_sys_flush(ndso)

    !--------------------------------------------------------------------
    ! Initialize run type
    !--------------------------------------------------------------------

    call NUOPC_CompAttributeGet(gcomp, name='start_type', value=cvalue, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    read(cvalue,*) starttype

    if (     trim(starttype) == trim('startup')) then
       runtype = "initial"
       write(ndso,*) 'starttype: initial'
    else if (trim(starttype) == trim('continue') ) then
       runtype = "continue"
       write(ndso,*) 'starttype: continue'
    else if (trim(starttype) == trim('branch')) then
       runtype = "branch"
       write(ndso,*) 'starttype: branch'
    end if

    if ( iaproc == napout) then
       write(ndso,*) trim(subname),' inst_name   = ',trim(inst_name)
       write(ndso,*) trim(subname),' inst_index  = ',inst_index
       write(ndso,*) trim(subname),' inst_suffix = ',trim(inst_suffix)
    endif

    !--------------------------------------------------------------------
    ! Define input fields
    !--------------------------------------------------------------------

    flags = .false.
    ! QL, 150525, flags for passing variables from coupler to ww3,
    !             lev, curr, wind, ice and mixing layer depth on
    flags(1:5) = .true.
    !      flags(1:4) = .true.   !changed by Adrean (lev,curr,wind,ice on)
    !      flags(3:4) = .true.   !changed by Adrean (wind,ice on)

    !--------------------------------------------------------------------
    ! Set time frame
    !--------------------------------------------------------------------

    ! TIME0 = from ESMF clock
    ! NOTE - are not setting TIMEN here

    if ( iaproc == napout ) write (ndso,930)
    call shr_sys_flush(ndso)

    ! Initial run or restart run
    if ( runtype == "initial") then
       call ESMF_ClockGet( clock, startTime=ETime, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
    else
       call ESMF_ClockGet( clock, currTime=ETime, rc=rc )
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
    endif
    call ESMF_TimeGet( ETime, yy=yy, mm=mm, dd=dd, s=start_tod, rc=rc )
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call shr_cal_ymd2date(yy, mm, dd, start_ymd)

    hh = start_tod/3600
    mm = (start_tod - (hh * 3600))/60
    ss = start_tod - (hh*3600) - (mm*60)

    time0(1) = start_ymd
    time0(2) = hh*10000 + mm*100 + ss

    call ESMF_ClockGet( clock, stopTime=ETime, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call ESMF_TimeGet( ETime, yy=yy, mm=mm, dd=dd, s=stop_tod, rc=rc )
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call shr_cal_ymd2date(yy, mm, dd, stop_ymd)

    hh = stop_tod/3600
    mm = (stop_tod - (hh * 3600))/60
    ss = stop_tod - (hh*3600) - (mm*60)

    timen(1) = stop_ymd
    timen(2) = hh*10000 + mm*100 + ss

    call stme21 ( time0 , dtme21 )
    if ( iaproc .eq. napout ) write (ndso,931) dtme21
    call shr_sys_flush(ndso)
    time = time0

    !--------------------------------------------------------------------
    ! Define output type and fields
    !--------------------------------------------------------------------

    iostyp = 1        ! gridded field
    write (ndso,940) 'no dedicated output process, any file system '
    call shr_sys_flush(ndso)

    ! Actually will need a new restart flag - since all of the ODAT
    ! should be set to 0 - since they are initializated in w3initmd
    ! ODAT    I.A.   I   Output data, five parameters per output type
    !                          1 YYYMMDD for first output.
    !                          2 HHMMSS for first output.
    !                          3 Output interval in seconds.
    !                          4 YYYMMDD for last output.
    !                          5 HHMMSS for last output.
    !                     1-5  Data for OTYPE = 1; gridded fields.
    !                     6-10 Id.  for OTYPE = 2; point output.
    !                    11-15 Id.  for OTYPE = 3; track point output.
    !                    16-20 Id.  for OTYPE = 4; restart files.
    !                    21-25 Id.  for OTYPE = 5; boundary data.
    ! FLGRD   L.A.   I   Flags for gridded output.
    ! NPT     Int.   I   Number of output points
    ! X/YPT   R.A.   I   Coordinates of output points.
    ! PNAMES  C.A.   I   Output point names.

    do j=1, 6
       odat(5*(j-1)+3) = 0
    end do

    ! get coupling interval
    call ESMF_ClockGet( clock, timeStep=timeStep, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call ESMF_TimeIntervalGet( timeStep, s=dtime_sync, rc=rc )
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! Hardwire gridded output for now
    ! first output time stamp is now read from file
    ! QL, 150525, 1-5 for history files, 16-20 for restart files
    !     150629, restart output interval is set to the total time of run
    !     150823, restart is taken over by rstwr
    !     160601, output interval is set to coupling interval, so that
    !             variables calculated in W3IOGO could be updated at
    !             every coupling interval
    odat(1) = time(1)     ! YYYYMMDD for first output
    odat(2) = time(2)     ! HHMMSS for first output
    odat(3) = dtime_sync  ! output interval in sec ! changed by Adrean
    odat(4) = 99990101    ! YYYYMMDD for last output
    odat(5) = 0           ! HHMMSS for last output
    odat(16) = time(1)    ! YYYYMMDD for first output
    odat(17) = time(2)    ! HHMMSS for first output
    odat(18) = dtime_sync ! output interval in sec
    odat(19) = 99990101   ! YYYYMMDD for last output
    odat(20) = 0          ! HHMMSS for last output

    ! Output Type 1: fields of mean wave parameters gridded output

    flgrd( 1) = .false. !   1. depth (m)
    flgrd( 2) = .false. !   2. mean current vel (vec, m/s)
    flgrd( 3) = .true.  !   3. mean wind vel (vec, m/s)
    flgrd( 4) = .false. !   4. air-sea temp diff (deg C)
    flgrd( 5) = .false. !   5. skin friction vel (scalar, m/s)
    flgrd( 6) = .true.  !   6. significant wave height (m)
    flgrd( 7) = .false. !   7. mean wave length (m)
    flgrd( 8) = .true.  !   8. mean wave period (Tn1, s)
    flgrd( 9) = .true.  !   9. mean wave dir (deg: met conv)
    flgrd(10) = .false. !  10. mean dir spread (deg: )
    flgrd(11) = .false. !  11. peak freq (Hz)
    flgrd(12) = .false. !  12. peak dir (deg: )
    flgrd(13) = .false. !  13. peak freq of wind-sea part
    flgrd(14) = .false. !  14. wind-sea dir (deg: met conv)
    flgrd(15) = .false. !  15. wave height of partitions
    flgrd(16) = .false. !  16. peak periods of partitions
    flgrd(17) = .false. !  17. peak wave length of partitions
    flgrd(18) = .false. !  18. mean dir of partitions
    flgrd(19) = .false. !  19. dir spread of partitions
    flgrd(20) = .false. !  20. wind-sea frac of partitions
    flgrd(21) = .false. !  21. wind-sea frac of entire spec
    flgrd(22) = .false. !  22. number of partitions
    flgrd(23) = .false. !  23. average time step (s)
    flgrd(24) = .false. !  24. cut-off freq (Hz)
    flgrd(25) = .false. !  25. ice concentration (frac)
    flgrd(26) = .false. !  26. water level (m?)
    flgrd(27) = .false. !  27. near-bottom rms exclusion amp
    flgrd(28) = .false. !  28. near-bottom rms orbital vel
    flgrd(29) = .false. !  29. radiation stresses
    flgrd(30) = .false. !  30. user defined (1)
    flgrd(31) = .false. !  31. user defined (2)

    ! QL, 150525, new output
    flgrd(32) = .false. !  32. Stokes drift at z=0
    flgrd(33) = .false. !  33. Turbulent Langmuir number (La_t)
    flgrd(34) = .false. !  34. Langmuir number (La_Proj)
    flgrd(35) = .false. !  35. Angle between wind and LC direction
    flgrd(36) = .false. !  36. Depth averaged Stokes drift (0-H_0.2ML)
    flgrd(37) = .false. !  37. Langmuir number (La_SL)
    flgrd(38) = .false. !  38. Langmuir number (La_SL,Proj)
    flgrd(39) = .false. !  39. Enhancement factor with La_SL,Proj

    if ( iaproc .eq. napout ) then
       flt = .true.
       do i=1, nogrd
          if ( flgrd(i) ) then
             if ( flt ) then
                write (ndso,1945) idout(i)
                flt    = .false.
             else
                write (ndso,1946) idout(i)
             end if
          end if
       end do
       if ( flt ) write (ndso,1945) 'no fields defined'
    end if
    call shr_sys_flush(ndso)

    !--------------------------------------------------------------------
    ! Wave model initializations
    !--------------------------------------------------------------------

    ! Notes on ww3 initialization:
    ! ww3 read initialization occurs in w3iors (which is called by initmd)
    ! For a startup (including hybrid) or branch run the initial datafile is
    ! set in namelist input 'initfile'
    ! For a continue run - the initfile vluae is created from the time(1:2)
    ! array set below

    if ( iaproc .eq. napout ) write (ndso,950)
    if ( iaproc .eq. napout ) write (ndso,951) 'wave model ...'
    call shr_sys_flush(ndso)

    ! Read namelist (set initfile in w3cesmmd)
    if ( iaproc .eq. napout ) then
       write(ndso,*) 'Read in ww3_inparm namelist from wav_in'//trim(inst_suffix)
       open(newunit=unitn, file='wav_in'//trim(inst_suffix), status='old')
       call shr_nl_find_group_name(unitn, 'ww3_inparm', status=ierr)
       if (ierr == 0) then
          read (unitn, ww3_inparm, iostat=ierr)
          if (ierr /= 0) then
             call shr_sys_abort('problem reading ww3_inparm namelist')
          end if
       end if
       close( unitn )
    end if
    call shr_mpi_bcast(initfile, mpi_comm)
    call shr_mpi_bcast(outfreq, mpi_comm)

    ! Set casename (in w3cesmmd)
    call NUOPC_CompAttributeGet(gcomp, name='case_name', value=cvalue, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    read(cvalue,*) casename

    ! Read in input data and initialize the model
    ! w3init calls w3iors which:
    ! - reads either the initfile if the run is startup or branch
    ! - constructs the filename from the casename variable and the time(:) array
    !   which is set above

    npts = 0
    allocate ( x(1), y(1), pnames(1) )
    pnames(1) = ' '

    call w3init ( 1, 'ww3', nds, ntrace, odat, flgrd, npts, x, y, pnames, iprt, prtfrm, mpi_comm )
    call shr_sys_flush(ndso)

    ! overwrite dt values with variables from coupler
    ! is this a problem with any things being set in w3init?
    dtmax  = real(dtime_sync)
    dtcfl  = real(dtime_sync) / 2. !checked by adrean
    dtcfli = real(dtime_sync)      !checked by adrean
    dtmin  = real(dtime_sync) / 12 !checked by adrean

    call mpi_barrier ( mpi_comm, ierr )

    !--------------------------------------------------------------------
    ! Mesh initialization
    !--------------------------------------------------------------------

    ! Note that nsea is the global number of sea points - and nseal is
    ! the local number of sea points

    !-------------
    ! create a  global index array for sea points
    !-------------

    allocate(gindex_sea(nseal))
    do jsea=1, nseal
       isea = iaproc + (jsea-1)*naproc
       ix = mapsf(isea,1)
       iy = mapsf(isea,2)
       gindex_sea(jsea) = ix + (iy-1)*nx
    end do

    !-------------
    ! create a global index array for non-sea (i.e. land points)
    !-------------

    allocate(mask_global(nx*ny), mask_local(nx*ny))
    mask_local(:) = 0
    do jsea=1, nseal
       isea = iaproc + (jsea-1)*naproc
       ix = mapsf(isea,1)
       iy = mapsf(isea,2)
       mask_local(ix + (iy-1)*nx) = 1
    end do
    call ESMF_VMAllReduce(vm, sendData=mask_local, recvData=mask_global, count=nx*ny, &
         reduceflag=ESMF_REDUCE_MAX, rc=rc)

    nlnd_global = nx*ny - nsea
    nlnd_local = nlnd_global / naproc
    my_lnd_start = nlnd_local*iam + min(iam, mod(nlnd_global, naproc)) + 1
    if (iam < mod(nlnd_global, naproc)) then
       nlnd_local = nlnd_local + 1
    end if
    my_lnd_end = my_lnd_start + nlnd_local - 1

    allocate(gindex_lnd(my_lnd_end - my_lnd_start + 1))
    ncnt = 0
    do n = 1,nx*ny
       if (mask_global(n) == 0) then ! this is a land pont
          ncnt = ncnt + 1
          if (ncnt >= my_lnd_start .and. ncnt <= my_lnd_end) then
             gindex_lnd(ncnt - my_lnd_start + 1) = n
          end if
       end if
    end do

    !-------------
    ! create a global index that includes both sea and land - but put land at the end
    !-------------
    nlnd = (my_lnd_end - my_lnd_start + 1) 
    allocate(gindex(nlnd + nseal))
    do ncnt = 1,nlnd + nseal
       if (ncnt <= nseal) then
          gindex(ncnt) = gindex_sea(ncnt)
       else
          gindex(ncnt) = gindex_lnd(ncnt-nseal)
       end if
    end do

    !-------------
    ! create distGrid from global index array
    !-------------
    DistGrid = ESMF_DistGridCreate(arbSeqIndexList=gindex, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    !-------------
    ! create the mesh 
    !-------------
    call NUOPC_CompAttributeGet(gcomp, name='mesh_wav', value=cvalue, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! read in the mesh with an auto-generated distGrid
    EMeshTemp = ESMF_MeshCreate(filename=trim(cvalue), fileformat=ESMF_FILEFORMAT_ESMFMESH, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    if (masterproc) then
       write(stdout,*)'mesh file for domain is ',trim(cvalue)
    end if

    ! recreate the mesh using the above distGrid
    EMesh = ESMF_MeshCreate(EMeshTemp, elementDistgrid=Distgrid, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    !--------------------------------------------------------------------
    ! Realize the actively coupled fields
    !--------------------------------------------------------------------

    call realize_fields(gcomp, mesh=Emesh, flds_scalar_name=flds_scalar_name, flds_scalar_num=flds_scalar_num, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    !--------------------------------------------------------------------
    ! end redirection of share output to wav log
    !--------------------------------------------------------------------
    call shr_sys_flush(ndso)
    call shr_file_setlogunit (shrlogunit)


900 FORMAT (/15X,'      *** WAVEWATCH III Program shell ***      '/ &
         15X,'==============================================='/)
901 FORMAT ( '  Comment character is ''',A,''''/)

930 FORMAT (/'  Time interval : '/                                  &
         ' --------------------------------------------------')
931 FORMAT ( '       Starting time : ',A)

940 FORMAT (/'  Output requests : '/                                &
         ' --------------------------------------------------'/ &
         '       ',A)

950 FORMAT (/'  Initializations :'/                                 &
         ' --------------------------------------------------')
951 FORMAT ( '       ',A)
1945 FORMAT ( '            Fields   : ',A)
1946 FORMAT ( '                       ',A)

  end subroutine InitializeRealize

  !===============================================================================

  subroutine DataInitialize(gcomp, rc)

    ! input/output variables
    type(ESMF_GridComp)  :: gcomp
    integer, intent(out) :: rc

    ! local variables
    type(ESMF_State)  :: exportState
    integer           :: jsea
    real(r8), pointer :: sw_lamult(:)
    real(r8), pointer :: sw_ustokes(:)
    real(r8), pointer :: sw_vstokes(:)
    character(len=*),parameter :: subname = '(wav_comp_nuopc:DataInitialize)'
    ! -------------------------------------------------------------------

    !--------------------------------------------------------------------
    ! Create export state 
    !--------------------------------------------------------------------
    call NUOPC_ModelGet(gcomp, exportState=exportState, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call state_getfldptr(exportState, 'Sw_lamult', sw_lamult, rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call state_getfldptr(exportState, 'Sw_ustokes', sw_ustokes, rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call state_getfldptr(exportState, 'Sw_vstokes', sw_vstokes, rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    do jsea=1, nseal
       sw_lamult(jsea)  = 1.
       sw_ustokes(jsea) = 0.
       sw_vstokes(jsea) = 0.
    enddo

    ! Set global grid size scalars in export state
    call State_SetScalar(dble(NX), flds_scalar_index_nx, exportState, &
         flds_scalar_name, flds_scalar_num, rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call State_SetScalar(dble(NY), flds_scalar_index_ny, exportState, &
         flds_scalar_name, flds_scalar_num, rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

  end subroutine DataInitialize

  !=====================================================================

  subroutine ModelAdvance(gcomp, rc)
  
    !------------------------
    ! Run WW3
    !------------------------

    ! arguments:
    type(ESMF_GridComp)  :: gcomp
    integer, intent(out) :: rc

    ! local variables
    type(ESMF_State) :: importState
    type(ESMF_State) :: exportState
    type(ESMF_Clock) :: clock
    type(ESMF_Alarm) :: alarm
    type(ESMF_TIME)  :: ETime
    integer          :: yy,mm,dd,hh,ss
    integer          :: ymd        ! current year-month-day
    integer          :: tod        ! current time of day (sec)
    integer          :: time0(2)
    integer          :: timen(2)
    integer          :: shrlogunit ! original log unit and level
    character(len=*),parameter :: subname = '(wav_comp_nuopc:ModelAdvance) ' 
    !-------------------------------------------------------

    !------------
    ! Reset shr logging to my log file
    !------------
    call shr_file_getLogUnit (shrlogunit)
    call shr_file_setLogUnit (stdout)

    !------------
    ! query the Component for its importState, exportState and clock
    !------------

    call ESMF_GridCompGet(gcomp, importState=importState, exportState=exportState, clock=clock, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    !------------
    ! Determine if time to write restart
    !------------

    call ESMF_ClockGetAlarm(clock, alarmname='alarm_restart', alarm=alarm, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    if (ESMF_AlarmIsRinging(alarm, rc=rc)) then
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       rstwr = .true.
       call ESMF_AlarmRingerOff( alarm, rc=rc )
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
    else
       rstwr = .false.
    endif

    !------------
    ! Determine if time to write restart
    !------------
    if (outfreq .gt. 0 .and. mod(hh, outfreq) .eq. 0 ) then
       ! output every outfreq hours
       histwr = .true.
    else
       call ESMF_ClockGetAlarm(clock, alarmname='alarm_history', alarm=alarm, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       if (ESMF_AlarmIsRinging(alarm, rc=rc)) then
          if (ChkErr(rc,__LINE__,u_FILE_u)) return
          histwr = .true.
          call ESMF_AlarmRingerOff( alarm, rc=rc )
          if (ChkErr(rc,__LINE__,u_FILE_u)) return
       else
          histwr = .false.
       endif
    end if

    !------------
    ! Determine time info 
    !------------

    ! use current time for next time step the NUOPC clock is not updated
    ! until the end of the time interval
    call ESMF_ClockGetNextTime(clock, nextTime=ETime, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call ESMF_TimeGet( ETime, yy=yy, mm=mm, dd=dd, s=tod, rc=rc )
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call shr_cal_ymd2date(yy, mm, dd, ymd)
    hh = tod/3600
    mm = (tod - (hh * 3600))/60
    ss = tod - (hh*3600) - (mm*60)
    timen(1) = ymd
    timen(2) = hh*10000 + mm*100 + ss

    call ESMF_ClockGetNextTime(clock, nextTime=ETime, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call ESMF_TimeGet( ETime, yy=yy, mm=mm, dd=dd, s=tod, rc=rc )
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call shr_cal_ymd2date(yy, mm, dd, ymd)
    hh = tod/3600
    mm = (tod - (hh * 3600))/60
    ss = tod - (hh*3600) - (mm*60)
    time0(1) = ymd
    time0(2) = hh*10000 + mm*100 + ss

    time = time0

    !------------
    ! Obtain import data from import state
    !------------
    call import_fields(gcomp, rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    !------------
    ! Run the wave model for the given interval
    !------------
    call w3wave ( 1, timen )

    !------------
    ! Create export state
    !------------

    call export_fields(gcomp, rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    !------------
    ! Reset shr logging to original values
    !------------

    call shr_file_setLogUnit (shrlogunit)
    call shr_sys_flush(stdout)

  end subroutine ModelAdvance

  !===============================================================================

  subroutine ModelSetRunClock(gcomp, rc)

    ! input/output variables
    type(ESMF_GridComp)  :: gcomp
    integer, intent(out) :: rc

    ! local variables
    type(ESMF_Clock)         :: mclock, dclock
    type(ESMF_Time)          :: mcurrtime, dcurrtime
    type(ESMF_Time)          :: mstoptime
    type(ESMF_TimeInterval)  :: mtimestep, dtimestep
    character(len=256)       :: cvalue
    character(len=256)       :: restart_option ! Restart option units
    integer                  :: restart_n      ! Number until restart interval
    integer                  :: restart_ymd    ! Restart date (YYYYMMDD)
    type(ESMF_ALARM)         :: restart_alarm
    character(len=256)       :: stop_option    ! Stop option units
    integer                  :: stop_n         ! Number until stop interval
    integer                  :: stop_ymd       ! Stop date (YYYYMMDD)
    type(ESMF_ALARM)         :: stop_alarm
    character(len=256)       :: history_option ! History option units
    integer                  :: history_n      ! Number until history interval
    integer                  :: history_ymd    ! History date (YYYYMMDD)
    type(ESMF_ALARM)         :: history_alarm
    character(len=128)       :: name
    integer                  :: alarmcount
    character(len=*),parameter :: subname=trim(modName)//':(ModelSetRunClock) '
    !-------------------------------------------------------------------------------

    rc = ESMF_SUCCESS
    call ESMF_LogWrite(subname//' called', ESMF_LOGMSG_INFO)

    ! query the Component for its clocks
    call NUOPC_ModelGet(gcomp, driverClock=dclock, modelClock=mclock, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call ESMF_ClockGet(dclock, currTime=dcurrtime, timeStep=dtimestep, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call ESMF_ClockGet(mclock, currTime=mcurrtime, timeStep=mtimestep, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    !--------------------------------
    ! force model clock currtime and timestep to match driver and set stoptime
    !--------------------------------

    mstoptime = mcurrtime + dtimestep
    call ESMF_ClockSet(mclock, currTime=dcurrtime, timeStep=dtimestep, stopTime=mstoptime, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    !--------------------------------
    ! set restart, stop and history alarms
    !--------------------------------

    call ESMF_ClockGetAlarmList(mclock, alarmlistflag=ESMF_ALARMLIST_ALL, alarmCount=alarmCount, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    if (alarmCount == 0) then

       call ESMF_GridCompGet(gcomp, name=name, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       call ESMF_LogWrite(subname//'setting alarms for' // trim(name), ESMF_LOGMSG_INFO)

       !----------------
       ! Restart alarm
       !----------------
       call NUOPC_CompAttributeGet(gcomp, name="restart_option", value=restart_option, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return

       call NUOPC_CompAttributeGet(gcomp, name="restart_n", value=cvalue, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       read(cvalue,*) restart_n

       call NUOPC_CompAttributeGet(gcomp, name="restart_ymd", value=cvalue, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       read(cvalue,*) restart_ymd

       call alarmInit(mclock, restart_alarm, restart_option, &
            opt_n   = restart_n,           &
            opt_ymd = restart_ymd,         &
            RefTime = mcurrTime,           &
            alarmname = 'alarm_restart', rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return

       call ESMF_AlarmSet(restart_alarm, clock=mclock, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return

       !----------------
       ! Stop alarm
       !----------------
       call NUOPC_CompAttributeGet(gcomp, name="stop_option", value=stop_option, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return

       call NUOPC_CompAttributeGet(gcomp, name="stop_n", value=cvalue, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       read(cvalue,*) stop_n

       call NUOPC_CompAttributeGet(gcomp, name="stop_ymd", value=cvalue, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       read(cvalue,*) stop_ymd

       call alarmInit(mclock, stop_alarm, stop_option, &
            opt_n   = stop_n,           &
            opt_ymd = stop_ymd,         &
            RefTime = mcurrTime,           &
            alarmname = 'alarm_stop', rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return

       call ESMF_AlarmSet(stop_alarm, clock=mclock, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return

    end if

    !----------------
    ! History alarm
    !----------------
    call NUOPC_CompAttributeGet(gcomp, name="history_option", value=history_option, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call NUOPC_CompAttributeGet(gcomp, name="history_n", value=cvalue, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    read(cvalue,*) history_n

    call NUOPC_CompAttributeGet(gcomp, name="history_ymd", value=cvalue, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    read(cvalue,*) history_ymd

    call alarmInit(mclock, history_alarm, history_option, &
         opt_n   = history_n,           &
         opt_ymd = history_ymd,         &
         RefTime = mcurrTime,           &
         alarmname = 'alarm_history', rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call ESMF_AlarmSet(history_alarm, clock=mclock, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    !--------------------------------
    ! Advance model clock to trigger alarms then reset model clock back to currtime
    !--------------------------------

    call ESMF_ClockAdvance(mclock,rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call ESMF_ClockSet(mclock, currTime=dcurrtime, timeStep=dtimestep, stopTime=mstoptime, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call ESMF_LogWrite(subname//' done', ESMF_LOGMSG_INFO)

  end subroutine ModelSetRunClock

  !===============================================================================

  subroutine ModelFinalize(gcomp, rc)

    ! input/output variables
    type(ESMF_GridComp)  :: gcomp
    integer, intent(out) :: rc

    ! local variables
    character(*), parameter :: F00   = "('(ww3_comp_nuopc) ',8a)"
    character(*), parameter :: F91   = "('(ww3_comp_nuopc) ',73('-'))"
    character(len=*),parameter  :: subname=trim(modName)//':(ModelFinalize) '
    !-------------------------------------------------------------------------------

    !--------------------------------
    ! Finalize routine
    !--------------------------------

    rc = ESMF_SUCCESS
    call ESMF_LogWrite(subname//' called', ESMF_LOGMSG_INFO)

    if (masterproc) then
       write(stdout,F91)
       write(stdout,F00) 'WW3: end of main integration loop'
       write(stdout,F91)
    end if

    call ESMF_LogWrite(subname//' done', ESMF_LOGMSG_INFO)

  end subroutine ModelFinalize

  !===============================================================================

end module wav_comp_nuopc
