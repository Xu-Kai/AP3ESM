module wav_import_export

  use ESMF
  use NUOPC
  use NUOPC_Model
  use shr_kind_mod    , only : r8 => shr_kind_r8
  use shr_cal_mod     , only : shr_cal_ymd2date
  use wav_shr_methods , only : chkerr

  implicit none
  private ! except

  public  :: advertise_fields
  public  :: realize_fields
  public  :: import_fields
  public  :: export_fields
  public  :: state_getfldptr

  private :: fldlist_add
  private :: fldlist_realize

  type fld_list_type
     character(len=128) :: stdname
  end type fld_list_type

  integer, parameter     :: fldsMax = 100
  integer                :: fldsToWav_num = 0
  integer                :: fldsFrWav_num = 0
  type (fld_list_type)   :: fldsToWav(fldsMax)
  type (fld_list_type)   :: fldsFrWav(fldsMax)

  integer     ,parameter :: dbug_flag = 0 ! internal debug level
  character(*),parameter :: u_FILE_u = &
       __FILE__

!===============================================================================
contains
!===============================================================================

  subroutine advertise_fields(importState, ExportState, flds_scalar_name, rc)

    ! input/output variables
    type(ESMF_State)               :: importState
    type(ESMF_State)               :: exportState
    character(len=*) , intent(in)  :: flds_scalar_name
    integer          , intent(out) :: rc

    ! local variables
    integer          :: n, num
    character(len=*), parameter :: subname='(wav_import_export:advertise_fields)'
    !-------------------------------------------------------------------------------

    rc = ESMF_SUCCESS

    !--------------------------------
    ! Advertise import fields
    !--------------------------------

    call fldlist_add(fldsToWav_num, fldsToWav, trim(flds_scalar_name))

    call fldlist_add(fldsToWav_num, fldsToWav, 'Sa_u'       )
    call fldlist_add(fldsToWav_num, fldsToWav, 'Sa_v'       )
    call fldlist_add(fldsToWav_num, fldsToWav, 'Sa_tbot'    )
    call fldlist_add(fldsToWav_num, fldsToWav, 'Si_ifrac'   )
    call fldlist_add(fldsToWav_num, fldsToWav, 'So_t'       )
    call fldlist_add(fldsToWav_num, fldsToWav, 'So_u'       )
    call fldlist_add(fldsToWav_num, fldsToWav, 'So_v'       )
    call fldlist_add(fldsToWav_num, fldsToWav, 'So_bldepth' )

    do n = 1,fldsToWav_num
       call NUOPC_Advertise(importState, standardName=fldsToWav(n)%stdname, &
            TransferOfferGeomObject='will provide', rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
    end do

    !--------------------------------
    ! Advertise export fields
    !--------------------------------

    call fldlist_add(fldsFrWav_num, fldsFrWav, trim(flds_scalar_name))

    call fldlist_add(fldsFrWav_num, fldsFrWav, 'Sw_lamult' )
    call fldlist_add(fldsFrWav_num, fldsFrWav, 'Sw_ustokes')
    call fldlist_add(fldsFrWav_num, fldsFrWav, 'Sw_vstokes')
   !call fldlist_add(fldsFrWav_num, fldsFrWav, 'Sw_hstokes')

    do n = 1,fldsFrWav_num
       call NUOPC_Advertise(exportState, standardName=fldsFrWav(n)%stdname, &
            TransferOfferGeomObject='will provide', rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
    end do

  end subroutine advertise_fields

!===============================================================================

  subroutine realize_fields(gcomp, mesh, flds_scalar_name, flds_scalar_num, rc)

    ! input/output variables
    type(ESMF_GridComp)            :: gcomp
    type(ESMF_Mesh)                :: mesh
    character(len=*) , intent(in)  :: flds_scalar_name
    integer          , intent(in)  :: flds_scalar_num
    integer          , intent(out) :: rc

    ! local variables
    type(ESMF_State)     :: importState
    type(ESMF_State)     :: exportState
    character(len=*), parameter :: subname='(wav_import_export:realize_fields)'
    !---------------------------------------------------------------------------

    rc = ESMF_SUCCESS

    call NUOPC_ModelGet(gcomp, importState=importState, exportState=exportState, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call fldlist_realize( &
         state=ExportState, &
         fldList=fldsFrWav, &
         numflds=fldsFrWav_num, &
         flds_scalar_name=flds_scalar_name, &
         flds_scalar_num=flds_scalar_num, &
         tag=subname//':WW3Export',&
         mesh=mesh, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call fldlist_realize( &
         state=importState, &
         fldList=fldsToWav, &
         numflds=fldsToWav_num, &
         flds_scalar_name=flds_scalar_name, &
         flds_scalar_num=flds_scalar_num, &
         tag=subname//':WW3Import',&
         mesh=mesh, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

  end subroutine realize_fields

!===============================================================================

  subroutine import_fields( gcomp, rc )

    !---------------------------------------------------------------------------
    ! Obtain the wave input from the mediator
    !---------------------------------------------------------------------------

    use w3gdatmd    , only: nseal, MAPSTA, MAPFS, MAPSF, NX, NY
    use w3idatmd    , only: CX0, CY0, CXN, CYN, DT0, DTN, ICEI, HML, WLEV, flags
    use w3idatmd    , only: TC0, TCN, TLN, TIN, TW0, TWN, WX0, WY0, WXN, WYN
    use w3odatmd    , only: naproc, iaproc, napout
    use w3wdatmd    , only: time
    
    ! input/output variables
    type(ESMF_GridComp)  :: gcomp
    integer, intent(out) :: rc

    ! Local variables
    type(ESMF_State)  :: importState
    type(ESMF_VM)     :: vm
    type(ESMF_Clock)  :: clock 
    type(ESMF_Time)   :: ETime
    type(ESMF_Field)  :: lfield
    integer           :: ymd, tod
    integer           :: yy, mm, dd, hh, ss
    real(r8)          :: so_u_global(nx*ny)
    real(r8)          :: so_v_global(nx*ny)
    real(r8)          :: so_t_global(nx*ny)
    real(r8)          :: so_bldepth_global(nx*ny)
    real(r8)          :: sa_u_global(nx*ny)
    real(r8)          :: sa_v_global(nx*ny)
    real(r8)          :: sa_tbot_global(nx*ny)
    real(r8)          :: si_ifrac_global(nx*ny)
    real(r8)          :: temp_global(nx*ny)
    real(r8), pointer :: so_u(:)
    real(r8), pointer :: so_v(:)
    real(r8), pointer :: so_t(:)
    real(r8), pointer :: so_bldepth(:)
    real(r8), pointer :: sa_u(:)
    real(r8), pointer :: sa_v(:)
    real(r8), pointer :: sa_tbot(:)
    real(r8), pointer :: si_ifrac(:)
    integer           :: n, ix, iy, isea
    real(r8)          :: def_value
    integer           :: time0(2) ! starting time of the run interval
    integer           :: timen(2) ! ending time of the run interval
    character(len=*), parameter :: subname='(wav_import_export:import_fields)'
    !---------------------------------------------------------------------------

    rc = ESMF_SUCCESS
    call ESMF_LogWrite(subname//' called', ESMF_LOGMSG_INFO)

    ! Get import state, clock and vm

    call ESMF_GridCompGet(gcomp, clock=clock, importState=importState, vm=vm, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! Determine time0 and timen - note have not advanced the model clock yet with nuopc

    call ESMF_ClockGetNextTime( clock, Etime, rc=rc )
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call ESMF_TimeGet( ETime, yy=yy, mm=mm, dd=dd, s=tod, rc=rc )
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call shr_cal_ymd2date(yy, mm, dd, ymd)
    hh = tod/3600
    mm = (tod - (hh * 3600))/60
    ss = tod - (hh*3600) - (mm*60)
    timen(1) = ymd
    timen(2) = hh*10000 + mm*100 + ss

    call ESMF_ClockGet( clock, currTime=ETime, rc=rc )
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

    ! Determine global data 

    call state_getfldptr(importState, 'So_u', so_u, rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call state_getfldptr(importState, 'So_v', so_v, rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call state_getfldptr(importState, 'So_t', so_t, rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call state_getfldptr(importState, 'So_bldepth', so_bldepth, rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call state_getfldptr(importState, 'Sa_u', sa_u, rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call state_getfldptr(importState, 'Sa_v', sa_v, rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call state_getfldptr(importState, 'Sa_tbot', sa_tbot, rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call state_getfldptr(importState, 'Si_ifrac', si_ifrac, rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    so_u_global(:)       = 0._r8
    so_v_global(:)       = 0._r8
    so_t_global(:)       = 0._r8
    so_bldepth_global(:) = 0._r8
    sa_u_global          = 0._r8
    sa_v_global          = 0._r8
    sa_tbot_global       = 0._r8
    si_ifrac_global      = 0._r8

    do n = 1, nseal
       isea = iaproc + (n-1)*naproc
       ix = mapsf(isea,1)
       iy = mapsf(isea,2)

       so_u_global       (ix + (iy-1)*nx) = so_u(n)
       so_v_global       (ix + (iy-1)*nx) = so_v(n)
       so_t_global       (ix + (iy-1)*nx) = so_t(n)
       so_bldepth_global (ix + (iy-1)*nx) = so_bldepth(n)
       sa_u_global       (ix + (iy-1)*nx) = sa_u(n)
       sa_v_global       (ix + (iy-1)*nx) = sa_v(n)
       sa_tbot_global    (ix + (iy-1)*nx) = sa_tbot(n)
       si_ifrac_global   (ix + (iy-1)*nx) = si_ifrac(n)
    end do

    call ESMF_VMAllReduce(vm, sendData=so_u_global, recvData=temp_global, count=nx*ny, reduceflag=ESMF_REDUCE_SUM, rc=rc)
    so_u_global(:) = temp_global(:)

    call ESMF_VMAllReduce(vm, sendData=so_v_global, recvData=temp_global, count=nx*ny, reduceflag=ESMF_REDUCE_SUM, rc=rc)
    so_v_global(:) = temp_global(:)

    call ESMF_VMAllReduce(vm, sendData=so_t_global, recvData=temp_global, count=nx*ny, reduceflag=ESMF_REDUCE_SUM, rc=rc)
    so_t_global(:) = temp_global(:)

    call ESMF_VMAllReduce(vm, sendData=so_bldepth_global, recvData=temp_global, count=nx*ny, reduceflag=ESMF_REDUCE_SUM, rc=rc)
    so_bldepth_global(:) = temp_global(:)

    call ESMF_VMAllReduce(vm, sendData=sa_u_global, recvData=temp_global, count=nx*ny, reduceflag=ESMF_REDUCE_SUM, rc=rc)
    sa_u_global(:) = temp_global(:)

    call ESMF_VMAllReduce(vm, sendData=sa_v_global, recvData=temp_global, count=nx*ny, reduceflag=ESMF_REDUCE_SUM, rc=rc)
    sa_v_global(:) = temp_global(:)

    call ESMF_VMAllReduce(vm, sendData=sa_tbot_global, recvData=temp_global, count=nx*ny, reduceflag=ESMF_REDUCE_SUM, rc=rc)
    sa_tbot_global(:) = temp_global(:)

    call ESMF_VMAllReduce(vm, sendData=si_ifrac_global, recvData=temp_global, count=nx*ny, reduceflag=ESMF_REDUCE_SUM, rc=rc)
    si_ifrac_global(:) = temp_global(:)

    ! input fields associated with W3FLDG calls in ww3_shel.ftn
    ! these arrays are global, just fill the local cells for use later
    ! fill both the lower (0) and upper (N) bound data with the same values
    ! fill with special values as default, these should not be used in practice
    ! set time for input data to time0 and timen (shouldn't matter)

    def_value = 0.0
    if (flags(1)) then
       TLN  = timen
       WLEV = def_value   ! ssh
    endif
    if (flags(2)) then
       TC0  = time0       ! times for ocn current fields
       TCN  = timen       
       CX0  = def_value   ! ocn u current
       CXN  = def_value
       CY0  = def_value   ! ocn v current
       CYN  = def_value
    endif
    if (flags(3)) then
       TW0  = time0       ! times for atm wind/temp fields.
       TWN  = timen
       WX0  = def_value   ! atm u wind
       WXN  = def_value
       WY0  = def_value   ! atm v wind
       WYN  = def_value
       DT0  = def_value   ! air temp - ocn temp
       DTN  = def_value
    endif
    if (flags(4)) then
       TIN  = timen       ! time for ice field
       ICEI = def_value   ! ice frac
    endif

    ! use these loops for global copy
    n = 0
    do iy = 1,NY
       do ix = 1,NX
          n = n + 1
          if (flags(1)) then
             WLEV(ix,iy) = 0.0
          endif
          if (flags(2)) then
             CX0(ix,iy)  = so_u_global(n) ! ocn u current
             CXN(ix,iy)  = so_u_global(n)
             CY0(ix,iy)  = so_v_global(n) ! ocn v current
             CYN(ix,iy)  = so_v_global(n)
          endif
          if (flags(3)) then
             WX0(ix,iy)  = sa_u_global(n)  ! atm u wind
             WXN(ix,iy)  = sa_u_global(n)
             WY0(ix,iy)  = sa_v_global(n)  ! atm v wind
             WYN(ix,iy)  = sa_v_global(n)
             DT0(ix,iy)  = sa_tbot_global(n) - so_t_global(n)  ! air temp - ocn temp
             DTN(ix,iy)  = sa_tbot_global(n) - so_t_global(n)
          endif
          if (flags(4)) then
             ICEI(ix,iy) = si_ifrac_global(n) ! ice frac
          endif
          if (flags(5)) then
             HML(ix,iy) = max(so_bldepth_global(n), 5.) ! ocn mixing layer depth
          endif
       enddo
    enddo
100 format(a,i6,2x,d21.14)

    call ESMF_LogWrite(trim(subname)//": done", ESMF_LOGMSG_INFO)

  end subroutine import_fields

!====================================================================================

  subroutine export_fields (gcomp, rc)

    !---------------------------------------------------------------------------
    ! Create the export state
    !---------------------------------------------------------------------------

    use shr_const_mod , only : fillvalue=>SHR_CONST_SPVAL
    use w3adatmd      , only : LAMULT, USSX, USSY
    use w3odatmd      , only : naproc, iaproc
    use w3gdatmd      , only : nseal, MAPSTA, MAPFS, MAPSF

    ! input/output/variables
    type(ESMF_GridComp)  :: gcomp
    integer, intent(out) :: rc

    ! Local variables
    type(ESMF_State)  :: exportState
    integer           :: n, jsea, isea, ix, iy, lsize
    real(r8), pointer :: sw_lamult(:)
    real(r8), pointer :: sw_ustokes(:)
    real(r8), pointer :: sw_vstokes(:)
    character(len=*), parameter :: subname='(wav_import_export:export_fields)'
    !---------------------------------------------------------------------------

    rc = ESMF_SUCCESS
    call ESMF_LogWrite(subname//' called', ESMF_LOGMSG_INFO)

    ! Get export state
    call NUOPC_ModelGet(gcomp, exportState=exportState, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! copy ww3 data to coupling datatype
    ! copy enhancement factor, uStokes, vStokes to coupler

    call state_getfldptr(exportState, 'Sw_lamult', sw_lamult, rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call state_getfldptr(exportState, 'Sw_ustokes', sw_ustokes, rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call state_getfldptr(exportState, 'Sw_vstokes', sw_vstokes, rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    do jsea=1, nseal
       isea = iaproc + (jsea-1)*naproc
       ix  = MAPSF(ISEA,1)
       iy  = MAPSF(ISEA,2)
       if (MAPSTA(iy,ix) .eq. 1) then
          ! QL, 160530, LAMULT now calculated in WW3 (w3iogomd.f90)
          sw_lamult(jsea)  = LAMULT(ISEA)
          sw_ustokes(jsea) = USSX(ISEA)
          sw_vstokes(jsea) = USSY(ISEA)
       else
          sw_lamult(jsea)  = 1.
          sw_ustokes(jsea) = 0.
          sw_vstokes(jsea) = 0.
       endif
       ! sw_hstokes(jsea) = ??
    enddo

    ! Fill in the local land points with fill value
    lsize = size(sw_lamult)
    do n = nseal+1,lsize
       sw_lamult(n)  = fillvalue
       sw_ustokes(n) = fillvalue
       sw_vstokes(n) = fillvalue
    end do

  end subroutine export_fields

  !===============================================================================

  subroutine fldlist_add(num, fldlist, stdname)
    integer,                    intent(inout) :: num
    type(fld_list_type),        intent(inout) :: fldlist(:)
    character(len=*),           intent(in)    :: stdname

    ! local variables
    integer :: rc
    character(len=*), parameter :: subname='(wav_import_export:fldlist_add)'
    !-------------------------------------------------------------------------------

    ! Set up a list of field information

    num = num + 1
    if (num > fldsMax) then
       call ESMF_LogWrite(trim(subname)//": ERROR num > fldsMax "//trim(stdname), &
            ESMF_LOGMSG_ERROR, line=__LINE__, file=__FILE__)
       return
    endif
    fldlist(num)%stdname = trim(stdname)

  end subroutine fldlist_add

  !===============================================================================

  subroutine fldlist_realize(state, fldList, numflds, flds_scalar_name, flds_scalar_num, mesh, tag, rc)

    use NUOPC , only : NUOPC_IsConnected, NUOPC_Realize
    use ESMF  , only : ESMF_MeshLoc_Element, ESMF_FieldCreate, ESMF_TYPEKIND_R8
    use ESMF  , only : ESMF_MAXSTR, ESMF_Field, ESMF_State, ESMF_Mesh, ESMF_StateRemove
    use ESMF  , only : ESMF_LogFoundError, ESMF_LOGMSG_INFO, ESMF_SUCCESS
    use ESMF  , only : ESMF_LogWrite, ESMF_LOGMSG_ERROR, ESMF_LOGERR_PASSTHRU

    ! input/output variables
    type(ESMF_State)    , intent(inout) :: state
    type(fld_list_type) , intent(in)    :: fldList(:)
    integer             , intent(in)    :: numflds
    character(len=*)    , intent(in)    :: flds_scalar_name
    integer             , intent(in)    :: flds_scalar_num
    character(len=*)    , intent(in)    :: tag
    type(ESMF_Mesh)     , intent(in)    :: mesh
    integer             , intent(inout) :: rc

    ! local variables
    integer                :: n
    type(ESMF_Field)       :: field
    character(len=80)      :: stdname
    character(len=*),parameter  :: subname='(wav_import_export:fldlist_realize)'
    ! ----------------------------------------------

    rc = ESMF_SUCCESS

    do n = 1, numflds
       stdname = fldList(n)%stdname
       if (NUOPC_IsConnected(state, fieldName=stdname)) then
          if (stdname == trim(flds_scalar_name)) then
             call ESMF_LogWrite(trim(subname)//trim(tag)//" Field = "//trim(stdname)//" is connected on root pe", &
                  ESMF_LOGMSG_INFO)

             ! Create the scalar field
             call SetScalarField(field, flds_scalar_name, flds_scalar_num, rc=rc)
             if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=u_FILE_u)) return

          else

             call ESMF_LogWrite(trim(subname)//trim(tag)//" Field = "//trim(stdname)//" is connected using mesh", &
                  ESMF_LOGMSG_INFO)
             ! Create the field
             field = ESMF_FieldCreate(mesh, ESMF_TYPEKIND_R8, name=stdname, meshloc=ESMF_MESHLOC_ELEMENT, rc=rc)
             if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=u_FILE_u)) return

          endif

          ! NOW call NUOPC_Realize
          call NUOPC_Realize(state, field=field, rc=rc)
          if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=u_FILE_u)) return

       else

          if (stdname /= trim(flds_scalar_name)) then
             call ESMF_LogWrite(subname // trim(tag) // " Field = "// trim(stdname) // " is not connected.", &
                  ESMF_LOGMSG_INFO)
             call ESMF_StateRemove(state, (/stdname/), rc=rc)
             if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=u_FILE_u)) return
          end if

       end if
    end do

  contains  !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    subroutine SetScalarField(field, flds_scalar_name, flds_scalar_num, rc)
      ! ----------------------------------------------
      ! create a field with scalar data on the root pe
      ! ----------------------------------------------
      use ESMF, only : ESMF_Field, ESMF_DistGrid, ESMF_Grid
      use ESMF, only : ESMF_DistGridCreate, ESMF_GridCreate, ESMF_LogFoundError, ESMF_LOGERR_PASSTHRU
      use ESMF, only : ESMF_FieldCreate, ESMF_GridCreate, ESMF_TYPEKIND_R8

      type(ESMF_Field) , intent(inout) :: field
      character(len=*) , intent(in)    :: flds_scalar_name
      integer          , intent(in)    :: flds_scalar_num
      integer          , intent(inout) :: rc

      ! local variables
      type(ESMF_Distgrid) :: distgrid
      type(ESMF_Grid)     :: grid
      character(len=*), parameter :: subname='(wav_import_export:SetScalarField)'
      ! ----------------------------------------------

      rc = ESMF_SUCCESS

      ! create a DistGrid with a single index space element, which gets mapped onto DE 0.
      distgrid = ESMF_DistGridCreate(minIndex=(/1/), maxIndex=(/1/), rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=u_FILE_u)) return

      grid = ESMF_GridCreate(distgrid, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=u_FILE_u)) return

      field = ESMF_FieldCreate(name=trim(flds_scalar_name), grid=grid, typekind=ESMF_TYPEKIND_R8, &
           ungriddedLBound=(/1/), ungriddedUBound=(/flds_scalar_num/), gridToFieldMap=(/2/), rc=rc) ! num of scalar values
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=u_FILE_u)) return

    end subroutine SetScalarField

  end subroutine fldlist_realize

  !===============================================================================

  subroutine state_getfldptr(State, fldname, fldptr, rc)
    ! ----------------------------------------------
    ! Get pointer to a state field
    ! ----------------------------------------------
    use ESMF , only : ESMF_State, ESMF_Field, ESMF_Mesh, ESMF_FieldStatus_Flag
    use ESMF , only : ESMF_StateGet, ESMF_FieldGet, ESMF_MeshGet
    use ESMF , only : ESMF_FIELDSTATUS_COMPLETE, ESMF_FAILURE

    ! input/output variables
    type(ESMF_State),  intent(in)    :: State
    character(len=*),  intent(in)    :: fldname
    real(R8), pointer, intent(out)   :: fldptr(:)
    integer,           intent(out)   :: rc

    ! local variables
    type(ESMF_FieldStatus_Flag) :: status
    type(ESMF_Field)            :: lfield
    type(ESMF_Mesh)             :: lmesh
    integer                     :: nnodes, nelements
    character(len=*), parameter :: subname='(wav_import_export:state_getfldptr)'
    ! ----------------------------------------------

    rc = ESMF_SUCCESS
    if (dbug_flag > 5) then
       call ESMF_LogWrite(trim(subname)//": called", ESMF_LOGMSG_INFO)
    end if

    call ESMF_StateGet(State, itemName=trim(fldname), field=lfield, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call ESMF_FieldGet(lfield, status=status, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    if (status /= ESMF_FIELDSTATUS_COMPLETE) then
       call ESMF_LogWrite(trim(subname)//": ERROR data not allocated ", ESMF_LOGMSG_INFO, rc=rc)
       rc = ESMF_FAILURE
       return
    else
       call ESMF_FieldGet(lfield, mesh=lmesh, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return

       call ESMF_MeshGet(lmesh, numOwnedNodes=nnodes, numOwnedElements=nelements, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return

       if (nnodes == 0 .and. nelements == 0) then
          call ESMF_LogWrite(trim(subname)//": no local nodes or elements ", ESMF_LOGMSG_INFO)
          rc = ESMF_FAILURE
          return
       end if

       call ESMF_FieldGet(lfield, farrayPtr=fldptr, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
    endif  ! status

    if (dbug_flag > 5) then
       call ESMF_LogWrite(trim(subname)//": done", ESMF_LOGMSG_INFO)
    end if

  end subroutine state_getfldptr

end module wav_import_export
