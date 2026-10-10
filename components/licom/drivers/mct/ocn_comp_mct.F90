!#define LOGMSG()
!write(mytid+600,'(a,i4)')"LICOM",__LINE__
!#define LOGMSGCarbon()
!write(600+mytid,'(a,i4)')"Licom CARBON",__LINE__

module ocn_comp_mct

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   !MODULE:    OCN_COMP_MCT
!   !AUTHOR:    Huimin Li
!   !Date:      2012/5/30
!
!   !DESCRIPTION:
!               This is the main driver for the licom coupled in CPL7
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   use mct_mod
   use ESMF_ClockMod, only : ESMF_Clock
   use seq_flds_mod
   use seq_cdata_mod
   use seq_infodata_mod
   use seq_timemgr_mod
   use seq_comm_mct,      only : seq_comm_suffix, seq_comm_inst, seq_comm_name
   use shr_file_mod
   use shr_cal_mod, only : shr_cal_date2ymd,shr_cal_numDaysInMonth
   use shr_const_mod,only:SHR_CONST_SPVAL !linpf 2012Jul26
   use perf_mod
   use fluxcpl
   use POP_CplIndices
   use shr_dmodel_mod
   use POP_CommMod
   use domain
   use grid
   use blocks
   use gather_scatter
   use distribution
   use perf_mod
   use ocn_import_export, only: ocn_import, ocn_export

#include <def-undef.h>
#ifdef SW_ATHREAD
!#include <athread.h>
use sunway_mod
use sunway_init
#endif
use param_mod
use pconst_mod
use shr_msg_mod
use licom_sys_mod
use control_mod
use constant_mod
use shr_cal_mod,       only: shr_cal_date2ymd,shr_cal_numDaysInMonth
use msg_mod, only: tag_1d,tag_2d,tag_3d,tag_4d,nproc,status,mpi_comm_ocn
use tracer_mod
use pmix_mod
use forc_mod
use canuto_mod !LPF20160729
use output_mod !LPF20170728
use work_mod
!
#ifdef USE_OCN_CARBON
use carbon_mod
use cforce_mod
#endif


  implicit none
#include <netcdf.inc>
  public :: ocn_init_mct
  public :: ocn_run_mct
  public :: ocn_final_mct
  SAVE
  private

  private :: ocn_SetGSMap_mct
  private :: ocn_domain_mct


  type(seq_infodata_type), pointer :: &
     infodata

!==========================================================================
  contains
!==========================================================================


!**************************************************************************
!   !ROUTINE:   ocn_init_mct
!   !AUTHOR:    Huimin Li
!   !Date:      2012/5/30
!
!   !INTERFACE:
  subroutine ocn_init_mct(EClock, cdata_o, x2o_o, o2x_o, NLFilename)
!
!   !DESCRIPTION:
!               This is the initialize routine for licom coupled in CPL7
!   !INPUT/OUTPUT PARAMETERS:
    type(ESMF_Clock)            , intent(inout)    :: EClock
    type(seq_cdata)             , intent(inout) :: cdata_o
    type(mct_aVect)             , intent(inout) :: x2o_o
    type(mct_aVect)             , intent(inout) :: o2x_o
    character(len=*), optional  , intent(in)    :: NLFilename ! Namelist filename

!   !LOCAL VARIABLES
    type(mct_gsMap), pointer :: &
       gsMap_o

    type(mct_gGrid), pointer :: &
       dom_o

    integer (kind(1)) :: &
       nThreads

    integer (kind(1)) :: &
       OCNID,     &
       lsize

    real (r8) ::  &
       precadj

    integer (kind(1)) :: iam,ierr
    character(len=32)  :: starttype          ! infodata start type

    integer :: i_temp, j_temp, i_comp, j_comp ! used in specifying the grid number in each process

    character*80 :: diri,diro,logfile
    namelist /modelio/ diri, diro, logfile
    integer(8) :: long_nx, long_ny, long_km

  !------------------------------------------------------
  ! ROUTINE BEGIN
  !------------------------------------------------------
  mpi_comm_ocn=0

  ! lihuimin 2012.7.16
  ! OCNID, errcode
  call seq_cdata_setptrs(cdata_o, ID=OCNID, mpicom=mpi_comm_ocn, &
       gsMap=gsMap_o, dom=dom_o, infodata=infodata)

  ! five parameters initialized in msg_pass('init') in CPL6 coupled version
  ! get infodata from drv
  cdate = 0
  sec = 0
  ierr = 0
  info_time = 0
  call seq_infodata_GetData( infodata, info_debug=info_dbug)

  ! send initial state to drv
  call seq_infodata_PutData( infodata, ocn_nx=imt_global, ocn_ny=jmt_global)
!  call seq_infodata_PutData( infodata, ocn_prognostic=.true.) !LPF 20120829
  ! lihuimin, 2012.7.25, not use ocnrof_p
  call seq_infodata_PutData( infodata, ocn_prognostic=.true.,ocnrof_prognostic=.true.)
!open ocnrof. !LPF 20120829

  !-----------------------------
  ! namelist & logfile setup
  !-----------------------------
  ! TODO: redirect the log file


  mytid=0

! He - 2010-10-08 | NDay: Number of days since the first day when the simulation begins.
! lihuimin cancle , 2012.6.14
!     NDAY = 0

!
! mpicom_o <-> mpi_comm_ocn
! TODO
      call mpi_comm_rank (mpi_comm_ocn, mytid, ierr)
      call mpi_comm_size (mpi_comm_ocn, nproc, ierr)
!YU

    if (mytid == 0 ) then
        open (15,file="ocn_modelio.nml", form="formatted")
        read (15,modelio)
        close(15)
        open (16,file=trim(logfile), form="formatted",STATUS='REPLACE',ACTION='WRITE')
    end if

    if (mytid == 0 ) then
      write(*,*) "MYTID=",mytid,"Number of Processors is",nproc
      write(*,*) "mpi_comm_ocn:",mpi_comm_ocn
    end if

    call init_ocn_com(mpi_comm_ocn)


      my_task = mytid
      master_task = 0
      POP_mytask = mytid
      POP_mastertask = 0
!YU


   ! lihuimin 2012.6.14, initial indices
   call POP_CplIndicesSet()


!---------------------------------------------------------------------
!     SET THE CONSTANTS USED IN THE MODEL
!---------------------------------------------------------------------
#ifdef COUP
      call shr_sys_flush(6)
#endif
      call grid_init_mem()
      long_nx = nx_block
      long_ny = ny_block
      long_km = km
      call write_init(mpi_comm_ocn,long_nx, long_ny, long_km)
      call pconst_init()
      call dyn_init()
      call work_init()
      call pmix_init()
      call tracer_init()
      CALL CONST
      call output_init()
      call cdf_init()
#ifdef COUP
      call shr_sys_flush(6)
#endif
#ifdef SHOW_TIME
      call run_time('CONST')
#endif

!***********************************************************************
!          SET SOME CONSTANTS FOR THE BOGCM
!***********************************************************************
#ifdef USE_OCN_CARBON
      LOGMSGCarbon()
      CALL CTRLC
      LOGMSGCarbon()
#ifdef SHOW_TIME
      call run_time('CTRLC')
#endif
#endif

!---------------------------------------------------------------------
!     SET MODEL'S RESOLUTION,TOPOGRAPHY AND THE CONSTANT
!     PARAMETERS RELATED TO LATITUDES (J)
!---------------------------------------------------------------------
   call init_domain_blocks

   call init_grid1
#ifdef READ_NCON_SPLIT
      call init_domain_distribution()
#else
      call init_domain_distribution(KMT_G)
#endif
   call init_grid2

   CALL GRIDS

   call calc_coeff

#ifdef SHOW_TIME
      call run_time('GRIDS')
#endif


!---------------------------------------------------------------------
!     SET SURFACE FORCING FIELDS (1: Annual mean; 0: Seasonal cycle)
!---------------------------------------------------------------------
     CALL RDRIVER !LPF20151006

#ifdef SHOW_TIME
      call run_time('RDRIVER')
#endif


!***********************************************************************
!      SET FORCING DATA USED IN CARBON CYCLYE
!***********************************************************************
#ifdef USE_OCN_CARBON
      LOGMSGCarbon()
       CALL CFORCE
      LOGMSGCarbon()
#ifdef SHOW_TIME
       call run_time('CFORCE')
#endif
#endif

!---------------------------------------------------------------------
!     INITIALIZATION
!---------------------------------------------------------------------
      CALL INIRUN

#ifdef SHOW_TIME
      call run_time('INRUN')
#endif
    if (NSTART==1.or.NSTART==2)then !initial run from observation T/S when nstart==1,or from fort22,nstart ==2
    !if (NSTART==1)then !initial run from observation T/S
        number_month=month_ref
        number_day=day_ref
        yy_licom=year_ref
        mm_licom=month_ref
        dd_licom=number_day

        iyfm = year_ref
        month = number_month
      tod_licom=0
      curr_ymd_licom=10000*yy_licom+100*mm_licom+dd_licom

       if (mytid==0) write (*,*) "yy_licom,mm_licom,dd_licom,tod_licom,curr_ymd_licom",yy_licom,mm_licom,dd_licom,tod_licom,curr_ymd_licom
       call seq_timemgr_EClockGetData(EClock,start_ymd=curr_ymd_cpl,start_tod=tod_cpl)
       if (mytid==0) write (*,*)"curr_ymd_cpl",curr_ymd_cpl,tod_cpl

!   IF (mytid.eq.0) THEN
!     if (curr_ymd_licom.ne.curr_ymd_cpl) then
!      if (mytid==0) write (*,*) "initial state mismacth! STOP Now"
!         stop 
!     endif
!   ENDIF
   else !nstart.ne.1.or.nstart.ne.2
        yy_licom = rst_year
        mm_licom = number_month
        dd_licom = number_day
        month = mm_licom  !number_month

      if (mytid==0) write (*,*) "restart mode yy_licom,mm_licom,dd_licom,tod_licom,curr_ymd_licom",yy_licom,mm_licom,dd_licom,tod_licom,curr_ymd_licom
   endif
      totday = 0
      Num_output=NSS*(NMONTH (mm_licom) - dd_licom + 1)
      if(mytid==0) write(*,*)'the days left in current months, nss is Num_output= ',NMONTH (mm_licom) - dd_licom + 1,nss, Num_output
!lhl20130419

!----------------------------------------------------------------------
!     Inialize mct attribute vectors
!     imitate CPL7/POP
!     lihuimin 2012.7.16
!----------------------------------------------------------------------

   call ocn_SetGSMap_mct(mpi_comm_ocn, OCNID, GSMap_o)


   lsize = mct_gsMap_lsize(gsMap_o, mpi_comm_ocn)
!   write(*,*) "mct_gsMap_lsize = ",lsize

   call ocn_domain_mct(lsize, gsMap_o, dom_o)

   call mct_aVect_init(x2o_o, rList=seq_flds_x2o_fields, lsize=lsize)
   call mct_aVect_zero(x2o_o)


   call mct_aVect_init(o2x_o, rList=seq_flds_o2x_fields, lsize=lsize)
   call mct_aVect_zero(o2x_o)



!**********************************************************************
!      INITIALIZATION CARBON MODEL
!**********************************************************************
#ifdef USE_OCN_CARBON
      LOGMSGCarbon()
      CALL INIRUN_PT
      LOGMSGCarbon()
#ifdef SHOW_TIME
      call run_time('INIRUN_PT')
#endif
#endif

!----------------------------------------------------------------------
!      INITIALIZATION Vertical Mixing  
!----------------------------------------------------------------------
#ifdef CANUTO
      call turb_ini
#endif



!---------------------------------------------------------------------
!     INITIALIZATION OF ISOPYCNAL MIXING
!---------------------------------------------------------------------
#ifdef ISO
      CALL ISOPYI
#ifdef SHOW_TIME
      call run_time('ISOPYI')
#endif
#endif

#ifdef COUP
         call shr_sys_flush(6)
#endif
!LPF20130420
        IF (num_cpl-int(num_cpl)>0) THEN
         if(mytid==0) write(6,*)'num_cpl is not integer! STOP Now'
          STOP
        ENDIF
        call seq_timemgr_EClockGetData(EClock, dtime=ocn_cpl_dt)
        licom_cpl_dt = 86400.0/num_cpl
        if (licom_cpl_dt /= ocn_cpl_dt) then
          if(mytid==0) write(*,*)'licom_cpl_dt= ',licom_cpl_dt, &
                     ' ocn_cpl_dt= ',ocn_cpl_dt
           if(mytid == 0) write(*,*)'licom_cpl_dt and &
                       ocn_cpl_dt is not the same! STOP Now'
             stop
        end if

      call flux_cpl !LPF20210802
      call ocn_export(o2x_o%rattr)


      isb=0
      isc=0
      ist=0
      first_fort = 0
#ifdef SW_ATHREAD
       call athread_init()
       CALL init_para !init sw slave info
       call readyt_sw_init
       call readyc_sw_init
       call barotr_sw_init
       call bclinc_sw_init
       call tracer_sw_init
#endif
     if(mytid == 0) write(*,*) 'ocn_init_mct finished'
  end subroutine ocn_init_mct


!**************************************************************************
!! ROUTINE:   ocn_run_mct
!! By Huimin Li, 2012/6/2
!! INTERFACE:
  subroutine ocn_run_mct( EClock, cdata_o, x2o_o, o2x_o)
!   !DESCRIPTION: RUN licom for a coupling interval

!   !INPUT/OUTPUT PARAMETERS:
!lhl20130419
    type(ESMF_Clock)            , intent(inout) :: EClock
!lhl20130419
    type(seq_cdata)             , intent(inout) :: cdata_o
    type(mct_aVect)             , intent(inout) :: x2o_o
    type(mct_aVect)             , intent(inout) :: o2x_o

!lhl20130419
    integer(kind(1))   :: curr_ymd     ! Current date YYYYMMDD
    integer(kind(1))   :: yy,mm,dd     ! year, month, day
    character (len=22) :: fnamet1,fnames1
    character (len=80) :: fnamett,fnamess
    character (len=80) :: datapath
    character (len=18) :: fname
    character (len=80) :: calendar  ! calendar type
!lhl20130419
!lhl20130913
    INTEGER :: IPT1,IPT2
    INTEGER :: IPT14_1,IPT14_2
    INTEGER :: iblock
    INTEGER :: iyfmnew
    REAL(r8):: FACTOR, ek0,t1,t2
!lhl20130913
!----------------------------------------------------------------------
!   get Eclock data information
!----------------------------------------------------------------------

        call seq_timemgr_EClockGetData( EClock, curr_ymd=ymd_sync, &
          curr_tod=tod_sync )
        call shr_cal_date2ymd(ymd_sync,yy,mm,dd)
        iyfm = yy_licom  ! nwmf = iyfm
        mon0 = mm_licom  ! month = (iyfm-1)*12 + mon0
        iday = dd_licom  ! number_day = iday + 1
        if(mytid == 0) write(*,*) 'ymd_sync=,tod_sync=',ymd_sync,tod_sync
        call seq_timemgr_EClockGetData( EClock, calendar=calendar)
        !calendar = 'gregorian'
        dd_cpl = shr_cal_numDaysInMonth(iyfm,mon0,trim(calendar))
    if(mytid == 0) then
        write(*,*) "From CPL7-Time yy=",yy,"mm=",mm,"dd=",dd,"dd_cpl=",dd_cpl,"calendar=",calendar
        write(*,*) "From Ocn-Time iyfm=",iyfm,"mon0=",mon0,"iday=",iday
        write(*,*) "current step=",num_step_per_day,"total=",NSS," num_cpl=",num_cpl 
    endif

!=====================time control for licom output===============
!lhl20120728
       month=(IYFM-1)*12+mon0
       !IMD = NMONTH (MON0)
       IMD = dd_cpl !use cpl time for gregorian or NO_LEAP by mxyu
!lhl20120728
!      month=mm
!      IY0 = (MONTH -1)/12
!      IYFM = IY0+1
!!     IYFM is the number of the current year
!----------------------------------------------------------------------
!   set time 
!----------------------------------------------------------------------
!LPF20170831
           if(first_step.eq.0) num_step_per_day=ocn_cpl_dt/DTS
!LPF20170831
           if(num_step_per_day==NSS) then
             num_step_per_day=0
             dd_licom=dd_licom+1   
           endif
!----------------------------------------------------------------------
!   set number of output 
!----------------------------------------------------------------------
   IF (num_cpl>=1) THEN !0
!LPF20170831
     if (number_month==1) then !1
       if (first_step==0) then !2-0
!           Num_output=NSS*NMONTH (MON0)-NSS/num_cpl
           Num_dailyoutput=NSS-NSS/num_cpl !20170831
          if(mytid==0) write(16,*)'Num_output= ',Num_output,NSS,num_cpl
       endif !2-0
!LPF20170831
     else
!          Num_output=NSS*NMONTH (MON0)
          Num_dailyoutput=NSS !*num_cpl
!        if(mytid==0.and.iday==1) write(*,*)'Num_output= ',Num_output,first_step
     endif !1
   ELSE   
      if(mytid==0) write(16,*)'Error due to num_cpl < 1 ! STOP Now!'
        stop
   ENDIF !0


         if ((mod(num_step_per_day,NSS/num_cpl).eq.0.and.first_step.eq.0)&
           .or.(iday>=1.and.num_step_per_day.eq.0)) then
         if ( simple_assm .eqv. .true.) then
          if(mytid==0)write(*,*) 'into call readtsforres,iday',iday,num_step_per_day
           !fname=./mdatalink//'Mercator_S_'//20140501//'.nc'
            iyfmnew=iyfm+iyfmfnew
            fnamet1(1:9)='Mercator_'
            fnamet1(10:11)='T_'
            fnamet1(20:22)='.nc'
            write(fnamet1(16:17),'(i2.2)')mon0
            write(fnamet1(12:15),'(i4.4)')iyfmnew
            write(fnamet1(18:19),'(i2.2)') iday !dd_licom?
            fnames1(1:9)='Mercator_'
            fnames1(10:11)='S_'
            fnames1(20:22)='.nc'
            write(fnames1(16:17),'(i2.2)')mon0
            write(fnames1(12:15),'(i4.4)')iyfmnew
            write(fnames1(18:19),'(i2.2)') iday !dd_licom?
             datapath='./mdatalink/'
             fnamett=trim(datapath)//fnamet1
             fnamess=trim(datapath)//fnames1
             !fnamett=fnamet1
             !fnamess=fnames1
            if (mytid==0) write(*,*) 'fnamett,fnamess',fnamett,fnamess
           !CALL READTSFORRES(fnamet1,fnames1)
           CALL READTSFORRES(fnamett,fnamess)
         endif
        endif

#ifdef ENABLE_GPLT
        call t_startf("ocn_import")
#endif
         if (mod(num_step_per_day,NSS/num_cpl).eq.0.or.first_step.eq.0) then
           ! if (mytid==0) write(*,*) 'mod(II,NSS/num_cpl)=',mod(num_step_per_day,NSS/num_cpl)
            if (mytid==0) write(16,*)'c2o,num_step_per_day=',num_step_per_day,dd_licom,iday,NSS,num_cpl
            call ocn_import(x2o_o%rAttr)
            call post_cpl
            first_step = 1
         else
!           if (mytid==0) write(*,*) 'mod(II,NSS/num_cpl)=',mod(num_step_per_day,NSS/num_cpl)
         endif

         if(first_fort .eq. 0) then
            if (mytid==0) write(*,*)'SSAVEINSSTEP write start field once'
#if (defined LOWRES)
             call SSAVEINS
#else
             call SSAVEINSSTEP(1)
#endif
             first_fort = 1
         endif
#ifdef ENABLE_GPLT
        call t_stopf("ocn_import")
#endif
!           if (mytid==0) write(*,*) 'NSS/num_cpl,num_step_per_day=',num_step_per_day,NSS/num_cpl,iday
!        if(num_step_per_day.eq.0) then !daily control   
!         if ( simple_assm==.true.) then
!          if(mytid==0)write(*,*) 'into call readtsforres,iday',iday,num_step_per_day
!          CALL READTSFORRES('temp_initial','salt_initial')
!         endif
!        endif
!----------------------------------------------------------------------
!   time interplate phc2 sss and CHLOROPHYLL
!----------------------------------------------------------------------
      IF ( IDAY <= 15) THEN
         IPT1 = MON0-1
         IPT14_1 = IPT1 + 1
         IF (IPT1 == 0 ) then 
            IPT1 = 12
            IPT14_1 = 1
         ENDIF
         IPT2 = MON0
         IPT14_2 = IPT2 + 1
         FACTOR = FLOAT (IDAY -15)/ FLOAT (NMONTH (IPT1)) + 1
      ELSE
         IPT1 = MON0
         IPT2 = MOD (MON0,12) +1
         IPT14_1 = IPT1 + 1
         IPT14_2 = IPT2 + 1
         if(MON0 == 12) IPT14_2 = 14
         FACTOR = FLOAT (IDAY -15)/ FLOAT (NMONTH (IPT1))
      END IF
      if (mytid==0) write(6,*)'MON0,IPT1,IPT2,IPT14_1,IPT14_2,FACTOR:',MON0,IPT1,IPT2,IPT14_1,IPT14_2,FACTOR
!!$OMP PARALLEL DO PRIVATE(J,I)
      DO iblock = 1, nblocks_clinic 
       DO J = 1,JMT
         DO I = 1,IMT
#ifdef SSS_restore
#if (defined LOWRES)
           SSS(I,J,iblock)=(SSS3(I,J,IPT2,iblock)-SSS3(I,J,IPT1,iblock))*FACTOR+SSS3(I,J,IPT1,iblock)
#else
           SSS(I,J,iblock)=(SSS3(I,J,IPT14_2,iblock)-SSS3(I,J,IPT14_1,iblock))*FACTOR+SSS3(I,J,IPT14_1,iblock)
#endif
#endif
           !RESTORE(I,J,1,2,iblock)=SSS(I,J,iblock) !LPF20151103 for restore SSSm to SSSo
#if (defined SOLARCHLORO) 
#if (defined LOWRES)
           chloro(I,J,iblock)= ( (chloro3 (I,J,IPT2,iblock) - chloro3 (I,J,IPT1,iblock))     &
                   * FACTOR + chloro3 (I,J,IPT1,iblock))
#else
           chloro(I,J,iblock)= ( (chloro3 (I,J,IPT14_2,iblock) - chloro3 (I,J,IPT14_1,iblock))     &
                  * FACTOR + chloro3 (I,J,IPT14_1,iblock))
#endif
#endif
#ifdef SST_restore
#if (defined LOWRES)

           SST(I,J,iblock)=(SST3(I,J,IPT2,iblock)-SST3(I,J,IPT1,iblock))*FACTOR+SST3(I,J,IPT1,iblock)
#else
           SST(I,J,iblock)=(SST3(I,J,IPT14_2,iblock)-SST3(I,J,IPT14_1,iblock))*FACTOR+SST3(I,J,IPT14_1,iblock)
                           
#endif
#endif
         END DO
       END DO
      END DO
!!$OMP END PARALLEL DO

#if (defined SOLARCHLORO) 
          CALL  SW_ABSOR
#endif
!LPF20131024


!---------------------------------------------------------------------
!     THERMAL CYCLE
!---------------------------------------------------------------------
#ifdef ENABLE_GPLT
        call t_startf("stepon")
#endif
       DO II = 1,NSS/num_cpl !thermal cycle in a couple interval
!       if (mytid.eq.0) write(16,*)'num_step_per_day=',num_step_per_day,NSS,num_cpl,II
         num_step_per_day=num_step_per_day+1
         
           tod_licom=(num_step_per_day-1)*DTS !+ocn_cpl_dt
           curr_ymd_licom=iyfm*10000+mon0*100+iday   
!          call seq_timemgr_EClockGetData(EClock,curr_ymd=curr_ymd_cpl,curr_tod=tod_cpl)
!           if (mytid==0) write(6,*)'curr_ymd_cpl,tod_cpl',curr_ymd_cpl,tod_cpl
       if (II==1) call energy

         !write fort file
         if(fortenable) then
           if(mod(num_step_per_day-1,nstep) .eq. 0) then
             NDAYS = NDAYS - 1
             if(NDAYS > 0) then
                 if (mytid.eq.0) write(6,*)'write fort interval left ',NDAYS,' days'
             else
                if (mytid.eq.0) write(6,*)'write fort num_step_per_day:',num_step_per_day,nstep,NSS
#if (defined LOWRES)
                 call SSAVEINS
#else
                CALL SSAVEINSSTEP(num_step_per_day)
#endif
                 NDAYS = nstep / NSS
             end if
           end if
         end if

!---------------------------------------------------------------------
!     COMPUTE DENSITY, BAROCLINIC PRESSURE AND THE RELAVANT VARIABLES
#ifdef ENABLE_GPLT
        call t_startf("READYT")
#endif
            CALL READYT
#ifdef ENABLE_GPLT
        call t_stopf("READYT")
#endif
!---------------------------------------------------------------------
!     BAROCLINIC & BAROTROPIC CYCLE
!---------------------------------------------------------------------

            DO JJ = 1,NCC

!     COMPUTE MOMENTUM ADVECTION, DIFFUSION & THEIR VERTICAL INTEGRALS
#ifdef ENABLE_GPLT
        call t_startf("READYC")
#endif
               CALL READYC
#ifdef ENABLE_GPLT
        call t_stopf("READYC")
#endif
!     PREDICTION OF BAROTROPIC MODE
#ifdef ENABLE_GPLT
        call t_startf("BAROTR")
#endif
               CALL BAROTR
#ifdef ENABLE_GPLT
        call t_stopf("BAROTR")
#endif

#ifdef ENABLE_GPLT
        call t_startf("BCLINC")
#endif
               CALL BCLINC
#ifdef ENABLE_GPLT
        call t_stopf("BCLINC")
#endif

            END DO !JJ

!--------------------------------------------------------------------
!     PREDICTION OF PASSIVE TRACER
!--------------------------------------------------------------------

#ifdef ENABLE_GPLT
        call t_startf("TRACER")
#endif
            CALL TRACER
#ifdef ENABLE_GPLT
        call t_stopf("TRACER")
#endif

#ifdef ENABLE_GPLT
        call t_startf("ICESNOW")
#endif
            CALL ICESNOW
#ifdef ENABLE_GPLT
        call t_stopf("ICESNOW")
#endif
          call energy

!--------------------------------------------------------------------
!     PERFORM CONVECTIVE ADJUSTMENT IF UNSTABLE STRATIFICATION OCCURS
!--------------------------------------------------------------------
!        IF (mod(IST,12)==0 ) CALL CONVADJ
#ifdef ENABLE_GPLT
        call t_startf("CONVADJ")
#endif
        CALL CONVADJ
#ifdef ENABLE_GPLT
        call t_stopf("CONVADJ")
#endif

!--------------------------------------------------------------------
!Accumulate the data in a T/S time step interval or in a day for output
#if (defined DAILYACC)
           ! accumulate monthly data
           if (dts_accum .eqv. .true.) CALL ACCUMM
           !accumlate daily data
             if (daily_accum .eqv. .true.) then
              CALL DAILYACCUMM
              Num_outputacc=Num_outputacc+1
             endif
#endif

!--------------------------------------------------------------------
!       Finish thermal cycle loop in a couple interval
   END DO
#ifdef ENABLE_GPLT
        call t_stopf("stepon")
#endif
!--------------------------------------------------------------------

           call energy
!--------------------------------------------------------------------
!          Diagnose if the energy is out of control
!--------------------------------------------------------------------
               iday=dd_licom
!--------------------------------------------------------------------

!--------------------------------------------------------------------
!             Sent ocean data to couple
!--------------------------------------------------------------------
#ifdef ENABLE_GPLT
        call t_startf("ocn_export")
#endif

    IF (mod(num_step_per_day,NSS/num_cpl).eq.0) then
         if (mytid==0) write(16,*)'o2c num_step_per_day=',num_step_per_day
         call flux_cpl
         call ocn_export(o2x_o%rattr)
    ELSE 
      if(mytid==0) write(16,*) 'No o2c num_step_per_day=',num_step_per_day
    ENDIF
#ifdef ENABLE_GPLT
        call t_stopf("ocn_export")
#endif
!--------------------------------------------------------------------
!
! check that internal clock is in sync with master clock
!
!--------------------------------------------------------------------

    IF (mod(num_step_per_day,NSS/num_cpl).eq.0) then
         ymd = yy_licom*10000 + mm_licom*100 + dd_licom 
         !tod = mod(num_step_per_day,NSS)*3600 + licom_cpl_dt*2 !iminute*60 + isecond
         tod = num_step_per_day*DTS !+ ocn_cpl_dt !iminute*60 + isecond
         if(tod>=86400) then
           tod=tod-86400
            if(dd_licom<imd) then
              dd_licom=dd_licom+1
              iday = dd_licom  ! number_day = iday + 1
            else
              mm_licom=mm_licom+1
              dd_licom=1
            endif
             IF (mm_licom>12) THEN 
              mm_licom=1
              yy_licom=yy_licom+1
             ENDIF
         endif
          ymd = yy_licom*10000 + mm_licom*100 + dd_licom
         if(mytid == 0) write(*,*) 'ocn run finish ymd=,imd,tod=',ymd,imd,tod,num_step_per_day
        !tod = ihour*3600 + iminute*60 + isecond
        !if ( .not. seq_timemgr_EClockDateInSync( EClock, ymd, tod ) )then
        call seq_timemgr_EClockGetData( EClock, curr_ymd=ymd_sync, &
          curr_tod=tod_sync )
        if (mytid==0) then
          if (ymd==ymd_sync.and.tod==tod_sync) then
            !  write(6,*)' Internal licom clock in sync with Sync Clock'
          else
           write(16,*)' licom clock',ymd,tod,'drv clock',ymd_sync,tod_sync
           write(16,*)' Internal licom clock not in sync with Sync Clock'
             call licom_shr_abort( 'stop'// &
             ":: Internal licom clock not in sync with Sync Clock")
          endif
        endif
    END IF
!LPF 20130811



!--------------------------------------------------------------------
      if(num_step_per_day==NSS) then 
           num_step_per_day=0
           CALL ADDPS
      end if
!--------------------------------------------------------------------

!     ACCUMULATE SOME VARIABLES FOR MONTHLY OUTPUT
         IF (dd_licom ==1.and.num_step_per_day==0.and.imt_global<3000) then
!       if(mytid==0)    call cpu_time(t1)
!          if(test_input) call SSAVEINPUT
!      if(mytid==0)    call cpu_time(t2)
!       if(mytid==0) write(16,*) 'test time1',t2-t1
!          if(diag_budget) call SSAVEBUDGET
!            if(mytid==0)    call cpu_time(t1)
!       if(mytid==0) write(16,*) 'test time2',t1-t2
       if(mytid==0)    call cpu_time(t1)
#if (defined LOWRES)
           CALL SSAVEMON !turnon onlyfor low
#endif
        if(mytid==0)    call cpu_time(t2)
        if(mytid==0) write(*,*) 'test time3',t2-t1
  
         ENDIF
!------------------------time control for licom output---------------
          iyfm = yy_licom  ! nwmf = iyfm
           mon0 = mm_licom  ! month = (iyfm-1)*12 + mon0
           iday = dd_licom  ! number_day = iday + 1
           month=(IYFM-1)*12+mon0
           !IMD = NMONTH (MON0)
           IMD = dd_cpl !use cpl time for gregorian or NO_LEAP by mxyu

!--------------------------------------------------------------------
!     COMPENSATE THE LOSS OF GROSS MASS
!--------------------------------------------------------------------

#ifdef COUP
         call shr_sys_flush(6)
#endif

!--------------------------------------------------------------------
!     SAVE FOR OUTPUT DAILY AVG
!--------------------------------------------------------------------
!LPF20170726
    if(num_step_per_day==0) then
#ifdef DAILYACC
    CALL DAILYAVG
#ifdef DAILYBUGDET     
    CALL SSAVEBUDGETDAILYAVG
#endif
!     CALL SSAVEDAILYAVG
!     CALL SSAVEMON ! for mon test

#endif

    endif
!LPF20170726

!--------------------------------------------------------------------
!     SAVE FOR OUTPUT DAILY instant value for restart
!--------------------------------------------------------------------
   if(mytid==0)    call cpu_time(t1)

     if ( num_step_per_day== 0) then
      number_day = iday
      if ( number_day  > imd ) then
          number_day = 1
      else if(number_day == imd) then
          if (mytid==0) write(6,*)'write month avg:',iyfm,mon0,number_day
#ifdef DAILYACC
          call SSAVEMONTHAVG
#endif
          Num_output=NSS*NMONTH (MON0)
          if (mytid==0) write(6,*)'new Num_output:',Num_output,mon0
      end if
      fname(1:8)='fort.22.'
      fname(13:13)='-'
      fname(16:16)='-'
      write(fname(9:12),'(i4.4)') iyfm
      write(fname(14:15),'(i2.2)') mon0
      write(fname(17:18),'(i2.2)') number_day
      if ( mytid == 0 ) then
            open (17, file="rpointer.ocn", form='formatted')  
            write(17,'(a18)') fname
            close(17)
!            write(*,*)'fname=',fname
      end if

      totday = totday + 1
#if (defined LOWRES)
         CALL SSAVEINS
#endif
!      CALL SSAVEINSSTEP(NSS)
      CALL NEXTSTEP
     end if


!--------------------------------------------------------------------

  end subroutine ocn_run_mct
!***********************************************************************
!BOP
!
! !IROUTINE: ocn_final_mct
!
! !INTERFACE:
  subroutine ocn_final_mct( EClock, cdata_o, x2o_o, o2x_o)
!
! !DESCRIPTION:
! Finalize LICOM
!
! !ARGUMENTS:
    type(ESMF_Clock)            , intent(inout) :: EClock
    type(seq_cdata)             , intent(inout) :: cdata_o
    type(mct_aVect)             , intent(inout) :: x2o_o
    type(mct_aVect)             , intent(inout) :: o2x_o
!
!EOP
!BOC
!-----------------------------------------------------------------------
!
    if (mytid.eq.0) write(6,*)'write final fort num_step_per_day:',num_step_per_day,NSS
#if (defined LOWRES)
      call SSAVEINS
#else
    CALL SSAVEINSSTEP(num_step_per_day)
#endif
    call cleanup

  end subroutine ocn_final_mct


!*************************************************************************
!ROUTINE: ocn_SetGSMap_mct
!         lihuimin 2012.7.16
!
  subroutine ocn_SetGSMap_mct(mpicom_ocn, OCNID, gsMap_ocn)

    implicit none
    integer        , intent(in)    :: mpicom_ocn
    integer        , intent(in)    :: OCNID
    type(mct_gsMap), intent(inout) :: gsMap_ocn

    integer,allocatable :: &
      gindex(:)

    integer (kind(1)) ::   &
      i,j, k, n, iblock, &
      lsize, gsize,   &
      ier

    type (block) :: this_block          ! block information for current block

!-----------------------------------------------------------------------
!  Build the LICOM grid numbering for MCT
!  NOTE:  Numbering scheme is: West to East and South to North starting
!  at the south pole.  Should be the same as what's used in SCRIP
!-----------------------------------------------------------------------

    n = 0
    do iblock = 1, nblocks_clinic
       this_block = get_block(blocks_clinic(iblock),iblock)
       do j=this_block%jb,this_block%je
       do i=this_block%ib,this_block%ie
          n=n+1
       enddo
       enddo
    enddo
    lsize = n

! not correct for padding, use "n" above
!    lsize = block_size_x*block_size_y*nblocks_clinic
    gsize = imt_global*jmt_global
    allocate(gindex(lsize),stat=ier)

    n = 0
    do iblock = 1, nblocks_clinic
       this_block = get_block(blocks_clinic(iblock),iblock)
       do j=this_block%jb,this_block%je
       do i=this_block%ib,this_block%ie
          n=n+1
          gindex(n) = (jmt_global- this_block%j_glob(j))*(imt_global) + this_block%i_glob(i)
       enddo
       enddo
    enddo


    call mct_gsMap_init( gsMap_ocn, gindex, mpicom_ocn, OCNID, lsize, gsize )

    deallocate(gindex)

  end subroutine ocn_SetGSMap_mct

!*************************************************************************
!ROUTINE: ocn_domain_mct
!         lihuimin 2012.7.16
!
  subroutine ocn_domain_mct(lsize, gsMap_o, dom_o)

    implicit none
    integer        , intent(in)    :: lsize
    type(mct_gsMap), intent(in)    :: gsMap_o
    type(mct_ggrid), intent(inout) :: dom_o
    ! local
    integer :: j_begin,j_end
    integer, pointer :: &
      idata(:)

    real(r8), pointer :: &
      data(:)

    integer (kind(1)) ::    i,j, k, n, iblock, ier

    type (block) :: this_block          ! block information for current block

!
  integer            :: fid    ! nc domain file ID
  integer            :: dimid  ! nc dimension id
  integer            :: vid    ! nc variable ID
  integer            :: rcode  ! nc return code
  integer            :: ntim   ! temporary
 
!
!     Define Variables.
      integer*4   :: ncid, iret
      integer*4,  dimension(4) :: start(4)
      integer*4,  dimension(4) :: count(4)
      character (len=18) :: fname
      real(r8), allocatable :: data1(:,:),data2(:,:), temp1(:,:,:),temp2(:,:,:),temp3(:,:,:)
      INTEGER :: NMFF

!
  if (trim(horiz_grid_opt) == 'lat_lon') then
  if (mytid==0 ) then
     call wrap_open ('domain_licom.nc', NF_NOWRITE, fid)

     ! obtain dimensions
     call wrap_inq_dimid  (fid, 'ni' , dimid)
     call wrap_inq_dimlen (fid, dimid, nx   )
     call wrap_inq_dimid  (fid, 'nj' , dimid)
     call wrap_inq_dimlen (fid, dimid, ny   )
     call wrap_inq_dimid  (fid, 'nv' , dimid)
     call wrap_inq_dimlen (fid, dimid, nv   )

     if (nx/=imt_global) then
        write(16,*)"nx=",nx,"imt=",imt_global
        stop
     end if
  endif


  call shr_mpi_bcast(nx,mpi_comm_ocn,"nx")
  call shr_mpi_bcast(ny,mpi_comm_ocn,"ny")

     ! obtain grid variables
     allocate(data1(nx,ny), data2(nx,ny))
     allocate(temp1(imt,jmt,max_blocks_clinic),temp2(imt,jmt,max_blocks_clinic),temp3(imt,jmt,max_blocks_clinic))

  if (mytid==0 ) then
     call wrap_inq_varid(fid, 'xc' ,vid)
     call wrap_get_var_realx(fid, vid, data1)
     do j=1, ny
     do i=1, nx
        data2(i,j) = data1(i,1+ny-j)
     end do
     end do
  end if
      call scatter_global(temp1, data2, master_task, distrb_clinic, &
                          field_loc_center, field_type_scalar)

  if (mytid==0 ) then
     call wrap_inq_varid(fid, 'yc' ,vid)
     call wrap_get_var_realx(fid, vid, data1)
     do j=1, ny
     do i=1, nx
        data2(i,j) = data1(i,1+ny-j)
     end do
     end do
  end if
      call scatter_global(temp2, data2, master_task, distrb_clinic, &
                          field_loc_center, field_type_scalar)
!
  if (mytid==0 ) then
     call wrap_inq_varid(fid, 'area' ,vid)
     call wrap_get_var_realx(fid, vid, data1)
     do j=1, ny
     do i=1, nx
        data2(i,j) = data1(i,1+ny-j)
     end do
     end do
     call wrap_close(fid)
  end if
      call scatter_global(temp3, data2, master_task, distrb_clinic, &
                          field_loc_center, field_type_scalar)

     deallocate(data1,data2)
  end if
!
    call mct_gGrid_init( GGrid=dom_o, CoordChars=trim(seq_flds_dom_coord), &
       OtherChars=trim(seq_flds_dom_other), lsize=lsize )
    call mct_aVect_zero(dom_o%data)
    allocate(data(lsize))


    call mct_gsMap_orderedPoints(gsMap_o, mytid, idata)
    call mct_gGrid_importIAttr(dom_o,'GlobGridNum',idata,lsize)



    data(:) = shr_const_spval !-9999.0_R8 !linpf 2012Jul27
    call mct_gGrid_importRAttr(dom_o,"lat"  ,data,lsize)
    call mct_gGrid_importRAttr(dom_o,"lon"  ,data,lsize)
    call mct_gGrid_importRAttr(dom_o,"area" ,data,lsize)
    call mct_gGrid_importRAttr(dom_o,"aream",data,lsize)
    data(:) = shr_const_spval !0.0_R8  !linpf 2012Jul27
    call mct_gGrid_importRAttr(dom_o,"mask",data,lsize)
    call mct_gGrid_importRAttr(dom_o,"frac",data,lsize)

    if (mytid==0 ) then
      write(16,*)"domain_mct  special, id =  ", mytid,"   imt=",imt,"   jmt=",jmt,"   nx=",nx,"   ny=",ny,"lsize=",lsize
    endif
!
!-------------------------------------------------------------------
!
! Fill in correct values for domain components
!
!-------------------------------------------------------------------

    n=0
    do iblock = 1, nblocks_clinic
       this_block = get_block(blocks_clinic(iblock),iblock)
       do j=this_block%jb,this_block%je
       do i=this_block%ib,this_block%ie
          n=n+1
          if (trim(horiz_grid_opt) == 'lat_lon') then
            data(n) = temp1(i,j,iblock)
          else
            data(n) = TLON(i,j,iblock)/DEGtoRAD
          end if
       enddo
       enddo
    enddo
    call mct_gGrid_importRattr(dom_o,"lon",data,lsize)
!
    n=0
    do iblock = 1, nblocks_clinic
       this_block = get_block(blocks_clinic(iblock),iblock)
       do j=this_block%jb,this_block%je
       do i=this_block%ib,this_block%ie
          n=n+1
          if (trim(horiz_grid_opt) == 'lat_lon') then
            data(n) = temp2(i,j,iblock)
          else
            data(n) = TLAT(i,j,iblock)/DEGtoRAD
          end if
       enddo
       enddo
    enddo
    call mct_gGrid_importRattr(dom_o,"lat",data,lsize)

    n=0
    do iblock = 1, nblocks_clinic
       this_block = get_block(blocks_clinic(iblock),iblock)
       do j=this_block%jb,this_block%je
       do i=this_block%ib,this_block%ie
          n=n+1
          if (trim(horiz_grid_opt) == 'lat_lon') then
            data(n) = temp3(i,j,iblock)
          else
            data(n) = TAREA(i,j,iblock)/(radius*radius)
          end if
       enddo
       enddo
    enddo
    call mct_gGrid_importRattr(dom_o,"area",data,lsize)
!
     n=0
    do iblock = 1, nblocks_clinic
       this_block = get_block(blocks_clinic(iblock),iblock)
       do j=this_block%jb,this_block%je
       do i=this_block%ib,this_block%ie
          n=n+1
!         data(n) = float(mask_r(i,j,iblock))
          data(n) = float(KMT(i,j,iblock))
          if (data(n) > 1.0_r8) data(n) = 1.0_r8
       enddo
       enddo
    enddo
    call mct_gGrid_importRattr(dom_o,"mask",data,lsize)
    call mct_gGrid_importRattr(dom_o,"frac",data,lsize)

    deallocate(data)
    deallocate(idata)
    if (trim(horiz_grid_opt) == 'lat_lon') deallocate (temp1, temp2, temp3)

  end subroutine ocn_domain_mct

!*****************************************************************************
!ROUTINE:   CLEANUP
!           lihuimin, 2012.7.22
!
  subroutine cleanup

    ! allocated in inirun.F90
    deallocate(t_cpl, s_cpl, u_cpl, v_cpl, dhdx, dhdy, Qheat)
    !deallocate(taux, tauy, netsw, lat1, sen, lwup, lwdn, melth, salt, prec, evap, meltw, roff, ifrac, patm, duu10n)
    deallocate(taux, tauy, netsw, lat1, sen, lwup, lwdn, melth, salt, prec, evap, meltw, roff, patm, duu10n)
    deallocate(iceoff, snow1)

    ! allocated in inirun.F90
    ! deallocate(h0, u, v, at)
    deallocate(h0, u, v)

! allocated in rdriver.F90
    !deallocate(su3,sv3,psa3,tsa3,qar3,uva3,swv3,cld3,sss3,sst3,nswv3,dqdt3)
    deallocate(sss3)
!    deallocate(seaice3,runoff3)
!    deallocate(wspd3,wspdu3,wspdv3,lwv3,rain3,snow3)
#if ( defined TIDEMIX )
    deallocate(wave_dis)
#endif
#if (defined SOLARCHLORO)
      deallocate(chloro3)
#endif
#ifdef DAILYBUGDET 
      deallocate(axdaily,aydaily,azdaily,dxdaily,dzdaily,dt_diffdaily)
      deallocate(tenddaily,dt_convdaily,pendaily)
  if ( simple_assm .eqv. .true.) deallocate(dt_restoredaily)
#endif
  if (allocated(restore_at)) then
    deallocate(restore_at) !LPF20170730
  endif

  end subroutine cleanup



 end module ocn_comp_mct
