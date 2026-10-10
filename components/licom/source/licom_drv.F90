#define LOGMSG()
!write(mytid+600,'(a,i4)')"LICOM",__LINE__
#define LOGMSGCarbon()
!write(600+mytid,'(a,i4)')"Licom CARBON",__LINE__

   use licom_sys_mod
   use work_mod
   use shr_const_mod,only:SHR_CONST_SPVAL !linpf 2012Jul26
   use POP_CommMod
   use domain
   use grid
   use blocks
   use gather_scatter
   use distribution

#include <def-undef.h>
#ifdef SW_ATHREAD
use sunway_mod
use sunway_init
#endif
use param_mod
use pconst_mod
use control_mod
!use control_licom_mod
use constant_mod
use msg_mod, only: tag_1d,tag_2d,tag_3d,tag_4d,nproc,status,mpi_comm_ocn
use tracer_mod
use pmix_mod
use forc_mod
!
#ifdef USE_OCN_CARBON
use carbon_mod
use cforce_mod
#endif

#ifdef GFS_FORCING
use gfs_forcing_mod
#endif

use output_mod !LPF20170728
use work_mod
use cdf_mod
use dyn_mod
!use sw_mod


  implicit none
#include <netcdf.inc>

#ifdef SW_ATHREAD
interface
   subroutine init_para_() bind(c,name='init_para_')
   end subroutine init_para_
end interface
#endif
!totalday: day interval from time ref
  integer :: mm,totalday
  REAL(r8):: t3a,t3b,t1,t2,clock0f,t3,t4,t4a,t4b,t4c,t4d
  real(r8):: jra_time,readyt_time,readyc_time,bclinc_time,barotr_time,tracer_time,ssave_time,nextstep_time
  real(r8):: tt1,tt2,tt3,tt4,tt5,tt6,tt7,tt8,tt9,tt10,tt11,tt12,tt13,tt14,tt15,tt16,tt17,tt18,tt19,tt20
  real(r8):: b1,b2,b3,b4,b5,b6,b7,b8,b9,b10,total_time,icesnow_time,energy_time
  integer(8) :: long_nx, long_ny, long_km

  INTEGER :: IPT1,IPT2
  INTEGER :: IPT14_1,IPT14_2
  INTEGER :: iblock
  INTEGER :: avg3counter
  REAL(r8):: FACTOR

  jra_time = 0
  readyt_time = 0
  readyc_time = 0
  bclinc_time = 0
  barotr_time = 0 
  tracer_time = 0
  ssave_time = 0
  nextstep_time = 0 
  mpi_comm_ocn=0
  icesnow_time = 0 
  energy_time = 0
 total_time=0
  cdate = 0
  sec = 0
  ierr = 0
  info_time = 0
totalday = 0
  yeartoggle = .false.

      call mpi_init (ierr)
      call mpi_comm_dup(mpi_comm_world, mpi_comm_ocn, ierr)
      call mpi_comm_rank (mpi_comm_ocn, mytid, ierr)
      if (mytid == 0 ) then
         open (16,file="ocn.log", form="formatted")
      end if
!
      if (mytid==0)  write(16,*)"End mpi_comm_rank"
      call mpi_comm_size (mpi_comm_ocn, nproc, ierr)
      if (mytid ==0) write(16,*) "MYTID=",mytid,"Number of Processors is",nproc
!YU

      my_task = mytid
      master_task = 0
      POP_mytask = mytid
      POP_mastertask = 0

#ifdef ENABLE_GPLT
    call gptl_initialize()
#endif

  !    call param_init()
      call cdf_init()
      call pconst_init()
      call grid_init_mem()
      call pmix_init()
      call work_init()
      call dyn_init()
      !call sw_init()
      call tracer_init()
    long_nx = nx_block
    long_ny = ny_block
    long_km = km
    call write_init(mpi_comm_ocn,long_nx, long_ny, long_km)
     ! call cesm_swlu_debug_init()
     ! call cesm_swlu_prof_init(0)

!---------------------------------------------------------------------
!     SET THE CONSTANTS USED IN THE MODEL
!---------------------------------------------------------------------

      tt1 = clock0f()
      CALL CONST
      tt2 = clock0f()
      b1 = tt2 -tt1
      call output_init()

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
   tt3 = clock0f()
   call init_domain_blocks
   tt4 = clock0f()  
   b2 = tt4 -tt3

#ifdef ENABLE_GPLT
    call t_startf("grid1_init"//CHAR(0))
#endif
   call init_grid1
#ifdef ENABLE_GPLT
    call t_stopf("grid1_init"//CHAR(0))
#endif



#ifdef READ_NCON_SPLIT
      call init_domain_distribution()
#else
call init_domain_distribution(KMT_G)
#endif

    call forc_init()

#ifdef JRA_FORCING
    call JRA_DAILY_INIT
#endif

#ifdef GFS_FORCING
    call gfs_forcing_init
#endif

#ifdef ENABLE_GPLT
    call t_startf("grid2_init"//CHAR(0))
#endif
   call init_grid2
#ifdef ENABLE_GPLT
    call t_stopf("grid2_init"//CHAR(0))
#endif
   
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
#ifdef ENABLE_GPLT
    call t_startf("inirun"//CHAR(0))
#endif
      CALL INIRUN
#ifdef ENABLE_GPLT
    call t_stopf("inirun"//CHAR(0))
#endif

#ifdef SHOW_TIME
      call run_time('INRUN')
#endif
    if (NSTART==1)then !initial run from observation T/S
      number_month=month_ref
      number_day=day_ref
      yy_licom=year_ref
      mm_licom=month_ref
      dd_licom=number_day
      iyfm = year_ref
      month = number_month
      if(number_month > 1) then
        do i=1,number_month-1
            imd  = nmonth(i)
            totalday = totalday+imd
        end do
        totalday = totalday + number_day - 1
      else
        totalday = 0
      endif
      r_day = r_day + number_day - 1
#ifdef GFS_FORCING
      number_day = 1
      r_day = NTIMES * int(DTS) / 86400
#endif
   else
      !mxyu 20241211 read directly from rpointer.ocn
      yy_licom = rst_year
      mm_licom = number_month
      dd_licom = number_day
        if(number_month > 1) then
            do i=1,number_month-1
                imd  = nmonth(i)
                totalday = totalday+imd
            end do
            totalday = totalday + number_day - 1
        else
            totalday = number_day - 1
        end if
      month = mm_licom  !number_month
      r_day = r_day + number_day - 1
#ifdef GFS_FORCING
      call calc_date_diff(year_ref, month_ref, day_ref, yy_licom, mm_licom,dd_licom, number_day)
      number_day = number_day+1
      r_day = number_day + NTIMES * int(DTS) / 86400 - 1
      totalday = number_day - 1
      if(mytid==0) write(*,*) "gfs restart timing:", number_day, r_day, totalday
#endif
   endif
      iyfm=yy_licom
      if (mytid==0) write (*,*) "iyfm,month,number_month,mm_licom,number_day,r_day,totalday: ",iyfm,month,number_month,mm_licom,number_day,r_day,totalday

    Num_output=NSS*(NMONTH (mm_licom) - dd_licom + 1)
    if(mytid==0) write(*,*)'the days left in current months, nss is Num_output= ',NMONTH (mm_licom) - dd_licom + 1,nss, Num_output

#ifdef USE_OCN_CARBON
      LOGMSGCarbon()
      CALL INIRUN_PT
      LOGMSGCarbon()
#ifdef SHOW_TIME
      call run_time('INIRUN_PT')
#endif
#endif


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

!---------------------------------------------------------------------
!     THERMAL CYCLE
!---------------------------------------------------------------------
       if(mytid==0) write(*,*) 'stepon begin, version: 05km (CPU)', nproc,'CORES'
#ifdef SW_ATHREAD
       call athread_init()
       CALL init_para_ !init sw slave info
       call readyt_sw_init
       call readyc_sw_init
       call barotr_sw_init
       call bclinc_sw_init
       call tracer_sw_init
#endif

!       DO MM = 1, number

#ifdef GFS_FORCING
        ! iau only run 1 day
        if(enoi_iau ) then
            number_day = 1
            r_day = 1
        endif
#endif

#ifdef IAU_FST
       number_day = 1
#endif
       DO iday =number_day, r_day

           mon0 = mod(month-1,12)+1
!           iyfm = iyfm + int((month-1)/12)
           imd  = nmonth(mon0)

           if (mytid==0) then
              write(6,*)'test-mm,iyfm,imd,number_month,mon0,month,number_day,r_day, dd_licom:', &
                        iyfm,imd,number_month,mon0,month,number_day,r_day,dd_licom
           end if
           t3a=clock0f()
           totalday =  totalday + 1
           !call COREIDEAL_DAILY
!            call mpi_barrier(mpi_comm_ocn,ierr)
            t4c=clock0f()
            if (mytid==0) write(6,*)'iday,number_day,mon0,dd_licom,totalday:',iday,number_day,mon0,dd_licom,totalday
#ifdef JRA_FORCING
            call read_daily_jra(totalday)
#endif
#ifdef GFS_FORCING
            call read_hour_gfs_forcing(iday)
#endif
!            call mpi_barrier(mpi_comm_ocn,ierr)
            t4a=clock0f()
           ! call cesm_swlu_prof_start(0)
            t1=clock0f()

!----------------------------------------------------------------------
!   time interplate  CHLOROPHYLL, chloro keep 12 month data
!   new sss/sst use 14 month data,the 1st is 12th month of last year, the 14th is 1st moth of next year
!----------------------------------------------------------------------
           IF ( dd_licom <= 15) THEN
             IPT1 = MON0-1
             IPT14_1 = IPT1 + 1
             IF (IPT1 == 0 ) then
                 IPT1 = 12
                 IPT14_1 = 1
             END IF
             IPT2 = MON0
             IPT14_2 = IPT2 + 1
             FACTOR = FLOAT (dd_licom -15)/ FLOAT (NMONTH (IPT1)) + 1
           ELSE
             IPT1 = MON0
             IPT2 = MOD (MON0,12) +1
             IPT14_1 = IPT1 + 1
             IPT14_2 = IPT2 + 1
             if(MON0 == 12) IPT14_2 = 14
             FACTOR = FLOAT (dd_licom -15)/ FLOAT (NMONTH (IPT1))
           END IF
           if (mytid==0) write(6,*)'MON0,IPT1,IPT2,IPT14_1,IPT14_2,FACTOR:',MON0,IPT1,IPT2,IPT14_1,IPT14_2,FACTOR

           DO iblock = 1, nblocks_clinic
             DO J = 1,JMT
                 DO I = 1,IMT
                     !RESTORE(I,J,1,2,iblock)=SSS(I,J,iblock) !LPF20151103 for restore SSSm to SSSo

#ifdef SSS_restore
                   SSS(I,J,iblock)=(SSS3(I,J,IPT14_2,iblock)-SSS3(I,J,IPT14_1,iblock))*FACTOR+SSS3(I,J,IPT14_1,iblock)
#endif
                   
#if (defined SOLARCHLORO)
                   chloro(I,J,iblock)= ( (chloro3 (I,J,IPT14_2,iblock) - chloro3 (I,J,IPT14_1,iblock))     &
                           * FACTOR + chloro3 (I,J,IPT14_1,iblock))
#endif
#ifdef SST_restore
                   SST(I,J,iblock)=(SST3(I,J,IPT14_2,iblock)-SST3(I,J,IPT14_1,iblock))*FACTOR+SST3(I,J,IPT14_1,iblock)
                           
#endif
                 END DO
               END DO
             END DO

#if (defined SOLARCHLORO)
       CALL  SW_ABSOR
#endif
#ifdef ENABLE_GPLT
    call t_startf("stepon"//CHAR(0))
#endif
       !NSS = 10
       DO II = 1, NSS

           if (II==1)call energy
           num_step_per_day=num_step_per_day+1
!           curr_ymd_licom=iyfm*10000+mon0*100+iday
!           if (mytid==0) write(6,*)'curr_ymd_cpl,tod_cpl',curr_ymd_cpl,tod_cpl
!---------------------------------------------------------------------
! set the ideal temperature and salinity for ideal simulation !LPF 20200304
!if the true value in the inirun, then this set will be undone!!
#if (defined TEST_IDEAL)
            if(.false.)  CALL setidealvalue
#endif

#ifdef JRA_FORCING
            ! 3hour jra forcing field input
            if(mod(II-1,fstep) .eq. 0) then
               frec = II / fstep + 1
               if (mytid .eq. 0) write(*,*) "frec, fstep, II: ", frec, fstep, II
               call jra_daily
            end if
#endif

#ifdef GFS_FORCING
            ! 6hour gfs forcing field input
            if(mod(II-1,f6hstep) .eq. 0) then
               F6hREC = II / f6hstep + 1
               if (mytid .eq. 0) write(*,*) "F6hREC, f6hstep, II: ", F6hREC, f6hstep, II
               call GFS_INTERPOLATE
            end if
#endif

#ifdef ENABLE_GPLT
            call t_startf("write_fort"//CHAR(0))
#endif

            !write fort file (NHIS=0 => fortenable=.false., nstep=0; skip to avoid mod/div by 0)
           if(fortenable .and. nstep > 0) then
              if(mod(II-1,nstep) .eq. 0) then
                  CALL SSAVEINSSTEP(II)
              end if
           end if

#ifdef ENABLE_GPLT
            call t_stopf("write_fort"//CHAR(0))
#endif

#if (defined AVG3HOUR || defined OBP_AVG3HOUR)
           avg3counter = II-fstep
           if(avg3counter > 0) then
             if(mod(avg3counter-1,fstep) .eq. 0) then
               if (mytid .eq. 0) write(*,*) "avg3counter,II,fstep: ", avg3counter,II,fstep
               CALL CALAVG3HOUR
               CALL SSAVE3HOURAVG(avg3counter)
             end if
           end if
#endif

!---------------------------------------------------------------------
!     COMPUTE DENSITY, BAROCLINIC PRESSURE AND THE RELAVANT VARIABLES
#ifdef ENABLE_GPLT
            call t_startf("READYT"//CHAR(0))
#endif
            CALL READYT
#ifdef ENABLE_GPLT
            call t_stopf("READYT"//CHAR(0))
#endif
!---------------------------------------------------------------------
!     BAROCLINIC & BAROTROPIC CYCLE
!---------------------------------------------------------------------
            DO JJ = 1,NCC

!     COMPUTE MOMENTUM ADVECTION, DIFFUSION & THEIR VERTICAL INTEGRALS
                if (mytid .eq. 0) then
                    if(mod(II-1,100) .eq. 0) write(6,*) "Ok 1, NSS is ", NSS,II
                end if
#ifdef ENABLE_GPLT
                call t_startf("READYC"//CHAR(0))
#endif
               CALL READYC
#ifdef ENABLE_GPLT
                call t_stopf("READYC"//CHAR(0))
#endif
!     PREDICTION OF BAROTROPIC MODE
#ifdef ENABLE_GPLT
                call t_startf("BAROTR"//CHAR(0))
#endif
               CALL BAROTR
#ifdef ENABLE_GPLT
                call t_stopf("BAROTR"//CHAR(0))
#endif
!     PREDICTION OF BAROCLINIC MODE
#ifdef ENABLE_GPLT
                call t_startf("BCLINC"//CHAR(0))
#endif
               CALL BCLINC
#ifdef ENABLE_GPLT
                call t_stopf("BCLINC"//CHAR(0))
#endif
            END DO

!--------------------------------------------------------------------
!     PREDICTION OF PASSIVE TRACER
!--------------------------------------------------------------------
#ifdef ENABLE_GPLT
            call t_startf("TRACER"//CHAR(0))
#endif
            CALL TRACER
#ifdef ENABLE_GPLT
            call t_stopf("TRACER"//CHAR(0))
#endif

#ifdef ENABLE_GPLT
            call t_startf("ICESNOW"//CHAR(0))
#endif
            CALL ICESNOW
#ifdef ENABLE_GPLT
            call t_stopf("ICESNOW"//CHAR(0))
#endif

#ifdef ENABLE_GPLT
            call t_startf("energy"//CHAR(0))
#endif
            call energy
#ifdef ENABLE_GPLT
            call t_stopf("energy"//CHAR(0))
#endif

!--------------------------------------------------------------------
!     PERFORM CONVECTIVE ADJUSTMENT IF UNSTABLE STRATIFICATION OCCURS
!--------------------------------------------------------------------
!       IF (mod(IST,(1800/IDTS))==0 ) CALL CONVADJ
#ifdef ENABLE_GPLT
            call t_startf("CONVADJ"//CHAR(0))
#endif
            IF (II==1) CALL CONVADJ
#ifdef ENABLE_GPLT
            call t_stopf("CONVADJ"//CHAR(0))
#endif
!--------------------------------------------------------------------
!Accumulate the data in a T/S time step interval or in a day for output
          ! if (mytid==0) write(6,*)'begin accumm'
#if (defined DAILYACC || defined MONTHACC)
             if (dts_accum .eqv. .true.) CALL ACCUMM
#endif
!LPF20170719
#if (defined DAILYACC)
              CALL DAILYACCUMM
              Num_outputacc=Num_outputacc+1
#endif
#if (defined AVG3HOUR || defined OBP_AVG3HOUR)
              CALL ACCUMM3HOUR
              Num_3houracc = Num_3houracc + 1
#endif
!LPF20170719
!--------------------------------------------------------------------
!       Finish thermal cycle loop in a couple interval
            if(II == NSS) then
#ifdef ENABLE_GPLT
                call t_stopf("stepon"//CHAR(0))
#endif
                t2=clock0f()
                if(mytid==0) write(*,*) ' (CPU) stepon test :',t2-t1, 'loops :',II,'/',NSS
            end if

   END DO

!     call cesm_swlu_prof_stop(0)
!     call cesm_swlu_prof_set_print_every(1)
!     call cesm_swlu_prof_print(0)


!--------------------------------------------------------------------

!--------------------------------------------------------------------
!          Diagnose if the energy is out of control
!--------------------------------------------------------------------

      if(num_step_per_day==NSS) then 
           num_step_per_day=0
           CALL ADDPS
      end if
!--------------------------------------------------------------------
!Accumulate the data in a T/S time step interval or in a day for output
          ! if (mytid==0) write(6,*)'begin accumm'
!            if ((.not.dts_accum).and.(num_step_per_day==0)) CALL ACCUMM
           !if (mytid==0) write(6,*)'end accumm'
!--------------------------------------------------------------------

       t4b=clock0f()

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
#ifdef ENABLE_GPLT
                call t_startf("write_daily"//CHAR(0))
#endif
                CALL SSAVEDAILYAVG
#ifdef ENABLE_GPLT
                call t_stopf("write_daily"//CHAR(0))
#endif
#endif
#if (defined AVG3HOUR || defined OBP_AVG3HOUR)
                avg3counter = NSS-fstep
                CALL CALAVG3HOUR
                CALL SSAVE3HOURAVG(avg3counter)
#endif
           endif

           !     ACCUMULATE SOME VARIABLES FOR MONTHLY OUTPUT
           IF (dd_licom ==imd .and.num_step_per_day==0) then
               month = month + 1
               if(month > 12) then
                   month = month -12
                   iyfm = iyfm + 1
                   yeartoggle = .true.
               end if
               mon0 = mod(month-1,12)+1
               !            number_day = 1
               dd_licom=0 !LPF20200315
               !          if (imt_global < 3000) CALL SSAVEMON !turnon onlyfor low
               if (mytid==0) write(6,*)'test-mm,month update to, mon0 is ',month,mon0
#if (defined DAILYACC || defined MONTHACC)
               call SSAVEMONTHAVG
#endif
               Num_output=NSS*NMONTH (MON0)
               yeartoggle = .false.
               if(mytid==0) write(*,*)'Num_output= ',Num_output
           ENDIF

           !mxyu 20241212 dd_licom now is used for the day in one month
           dd_licom = dd_licom + 1
!LPF20170726
!--------------------------------------------------------------------
!     SAVE FOR OUTPUT DAILY
     if ( num_step_per_day== 0) then
!         CALL SSAVEINS
!         CALL SSAVEINSSTEP(II)
       t4d=clock0f()
!       if(mytid==0) write(*,*) 'write_file time or ssaveins is  :',t4d-t4b,"s"
         t3=clock0f()
         CALL NEXTSTEP
         t4=clock0f()
         nextstep_time = nextstep_time +  t4 - t3
     end if
!--------------------------------------------------------------------
       if(mytid==0) then
           write(*,FMT='(A6,F20.15,A10,F20.15,A6, I2,I3)') ' day:',t4b-t4a,' SYPD:',86400./((t4b-t4a)*365.),'   ',MM,iday
           write(*,FMT='(A6,F20.15,A10,F20.15,A6, I2,I3)') ' day:',t4d-t4c,' SYPD:',86400./((t4d-t4c)*365.),'   ',MM,iday,"io"
           total_time = total_time + (t4b-t4a)
        endif


   end do
       !save the last fort file
       if(mytid==0) write(*,*) "save the last fort file"
       CALL SSAVEINSSTEP(0)
       call energy
       if(mytid==0) then
           write(*,*)"total time is ", total_time,"s"
       endif
 !      exit
!    end do

#ifdef ENABLE_GPLT
      call gptl_finalize(mpi_comm_ocn)
#endif

!#endif 
!--------------------------------------------------------------------
     call mpi_finalize(ierr)
     end 
