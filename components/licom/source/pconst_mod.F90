!  CVS: $Id: pconst_mod.F90,v 1.1.1.1 2004/04/29 06:22:39 lhl Exp $
module pconst_mod
!
#include <def-undef.h>
use precision_mod
use param_mod
!     -----------------------------------------------------------
!     Index Fields
!     -----------------------------------------------------------
!YU   real,dimension(:),allocatable:: dyr_global
      !integer,dimension(jmt):: j_global
      !integer,dimension(imt):: i_global

      integer,dimension(:), allocatable:: j_global
      integer,dimension(:), allocatable:: i_global
      integer :: ix,iy
       real(r8),dimension(imt,jmt,km,max_blocks_clinic),target :: vit,viv
      !real(r8),dimension(:,:,:,:), allocatable, target:: vit,viv

!     real(r8),dimension(imt_global,jmt,km):: vit_1d,viv_1d
      !real(r8),dimension(jmt_global):: ahv_back
      real(r8),dimension(:), allocatable:: ahv_back
!      integer,dimension(imt,jmt,max_blocks_clinic):: basin
!     integer,dimension(imt,jmt,max_blocks_clinic):: itnu
!lhl1204
      !integer,dimension(imt,jmt,max_blocks_clinic):: na
      integer,dimension(:,:,:), allocatable:: na
      real(r8) :: dfricmx,dwndmix
      real(r8) :: Pr_number  
#ifdef BCKMEX
      real(r8) :: diff_back_eq,diff_back_coef,diff_back_coef_max  
#endif

      ! lihuimin, 2012.7.15
!     integer :: i_num   ! actual grid number in i direction in this process
!     integer :: j_num   ! actual grid number in j direction in this process
!     integer :: i_f_num ! formal grid number in i dircetion 
!     integer :: j_f_num ! formal grid number in j direction !lhl1204
!
!
!     -----------------------------------------------------------
!     Grids
!     -----------------------------------------------------------
#if (defined NETCDF) || (defined ALL)
      !real(r4),dimension(imt_global):: lon
      !real(r4),dimension(jmt_global):: lat
      !real(r4),dimension(imt_global,jmt_global):: lon_o,lat_o
      !real(r4),dimension(imt_global,jmt_global):: ulon_o,ulat_o
      real(r4),dimension(km):: lev
      real(r4),dimension(km+1):: lev1

      real(r4),dimension(:), allocatable:: lon
      real(r4),dimension(:), allocatable:: lat

      real(r4),dimension(:,:),allocatable:: lon_o,lat_o
      real(r4),dimension(:,:),allocatable:: ulon_o,ulat_o

     ! real(r4),dimension(:), allocatable:: lev
     ! real(r4),dimension(:), allocatable:: lev1

#endif
!lhl090729
      !real(r8),dimension(s_imt,s_jmt):: s_lon
      !real(r8),dimension(s_imt,s_jmt):: s_lat
      real(r8),dimension(:,:), allocatable:: s_lon
      real(r8),dimension(:,:), allocatable:: s_lat
      real(r8),dimension(:,:), allocatable:: r_lon
      real(r8),dimension(:,:), allocatable:: r_lat
!lhl090729
       real(r8),dimension(km),target :: zkt,dzp,odzp,odzt
       real(r8),dimension(kmp1),target :: zkp
      ! real(r8),dimension(imt,jmt,max_blocks_clinic):: EBEA,EBEB,EBLA,EBLB,EPEA,EPEB,EPLA,EPLB
      ! real(r8),dimension(imt,jmt,max_blocks_clinic):: RRD1,RRD2
      
      !real(r8),dimension(:), allocatable,target:: zkt,dzp,odzp,odzt
      !real(r8),dimension(:), allocatable,target :: zkp
      real(r8),dimension(:,:,:), allocatable,target:: EBEA,EBEB,EBLA,EBLB,EPEA,EPEB,EPLA,EPLB
      real(r8),dimension(:,:,:), allocatable:: RRD1,RRD2
!lhl060506
       integer(4), dimension(600000) :: remap_ocn_blocks

#if ( defined SMAG)

#endif

!Yu

      real(r8),dimension(imt,jmt,max_blocks_clinic),target::ohbt,ohbu,dzph,hbx,hby
      real(r8),dimension(imt,jmt,max_blocks_clinic),target:: SNLAT
!     real(r8),dimension(jmt):: COSU,COST
!     real(r8),dimension(:),allocatable::COSU_global,COST_global
      integer,dimension(:),allocatable:: i_start,j_start
!     real(r8),dimension(imt)::CF1,CF2,SF1,SF2
!
!
!     -----------------------------------------------------------
!     Reference T S & coefficients for calculation of d(density)
!     -----------------------------------------------------------

      ! REAL(r8),target :: TO(KM),SO(KM),C(KM,9),PO(KM)
      REAL(r8), allocatable,target :: TO(:),SO(:),C(:,:),PO(:)

!lhl1204

!YU
!
!
!     -----------------------------------------------------------
!     Control Parameter
!     -----------------------------------------------------------
      INTEGER:: ISOP
!
!
!     -----------------------------------------------------------
!     Phycical Parameter
!     -----------------------------------------------------------
      real(r8)::amv,ahv,ahice
!lhl1204
#if (defined TIDEFORCE)
      real(r8),dimension(imt,jmt,max_blocks_clinic)::pboz0_save 
      real(r8):: CD_COASTAL,CD_INDIAN,CD_PACIFIC,CD_ATLANTIC,BETA_TIDE 
      real(r8),dimension(imt,jmt,max_blocks_clinic),target::CD_TIDE
#endif
      !real(r8),dimension(imt,jmt,km,max_blocks_clinic)::akmu,akmt
      !real(r8),dimension(imt,jmt,km,NTRA,max_blocks_clinic)::akt
      real(r8),dimension(:,:,:,:), allocatable,target ::akmu,akmt
      real(r8),dimension(:,:,:,:,:), allocatable,target::akt
!lhl1204
      real(r8),dimension(jmt)::AM,AH
      !real(r8),dimension(imt,jmt,km,max_blocks_clinic)::am3,ah3
!lhl
      !real(r8),dimension(imt,jmt,km,max_blocks_clinic)::amx,amy
      real(r8),dimension(:,:,:,:), allocatable::am3,ah3
! !lhl
      real(r8),dimension(:,:,:,:),allocatable::amx,amy
!lhl
      real(r8)::gamma
      real(r8)::iau_gamma
#if ( defined SMAG)
      real(r8):: D0,CP,C0F,TBICE,OD0,SAG,CAG,OD0CP,ASEA, &
                    VSEA,AFB1,AFB2,AFC1,AFC2,AFT1,AFT2,KARMAN,RR
#else
      real(r8)::  D0,CP,C0F,TBICE,OD0,SAG,CAG,OD0CP,ASEA, &
                    VSEA,AFB1,AFB2,AFC1,AFC2,AFT1,AFT2
#endif
      logical :: diag_msf, diag_bsf, diag_mth, diag_budget,test_input
!
!
      CHARACTER (LEN=3):: ABMON(12)
      CHARACTER (LEN=3):: ABMON1(12)
      CHARACTER (LEN=80):: out_dir
!
      REAL(r8):: DTB,DTC,DTS,DTB2,DTC2,ONBB,ONBC,ONCC
      INTEGER:: NBB,NCC,NSS,ISB,ISC,IST,MONTH
      INTEGER:: NMONTH(12),NNMONTH(12)
      INTEGER:: number_day, number_month
      !mxyu20241211 read restart year time
      INTEGER:: rst_year
      ! running time(s): time = ntimes * IDTS
      INTEGER :: NTIMES, r_day
      !mxyu241012 frequency of write fort file
      INTEGER :: NHIS, NSTEP,NDAYS
      INTEGER :: FREC, FSTEP ! steps interval for 3 hours forcing input
      logical :: fortenable
      !mxyu241012
      INTEGER :: F6hREC, F6hSTEP ! steps interval for 6 hours gfs forcing input
!
      ! lihuimin 2012.6.18, add REFDATE
      INTEGER :: NUMBER,NSTART,IY0,IYFM,MON0,MEND,IMD,IDAY,II,JJ,IO_HIST,IO_REST,rest_freq,hist_freq,REFDATE,boundary_restore
      INTEGER :: timeref
      integer :: klv
      integer :: kvt !use the layer for simple assimilation
!
      character (len=80) :: adv_momentum, adv_tracer
      integer(kind(1))   :: curr_ymd_licom     ! Current date YYYYMMDD
      integer(kind(1))   :: curr_ymd_cpl     ! Current date YYYYMMDD
      integer(kind(1))   :: yy_licom,mm_licom,dd_licom,tod_licom     ! year, month, day
      integer(kind(1))   :: yy_cpl,mm_cpl,dd_cpl,tod_cpl     ! year, month, day
      integer(kind(1))   :: year_ref,month_ref,day_ref,hour_ref     ! reference year,month,day,hour
      integer(kind(1))   :: imonth,ihour,iminute,isecond,tod,ymd,ymd_sync,tod_sync     ! year, month, day
!lhl20130419
      integer:: ocn_cpl_dt,licom_cpl_dt
      integer:: totday !total number days simulation counting
      integer:: yearadd,dayadd !number year and day will be added to the restart file 
      integer:: iyfmfnew !number forcing yr will be added to the restart file 
      integer:: first_step !test if first step 
      integer:: first_fort !output forcing field before calculation 
      integer:: num_step_per_day !test if first step 
      integer:: Num_output,Num_outputacc,Num_3houracc !number of accumulation for output
      integer:: Num_dailyoutput !number of accumulation for daily output 
      logical:: dts_accum !test if daily accumulation or T/S step accumula
      logical:: daily_accum !test if daily accumulation or T/S step accumula
      logical:: dailybudget_accum !test if daily accumulation or T/S step accumula
      logical:: simple_assm !using the simple assimilation(restoring) or not 
      logical :: yeartoggle
      logical:: enoi_iau !using the IAU  or not
#if (defined TIDEFORCEGW)
      real(r8),parameter  ::tidetune=0.5 
#endif
#if ( defined TIDEMIX )
       real(r8),parameter  ::local_mixing_fraction=0.33, mixing_ef=0.20, shelf_cutoff=1500., decay_scale=1000.
       real(r8),parameter  ::back_tidalmixing=1.0d-6,max_tidalmixing=1.0D-2, max_wavedis=0.1     !yuzp-2016/12/8
      !  real(r8),dimension(imt,jmt,km,max_blocks_clinic) :: fz_tide,ak_tide
      !  real(r8),dimension(imt,jmt,km,max_blocks_clinic) :: fztidal,richardson,wp3_tidal     !yuzp-2016/11/13
      !  real(r8),dimension(imt,jmt,km,max_blocks_clinic) :: ak_tide1     !yuzp-2016/11/19

       real(r8),dimension(:,:,:,:), allocatable,target :: fz_tide,ak_tide
       real(r8),dimension(:,:,:,:), allocatable,target :: fztidal,richardson,wp3_tidal     !yuzp-2016/11/13
       real(r8),dimension(:,:,:,:), allocatable :: ak_tide1     !yuzp-2016/11/19
#endif
#if ( defined CANUTOMIXOUT )
       real(r8),dimension(imt,jmt,km,max_blocks_clinic) ::wp1_canuto,wp2_canuto,wp3_canuto,wp4_canuto,&
                                                          wp5_canuto,wp6_canuto,wp7_canuto,wp8_canuto     !yuzp-2016/12/4

       real(r8),dimension(imt,jmt,km,max_blocks_clinic) ::wk1_canuto,wk2_canuto,wk3_canuto,wk4_canuto     !yuzp-2016/12/4

       real(r8),dimension(imt,jmt,max_blocks_clinic) ::fcor_canuto,fcort_canuto     !yuzp-2016/12/4

       real(r8),dimension(imt,jmt,km,max_blocks_clinic) ::wp12_canuto,wp13_canuto     !yuzp-2016/12/4
       real(r8),dimension(imt,jmt,max_blocks_clinic) ::wp10_canuto,wp11_canuto     !yuzp-2016/12/4
       real(r8),dimension(imt,jmt,km,max_blocks_clinic) ::alpha_canuto,beta_canuto     !yuzp-2016/12/4
#endif
contains 

subroutine pconst_init()
      implicit none
      ! integer,dimension(jmt):: j_global
      ! integer,dimension(imt):: i_global
      ! real(r8),dimension(imt,jmt,km,max_blocks_clinic):: vit,viv 
      ! real(r8),dimension(jmt_global):: ahv_back
      ! integer,dimension(imt,jmt,max_blocks_clinic):: na
      ! real(r4),dimension(imt_global):: lon
      ! real(r4),dimension(jmt_global):: lat
      !real(r4),dimension(imt_global,jmt_global):: lon_o,lat_o
      !real(r4),dimension(imt_global,jmt_global):: ulon_o,ulat_o
      ! real(r4),dimension(km):: lev
      ! real(r4),dimension(km+1):: lev1
      ! real(r8),dimension(km):: zkt,dzp,odzp,odzt
      ! real(r8),dimension(kmp1):: zkp
      ! real(r8),dimension(imt,jmt,max_blocks_clinic):: EBEA,EBEB,EBLA,EBLB,EPEA,EPEB,EPLA,EPLB
      ! real(r8),dimension(imt,jmt,max_blocks_clinic):: RRD1,RRD2 
      ! REAL(r8):: TO(KM),SO(KM),C(KM,9),PO(KM)
       ! real(r8),dimension(imt,jmt,km,max_blocks_clinic)::akmu,akmt
      ! real(r8),dimension(imt,jmt,km,NTRA,max_blocks_clinic)::akt
      ! real(r8),dimension(imt,jmt,km,max_blocks_clinic)::am3,ah3
      ! real(r8),dimension(imt,jmt,km,max_blocks_clinic)::amx,amy
      !  real(r8),dimension(imt,jmt,km,max_blocks_clinic) :: fz_tide,ak_tide
      !  real(r8),dimension(imt,jmt,km,max_blocks_clinic) :: ak_tide1     !yuzp-2016/11/19

      allocate(i_global(jmt))
      allocate(j_global(imt))
      !allocate(vit(imt,jmt,km, max_blocks_clinic))
      !allocate(viv(imt,jmt,km, max_blocks_clinic))
      allocate(ahv_back(jmt_global))
      allocate(na(imt,jmt, max_blocks_clinic))
!      allocate(lon(imt_global))
!      allocate(lat(jmt_global))
#if (defined LOWRES)
       allocate(lon_o(imt_global, jmt_global))
       allocate(lat_o(imt_global, jmt_global))
       allocate(ulon_o(imt_global, jmt_global))
       allocate(ulat_o(imt_global, jmt_global))
#endif
 !     allocate(lev(km))
 !     allocate(lev1(km+1))
 !     allocate(zkt(km))
 !     allocate(dzp(km))
 !     allocate(odzp(km))
 !     allocate(odzt(km))
 !     allocate(zkp(kmp1))
      allocate(EBEA(imt,jmt,max_blocks_clinic))
      allocate(EBEB(imt,jmt,max_blocks_clinic))
      allocate(EBLA(imt,jmt,max_blocks_clinic))
      allocate(EBLB(imt,jmt,max_blocks_clinic))
      allocate(EPEA(imt,jmt,max_blocks_clinic))
      allocate(EPEB(imt,jmt,max_blocks_clinic))
      allocate(EPLA(imt,jmt,max_blocks_clinic))
      allocate(EPLB(imt,jmt,max_blocks_clinic))
      allocate(RRD1(imt,jmt,max_blocks_clinic))      
      allocate(RRD2(imt,jmt,max_blocks_clinic))      
      allocate(TO(KM))
      allocate(SO(KM))
      allocate(C(KM, 9))
      allocate(PO(KM))
      allocate(akmu(imt,jmt,km, max_blocks_clinic))
      allocate(akmt(imt,jmt,km,max_blocks_clinic))
      allocate(akt(imt,jmt,km,NTRA,max_blocks_clinic))
#if ( defined SMAG)
      allocate(am3(imt,jmt,km,max_blocks_clinic))
      allocate(ah3(imt,jmt,km,max_blocks_clinic))
      allocate(amx(imt,jmt,km,max_blocks_clinic))
      allocate(amy(imt,jmt,km,max_blocks_clinic))
#endif     
      allocate(fz_tide(imt,jmt,km,max_blocks_clinic))
      allocate(ak_tide(imt,jmt,km,max_blocks_clinic))
!      allocate(fztidal(imt,jmt,km,max_blocks_clinic))
!      allocate(richardson(imt,jmt,km,max_blocks_clinic))
!      allocate(wp3_tidal(imt,jmt,km,max_blocks_clinic))
!      allocate(ak_tide1(imt,jmt,km,max_blocks_clinic))

      allocate(s_lon(s_imt, s_jmt))
      allocate(s_lat(s_imt, s_jmt))
      allocate(r_lon(r_imt, r_jmt))
      allocate(r_lat(r_imt, r_jmt))
      !real(r8),dimension(s_imt,s_jmt):: s_lon
      !real(r8),dimension(s_imt,s_jmt):: s_lat
end subroutine



end module pconst_mod
