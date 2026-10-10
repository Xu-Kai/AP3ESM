!  CVS: $Id: rdriver.F90,v 1.7 2003/08/25 07:47:52 lhl Exp $
!     ======================
      SUBROUTINE RDRIVER
!     ======================
 
#include <def-undef.h>
use param_mod
use pconst_mod
!use POP_HaloMod
!use POP_GridHorzMod
use domain
use constant_mod
use gather_scatter
use forc_mod
use work_mod
use dyn_mod, only : buffer
#ifdef SPMD
use msg_mod
#endif

#ifdef COUP
use shr_sys_mod
#endif 

      IMPLICIT NONE
#include <netcdf.inc>
!
!     Define Variables.
      integer*4   :: ncid1, iret,iblock
      integer*4,  dimension(4) :: start
      integer*4,  dimension(4) :: count
!lhl20130913
      integer*4,  dimension(3) :: start1
      integer*4,  dimension(3) :: count1
!LPF20161112
      integer*4,  dimension(2) :: start2
      integer*4,  dimension(2) :: count2
!lhl20130913
      real*4 xx_r4,yy_r4 !LPF20151005
      real*8 xx_r8,yy_r8 !LPF20151005
      real(r4),dimension(imt,jmt,max_blocks_clinic):: vars_r4 
!      REAL    :: WCOE (JMT),ABC
      integer*8 :: data_size
#ifndef COUP 
      allocate(su3(imt,jmt,12,max_blocks_clinic),&
               sv3(imt,jmt,12,max_blocks_clinic),&
               psa3(imt,jmt,12,max_blocks_clinic),&
               tsa3(imt,jmt,12,max_blocks_clinic),&
               qar3(imt,jmt,12,max_blocks_clinic),&
               uva3(imt,jmt,12,max_blocks_clinic))
      ! if (mytid==0) write(*,*)'rd-ok1'
      allocate(swv3(imt,jmt,12,max_blocks_clinic),&
               cld3(imt,jmt,12,max_blocks_clinic),&
               nswv3(imt,jmt,12,max_blocks_clinic),&
               dqdt3(imt,jmt,12,max_blocks_clinic))
      ! if (mytid==0) write(*,*)'rd-ok2'
      allocate(seaice3(imt,jmt,12,max_blocks_clinic),&
               runoff3(imt,jmt,12,max_blocks_clinic))
      ! if (mytid==0) write(*,*)'rd-ok3'
      allocate(wspd3(imt,jmt,12,max_blocks_clinic),&
               wspdu3(imt,jmt,12,max_blocks_clinic),&
               wspdv3(imt,jmt,12,max_blocks_clinic),&
               lwv3(imt,jmt,12,max_blocks_clinic),&
               rain3(imt,jmt,12,max_blocks_clinic),&
               snow3(imt,jmt,12,max_blocks_clinic))
#endif

#ifdef LOWRES
      allocate(sss3(imt,jmt,12,max_blocks_clinic),&
               sst3(imt,jmt,12,max_blocks_clinic))
      allocate(chloro3(imt,jmt,12,max_blocks_clinic))

#else
      allocate(sss3(imt,jmt,14,max_blocks_clinic),&
               sst3(imt,jmt,14,max_blocks_clinic))
      allocate(chloro3(imt,jmt,14,max_blocks_clinic))
#endif

#ifndef READ_SPLIT
      allocate(buffer_real4(imt_global,jmt_global))
      allocate(buffer(imt_global,jmt_global))
#endif

#if ( defined TIDEFORCEGW )
      allocate(depthstd(imt,jmt,max_blocks_clinic)) 
      allocate(CD_InWave(imt,jmt,max_blocks_clinic))
      CD_InWave=0.0D0
#endif
#if ( defined TIDEMIX )
      allocate(wave_dis(imt,jmt,max_blocks_clinic))
#endif

      if (mytid==0) write(16,*)"Beginning------RDRIVER! "

!!$OMP PARALLEL DO PRIVATE (K,J,I) collapse(4)
      DO iblock = 1, nblocks_clinic
       DO K = 1,KM
         DO J = 1,JMT
            DO I = 1,IMT
               WKA (I,J,K,iblock)= 0.0
            END DO
         END DO
       END DO
      END DO
!!$OMP END PARALLEL DO
 

!#undef READ_SPLIT
#ifdef READ_SPLIT
      data_size = imt * jmt * 8 * 43
      call c_fopen_read("ocn"//CHAR(0), "rdriver.dir/rdriver"//CHAR(0), data_size)
#endif

!-----------------------------------------------------------------------
!     SU3  : Sea surface zonal wind stress         (N/M**2)
!     SV3  : Sea surface meridional wind stress    (N/M**2)
!     PSA3 : Sea surface air pressure              (Pa)
!     SWV3 : Total net downward solar radiation    (W/M**2)
!    NSWV3 : None Solar flux                       (Wm-2)
!    DQDT3 : Dq/Dt                                 (WK-1m-2)
!     SST3 : Sea surface temperature               (Celsius)
!     SSS3 : Sea surface salinity                  (psu)
!   chloro3:chlorophll concentration               (mg m-3)
!-----------------------------------------------------------------------
 
!     READ FORCING FIELD
#ifdef SPMD
!
!===============================================
!input phc sss to substitute woa to use for the restoring
!===============================================
#if (defined SSSNORM) || (defined BOUNDNORTH)
#ifdef READ_SPLIT
      data_size = imt * jmt * 14 * 8
      call c_fread(sss3(1,1,1,1), data_size)
      do k = 1, 14
            sss3(:,:,k,:) = sss3(:,:,k,:) * VIT(:,:,1,:)
      end do

      call c_fread(sst3(1,1,1,1), data_size)
      do k = 1, 14
            sst3(:,:,k,:) = sst3(:,:,k,:) * VIT(:,:,1,:)
      end do
#else
! for not READ_SPLIT

#ifdef SSS_restore
      if(mytid==0)then
        iret=nf_open('sss_mercator_monthly.nc',nf_nowrite,ncid1)
        call check_err (iret)
      end if
#if (defined LOWRES)
      do k=1, 12
#else
      do k=1, 14
#endif
        start(1)=1 ; count(1)=imt_global
        start(2)=1 ; count(2)=jmt_global
        start(3)=1 ; count(3)=1
        start(4)=k ; count(4)=1
       
      if (mytid == 0 ) then
#if (defined LOWRES)
         iret=nf_get_vara_double(ncid1,  7,start,count,buffer)
         call check_err (iret)

         do j=1 ,jmt_global/2
         do i=1, imt_global
           xx_r8 = buffer(i,j)
           yy_r8 = buffer(i,jmt_global+1-j)
           buffer(i,jmt_global+1-j) = xx_r8
           buffer(i,j) = yy_r8
        end do
        end do
#else
         iret=nf_get_vara_real(ncid1,  7,start,count,buffer_real4)
         call check_err (iret)
      do j=1 ,jmt_global/2
        do i=1, imt_global
           xx_r8 = buffer_real4(i,j)
           yy_r8 = buffer_real4(i,jmt_global+1-j)
           buffer(i,jmt_global+1-j) = xx_r8
           buffer(i,j) = yy_r8
        end do
        end do

        if(mod(jmt_global,2) .eq. 1) then
            write(6,*) "sss3 the middle j value need to fill,k:",k
            do i = 1, imt_global
                buffer(i, jmt_global/2 + 1) = buffer_real4(i, jmt_global/2 + 1)
            end do
        endif

#endif
      end if

      call scatter_global(sss3(:,:,k,:), buffer, master_task, distrb_clinic, &
                          field_loc_center, field_type_scalar)
      !sss3(:,:,k,:)=vars_r4(:,:,:)*1.0D0
      !sss3(:,:,k,:)=vars_r4(:,:,:)*1.0D0
      sss3(:,:,k,:)=sss3(:,:,k,:)*VIT(:,:,1,:)
!      if(mytid == 0) write(*,*) "sss3 max/min:",k, maxval(sss3), minval(sss3)

      end do
       if (mytid == 0 ) then
         iret = nf_close (ncid1)
         call check_err (iret)
        write(6,*)'ok111 rdriver sss3'
       end if
#endif

#ifdef SST_restore
     if(mytid==0)then
             iret=nf_open('sst_ostia_monthly.nc',nf_nowrite,ncid1)
             call check_err (iret)
     end if
#if (defined LOWRES)
      do k=1, 12
#else
     do k=1, 14
#endif

        if (mytid == 0 ) then
#if (defined LOWRES)
        start(1)=1 ; count(1)=imt_global
        start(2)=1 ; count(2)=jmt_global
        start(3)=1 ; count(3)=1
        start(4)=k ; count(4)=1
         iret=nf_get_vara_double(ncid1,  7,start,count,buffer)
         call check_err (iret)

         do j=1 ,jmt_global/2
         do i=1, imt_global
           xx_r8 = buffer(i,j)
           yy_r8 = buffer(i,jmt_global+1-j)
           buffer(i,jmt_global+1-j) = xx_r8
           buffer(i,j) = yy_r8
        end do
        end do
#else

       start1(1)=1 ; count1(1)=imt_global
       start1(2)=1 ; count1(2)=jmt_global
       start1(3)=k ; count1(3)=1
            iret=nf_get_vara_real(ncid1,  6,start1,count1,buffer_real4)
            call check_err (iret)
            do j=1 ,jmt_global/2
               do i=1, imt_global
                  xx_r8 = buffer_real4(i,j)
                  yy_r8 = buffer_real4(i,jmt_global+1-j)
                  buffer(i,jmt_global+1-j) = xx_r8
                  buffer(i,j) = yy_r8
               end do
            end do
#endif
        end if

        call scatter_global(sst3(:,:,k,:), buffer, master_task, distrb_clinic, &
                           field_loc_center, field_type_scalar)
      sst3(:,:,k,:)=sst3(:,:,k,:)*VIT(:,:,1,:)
     end do
     if (mytid == 0 ) then
        iret = nf_close (ncid1)
        call check_err (iret)
        write(6,*)'ok222 rdriver sst3'
     end if
#endif

#endif     
! for READ_SPLIT

#endif

!===============================================
!input seaice
#ifdef FRC_CORE
       if(.false.) then !LPF20151026
      if (mytid == 0 ) then
      iret=nf_open('seaice.db.NSIDC.1979-2006.clim.monthmean.modelgrid.01x01.nc',nf_nowrite,ncid1)
      call check_err (iret)
      end if

      do k=1, 12
      start(1)=1 ; count(1)=imt_global
      start(2)=1 ; count(2)=jmt_global
      start(3)=1 ; count(3)=1
      start(4)=k ; count(4)=1
      if (mytid == 0 ) then
      iret=nf_get_vara_double(ncid1, 5,start,count,buffer)
      call check_err (iret)

       do j=1 ,jmt_global/2
        do i=1, imt_global
           xx_r8 = buffer(i,j)
           yy_r8 = buffer(i,jmt_global+1-j)
           buffer(i,jmt_global+1-j) = xx_r8
           buffer(i,j) = yy_r8
        end do
       end do
      end if

    call scatter_global(seaice3(:,:,k,:), buffer, master_task, distrb_clinic, &
                          field_loc_center, field_type_scalar)
      end do 

      if (mytid == 0 ) then
      iret = nf_close (ncid1)
      call check_err (iret)
      end if

!===============================================
!input runoff
      if (mytid == 0 ) then
      iret=nf_open('runoff.db.clim.modelgrid.01x01.nc',nf_nowrite,ncid1)
      call check_err (iret)
      end if

      do k=1, 1
      start(1)=1 ; count(1)=imt_global
      start(2)=1 ; count(2)=jmt_global
      start(3)=1 ; count(3)=1
      start(4)=k ; count(4)=1
      if (mytid == 0 ) then
      iret=nf_get_vara_double(ncid1, 4,start,count,buffer)
      call check_err (iret)
       
      do j=1 ,jmt_global/2
        do i=1, imt_global
           xx_r8 = buffer(i,j)
           yy_r8 = buffer(i,jmt_global+1-j)
           buffer(i,jmt_global+1-j) = xx_r8
           buffer(i,j) = yy_r8
        end do
       end do
      end if

      call scatter_global(runoff3(:,:,k,:), buffer, master_task, distrb_clinic, &
                          field_loc_center, field_type_scalar)
      end do 

      if (mytid == 0 ) then
      iret = nf_close (ncid1)
      call check_err (iret)
      end if
      end if !LPF20151026
#endif

#if ( defined TIDEMIX )
#ifdef READ_SPLIT
      data_size = imt * jmt * 8
      call c_fread(wave_dis(1,1,1), data_size)
      wave_dis(:,:,:)=wave_dis(:,:,:)*VIT(:,:,1,:)
#else 
!==============================================
!input annual mean wave dissipation
      !if(.false.) then !LPF20200307
      if (mytid == 0 ) then
      iret=nf_open('tidal_energy.nc',nf_nowrite,ncid1) !LPF20140424
      call check_err (iret)
      endif

      do k=1,1
      start2(1)=1 ; count2(1)=imt_global
      start2(2)=1 ; count2(2)=jmt_global
      !start(3)=1 ; count(3)=1
      !start(4)=k ; count(4)=1
      if (mytid == 0 ) then
        iret=nf_get_vara_double(ncid1, 5,start2,count2,buffer) !LPF201310026
        !iret=nf_get_vara_real(ncid1, 3,start,count,buffer_real4) !LPF20131026
        call check_err (iret)
        !buffer=buffer_real4 !LPF20131026
        do j=1 ,jmt_global/2
        do i=1, imt_global
           xx_r8 = buffer(i,j)
           yy_r8 = buffer(i,jmt_global+1-j)
           buffer(i,jmt_global+1-j) = xx_r8
           buffer(i,j) = yy_r8
        end do
        end do
      endif

    call scatter_global(wave_dis(:,:,:), buffer, master_task, distrb_clinic, &
                          field_loc_center, field_type_scalar)
      !wave_dis(:,:,:)=vars_r4(:,:,:)*1.0D0

      !call global_distribute(buffer,wave_dis)
      !write(500+mytid,*) wave_dis
      wave_dis(:,:,:)=wave_dis(:,:,:)*VIT(:,:,1,:)
      end do 

      if (mytid == 0 ) then
      iret = nf_close (ncid1)
      call check_err (iret)
       write(6,*)'ok333 rdriver wave_dis'
      endif
!      stop

      !endif !LPF20200307
       if(.false.)then
      if(mytid==0) then
      buffer=1.0D-4
       write(6,*)'ok111 rdriver wave_dis'
      call flush(6)
      endif
    call scatter_global(wave_dis(:,:,:), buffer, master_task, distrb_clinic, &
                          field_loc_center, field_type_scalar)
       endif
#endif
#endif

#if ( defined TIDEFORCEGW )
      if (mytid == 0 ) then 
      iret=nf_open('topography_std_remapbil.nc',nf_nowrite,ncid1) 
      call check_err (iret) 
      endif 

      do k=1,1 
      start2(1)=1 ; count2(1)=imt_global 
      start2(2)=1 ; count2(2)=jmt_global 
      if (mytid == 0 ) then 
      iret=nf_get_vara_double(ncid1, 5,start2,count2,buffer) 
      call check_err (iret) 

      do j=1 ,jmt_global/2 
        do i=1, imt_global 
           xx_r8 = buffer(i,j) 
           yy_r8 = buffer(i,jmt_global+1-j) 
           buffer(i,jmt_global+1-j) = xx_r8 
           buffer(i,j) = yy_r8 
        end do 
        end do 
      endif 

    call scatter_global(depthstd(:,:,:), buffer, master_task, distrb_clinic, &
                          field_loc_center, field_type_scalar) 
      end do 
      depthstd(:,:,:)=depthstd(:,:,:)*VIT(:,:,1,:)
      if (mytid == 0 ) then 
      iret = nf_close (ncid1) 
      call check_err (iret) 
      write(6,*)'ok222 rdriver depthstd'
      endif 
#endif
!===============================================
!input the chlorophyll concentration
!===============================================
#if (defined SOLARCHLORO)
#ifdef READ_SPLIT

      data_size = imt * jmt * 8 * 14
      call c_fread(chloro3(1,1,1,1), data_size)
      do k = 1, 14
            chloro3(:,:,k,:)=chloro3(:,:,k,:)*VIT(:,:,1,:)
      end do
#else 


      if (mytid == 0)  then
         write(*,*)'open chl data'
         iret=nf_open('MODEL_CHLFRC',nf_nowrite,ncid1)
        call check_err (iret)
      endif
#if (defined LOWRES)
      do k=1, 12
#else
      do k=1, 14
#endif
       start1(1)=1 ; count1(1)=imt_global
       start1(2)=1 ; count1(2)=jmt_global
       !start1(3)=1 ; count(3)=1
       start1(3)=k ; count1(3)=1
      if (mytid == 0)  then
        iret=nf_get_vara_double(ncid1,      6,start1,count1,buffer)
        !iret=nf_get_vara_real(ncid1,      5,start,count,buffer_real4)
        call check_err (iret)

        do j=1 ,jmt_global/2
        do i=1, imt_global
           xx_r8 = buffer(i,j)
           yy_r8 = buffer(i,jmt_global+1-j)
           buffer(i,jmt_global+1-j) = xx_r8
           buffer(i,j) = yy_r8
        end do
        end do
      endif

      call scatter_global(chloro3(:,:,k,:), buffer, master_task, distrb_clinic, &
                          field_loc_center, field_type_scalar)
      !chloro3(:,:,k,:)=vars_r4(:,:,:)*1.0D0 
      chloro3(:,:,k,:)=chloro3(:,:,k,:)*VIT(:,:,1,:)
      enddo 

      if (mytid == 0 ) then
      iret = nf_close (ncid1)
        call check_err (iret)
       write(6,*)'ok444 rdriver chloro3'
      endif
#endif
#endif
#ifdef READ_SPLIT

      call c_fclose()
#endif

!==============================================
!
!        write(*,'(i4,11f8.2)') mytid,((su3(i,j,1),i=190,200),j=10,20)
!
#else
!Yu
#if (!defined CDFIN)
      OPEN (90,FILE ='MODEL.FRC',STATUS ='OLD',FORM ='UNFORMATTED')
      READ (90) SWV3_io,NSWV3_io,DQDT3_io,SU3_io,SV3_io,SST3_io,SSS3_io
      CLOSE (90)
      OPEN (91,FILE ='MODEL_CHLFRC',STATUS ='OLD',FORM ='UNFORMATTED')
      READ (91) chloro3_io
      CLOSE (91)
!
#else
!----------------------------------------------------
! Open netCDF file.
!----------------------------------------------------
      iret=nf_open('MODEL.FRC',nf_nowrite,ncid1)
      call check_err (iret)

!----------------------------------------------------
!   Retrieve data
!----------------------------------------------------
      start(1)=1 ; count(1)=imt_global
      start(2)=1 ; count(2)=jmt_global
      start(3)=1 ; count(3)=1
      start(4)=1 ; count(4)=12

      iret=nf_get_vara_real(ncid1,   5,start,count,swv3_io)
      call check_err (iret)
      iret=nf_get_vara_real(ncid1,   6,start,count,nswv3_io)
      call check_err (iret)
      iret=nf_get_vara_real(ncid1,   7,start,count,dqdt3_io)
      call check_err (iret)
      iret=nf_get_vara_real(ncid1,   8,start,count,su3_io)
      call check_err (iret)
      iret=nf_get_vara_real(ncid1,   9,start,count,sv3_io)
      call check_err (iret)
      iret=nf_get_vara_real(ncid1,  10,start,count,sst3_io)
      call check_err (iret)
      iret=nf_get_vara_real(ncid1,  11,start,count,sss3_io)
      call check_err (iret)
!
      iret = nf_close (ncid1)
      call check_err (iret)

 
#if (defined SOLARCHLORO)
      iret=nf_open('MODEL_CHLFRC',nf_nowrite,ncid1)
      call check_err (iret)
      iret=nf_get_vara_real(ncid1,      5,start,count,chloro3_io)
      call check_err (iret)
      iret = nf_close (ncid1)
      call check_err (iret)
#endif
!
#endif
!$OMP PARALLEL
!$OMP WORKSHARE
      SWV3=SWV3_io
      NSWV3=NSWV3_io
      DQDT3=DQDT3_io
      SU3=SU3_io
      SV3=SV3_io
      SST3=SST3_io
      SSS3=SSS3_io
      chloro3=chloro3_io
!$OMP END WORKSHARE
!$OMP END PARALLEL
#endif
!SPMD
!Yu
      if (mytid==0) then
#ifdef COUP
      call shr_sys_flush(6)
#endif
      end if
 
!-----------------------------------------------------------------------
!     land grids of the forcing fields assigned to 0 
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!     salinity = (psu-35)*0.001
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!     reverse VS (southward is positive)
!     notice: the former program VS is reversed during preparing forcing field
!-----------------------------------------------------------------------
#ifndef COUP
!!$OMP PARALLEL DO PRIVATE (K,J,I) collapse(4)
      DO iblock = 1, nblocks_clinic
       DO K = 1,12
         DO J = 1,JMT
            DO I = 1,IMT
                SWV3(I,J,K,iblock)= SWV3(I,J,K,iblock)*VIT(I,J,1,iblock)
               NSWV3(I,J,K,iblock)=NSWV3(I,J,K,iblock)*VIT(I,J,1,iblock)
               DQDT3(I,J,K,iblock)=DQDT3(I,J,K,iblock)*VIT(I,J,1,iblock)
                 SU3(I,J,K,iblock)=  SU3(I,J,K,iblock)*VIV(I,J,1,iblock)
                 SV3(I,J,K,iblock)=  (-1)*SV3(I,J,K,iblock)*VIV(I,J,1,iblock)
                !SST3(I,J,K,iblock)= SST3(I,J,K,iblock) !*VIT(I,J,1) !LPF20140227
                !SSS3(I,J,K,iblock)= SSS3(I,J,K,iblock) !*VIT(I,J,1) !LPF20140227
                !SSS3 (I,J,K,iblock) = (SSS3 (I,J,K,iblock) -35.0D0)*0.001D0
                seaice3(I,J,K,iblock)= seaice3(I,J,K,iblock)*VIT(I,J,1,iblock)
                runoff3(I,J,K,iblock)= runoff3(I,J,K,iblock)*VIT(I,J,1,iblock)
! #if (defined SOLARCHLORO)
!         chloro3(I,J,K,iblock)= chloro3(I,J,K,iblock)*VIT(I,J,1,iblock) 
! #endif
                END DO
           END DO
         END DO
        END DO
#endif
        DO iblock = 1, nblocks_clinic
#if (defined LOWRES)
        do k=1, 12
#else
       DO K = 1,14
#endif
         DO J = 1,JMT
            DO I = 1,IMT
                SST3(I,J,K,iblock)= SST3(I,J,K,iblock) !*VIT(I,J,1) !LPF20140227
                SSS3 (I,J,K,iblock) = (SSS3 (I,J,K,iblock) -35.0D0)*0.001D0
#if (defined SOLARCHLORO)
        chloro3(I,J,K,iblock)= chloro3(I,J,K,iblock)*VIT(I,J,1,iblock) 
#endif
                END DO
           END DO
         END DO
        END DO
!!!$OMP END PARALLEL DO
!-----------------------------------------------------------------------
!     reverse VS (southward is positive)
!     notice: the former program VS is reversed during preparing forcing field
!-----------------------------------------------------------------------
 
!!!$OMP PARALLEL DO PRIVATE (K,J,I)
!      DO K = 1,12
!         DO J = 1,JMT
!            DO I = 1,IMT
!               SV3 (I,J,K,iblock) = -SV3 (I,J,K,iblock)
!            END DO
!         END DO
!      END DO
!!!$OMP END PARALLEL DO
 
 
!-----------------------------------------------------------------------
!     CALCULATING THE ANNUAL MEAN FORCING FIELD 
!-----------------------------------------------------------------------
 
#if (defined FRC_ANN)
      DO iblock = 1, nblocks_clinic 
      DO M = 1,12
!!$OMP PARALLEL DO PRIVATE (J,I) collapse(2)
         DO J = 1,JMT
            DO I = 1,IMT
               WKA (I,J,1,iblock)= WKA (I,J,1,iblock) + SU3 (I,J,M,iblock)
               WKA (I,J,2,iblock)= WKA (I,J,2,iblock) + SV3 (I,J,M,iblock)
               WKA (I,J,3,iblock)= WKA (I,J,3,iblock) + SSS3 (I,J,M,iblock)
               WKA (I,J,4,iblock)= WKA (I,J,4,iblock) + SWV3 (I,J,M,iblock)
               WKA (I,J,5,iblock)= WKA (I,J,5,iblock) + SST3 (I,J,M,iblock)
               WKA (I,J,6,iblock)= WKA (I,J,6,iblock) + NSWV3 (I,J,M,iblock)
               WKA (I,J,7,iblock)= WKA (I,J,7,iblock) + DQDT3 (I,J,M,iblock)
#if (defined SOLARCHLORO)
               WKA (I,J,8,iblock)= WKA (I,J,8,iblock) + chloro3 (I,J,M,iblock)
#endif
            END DO
         END DO
!!$OMP END PARALLEL DO
      END DO
      END DO
 
!!$OMP PARALLEL DO PRIVATE (M,J,I) collapse(4)
      DO iblock = 1, nblocks_clinic 
       DO M = 1,12
         DO J = 1,JMT
            DO I = 1,IMT
               SU3 (I,J,M,iblock) = WKA (I,J,1,iblock)/12.0D0
               SV3 (I,J,M,iblock) = WKA (I,J,2,iblock)/12.0D0
              SSS3 (I,J,M,iblock) = WKA (I,J,3,iblock)/12.0D0
              SWV3 (I,J,M,iblock) = WKA (I,J,4,iblock)/12.0D0
              SST3 (I,J,M,iblock) = WKA (I,J,5,iblock)/12.0D0
             NSWV3 (I,J,M,iblock) = WKA (I,J,6,iblock)/12.0D0
             DQDT3 (I,J,M,iblock) = WKA (I,J,7,iblock)/12.0D0
#if (defined SOLARCHLORO)
             chloro3 (I,J,M,iblock) = WKA (I,J,8,iblock)/12.0D0
#endif
            END DO
         END DO
       END DO
      END DO
!!$OMP END PARALLEL DO
 
#endif
!
      if (mytid==0)then
      write(16,*)"END-----------RDRIVER!"
!#ifdef COUP
!      call shr_sys_flush(6)
!#endif
      endif
#ifndef READ_SPLIT
    deallocate(buffer_real4,buffer)
#endif
      RETURN
      END SUBROUTINE RDRIVER
 
