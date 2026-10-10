!  CVS: $Id: inirun.F90,v 1.1.1.1 2004/04/29 06:22:39 lhl Exp $
!     =================
      SUBROUTINE INIRUN
!     =================
!     INITIALIZING FOR ALL PHYSICAL FIELDS
 
#include <def-undef.h>
use param_mod
use pconst_mod
use dyn_mod
use tracer_mod
use forc_mod
use work_mod
use msg_mod, only: tag_1d,tag_2d,tag_3d,tag_4d,nproc,status,mpi_comm_ocn
!use shr_msg_mod
use shr_mpi_mod
use licom_sys_mod
use buf_mod
use control_mod
use domain
use constant_mod
use gather_scatter
use POP_HaloMod
use POP_GridHorzMod
use operators
use hmix_del2
use hmix_del4
use output_mod !LPF20170831
use diag_mod
use grid

      IMPLICIT NONE

  !----- local  ------
  integer            :: fid    ! nc domain file ID
  integer            :: dimid  ! nc dimension id
  integer            :: vid    ! nc variable ID
  integer            :: rcode  ! nc return code
  integer            :: ntim   ! temporary
  integer            :: iblock   ! temporary
  integer            :: ErrorCode   ! temporary
  real               :: abs_lat    ! xjw20250724
 
#include <netcdf.inc>
!
!     Define Variables.
      integer*4,  dimension(4) :: start(4)
      integer*4,  dimension(4) :: count(4)
      integer*4,  dimension(3) :: startiau(3)
      integer*4,  dimension(3) :: countiau(3)
      character (len=18) :: fname
      INTEGER :: NMFF
      real(r8) :: xx,yy
      integer*4   :: ncid1, iret1
#ifndef LOWRES
      integer*4 :: ncid,iret
#endif
        character(256) :: seq_no, file_name
        integer :: ioerr
        integer*8 :: io_size
        real(r8), allocatable:: ini_u(:, :, :, :)
        real(r8), allocatable:: ini_v(:, :, :, :)
        real(r8), allocatable:: uu_data(:, :, :)
        real(r8), allocatable:: vv_data(:, :, :)
       allocate(h0(imt,jmt,max_blocks_clinic),u(imt,jmt,km,max_blocks_clinic),v(imt,jmt,km,max_blocks_clinic))
!      allocate(h0(imt,jmt,max_blocks_clinic),u(imt,jmt,km,max_blocks_clinic),v(imt,jmt,km,max_blocks_clinic), &
!               at(imt,jmt,km,ntra,max_blocks_clinic))
      !allocate(atzwp(imt,jmt,km,ntra,max_blocks_clinic))
#ifdef TIDEFORCE
       allocate(elevat(imt,jmt,max_blocks_clinic)) 
       allocate(h0_tide(imt,jmt,max_blocks_clinic)) 
#endif

    if (mytid==0)then
        write(6,*)"BEGINNING-----INIRUN !"
    endif 

    ! xjw20250724
    if (mytid==0)then
        write(6,*)"Set the prescribed ratio profile for ss/tt/ssh diffusivity !"
    endif 

    DO iblock = 1, nblocks_clinic
      do j=1,jmt
        do i=1,imt
!           factor_vis_tsz0(i,j,iblock) = 1.0d0 - 0.9d0 * exp(-(((tlat(i,j,iblock)/DEGtoRAD)-35.0)**2/200.0)) &
!                                               - 0.9d0 * exp(-(((tlat(i,j,iblock)/DEGtoRAD)+40.0)**2/200.0))
          abs_lat = abs(tlat(i,j,iblock)/DEGtoRAD)

          if(abs_lat.lt.2)then
              factor_vis_tsz0(i,j,iblock) = 5.0
          else if((abs_lat.ge.2).and.(abs_lat.lt.5))then
              factor_vis_tsz0(i,j,iblock) = 3.0
          else if((abs_lat.ge.5).and.(abs_lat.lt.10))then
              factor_vis_tsz0(i,j,iblock) = 1.5
          else if((abs_lat.ge.10).and.(abs_lat.lt.15))then
              factor_vis_tsz0(i,j,iblock) = 0.5
          else if((abs_lat.ge.15).and.(abs_lat.lt.20))then
              factor_vis_tsz0(i,j,iblock) = 0.15
          else if((abs_lat.ge.20).and.(abs_lat.lt.25))then
              factor_vis_tsz0(i,j,iblock) = 0.05
          else if((abs_lat.ge.25).and.(abs_lat.lt.35))then
              factor_vis_tsz0(i,j,iblock) = 0.02
          else if((abs_lat.ge.35).and.(abs_lat.lt.55))then
              factor_vis_tsz0(i,j,iblock) = 0.015
          else
              factor_vis_tsz0(i,j,iblock) = 0.02
          end if

          if(abs_lat.lt.1)then
              factor_vis_uv(i,j,iblock) = 5.0
          else if((abs_lat.ge.1).and.(abs_lat.lt.2))then
              factor_vis_uv(i,j,iblock) = 4.0
          else if((abs_lat.ge.2).and.(abs_lat.lt.3))then
              factor_vis_uv(i,j,iblock) = 2.5
          else if((abs_lat.ge.3).and.(abs_lat.lt.4))then
              factor_vis_uv(i,j,iblock) = 2.0
          else if((abs_lat.ge.4).and.(abs_lat.lt.5))then
              factor_vis_uv(i,j,iblock) = 1.5
          else if((abs_lat.ge.5).and.(abs_lat.lt.6))then
              factor_vis_uv(i,j,iblock) = 1.2
          else
              factor_vis_uv(i,j,iblock) = 1.0
          end if

        end do 
      end do
    END DO

!    if (mytid==0)then
!        open(22,file='factor_vis_tsz0.dat',form='unformatted')
!    end if
!
!    call gather_global(factor_vis_tsz0_global,tlat, 0,distrb_clinic)
!    if (mytid==0)then
!        write(22) factor_vis_tsz0_global
!    end if
!
!    call gather_global(factor_vis_tsz0_global,tlon, 0,distrb_clinic)
!    if (mytid==0)then
!        write(22) factor_vis_tsz0_global
!    end if
!
!    call gather_global(factor_vis_tsz0_global,factor_vis_tsz0, 0,distrb_clinic)
!    if (mytid==0)then
!        write(22) factor_vis_tsz0_global
!        !close(22)
!    end if
!
!    call gather_global(factor_vis_tsz0_global,DXU, 0,distrb_clinic)
!    if (mytid==0)then
!        write(22) factor_vis_tsz0_global
!    end if
!
!    call gather_global(factor_vis_tsz0_global,DYU, 0,distrb_clinic)
!    if (mytid==0)then
!        write(22) factor_vis_tsz0_global
!        close(22)
!    end if
!
!    if (mytid==0)then
!        write(6,*)"Finish setting the prescribed ratio profile !"
!    endif

 
#ifdef COUP
  ! lihuimin, 2012.7.17, nx,ny --> imt,jmt
  allocate (t_cpl(imt,jmt,max_blocks_clinic))
  allocate (s_cpl(imt,jmt,max_blocks_clinic))
  allocate (u_cpl(imt,jmt,max_blocks_clinic))
  allocate (v_cpl(imt,jmt,max_blocks_clinic))
  allocate (dhdx(imt,jmt,max_blocks_clinic))
  allocate (dhdy(imt,jmt,max_blocks_clinic))
  allocate (Qheat   (imt,jmt,max_blocks_clinic))

  allocate (taux (imt,jmt,max_blocks_clinic))
  allocate (tauy (imt,jmt,max_blocks_clinic))
  allocate (netsw(imt,jmt,max_blocks_clinic))
  allocate (lat1 (imt,jmt,max_blocks_clinic))
  allocate (sen  (imt,jmt,max_blocks_clinic))
  allocate (lwup (imt,jmt,max_blocks_clinic))
  allocate (lwdn (imt,jmt,max_blocks_clinic))
  allocate (melth(imt,jmt,max_blocks_clinic))
  allocate (salt (imt,jmt,max_blocks_clinic))
  allocate (prec (imt,jmt,max_blocks_clinic))
  allocate (evap (imt,jmt,max_blocks_clinic))
  allocate (meltw(imt,jmt,max_blocks_clinic))
  allocate (roff (imt,jmt,max_blocks_clinic))
  allocate (iceoff (imt,jmt,max_blocks_clinic))
  allocate (snow1 (imt,jmt,max_blocks_clinic))
  !allocate (ifrac(imt,jmt,max_blocks_clinic)) 
  allocate (patm (imt,jmt,max_blocks_clinic)) 
  allocate (duu10n(imt,jmt,max_blocks_clinic))

   ! if(nstart==1) then
    if(nstart>=1) then !LPF20140423 nstart==1 initilized from OBS T/S nstart=2,initialized from any ocean run
         first_step=0
    else
         first_step=1
    endif
    num_step_per_day=0
    !num_cpl=1 !6

#endif

! if (mytid==0 ) then
!    call wrap_open ('domain_licom.nc', NF_NOWRITE, fid)

!    write(6,*) 'read domain data...'
!    call shr_sys_flush(6)
!    call wrap_inq_varid(fid, 'mask', vid)
!    call wrap_get_var_int(fid,vid,mask_g)
!    do j=1,jmt_global
!    do i=1,imt_global
!       mmm(i,j) = mask_g(i,jmt_global+1-j)
!    end do
!    end do
! end if
! call scatter_global(mask_r, mmm, master_task, distrb_clinic, &
!                         field_loc_center, field_type_scalar)

     ub(:,:,:) = 0.0D0
     vb(:,:,:) = 0.0D0
     h0(:,:,:) = 0.0D0
     ubp = 0.0D0
     vbp(:,:,:) = 0.0D0
     h0p = 0.0D0
#ifdef TIDEFORCE
     elevat = 0.0D0 
     h0_tide = 0.0D0 
#endif    
     U    = 0.0D0
     V    = 0.0D0
     UP   = 0.0D0
     VP   = 0.0D0
     WS   = 0.0D0
    wsbak = 0.0D0
     H0L  = 0.0D0
     H0F  = 0.0D0
     H0BL = 0.0D0
     H0BF = 0.0D0
     UTL  = 0.0D0
     UTF  = 0.0D0
     VTL  = 0.0D0
     VTF  = 0.0D0
     NET  = 0.0D0
     PXB  = 0.0D0
     PYB  = 0.0D0
     PAX  = 0.0D0
     PAY  = 0.0D0
     WHX  = 0.0D0
     WHY  = 0.0D0
     WGP  = 0.0D0
!     ------------------------------------------------------------------
!     Output Arrays
!     ------------------------------------------------------------------
      CALL YY00

#if (defined DAILYACC)
        call DAILYAVGINI
        Num_outputacc=0
#endif
#if (defined DAILYACC || defined MONTHACC)
        call MM00
#endif
#if (defined AVG3HOUR || defined OBP_AVG3HOUR)
        call AVG3HOURINI
        Num_3houracc = 0
#endif
!
      if (mytid==0) then
          write(*,*)"wjl debug 1"
      end if
 
 
      IF (NSTART == 1 .or. boundary_restore > 0 ) THEN
!
      MONTH = 1
      allocate(buffer(imt_global,jmt_global))
 
     
!

#if (defined LOWRES)
      if (mytid==0) then
          iret=nf_open('temp_initial',nf_nowrite,ncid)
          call check_err (iret)

          iret1=nf_open('salt_initial',nf_nowrite,ncid1)
          call check_err (iret1)
      end if
!------------------------------------------------------------------
!     READ LEVITUS ANNUAL MEAN TEMPERATURE AND SALINITY
!----------------------------------------------------
! Open netCDF file.
!----------------------------------------------------
!   Retrieve data
!----------------------------------------------------
     
      do k=1,km !km
!
       if (mytid == 0) then
        start(1)=1 ; count(1)=imt_global
        start(2)=1 ; count(2)=jmt_global
        start(3)=k ; count(3)=1
        start(4)=1 ; count(4)=1

        iret=nf_get_vara_double(ncid,   7,start,count, buffer)
        call check_err (iret)
        do j=1 ,jmt_global/2
        do i=1, imt_global
           xx = buffer(i,j)
           yy = buffer(i,jmt_global+1-j)
           buffer(i,jmt_global+1-j) = xx
           buffer(i,j) = yy
        end do
        end do
!
!       do j=1 ,jmt_global
!       do i=1, imt_global/2
!          xx = buffer(i,j)
!          yy = buffer(imt_global/2+i,j)
!          buffer(imt_global/2+i,j) = xx
!          buffer(i,j) = yy
!       end do
!       end do
       end if
!
      call scatter_global(at(:,:,k,1,:), buffer, master_task, distrb_clinic, &
                          field_loc_center, field_type_scalar)

!
       if (mytid == 0 ) then
        iret=nf_get_vara_double(ncid1,   7,start,count, buffer)
        call check_err (iret)
        do j=1 ,jmt_global/2
        do i=1, imt_global
           xx = buffer(i,j)
           yy = buffer(i,jmt_global+1-j)
           buffer(i,jmt_global+1-j) = xx
           buffer(i,j) = yy
        end do
        end do
!
!       do j=1 ,jmt_global
!       do i=1, imt_global/2
!          xx = buffer(i,j)
!          yy = buffer(imt_global/2+i,j)
!          buffer(imt_global/2+i,j) = xx
!          buffer(i,j) = yy
!       end do
!       end do
       end if
!
       call scatter_global(at(:,:,k,2,:), buffer,master_task, distrb_clinic, &
                          field_loc_center, field_type_scalar)
!
       end do !km
       if (mytid == 0 ) then
         iret = nf_CLOSE (ncid)
         iret = nf_CLOSE (ncid1)
       end if
#endif

#ifdef IAU_FST_SPLIT
       io_size = nx_block * ny_block * km * 8

       call c_fopen_read("ocn"//CHAR(0), "temp_iau.dir/temp_initial"//CHAR(0), io_size)
       call c_fread(at_iau(1,1,1,1,1), io_size)
       call c_fclose()

       call c_fopen_read("ocn"//CHAR(0), "salt_iau.dir/salt_initial"//CHAR(0), io_size)
       call c_fread(at_iau(1,1,1,2,1), io_size)
       call c_fclose()

      call c_fopen_read("ocn"//CHAR(0), "u_iau.dir/u"//CHAR(0), io_size)
      call c_fread(uu_iau(1,1,1,1), io_size)
      call c_fclose()

      call c_fopen_read("ocn"//CHAR(0), "v_iau.dir/v"//CHAR(0), io_size)
      call c_fread(vv_iau(1,1,1,1), io_size)
      call c_fclose()


      io_size = nx_block * ny_block  * 8
      call c_fopen_read("ocn"//CHAR(0), "h0_iau.dir/h0"//CHAR(0), io_size)
      call c_fread(h0_iau(1,1,1), io_size)
      call c_fclose()

    call POP_HaloUpdate(h0_iau , POP_haloClinic, POP_gridHorzLocCenter,&
            POP_fieldKindScalar, errorCode, fillValue = 0.0_r8)
    call POP_HaloUpdate(uu_iau , POP_haloClinic, POP_gridHorzLocSWcorner , &
            POP_fieldKindVector, errorCode, fillValue = 0.0_r8)
    call POP_HaloUpdate(vv_iau , POP_haloClinic, POP_gridHorzLocSWcorner , &
            POP_fieldKindVector, errorCode, fillValue = 0.0_r8)
    call POP_HaloUpdate(at_iau(:,:,:,1,:) , POP_haloClinic, POP_gridHorzLocCenter , &
            POP_fieldKindScalar, errorCode, fillValue = 0.0_r8)
    call POP_HaloUpdate(at_iau(:,:,:,2,:) , POP_haloClinic, POP_gridHorzLocCenter , &
            POP_fieldKindScalar, errorCode, fillValue = 0.0_r8)

#else
      if(enoi_iau ) then
          if (mytid==0) then
              write(*,*)"read iau"
          end if

          allocate(uu_data(imt,jmt,max_blocks_clinic),vv_data(imt,jmt,max_blocks_clinic))
          if (mytid==0) then
           iret=nf_open('salt_iau',nf_nowrite,ncid)
           call check_err (iret)

           iret1=nf_open('temp_iau',nf_nowrite,ncid1)
           call check_err (iret1)
          end if
!
           allocate(buffer_real4(imt_global,jmt_global))
          do k=1,km !km
!
           if (mytid == 0) then
            start(1)=1 ; count(1)=imt_global
            start(2)=1 ; count(2)=jmt_global
            start(3)=k ; count(3)=1
            start(4)=1 ; count(4)=1

            iret=nf_get_vara_real(ncid,   7,start,count, buffer_real4)
            call check_err (iret)
            do j=1 ,jmt_global/2
            do i=1, imt_global
               xx = buffer_real4(i,j)
               yy = buffer_real4(i,jmt_global+1-j)
               buffer(i,jmt_global+1-j) = xx
               buffer(i,j) = yy
            end do
            end do
            !buffer(:,:)=buffer_real4(:,:)*1000+35.0
            write(*,*) "salt_iau, min/max:", minval(buffer),maxval(buffer),k
          end if

          call scatter_global(at_iau(:,:,k,2,:), buffer,master_task, distrb_clinic, &
                              field_loc_center, field_type_scalar)

           if (mytid == 0) then
            start(1)=1 ; count(1)=imt_global
            start(2)=1 ; count(2)=jmt_global
            start(3)=k ; count(3)=1
            start(4)=1 ; count(4)=1

            iret1=nf_get_vara_real(ncid1,   7,start,count, buffer_real4)
            call check_err (iret)
            do j=1 ,jmt_global/2
            do i=1, imt_global
               xx = buffer_real4(i,j)
               yy = buffer_real4(i,jmt_global+1-j)
               buffer(i,jmt_global+1-j) = xx
               buffer(i,j) = yy
            end do
            end do
            write(*,*) "temp_iau, min/max:", minval(buffer),maxval(buffer),k
           ! buffer(:,:)=buffer_real4(:,:)
          end if

          call scatter_global(at_iau(:,:,k,1,:), buffer,master_task, distrb_clinic, &
                          field_loc_center, field_type_scalar)


         end do !K
         call POP_HaloUpdate(at_iau(:,:,:,1,:), POP_haloClinic, POP_gridHorzLocSWcorner , &
                           POP_fieldKindVector, errorCode, fillValue = 0.0_r8)
         call POP_HaloUpdate(at_iau(:,:,:,2,:), POP_haloClinic, POP_gridHorzLocSWcorner , &
                           POP_fieldKindVector, errorCode, fillValue = 0.0_r8)
          if (mytid==0) then
              iret = nf_CLOSE (ncid)
              iret1 = nf_CLOSE (ncid1)
          end if

          if (mytid==0) then
           iret=nf_open('uu_iau',nf_nowrite,ncid)
           call check_err (iret)
           iret=nf_open('vv_iau',nf_nowrite,ncid1)
           call check_err (iret)
          end if
         do k=1,km !km
!
           if (mytid == 0) then
            start(1)=1 ; count(1)=imt_global
            start(2)=1 ; count(2)=jmt_global
            start(3)=k ; count(3)=1
            start(4)=1 ; count(4)=1

            iret=nf_get_vara_real(ncid,   7,start,count, buffer_real4)
            call check_err (iret)
            do j=1 ,jmt_global/2
            do i=1, imt_global
               xx = buffer_real4(i,j)
               yy = buffer_real4(i,jmt_global+1-j)
               buffer(i,jmt_global+1-j) = xx
               buffer(i,j) = yy
            end do
            end do
            !buffer(:,:)=buffer_real4(:,:)*1000+35.0
            write(*,*) "uu_iau, min/max:", minval(buffer),maxval(buffer),k
          end if

           call scatter_global(uu_data(:,:,:), buffer,master_task, distrb_clinic, &
                              field_loc_center, field_type_scalar)

           if (mytid == 0) then
               start(1)=1 ; count(1)=imt_global
               start(2)=1 ; count(2)=jmt_global
               start(3)=k ; count(3)=1
               start(4)=1 ; count(4)=1

               iret=nf_get_vara_real(ncid1,   7,start,count, buffer_real4)
               call check_err (iret)
               do j=1 ,jmt_global/2
                   do i=1, imt_global
                       xx = buffer_real4(i,j)
                       yy = buffer_real4(i,jmt_global+1-j)
                       buffer(i,jmt_global+1-j) = xx
                       buffer(i,j) = yy
                   end do
               end do
               !buffer(:,:)=buffer_real4(:,:)*1000+35.0
               write(*,*) "vv_iau, min/max:", minval(buffer),maxval(buffer),k
           end if

           call scatter_global(vv_data(:,:,:), buffer,master_task, distrb_clinic, &
                   field_loc_center, field_type_scalar)

           uu_iau(:,:,k,1) = uu_data(:,:,1)*cos(local_rho(:,:)) + vv_data(:,:,1)*sin(local_rho(:,:))
           vv_iau(:,:,k,1) = -uu_data(:,:,1)*sin(local_rho(:,:)) + vv_data(:,:,1)*cos(local_rho(:,:))

         end do !K
         call POP_HaloUpdate(uu_iau, POP_haloClinic, POP_gridHorzLocSWcorner , &
                           POP_fieldKindVector, errorCode, fillValue = 0.0_r8)

          call POP_HaloUpdate(vv_iau, POP_haloClinic, POP_gridHorzLocSWcorner , &
                  POP_fieldKindVector, errorCode, fillValue = 0.0_r8)

          if (mytid==0) then
            iret = nf_CLOSE (ncid)
            iret = nf_CLOSE (ncid)
          end if


          if (mytid==0) then
           iret=nf_open('ssh_iau',nf_nowrite,ncid)
           call check_err (iret)
          end if
!
           if (mytid == 0) then
            startiau(1)=1 ; countiau(1)=imt_global
            startiau(2)=1 ; countiau(2)=jmt_global
            startiau(3)=1 ; countiau(3)=1

            iret=nf_get_vara_real(ncid,   6,startiau,countiau, buffer_real4)
            call check_err (iret)
            do j=1 ,jmt_global/2
            do i=1, imt_global
               xx = buffer_real4(i,j)
               yy = buffer_real4(i,jmt_global+1-j)
               buffer(i,jmt_global+1-j) = xx
               buffer(i,j) = yy
            end do
            end do
            !buffer(:,:)=buffer_real4(:,:)*1000+35.0
            write(*,*) "ssh_iau, min/max:", minval(buffer),maxval(buffer)
          end if

           call scatter_global(h0_iau, buffer,master_task, distrb_clinic, &
                              field_loc_center, field_type_scalar)

         call POP_HaloUpdate(h0_iau, POP_haloClinic, POP_gridHorzLocSWcorner , &
                           POP_fieldKindVector, errorCode, fillValue = 0.0_r8)
          if (mytid==0) then
            iret = nf_CLOSE (ncid)
          end if
          deallocate(buffer_real4)
          deallocate(uu_data)
          deallocate(vv_data)
  endif

#endif

#if (defined HIGHRES)
#define READ_SPLIT
#define U_V_INITIAL
#ifdef READ_SPLIT

#ifdef U_V_INITIAL
        allocate(ini_u(imt,jmt,km,max_blocks_clinic))
        allocate(ini_v(imt,jmt,km,max_blocks_clinic))
#endif
       io_size = nx_block * ny_block * km * 8
      
       call c_fopen_read("ocn"//CHAR(0), "temp.dir/temp_initial"//CHAR(0), io_size)
       call c_fread(at(1,1,1,1,1), io_size)
       call c_fclose()

       call c_fopen_read("ocn"//CHAR(0), "salt.dir/salt_initial"//CHAR(0), io_size)
       call c_fread(at(1,1,1,2,1), io_size)
       call c_fclose()

#ifdef U_V_INITIAL
      call c_fopen_read("ocn"//CHAR(0), "u.dir/u"//CHAR(0), io_size)
#ifdef COUP
      call c_fread(u(1,1,1,1), io_size)
#else
      call c_fread(ini_u(1,1,1,1), io_size)
#endif
      call c_fclose()

      call c_fopen_read("ocn"//CHAR(0), "v.dir/v"//CHAR(0), io_size)
#ifdef COUP
      call c_fread(v(1,1,1,1), io_size)
#else
      call c_fread(ini_v(1,1,1,1), io_size)
#endif
      call c_fclose()


      io_size = nx_block * ny_block  * 8
      call c_fopen_read("ocn"//CHAR(0), "h0.dir/h0"//CHAR(0), io_size)
      call c_fread(h0(1,1,1), io_size)
      call c_fclose()

#ifndef COUP
      do k=1,km
        u(:,:,k,1) = ini_u(:,:,k, 1)*cos(local_rho(:,:)) + ini_v(:,:,k,1)*sin(local_rho(:,:))
        v(:,:,k,1) = -ini_u(:,:,k, 1)*sin(local_rho(:,:)) + ini_v(:,:,k,1)*cos(local_rho(:,:))
      end do
#endif
      deallocate(ini_u)
      deallocate(ini_v)
#endif
!U_V_INITIAL

#else
!READ_SPLIT


      if (mytid==0) then
          write(*,*)"wjl debug 2"
      end if
      if (mytid==0) then
       iret=nf_open('salt_initial',nf_nowrite,ncid)
       call check_err (iret)

       iret1=nf_open('temp_initial',nf_nowrite,ncid1)
       call check_err (iret1)
      end if
!
      if(.true.)then !work-20160501
       allocate(buffer_real4(imt_global,jmt_global))
      do k=1,km !km
!
       if (mytid == 0) then
        start(1)=1 ; count(1)=imt_global
        start(2)=1 ; count(2)=jmt_global
        start(3)=k ; count(3)=1
        start(4)=1 ; count(4)=1

        iret=nf_get_vara_real(ncid,   7,start,count, buffer_real4)
        call check_err (iret)
        do j=1 ,jmt_global/2
        do i=1, imt_global
           xx = buffer_real4(i,j)
           yy = buffer_real4(i,jmt_global+1-j)
           buffer(i,jmt_global+1-j) = xx
           buffer(i,j) = yy
        end do
        end do
        !buffer(:,:)=buffer_real4(:,:)*1000+35.0
      end if

       call scatter_global(at(:,:,k,2,:), buffer,master_task, distrb_clinic, &
                          field_loc_center, field_type_scalar)
      
       if (mytid == 0) then
        start(1)=1 ; count(1)=imt_global
        start(2)=1 ; count(2)=jmt_global
        start(3)=k ; count(3)=1
        start(4)=1 ; count(4)=1

        iret1=nf_get_vara_real(ncid1,   7,start,count, buffer_real4)
        call check_err (iret)
        do j=1 ,jmt_global/2
        do i=1, imt_global
           xx = buffer_real4(i,j)
           yy = buffer_real4(i,jmt_global+1-j)
           buffer(i,jmt_global+1-j) = xx
           buffer(i,j) = yy
        end do
        end do

       ! buffer(:,:)=buffer_real4(:,:)
      end if

       call scatter_global(at(:,:,k,1,:), buffer,master_task, distrb_clinic, &
                          field_loc_center, field_type_scalar)


     end do !K
      deallocate(buffer_real4)

     endif !work-20160501

!
      if(.false.)then !nowork
      do k=1,km !km
!
       if (mytid == 0) then
        start(1)=1 ; count(1)=imt_global
        start(2)=1 ; count(2)=jmt_global
        start(3)=k ; count(3)=1
        start(4)=1 ; count(4)=1

        iret=nf_get_vara_double(ncid,   7,start,count, buffer)
        call check_err (iret)
        do j=1 ,jmt_global/2
        do i=1, imt_global
           xx = buffer(i,j)
           yy = buffer(i,jmt_global+1-j)
           buffer(i,jmt_global+1-j) = xx
           buffer(i,j) = yy
        end do
        end do

      end if

       call scatter_global(at(:,:,k,2,:), buffer,master_task, distrb_clinic, &
                          field_loc_center, field_type_scalar)
      
       if (mytid == 0) then
        start(1)=1 ; count(1)=imt_global
        start(2)=1 ; count(2)=jmt_global
        start(3)=k ; count(3)=1
        start(4)=1 ; count(4)=1

        iret1=nf_get_vara_double(ncid1,   7,start,count, buffer)
        call check_err (iret)
        do j=1 ,jmt_global/2
        do i=1, imt_global
           xx = buffer(i,j)
           yy = buffer(i,jmt_global+1-j)
           buffer(i,jmt_global+1-j) = xx
           buffer(i,j) = yy
        end do
        end do

      end if

       call scatter_global(at(:,:,k,1,:), buffer,master_task, distrb_clinic, &
                          field_loc_center, field_type_scalar)


     end do

     endif !nowork
      iret = nf_CLOSE (ncid)
      iret1 = nf_CLOSE (ncid1)
#endif
!READ_SPLIT
#endif
!HIGHRES
      if (mytid==0) then
          write(*,*)"wjl debug 3"
      end if

#ifdef U_V_INITIAL
    call POP_HaloUpdate(h0 , POP_haloClinic, POP_gridHorzLocCenter,&
            POP_fieldKindScalar, errorCode, fillValue = 0.0_r8)
    call POP_HaloUpdate(u , POP_haloClinic, POP_gridHorzLocSWcorner , &
            POP_fieldKindVector, errorCode, fillValue = 0.0_r8)
    call POP_HaloUpdate(v , POP_haloClinic, POP_gridHorzLocSWcorner , &
            POP_fieldKindVector, errorCode, fillValue = 0.0_r8)
#endif

    call POP_HaloUpdate(at(:,:,:,1,:) , POP_haloClinic, POP_gridHorzLocCenter , &
            POP_fieldKindScalar, errorCode, fillValue = 0.0_r8)
    call POP_HaloUpdate(at(:,:,:,2,:) , POP_haloClinic, POP_gridHorzLocCenter , &
            POP_fieldKindScalar, errorCode, fillValue = 0.0_r8)

      if (mytid==0) then
          write(*,*)"wjl debug 4"
      end if


!----------------------------------------------------
!   assign 0 to land grids of TSinital
!----------------------------------------------------
      DO IBLOCK= 1, NBLOCKS_CLINIC
         DO K = 1,KM
            DO J = 1,JMT
               DO I = 1,IMT
                   if(VIT(I,J,K,IBLOCK).gt.0.5) then
                   if( AT(I,J,K,1,IBLOCK).gt.50.or.AT (I,J,K,1,IBLOCK).lt.-5.0  ) &
                     write(16,*) 'AT (I,J,K,1,IBLOCK)=',i,j,k,AT (I,J,K,1,IBLOCK) 
                   if( AT(I,J,K,2,IBLOCK).gt.70.or.AT (I,J,K,2,IBLOCK).lt.0  ) &
                     write(16,*) 'AT (I,J,K,2,IBLOCK)=',i,j,k,AT (I,J,K,2,IBLOCK) 
                   !if( AT(I,J,K,1,IBLOCK).gt.50.) AT (I,J,K,1,IBLOCK) =0.0
                   !if( AT(I,J,K,1,IBLOCK).lt.-5.) AT (I,J,K,1,IBLOCK) =0.0
                   !if( AT(I,J,K,2,IBLOCK).gt.70.) AT (I,J,K,2,IBLOCK) =35.0
                   !if( AT(I,J,K,2,IBLOCK).lt.0) AT (I,J,K,2,IBLOCK) =35.0
                   endif
                  AT (I,J,K,1,IBLOCK) = AT (I,J,K,1,IBLOCK)*VIT(I,J,K,IBLOCK)
                  AT (I,J,K,2,IBLOCK) = (AT (I,J,K,2,IBLOCK)- 35.0D0)*0.001D0*VIT(I,J,K,IBLOCK)
#ifdef U_V_INITIAL
                   U(I,J,K,IBLOCK) = U(I,J,K,IBLOCK) * VIV(I,J,K,IBLOCK)
                   V(I,J,K,IBLOCK) = V(I,J,K,IBLOCK) * VIV(I,J,K,IBLOCK)
#endif
               END DO
            END DO
         END DO
      END DO
!

!ZWP20170218
!
!      DO IBLOCK= 1, NBLOCKS_CLINIC
!         DO K = 1,KM
!            DO J = 1,JMT
!               DO I = 1,IMT
!                  ATZWP (I,J,K,1,IBLOCK) = ATZWP (I,J,K,1,IBLOCK)*VIT(I,J,K,IBLOCK)
!                  ATZWP (I,J,K,2,IBLOCK) = (ATZWP (I,J,K,2,IBLOCK)- 35.0D0)*0.001D0*VIT(I,J,K,IBLOCK)
!               END DO
!            END DO
!         END DO
!      END DO
!ZWP20170218

      RESTORE = AT
!      RESTORE = ATZWP
!      deallocate(buffer)
!      deallocate(ATZWP) !ZWP20170218
!
      end if
! above codes read initial data no matter is startup or restart
      IF ( NSTART == 1) THEN

        number_day = 1

#ifndef U_V_INITIAL
         DO K = 1,KM
               ATB (:,:,K,:,:) = AT (:,:,K,:,:)
         END DO
         ATB (:,:,0,:,:) = 0.0D0
         u = 0.0d0
         v = 0.0d0
         h0= 0.0d0
#else

        CALL VINTEG (U,UB)
        CALL VINTEG (V,VB)

        DO IBLOCK = 1, NBLOCKS_CLINIC
        DO K = 1,KM
        DO J = 1,jmt
        DO I = 1,IMT
            UP (I,J,K,IBLOCK) = U (I,J,K,IBLOCK)
            VP (I,J,K,IBLOCK) = V (I,J,K,IBLOCK)
            UTF (I,J,K,IBLOCK) = U (I,J,K,IBLOCK)
            VTF (I,J,K,IBLOCK) = V (I,J,K,IBLOCK)
            ATB (I,J,K,1,IBLOCK) = AT (I,J,K,1,IBLOCK)
            ATB (I,J,K,2,IBLOCK) = AT (I,J,K,2,IBLOCK)
        END DO
        END DO
        END DO
        END DO


        DO IBLOCK = 1, NBLOCKS_CLINIC
        DO J = 1,jmt
        DO I = 1,IMT
            H0P (I,J,IBLOCK)= H0 (I,J,IBLOCK)
            UBP (I,J,IBLOCK)= UB (I,J,IBLOCK)
            VBP (I,J,IBLOCK)= VB (I,J,IBLOCK)
            H0F (I,J,IBLOCK)= H0 (I,J,IBLOCK)
            H0BF (I,J,IBLOCK)= H0 (I,J,IBLOCK)
        END DO
        END DO
        END DO
#endif

      if (mytid==0) then
          write(*,*)"wjl debug 5"
      end if
!M
#ifdef TIDEFORCE
      elevat = 0.0D0 
      h0_tide = 0.0D0 
#endif 
#ifdef COUP
            do iblock = 1, nblocks_clinic
            do j=1,jmt
            do i=1,imt
               t_cpl (i,j,iblock)  = 273.15+at(i,j,1,1,iblock)
               s_cpl (i,j,iblock)  = at(i,j,1,2,iblock)*1000.+35.
               qheat     (i,j,iblock)  = 0.0
               u_cpl (i,j,iblock)  = 0.0
               v_cpl (i,j,iblock)  = 0.0
               dhdx  (i,j,iblock)  = 0.0
               dhdy  (i,j,iblock)  = 0.0
            end do
            end do
            end do
#endif
      ELSE
 
!     ------------------------------------------------------------------
!     READ INTERMEDIATE RESULTS (fort.22/fort.21)
!     ------------------------------------------------------------------
 

!xk add 20240626
#ifdef SPLIT_RPT
!      io_size = nx_block * ny_block * (4 * km  + 1)* 8
!      call c_fopen_read("ocn"//CHAR(0), "rpt.dir/rpt"//CHAR(0), io_size)
       io_size = nx_block * ny_block * km * 8
       call c_fopen_restart_read("restart"//CHAR(0), "at1_temp.dir/at1_temp"//CHAR(0), io_size)
       call c_fread(at(1,1,1,1,1), io_size)
       call c_fclose()

       call c_fopen_restart_read("restart"//CHAR(0), "at2_salt.dir/at2_salt"//CHAR(0), io_size)
       call c_fread(at(1,1,1,2,1), io_size)
       call c_fclose()

       call c_fopen_restart_read("restart"//CHAR(0), "u.dir/u"//CHAR(0), io_size)
       call c_fread(u(1,1,1,1), io_size)
       call c_fclose()

       call c_fopen_restart_read("restart"//CHAR(0), "v.dir/v"//CHAR(0), io_size)
       call c_fread(v(1,1,1,1), io_size)
       call c_fclose()

       io_size = nx_block * ny_block * 12 * 8
       call c_fopen_restart_read("restart"//CHAR(0), "other_data.dir/other_data"//CHAR(0), io_size)
       io_size = nx_block * ny_block  * 8
       call c_fread(h0(1,1,1), io_size)
       call c_fclose()


#ifdef COUP
       io_size = nx_block * ny_block * 7 * 8
       call c_fopen_restart_read("restart"//CHAR(0), "coup_data.dir/coup_data"//CHAR(0), io_size)
       io_size = nx_block * ny_block  * 8
       call c_fread(t_cpl(1,1,1), io_size)
       call c_fread(s_cpl(1,1,1), io_size)
       call c_fread(u_cpl(1,1,1), io_size)
       call c_fread(v_cpl(1,1,1), io_size)
       call c_fread(dhdx(1,1,1), io_size)
       call c_fread(dhdy(1,1,1), io_size)
       call c_fread(qheat(1,1,1), io_size)
       call c_fclose()
#endif

!      io_size = 4 * 2
!      call c_fopen_read("restart"//CHAR(0), "date.dir/date"//CHAR(0), io_size)
!
!      io_size = 4
!      call c_fread(r_month, io_size)
!      call c_fread(r_day, io_size)
!      call c_fclose()

      if (mytid==0) then
          open (17,file='rpointer.ocn',form='formatted')
          read(17,'(a18)') fname
          close(17)
          read(fname(9:12),"(i4)") rst_year
          read(fname(14:15),"(i2)") number_month
          read(fname(17:18),"(i2)") number_day
          write(*,*) "rpointer year month day: ",rst_year,number_month,number_day
       end if
#ifdef TIDEFORCE
       elevat = 0.0D0 
       h0_tide = 0.0D0 
#endif      
!      number_month = 2
!      number_day = 1
      if(nstart==3) then !for branch run !LPF20170621
            number_month= yearadd*12+number_month 
            number_day= dayadd+number_day
      endif
      month= number_month

       if(mytid == 0)     write(*,*) 'number_month =',number_month,'mon0=',mon0,&
                       'number_day=',number_day,'iday=',iday
       call POP_HaloUpdate(h0 , POP_haloClinic, POP_gridHorzLocCenter,&
                       POP_fieldKindScalar, errorCode, fillValue = 0.0_r8)
      call POP_HaloUpdate(u , POP_haloClinic, POP_gridHorzLocSWcorner , &
                  POP_fieldKindVector, errorCode, fillValue = 0.0_r8)
      call POP_HaloUpdate(v , POP_haloClinic, POP_gridHorzLocSWcorner , &
                       POP_fieldKindVector, errorCode, fillValue = 0.0_r8)
      call POP_HaloUpdate(at(:,:,:,1,:) , POP_haloClinic, POP_gridHorzLocCenter , &
                  POP_fieldKindScalar, errorCode, fillValue = 0.0_r8)
      call POP_HaloUpdate(at(:,:,:,2,:) , POP_haloClinic, POP_gridHorzLocCenter , &
                       POP_fieldKindScalar, errorCode, fillValue = 0.0_r8)

      call mpi_bcast(month,1,mpi_integer,0,mpi_comm_ocn,ierr)
      call mpi_bcast(rst_year,1,mpi_integer,0,mpi_comm_ocn,ierr)
      call mpi_bcast(number_month,1,mpi_integer,0,mpi_comm_ocn,ierr)
      call mpi_bcast(number_day,1,mpi_integer,0,mpi_comm_ocn,ierr)

#else
       ! wjl 20240626
        if (mytid==0) then
          open (17,file='rpointer.ocn',form='formatted')
          read(17,'(a18)') fname
          close(17)
          open(22,file=trim(out_dir)//fname,form='unformatted')
       end if
      if (mytid==0) then
          write(*,*)"wjl debug 42"
      end if
!
         allocate (buffer(imt_global,jmt_global))
!
         if (mytid==0) then
         READ (22)buffer
         end if

       !TODO 
!        call scatter_global(h0,buffer(2:imt_global+1,:), master_task, distrb_clinic, &
         call scatter_global(h0,buffer, master_task, distrb_clinic, &
                          field_loc_center, field_type_scalar)
         call POP_HaloUpdate(h0 , POP_haloClinic, POP_gridHorzLocCenter,&
                       POP_fieldKindScalar, errorCode, fillValue = 0.0_r8)
!
         do k=1,km
         if (mytid==0) then
         READ (22)buffer
         end if
!        call scatter_global(u(:,:,k,:), buffer(2:imt_global+1,:), master_task, distrb_clinic, &
         call scatter_global(u(:,:,k,:), buffer, master_task, distrb_clinic, &
                          field_loc_swcorner, field_type_vector)
         end do
         call POP_HaloUpdate(u , POP_haloClinic, POP_gridHorzLocSWcorner , &
                       POP_fieldKindVector, errorCode, fillValue = 0.0_r8)
!
         do k=1,km
         if (mytid==0) then
         READ (22)buffer
         end if
!        call scatter_global(v(:,:,k,:),buffer(2:imt_global+1,:), master_task, distrb_clinic, &
         call scatter_global(v(:,:,k,:),buffer, master_task, distrb_clinic, &
                          field_loc_swcorner, field_type_vector)
         end do
         call POP_HaloUpdate(v , POP_haloClinic, POP_gridHorzLocSWcorner , &
                       POP_fieldKindVector, errorCode, fillValue = 0.0_r8)
!
         do k=1,km
         if (mytid==0) then
         READ (22)buffer
         end if
!        call scatter_global(at(:,:,k,1,:),buffer(2:imt_global+1,:), master_task, distrb_clinic, &
         call scatter_global(at(:,:,k,1,:),buffer, master_task, distrb_clinic, &
                          field_loc_center, field_type_scalar)
         end do
         call POP_HaloUpdate(at(:,:,:,1,:) , POP_haloClinic, POP_gridHorzLocCenter , &
                       POP_fieldKindScalar, errorCode, fillValue = 0.0_r8)
!
         do k=1,km
         if (mytid==0) then
         READ (22)buffer
         end if
!        call scatter_global(at(:,:,k,2,:),buffer(2:imt_global+1,:), master_task, distrb_clinic, &
         call scatter_global(at(:,:,k,2,:),buffer, master_task, distrb_clinic, &
                          field_loc_center, field_type_scalar)
         end do
         call POP_HaloUpdate(at(:,:,:,2,:) , POP_haloClinic, POP_gridHorzLocCenter , &
                       POP_fieldKindScalar, errorCode, fillValue = 0.0_r8)
!lhl20110728 for ws
         do k=1,km
         if (mytid==0) READ (22)buffer
         end do
         if (mytid==0) READ (22)buffer !su
         if (mytid==0) READ (22)buffer !sv
         if (mytid==0) READ (22)buffer !swv
         if (mytid==0) READ (22)buffer !sshf
         if (mytid==0) READ (22)buffer !lthf
         if (mytid==0) READ (22)buffer  !fresh
        ! if (mytid==0) READ (22)buffer
        !if (mytid==0) READ (22)buffer !add lwv !noneed for new
        ! if (mytid==0) READ (22)buffer
        ! if (mytid==0) READ (22)buffer
!lhl20110728
         if (mytid == 0) then
          read(22)number_month,number_day
           if(nstart==3) then !for branch run !LPF20170621
             number_month= yearadd*12+number_month 
             number_day= dayadd+number_day
           endif
           month= number_month
         endif
      call mpi_bcast(month,1,mpi_integer,0,mpi_comm_ocn,ierr)
      call mpi_bcast(number_month,1,mpi_integer,0,mpi_comm_ocn,ierr)
      call mpi_bcast(number_day,1,mpi_integer,0,mpi_comm_ocn,ierr)
            write(*,*) 'number_month =',number_month,'mon0=',mon0,&
                       'number_day=',number_day,'iday=',iday
#ifdef COUP
         if (nstart==2) then
            month=(cdate/10000-1)*12+mod(cdate,10000)/100
!M
            do iblock = 1, nblocks_clinic
            do j=1,jmt
            do i=1,imt
               ! lihuimin, TODO, to be considered
               ! lihuimin, 2012.7.23, coordinate with flux_cpl, ft. yu
               t_cpl (i,j,iblock)  = 273.15+at(i,j,1,1,iblock)
               s_cpl (i,j,iblock)  = at(i,j,1,2,iblock)*1000.+35.
               ! modi end
               qheat     (i,j,iblock)  = 0.0
               u_cpl (i,j,iblock)  = 0.0
               v_cpl (i,j,iblock)  = 0.0
               dhdx  (i,j,iblock)  = 0.0
               dhdy  (i,j,iblock)  = 0.0
            end do
            end do
            end do
         else
!LPF 20120815
!for t_cpl
          if (mytid==0) then
           READ (22)buffer
          end if
!         call scatter_global(t_cpl,buffer(2:imt_global+1,:), master_task, distrb_clinic, &
          call scatter_global(t_cpl,buffer, master_task, distrb_clinic, &
                          field_loc_center, field_type_scalar)
!for s_cpl
          if (mytid==0) then
           READ (22)buffer
          end if
!         call scatter_global(s_cpl, buffer(2:imt_global+1,:), master_task, distrb_clinic, &
          call scatter_global(s_cpl, buffer, master_task, distrb_clinic, &
                          field_loc_center, field_type_scalar)
!for u_cpl
          if (mytid==0) then
           READ (22)buffer
          end if
!for v_cpl
!         call scatter_global(u_cpl, buffer(2:imt_global+1,:), master_task, distrb_clinic, &
          call scatter_global(u_cpl, buffer, master_task, distrb_clinic, &
                          field_loc_center, field_type_vector)
          if (mytid==0) then
           READ (22)buffer
          end if
!         call scatter_global(v_cpl, buffer(2:imt_global+1,:), master_task, distrb_clinic, &
          call scatter_global(v_cpl, buffer, master_task, distrb_clinic, &
                          field_loc_center, field_type_vector)

!for dhdx 
          if (mytid==0) then
           READ (22)buffer
          end if
!         call scatter_global(dhdx,buffer(2:imt_global+1,:), master_task, distrb_clinic, &
          call scatter_global(dhdx,buffer, master_task, distrb_clinic, &
                          field_loc_center, field_type_scalar)
!for dhdy 
          if (mytid==0) then
           READ (22)buffer
          end if
!         call scatter_global(dhdy, buffer(2:imt_global+1,:), master_task, distrb_clinic, &
          call scatter_global(dhdy, buffer, master_task, distrb_clinic, &
                          field_loc_center, field_type_scalar)
!for q
          if (mytid==0) then
           READ (22)buffer
          end if
!         call scatter_global(q ,buffer(2:imt_global+1,:),  master_task, distrb_clinic, &
          call scatter_global(qheat ,buffer,  master_task, distrb_clinic, &
                          field_loc_center, field_type_scalar)
!
         end if
!LPF 20120815
!Yu

#endif
          deallocate(buffer)
          if (mytid == 0) CLOSE(22)
!xk add 20240626
#endif

       if (simple_assm.and.nstart==1) then !LPF20170821 
             number_month= yearadd*12+number_month
             number_day= dayadd+number_day
       endif
 
         NMFF = MOD (MONTH -1,12)
 
         CALL VINTEG (U,UB)
         CALL VINTEG (V,VB)
 
      DO IBLOCK = 1, NBLOCKS_CLINIC
         DO K = 1,KM
            DO J = 1,jmt
               DO I = 1,IMT
                  UP (I,J,K,IBLOCK) = U (I,J,K,IBLOCK)
                  VP (I,J,K,IBLOCK) = V (I,J,K,IBLOCK)
                  UTF (I,J,K,IBLOCK) = U (I,J,K,IBLOCK)
                  VTF (I,J,K,IBLOCK) = V (I,J,K,IBLOCK)
                  ATB (I,J,K,1,IBLOCK) = AT (I,J,K,1,IBLOCK)
                  ATB (I,J,K,2,IBLOCK) = AT (I,J,K,2,IBLOCK)
               END DO
            END DO
         END DO
       END DO
 
 
      DO IBLOCK = 1, NBLOCKS_CLINIC
         DO J = 1,jmt
            DO I = 1,IMT
               H0P (I,J,IBLOCK)= H0 (I,J,IBLOCK)
               UBP (I,J,IBLOCK)= UB (I,J,IBLOCK)
               VBP (I,J,IBLOCK)= VB (I,J,IBLOCK)
               H0F (I,J,IBLOCK)= H0 (I,J,IBLOCK)
               H0BF (I,J,IBLOCK)= H0 (I,J,IBLOCK)
            END DO
         END DO
      END DO
      END IF
!

#ifdef BIHAR
      call init_del4t
      call init_del4u
#else
      call init_del2t
      call init_del2u
#endif
!
#if (defined LOWRES)
      if (diag_mth .or. diag_msf)  then
          call init_diagnostics
      end if
#endif
!
      if (mytid==0)then
          write(6,*)"END-----------INIRUN !"
#ifdef COUP
          call licom_shr_flush(6)
#endif
      endif 

      RETURN

      END SUBROUTINE INIRUN
 

