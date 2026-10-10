!  CVS: $Id: ssave-cdf.F90,v 1.1.1.1 2004/04/29 06:22:39 lhl Exp $
!     =================
      subroutine SSAVEINSSTEP(IIS)
#include <def-undef.h>
!     =================
!     output in NETcdf format
!     written by liu hai long 2001 jun
use param_mod
use pconst_mod
use dyn_mod
use tracer_mod
use forc_mod
use work_mod
use buf_mod, only:t_cpl,s_cpl,u_cpl,v_cpl,dhdx,dhdy,qheat      
use domain
use gather_scatter
use distribution
use msg_mod
use blocks
!use mct_mod
!use esmf
!use seq_flds_mod
!use seq_cdata_mod
!use seq_infodata_mod
!use seq_timemgr_mod

  implicit none


      logical       :: write_restart     ! restart now
      character (len=22) :: fname
      integer :: klevel, iiday, nwmf,hour,last
      integer :: IIS 
      integer(kind(1))   :: curr_ymd     ! Current date YYYYMMDD
      integer :: wid3,wid4,wid5,wid6,wid7 ! mpi id for parallel write
      integer * 8 :: c_nx, c_ny, c_km, c_month, c_day
! 
       if(nproc > 50) then
        wid3=50
        wid4=20
        wid5=30
        wid6=40
       else
        wid3=12
        wid4=13
        wid5=14
        wid6=15
       endif
       if(nproc > 600) then
        wid3=300
        wid4=400
        wid5=500
        wid6=600
       endif

#ifdef SPLIT_WRITE


  nwmf= iyfm
  !DO LPF 20200401
  !mon0=number_month
!  mon0=month-(iyfm-1)*12
  !DO LPF 20200401
  !mon0 = mod(number_month-1,12) + 1
  !mon0=number_month
  iiday = nnmonth(mon0) + iday + 1+ (nwmf-1)*365

  !
  write_restart = .false.
  if ( rest_freq < 0) then
  if (number_day == 1) write_restart = .true.
  else if ( mod(iiday-1,rest_freq) == 0  ) then
  write_restart = .true.
  end if
  if(mytid==0) write(*,*)'in instant,iday=,imd=,mon0=, nwmf=, rest_freq= ',iday,imd,mon0,nwmf,rest_freq
  if(mytid==0) write(*,*)'in instant,write_restart', write_restart
!

  where (vit(:,:,1,:) < 0.5D0) h0 = 0.0_r8

  where (viv < 0.5D0) u = 0.0_r8
  where (viv < 0.5D0) v = 0.0_r8
  where (vit < 0.5D0) at(:,:,:,1,:) = 0.0_r8
  where (vit < 0.5D0) at(:,:,:,2,:) = 0.0_r8
  do klevel=1,km
     where (vit(:,:,klevel,:) < 0.5D0) ws(:,:,klevel,:) = 0.0_r8
  end do

!    do j=1,jmt
!    do i=1,imt
!    if(ws(i,j,1,1) .eq. 1.0D0) write(*,*) "ws fort is 1, i,j",i,j,mytid
!    end do
!    end do

  ! write(fname(1:2),'(i2.2)')mon0
  ! write(fname(3:6),'(i4.4)')nwmf
  ! write(fname(7:9),'(i2.2)') number_day
  ! write(fname(10:14),'(i5.5)') IIS
  if (nstep > 0) then
    hour = IIS / nstep * (NHIS/60)
  else
    hour = 0
  end if
  last = 0
  fname(1:8)='fort.22.'
  fname(13:13)='-'
  fname(16:16)='-'
  write(fname(14:15),'(i2.2)')mon0
  write(fname(9:12),'(i4.4)')nwmf
  write(fname(17:18),'(i2.2)') dd_licom
  write(fname(19:20),'(i2.2)') hour
  write(fname(21:22),'(i2.2)') last

  if(mytid==0) write(*,*)'in instant,write_restart,hour ',fname,nwmf,mon0,dd_licom,hour

  c_nx = imt 
  c_ny = jmt
  c_km = km
  call write_init(mpi_comm_ocn,c_nx, c_ny, c_km)

  call write_km_data(u(1,1,1,1), trim(fname)//CHAR(0), "u"//CHAR(0))
  call write_km_data(v(1,1,1,1), trim(fname)//CHAR(0), "v"//CHAR(0))
  call write_km_data(at(1,1,1,1,1), trim(fname)//CHAR(0), "at1_temp"//CHAR(0))
  call write_km_data(at(1,1,1,2,1), trim(fname)//CHAR(0), "at2_salt"//CHAR(0))
  call write_km_data(PDENSITY(1,1,1,1), trim(fname)//CHAR(0), "dens"//CHAR(0))
  call write_kmp1_data(ws(1,1,1,1), trim(fname)//CHAR(0), "ws"//CHAR(0))
!  if(mytid==0) write(*,*)'write u v at ws done'

  where (vit(:,:,1,:) < 0.5D0) h0 = 0.0_r8
  where (viv(:,:,1,:) < 0.5D0) su = 0.0_r8
  where (viv(:,:,1,:) < 0.5D0) sv = 0.0_r8
  where (vit(:,:,1,:) < 0.5D0) swv = 0.0_r8
  where (vit(:,:,1,:) < 0.5D0) lwv = 0.0_r8
  where (vit(:,:,1,:) < 0.5D0) sshf = 0.0_r8
  where (vit(:,:,1,:) < 0.5D0) lthf = 0.0_r8
  where (vit(:,:,1,:) < 0.5D0) fresh = 0.0_r8
  where (vit(:,:,1,:) < 0.5D0) seaice = 0.0_r8

  call write_other_data(h0(1,1,1), su(1,1,1), sv(1,1,1), swv(1,1,1),lwv(1,1,1), sshf(1,1,1), &
          lthf(1,1,1), fresh(1,1,1),seaice(1,1,1),runoff(1,1,1),rain(1,1,1),psa(1,1,1), trim(fname)//CHAR(0), "other_data"//CHAR(0))
!  if(mytid==0) write(*,*)'write h0 su sv swv sshf lthf fresh done'

#ifdef COUP
  where (vit(:,:,1,:) < 0.5D0) t_cpl = 0.0_r8
  where (viv(:,:,1,:) < 0.5D0) s_cpl = 0.0_r8
  where (viv(:,:,1,:) < 0.5D0) u_cpl = 0.0_r8
  where (vit(:,:,1,:) < 0.5D0) v_cpl = 0.0_r8
  where (vit(:,:,1,:) < 0.5D0) dhdx = 0.0_r8
  where (vit(:,:,1,:) < 0.5D0) dhdy = 0.0_r8
  where (vit(:,:,1,:) < 0.5D0) qheat = 0.0_r8
   call write_coup_data(t_cpl(1,1,1), s_cpl(1,1,1), u_cpl(1,1,1), v_cpl(1,1,1), &
          dhdx(1,1,1), dhdy(1,1,1),qheat(1,1,1), trim(fname)//CHAR(0), "coup_data"//CHAR(0))

#endif
  c_month = month
  c_day = number_day

  call write_date_data(c_month, c_day, trim(fname)//CHAR(0), "date"//CHAR(0))



#else


!	if(mytid==master_task .or. (mytid==master_task+wid3) .or.(mytid==master_task+wid4) ) then
!        if (.not. allocated(buffer3d))allocate(buffer3d(imt,jmt,km,nblocks_tot))
!        endif
!        if((mytid==master_task+wid5) .or. (mytid==master_task+wid6) ) then
!        if (.not. allocated(buffer3d)) allocate(buffer3d(imt,jmt,km,nblocks_tot))
!        endif
	if(mytid==master_task .or. (mytid==master_task+wid3) .or.(mytid==master_task+wid4) ) then
        if (.not. allocated(buffer2d))allocate(buffer2d(imt,jmt,nblocks_tot))
        endif
        if((mytid==master_task+wid5) .or. (mytid==master_task+wid6) ) then
        if (.not. allocated(buffer2d)) allocate(buffer2d(imt,jmt,nblocks_tot))
        endif

         number_day = iday +1 !LPF 20120816
         !number_month = month
      !if (mod(number_month,12)==0) then
         if ( number_day  > imd ) then
           if (mod(number_month,12)==0) then
            !DO LPF 20200401
            iyfm=int(number_month/12)+1
            !iyfm=iyfm+1 !20200401
            !number_month=1
            number_month=number_month+1
            !DO LPF 20200401 
            number_day=1
           else
             number_day = 1
             number_month = number_month + 1
           !if (mod(number_month,12)==0) then
           ! iyfm=iyfm+1
           end if
         end if
         nwmf= iyfm
         !DO LPF 20200401
         !mon0=number_month
         mon0=number_month-(iyfm-1)*12
         !DO LPF 20200401
         !mon0 = mod(number_month-1,12) + 1
         !mon0=number_month
         iiday = nnmonth(mon0) + iday + 1+ (nwmf-1)*365

!
     write_restart = .false.
     if ( rest_freq < 0) then
        if (number_day == 1) write_restart = .true.
     else if ( mod(iiday-1,rest_freq) == 0  ) then
        write_restart = .true.
     end if
     if(mytid==0) write(*,*)'in instant,iday=,imd=,mon0=, nwmf=, rest_freq= ',iday,imd,mon0,nwmf,rest_freq
     if(mytid==0) write(*,*)'in instant,write_restart', write_restart
!
!    if ( mod(iday,rest_freq) == 0 .or. iday == imd .or. iday == 10 .or. iday ==20) then
!    if ( mod(iday,rest_freq) == 0 .or. iday == 1 ) then
!    if ( mod(iiday,rest_freq) == 1  ) then
     if ( write_restart  ) then
!
         fname(1:8)='fort.22.'
         fname(13:13)='-'
         fname(16:16)='-'
         write(fname(14:15),'(i2.2)')mon0
         write(fname(9:12),'(i4.4)')nwmf
         write(fname(17:18),'(i2.2)') number_day
         write(fname(19:22),'(i5.5)') IIS
 
!        if ( number_day ==1 .and. mon0 < 12) then
!            write(fname(14:15),'(i2.2)')mon0+1
!        end if
!
!        if ( number_day ==1 .and. mon0 == 12) then
!            write(fname(14:15),'(i2.2)')mon0-11
!            write(fname(9:12),'(i4.4)')nwmf+1
!        end if

       if ( mytid == 0 ) then
!         open (17, file="rpointer.ocn", form='formatted')
!         write(17,'(a18)') fname
!         close(17)
         write(*,*)'fname=',fname
         open(22,file=trim(out_dir)//fname,form='unformatted')
       end if

       if ( mytid == wid3 ) then
         fname(1:8)='fort.23.'
         open(23,file=trim(out_dir)//fname,form='unformatted')
       end if
       if ( mytid == wid4 ) then
         fname(1:8)='fort.24.'
         open(24,file=trim(out_dir)//fname,form='unformatted')
       end if
       if ( mytid == wid5 ) then
         fname(1:8)='fort.25.'
         open(25,file=trim(out_dir)//fname,form='unformatted')
       end if
       if ( mytid == wid6 ) then
         fname(1:8)='fort.26.'
         open(26,file=trim(out_dir)//fname,form='unformatted')
       end if
       if ( mytid == 0 ) then
         fname(1:8)='fort.27.'
         open(27,file=trim(out_dir)//fname,form='unformatted')
       end if
!
       allocate(buffer(imt_global,jmt_global))
       where (vit(:,:,1,:) < 0.5D0) h0 = 0.0_r8
       call gather_global(buffer, h0, master_task,distrb_clinic)
!
       if (mytid==0) then
          WRITE (22)buffer
       end if
!
      where (viv < 0.5D0) u = 0.0_r8
      where (viv < 0.5D0) v = 0.0_r8
      where (vit < 0.5D0) at(:,:,:,1,:) = 0.0_r8
      where (vit < 0.5D0) at(:,:,:,2,:) = 0.0_r8
      do klevel=1,km
         where (vit(:,:,klevel,:) < 0.5D0) ws(:,:,klevel,:) = 0.0_r8
      end do

!start to write 3D file, fort.22 - 26
      do klevel=1,km
      call MPI_gather(u (:,:,klevel,1),  imt*jmt,MPI_PR,buffer2d,imt*jmt,MPI_PR,master_task,mpi_comm_ocn,ierr)
      call MPI_gather(v (:,:,klevel,1),  imt*jmt,MPI_PR,buffer2d,imt*jmt,MPI_PR,master_task+wid3,mpi_comm_ocn,ierr)
      call MPI_gather(at(:,:,klevel,1,1),imt*jmt,MPI_PR,buffer2d,imt*jmt,MPI_PR,master_task+wid4,mpi_comm_ocn,ierr)
      call MPI_gather(at(:,:,klevel,2,1),imt*jmt,MPI_PR,buffer2d,imt*jmt,MPI_PR,master_task+wid5,mpi_comm_ocn,ierr)
      call MPI_gather(ws(:,:,klevel,1),  imt*jmt,MPI_PR,buffer2d,imt*jmt,MPI_PR,master_task+wid6,mpi_comm_ocn,ierr)

      if (mytid==0) then
             call gather_global2(buffer, buffer2d, master_task,distrb_clinic)
             WRITE (22) buffer
      end if

      if (mytid==wid3) then
             call gather_global2(buffer, buffer2d, master_task,distrb_clinic)
             WRITE (23) buffer
      end if

      if (mytid==wid4) then
             call gather_global2(buffer, buffer2d, master_task,distrb_clinic)
             WRITE (24) buffer
      end if

      if (mytid==wid5) then
             call gather_global2(buffer, buffer2d, master_task,distrb_clinic)
             WRITE (25) buffer
      end if

      if (mytid==wid6) then
             call gather_global2(buffer, buffer2d, master_task,distrb_clinic)
             WRITE (26) buffer
      end if
      end do

!start to write fort.27

      where (viv(:,:,1,:) < 0.5D0) su = 0.0_r8
      call gather_global(buffer, su, master_task,distrb_clinic)
      if (mytid==0) then
         WRITE (27)buffer
      end if
!         write(*,*)'finish su'

      where (viv(:,:,1,:) < 0.5D0) sv = 0.0_r8
      call gather_global(buffer, sv, master_task,distrb_clinic)
      if (mytid==0) then
         WRITE (27)buffer
      end if

!         write(*,*)'finish sv'
      where (vit(:,:,1,:) < 0.5D0) swv = 0.0_r8
     call gather_global(buffer, swv, master_task,distrb_clinic)
     if (mytid==0) then
         WRITE (27)buffer
     end if
!         write(*,*)'finish lwv'
      where (vit(:,:,1,:) < 0.5D0) sshf = 0.0_r8
     call gather_global(buffer, sshf, master_task,distrb_clinic)
     if (mytid==0) then
         WRITE (27)buffer
     end if
!         write(*,*)'finish sshf'
      where (vit(:,:,1,:) < 0.5D0) lthf = 0.0_r8
     call gather_global(buffer, lthf, master_task,distrb_clinic)
     if (mytid==0) then
         WRITE (27)buffer
     end if
!         write(*,*)'finish lthf'
      where (vit(:,:,1,:) < 0.5D0) fresh = 0.0_r8
     call gather_global(buffer, fresh, master_task,distrb_clinic)
     if (mytid==0) then
         WRITE (27)buffer
     end if
!         write(*,*)'finish fresh'
     if (mytid==0) then
         write(27) number_month, number_day
     end if
!
#ifdef COUP
      where (vit(:,:,1,:) < 0.5D0) t_cpl = 0.0_r8
     call gather_global(buffer, t_cpl, master_task,distrb_clinic)
     if (mytid==0) then
         WRITE (22)buffer
     end if
!         write(*,*)'finish t_cpl'
      where (vit(:,:,1,:) < 0.5D0) s_cpl = 0.0_r8
     call gather_global(buffer, s_cpl, master_task,distrb_clinic)
     if (mytid==0) then
         WRITE (22)buffer
     end if
!         write(*,*)'finish s_cpl'
      where (vit(:,:,1,:) < 0.5D0) u_cpl = 0.0_r8
     call gather_global(buffer, u_cpl, master_task,distrb_clinic)
     if (mytid==0) then
         WRITE (22)buffer
     end if
!         write(*,*)'finish u_cpl'
      where (vit(:,:,1,:) < 0.5D0) v_cpl = 0.0_r8
     call gather_global(buffer, v_cpl, master_task,distrb_clinic)
     if (mytid==0) then
         WRITE (22)buffer
     end if
!         write(*,*)'finish v_cpl'
      where (vit(:,:,1,:) < 0.5D0) dhdx = 0.0_r8
     call gather_global(buffer, dhdx, master_task,distrb_clinic)
     if (mytid==0) then
         WRITE (22)buffer
     end if
!         write(*,*)'finish dhdx'
      where (vit(:,:,1,:) < 0.5D0) dhdy = 0.0_r8
     call gather_global(buffer, dhdy, master_task,distrb_clinic)
     if (mytid==0) then
         WRITE (22)buffer
     end if
!         write(*,*)'finish dhdy'
      where (vit(:,:,1,:) < 0.5D0) qheat = 0.0_r8
     call gather_global(buffer, qheat, master_task,distrb_clinic)
     if (mytid==0) then
         WRITE (22)buffer
     end if
!         write(*,*)'finish q'
#endif

     deallocate( buffer)
      if(mytid==0) then
        close(22)
        close(27)
      end if
      if(mytid==wid3) then
        close(23)
      end if
      if(mytid==wid4) then
        close(24)
      end if
      if(mytid==wid5) then
        close(25)
      end if
      if(mytid==wid6) then
        close(26)
      end if
!
  end if
#endif
      return
      end

!#if (defined NETCDF) || (defined ALL)
!      SUBROUTINE check_err (iret)
!#include <netcdf.inc>
!      INTEGER :: iret
!      IF (iret /= NF_NOERR) THEN
!         PRINT *, nf_strerror (iret)
!         STOP
!      END IF
!      END SUBROUTINE check_err
!#endif

