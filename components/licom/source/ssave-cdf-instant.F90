!  CVS: $Id: ssave-cdf.F90,v 1.1.1.1 2004/04/29 06:22:39 lhl Exp $
!     =================
      subroutine SSAVEINS
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
use blocks
!use mct_mod
!use esmf
!use seq_flds_mod
!use seq_cdata_mod
!use seq_infodata_mod
!use seq_timemgr_mod
use msg_mod
use MPIIO
use module_add
use blocks
  implicit none


      logical       :: write_restart     ! restart now
      character (len=18) :: fname
      integer :: klevel, iiday, nwmf
      integer(kind(1))   :: curr_ymd     ! Current date YYYYMMDD
      integer::fh,Byte_size_time,date_count
      integer*8 all_offset
      type (block) :: this_block  ! block info for current block
	  integer :: global_sizes(2),local_sizes(2),start_point(2),subarray_type
	   
!
         number_day = iday !+1 !LPF 20120816
         number_month = month
         if ( number_day  > imd ) then
             number_day = 1
             number_month = month + 1
         end if
         nwmf= iyfm
         iiday = nnmonth(mon0) + iday + (nwmf-1)*365
!
!     write_restart = .true.
     if ( rest_freq < 0) then
        if (number_day == 1) write_restart = .true.
     else if ( mod(totday,rest_freq) == 0  ) then
        write_restart = .true.
     end if
     if(mytid==0) write(*,*)'in instant,iday=,imd=,mon0=, nwmf=, rest_freq= ',iday,imd,mon0,nwmf,rest_freq
     if(mytid==0) write(*,*)'in instant,write_restart:', totday, rest_freq, write_restart
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
 
!        if ( number_day ==1 .and. mon0 < 12) then
!            write(fname(14:15),'(i2.2)')mon0+1
!        end if
!
!        if ( number_day ==1 .and. mon0 == 12) then
!            write(fname(14:15),'(i2.2)')mon0-11
!            write(fname(9:12),'(i4.4)')nwmf+1
!        end if
    if ( mytid == 0 ) then
         open (17, file="rpointer.ocn", form='formatted')  
         write(17,'(a18)') fname
         close(17)
         write(*,*)'fname=',fname
    end if	

all_offset=0
	
       call mpi_file_open(mpi_comm_ocn,fname,mpi_mode_wronly+mpi_mode_create,mpi_info_null,fh,ierror)    

	do n=1,nblocks_tot
	if (distrb_clinic%proc(n) == my_task+1) then
        this_block = get_block(n,n)
        date_count=(this_block%ie-this_block%ib+1)*(this_block%je-this_block%jb+1)
		
		global_sizes(1)=imt_global
		global_sizes(2)=jmt_global
		local_sizes(1)=this_block%ie-this_block%ib+1
		local_sizes(2)=this_block%je-this_block%jb+1
		start_point(1)=this_block%i_glob(this_block%ib)-1 !i=3 to G_i=1
		start_point(2)=this_block%j_glob(this_block%jb)-1
		call MPI_Type_create_subarray(2,global_sizes,local_sizes,start_point, &
		MPI_ORDER_FORTRAN, MPI_dbl, subarray_type,ierror)
		call MPI_Type_commit(subarray_type,ierror)
		exit
	end if
	end do
	
		where (vit(:,:,1,:) < 0.5D0) h0 = 0.0_r8 
        call mpiio_write(fh,h0,this_block,distrb_clinic,subarray_type,date_count,all_offset)

		 
     where (viv < 0.5D0) u = 0.0_r8
      do klevel=1,km
		call mpiio_write(fh,u(:,:,klevel,:),this_block,distrb_clinic,subarray_type,date_count,all_offset)		  
       end do

 !        write(*,*)'finish U'
      where (viv < 0.5D0) v = 0.0_r8
      do klevel=1,km
		  call mpiio_write(fh,v(:,:,klevel,:),this_block,distrb_clinic,subarray_type,date_count,all_offset)
       end do
 !        write(*,*)'finish V'
      where (vit < 0.5D0) at(:,:,:,1,:) = 0.0_r8
      do klevel=1,km
		 call mpiio_write(fh,at(:,:,klevel,1,:),this_block,distrb_clinic,subarray_type,date_count,all_offset)
      end do
 !        write(*,*)'finish at1'
      where (vit < 0.5D0) at(:,:,:,2,:) = 0.0_r8
      do klevel=1,km
		call mpiio_write(fh,at(:,:,klevel,2,:),this_block,distrb_clinic,subarray_type,date_count,all_offset)		 
      end do
!         write(*,*)'finish at2'
!lhl20110728
      do klevel=1,km
         where (vit(:,:,klevel,:) < 0.5D0) ws(:,:,klevel,:) = 0.0_r8
		call mpiio_write(fh,ws(:,:,klevel,:),this_block,distrb_clinic,subarray_type,date_count,all_offset)		 
      end do


      where (viv(:,:,1,:) < 0.5D0) su = 0.0_r8
         call mpiio_write(fh,su,this_block,distrb_clinic,subarray_type,date_count,all_offset)
!         write(*,*)'finish su'
      where (viv(:,:,1,:) < 0.5D0) sv = 0.0_r8
         call mpiio_write(fh,sv,this_block,distrb_clinic,subarray_type,date_count,all_offset)
!         write(*,*)'finish sv'
      where (vit(:,:,1,:) < 0.5D0) swv = 0.0_r8
         call mpiio_write(fh,swv,this_block,distrb_clinic,subarray_type,date_count,all_offset)
 !        write(*,*)'finish swv'
      where (vit(:,:,1,:) < 0.5D0) lwv = 0.0_r8
         call mpiio_write(fh,lwv,this_block,distrb_clinic,subarray_type,date_count,all_offset)
 !        write(*,*)'finish lwv'
      where (vit(:,:,1,:) < 0.5D0) sshf = 0.0_r8
         call mpiio_write(fh,sshf,this_block,distrb_clinic,subarray_type,date_count,all_offset)
 !        write(*,*)'finish sshf'
      where (vit(:,:,1,:) < 0.5D0) lthf = 0.0_r8
         call mpiio_write(fh,lthf,this_block,distrb_clinic,subarray_type,date_count,all_offset)
 !        write(*,*)'finish lthf'
      where (vit(:,:,1,:) < 0.5D0) fresh = 0.0_r8
         call mpiio_write(fh,fresh,this_block,distrb_clinic,subarray_type,date_count,all_offset)  
 !      write(*,*)'finish fresh'  


      if(mytid==0)then
      open(57,file=fname,form='unformatted',status='old',access='direct',recl=4,action='write')
      write(57,rec=(all_offset/4+1)) number_month
      write(57,rec=(all_offset/4+2)) number_day
      close(57)
      end if
      all_offset=all_offset+8


#ifdef COUP
      where (vit(:,:,1,:) < 0.5D0) t_cpl = 0.0_r8
      call mpiio_write(fh,t_cpl,this_block,distrb_clinic,subarray_type,date_count,all_offset) 
!         write(*,*)'finish t_cpl'
      where (vit(:,:,1,:) < 0.5D0) s_cpl = 0.0_r8
      call mpiio_write(fh,s_cpl,this_block,distrb_clinic,subarray_type,date_count,all_offset)
 !        write(*,*)'finish s_cpl'
      where (vit(:,:,1,:) < 0.5D0) u_cpl = 0.0_r8
      call mpiio_write(fh,u_cpl,this_block,distrb_clinic,subarray_type,date_count,all_offset)
 !        write(*,*)'finish u_cpl'
      where (vit(:,:,1,:) < 0.5D0) v_cpl = 0.0_r8
      call mpiio_write(fh,v_cpl,this_block,distrb_clinic,subarray_type,date_count,all_offset)
 !        write(*,*)'finish v_cpl'
      where (vit(:,:,1,:) < 0.5D0) dhdx = 0.0_r8
      call mpiio_write(fh,dhdx,this_block,distrb_clinic,subarray_type,date_count,all_offset)
 !        write(*,*)'finish dhdx'
      where (vit(:,:,1,:) < 0.5D0) dhdy = 0.0_r8
      call mpiio_write(fh,dhdy,this_block,distrb_clinic,subarray_type,date_count,all_offset)
 !        write(*,*)'finish dhdy'
      where (vit(:,:,1,:) < 0.5D0) qheat = 0.0_r8
      call mpiio_write(fh,qheat,this_block,distrb_clinic,subarray_type,date_count,all_offset)
 !        write(*,*)'finish qheat'
#endif

  call mpi_file_close(fh,ierror)

!write(*,*)'finish close 3'
  end if

      return
      end

#if (defined NETCDF) || (defined ALL)
      SUBROUTINE check_err (iret)
#include <netcdf.inc>
      INTEGER :: iret
      IF (iret /= NF_NOERR) THEN
         PRINT *, nf_strerror (iret)
         STOP
      END IF
      END SUBROUTINE check_err
#endif

