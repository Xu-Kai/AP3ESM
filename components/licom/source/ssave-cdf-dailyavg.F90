!  CVS: $Id: ssave-cdf.F90,v 1.1.1.1 2004/04/29 06:62:39 lhl Exp $
!     =================
      subroutine SSAVEDAILYAVG
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
!use mct_mod
!use esmf_mod
!use seq_flds_mod
!use seq_cdata_mod
!use seq_infodata_mod
!use seq_timemgr_mod
use output_mod
use msg_mod
use MPIIO
use module_add
use blocks
  implicit none


      logical       :: write_restart     ! restart now
      character (len=18) :: fname
      integer :: klevel, iiday, nwmf1,mon01,number_day1
      integer :: imd11
      integer(kind(1))   :: curr_ymd     ! Current date YYYYMMDD

      integer::fh,Byte_size_time,date_count
      integer*8 all_offset
      type (block) :: this_block  ! block info for current block
	  integer :: global_sizes(2),local_sizes(2),start_point(2),subarray_type
      integer * 8 :: c_nx, c_ny, c_km, c_month, c_day, k_n

         number_day1 = dd_licom !-1 !LPF 20170725
!         number_month = month
         !if ( number_day1  > imd ) then
         !    number_day1 = 1
         !    number_month = month + 1
         !end if
        
             mon01=mon0
             nwmf1= iyfm
!             number_day1 = iday-1 !nnmonth(mon01)
!             number_day1 = nnmonth(mon01)
             imd11=NMONTH (MON01) 
!        if ( iday ==1 ) then
!         if(mon0>1) then
!             nwmf1= iyfm
!             mon01=mon0-1
!          elseif(mon0==1) then
!             nwmf1= iyfm-1
!             mon01=12
!         endif
!             number_day1 = nmonth(mon01)
!         endif

!        if( iday  > imd11 ) then
!             number_day1 = 1
!             number_month = month + 1
!        !else
!             !mon01=mon0
!             !nwmf1= iyfm
!             !number_day1 = nnmonth(mon01)
!         !    number_day1 = iday-1 !nnmonth(mon01)
!        endif
         
#if (defined DAILYACC)
     write_restart = .false.
     if ( daily_accum ) then !LPF20170725
     !   if (number_day == 1) write_restart = .true.
     !else if ( mod(iiday-1,rest_freq) == 0  ) then
        write_restart = .true.
     end if
     if(mytid==0) write(*,*)'in dailyavg,iday=,imd=,mon0=, nwmf= ',number_day1,imd,mon01,nwmf1
     if(mytid==0) write(*,*)'in dailyavg,write_restart', write_restart
!
!    if ( mod(iday,rest_freq) == 0 .or. iday == imd .or. iday == 10 .or. iday ==20) then
!    if ( mod(iday,rest_freq) == 0 .or. iday == 1 ) then
!    if ( mod(iiday,rest_freq) == 1  ) then
     if ( write_restart  ) then
!

         fname(1:8)='dailyavg'
         fname(13:13)='-'
         fname(16:16)='-'
         write(fname(14:15),'(i2.2)')mon01
         write(fname(9:12),'(i4.4)')nwmf1
         write(fname(17:18),'(i2.2)') number_day1
 
!        if ( number_day ==1 .and. mon0 < 12) then
!            write(fname(14:15),'(i2.2)')mon0+1
!        end if
!
!        if ( number_day ==1 .and. mon0 == 12) then
!            write(fname(14:15),'(i2.2)')mon0-11
!            write(fname(9:12),'(i4.4)')nwmf+1
!        end if
       if ( mytid == 0 ) then
         write(*,*)'fname=',fname
       end if


#ifdef SPLIT_WRITE
       c_nx = imt
       c_ny = jmt
       c_km = km
       call write_init(mpi_comm_ocn,c_nx, c_ny, c_km)

       where (viv < 0.5D0) usdaily = 0.0_r4
       k_n = km
       call write_data_r4(usdaily(1,1,1,1), trim(fname)//CHAR(0), "usdaily"//CHAR(0), k_n)

 !        write(*,*)'finish U'
       where (viv < 0.5D0) vsdaily = 0.0_r4
       k_n = km
       call write_data_r4(vsdaily(1,1,1,1), trim(fname)//CHAR(0), "vsdaily"//CHAR(0), k_n)
 !        write(*,*)'finish V'
       where (vit < 0.5D0) tsdaily = 0.0_r4
       k_n = km
       call write_data_r4(tsdaily(1,1,1,1), trim(fname)//CHAR(0), "tsdaily"//CHAR(0), k_n)
 !        write(*,*)'finish at1'
       where (vit < 0.5D0) ssdaily = 0.0_r4
       k_n = km
       call write_data_r4(ssdaily(1,1,1,1), trim(fname)//CHAR(0), "ssdaily"//CHAR(0), k_n)
!         write(*,*)'finish at2'
       k_n = km
       where (vit < 0.5D0) densdaily = 0.0_r4
       call write_data_r4(densdaily(1,1,1,1), trim(fname)//CHAR(0), "densdaily"//CHAR(0), k_n)
!lhl20110728
        do klevel=1,km+1
            where (vit(:,:,klevel,:) < 0.5D0) wsdaily(:,:,klevel,:) = 0.0_r4
        end do
        k_n = km+1
!        do j=1,jmt
!        do i=1,imt
!        if(wsdaily(i,j,1,1) .eq. 1.0D0) write(*,*) "wsdaily is 1, i,j",i,j,mytid
!        end do
!        end do
        call write_data_kmp1_r4(wsdaily(1,1,1,1), trim(fname)//CHAR(0), "wsdaily"//CHAR(0), k_n)

       where (vit(:,:,1,:) < 0.5D0) z0daily = 0.0_r4
       k_n = 1
       call write_data_r4(z0daily(1,1,1), trim(fname)//CHAR(0), "z0daily"//CHAR(0), k_n)

       where (viv(:,:,1,:) < 0.5D0) sudaily = 0.0_r4
       k_n = 1
       call write_data_r4(sudaily(1,1,1), trim(fname)//CHAR(0), "sudaily"//CHAR(0), k_n)
!         write(*,*)'finish su'

       where (viv(:,:,1,:) < 0.5D0) svdaily = 0.0_r4
       k_n = 1
       call write_data_r4(svdaily(1,1,1), trim(fname)//CHAR(0), "svdaily"//CHAR(0), k_n)
!         write(*,*)'finish sv'

       where (vit(:,:,1,:) < 0.5D0) swvdaily = 0.0_r4
       k_n = 1
       call write_data_r4(swvdaily(1,1,1), trim(fname)//CHAR(0), "swvdaily"//CHAR(0), k_n)
 !        write(*,*)'finish swv'

       where (vit(:,:,1,:) < 0.5D0) lwvdaily = 0.0_r4
       k_n = 1
       call write_data_r4(lwvdaily(1,1,1), trim(fname)//CHAR(0), "lwvdaily"//CHAR(0), k_n)
 !        write(*,*)'finish lwv'

       where (vit(:,:,1,:) < 0.5D0) sshfdaily = 0.0_r4
       k_n = 1
       call write_data_r4(sshfdaily(1,1,1), trim(fname)//CHAR(0), "sshfdaily"//CHAR(0), k_n)
 !        write(*,*)'finish sshf'

       where (vit(:,:,1,:) < 0.5D0) lthfdaily = 0.0_r4
       k_n = 1
       call write_data_r4(lthfdaily(1,1,1), trim(fname)//CHAR(0), "lthfdaily"//CHAR(0), k_n)

       where (vit(:,:,1,:) < 0.5D0) freshdaily = 0.0_r4
       k_n = 1
       call write_data_r4(freshdaily(1,1,1), trim(fname)//CHAR(0), "freshdaily"//CHAR(0), k_n)

 !        write(*,*)'finish lthf'
!      where (vit(:,:,1,:) < 0.5D0) mlddaily = 0.0_r4
!      k_n = 1
!      call write_data_r4_(mlddaily(1,1,1), trim(fname)//CHAR(0), "mlddaily"//CHAR(0), k_n)
!
 !      write(*,*)'finish fresh'


!      where (vit(:,:,1,:) < 0.5D0) precdaily = 0.0_r4
!         call mpiio_write(fh,precdaily,this_block,distrb_clinic,subarray_type,date_count,all_offset)
!!         write(*,*)'finish prec'

!      where (vit(:,:,1,:) < 0.5D0) evapdaily = 0.0_r4
!         call mpiio_write(fh,evapdaily,this_block,distrb_clinic,subarray_type,date_count,all_offset)
!!         write(*,*)'finish evap'

!      where (vit(:,:,1,:) < 0.5D0) roffdaily = 0.0_r4
!         call mpiio_write(fh,roffdaily,this_block,distrb_clinic,subarray_type,date_count,all_offset)
!!         write(*,*)'finish roffincluding iceoff'



  if(mytid==0) write(*,*)'write u v at ws done'


!  c_month = number_month
!  c_day = number_day
!
!  call write_date_data(c_month, c_day, trim(fname)//CHAR(0), "date"//CHAR(0))


#else
!
!==================================================================liyl
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
		MPI_ORDER_FORTRAN, MPI_REAL, subarray_type,ierror)
		call MPI_Type_commit(subarray_type,ierror)
		exit
	end if
	end do
	
       where (vit(:,:,1,:) < 0.5D0) z0daily = 0.0_r4
        call mpiio_write(fh,z0daily,this_block,distrb_clinic,subarray_type,date_count,all_offset)
 
      where (viv < 0.5D0) usdaily = 0.0_r4
      do klevel=1,km
		call mpiio_write(fh,usdaily(:,:,klevel,:),this_block,distrb_clinic,subarray_type,date_count,all_offset)		  
       end do

 !        write(*,*)'finish U'
      where (viv < 0.5D0) vsdaily = 0.0_r4
      do klevel=1,km
		  call mpiio_write(fh,vsdaily(:,:,klevel,:),this_block,distrb_clinic,subarray_type,date_count,all_offset)
       end do
 !        write(*,*)'finish V'
      where (vit < 0.5D0) tsdaily(:,:,:,:) = 0.0_r4
      do klevel=1,km
		 call mpiio_write(fh,tsdaily(:,:,klevel,:),this_block,distrb_clinic,subarray_type,date_count,all_offset)
      end do
 !        write(*,*)'finish at1'
      where (vit < 0.5D0) ssdaily(:,:,:,:) = 0.0_r4
      do klevel=1,km
		call mpiio_write(fh,ssdaily(:,:,klevel,:),this_block,distrb_clinic,subarray_type,date_count,all_offset)		 
      end do
!         write(*,*)'finish at2'
!lhl20110728
      do klevel=1,km
        where (vit(:,:,klevel,:) < 0.5D0) wsdaily(:,:,klevel,:) = 0.0_r4
		call mpiio_write(fh,wsdaily(:,:,klevel,:),this_block,distrb_clinic,subarray_type,date_count,all_offset)		 
      end do

	  
      where (viv(:,:,1,:) < 0.5D0) sudaily = 0.0_r4
         call mpiio_write(fh,sudaily,this_block,distrb_clinic,subarray_type,date_count,all_offset)
!         write(*,*)'finish su'
      where (viv(:,:,1,:) < 0.5D0) svdaily = 0.0_r4
         call mpiio_write(fh,svdaily,this_block,distrb_clinic,subarray_type,date_count,all_offset)
!         write(*,*)'finish sv'
      where (vit(:,:,1,:) < 0.5D0) swvdaily = 0.0_r4
         call mpiio_write(fh,swvdaily,this_block,distrb_clinic,subarray_type,date_count,all_offset)
 !        write(*,*)'finish swv'
      where (vit(:,:,1,:) < 0.5D0) lwvdaily = 0.0_r4
         call mpiio_write(fh,lwvdaily,this_block,distrb_clinic,subarray_type,date_count,all_offset)
 !        write(*,*)'finish lwv'
      where (vit(:,:,1,:) < 0.5D0) sshfdaily = 0.0_r4
         call mpiio_write(fh,sshfdaily,this_block,distrb_clinic,subarray_type,date_count,all_offset)
 !        write(*,*)'finish sshf'
      where (vit(:,:,1,:) < 0.5D0) lthfdaily = 0.0_r4
         call mpiio_write(fh,lthfdaily,this_block,distrb_clinic,subarray_type,date_count,all_offset)
 !        write(*,*)'finish lthf'
      where (vit(:,:,1,:) < 0.5D0) mlddaily = 0.0_r4
         call mpiio_write(fh,mlddaily,this_block,distrb_clinic,subarray_type,date_count,all_offset)  
 !      write(*,*)'finish fresh'  


      where (vit(:,:,1,:) < 0.5D0) precdaily = 0.0_r4
         call mpiio_write(fh,precdaily,this_block,distrb_clinic,subarray_type,date_count,all_offset)  
!         write(*,*)'finish prec'

      where (vit(:,:,1,:) < 0.5D0) evapdaily = 0.0_r4
         call mpiio_write(fh,evapdaily,this_block,distrb_clinic,subarray_type,date_count,all_offset)  
!         write(*,*)'finish evap'		 

      where (vit(:,:,1,:) < 0.5D0) roffdaily = 0.0_r4
         call mpiio_write(fh,roffdaily,this_block,distrb_clinic,subarray_type,date_count,all_offset) 
!         write(*,*)'finish roffincluding iceoff'		 


      if(mytid==0)then
      open(62,file=fname,form='unformatted',status='old',access='direct',recl=4,action='write')
      write(62,rec=(all_offset/4+1)) number_month
      write(62,rec=(all_offset/4+2)) number_day1
      close(62)
      end if

	  call mpi_file_close(fh,ierror)

#endif
end if	  
!==================================================================
#endif

!
  !end if

!----------put the initial value 0.0------
      call DAILYAVGINI
      Num_outputacc=0
!----------put the initial value 0.0------

      return
       END SUBROUTINE SSAVEDAILYAVG

