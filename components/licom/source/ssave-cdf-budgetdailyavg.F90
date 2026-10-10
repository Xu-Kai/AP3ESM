!  CVS: $Id: ssave-cdf.F90,v 1.1.1.1 2004/04/29 06:62:39 lhl Exp $
!     =================
      subroutine SSAVEBUDGETDAILYAVG
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
use distribution
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
	  
         number_day1 = iday !-1 !LPF 20170725
         number_month = month
         !if ( number_day1  > imd ) then
         !    number_day1 = 1
         !    number_month = month + 1
         !end if
        
             mon01=mon0
             nwmf1= iyfm
             number_day1 = iday-1 !nnmonth(mon01)
!             number_day1 = nnmonth(mon01)
             imd11=NMONTH (MON01) 
        if ( iday ==1 ) then
         if(mon0>1) then
             nwmf1= iyfm
             mon01=mon0-1
          elseif(mon0==1) then
             nwmf1= iyfm-1
             mon01=12
         endif
             number_day1 = nmonth(mon01)
         endif

        if( iday  > imd11 ) then
             number_day1 = 1 
             number_month = month + 1
        !else
             !mon01=mon0
             !nwmf1= iyfm
             !number_day1 = nnmonth(mon01)
         !    number_day1 = iday-1 !nnmonth(mon01)
        endif
         
#if (defined DAILYBUGDET)
     write_restart = .false.
     if ( dailybudget_accum ) then !LPF20170725
     !   if (number_day == 1) write_restart = .true.
     !else if ( mod(iiday-1,rest_freq) == 0  ) then
        write_restart = .true.
     end if
!
!    if ( mod(iday,rest_freq) == 0 .or. iday == imd .or. iday == 10 .or. iday ==20) then
!    if ( mod(iiday,rest_freq) == 1  ) then
     if ( write_restart  ) then
!

         fname(1:8)='budgdyvg'
         fname(13:13)='-'
         fname(16:16)='-'
         write(fname(14:15),'(i2.2)')mon01
         write(fname(9:12),'(i4.4)')nwmf1
         write(fname(17:18),'(i2.2)') number_day1
 
      !real(r4),dimension(imt,jmt,km,NTRA,max_blocks_clinic)::tenddaily
      !real(r4),dimension(imt,jmt,km,NTRA,max_blocks_clinic)::dt_diffdaily
      !real(r4),dimension(imt,jmt,km,NTRA,max_blocks_clinic)::axdaily,aydaily,azdaily
      !real(r4),dimension(imt,jmt,km,NTRA,max_blocks_clinic)::dxdaily,dydaily,dzdaily
      !real(r4),dimension(imt,jmt,km,max_blocks_clinic)::pendaily
      !real(r4),dimension(imt,jmt,km,NTRA,max_blocks_clinic)::dt_convdaily
      if ( mytid==0 ) then
       write(*,*)'fname=',fname
      end if !mytid==0
	  
!=================================================================liyl
	  
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

!klv=45 for upper1500m !LPF20170729 !to save the memorysize

     DO N=1,NTRA


		 
      where (vit < 0.5D0) axdaily(:,:,:,N,:) = 0.0_r4
      do klevel=1,klv
		call mpiio_write(fh,axdaily(:,:,klevel,N,:),this_block,distrb_clinic,subarray_type,date_count,all_offset)		  
       end do
!         write(*,*)'adv u'

      where (vit < 0.5D0) aydaily(:,:,:,N,:) = 0.0_r4
      do klevel=1,klv
		call mpiio_write(fh,aydaily(:,:,klevel,N,:),this_block,distrb_clinic,subarray_type,date_count,all_offset)		  
       end do
!         write(*,*)'adv y'

      where (vit < 0.5D0) azdaily(:,:,:,N,:) = 0.0_r4
      do klevel=1,klv
		call mpiio_write(fh,azdaily(:,:,klevel,N,:),this_block,distrb_clinic,subarray_type,date_count,all_offset)		  
       end do
!         write(*,*)'adv z'

      where (vit < 0.5D0) dxdaily(:,:,:,N,:) = 0.0_r4
      do klevel=1,klv
		call mpiio_write(fh,dxdaily(:,:,klevel,N,:),this_block,distrb_clinic,subarray_type,date_count,all_offset)		  
       end do
!         write(*,*)'diff dx'


      where (vit < 0.5D0) dzdaily(:,:,:,N,:) = 0.0_r4
      do klevel=1,klv
		call mpiio_write(fh,dzdaily(:,:,klevel,N,:),this_block,distrb_clinic,subarray_type,date_count,all_offset)		  
       end do
!         write(*,*)'diff dz'

      where (vit < 0.5D0) dt_diffdaily(:,:,:,N,:) = 0.0_r4
      do klevel=1,klv
		call mpiio_write(fh,dt_diffdaily(:,:,klevel,N,:),this_block,distrb_clinic,subarray_type,date_count,all_offset)		  
       end do
!         write(*,*)'diff dt impilict'


      where (vit < 0.5D0) dt_convdaily(:,:,:,N,:) = 0.0_r4
      do klevel=1,klv
		call mpiio_write(fh,dt_convdaily(:,:,klevel,N,:),this_block,distrb_clinic,subarray_type,date_count,all_offset)		  
       end do
!         write(*,*)'diff dt due to conv'



      where (vit < 0.5D0) tenddaily(:,:,:,N,:) = 0.0_r4
      do klevel=1,klv
		call mpiio_write(fh,tenddaily(:,:,klevel,N,:),this_block,distrb_clinic,subarray_type,date_count,all_offset)		  
       end do
!         write(*,*)'tendency'

     if ( simple_assm ) then 
     !if ( simple_assm==.true.) then 
      where (vit < 0.5D0) dt_restoredaily(:,:,:,N,:) = 0.0_r4
      do klevel=1,klv
		call mpiio_write(fh,dt_restoredaily(:,:,klevel,N,:),this_block,distrb_clinic,subarray_type,date_count,all_offset)		  
       end do
!         write(*,*)'restore'
     endif !assm

	   where (vit(:,:,1,:) < 0.5D0) net_daily(:,:,N,:) = 0.0_r4
		call mpiio_write(fh,net_daily(:,:,N,:),this_block,distrb_clinic,subarray_type,date_count,all_offset)		  
!         write(*,*)'net'

    ENDDO !for tracer number
    
	 where (vit < 0.5D0) pendaily = 0.0_r4
     do klevel=1,klv
		call mpiio_write(fh,pendaily(:,:,klevel,:),this_block,distrb_clinic,subarray_type,date_count,all_offset)		  
     end do
      
!      where (vit(:,:,1,:) < 0.5D0) mlddaily = 0.0_r4
!     call gather_global(buffer_r4_global, mlddaily, master_task,distrb_clinic)
!     if (mytid==0) then
!         WRITE (72)buffer_r4_global
!     end if
!         write(*,*)'finish mld'

      if(mytid==0)then
      open(72,file=fname,form='unformatted',status='old',access='direct',recl=4,action='write')
      write(72,rec=(all_offset/4+1)) number_month
      write(72,rec=(all_offset/4+2)) number_day1
      close(72)
      end if

	call mpi_file_close(fh,ierror)
!
  end if !write or not
!      deallocate(axdaily,aydaily,azdaily,dxdaily,dzdaily,dt_diffdaily)
!      deallocate(tenddaily,dt_convdaily,pendaily)
#endif

!----------put the initial value 0.0------
!      call DAILYAVGINI
!----------put the initial value 0.0------

      return
       END SUBROUTINE SSAVEBUDGETDAILYAVG

