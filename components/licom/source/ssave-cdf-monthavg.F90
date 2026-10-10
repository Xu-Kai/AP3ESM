!  CVS: $Id: ssave-cdf.F90,v 1.1.1.1 2004/04/29 06:62:39 lhl Exp $
!     =================
      subroutine SSAVEMONTHAVG
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

         number_day1 = 0 !-1 !LPF 20170725

!         if(yeartoggle .eqv. .true.) then
!             mon01 = 12
!             nwmf1= iyfm - 1
!         else
#ifdef COUP
             mon01=mon0
#else
             mon01=mon0 - 1
#endif
             nwmf1= iyfm
!         end if

!             number_day1 = iday-1 !nnmonth(mon01)
!             number_day1 = nnmonth(mon01)
             imd11=NMONTH (MON01) 

         
#if (defined DAILYACC || defined MONTHACC)
     write_restart = .false.
     if ( daily_accum ) then !LPF20170725
     !   if (number_day == 1) write_restart = .true.
     !else if ( mod(iiday-1,rest_freq) == 0  ) then
        write_restart = .true.
     end if
     if(mytid==0) write(*,*)'in monthavg,iday=,imd=,mon0=, nwmf= ',number_day1,imd,mon01,nwmf1
     if(mytid==0) write(*,*)'in monthavg,write_restart', write_restart

     if ( write_restart  ) then
!

         fname(1:8)='monthavg'
         fname(13:13)='-'
         fname(16:16)='-'
         write(fname(14:15),'(i2.2)')mon01
         write(fname(9:12),'(i4.4)')nwmf1
         write(fname(17:18),'(i2.2)') number_day1
 

       if ( mytid == 0 ) then
         write(*,*)'fname=',fname
       end if

       c_nx = imt
       c_ny = jmt
       c_km = km
       call write_init(mpi_comm_ocn,c_nx, c_ny, c_km)

       where (viv < 0.5D0) usmon = 0.0_r4
       k_n = km
       call write_data_r4(usmon(1,1,1,1), trim(fname)//CHAR(0), "usmon"//CHAR(0), k_n)

 !        write(*,*)'finish U'
       where (viv < 0.5D0) vsmon = 0.0_r4
       k_n = km
       call write_data_r4(vsmon(1,1,1,1), trim(fname)//CHAR(0), "vsmon"//CHAR(0), k_n)
 !        write(*,*)'finish V'
       where (vit < 0.5D0) tsmon = 0.0_r4
       k_n = km
       call write_data_r4(tsmon(1,1,1,1), trim(fname)//CHAR(0), "tsmon"//CHAR(0), k_n)
 !        write(*,*)'finish at1'
       where (vit < 0.5D0) ssmon = 0.0_r4
       k_n = km
       call write_data_r4(ssmon(1,1,1,1), trim(fname)//CHAR(0), "ssmon"//CHAR(0), k_n)

       where (vit < 0.5D0) densmon = 0.0_r4
       k_n = km
       call write_data_r4(densmon(1,1,1,1), trim(fname)//CHAR(0), "densmon"//CHAR(0), k_n)
!         write(*,*)'finish at2'
!lhl20110728
        do klevel=1,km+1
            where (vit(:,:,klevel,:) < 0.5D0) wsmon(:,:,klevel,:) = 0.0_r4
        end do
        k_n = km+1
        call write_data_kmp1_r4(wsmon(1,1,1,1), trim(fname)//CHAR(0), "wsmon"//CHAR(0), k_n)

       where (vit(:,:,1,:) < 0.5D0) z0mon = 0.0_r4
       k_n = 1
       call write_data_r4(z0mon(1,1,1), trim(fname)//CHAR(0), "z0mon"//CHAR(0), k_n)

       where (viv(:,:,1,:) < 0.5D0) sumon = 0.0_r4
       k_n = 1
       call write_data_r4(sumon(1,1,1), trim(fname)//CHAR(0), "sumon"//CHAR(0), k_n)
!         write(*,*)'finish su'

       where (viv(:,:,1,:) < 0.5D0) svmon = 0.0_r4
       k_n = 1
       call write_data_r4(svmon(1,1,1), trim(fname)//CHAR(0), "svmon"//CHAR(0), k_n)
!         write(*,*)'finish sv'

       where (vit(:,:,1,:) < 0.5D0) swvmon = 0.0_r4
       k_n = 1
       call write_data_r4(swvmon(1,1,1), trim(fname)//CHAR(0), "swvmon"//CHAR(0), k_n)
 !        write(*,*)'finish swv'

       where (vit(:,:,1,:) < 0.5D0) lwvmon = 0.0_r4
       k_n = 1
       call write_data_r4(lwvmon(1,1,1), trim(fname)//CHAR(0), "lwvmon"//CHAR(0), k_n)
 !        write(*,*)'finish lwv'

       where (vit(:,:,1,:) < 0.5D0) sshfmon = 0.0_r4
       k_n = 1
       call write_data_r4(sshfmon(1,1,1), trim(fname)//CHAR(0), "sshfmon"//CHAR(0), k_n)
 !        write(*,*)'finish sshf'

       where (vit(:,:,1,:) < 0.5D0) lthfmon = 0.0_r4
       k_n = 1
       call write_data_r4(lthfmon(1,1,1), trim(fname)//CHAR(0), "lthfmon"//CHAR(0), k_n)

       where (vit(:,:,1,:) < 0.5D0) freshmon = 0.0_r4
       k_n = 1
       call write_data_r4(freshmon(1,1,1), trim(fname)//CHAR(0), "freshmon"//CHAR(0), k_n)


!      where (vit(:,:,1,:) < 0.5D0) precdaily = 0.0_r4
!         call mpiio_write(fh,precdaily,this_block,distrb_clinic,subarray_type,date_count,all_offset)
!!         write(*,*)'finish prec'

!      where (vit(:,:,1,:) < 0.5D0) evapdaily = 0.0_r4
!         call mpiio_write(fh,evapdaily,this_block,distrb_clinic,subarray_type,date_count,all_offset)
!!         write(*,*)'finish evap'

!      where (vit(:,:,1,:) < 0.5D0) roffdaily = 0.0_r4
!         call mpiio_write(fh,roffdaily,this_block,distrb_clinic,subarray_type,date_count,all_offset)
!!         write(*,*)'finish roffincluding iceoff'



  if(mytid==0) write(*,*)'write month average done'



end if	  
!==================================================================
#endif


!----------put the initial value 0.0------
      call MM00
!----------put the initial value 0.0------

      return
       END SUBROUTINE SSAVEMONTHAVG

