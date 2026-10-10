!  CVS: $Id: ssave-cdf.F90,v 1.1.1.1 2004/04/29 06:62:39 lhl Exp $
!     =================
subroutine SSAVE3HOURAVG(IIS)

#include <def-undef.h>
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

      integer, intent(in) :: IIS
      logical       :: write_restart     ! restart now
      character (len=22) :: fname
      integer :: klevel, iiday, nwmf1,mon01,number_day1,hour,last
      integer :: imd11
      integer(kind(1))   :: curr_ymd     ! Current date YYYYMMDD

      integer::fh,Byte_size_time,date_count
      integer*8 all_offset
      type (block) :: this_block  ! block info for current block
	  integer :: global_sizes(2),local_sizes(2),start_point(2),subarray_type
      integer * 8 :: c_nx, c_ny, c_km, c_month, c_day, k_n

         number_day1 = dd_licom !-1 !LPF 20170725

        
             mon01=mon0
             nwmf1= iyfm
             imd11=NMONTH (MON01)

         
#if (defined AVG3HOUR || defined OBP_AVG3HOUR)
     write_restart = .false.
     if ( daily_accum ) then !LPF20170725
     !   if (number_day == 1) write_restart = .true.
     !else if ( mod(iiday-1,rest_freq) == 0  ) then
        write_restart = .true.
     end if
     hour = IIS / (FSTEP/3)
     if(mytid==0) write(*,*)'in 3houravg,iday=,imd=,mon0=, nwmf=,hour=',number_day1,imd,mon01,nwmf1,hour
     if(mytid==0) write(*,*)'in 3houravg,write_restart', write_restart
!
     if ( write_restart  ) then

       c_nx = imt
       c_ny = jmt
       c_km = km
       call write_init(mpi_comm_ocn,c_nx, c_ny, c_km)
#if (defined AVG3HOUR)
       last = 0
       fname(1:8)='avg3hour'
       fname(13:13)='-'
       fname(16:16)='-'
       write(fname(14:15),'(i2.2)')mon01
       write(fname(9:12),'(i4.4)')nwmf1
       write(fname(17:18),'(i2.2)') number_day1
       write(fname(19:20),'(i2.2)') hour
       write(fname(21:22),'(i2.2)') last

       if ( mytid == 0 ) write(*,*)'avg3hour fname=',fname

       where (viv < 0.5D0) us3hour = 0.0_r4
       k_n = km
       call write_data_r4(us3hour(1,1,1,1), trim(fname)//CHAR(0), "us3hour"//CHAR(0), k_n)

 !        write(*,*)'finish U'
       where (viv < 0.5D0) vs3hour = 0.0_r4
       k_n = km
       call write_data_r4(vs3hour(1,1,1,1), trim(fname)//CHAR(0), "vs3hour"//CHAR(0), k_n)
 !        write(*,*)'finish V'
       where (vit < 0.5D0) ts3hour = 0.0_r4
       k_n = km
       call write_data_r4(ts3hour(1,1,1,1), trim(fname)//CHAR(0), "ts3hour"//CHAR(0), k_n)
 !        write(*,*)'finish at1'
       where (vit < 0.5D0) ss3hour = 0.0_r4
       k_n = km
       call write_data_r4(ss3hour(1,1,1,1), trim(fname)//CHAR(0), "ss3hour"//CHAR(0), k_n)
         k_n = km
         call write_data_r4(dens3hour(1,1,1,1), trim(fname)//CHAR(0), "dens3hour"//CHAR(0), k_n)
!         write(*,*)'finish at2'
!lhl20110728
        do klevel=1,km+1
            where (vit(:,:,klevel,:) < 0.5D0) ws3hour(:,:,klevel,:) = 0.0_r4
        end do
        k_n = km+1

        call write_data_kmp1_r4(ws3hour(1,1,1,1), trim(fname)//CHAR(0), "ws3hour"//CHAR(0), k_n)

       where (vit(:,:,1,:) < 0.5D0) z03hour = 0.0_r4
       k_n = 1
       call write_data_r4(z03hour(1,1,1), trim(fname)//CHAR(0), "z03hour"//CHAR(0), k_n)

       where (viv(:,:,1,:) < 0.5D0) su3hour = 0.0_r4
       k_n = 1
       call write_data_r4(su3hour(1,1,1), trim(fname)//CHAR(0), "su3hour"//CHAR(0), k_n)
!         write(*,*)'finish su'

       where (viv(:,:,1,:) < 0.5D0) sv3hour = 0.0_r4
       k_n = 1
       call write_data_r4(sv3hour(1,1,1), trim(fname)//CHAR(0), "sv3hour"//CHAR(0), k_n)
!         write(*,*)'finish sv'

       where (vit(:,:,1,:) < 0.5D0) swv3hour = 0.0_r4
       k_n = 1
       call write_data_r4(swv3hour(1,1,1), trim(fname)//CHAR(0), "swv3hour"//CHAR(0), k_n)
 !        write(*,*)'finish swv'

       where (vit(:,:,1,:) < 0.5D0) lwv3hour = 0.0_r4
       k_n = 1
       call write_data_r4(lwv3hour(1,1,1), trim(fname)//CHAR(0), "lwv3hour"//CHAR(0), k_n)
 !        write(*,*)'finish lwv'

       where (vit(:,:,1,:) < 0.5D0) sshf3hour = 0.0_r4
       k_n = 1
       call write_data_r4(sshf3hour(1,1,1), trim(fname)//CHAR(0), "sshf3hour"//CHAR(0), k_n)
 !        write(*,*)'finish sshf'

       where (vit(:,:,1,:) < 0.5D0) lthf3hour = 0.0_r4
       k_n = 1
       call write_data_r4(lthf3hour(1,1,1), trim(fname)//CHAR(0), "lthf3hour"//CHAR(0), k_n)

       where (vit(:,:,1,:) < 0.5D0) fresh3hour = 0.0_r4
       k_n = 1
       call write_data_r4(fresh3hour(1,1,1), trim(fname)//CHAR(0), "fresh3hour"//CHAR(0), k_n)


        if(mytid==0) write(*,*)'write u v at ws done'
#endif
#if (defined OBP_AVG3HOUR)
       last = 0
       fname(1:8)='obpavg3h'
       fname(13:13)='-'
       fname(16:16)='-'
       write(fname(14:15),'(i2.2)')mon01
       write(fname(9:12),'(i4.4)')nwmf1
       write(fname(17:18),'(i2.2)') number_day1
       write(fname(19:20),'(i2.2)') hour
       write(fname(21:22),'(i2.2)') last

       if ( mytid == 0 ) write(*,*)'obpavg3h fname=',fname
       where (vit(:,:,1,:) < 0.5D0) pbo3hour = 0.0_r4
       k_n = 1
       call write_data_r4(pbo3hour(1,1,1), trim(fname)//CHAR(0), "pbo3hour"//CHAR(0), k_n)

       where (vit(:,:,1,:) < 0.5D0) pso3hour = 0.0_r4
       k_n = 1
       call write_data_r4(pso3hour(1,1,1), trim(fname)//CHAR(0), "pso3hour"//CHAR(0), k_n)

       where (vit(:,:,1,:) < 0.5D0) ssh3hour = 0.0_r4
       k_n = 1
       call write_data_r4(ssh3hour(1,1,1), trim(fname)//CHAR(0), "ssh3hour"//CHAR(0), k_n)       
       if(mytid==0) write(*,*)'write obp done'
#endif

     end if
!==================================================================
#endif


!----------put the initial value 0.0------
      call AVG3HOURINI
        Num_3houracc=0
!----------put the initial value 0.0------

      return
END SUBROUTINE SSAVE3HOURAVG

