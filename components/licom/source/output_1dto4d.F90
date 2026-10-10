
module output_netcdf

#include <def-undef.h>
use precision_mod
use param_mod
use pconst_mod
use output_mod
use cdf_mod !, only: time_dim,y_dim,x_dim,buffer_r4_local,buffer_r4_global
use domain
use distribution
use gather_scatter
      IMPLICIT NONE
#include <netcdf.inc>
  public :: output1d
  public :: output_define
  public :: output_attribute
  public :: output2d
!  public :: output_writedata
  SAVE
  private

!==========================================================================
  contains
!==========================================================================


!  CVS: $Id: output1d.F90,v 1.1.1.1 2004/04/29 06:22:39 lhl Exp $
!     ============================================
      subroutine output1d()
!     ============================================
! output 1D data by LPF 
!      IMPLICIT NONE
!      logical :: hist_output,rest_output 
          return
          end subroutine output1d 

!     ================================================================================
     subroutine output_define(x2d_id,x2d_rank,x2d_dims,var2d_shortname)
     !subroutine output2d(var2d,ixx,jyy,kk,lat_len,lon_len,x2d_id,ncid1,iret1,start3,count3)
!     ================================================================================
      !IMPLICIT NONE
      integer :: x2d_id !,ncid1,iret1
      integer:: x2d_rank 
      integer:: x2d_dims(x2d_rank)
      character(len=*):: var2d_shortname !,var2d_unit,var2d_coord

       if(mytid==0) then
        if(x2d_rank > 3) then
         x2d_dims (4) = time_dim
         x2d_dims (3) = lev_dim
         x2d_dims (2) = y_dim
         x2d_dims (1) = x_dim
        else
         x2d_dims (3) = time_dim
         x2d_dims (2) = y_dim
         x2d_dims (1) = x_dim
        endif
         iret = nf_def_var (ncid, var2d_shortname, NF_REAL, x2d_rank, x2d_dims, x2d_id)
         CALL check_err (iret)
       endif
             return
        end subroutine output_define 

!     ================================================================================
     subroutine output_attribute(x2d_id,x2d_rank,x2d_dims,var2d_longname,var2d_unit,var2d_coord)
!     ================================================================================
      integer :: x2d_id !,ncid1,iret1
      integer:: x2d_rank,len_longname,len_unit,len_coord 
      integer:: x2d_dims(x2d_rank)
      character(len=*):: var2d_longname,var2d_unit,var2d_coord


       if(mytid==0) then
        len_unit=len_trim(var2d_unit)
        len_longname=len_trim(var2d_longname)
        len_coord=len_trim(var2d_coord)

       !if(mytid==0) then
         iret = nf_put_att_text (ncid, x2d_id, 'long_name', len_longname, var2d_longname)
         !iret = nf_put_att_text (ncid, x2d_id, 'long_name', 18, 'sea surface height')
         CALL check_err (iret)
         iret = nf_put_att_text (ncid, x2d_id, 'units', len_unit, var2d_unit)
         !iret = nf_put_att_text (ncid, x2d_id, 'units', 5, 'meter')
         CALL check_err (iret)
         iret = nf_put_att_real (ncid, x2d_id, '_FillValue', NF_REAL, 1, spval)
         CALL check_err (iret)
         iret = nf_put_att_text (ncid, x2d_id, 'coordinates', len_coord, var2d_coord)
         !iret = nf_put_att_text (ncid, x2d_id, 'coordinates', 7,  'lon lat')
         CALL check_err (iret)
       
       endif
             return
        end subroutine output_attribute 
      

!     ================================================================================
     subroutine output2d(x2d_id,x2d_rank,x2d_dims,starttmp,counttmp)
     !subroutine output2d(var2d,ixx,jyy,kk,lat_len,lon_len,x2d_id,ncid1,iret1,start3,count3)
!     ================================================================================

      integer :: x2d_id !,ncid1,iret1
      !integer :: kk,x2d_id,ncid1,iret1,ixx,jyy
      integer :: x2d_rank 
      !real(r4)    :: var2d(ixx,jyy,kk) !,x2d(ixx,jyy)
      integer:: starttmp(x2d_rank),counttmp(x2d_rank)
      integer :: x2d_dims(x2d_rank)
      !integer:: lon_len,lat_len
          

       if(mytid==0) then
!        if (x2d_rank >3 ) then   
!         starttmp (1)= 1
!         starttmp (2)= 1
!         starttmp (4)= 1
!         counttmp (1)= lon_len
!         counttmp (2)= lat_len
!         counttmp (4)= 1
!         starttmp (3)= k
!         counttmp (3)= 1
!        else
!         starttmp (1)= 1
!         starttmp (2)= 1
!         starttmp (3)= 1
!         counttmp (1)= lon_len
!         counttmp (2)= lat_len
!         counttmp (3)= 1
!        endif 



!         allocate(buffer_r4_global(imt_global,jmt_global), buffer_r4_local(imt,jmt,max_blocks_clinic))
!
         !where ( vit(:,:,1,:) > 0.5D0 )
         !   buffer_r4_local = var2d 
         !elsewhere
         !   buffer_r4_local = spval
         !endwhere 
         !call gather_global(buffer_r4_global,buffer_r4_local, master_task,distrb_clinic) 
         
            iret = nf_put_vara_real (ncid,x2d_id,starttmp, counttmp, buffer_r4_global)
            CALL check_err (iret)
        endif !mytid==0


             return
        end subroutine output2d 

end module output_netcdf
