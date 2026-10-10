!  CVS: $Id: inirun.F90,v 1.1.1.1 2004/04/29 06:22:39 lhl Exp $
!     =================
      SUBROUTINE READTSFORRES (fnamet,fnames)
!     =================
!     INITIALIZING FOR ALL PHYSICAL FIELDS
 
#include <def-undef.h>
use param_mod
use pconst_mod
use dyn_mod
use tracer_mod
use cdf_mod,only:buffer_r4_global,buffer_r4_local 
use msg_mod, only: tag_1d,tag_2d,tag_3d,tag_4d,nproc,status,mpi_comm_ocn
use shr_msg_mod
use shr_mpi_mod
use shr_sys_mod
use domain
use constant_mod
use gather_scatter
use POP_HaloMod
use POP_GridHorzMod
use operators

      IMPLICIT NONE

  !----- local  ------
  integer            :: fid    ! nc domain file ID
  integer            :: dimid  ! nc dimension id
  integer            :: vid    ! nc variable ID
  integer            :: rcode  ! nc return code
  integer            :: ntim   ! temporary
  integer            :: iblock   ! temporary
  integer            :: ErrorCode   ! temporary
 
#include <netcdf.inc>
!
!     Define Variables.
      integer*4,  dimension(4) :: start(4)
      integer*4,  dimension(4) :: count(4)
      character (len=*) :: fnamet,fnames
      INTEGER :: NMFF
      real(r4) :: xx,yy
      integer*4   :: ncid1, iret1
      integer*4 :: ncid,iret

      !allocate(h0(imt,jmt,max_blocks_clinic),u(imt,jmt,km,max_blocks_clinic),v(imt,jmt,km,max_blocks_clinic), &
      !         at(imt,jmt,km,ntra,max_blocks_clinic))
      if (.not. allocated(restore_at)) then
       allocate(restore_at(imt,jmt,km,ntra,max_blocks_clinic))
         restore_at=0.0
      endif

      if (mytid==0)then
          write(6,*)"BEGINNING-----read_tsforres !"
      endif 

 
     if ( simple_assm ) then
     !if ( simple_assm==.true.) then
      !if (mytid==0) then
      if (.not. allocated(buffer_r4_global)) then
        allocate(buffer_r4_global(imt_global,jmt_global))
        !allocate(buffer(imt_global,jmt_global))
      endif
      if (.not. allocated(buffer_r4_local)) then
        allocate(buffer_r4_local(imt,jmt,max_blocks_clinic))
        !allocate(buffer(imt_global,jmt_global))
      endif
     
#if (defined HIGHRES)
!----------------------------------------------------
! Open netCDF file.
!----------------------------------------------------
      if (mytid==0) then
       iret=nf_open(fnames,nf_nowrite,ncid)
       call check_err (iret)

       iret1=nf_open(fnamet,nf_nowrite,ncid1)
       call check_err (iret1)
      end if
!

!----------------------------------------------------
!   Retrieve data
!----------------------------------------------------
     

      do k=1,km
!
       if (mytid==0) then
        start(1)=1 ; count(1)=imt_global
        start(2)=1 ; count(2)=jmt_global
        start(3)=k ; count(3)=1
        start(4)=1 ; count(4)=1

        iret=nf_get_vara_real(ncid,   7,start,count, buffer_r4_global)
        !iret=nf_get_vara_double(ncid,   7,start,count, buffer)
        call check_err (iret)
        do j=1 ,jmt_global/2
        do i=1, imt_global
           xx = buffer_r4_global(i,j)
           yy = buffer_r4_global(i,jmt_global+1-j)
           buffer_r4_global(i,jmt_global+1-j) = xx
           buffer_r4_global(i,j) = yy
           !xx = buffer(i,j)
           !yy = buffer(i,jmt_global+1-j)
           !buffer(i,jmt_global+1-j) = xx
           !buffer(i,j) = yy
        end do
        end do

       
      end if

       call scatter_global(buffer_r4_local, buffer_r4_global,master_task, distrb_clinic, &
                          field_loc_center, field_type_scalar)
         restore_at(:,:,k,2,:)=buffer_r4_local(:,:,:)*1.0D0
      


       if (mytid==0) then
        start(1)=1 ; count(1)=imt_global
        start(2)=1 ; count(2)=jmt_global
        start(3)=k ; count(3)=1
        start(4)=1 ; count(4)=1

        iret1=nf_get_vara_real(ncid1,   7,start,count, buffer_r4_global)
        !iret1=nf_get_vara_double(ncid1,   7,start,count, buffer)
        call check_err (iret)
        do j=1 ,jmt_global/2
        do i=1, imt_global
           xx = buffer_r4_global(i,j)
           yy = buffer_r4_global(i,jmt_global+1-j)
           buffer_r4_global(i,jmt_global+1-j) = xx
           buffer_r4_global(i,j) = yy
           !xx = buffer(i,j)
           !yy = buffer(i,jmt_global+1-j)
           !buffer(i,jmt_global+1-j) = xx
           !buffer(i,j) = yy
        end do
        end do

      end if


       call scatter_global(buffer_r4_local, buffer_r4_global,master_task, distrb_clinic, &
                          field_loc_center, field_type_scalar)
         restore_at(:,:,k,1,:)=buffer_r4_local(:,:,:)*1.0D0

     end do !kvt
      iret = nf_CLOSE (ncid)
      iret1 = nf_CLOSE (ncid1)
#endif

 
       call POP_HaloUpdate(restore_at(:,:,:,1,:) , POP_haloClinic, POP_gridHorzLocCenter,&
                       POP_fieldKindScalar, errorCode, fillValue = 0.0_r8)
       call POP_HaloUpdate(restore_at(:,:,:,2,:) , POP_haloClinic, POP_gridHorzLocCenter,&
                       POP_fieldKindScalar, errorCode, fillValue = 0.0_r8)

!----------------------------------------------------
!   assign 0 to land grids of TSinital
!----------------------------------------------------
!$OMP PARALLEL DO PRIVATE (IBLOCK,K,J,I)
      DO IBLOCK= 1, NBLOCKS_CLINIC
         DO K = 1,KVT
            DO J = 1,JMT
               DO I = 1,IMT
                   if(VIT(I,J,K,IBLOCK).gt.0.5) then
                   if( RESTORE_AT(I,J,K,1,IBLOCK).gt.50.or.RESTORE_AT (I,J,K,1,IBLOCK).lt.-5.0  ) &
                     write(16,*) 'resAT (I,J,K,1,IBLOCK)=',i,j,k,RESTORE_AT (I,J,K,1,IBLOCK) 
                   if( RESTORE_AT(I,J,K,2,IBLOCK).gt.70.or.RESTORE_AT (I,J,K,2,IBLOCK).lt.0  ) &
                     write(16,*) 'resAT (I,J,K,2,IBLOCK)=',i,j,k,RESTORE_AT (I,J,K,2,IBLOCK) 
                   endif
                  RESTORE_AT (I,J,K,1,IBLOCK) = RESTORE_AT (I,J,K,1,IBLOCK)*VIT(I,J,K,IBLOCK)
                  RESTORE_AT (I,J,K,2,IBLOCK) = (RESTORE_AT (I,J,K,2,IBLOCK)- 35.0D0)*0.001D0*VIT(I,J,K,IBLOCK)
               END DO
            END DO
         END DO
      END DO

     !if (mytid==0) then
      deallocate(buffer_r4_global)
      deallocate(buffer_r4_local)
      !deallocate(buffer)
     !endif
 
     end if
   
    END SUBROUTINE READTSFORRES 
