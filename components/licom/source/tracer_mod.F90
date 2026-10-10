!  CVS: $Id: tracer_mod.F90,v 1.1.1.1 2004/04/29 06:22:39 lhl Exp $
module tracer_mod
#include <def-undef.h>
!
use precision_mod
use param_mod
use constant_mod
!     ------------------------------------------------------------------
!     U V T S H0 W RHO
!     ------------------------------------------------------------------

implicit none

       real(r8),dimension(imt,jmt,0:km,NTRA,max_blocks_clinic),target ::atb
      !real(r8),dimension(:,:,:,:,:), allocatable, target ::atb

      real(r8),dimension(imt,jmt,NTRA,max_blocks_clinic),target ::net
!mohr
      real(r8),dimension(imt,jmt,km,NTRA,max_blocks_clinic),target::at
!      real(r8),allocatable,dimension(:,:,:,:,:),target ::at
      real(r8),allocatable,dimension(:,:,:,:,:)::restore_at !LPF20170729
!      real(r8),allocatable,dimension(:,:,:,:,:)::atzwp !ZWP20170218
!
      ! real(r8),dimension(imt,jmt,km,max_blocks_clinic),target ::pdensity
      ! real(r8),dimension(imt,jmt,max_blocks_clinic)::amld
      real(r8),dimension(:,:,:,:), allocatable,target::pdensity
      real(r8),dimension(:,:,:), allocatable,target::amld
!lhl1204
! !
!       real(r8),dimension(imt,jmt,km,NTRA,max_blocks_clinic)::tend
!       real(r8),dimension(imt,jmt,km,NTRA,max_blocks_clinic)::ax,ay,az
!       real(r8),dimension(imt,jmt,km,NTRA,max_blocks_clinic)::dx,dy,dz
       real(r8),dimension(imt,jmt,km,max_blocks_clinic)::penetrate
!       real(r8),dimension(imt,jmt,km,NTRA,max_blocks_clinic)::dt_diff
! !
!       real(r8),dimension(imt,jmt,km,NTRA,max_blocks_clinic)::ddy
!       real(r8),dimension(imt,jmt,km,NTRA,max_blocks_clinic)::dt_conv
!       real(r8),dimension(imt,jmt,km,NTRA,max_blocks_clinic)::dt_restore !20170802
      real(r8),dimension(:,:,:,:,:), allocatable,target::tend
      real(r8),dimension(:,:,:,:,:), allocatable,target::ax,ay,az
      real(r8),dimension(:,:,:,:,:), allocatable,target::dx,dy,dz
!      real(r8),dimension(:,:,:,:), allocatable::penetrate
      real(r8),dimension(:,:,:,:,:), allocatable,target::dt_diff
      real(r8),dimension(:,:,:,:,:), allocatable::ddy
      real(r8),dimension(:,:,:,:,:), allocatable,target::dt_conv
      real(r8),dimension(:,:,:,:,:), allocatable::dt_restore
      real(r8),dimension(:,:,:),allocatable,target :: atmax,atmin,at00
      real(r8),dimension(:,:,:),allocatable,target :: adv_x,adv_y,adv_z
      real(r8),dimension(:,:,:),allocatable,target :: flxqRtmp,flxqLtmp
!
#ifdef ISO
      real(r8),dimension(imt,jmt,km,NTRA,max_blocks_clinic)::aay_iso,ddy_iso
      real(r8),dimension(imt,jmt,km,NTRA,max_blocks_clinic)::ax_iso,ay_iso,az_iso
      real(r8),dimension(imt,jmt,km,NTRA,max_blocks_clinic)::dx_iso,dy_iso,dz_iso
#endif

#ifdef SST_restore
      real(r8),dimension(imt,jmt,max_blocks_clinic),target::sst_restore
#endif

#ifdef SSS_restore
      real(r8),dimension(imt,jmt,max_blocks_clinic),target::sss_restore
#endif

    ! IAU_FST declare
      real(r8),dimension(imt,jmt,km,NTRA,max_blocks_clinic),target::at_iau
     ! real(r8),allocatable,dimension(:,:,:,:,:),target::at_iau !LPF2020805
      real(r8),allocatable,dimension(:,:,:,:),target::uu_iau,vv_iau !LPF2020805
      real(r8),allocatable,dimension(:,:,:),target::h0_iau !LPF2020805
!
!     ------------------------------------------------------------------
!     Sea Ice
!     ------------------------------------------------------------------
!
      real(r8),dimension(imt,jmt,max_blocks_clinic),public:: licomqice
      real(r8) FW_norm2
      real(r8),dimension(imt,jmt,km,max_blocks_clinic),target ::VTL_ori !for output dt diffusion
      REAL(r8),target :: HDTK(imt,jmt), adv_tt(imt,jmt,km)
      REAL(r8),target    :: h_d2tk(imt,jmt,km)
      real(r8),public :: C2DTTS
      real(r8),dimension(3),target ::d_err_norm2
 !
  contains
  subroutine tracer_init()
      implicit none
      ! real(r8), allocatable :: test(:,:)
      ! integer :: ier
      ! allocate(test(100, 0), stat=ier)

      ! if(ier/=0) write(6,*) "allocate wrong"

!      allocate(atb(imt,jmt,0:km,NTRA,max_blocks_clinic))
      allocate(pdensity(imt,jmt,km,max_blocks_clinic))
      allocate(amld(imt,jmt,max_blocks_clinic))


!       real(r8),dimension(imt,jmt,km,max_blocks_clinic)::penetrate
!       real(r8),dimension(imt,jmt,km,NTRA,max_blocks_clinic)::dt_diff
! !
!       real(r8),dimension(imt,jmt,km,NTRA,max_blocks_clinic)::ddy
!       real(r8),dimension(imt,jmt,km,NTRA,max_blocks_clinic)::dt_conv
!       real(r8),dimension(imt,jmt,km,NTRA,max_blocks_clinic)::dt_restore !20170802
#if (defined LOWRES)
      allocate(tend(imt,jmt,km,NTRA,max_blocks_clinic))
      allocate(ax(imt,jmt,km,NTRA,max_blocks_clinic))
      allocate(ay(imt,jmt,km,NTRA,max_blocks_clinic))
      allocate(az(imt,jmt,km,NTRA,max_blocks_clinic))
      allocate(dx(imt,jmt,km,NTRA,max_blocks_clinic))
!      allocate(dy(imt,jmt,km,NTRA,max_blocks_clinic))
      allocate(dz(imt,jmt,km,NTRA,max_blocks_clinic))
      allocate(dt_diff(imt,jmt,km,NTRA,max_blocks_clinic))
      allocate(dt_conv(imt,jmt,km,NTRA,max_blocks_clinic))
      tend = c0
      ax = c0
      ay = c0
      az = c0
      dx = c0
!      dy = c0
      dt_diff = c0
      dt_conv = c0
#endif

!      allocate(penetrate(imt,jmt,km,max_blocks_clinic))
!      allocate(ddy(imt,jmt,km,NTRA,max_blocks_clinic))
!      allocate(dt_restore(imt,jmt,km,NTRA,max_blocks_clinic))
      allocate(at00(imt,jmt,km), atmax(imt,jmt,km), atmin(imt,jmt,km) )
      allocate(adv_x(imt,jmt,km), adv_y(imt,jmt,km), adv_z(imt,jmt,km) )

      adv_x = c0
      adv_y = c0
      adv_z = c0
#ifdef IAU_FST
!      allocate(at_iau(imt,jmt,km,ntra,max_blocks_clinic))
      allocate(uu_iau(imt,jmt,km,max_blocks_clinic))
      allocate(vv_iau(imt,jmt,km,max_blocks_clinic))
      allocate(h0_iau(imt,jmt,max_blocks_clinic))

      at_iau = c0
      uu_iau = c0
      vv_iau = c0
      h0_iau = c0

#endif
  end subroutine
end module tracer_mod
