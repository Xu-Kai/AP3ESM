!  CVS: $Id: work_mod.F90,v 1.1.1.1 2004/04/29 06:22:39 lhl Exp $
module work_mod
#include <def-undef.h>
!
use precision_mod
use param_mod
!
!
!     ------------------------------------------------------------------
!     Working Arrays
!     ------------------------------------------------------------------
      !real(r8),dimension(imt,jmt,max_blocks_clinic):: PXB,PYB,PAX,PAY,WHX,WHY,WGP
      !real(r8),dimension(imt,jmt,km,max_blocks_clinic)::wka
      real(r8),dimension(:,:,:),allocatable,target:: PXB,PYB,PAX,PAY,WHX,WHY,WGP
      real(r8),dimension(:,:,:,:),allocatable,target::wka,wkz
      real(r8),dimension(:,:,:,:),allocatable:: work_1,work_2,work_3,temp
      real(r8),dimension(:,:,:,:),allocatable,target:: tmp1,tmp2,uk,vk
      real(r8),dimension(imt,jmt,max_blocks_clinic),target :: work
!      real(r8),dimension(:,:,:), allocatable,target :: work
      !real(r8):: wkk(kmp1)
      real(r8), allocatable:: wkk(:)
      real(r8),dimension(imt,jmt,km,max_blocks_clinic),target:: WKB,WKC,WKD,TF
!      real(r8),dimension(:,:,:,:),allocatable,target:: WKB,WKC,WKD,TF
      real(r8),dimension(imt,jmt,max_blocks_clinic),target::stf
!      real(r8),dimension(:,:,:),allocatable,target::stf
      real(r4),dimension(:,:),allocatable:: buffer_real4
      real(r8),dimension(:,:), allocatable :: work1_g, work2_g, work3_g
      
      real (r8), dimension (:,:,:), allocatable,target :: &
      DTN,DTS_s,DTE,DTW,                          &
      DUC,DUN,DUS,DUE,DUW,                      &
      DMC,DMN,DMS,DME,DMW,                      &
      DUM,                                      &
      AHF,              &! variable mixing factor for tracer   mixing
      AMF,              &! variable mixing factor for momentum mixing
      ratio_dxy         !changing mixing !YYQ add
      
      real (r8), dimension (nx_block,ny_block),target :: D2UK,D2VK        ! intermediate Del**2 results
      
      real(r8), target  :: z_curl(nx_block,ny_block,max_blocks_clinic)
      real(r8), target  :: diverg(nx_block,ny_block,max_blocks_clinic)
      real(r8), target  :: am_factor(nx_block,ny_block,max_blocks_clinic)
      real(r8), target  :: am_factor_d(nx_block,ny_block,max_blocks_clinic)
      real(r8), target  :: factor_vis_tsz0(imt,jmt,max_blocks_clinic)   ! lat-profie for ss/tt/ssh diffusivity, yjf&xjw 20250724
      real(r8), target  :: factor_vis_uv(imt,jmt,max_blocks_clinic)   ! lat-profie for horizontal viscosity, yjf&xjw 20251029
!      real(r8) :: factor_vis_tsz0_global(imt_global,jmt_global)

contains
subroutine work_init()
      implicit none
      allocate(PXB(imt,jmt,max_blocks_clinic))
      allocate(PYB(imt,jmt,max_blocks_clinic))
      allocate(PAX(imt,jmt,max_blocks_clinic))
      allocate(PAY(imt,jmt,max_blocks_clinic))
      allocate(WHX(imt,jmt,max_blocks_clinic))
      allocate(WHY(imt,jmt,max_blocks_clinic))
      allocate(WGP(imt,jmt,max_blocks_clinic))
      allocate(wka(imt,jmt,km,max_blocks_clinic))
      allocate(wkz(imt,jmt,km,max_blocks_clinic))
!      allocate(work(imt,jmt,max_blocks_clinic))
      allocate(wkk(kmp1))
      allocate(uk(imt,jmt,km,max_blocks_clinic))
      allocate(vk(imt,jmt,km,max_blocks_clinic))
      z_curl = 0.0
      diverg = 0.0
      am_factor = 0.0
      am_factor_d = 0.0
      wkz = 0.0
!      allocate(stf(imt,jmt,max_blocks_clinic))
 !     allocate(tf(imt,jmt,km,max_blocks_clinic))
 !     allocate(wkb(imt,jmt,km,max_blocks_clinic))
 !     allocate(wkc(imt,jmt,km,max_blocks_clinic))
 !     allocate(wkd(imt,jmt,km,max_blocks_clinic))
end subroutine


!
end module work_mod
