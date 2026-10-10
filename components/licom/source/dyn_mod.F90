!  CVS: $Id: dyn_mod.F90,v 1.1.1.1 2004/04/29 06:22:39 lhl Exp $
module dyn_mod
#include <def-undef.h>
use precision_mod
use param_mod
!     ------------------------------------------------------------------
!     U V T S H0 W RHO
!     ------------------------------------------------------------------
!       real(r8),dimension(imt,jmt,max_blocks_clinic)::ub,vb,ubp,vbp,h0p
!       real(r8),dimension(imt,jmt,km,max_blocks_clinic)::up,vp
       real(r8),dimension(imt,jmt,kmp1,max_blocks_clinic),target ::ws
       real(r8),dimension(imt,jmt,max_blocks_clinic),target ::h0l,h0f,h0bl,h0bf
       real(r8),dimension(imt,jmt,km,max_blocks_clinic),target ::utl,utf,vtl,vtf
       REAL(r8),dimension(IMT,JMT,KM,MAX_BLOCKS_CLINIC),target :: WP12,WP13,riv1,riv2
       real(r8),target :: u_wface(imt,jmt,km), v_sface(imt,jmt,km)
       real(r8),target :: d2ukm(imt,jmt,km), d2vkm(imt,jmt,km), am_factor_m(imt,jmt,km)
      !real(r8),dimension(:,:,:),allocatable,target ::hdukkm,hdvkkm,am_f_km
! !
!       REAL(r8)    :: SBCX(imt,jmt,max_blocks_clinic),BBCX(imt,jmt,max_blocks_clinic)
!       REAL(r8)    :: SBCY(imt,jmt,max_blocks_clinic),BBCY(imt,jmt,max_blocks_clinic)
      real(r8),dimension(:,:,:), allocatable,target::ub,vb,ubp,vbp,h0p
      real(r8),dimension(:,:,:,:), allocatable,target::up,vp
      real(r8),dimension(:,:,:,:), allocatable,target::wsbak,wp3tmp,wk1tmp,wk2tmp,wk3tmp
!      real(r8),dimension(:,:,:), allocatable,target::h0l,h0f,h0bl,h0bf
!      real(r8),dimension(:,:,:,:), allocatable ,target::utl,utf,vtl,vtf
!
      integer,dimension(:,:,:), allocatable, target::mldlev

      REAL(r8), allocatable,target :: SBCX(:,:,:),BBCX(:,:,:)
      REAL(r8), allocatable,target :: SBCY(:,:,:),BBCY(:,:,:)
#ifdef ISOOUT
      real(r8),dimension(imt,jmt,km,max_blocks_clinic)::vetiso,vntiso,vbtiso
#endif      
!
      real(r8),allocatable,dimension(:,:) :: buffer
!
      real(r8),allocatable,dimension(:,:,:),target ::h0
      real(r8),allocatable,dimension(:,:,:,:),target ::u,v
#ifdef TIDEFORCE
      real(r8),allocatable,dimension(:,:,:),target :: elevat 
      real(r8),allocatable,dimension(:,:,:),target :: h0_tide 
#endif
!
!
!     ------------------------------------------------------------------
!     Pressure gradient
!     ------------------------------------------------------------------
      real(r8),dimension(:,:,:,:),allocatable,target ::gg,dlu,dlv
      real(r8),dimension(:,:,:),allocatable,target ::dlub,dlvb

      real(r8),dimension(:,:,:,:),allocatable,target :: alpha, beta, pp, ppa, ppb, ppc
      real(r8),dimension(:,:,:),allocatable,target:: work1, work2
      REAL (r8),target :: dt2k_b(imt,jmt)
      real(r8),target, dimension(imt,jmt,km,max_blocks_clinic)::wka1, wka2
      real(r8),target, dimension(imt,jmt,km+1,max_blocks_clinic)::wktmp
      REAL(r8),target :: ucheck(imt,jmt,km), vcheck(imt,jmt,km)!ZWP20170212

      !integer,dimension(:),allocatable,target :: bit_3d_ocean_point
      integer :: num_ocean_point
contains
      subroutine dyn_init()

      implicit none

!       real(r8),dimension(imt,jmt,max_blocks_clinic)::ub,vb,ubp,vbp,h0p
!       real(r8),dimension(imt,jmt,km,max_blocks_clinic)::up,vp
!       real(r8),dimension(imt,jmt,kmp1,max_blocks_clinic)::ws
!       real(r8),dimension(imt,jmt,max_blocks_clinic)::h0l,h0f,h0bl,h0bf
!       real(r8),dimension(imt,jmt,km,max_blocks_clinic)::utl,utf,vtl,vtf
! !
!       REAL(r8)    :: SBCX(imt,jmt,max_blocks_clinic),BBCX(imt,jmt,max_blocks_clinic)
!       REAL(r8)    :: SBCY(imt,jmt,max_blocks_clinic),BBCY(imt,jmt,max_blocks_clinic)

      allocate(ub(imt, jmt, max_blocks_clinic))
      allocate(vb(imt, jmt, max_blocks_clinic))
      allocate(ubp(imt, jmt, max_blocks_clinic))
      allocate(vbp(imt, jmt, max_blocks_clinic))
      allocate(h0p(imt, jmt, max_blocks_clinic))
      allocate(up(imt, jmt, km, max_blocks_clinic))
      allocate(vp(imt, jmt, km, max_blocks_clinic))
!      allocate(ws(imt, jmt, kmp1, max_blocks_clinic))
      allocate(wsbak(imt, jmt, kmp1, max_blocks_clinic))
!      allocate(wp3tmp(imt, jmt, km, max_blocks_clinic))
!      allocate(wk1tmp(imt, jmt, km, max_blocks_clinic))
!      allocate(wk2tmp(imt, jmt, km, max_blocks_clinic))
!      allocate(wk3tmp(imt, jmt, km, max_blocks_clinic))
      allocate(mldlev(imt,jmt,max_blocks_clinic))
!      allocate(h0l(imt, jmt, max_blocks_clinic))
 !     allocate(h0f(imt, jmt, max_blocks_clinic))
 !     allocate(h0bl(imt, jmt, max_blocks_clinic))
 !     allocate(h0bf(imt, jmt, max_blocks_clinic))
 !     allocate(utl(imt, jmt, km, max_blocks_clinic))
 !     allocate(utf(imt, jmt, km, max_blocks_clinic))
 !     allocate(vtl(imt, jmt, km, max_blocks_clinic))
 !     allocate(vtf(imt, jmt, km, max_blocks_clinic))

      allocate(SBCX(imt, jmt, max_blocks_clinic))
      allocate(BBCX(imt, jmt, max_blocks_clinic))
      allocate(SBCY(imt, jmt, max_blocks_clinic))
      allocate(BBCY(imt, jmt, max_blocks_clinic))

      allocate(dlu(imt,jmt,km,max_blocks_clinic),dlv(imt,jmt,km,max_blocks_clinic),gg(imt,jmt,km,max_blocks_clinic))
      SBCX = 0.0D0
      BBCX = 0.0D0
      SBCY = 0.0D0
      BBCY = 0.0D0

      allocate ( alpha(imt,jmt,km,max_blocks_clinic), beta(imt,jmt,km,max_blocks_clinic))
      allocate ( pp(imt,jmt,km,max_blocks_clinic), ppa(imt,jmt,km,max_blocks_clinic))
      allocate ( ppb(imt,jmt,km,max_blocks_clinic), ppc(imt,jmt,km,max_blocks_clinic))
      allocate ( work1(imt,jmt,max_blocks_clinic), work2(imt,jmt,max_blocks_clinic))
      allocate(dlub(imt,jmt,max_blocks_clinic),dlvb(imt,jmt,max_blocks_clinic),stat=ierr)
      !allocate(bit_3d_ocean_point(imt*jmt))

      !allocate(hdukkm(imt,jmt,km),hdvkkm(imt,jmt,km),am_f_km(imt,jmt,km))
      !hdukkm = 0.0D0
      !hdvkkm = 0.0D0
      !am_f_km = 0.0D0
      wsbak = 0.0D0
      end subroutine
end module dyn_mod

