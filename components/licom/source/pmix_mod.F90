!  CVS: $Id: pmix_mod.F90,v 1.1.1.1 2004/04/29 06:22:39 lhl Exp $
module pmix_mod
#include <def-undef.h>
!
use precision_mod
use param_mod
!     Calculating Richardson number riu (at U/V-point);
!                                   rit (at T-point)
!     in order to computing vertical mixing coefficients based on
!     Pacanowski & Philander (JPO vol 11, #11, 1981).
!
!     Note: this parameterization was designed for equatorial models
!     and may not do a good job in mid or high latitudes. Simulations
!     in these regions (where vertical shear is small) are improved with
!     the addition of solar short wave penetration into the ocean which
!     reduces buoyancy and enhances vertical mixing.
!
      INTEGER:: RTST,RTEND,RUST,RUEND
!
      real(r8):: wndmix, fricmx, diff_cbt_back, diff_cbt_limit,&
                     visc_cbu_back, visc_cbu_limit
!lhl1204 real(r8),dimension(:,:,:),allocatable:: ric,rit
      real(r8),dimension(:,:,:,:),allocatable,target :: ric,rict,rict_replace
      real(r8),dimension(:,:,:),allocatable:: rict_ref
      real(r8),dimension(:,:,:,:),allocatable,target:: rit,riu
!      real(r8),dimension(imt,jmt,kmm1):: ricdt,ricdu
      !real(r8),dimension(imt,jmt,kmm1,max_blocks_clinic):: ricdt
      !real(r8),dimension(imt,jmt,kmm1,max_blocks_clinic):: ricdtTmS !LPF20160728
      real(r8),dimension(:,:,:,:), allocatable,target:: ricdt
      real(r8),dimension(:,:,:,:), allocatable,target:: ricdtTmS !LPF20160728
!      real(r8),dimension(imt,jmt,kmm1):: ridu,ridt,s2u,s2t
      !real(r8),dimension(imt,jmt,kmm1,max_blocks_clinic):: ridt,s2u,s2t
      real(r8),dimension(:,:,:,:), allocatable,target:: ridt,s2u,s2t
!lhl1204
!
!
!
!
!
!
!-----------------------------------------------------------------------
!     Solar Shortwave energy penetrates below the ocean surface. Clear
!     water assumes energy partitions between two exponentials as
!     follows:
!
!     58% of the energy decays with a 35 cm e-folding scale
!     42% of the energy decays with a 23 m e-folding scale
!
!     if the thickness of the first ocean level "dzt(1)" is 50 meters,
!     then shortwave penetration wouldn't matter. however, for
!     dzt(1) = 10 meters, the effect can be significant. this can be
!     particularly noticable in the summer hemisphere.
!
!     Paulson and Simpson (1977 Irradiance measurements in the upper
!                               ocean JPO 7, 952-956)
!     Also see ... Jerlov (1968 Optical oceanography. Elsevier)
!                  A General Circulation Model for Upper Ocean
!                  Simulaton (Rosati and Miyakoda JPO vol 18,Nov 1988)
!   based on chlorophyll concentration shortwave radiation
!   can see Ohlmann Ocean Radiant Heating in Climate models( JC 5 2003 )
!-----------------------------------------------------------------------
!
#if (defined SOLAR)
      REAL(r8):: pen(kmm1)
     ! REAL(r8), allocatable:: pen(:)
#endif

#if (defined SOLARCHLORO)
      real(r8),dimension(imt,jmt,km,max_blocks_clinic),target:: pen_chl   !  be different from pen(kmm1)
      !real(r8),dimension(:,:,:,:), allocatable:: pen_chl   !  be different from pen(kmm1)
#endif

contains
      subroutine pmix_init()
            implicit none
            !real(r8),dimension(imt,jmt,kmm1,max_blocks_clinic):: ricdt
            !real(r8),dimension(imt,jmt,kmm1,max_blocks_clinic):: ricdtTmS !LPF20160728
            !real(r8),dimension(imt,jmt,kmm1,max_blocks_clinic):: ridt,s2u,s2t
            !real(r8),dimension(imt,jmt,km,max_blocks_clinic):: pen_chl   !  be different from pen(kmm1)
            allocate(ricdt(imt,jmt,kmm1,max_blocks_clinic))
            allocate(ricdtTmS(imt,jmt,kmm1,max_blocks_clinic))
!            allocate(ridt(imt,jmt,kmm1,max_blocks_clinic))
            allocate(s2u(imt,jmt,kmm1,max_blocks_clinic))
            allocate(s2t(imt,jmt,kmm1,max_blocks_clinic))

            allocate(rit(imt,jmt,kmm1,max_blocks_clinic),ric(imt,jmt,kmm1,max_blocks_clinic))
            allocate(rict(imt,jmt,kmm1,max_blocks_clinic))
            allocate(riu(imt,jmt,0:km,max_blocks_clinic),stat=ierr)

            rit   = 0.0_r8
            ric   = 0.0_r8
            rict  = 0.0_r8
            ricdtTmS  = 0.0_r8
            ricdt = 0.0_r8
            riu = 0.0_r8
#if (defined LOWRES)
            allocate(rict_replace(imt,jmt,kmm1,max_blocks_clinic))
            allocate(rict_ref(imt,jmt,max_blocks_clinic))
            rict_replace = 0.0_r8
            rict_ref = 0.0_r8
#endif
!#if (defined SOLAR)
!        allocate(pen(kmm1))
!#endif
!#if (defined SOLARCHLORO)
!            allocate(pen_chl(imt,jmt,km,max_blocks_clinic))
!#endif
      end subroutine



end module pmix_mod
