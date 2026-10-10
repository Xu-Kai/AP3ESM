!  CVS: $Id: output_mod.F90,v 1.1.1.1 2004/04/29 06:22:39 lhl Exp $
module output_mod
#include <def-undef.h>
!
use precision_mod
use param_mod
!     ------------------------------------------------------------------
!     Output Arrays
!     ------------------------------------------------------------------
#if (defined DAILYACC)


      real(r4),dimension(:,:,:),allocatable::z0daily
      real(r4),dimension(:,:,:),allocatable::mlddaily
      real(r4),dimension(:,:,:),allocatable::lthfdaily,sshfdaily,lwvdaily,swvdaily,freshdaily
      real(r4),dimension(:,:,:),allocatable::precdaily,evapdaily,roffdaily
      real(r4),dimension(:,:,:),allocatable::sudaily,svdaily
      real(r4),dimension(:,:,:,:),allocatable::tsdaily,ssdaily,densdaily
      real(r4),dimension(:,:,:,:),allocatable::usdaily,vsdaily,wsdaily
#if (defined DAILYBUGDET)
      real(r4),allocatable,dimension(:,:,:,:,:)::tenddaily
      real(r4),allocatable,dimension(:,:,:,:,:)::dt_diffdaily,axdaily,aydaily,azdaily
      real(r4),allocatable,dimension(:,:,:,:,:)::dx_diffdaily,dxdaily,dzdaily,dt_convdaily
      !real(r4),dimension(imt,jmt,km,NTRA,max_blocks_clinic)::dxdaily,dzdaily !,dydaily
      real(r4),allocatable,dimension(:,:,:,:)::pendaily
      !real(r4),dimension(imt,jmt,km,max_blocks_clinic)::pendaily
      real(r4),allocatable,dimension(:,:,:,:,:)::dt_restoredaily
      real(r4),allocatable,dimension(:,:,:,:)::net_daily
#endif
#endif

      real(r4),dimension(:,:,:),allocatable::z0mon,sumon,svmon,lthfmon,sshfmon,lwvmon,swvmon,freshmon
      real(r4),dimension(:,:,:,:),allocatable::tsmon,ssmon
      real(r4),dimension(:,:,:,:),allocatable::usmon,vsmon,wsmon,densmon
#ifdef AVG3HOUR
      real(r4),dimension(:,:,:),allocatable::z03hour,su3hour,sv3hour,lthf3hour,sshf3hour,lwv3hour,swv3hour,fresh3hour
      real(r4),dimension(:,:,:,:),allocatable::ts3hour,ss3hour,dens3hour
      real(r4),dimension(:,:,:,:),allocatable::us3hour,vs3hour,ws3hour
#endif
 
#ifdef OBP_AVG3HOUR
      real(r8),dimension(imt,jmt,max_blocks_clinic)::pbo_save
      real(r4),dimension(:,:,:),allocatable::pbo3hour !yuzp20200710
      real(r4),dimension(:,:,:),allocatable::pso3hour !yuzp20200710
      real(r4),dimension(:,:,:),allocatable::ssh3hour
#endif
      real(r4), parameter :: spval =1.0e+35
#if (defined LOWRES)
      real(r4),dimension(imt,jmt,max_blocks_clinic)::himon,hdmon,qicemon
      real(r4),dimension(imt,jmt,2,max_blocks_clinic),target:: icmon
      real(r4),dimension(imt,jmt,NTRA,max_blocks_clinic)::netmon
      real(r4),dimension(imt,jmt,max_blocks_clinic)::runoffmon
      real(r4),dimension(imt,jmt,max_blocks_clinic)::mldmon
      real(r4),dimension(imt,jmt,max_blocks_clinic)::ifracmon !LPF20160817
      real(r4),dimension(imt,jmt,km,max_blocks_clinic)::akmmon,aktmon,aksmon
      real(r4),dimension(imt,jmt,km,NTRA,max_blocks_clinic)::tendmon
      real(r4),dimension(imt,jmt,km,NTRA,max_blocks_clinic)::dt_diffmon
      real(r4),dimension(imt,jmt,km,NTRA,max_blocks_clinic)::axmon,aymon,azmon
      real(r4),dimension(imt,jmt,km,NTRA,max_blocks_clinic)::dxmon,dymon,dzmon
      real(r4),dimension(imt,jmt,km,max_blocks_clinic)::penmon
      real(r4),dimension(imt,jmt,km,NTRA,max_blocks_clinic)::dt_convmon
      !real(r4),dimension(imt,jmt,km,NTRA,max_blocks_clinic)::hdifmon !horizontal diffusion
      real(r4),allocatable, dimension(:,:,:):: psi_euler, psi_eddy, psi
      real(r4),allocatable, dimension(:,:):: mth, mth_adv, mth_adv_iso, mth_dif
#if ( defined TIDEMIX )
      real(r4),dimension(imt,jmt,km,max_blocks_clinic)::aktidemon
      real(r4),dimension(imt,jmt,max_blocks_clinic)::wavedismon
!
      real(r4),dimension(imt,jmt,km,max_blocks_clinic)::richardsonmon     !yuzp-2016/11/13
      real(r4),dimension(imt,jmt,km,max_blocks_clinic)::fztidalmon     !yuzp-2016/11/13
      real(r4),dimension(imt,jmt,km,max_blocks_clinic)::wp3_tidalmon     !yuzp-2016/11/13
      real(r4),dimension(imt,jmt,km,max_blocks_clinic)::ak_tide1mon     !yuzp-2016/11/19
#endif
!#if ( defined CANUTOMIXOUT )
!      real(r4),dimension(imt,jmt,km,max_blocks_clinic)::wp1_canutomon     !yuzp-2016/12/4
!      real(r4),dimension(imt,jmt,km,max_blocks_clinic)::wp2_canutomon     !yuzp-2016/12/4
!      real(r4),dimension(imt,jmt,km,max_blocks_clinic)::wp3_canutomon     !yuzp-2016/12/4
!      real(r4),dimension(imt,jmt,km,max_blocks_clinic)::wp4_canutomon     !yuzp-2016/12/4
!      real(r4),dimension(imt,jmt,km,max_blocks_clinic)::wp5_canutomon     !yuzp-2016/12/4
!      real(r4),dimension(imt,jmt,km,max_blocks_clinic)::wp6_canutomon     !yuzp-2016/12/4
!      real(r4),dimension(imt,jmt,km,max_blocks_clinic)::wp7_canutomon     !yuzp-2016/12/4
!      real(r4),dimension(imt,jmt,km,max_blocks_clinic)::wp8_canutomon     !yuzp-2016/12/4
!
!      real(r4),dimension(imt,jmt,km,max_blocks_clinic)::wk1_canutomon     !yuzp-2016/12/4
!      real(r4),dimension(imt,jmt,km,max_blocks_clinic)::wk2_canutomon     !yuzp-2016/12/4
!      real(r4),dimension(imt,jmt,km,max_blocks_clinic)::wk3_canutomon     !yuzp-2016/12/4
!      real(r4),dimension(imt,jmt,km,max_blocks_clinic)::wk4_canutomon     !yuzp-2016/12/4
!
!      real(r4),dimension(imt,jmt,max_blocks_clinic)::fcor_canutomon     !yuzp-2016/12/4
!      real(r4),dimension(imt,jmt,max_blocks_clinic)::fcort_canutomon     !yuzp-2016/12/4
!
!      real(r4),dimension(imt,jmt,max_blocks_clinic)::wp10_canutomon     !yuzp-2016/12/8
!      real(r4),dimension(imt,jmt,max_blocks_clinic)::wp11_canutomon     !yuzp-2016/12/8
!
!      real(r4),dimension(imt,jmt,km,max_blocks_clinic)::wp12_canutomon     !yuzp-2016/12/4
!      real(r4),dimension(imt,jmt,km,max_blocks_clinic)::wp13_canutomon     !yuzp-2016/12/4
!
!      real(r4),dimension(imt,jmt,km,max_blocks_clinic)::alpha_canutomon     !yuzp-2016/12/4
!      real(r4),dimension(imt,jmt,km,max_blocks_clinic)::beta_canutomon     !yuzp-2016/12/4
!#endif
!
#if (defined ISO_TYPE_BF)
      real(r4),dimension(imt,jmt,km,max_blocks_clinic)::athkdfmon
#endif
#ifdef ISO
      real(r4),dimension(imt,jmt,km,NTRA,max_blocks_clinic)::axmon_iso,aymon_iso,azmon_iso
      real(r4),dimension(imt,jmt,km,NTRA,max_blocks_clinic)::dxmon_iso,dymon_iso,dzmon_iso
      real(r4),dimension(imt,jmt,km,max_blocks_clinic):: wsmon_iso,vsmon_iso
#if ( defined ISOOUT )
      real(r4),dimension(imt,jmt,km,max_blocks_clinic)::vetisomon,vntisomon,vbtisomon
#endif
#endif
!
!#if (defined SMAG_OUT)
!      real(r4),dimension(imt,jmt,km)::a3mon
!#endif
!
      real(r8) ERR_norm2mon !LPF20160823
#endif
!for LOWREN
contains
subroutine output_init()
      implicit none


#if (defined DAILYACC)
      allocate(z0daily(imt,jmt,max_blocks_clinic))
      allocate(mlddaily(imt,jmt,max_blocks_clinic))
      allocate(lthfdaily(imt,jmt,max_blocks_clinic))
      allocate(sshfdaily(imt,jmt,max_blocks_clinic))
      allocate(lwvdaily(imt,jmt,max_blocks_clinic))
      allocate(swvdaily(imt,jmt,max_blocks_clinic))
      allocate(freshdaily(imt,jmt,max_blocks_clinic))
      allocate(precdaily(imt,jmt,max_blocks_clinic))
      allocate(evapdaily(imt,jmt,max_blocks_clinic))
      allocate(roffdaily(imt,jmt,max_blocks_clinic))
      allocate(sudaily(imt,jmt,max_blocks_clinic))
      allocate(svdaily(imt,jmt,max_blocks_clinic))      
      allocate(tsdaily(imt,jmt,km,max_blocks_clinic))
      allocate(ssdaily(imt,jmt,km,max_blocks_clinic))
      allocate(usdaily(imt,jmt,km,max_blocks_clinic))
      allocate(vsdaily(imt,jmt,km,max_blocks_clinic))
      allocate(densdaily(imt,jmt,km,max_blocks_clinic))
      allocate(wsdaily(imt,jmt,kmp1,max_blocks_clinic))
#endif
#if (defined DAILYACC || defined MONTHACC)
! for monthly data output
      allocate(z0mon(imt,jmt,max_blocks_clinic))
      allocate(sumon(imt,jmt,max_blocks_clinic))
      allocate(svmon(imt,jmt,max_blocks_clinic))
      allocate(lthfmon(imt,jmt,max_blocks_clinic))
      allocate(sshfmon(imt,jmt,max_blocks_clinic))
      allocate(lwvmon(imt,jmt,max_blocks_clinic))
      allocate(swvmon(imt,jmt,max_blocks_clinic))
      allocate(freshmon(imt,jmt,max_blocks_clinic))

      allocate(tsmon(imt,jmt,km,max_blocks_clinic))
      allocate(ssmon(imt,jmt,km,max_blocks_clinic))
      allocate(usmon(imt,jmt,km,max_blocks_clinic))
      allocate(vsmon(imt,jmt,km,max_blocks_clinic))
      allocate(wsmon(imt,jmt,kmp1,max_blocks_clinic))
      allocate(densmon(imt,jmt,km,max_blocks_clinic))

#endif
#ifdef AVG3HOUR
! for 3hour average data output
      allocate(z03hour(imt,jmt,max_blocks_clinic))
      allocate(su3hour(imt,jmt,max_blocks_clinic))
      allocate(sv3hour(imt,jmt,max_blocks_clinic))
      allocate(lthf3hour(imt,jmt,max_blocks_clinic))
      allocate(sshf3hour(imt,jmt,max_blocks_clinic))
      allocate(lwv3hour(imt,jmt,max_blocks_clinic))
      allocate(swv3hour(imt,jmt,max_blocks_clinic))
      allocate(fresh3hour(imt,jmt,max_blocks_clinic))

      allocate(ts3hour(imt,jmt,km,max_blocks_clinic))
      allocate(ss3hour(imt,jmt,km,max_blocks_clinic))
      allocate(us3hour(imt,jmt,km,max_blocks_clinic))
      allocate(vs3hour(imt,jmt,km,max_blocks_clinic))
      allocate(dens3hour(imt,jmt,km,max_blocks_clinic))
      allocate(ws3hour(imt,jmt,kmp1,max_blocks_clinic))
#endif
#ifdef OBP_AVG3HOUR
      allocate(pbo3hour(imt,jmt,max_blocks_clinic))
      allocate(pso3hour(imt,jmt,max_blocks_clinic))
      allocate(ssh3hour(imt,jmt,max_blocks_clinic))
#endif
end subroutine


end module output_mod
