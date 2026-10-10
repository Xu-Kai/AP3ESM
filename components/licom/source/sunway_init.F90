!  CVS: $Id: param_mod.F90,v 1.1.1.1 2004/04/29 06:22:39 lhl Exp $
module sunway_init
#include <def-undef.h>
use precision_mod
use param_mod,only:imt,jmt,km,mytid
use sunway_mod
!-------------------------------------------------------------------------------
! Author: Maoxue YU  ( 12 July, 2024) !
!-------------------------------------------------------------------------------


    real(r8),dimension(11),target::turb_param
contains

subroutine readyt_sw_init

    use pconst_mod
    use tracer_mod
    use dyn_mod
    use work_mod
    use pmix_mod
    use constant_mod
    use forc_mod, only: psa,USTAR,BUOYTUR, BUOYSOL,NSWV,SWV
    use grid
    implicit none

    readyt_cp_param%d_dlu        => DLU(1,1,1,1)
    readyt_cp_param%d_dlv        => DLV(1,1,1,1)
    readyt_cp_param%d_rit        => rit(1,1,1,1)
    readyt_cp_param%d_ric        => RIC(1,1,1,1)
    readyt_cp_param%d_rict       => RICT(1,1,1,1)
    readyt_cp_param%d_ricdtTms   => ricdtTmS(1,1,1,1)
    readyt_cp_param%d_ricdt      => RICDT(1,1,1,1)
    readyt_cp_param%d_akt        => AKT(1,1,1,1,1)
    readyt_cp_param%rank         = my_task 

    mpe_param%cuBlock(1)     = imt
    mpe_param%cuBlock(2)     = 1
    mpe_param%cuBlock(3)     = 1
    mpe_param%cuGrid(1)      = 1
    mpe_param%cuGrid(2)      = jmt
    mpe_param%cuGrid(3)      = km
    mpe_param%rank           =  my_task
    mpe_param%d_h0l         => H0L(1,1,1)
    mpe_param%d_h0f         => H0F(1,1,1)
    mpe_param%d_h0          => H0(1,1,1)
    mpe_param%d_utl         => UTL(1,1,1,1)
    mpe_param%d_vtl         => VTL(1,1,1,1)
    mpe_param%d_utf         => UTF(1,1,1,1)
    mpe_param%d_vtf         => VTF(1,1,1,1)
    mpe_param%d_u           => U(1,1,1,1)
    mpe_param%d_v           => V(1,1,1,1)
    mpe_param%d_dlu         => DLU(1,1,1,1)
    mpe_param%d_dlv         => DLV(1,1,1,1)
    mpe_param%d_atb         => ATB(1,1,0,1,1)
!    mpe_param%d_au0         => AU0(1,1,1)
!    mpe_param%d_aus         => AUS(1,1,1)
!    mpe_param%d_auw         => AUW(1,1,1)
!    mpe_param%d_ausw        => AUSW(1,1,1)


    dens_para%cuBlock(1)     = imt
    dens_para%cuBlock(2)     = 1
    dens_para%cuBlock(3)     = 1
    dens_para%cuGrid(1)      = 1
    dens_para%cuGrid(2)      = jmt
    dens_para%cuGrid(3)      = kmm1
    dens_para%rank           =  my_task
    dens_para%d_km           =  km
    dens_para%d_kmm1         =  kmm1
    dens_para%d_jmt          =  jmt
    dens_para%d_imt          =  imt
    dens_para%d_g            =  g
    dens_para%d_od0          =  od0
    dens_para%d_dlu          => DLU(1,1,1,1)
    dens_para%d_dlv          => DLV(1,1,1,1)
    dens_para%d_to           => TO(1)
    dens_para%d_so           => SO(1)
    dens_para%d_ric          => RIC(1,1,1,1)
    dens_para%d_viv          => VIV(1,1,1,1)
    dens_para%d_odzt         => ODZT(1)
    dens_para%d_c            => C(1,1)
    dens_para%d_rict         => RICT(1,1,1,1)
    dens_para%d_vit          => VIT(1,1,1,1)
    dens_para%d_atb          => ATB(1,1,0,1,1)

 mpe_2_param%cuBlock(1)     = imt
    mpe_2_param%cuBlock(2)     = 1
    mpe_2_param%cuBlock(3)     = 1
    mpe_2_param%cuGrid(1)      = 1
    mpe_2_param%cuGrid(2)      = jmt
    mpe_2_param%cuGrid(3)      = km
    mpe_2_param%d_km           =  km
    mpe_2_param%d_jmt          =  jmt
    mpe_2_param%d_imt          =  imt
    mpe_2_param%rank           =  my_task
    mpe_2_param%d_g            =  g
    mpe_2_param%d_od0          =  od0
    mpe_2_param%d_to           => TO(1)
    mpe_2_param%d_so           => SO(1)
    mpe_2_param%d_c            => C(1,1)
    mpe_2_param%d_po           => PO(1)
    mpe_2_param%d_gg           => GG(1,1,1,1)
    mpe_2_param%d_pdensity     => PDENSITY(1,1,1,1)
    mpe_2_param%d_vit          => VIT(1,1,1,1)
    mpe_2_param%d_ppb          => PPB(1,1,1,1)
    mpe_2_param%d_ppc          => PPC(1,1,1,1)
    mpe_2_param%d_at           => AT(1,1,1,1,1)
    mpe_2_param%d_dzp          => DZP(1)

    mpe_2_2_param%cuBlock(1)     = imt
    mpe_2_2_param%cuBlock(2)     = 1
    mpe_2_2_param%cuBlock(3)     = 1
    mpe_2_2_param%cuGrid(1)      = 1
    mpe_2_2_param%cuGrid(2)      = jmt
    mpe_2_2_param%cuGrid(3)      = 1
    mpe_2_2_param%d_km           =  km
    mpe_2_2_param%d_jmt          =  jmt
    mpe_2_2_param%d_imt          =  imt
    mpe_2_2_param%rank           =  my_task
    mpe_2_2_param%d_od0          =  od0
    mpe_2_2_param%d_psa          => PSA(1,1,1)
    mpe_2_2_param%d_ppb          => PPB(1,1,1,1)
    mpe_2_2_param%d_ppc          => PPC(1,1,1,1)
    mpe_2_2_param%d_at           => AT(1,1,1,1,1)
    mpe_2_2_param%d_gg           => GG(1,1,1,1)
    mpe_2_2_param%d_vit          => VIT(1,1,1,1)
    mpe_2_2_param%d_pp           => PP(1,1,1,1)
    mpe_2_2_param%d_ppa          => PPA(1,1,1,1)
    mpe_2_2_param%d_dzp          => DZP(1)

    thermal_param%cuBlock(1)     = imt
    thermal_param%cuBlock(2)     = 1
    thermal_param%cuBlock(3)     = 1
    thermal_param%cuGrid(1)      = 1
    thermal_param%cuGrid(2)      = jmt
    thermal_param%cuGrid(3)      = km
    thermal_param%d_km           =  km
    thermal_param%d_jmt          =  jmt
    thermal_param%d_imt          =  imt
    thermal_param%rank           =  my_task
    thermal_param%d_od0          =  od0
    thermal_param%d_ppb          => PPB(1,1,1,1)
    thermal_param%d_ppc          => PPC(1,1,1,1)
    thermal_param%d_ppa          => PPA(1,1,1,1)
    thermal_param%d_alpha        => alpha(1,1,1,1)
    thermal_param%d_beta         => beta(1,1,1,1)
    thermal_param%d_vit          => VIT(1,1,1,1)


    readyt_3_parm%cuBlock(1)     = imt
    readyt_3_parm%cuBlock(2)     = 1
    readyt_3_parm%cuBlock(3)     = 1
    readyt_3_parm%cuGrid(1)      = 1
    readyt_3_parm%cuGrid(2)      = jmt
    readyt_3_parm%cuGrid(3)      = km
    readyt_3_parm%d_kmm1         =  kmm1
    readyt_3_parm%d_km           =  km
    readyt_3_parm%d_jmt          =  jmt
    readyt_3_parm%d_imt          =  imt
    readyt_3_parm%rank           =  my_task
    readyt_3_parm%d_kmu          =>  KMU(1,1,1)
    readyt_3_parm%d_od0cp        =  od0cp
    readyt_3_parm%d_od0          =  od0
    readyt_3_parm%d_g            =  g
    readyt_3_parm%d_ricdttms     => ricdtTmS(1,1,1,1)
    readyt_3_parm%d_ricdt        => RICDT(1,1,1,1)
    readyt_3_parm%d_gg           => GG(1,1,1,1)
    readyt_3_parm%d_pdensity     => PDENSITY(1,1,1,1)
    readyt_3_parm%d_vit          => VIT(1,1,1,1)
    readyt_3_parm%d_pp           => PP(1,1,1,1)
    readyt_3_parm%d_at           => AT(1,1,1,1,1)
    readyt_3_parm%d_po           => PO(1)
    readyt_3_parm%d_dlu          => DLU(1,1,1,1)
    readyt_3_parm%d_dlv          => DLV(1,1,1,1)
    readyt_3_parm%d_dxur         => DXUR(1,1,1)
    readyt_3_parm%d_dyur         => DYUR(1,1,1)
    readyt_3_parm%d_buoytur      => BUOYTUR(1,1,1)
    readyt_3_parm%d_buoysol      => BUOYSOL(1,1,1)
    readyt_3_parm%d_nswv         => NSWV(1,1,1)
    readyt_3_parm%d_swv          => SWV(1,1,1)
    readyt_3_parm%d_alpha        => ALPHA(1,1,1,1)
    readyt_3_parm%d_beta         => BETA(1,1,1,1)
    readyt_3_parm%d_odzt         => ODZT(1)
    readyt_3_parm%d_akt          =>AKT(1,1,1,1,1)
#if (defined CANUTOMIXOUT)
    readyt_3_parm%d_alpha_canuto => alpha_canuto(1,1,1,1)
    readyt_3_parm%d_beta_canuto  => beta_canuto(1,1,1,1)
#endif

    readyt_3_2_parm%cuBlock(1)     = imt
    readyt_3_2_parm%cuBlock(2)     = 1
    readyt_3_2_parm%cuBlock(3)     = 1
    readyt_3_2_parm%cuGrid(1)      = 1
    readyt_3_2_parm%cuGrid(2)      = jmt
    readyt_3_2_parm%cuGrid(3)      = 1
    readyt_3_2_parm%d_km           =  km
    readyt_3_2_parm%d_jmt          =  jmt
    readyt_3_2_parm%d_imt          =  imt
    readyt_3_2_parm%rank           =  my_task
    readyt_3_2_parm%d_gg           => GG(1,1,1,1)
    readyt_3_2_parm%d_dzp          => DZP(1)
    readyt_3_2_parm%d_dlu          => DLU(1,1,1,1)
    readyt_3_2_parm%d_dlv          => DLV(1,1,1,1)
    readyt_3_2_parm%d_ohbu         => OHBU(1,1,1)
    readyt_3_2_parm%d_ohbt         => OHBT(1,1,1)
    readyt_3_2_parm%d_zkt          => ZKT(1)
    readyt_3_2_parm%d_viv          => VIV(1,1,1,1)
    readyt_3_2_parm%d_pxb          => PXB(1,1,1)
    readyt_3_2_parm%d_pyb          => PYB(1,1,1)

    readyt_4_parm%cuBlock(1)     = imt
    readyt_4_parm%cuBlock(2)     = 1
    readyt_4_parm%cuBlock(3)     = 1
    readyt_4_parm%cuGrid(1)      = 1
    readyt_4_parm%cuGrid(2)      = jmt
    readyt_4_parm%cuGrid(3)      = 1
    readyt_4_parm%d_km           =  km
    readyt_4_parm%d_jmt          =  jmt
    readyt_4_parm%d_imt          =  imt
    readyt_4_parm%rank           =  my_task
    readyt_4_parm%d_g            =  g
    readyt_4_parm%d_dlu          => DLU(1,1,1,1)
    readyt_4_parm%d_dlv          => DLV(1,1,1,1)
    readyt_4_parm%d_ohbt         => OHBT(1,1,1)


    readyt_5_parm%cuBlock(1)     = imt
    readyt_5_parm%cuBlock(2)     = 1
    readyt_5_parm%cuBlock(3)     = 1
    readyt_5_parm%cuGrid(1)      = 1
    readyt_5_parm%cuGrid(2)      = jmt
    readyt_5_parm%cuGrid(3)      = 1
    readyt_5_parm%d_jmt          =  jmt
    readyt_5_parm%d_imt          =  imt
    readyt_5_parm%rank           =  my_task
    readyt_5_parm%d_kmu          =>  KMU(1,1,1)
    readyt_5_parm%d_od0          =  od0
    readyt_5_parm%d_dlv          => DLV(1,1,1,1)
    readyt_5_parm%d_au0          => AU0(1,1,1)
    readyt_5_parm%d_aus          => AUS(1,1,1)
    readyt_5_parm%d_auw          => AUW(1,1,1)
    readyt_5_parm%d_ausw         => AUSW(1,1,1)
    readyt_5_parm%d_wgp          => WGP(1,1,1)
    readyt_5_parm%d_work         => WORK(1,1,1)
    readyt_5_parm%d_hbx          => HBX(1,1,1)
    readyt_5_parm%d_hby          => HBY(1,1,1)
    readyt_5_parm%d_viv          => VIV(1,1,1,1)
    readyt_5_parm%d_whx          => WHX(1,1,1)
    readyt_5_parm%d_why          => WHY(1,1,1)
    readyt_5_parm%d_work1        => WORK1(1,1,1)
    readyt_5_parm%d_work2        => WORK2(1,1,1)
    readyt_5_parm%d_psa          => PSA(1,1,1)
    readyt_5_parm%d_dxur         => DXUR(1,1,1)
    readyt_5_parm%d_dyur         => DYUR(1,1,1)
    readyt_5_parm%d_pax          => PAX(1,1,1)
    readyt_5_parm%d_pay          => PAY(1,1,1)

end subroutine readyt_sw_init


subroutine readyc_sw_init

    use pconst_mod
    use tracer_mod
    use dyn_mod
    use work_mod
    use pmix_mod
    use constant_mod
!    use forc_mod, only: psa,USTAR,BUOYTUR, BUOYSOL,NSWV,SWV
    use forc_mod
    use grid
    use blocks
    use domain
    implicit none
    type(block):: this_block
    integer :: iblock,biti,bitj,bitjj


    readyc1_param%cuBlock(1)     = imt
    readyc1_param%cuBlock(2)     = 1
    readyc1_param%cuBlock(3)     = 1
    readyc1_param%cuGrid(1)      = 1
    readyc1_param%cuGrid(2)      = jmt
    readyc1_param%cuGrid(3)      = km
    readyc1_param%d_h0bl        => H0BL(1,1,1)
    readyc1_param%d_h0bf        => H0BF(1,1,1)
    readyc1_param%d_h0          => H0(1,1,1)
    readyc1_param%d_up          => UP(1,1,1,1)
    readyc1_param%d_vp          => VP(1,1,1,1)
    readyc1_param%d_viv         => VIV(1,1,1,1)
    readyc1_param%d_vit         => VIT(1,1,1,1)
    readyc1_param%d_wp12        => WP12(1,1,1,1)
    readyc1_param%d_wp13        => WP13(1,1,1,1)
    readyc1_param%d_at0         => AT0(1,1,1)
    readyc1_param%d_atn         => ATN(1,1,1)
    readyc1_param%d_ate         => ATE(1,1,1)
    readyc1_param%d_atne        => ATNE(1,1,1)

    !readyc1_2_param
    readyc1_2_param%cuBlock(1)     = imt
    readyc1_2_param%cuBlock(2)     = 1
    readyc1_2_param%cuBlock(3)     = 1
    readyc1_2_param%cuGrid(1)      = 1
    readyc1_2_param%cuGrid(2)      = jmt
    readyc1_2_param%cuGrid(3)      = kmm1
    readyc1_2_param%d_wp12        => WP12(1,1,1,1)
    readyc1_2_param%d_wp13        => WP13(1,1,1,1)
    readyc1_2_param%d_vit         => VIT(1,1,1,1)
    readyc1_2_param%d_s2t         => s2t(1,1,1,1)
    readyc1_2_param%d_riv1        => riv1(1,1,1,1)
    readyc1_2_param%d_riv2        => riv2(1,1,1,1)
    readyc1_2_param%d_odzt        => ODZT(1)
    readyc1_2_param%d_rit         => rit(1,1,1,1)
    readyc1_2_param%d_rict        => rict(1,1,1,1)
    readyc1_2_param%d_up          => UP(1,1,1,1)
    readyc1_2_param%d_vp          => VP(1,1,1,1)
    readyc1_2_param%d_s2u         => s2u(1,1,1,1)
    readyc1_2_param%d_viv         => viv(1,1,1,1)
    readyc1_2_param%d_riu         => riu(1,1,0,1)
    readyc1_2_param%d_ric         => ric(1,1,1,1)
    !readyc1_2_param%rank          = my_task

#ifdef BCKMEX
!readyc1_3_param
    readyc1_3_param%cuBlock(1)     = imt
    readyc1_3_param%cuBlock(2)     = 1
    readyc1_3_param%cuBlock(3)     = 1
    readyc1_3_param%cuGrid(1)      = 1
    readyc1_3_param%cuGrid(2)      = jmt
    readyc1_3_param%cuGrid(3)      = 1
    readyc1_3_param%d_diff_back    => diff_back(1,1,1)
    readyc1_3_param%d_diff_back_sh => diff_back_sh(1,1,1)
    readyc1_3_param%d_diff_back_nh => diff_back_nh(1,1,1)
    readyc1_3_param%d_tlat         => tlat(1,1,1)
    readyc1_3_param%d_diff_back_coef_max        = diff_back_coef_max
    readyc1_3_param%d_diff_back_eq              = diff_back_eq
    readyc1_3_param%d_diff_back_coef            = diff_back_coef
#endif


    !readyc2_parm
    readyc2_parm%cuBlock(1)     = imt
    readyc2_parm%cuBlock(2)     = 1
    readyc2_parm%cuBlock(3)     = 1
    readyc2_parm%cuGrid(1)      = 1
    readyc2_parm%cuGrid(2)      = jmt
    readyc2_parm%cuGrid(3)      = kmm1
    readyc2_parm%d_akmu            =>akmu(1,1,1,1)
    readyc2_parm%d_au0             =>AU0(1,1,1)
    readyc2_parm%d_akmt            =>akmt(1,1,1,1)
    readyc2_parm%d_aus             =>aus(1,1,1)
    readyc2_parm%d_auw             =>AUW(1,1,1)
    readyc2_parm%d_ausw            =>AUSW(1,1,1)
    readyc2_parm%d_viv             =>viv(1,1,1,1)

    !readyc_upwell_0
    readyc_upwell_0_param%cuBlock(1)     = imt
    readyc_upwell_0_param%cuBlock(2)     = 1
    readyc_upwell_0_param%cuBlock(3)     = 1
    readyc_upwell_0_param%cuGrid(1)      = 1
    readyc_upwell_0_param%cuGrid(2)      = jmt
    readyc_upwell_0_param%cuGrid(3)      = km
    readyc_upwell_0_param%rank           = my_task
    readyc_upwell_0_param%d_au0            =>AU0(1,1,1)
    readyc_upwell_0_param%d_aus            =>aus(1,1,1)
    readyc_upwell_0_param%d_auw            =>AUW(1,1,1)
    readyc_upwell_0_param%d_ausw           =>AUSW(1,1,1)
    readyc_upwell_0_param%d_h0             =>H0(1,1,1)
    readyc_upwell_0_param%d_ohbu           =>ohbu(1,1,1)
    readyc_upwell_0_param%d_u              =>u(1,1,1,1)
    readyc_upwell_0_param%d_uk             =>UK(1,1,1,1)
    readyc_upwell_0_param%d_v              =>v(1,1,1,1)
    readyc_upwell_0_param%d_vk             =>VK(1,1,1,1)
    readyc_upwell_0_param%d_work           =>WORK(1,1,1)

    !readyc_upwell
    readyc_upwell_param%cuBlock(1)     = imt
    readyc_upwell_param%cuBlock(2)     = 1
    readyc_upwell_param%cuBlock(3)     = 1
    readyc_upwell_param%cuGrid(1)      = 1
    readyc_upwell_param%cuGrid(2)      = jmt
    readyc_upwell_param%cuGrid(3)      = km
    readyc_upwell_param%rank           = my_task
    readyc_upwell_param%d_dzp            =>DZP(1)
    readyc_upwell_param%d_h0             =>H0(1,1,1)
    readyc_upwell_param%d_hts            =>hts(1,1,1)
    readyc_upwell_param%d_htw            =>HTW(1,1,1)
    readyc_upwell_param%d_ohbt           =>OHBT(1,1,1)
    readyc_upwell_param%d_tarea_r        =>TAREA_R(1,1,1)
    readyc_upwell_param%d_uk             =>uk(1,1,1,1)
    readyc_upwell_param%d_vit            =>VIT(1,1,1,1)
    readyc_upwell_param%d_vk             =>VK(1,1,1,1)
    readyc_upwell_param%d_wka            =>wka(1,1,1,1)
    readyc_upwell_param%d_work           =>WORK(1,1,1)
    readyc_upwell_param%d_ws             =>ws(1,1,1,1)
    readyc_upwell_param%d_kmt            =>KMT(1,1,1)

    readyc_upwell_1_param%cuBlock(1)     = imt
    readyc_upwell_1_param%cuBlock(2)     = 1
    readyc_upwell_1_param%cuBlock(3)     = 1
    readyc_upwell_1_param%cuGrid(1)      = (imt + readyc_upwell_1_param%cuBlock(1) - 1) / readyc_upwell_1_param%cuBlock(1)
    readyc_upwell_1_param%cuGrid(2)      = (jmt + readyc_upwell_1_param%cuBlock(2) - 1) / readyc_upwell_1_param%cuBlock(2)
    readyc_upwell_1_param%cuGrid(3)      = 1
    readyc_upwell_1_param%rank           = my_task
    readyc_upwell_1_param%d_vit            =>VIT(1,1,1,1)
    readyc_upwell_1_param%d_dzp            =>DZP(1)
    readyc_upwell_1_param%d_wka            =>wka(1,1,1,1)
    readyc_upwell_1_param%d_work           =>WORK(1,1,1)

    readyc_upwell_2_param%cuBlock(1)     = imt
    readyc_upwell_2_param%cuBlock(2)     = 1
    readyc_upwell_2_param%cuBlock(3)     = 1
    readyc_upwell_2_param%cuGrid(1)      = (imt + readyc_upwell_2_param%cuBlock(1) - 1) / readyc_upwell_2_param%cuBlock(1)
    readyc_upwell_2_param%cuGrid(2)      = (jmt + readyc_upwell_2_param%cuBlock(2) - 1) / readyc_upwell_2_param%cuBlock(2)
    readyc_upwell_2_param%cuGrid(3)      = 1
    readyc_upwell_2_param%rank           = my_task
    readyc_upwell_2_param%d_ws             =>ws(1,1,1,1)
    readyc_upwell_2_param%d_vit            =>VIT(1,1,1,1)
    readyc_upwell_2_param%d_h0             =>H0(1,1,1)
    readyc_upwell_2_param%d_dzp            =>DZP(1)
    readyc_upwell_2_param%d_ohbt           =>OHBT(1,1,1)
    readyc_upwell_2_param%d_wka            =>wka(1,1,1,1)
    readyc_upwell_2_param%d_work           =>WORK(1,1,1)
    readyc_upwell_2_param%d_wsbak          =>wsbak(1,1,1,1)

    !readyc3
    readyc3_param%cuBlock(1)     = imt
    readyc3_param%cuBlock(2)     = 1
    readyc3_param%cuBlock(3)     = 1
    readyc3_param%cuGrid(1)      = 1
    readyc3_param%cuGrid(2)      = jmt
    readyc3_param%cuGrid(3)      = km
    readyc3_param%d_wka         =>wka(1,1,1,1)
    readyc3_param%d_ws          =>ws(1,1,1,1)
    readyc3_param%d_au0         =>AU0(1,1,1)
    readyc3_param%d_ausw        =>AUSW(1,1,1)
    readyc3_param%d_aus         =>aus(1,1,1)
    readyc3_param%d_auw         =>AUW(1,1,1)
    readyc3_param%d_u_wface     =>u_wface(1,1,1)
    readyc3_param%d_v_sface     =>v_sface(1,1,1)
    readyc3_param%d_u           =>u(1,1,1,1)
    readyc3_param%d_v           =>v(1,1,1,1)
    readyc3_param%d_hue         =>HUE(1,1,1)
    readyc3_param%d_hun         =>HUN(1,1,1)

    !advection_momentum
    advection_momentum_param%cuBlock(1)     = imt
    advection_momentum_param%cuBlock(2)     = 1
    advection_momentum_param%cuBlock(3)     = 1
    advection_momentum_param%cuGrid(1)      = 1
    advection_momentum_param%cuGrid(2)      = jmt
    advection_momentum_param%cuGrid(3)      = km
    advection_momentum_param%d_dlu         =>dlu(1,1,1,1)
    advection_momentum_param%d_dlv         =>dlv(1,1,1,1)
    advection_momentum_param%d_odzp        =>odzp(1)
    advection_momentum_param%d_u           =>u(1,1,1,1)
    advection_momentum_param%d_uarea_r     =>UAREA_R(1,1,1)
    advection_momentum_param%d_u_wface     =>u_wface(1,1,1)
    advection_momentum_param%d_v           =>v(1,1,1,1)
    advection_momentum_param%d_v_sface     =>v_sface(1,1,1)
    advection_momentum_param%d_wka         =>wka(1,1,1,1)
    !    advection_momentum_param%d_utest       =>v(1,1,1,1)

    readyc_vinteg_param%cuBlock(1)     = imt
    readyc_vinteg_param%cuBlock(2)     = 1
    readyc_vinteg_param%cuBlock(3)     = 1
    readyc_vinteg_param%cuGrid(1)      = 1
    readyc_vinteg_param%cuGrid(2)      = jmt
    readyc_vinteg_param%cuGrid(3)      = km
    readyc_vinteg_param%d_dlu         =>dlu(1,1,1,1)
    readyc_vinteg_param%d_dlv         =>dlv(1,1,1,1)
    readyc_vinteg_param%d_dlub      =>dlub(1,1,1)  !to allocate
    readyc_vinteg_param%d_dlvb      =>dlvb(1,1,1)  ! to allocate
    readyc_vinteg_param%d_dzp       =>DZP(1)
    readyc_vinteg_param%d_ohbu      =>ohbu(1,1,1)
    readyc_vinteg_param%d_viv       =>VIV(1,1,1,1)

    !readyc4
    readyc4_param%cuBlock(1)     = imt
    readyc4_param%cuBlock(2)     = 1
    readyc4_param%cuBlock(3)     = 1
    readyc4_param%cuGrid(1)      = 1
    readyc4_param%cuGrid(2)      = jmt
    readyc4_param%cuGrid(3)      = 1
    readyc4_param%d_kmu       =>KMU(1,1,1)
    readyc4_param%d_sbcx      =>SBCX(1,1,1)
    readyc4_param%d_sbcy      =>SBCY(1,1,1)
    readyc4_param%d_su        =>su(1,1,1)
    readyc4_param%d_sv        =>sv(1,1,1)
    readyc4_param%d_bbcx      =>BBCX(1,1,1)
    readyc4_param%d_od0       =od0
    readyc4_param%d_c0f       =C0F   !to check
    readyc4_param%d_up        =>UP(1,1,1,1)
    readyc4_param%d_vp        =>VP(1,1,1,1)
    readyc4_param%d_snlat     =>SNLAT(1,1,1)
    readyc4_param%d_cag       =CAG
    readyc4_param%d_sag       =SAG
    readyc4_param%d_dlub      =>dlub(1,1,1)  !to allocate
    readyc4_param%d_dlvb      =>dlvb(1,1,1)  ! to allocate
    readyc4_param%d_dzp       =>DZP(1)
    readyc4_param%d_ohbu      =>ohbu(1,1,1)
    readyc4_param%d_dlu       =>dlu(1,1,1,1)
    readyc4_param%d_viv       =>VIV(1,1,1,1)
    readyc4_param%d_dlv       =>dlv(1,1,1,1)
    readyc4_param%d_bbcy      =>BBCY(1,1,1)

    !readyc4_2
    readyc4_2_param%cuBlock(1)     = imt
    readyc4_2_param%cuBlock(2)     = 1
    readyc4_2_param%cuBlock(3)     = 1
    readyc4_2_param%cuGrid(1)      = 1
    readyc4_2_param%cuGrid(2)      = jmt
    readyc4_2_param%cuGrid(3)      = km
    readyc4_2_param%d_od0       =od0
    readyc4_2_param%d_sag       =SAG
    readyc4_2_param%d_cag       =CAG
    readyc4_2_param%d_su        =>su(1,1,1)
    readyc4_2_param%d_sv        =>sv(1,1,1)
    readyc4_2_param%d_up        =>UP(1,1,1,1)
    readyc4_2_param%d_vp        =>VP(1,1,1,1)
    readyc4_2_param%d_akmu      =>akmu(1,1,1,1)
    readyc4_2_param%d_odzt      =>ODZT(1)
    readyc4_2_param%d_viv       =>VIV(1,1,1,1)
    readyc4_2_param%d_wka       =>wka(1,1,1,1)
    readyc4_2_param%d_snlat     =>SNLAT(1,1,1)
    readyc4_2_param%d_dlv       =>dlv(1,1,1,1)
    readyc4_2_param%d_odzp      =>odzp(1)
    readyc4_2_param%d_dlu       =>dlu(1,1,1,1)

    !hdel_div
    hdel4_div_param%cuBlock(1)     = imt
    hdel4_div_param%cuBlock(2)     = 1
    hdel4_div_param%cuBlock(3)     = 1
    hdel4_div_param%cuGrid(1)      = 1
    hdel4_div_param%cuGrid(2)      = jmt
    hdel4_div_param%cuGrid(3)      = km
    hdel4_div_param%d_kmt          =>KMT(1,1,1)
    hdel4_div_param%d_hts          =>hts(1,1,1)
    hdel4_div_param%d_htw          =>HTW(1,1,1)
    hdel4_div_param%d_tarea_r      =>TAREA_R(1,1,1)
    hdel4_div_param%d_up           =>UP(1,1,1,1)
    hdel4_div_param%d_vp           =>VP(1,1,1,1)
    hdel4_div_param%d_wka          =>wka(1,1,1,1)

    hdel4_grad_param%cuBlock(1)     = imt
    hdel4_grad_param%cuBlock(2)     = 1
    hdel4_grad_param%cuBlock(3)     = 1
    hdel4_grad_param%cuGrid(1)      = 1
    hdel4_grad_param%cuGrid(2)      = jmt
    hdel4_grad_param%cuGrid(3)      = km
    hdel4_grad_param%rank           = my_task
    hdel4_grad_param%d_kmu          => KMU(1,1,1)
    hdel4_grad_param%d_wka          => WKA(1,1,1,1)
    hdel4_grad_param%d_dxur         => DXUR(1,1,1)
    hdel4_grad_param%d_dyur         => DYUR(1,1,1)
    hdel4_grad_param%d_am_factor    => am_factor_m(1,1,1)
    hdel4_grad_param%d_wkz          =>wkz(1,1,1,1)

    hdel4_zcurl_param%cuBlock(1)     = imt
    hdel4_zcurl_param%cuBlock(2)     = 1
    hdel4_zcurl_param%cuBlock(3)     = 1
    hdel4_zcurl_param%cuGrid(1)      = 1
    hdel4_zcurl_param%cuGrid(2)      = jmt
    hdel4_zcurl_param%cuGrid(3)      = km
    hdel4_zcurl_param%rank           = my_task
    hdel4_zcurl_param%d_kmt          => KMT(1,1,1)
    hdel4_zcurl_param%d_dxu          => DXU(1,1,1)
    hdel4_zcurl_param%d_dyu          => DYU(1,1,1)
    hdel4_zcurl_param%d_up           => UP(1,1,1,1)
    hdel4_zcurl_param%d_vp           => VP(1,1,1,1)
    hdel4_zcurl_param%d_wkz          => WKZ(1,1,1,1)
    hdel4_zcurl_param%d_tarea_r      =>TAREA_R(1,1,1)

    hdiffu_del4_0_param%cuBlock(1)     = imt
    hdiffu_del4_0_param%cuBlock(2)     = 1
    hdiffu_del4_0_param%cuBlock(3)     = 1
    hdiffu_del4_0_param%cuGrid(1)      = 1
    hdiffu_del4_0_param%cuGrid(2)      = jmt
    hdiffu_del4_0_param%cuGrid(3)      = km
    hdiffu_del4_0_param%rank           = my_task
    hdiffu_del4_0_param%d_kmu          => KMU(1,1,1)
    hdiffu_del4_0_param%d_am           = AM_BIHAR
    hdiffu_del4_0_param%d_uarea_r      => UAREA_R(1,1,1)
    hdiffu_del4_0_param%d_up           => UP(1,1,1,1)
    hdiffu_del4_0_param%d_vp           => VP(1,1,1,1)
    hdiffu_del4_0_param%d_duc          => DUC(1,1,1)
    hdiffu_del4_0_param%d_dum          => DUM(1,1,1)
    hdiffu_del4_0_param%d_dun          => DUN(1,1,1)
    hdiffu_del4_0_param%d_dus          => DUS(1,1,1)
    hdiffu_del4_0_param%d_due          => DUE(1,1,1)
    hdiffu_del4_0_param%d_duw          => DUW(1,1,1)
    hdiffu_del4_0_param%d_dmc          => DMC(1,1,1)
    hdiffu_del4_0_param%d_dmn          => DMN(1,1,1)
    hdiffu_del4_0_param%d_dms          => DMS(1,1,1)
    hdiffu_del4_0_param%d_dme          => DME(1,1,1)
    hdiffu_del4_0_param%d_dmw          => DMW(1,1,1)
    hdiffu_del4_0_param%d_amf          => AMF(1,1,1)
    hdiffu_del4_0_param%d_am_factor    => am_factor_m(1,1,1)
    hdiffu_del4_0_param%d_d2uk         => d2ukm(1,1,1)
    hdiffu_del4_0_param%d_d2vk         => d2vkm(1,1,1)

    hdiffu_del4_1_param%cuBlock(1)     = imt
    hdiffu_del4_1_param%cuBlock(2)     = 1
    hdiffu_del4_1_param%cuBlock(3)     = 1
    hdiffu_del4_1_param%cuGrid(1)      = 1
    hdiffu_del4_1_param%cuGrid(2)      = jmt
    hdiffu_del4_1_param%cuGrid(3)      = km
    hdiffu_del4_1_param%rank           = my_task
    hdiffu_del4_1_param%d_kmu          => KMU(1,1,1)
    hdiffu_del4_1_param%d_am           = AM_BIHAR
    hdiffu_del4_1_param%d_c0f          = C0F
    hdiffu_del4_1_param%d_wka          => WKA(1,1,1,1)
    hdiffu_del4_1_param%d_up          => UP(1,1,1,1)
    hdiffu_del4_1_param%d_vp          => VP(1,1,1,1)
    hdiffu_del4_1_param%d_duc          => DUC(1,1,1)
    hdiffu_del4_1_param%d_dum          => DUM(1,1,1)
    hdiffu_del4_1_param%d_dun          => DUN(1,1,1)
    hdiffu_del4_1_param%d_dus          => DUS(1,1,1)
    hdiffu_del4_1_param%d_due          => DUE(1,1,1)
    hdiffu_del4_1_param%d_duw          => DUW(1,1,1)
    hdiffu_del4_1_param%d_dmc          => DMC(1,1,1)
    hdiffu_del4_1_param%d_dmn          => DMN(1,1,1)
    hdiffu_del4_1_param%d_dms          => DMS(1,1,1)
    hdiffu_del4_1_param%d_dme          => DME(1,1,1)
    hdiffu_del4_1_param%d_dmw          => DMW(1,1,1)
    hdiffu_del4_1_param%d_dlu          => dlu(1,1,1,1)
    hdiffu_del4_1_param%d_dlv          => dlv(1,1,1,1)
    hdiffu_del4_1_param%d_d2uk         => d2ukm(1,1,1)
    hdiffu_del4_1_param%d_d2vk         => d2vkm(1,1,1)
    !hdiffu_del4_1_param%d_hdukkm         => hdukkm(1,1,1)
    !hdiffu_del4_1_param%d_hdvkkm         => hdvkkm(1,1,1)
    !hdiffu_del4_1_param%d_am_f_km         => am_f_km(1,1,1)
    hdiffu_del4_1_param%d_dts          = DTS
    hdiffu_del4_1_param%d_value        = 1.7*(0.25*3.6D7/real(IMT_GLOBAL, kind=8))**2.5
    hdiffu_del4_1_param%d_value2       = (real(IMT_GLOBAL, kind=8)/1.8D4)**1.3
    hdiffu_del4_1_param%d_am_factor    => am_factor_m(1,1,1)
    hdiffu_del4_1_param%d_uarea_r      => uarea(1,1,1)
    hdiffu_del4_1_param%d_fac_vis_uv   => factor_vis_uv(1,1,1)


    hdiffu_del2_param%cuBlock(1)     = imt
    hdiffu_del2_param%cuBlock(2)     = 1
    hdiffu_del2_param%cuBlock(3)     = 1
    hdiffu_del2_param%cuGrid(1)      = (imt + hdiffu_del2_param%cuBlock(1) - 1) / hdiffu_del2_param%cuBlock(1)
    hdiffu_del2_param%cuGrid(2)      = (jmt + hdiffu_del2_param%cuBlock(2) - 1) / hdiffu_del2_param%cuBlock(2)
    hdiffu_del2_param%cuGrid(3)      = (km +  hdiffu_del2_param%cuBlock(3) - 1) / hdiffu_del2_param%cuBlock(3)
    hdiffu_del2_param%rank           = my_task
    hdiffu_del2_param%d_kmu          => KMU(1,1,1)
    hdiffu_del2_param%d_am           = AM_HOR
    hdiffu_del2_param%d_c0f          = C0F
    hdiffu_del2_param%d_wka          => WKA(1,1,1,1)
    hdiffu_del2_param%d_up          => UP(1,1,1,1)
    hdiffu_del2_param%d_vp          => VP(1,1,1,1)
    hdiffu_del2_param%d_duc          => DUC(1,1,1)
    hdiffu_del2_param%d_dum          => DUM(1,1,1)
    hdiffu_del2_param%d_dun          => DUN(1,1,1)
    hdiffu_del2_param%d_dus          => DUS(1,1,1)
    hdiffu_del2_param%d_due          => DUE(1,1,1)
    hdiffu_del2_param%d_duw          => DUW(1,1,1)
    hdiffu_del2_param%d_dmc          => DMC(1,1,1)
    hdiffu_del2_param%d_dmn          => DMN(1,1,1)
    hdiffu_del2_param%d_dms          => DMS(1,1,1)
    hdiffu_del2_param%d_dme          => DME(1,1,1)
    hdiffu_del2_param%d_dmw          => DMW(1,1,1)
    hdiffu_del2_param%d_dlu          => dlu(1,1,1,1)
    hdiffu_del2_param%d_dlv          => dlv(1,1,1,1)
    hdiffu_del2_param%d_viv          => VIV(1,1,1,1)
#ifdef CANUTO
    turb_param = &
    (/rri,rnd2on2,dri,deltheta_r,b1,theta_rcrp,theta_rcrn,real(ncc),dfricmx,dwndmix,dand2on2/)
   if(my_task .eq. 0) then
       write(*,*) "turb_param:",turb_param(1),turb_param(2),turb_param(3),turb_param(4),turb_param(5),turb_param(6), &
            turb_param(7),turb_param(8),turb_param(9),turb_param(10),turb_param(11)
   end if
!readyc_turb_param
    readyc_turbo_parm%cuBlock(1)     = imt
    readyc_turbo_parm%cuBlock(2)     = 1
    readyc_turbo_parm%cuBlock(3)     = 1
    readyc_turbo_parm%cuGrid(1)      = 1
    readyc_turbo_parm%cuGrid(2)      = jmt
    readyc_turbo_parm%cuGrid(3)      = 1
    readyc_turbo_parm%d_km              = km
    readyc_turbo_parm%d_imt             = imt
    readyc_turbo_parm%d_jmt             = jmt
    readyc_turbo_parm%rank              = my_task
    readyc_turbo_parm%d_akmt            =>akmt(1,1,1,1)
    readyc_turbo_parm%d_akt             =>AKT(1,1,1,1,1)
    readyc_turbo_parm%d_ak_tide         =>AK_TIDE(1,1,1,1)
    readyc_turbo_parm%d_amld            =>amld(1,1,1)
    readyc_turbo_parm%d_amtaun2a1       =>amtaun2a1(-mt)
    readyc_turbo_parm%d_and2on2a1       =>and2on2a1(-mt)
    readyc_turbo_parm%d_at              =>AT(1,1,1,1,1)
    readyc_turbo_parm%d_back_ra_r       =>back_ra_r(-n_theta_r_oct)
    readyc_turbo_parm%d_buoysol         =>BUOYSOL(1,1,1)
    readyc_turbo_parm%d_buoytur         =>BUOYTUR(1,1,1)
    readyc_turbo_parm%d_dzp             =>DZP(1)
    readyc_turbo_parm%d_fcort           =>FCORT(1,1,1)
    readyc_turbo_parm%d_fztidal         =>fztidal(1,1,1,1)
    readyc_turbo_parm%d_fz_tide         =>FZ_TIDE(1,1,1,1)
    readyc_turbo_parm%d_pdensity        =>pdensity(1,1,1,1)
    readyc_turbo_parm%d_rib             =>rib(-mt)
    readyc_turbo_parm%d_ridb            =>ridb(-mt)
    readyc_turbo_parm%d_ricdttms        =>ricdtTmS(1,1,1,1)
    readyc_turbo_parm%d_richardson      =>richardson(1,1,1,1)
    readyc_turbo_parm%d_rict            =>rict(1,1,1,1)
    readyc_turbo_parm%d_rit             =>rit(1,1,1,1)
    readyc_turbo_parm%d_s2t             =>s2t(1,1,1,1)
    readyc_turbo_parm%d_sha1            =>sha1(-mt)
    readyc_turbo_parm%d_shb             =>shb(-mt,-mt)
    readyc_turbo_parm%d_sh_r1           =>sh_r1(-n_theta_r_oct)
    readyc_turbo_parm%d_slq2b           =>slq2b(-mt,-mt)
    readyc_turbo_parm%d_slq2_r1         =>slq2_r1(-n_theta_r_oct)
    readyc_turbo_parm%d_sma1            =>sma1(-mt)
    readyc_turbo_parm%d_smb             =>smb(-mt,-mt)
    readyc_turbo_parm%d_sm_r1           =>sm_r1(-n_theta_r_oct)
    readyc_turbo_parm%d_ssa1            =>ssa1(-mt)
    readyc_turbo_parm%d_ssb             =>ssb(-mt,-mt)
    readyc_turbo_parm%d_ss_r1           =>ss_r1(-n_theta_r_oct)
    readyc_turbo_parm%d_turb_param      =>turb_param(1)
    readyc_turbo_parm%d_vit             =>VIT(1,1,1,1)
    readyc_turbo_parm%d_wave_dis        =>WAVE_DIS(1,1,1)
    readyc_turbo_parm%d_wp3_tidal       =>wp3_tidal(1,1,1,1)
    readyc_turbo_parm%d_zkp             =>ZKP(1)
    readyc_turbo_parm%d_irimax          =>irimax(-mt)
    readyc_turbo_parm%d_kmt             =>KMT(1,1,1)
#ifdef BCKMEX
    readyc_turbo_parm%d_diff_back        =>diff_back(1,1,1)
#endif

#endif

    canuto2010_parm%cuBlock(1)     = imt
    canuto2010_parm%cuBlock(2)     = 1
    canuto2010_parm%cuBlock(3)     = 1
    canuto2010_parm%cuGrid(1)      = 1
    canuto2010_parm%cuGrid(2)      = jmt
    canuto2010_parm%cuGrid(3)      = 1
    canuto2010_parm%d_kmt          => KMT(1,1,1)
    canuto2010_parm%d_dzp          => DZP(1)
    canuto2010_parm%d_zkp          => ZKP(1)
    canuto2010_parm%d_ulat         => ulat(1,1,1)
    canuto2010_parm%d_rit          => rit(1,1,1,1)
    canuto2010_parm%d_rict         => rict(1,1,1,1)
    canuto2010_parm%d_ricdt        => ricdt(1,1,1,1)
    canuto2010_parm%d_s2t          => s2t(1,1,1,1)
    canuto2010_parm%d_vit          => VIT(1,1,1,1)
    canuto2010_parm%d_akt          => AKT(1,1,1,1,1)
    canuto2010_parm%d_at           => AT(1,1,1,1,1)
    canuto2010_parm%d_to           => TO(1)
    canuto2010_parm%d_so           => SO(1)
    canuto2010_parm%d_c            => C(1,1)
    canuto2010_parm%d_wave_dis     => WAVE_DIS(1,1,1)
    canuto2010_parm%d_fz_tide      => FZ_TIDE(1,1,1,1)
    canuto2010_parm%d_akmt         => akmt(1,1,1,1)
    !canuto2010_parm%d_wk1         => wk1tmp(1,1,1,1)
    !canuto2010_parm%d_wk2         => wk2tmp(1,1,1,1)
    !canuto2010_parm%d_wk3         => wk3tmp(1,1,1,1)
    !canuto2010_parm%d_wp3         => wp3tmp(1,1,1,1)
    canuto2010_parm%d_mldlev      => mldlev(1,1,1)
    canuto2010_parm%d_amld        => amld(1,1,1)
    canuto2010_parm%rank          = my_task

    canuto2010_parm%d_pi0_1        =  1.0d0/(sqrt(27.0d0/5.0d0*1.66d0**3)*(1+1.0d0/0.72d0))
 !   canuto2010_parm%d_pi0_2        =  1.0e0/3.0e0


!    num_ocean_point = 1
!    DO IBLOCK = 1, NBLOCKS_CLINIC
!            DO J = 1, JMT
!                DO I = 1, IMT
!                   if(vit(i,j,1,iblock) .gt. 0.5) then
!                       biti = and(i,1023) !0x3ff
!                       bitj = and(j,1023) !0x3ff
!                       bitjj = shiftl(bitj,10)
!                       bit_3d_ocean_point(num_ocean_point) = ior(bitjj, biti)
!!                       if(my_task .eq. 71) then
!!                           write(*,*) "i,j,bitjj,ocean_point",i,j,bitjj,num_ocean_point,bit_3d_ocean_point(num_ocean_point)
!!                       end if
!                       num_ocean_point = num_ocean_point + 1
!                   end if
!                END DO
!            END DO
!    END DO
!    canuto2010_parm%num_ocean_point        =  num_ocean_point
!    canuto2010_parm%bit_3d_ocean_point     =>  bit_3d_ocean_point(1)
end subroutine readyc_sw_init

subroutine barotr_sw_init

        use pconst_mod
        use tracer_mod
        use dyn_mod
        use work_mod
        use pmix_mod
        use constant_mod
        !    use forc_mod, only: psa,USTAR,BUOYTUR, BUOYSOL,NSWV,SWV
        use forc_mod
        use grid
        use blocks
        use domain
        implicit none
        type(block):: this_block
        integer :: iblock

        barotr1_param%cuBlock(1)     = imt
        barotr1_param%cuBlock(2)     = 1
        barotr1_param%cuBlock(3)     = 1
        barotr1_param%cuGrid(1)      = 1
        barotr1_param%cuGrid(2)      = jmt
        barotr1_param%cuGrid(3)      = 1
        barotr1_param%d_km           = km
        barotr1_param%d_imt          = imt
        barotr1_param%d_jmt          = jmt
        barotr1_param%rank           = my_task
        barotr1_param%d_work         => WORK(1,1,1)
        barotr1_param%d_wka          => WKA(1,1,1,1)
        barotr1_param%d_dzph         => DZPH(1,1,1)
        barotr1_param%d_h0           => H0(1,1,1)
        barotr1_param%d_ub           => UB(1,1,1)
        barotr1_param%d_vb           => VB(1,1,1)
        barotr1_param%d_au0          => AU0(1,1,1)
        barotr1_param%d_ausw         => AUSW(1,1,1)
        barotr1_param%d_aus          => AUS(1,1,1)
        barotr1_param%d_auw          => AUW(1,1,1)


        barotr1_1_param%cuBlock(1)     = imt
        barotr1_1_param%cuBlock(2)     = 1
        barotr1_1_param%cuBlock(3)     = 1
        barotr1_1_param%cuGrid(1)      = 1
        barotr1_1_param%cuGrid(2)      = jmt
        barotr1_1_param%cuGrid(3)      = 1
        barotr1_1_param%d_imt          = imt
        barotr1_1_param%d_jmt          = jmt
        barotr1_1_param%rank           = my_task
        barotr1_1_param%d_kmt          => KMT(1,1,1)
        barotr1_1_param%d_work         => WORK(1,1,1)
        barotr1_1_param%d_wka          => WKA(1,1,1,1)
        barotr1_1_param%d_vit          => VIT(1,1,1,1)
        barotr1_1_param%d_htw          => HTW(1,1,1)
        barotr1_1_param%d_hts          => HTS(1,1,1)
        barotr1_1_param%d_tarea_r      => TAREA_R(1,1,1)


        barotr2_1_param%cuBlock(1)     = imt
        barotr2_1_param%cuBlock(2)     = 1
        barotr2_1_param%cuBlock(3)     = 1
        barotr2_1_param%cuGrid(1)      = 1
        barotr2_1_param%cuGrid(2)      = jmt
        barotr2_1_param%cuGrid(3)      = 1
        barotr2_1_param%d_imt          = imt
        barotr2_1_param%d_jmt          = jmt
        barotr2_1_param%d_dtb          = DTB
        barotr2_1_param%d_h0           => H0(1,1,1)
        barotr2_1_param%d_h0p          => H0P(1,1,1)
        barotr2_1_param%d_work         => WORK(1,1,1)
#ifdef TIDEFORCE
        barotr2_1_param%BETA_TIDE      = BETA_TIDE
        barotr2_1_param%d_h0_tide      => h0_tide(1,1,1)
        barotr2_1_param%d_elev         => elevat(1,1,1)
#endif

#ifdef BIHAR

        bhdel4_div_param%cuBlock(1)     = imt
        bhdel4_div_param%cuBlock(2)     = 1
        bhdel4_div_param%cuBlock(3)     = 1
        bhdel4_div_param%cuGrid(1)      = 1
        bhdel4_div_param%cuGrid(2)      = jmt
        bhdel4_div_param%cuGrid(3)      = 1
        bhdel4_div_param%rank           = my_task
        bhdel4_div_param%d_kmt          => KMT(1,1,1)
        bhdel4_div_param%d_hts          => HTS(1,1,1)
        bhdel4_div_param%d_htw          => HTW(1,1,1)
        bhdel4_div_param%d_tarea_r      => TAREA_R(1,1,1)
        bhdel4_div_param%d_ubp          => UBP(1,1,1)
        bhdel4_div_param%d_vbp          => VBP(1,1,1)
        bhdel4_div_param%d_wka          => WKA(1,1,1,1)
        bhdel4_div_param%d_diverg       => diverg(1,1,1)

    bhdel4_grad_param%cuBlock(1)     = imt
    bhdel4_grad_param%cuBlock(2)     = 1
    bhdel4_grad_param%cuBlock(3)     = 1
    bhdel4_grad_param%cuGrid(1)      = 1
    bhdel4_grad_param%cuGrid(2)      = jmt
    bhdel4_grad_param%cuGrid(3)      = 1
    bhdel4_grad_param%rank           = my_task
    bhdel4_grad_param%d_kmu          => KMU(1,1,1)
    bhdel4_grad_param%d_wka          => WKA(1,1,1,1)
    bhdel4_grad_param%d_dxur         => DXUR(1,1,1)
    bhdel4_grad_param%d_dyur         => DYUR(1,1,1)
    bhdel4_grad_param%d_am_factor    => am_factor_d(1,1,1)
    bhdel4_grad_param%d_z_curl       => z_curl(1,1,1)
    bhdel4_grad_param%d_diverg       => diverg(1,1,1)

    bhdel4_zcurl_param%cuBlock(1)     = imt
    bhdel4_zcurl_param%cuBlock(2)     = 1
    bhdel4_zcurl_param%cuBlock(3)     = 1
    bhdel4_zcurl_param%cuGrid(1)      = 1
    bhdel4_zcurl_param%cuGrid(2)      = jmt
    bhdel4_zcurl_param%cuGrid(3)      = 1
    bhdel4_zcurl_param%rank           = my_task
    bhdel4_zcurl_param%d_kmt          => KMT(1,1,1)
    bhdel4_zcurl_param%d_dxu          => DXU(1,1,1)
    bhdel4_zcurl_param%d_dyu          => DYU(1,1,1)
    bhdel4_zcurl_param%d_ubp          => UBP(1,1,1)
    bhdel4_zcurl_param%d_vbp          => VBP(1,1,1)
    bhdel4_zcurl_param%d_wka          => WKA(1,1,1,1)
    bhdel4_zcurl_param%d_tarea_r      => TAREA_R(1,1,1)
    bhdel4_zcurl_param%d_z_curl       => z_curl(1,1,1)


    DO IBLOCK = 1, NBLOCKS_CLINIC
    this_block = get_block(blocks_clinic(iblock),iblock)

    barotr2_udel4_0_param%cuBlock(1)     = imt
    barotr2_udel4_0_param%cuBlock(2)     = 1
    barotr2_udel4_0_param%cuBlock(3)     = 1
    barotr2_udel4_0_param%cuGrid(1)      = 1
    barotr2_udel4_0_param%cuGrid(2)      = jmt
    barotr2_udel4_0_param%cuGrid(3)      = 1
    barotr2_udel4_0_param%d_jb           = this_block%JB
    barotr2_udel4_0_param%d_je           = this_block%JE
    barotr2_udel4_0_param%d_ib           = this_block%IB
    barotr2_udel4_0_param%d_ie           = this_block%IE
    barotr2_udel4_0_param%d_imt          = imt
    barotr2_udel4_0_param%d_jmt          = jmt
    barotr2_udel4_0_param%rank           = my_task
    barotr2_udel4_0_param%d_kmu          => KMU(1,1,1)
    barotr2_udel4_0_param%d_am           = AM_BIHAR
    barotr2_udel4_0_param%d_ubp          => UBP(1,1,1)
    barotr2_udel4_0_param%d_vbp          => VBP(1,1,1)
    barotr2_udel4_0_param%d_uarea_r      => UAREA_R(1,1,1)
    barotr2_udel4_0_param%d_duc          => DUC(1,1,1)
    barotr2_udel4_0_param%d_dum          => DUM(1,1,1)
    barotr2_udel4_0_param%d_dun          => DUN(1,1,1)
    barotr2_udel4_0_param%d_dus          => DUS(1,1,1)
    barotr2_udel4_0_param%d_due          => DUE(1,1,1)
    barotr2_udel4_0_param%d_duw          => DUW(1,1,1)
    barotr2_udel4_0_param%d_dmc          => DMC(1,1,1)
    barotr2_udel4_0_param%d_dmn          => DMN(1,1,1)
    barotr2_udel4_0_param%d_dms          => DMS(1,1,1)
    barotr2_udel4_0_param%d_dme          => DME(1,1,1)
    barotr2_udel4_0_param%d_dmw          => DMW(1,1,1)
    barotr2_udel4_0_param%d_amf          => AMF(1,1,1)
    barotr2_udel4_0_param%d_am_factor    => am_factor_d(1,1,1)
    barotr2_udel4_0_param%d_d2uk         => D2UK(1,1)
    barotr2_udel4_0_param%d_d2vk         => D2VK(1,1)

    barotr2_udel4_1_param%cuBlock(1)     = imt
    barotr2_udel4_1_param%cuBlock(2)     = 1
    barotr2_udel4_1_param%cuBlock(3)     = 1
    barotr2_udel4_1_param%cuGrid(1)      = 1
    barotr2_udel4_1_param%cuGrid(2)      = jmt
    barotr2_udel4_1_param%cuGrid(3)      = 1
    barotr2_udel4_1_param%d_jb           = this_block%JB
    barotr2_udel4_1_param%d_je           = this_block%JE
    barotr2_udel4_1_param%d_ib           = this_block%IB
    barotr2_udel4_1_param%d_ie           = this_block%IE
    barotr2_udel4_1_param%d_imt          = imt
    barotr2_udel4_1_param%d_jmt          = jmt
    barotr2_udel4_1_param%d_km           = km
    barotr2_udel4_1_param%rank           = my_task
    barotr2_udel4_1_param%d_kmu          => KMU(1,1,1)
    barotr2_udel4_1_param%d_am           = AM_BIHAR
    barotr2_udel4_1_param%d_dts          = DTS
    barotr2_udel4_1_param%d_value        = 1.7*(0.25*3.6D7/real(IMT_GLOBAL, kind=8))**2.5
    barotr2_udel4_1_param%d_viv          => VIV(1,1,1,1)
    barotr2_udel4_1_param%d_ebea         => EBEA(1,1,1)
    barotr2_udel4_1_param%d_ebeb         => EBEB(1,1,1)
    barotr2_udel4_1_param%d_fcor         => FCOR(1,1,1)
#ifdef TIDEFORCE
    barotr2_udel4_1_param%d_h0_tide      => h0_tide(1,1,1)
#endif 
    barotr2_udel4_1_param%d_h0           => H0(1,1,1)
    barotr2_udel4_1_param%d_ubp          => UBP(1,1,1)
    barotr2_udel4_1_param%d_vbp          => VBP(1,1,1)
    barotr2_udel4_1_param%d_dlub         => DLUB(1,1,1)
    barotr2_udel4_1_param%d_dlvb         => DLVB(1,1,1)
    barotr2_udel4_1_param%d_wka          => WKA(1,1,1,1)
    barotr2_udel4_1_param%d_work         => WORK(1,1,1)
    barotr2_udel4_1_param%d_wgp          => WGP(1,1,1)
    barotr2_udel4_1_param%d_pax          => PAX(1,1,1)
    barotr2_udel4_1_param%d_pay          => PAY(1,1,1)
    barotr2_udel4_1_param%d_pxb          => PXB(1,1,1)
    barotr2_udel4_1_param%d_pyb          => PYB(1,1,1)
    barotr2_udel4_1_param%d_whx          => WHX(1,1,1)
    barotr2_udel4_1_param%d_why          => WHY(1,1,1)
    barotr2_udel4_1_param%d_duc          => DUC(1,1,1)
    barotr2_udel4_1_param%d_dum          => DUM(1,1,1)
    barotr2_udel4_1_param%d_dun          => DUN(1,1,1)
    barotr2_udel4_1_param%d_dus          => DUS(1,1,1)
    barotr2_udel4_1_param%d_due          => DUE(1,1,1)
    barotr2_udel4_1_param%d_duw          => DUW(1,1,1)
    barotr2_udel4_1_param%d_dmc          => DMC(1,1,1)
    barotr2_udel4_1_param%d_dmn          => DMN(1,1,1)
    barotr2_udel4_1_param%d_dms          => DMS(1,1,1)
    barotr2_udel4_1_param%d_dme          => DME(1,1,1)
    barotr2_udel4_1_param%d_dmw          => DMW(1,1,1)
    barotr2_udel4_1_param%d_dxur         => DXUR(1,1,1)
    barotr2_udel4_1_param%d_dyur         => DYUR(1,1,1)
    barotr2_udel4_1_param%d_au0          => AU0(1,1,1)
    barotr2_udel4_1_param%d_ausw         => AUSW(1,1,1)
    barotr2_udel4_1_param%d_aus          => AUS(1,1,1)
    barotr2_udel4_1_param%d_auw          => AUW(1,1,1)
    barotr2_udel4_1_param%d_d2uk         => D2UK(1,1)
    barotr2_udel4_1_param%d_d2vk         => D2VK(1,1)
    barotr2_udel4_1_param%d_am_factor    => am_factor_d(1,1,1)
    barotr2_udel4_1_param%d_uarea_r      => uarea(1,1,1)
    barotr2_udel4_1_param%d_factor       => am_factor_d(1,1,1)
    barotr2_udel4_1_param%d_fac_vis_uv     => factor_vis_uv(1,1,1)
    barotr2_udel4_1_param%d_value2        =(real(IMT_GLOBAL, kind=8)/1.8D4)**1.3
     END DO
#else
    DO IBLOCK = 1, NBLOCKS_CLINIC
        this_block = get_block(blocks_clinic(iblock),iblock)
    barotr2_param%cuBlock(1)     = imt
    barotr2_param%cuBlock(2)     = 1
    barotr2_param%cuBlock(3)     = 1
    barotr2_param%cuGrid(1)      = (imt + barotr2_param%cuBlock(1) - 1) / barotr2_param%cuBlock(1)
    barotr2_param%cuGrid(2)      = (jmt + barotr2_param%cuBlock(2) - 1) / barotr2_param%cuBlock(2)
    barotr2_param%cuGrid(3)      = 1
    barotr2_param%d_jb           = this_block%JB
    barotr2_param%d_je           = this_block%JE
    barotr2_param%d_ib           = this_block%IB
    barotr2_param%d_ie           = this_block%IE
    barotr2_param%d_imt          = imt
    barotr2_param%d_jmt          = jmt
    barotr2_param%d_km           = km
    barotr2_param%rank           = my_task
    barotr2_param%d_kmu          => KMU(1,1,1)
    barotr2_param%d_am           = AM_HOR
    barotr2_param%d_viv          => VIV(1,1,1,1)
    barotr2_param%d_ebea         => EBEA(1,1,1)
    barotr2_param%d_ebeb         => EBEB(1,1,1)
    barotr2_param%d_fcor         => FCOR(1,1,1)
#ifdef TIDEFORCE
    barotr2_param%d_h0_tide      => h0_tide(1,1,1)
#endif
    barotr2_param%d_h0           => H0(1,1,1)
    barotr2_param%d_ubp          => UBP(1,1,1)
    barotr2_param%d_vbp          => VBP(1,1,1)
    barotr2_param%d_dlub         => DLUB(1,1,1)
    barotr2_param%d_dlvb         => DLVB(1,1,1)
    barotr2_param%d_wka          => WKA(1,1,1,1)
    barotr2_param%d_work         => WORK(1,1,1)
    barotr2_param%d_wgp          => WGP(1,1,1)
    barotr2_param%d_pax          => PAX(1,1,1)
    barotr2_param%d_pay          => PAY(1,1,1)
    barotr2_param%d_pxb          => PXB(1,1,1)
    barotr2_param%d_pyb          => PYB(1,1,1)
    barotr2_param%d_whx          => WHX(1,1,1)
    barotr2_param%d_why          => WHY(1,1,1)
    barotr2_param%d_duc          => DUC(1,1,1)
    barotr2_param%d_dum          => DUM(1,1,1)
    barotr2_param%d_dun          => DUN(1,1,1)
    barotr2_param%d_dus          => DUS(1,1,1)
    barotr2_param%d_due          => DUE(1,1,1)
    barotr2_param%d_duw          => DUW(1,1,1)
    barotr2_param%d_dmc          => DMC(1,1,1)
    barotr2_param%d_dmn          => DMN(1,1,1)
    barotr2_param%d_dms          => DMS(1,1,1)
    barotr2_param%d_dme          => DME(1,1,1)
    barotr2_param%d_dmw          => DMW(1,1,1)
    barotr2_param%d_dxur         => DXUR(1,1,1)
    barotr2_param%d_dyur         => DYUR(1,1,1)
    barotr2_param%d_au0          => AU0(1,1,1)
    barotr2_param%d_ausw         => AUSW(1,1,1)
    barotr2_param%d_aus          => AUS(1,1,1)
    barotr2_param%d_auw          => AUW(1,1,1)

        END DO
#endif

    DO IBLOCK = 1, NBLOCKS_CLINIC
        barotr2_tide_param%cuBlock(1)     = imt
        barotr2_tide_param%cuBlock(2)     = 1
        barotr2_tide_param%cuBlock(3)     = 1
        barotr2_tide_param%cuGrid(1)      = (imt + barotr2_tide_param%cuBlock(1) - 1) / barotr2_tide_param%cuBlock(1)
        barotr2_tide_param%cuGrid(2)      = (jmt + barotr2_tide_param%cuBlock(2) - 1) / barotr2_tide_param%cuBlock(2)
        barotr2_tide_param%cuGrid(3)      = 1

        barotr2_tide_param%d_d0       = D0
        barotr2_tide_param%d_wka      => WKA(1,1,1,1)
        barotr2_tide_param%d_ubp      => UBP(1,1,1)
        barotr2_tide_param%d_vbp      => VBP(1,1,1)
        barotr2_tide_param%d_ohbu     => OHBU(1,1,1)
        barotr2_tide_param%d_h0       => H0(1,1,1)
#ifdef TIDEFORCE
        barotr2_tide_param%d_cd_tide  => CD_TIDE(1,1,1)
#endif
        barotr2_tide_param%d_viv      => VIV(1,1,1,1)
        barotr2_tide_param%d_kmu      => KMU(1,1,1)
        barotr2_tide_param%d_ric      => ric(1,1,1,1)
#ifdef TIDEFORCEGW
        barotr2_tide_param%d_cd_inwave=> CD_InWave(1,1,1)
        barotr2_tide_param%d_depthstd => depthstd(1,1,1)
#endif 
    END DO



    barotr3_param%cuBlock(1)     = imt
    barotr3_param%cuBlock(2)     = 1
    barotr3_param%cuBlock(3)     = 1
    barotr3_param%cuGrid(1)      = 1
    barotr3_param%cuGrid(2)      = jmt
    barotr3_param%cuGrid(3)      = 1
    barotr3_param%d_imt          = imt
    barotr3_param%d_jmt          = jmt
    barotr3_param%rank           = my_task
    barotr3_param%d_dtb          = DTB
    barotr3_param%d_dzph         => DZPH(1,1,1)
    barotr3_param%d_ub           => UB(1,1,1)
    barotr3_param%d_vb           => VB(1,1,1)
    barotr3_param%d_ubp          => UBP(1,1,1)
    barotr3_param%d_vbp          => VBP(1,1,1)
    barotr3_param%d_wka          => WKA(1,1,1,1)
    barotr3_param%d_work         => WORK(1,1,1)

    DO IBLOCK = 1, NBLOCKS_CLINIC
        this_block = get_block(blocks_clinic(iblock),iblock)
#ifdef BIHAR

    barotr3_tdel4_1_param%cuBlock(1)     = imt
    barotr3_tdel4_1_param%cuBlock(2)     = 1
    barotr3_tdel4_1_param%cuBlock(3)     = 1
    barotr3_tdel4_1_param%cuGrid(1)      = 1
    barotr3_tdel4_1_param%cuGrid(2)      = jmt
    barotr3_tdel4_1_param%cuGrid(3)      = 1
    barotr3_tdel4_1_param%d_jb           = this_block%JB
    barotr3_tdel4_1_param%d_je           = this_block%JE
    barotr3_tdel4_1_param%d_ib           = this_block%IB
    barotr3_tdel4_1_param%d_ie           = this_block%IE
    barotr3_tdel4_1_param%d_imt          = imt
    barotr3_tdel4_1_param%d_jmt          = jmt
    barotr3_tdel4_1_param%d_km           = km
    barotr3_tdel4_1_param%rank           = my_task
    barotr3_tdel4_1_param%d_kmt          => KMT(1,1,1)
    barotr3_tdel4_1_param%d_kmtn         => KMTN(1,1,1)
    barotr3_tdel4_1_param%d_kmts         => KMTS(1,1,1)
    barotr3_tdel4_1_param%d_kmte         => KMTE(1,1,1)
    barotr3_tdel4_1_param%d_kmtw         => KMTW(1,1,1)
    barotr3_tdel4_1_param%d_h0p          => H0P(1,1,1)
    barotr3_tdel4_1_param%d_dtn          => DTN(1,1,1)
    barotr3_tdel4_1_param%d_dts          => DTS_s(1,1,1)
    barotr3_tdel4_1_param%d_dte          => DTE(1,1,1)
    barotr3_tdel4_1_param%d_dtw          => DTW(1,1,1)
    barotr3_tdel4_1_param%d_d2tk         => dt2k_b(1,1)
    barotr3_tdel4_1_param%d_ahf          => AHF(1,1,1)


    barotr3_1_param%d_ah           = AH_BIHAR
#else

    barotr3_1_param%d_ah           = AH_HOR
#endif

    barotr3_1_param%cuBlock(1)     = imt
    barotr3_1_param%cuBlock(2)     = 1
    barotr3_1_param%cuBlock(3)     = 1
    barotr3_1_param%cuGrid(1)      = 1
    barotr3_1_param%cuGrid(2)      = jmt
    barotr3_1_param%cuGrid(3)      = 1
    barotr3_1_param%d_jb           = this_block%JB
    barotr3_1_param%d_je           = this_block%JE
    barotr3_1_param%d_ib           = this_block%IB
    barotr3_1_param%d_ie           = this_block%IE
    barotr3_1_param%d_imt          = imt
    barotr3_1_param%d_jmt          = jmt
    barotr3_1_param%d_km           = km
    barotr3_1_param%rank           = my_task
    barotr3_1_param%d_kmt          => KMT(1,1,1)
    barotr3_1_param%d_kmtn         => KMTN(1,1,1)
    barotr3_1_param%d_kmts         => KMTS(1,1,1)
    barotr3_1_param%d_kmte         => KMTE(1,1,1)
    barotr3_1_param%d_kmtw         => KMTW(1,1,1)

    barotr3_1_param%d_vit          => VIT(1,1,1,1)
    barotr3_1_param%d_h0p          => H0P(1,1,1)
    barotr3_1_param%d_wka          => WKA(1,1,1,1)
    barotr3_1_param%d_work         => WORK(1,1,1)
    barotr3_1_param%d_htw          => HTW(1,1,1)
    barotr3_1_param%d_hts          => HTS(1,1,1)
    barotr3_1_param%d_tarea_r      => TAREA_R(1,1,1)
    barotr3_1_param%d_dtn          => DTN(1,1,1)
    barotr3_1_param%d_dts          => DTS_s(1,1,1)
    barotr3_1_param%d_dte          => DTE(1,1,1)
    barotr3_1_param%d_dtw          => DTW(1,1,1)
    barotr3_1_param%d_d2tk         => dt2k_b(1,1)
    barotr3_1_param%d_factor_vis_tsz0          => factor_vis_tsz0(1,1,1)
    end do
    barotr4_param%cuBlock(1)     = imt
    barotr4_param%cuBlock(2)     = 1
    barotr4_param%cuBlock(3)     = 1
    barotr4_param%cuGrid(1)      = 1
    barotr4_param%cuGrid(2)      = jmt
    barotr4_param%cuGrid(3)      = 1
    barotr4_param%d_imt          = imt
    barotr4_param%d_jmt          = jmt
    barotr4_param%d_dtb          = DTB
    barotr4_param%d_iau          = enoi_iau
    barotr4_param%d_lamda        = 1.0D0/(iau_gamma*86400.D0)
    barotr4_param%d_h0           => H0(1,1,1)
    barotr4_param%d_h0p          => H0P(1,1,1)
    barotr4_param%d_ub           => UB(1,1,1)
    barotr4_param%d_vb           => VB(1,1,1)
    barotr4_param%d_ubp          => UBP(1,1,1)
    barotr4_param%d_vbp          => VBP(1,1,1)
    barotr4_param%d_h0f          => H0F(1,1,1)
    barotr4_param%d_h0bf         => H0BF(1,1,1)
    barotr4_param%d_work         => WORK(1,1,1)
    barotr4_param%d_h0_iau       => h0_iau(1,1,1)
    barotr4_param%d_vit          => VIT(1,1,1,1)
end subroutine barotr_sw_init

subroutine bclinc_sw_init
    use pconst_mod
    use tracer_mod
    use dyn_mod
    use work_mod
    use pmix_mod
    use constant_mod
    !    use forc_mod, only: psa,USTAR,BUOYTUR, BUOYSOL,NSWV,SWV
    use forc_mod
    use grid
    use blocks
    use domain
    implicit none


#if (defined CANUTO) ||(defined CANUTO2010)
    bclinc12_cu_param%cuBlock(1) = imt
    bclinc12_cu_param%cuBlock(2) = 1
    bclinc12_cu_param%cuBlock(3) = 1
    bclinc12_cu_param%cuGrid(1)  = 1
    bclinc12_cu_param%cuGrid(2)  = jmt
    bclinc12_cu_param%cuGrid(3)  = 1
    bclinc12_cu_param%d_od0 = od0
    bclinc12_cu_param%d_c0f = c0f
    bclinc12_cu_param%d_cag = cag
    bclinc12_cu_param%d_snlat => snlat(1,1,1)
    bclinc12_cu_param%d_sag = sag
    bclinc12_cu_param%d_up => up(1,1,1,1)
    bclinc12_cu_param%d_vp => vp(1,1,1,1)
    bclinc12_cu_param%d_bbcy => bbcy(1,1,1)
    bclinc12_cu_param%d_sbcy => sbcy(1,1,1)
    bclinc12_cu_param%d_sbcx => sbcx(1,1,1)
    bclinc12_cu_param%d_bbcx => bbcx(1,1,1)
    bclinc12_cu_param%d_su => su(1,1,1)
    bclinc12_cu_param%d_sv => sv(1,1,1)
    bclinc12_cu_param%d_kmu => kmu(1,1,1)
#endif


    bclinc11_cu_param%cuBlock(1) = imt
    bclinc11_cu_param%cuBlock(2) = 1
    bclinc11_cu_param%cuBlock(3) = 1
    bclinc11_cu_param%cuGrid(1)  = (imt + bclinc11_cu_param%cuBlock(1) - 1) / bclinc11_cu_param%cuBlock(1)
    bclinc11_cu_param%cuGrid(2)  = (jmt + bclinc11_cu_param%cuBlock(2) - 1) / bclinc11_cu_param%cuBlock(2)
    bclinc11_cu_param%cuGrid(3)  = 1
    bclinc11_cu_param%d_od0      = od0
    bclinc11_cu_param%d_onbb     = onbb
    bclinc11_cu_param%d_ohbt     => ohbt(1,1,1)
    bclinc11_cu_param%d_vit      => vit(1, 1, 1, 1)
    bclinc11_cu_param%d_dzp      => dzp(1)
    bclinc11_cu_param%d_h0bf     => H0BF(1,1,1)
    bclinc11_cu_param%d_h0bl     => h0bl(1,1,1)
    bclinc11_cu_param%d_gg       => gg(1,1,1,1)
    bclinc11_cu_param%d_psa      => psa(1,1,1)
    bclinc11_cu_param%d_wka      => wka(1,1,1,1)
    bclinc11_cu_param%d_work     => work(1,1,1)
    bclinc11_cu_param%d_wka1     => wka1(1,1,1,1)
    bclinc11_cu_param%d_wka2     => wka2(1,1,1,1)
    bclinc11_cu_param%d_zkt      => zkt(1)

    bclinc6_cu_param%cuBlock(1) = imt
    bclinc6_cu_param%cuBlock(2) = 1
    bclinc6_cu_param%cuBlock(3) = 1
    bclinc6_cu_param%cuGrid(1)  = 1
    bclinc6_cu_param%cuGrid(2)  = jmt
    bclinc6_cu_param%cuGrid(3)  = km
    bclinc6_cu_param%d_viv      => viv(1, 1, 1, 1)
    bclinc6_cu_param%d_epea     => epea(1, 1, 1)
    bclinc6_cu_param%d_epeb     => epeb(1, 1, 1)
    bclinc6_cu_param%d_epla     => epla(1, 1, 1)
    bclinc6_cu_param%d_eplb     => eplb(1, 1, 1)
    bclinc6_cu_param%d_up       => up(1, 1, 1, 1)
    bclinc6_cu_param%d_vp       => vp(1, 1, 1, 1)
    bclinc6_cu_param%d_dlu      => dlu(1, 1, 1, 1)
    bclinc6_cu_param%d_dlv      => dlv(1, 1, 1, 1)
    bclinc6_cu_param%d_gg       => gg(1, 1, 1, 1)
    bclinc6_cu_param%d_kmu      => kmu(1, 1, 1)
    bclinc6_cu_param%d_fcor     => fcor(1, 1, 1)
    bclinc6_cu_param%d_wka1     => wka1(1,1, 1, 1)
    bclinc6_cu_param%d_wka2     => wka(1,1, 1, 1)
    bclinc6_cu_param%d_dxur     => dxur(1, 1, 1)
    bclinc6_cu_param%d_dyur     => dyur(1, 1, 1)



    bclinc71_cu_param%cuBlock(1) = imt
    bclinc71_cu_param%cuBlock(2) = 1
    bclinc71_cu_param%cuBlock(3) = 1
    bclinc71_cu_param%cuGrid(1)  = 1
    bclinc71_cu_param%cuGrid(2)  = jmt
    bclinc71_cu_param%cuGrid(3)  = km
    bclinc71_cu_param%d_viv      => viv(1, 1, 1, 1)
    bclinc71_cu_param%d_odzt     => odzt(1)
    bclinc71_cu_param%d_odzp     => odzp(1)
    bclinc71_cu_param%d_akmu     => akmu(1, 1, 1, 1)
    bclinc71_cu_param%d_dtc      =  DTC
    bclinc71_cu_param%rank       =  my_task
    bclinc71_cu_param%d_u        => u(1, 1, 1, 1)
    bclinc71_cu_param%d_up       => up(1, 1, 1, 1)
    bclinc71_cu_param%d_dlu      => dlu(1, 1, 1, 1)
    bclinc71_cu_param%d_v        => v(1, 1, 1, 1)
    bclinc71_cu_param%d_vp       => vp(1, 1, 1, 1)
    bclinc71_cu_param%d_dlv      => dlv(1, 1, 1, 1)
    bclinc71_cu_param%d_sbcx     => sbcx(1, 1, 1)
    bclinc71_cu_param%d_bbcx     => bbcx(1, 1, 1)
    bclinc71_cu_param%d_kmu      => kmu(1, 1, 1)

    !      bclinc72_cu_param%cuBlock(1) = TWOX
    !      bclinc72_cu_param%cuBlock(2) = TWOY
    !      bclinc72_cu_param%cuBlock(3) = 1
    !      bclinc72_cu_param%cuGrid(1)  = (imt + bclinc72_cu_param%cuBlock(1) - 1) / bclinc72_cu_param%cuBlock(1)
    !      bclinc72_cu_param%cuGrid(2)  = (jmt + bclinc72_cu_param%cuBlock(2) - 1) / bclinc72_cu_param%cuBlock(2)
    !      bclinc72_cu_param%cuGrid(3)  = 1
    !      bclinc72_cu_param%d_viv      => viv(1, 1, 1, 1)
    !      bclinc72_cu_param%d_odzt     => odzt(1)
    !      bclinc72_cu_param%d_odzp     => odzp(1)
    !      bclinc72_cu_param%d_akmu     => akmu(1, 1, 1, 1)
    !      bclinc72_cu_param%d_dtc      =  DTC
    !      bclinc72_cu_param%d_v        => v(1, 1, 1, 1)
    !      bclinc72_cu_param%d_vp       => vp(1, 1, 1, 1)
    !      bclinc72_cu_param%d_dlv      => dlv(1, 1, 1, 1)
    !      bclinc72_cu_param%d_sbcy     => sbcy(1, 1, 1)
    !      bclinc72_cu_param%d_bbcy     => bbcy(1, 1, 1)
    !      bclinc72_cu_param%d_kmu      => kmu(1, 1, 1)
    !      bclinc72_cu_param%d_a8       => ta8(1,1,1)
    !      bclinc72_cu_param%d_b8       => tb8(1,1,1)
    !      bclinc72_cu_param%d_c8       => tc8(1,1,1)
    !      bclinc72_cu_param%d_d8       => td8(1,1,1)
    !      bclinc72_cu_param%d_e8       => te8(1,1,1)
    !      bclinc72_cu_param%d_f8       => tf8(1,1,1)
    !      bclinc72_cu_param%d_wka2     => wka2(1,1,1,1)


    bclinc81_cu_param%cuBlock(1) = imt
    bclinc81_cu_param%cuBlock(2) = 1
    bclinc81_cu_param%cuBlock(3) = 1
    bclinc81_cu_param%cuGrid(1)  = 1
    bclinc81_cu_param%cuGrid(2)  = jmt
    bclinc81_cu_param%cuGrid(3)  = km
    bclinc81_cu_param%d_viv      => viv(1, 1, 1, 1)
    bclinc81_cu_param%d_dzp      => dzp(1)
    bclinc81_cu_param%d_u1       => wka1(1,1,1,1) !jenny
    bclinc81_cu_param%d_u        => u(1, 1, 1, 1)
    bclinc81_cu_param%d_ub       => ub(1, 1, 1)
    bclinc81_cu_param%d_utf      => utf(1, 1, 1, 1)
    bclinc81_cu_param%d_ohbu     => ohbu(1, 1, 1)
    bclinc81_cu_param%d_work      =>work(1,1,1)

    bclinc82_cu_param%cuBlock(1) = imt
    bclinc82_cu_param%cuBlock(2) = 1
    bclinc82_cu_param%cuBlock(3) = 1
    bclinc82_cu_param%cuGrid(1)  = 1
    bclinc82_cu_param%cuGrid(2)  = jmt
    bclinc82_cu_param%cuGrid(3)  = km
    bclinc82_cu_param%d_viv      => viv(1, 1, 1, 1)
    bclinc82_cu_param%d_dzp      => dzp(1)
    bclinc82_cu_param%d_v1     => wka2(1,1,1,1) !jenny
    bclinc82_cu_param%d_v        => v(1, 1, 1, 1)
    bclinc82_cu_param%d_vb       => vb(1, 1, 1)
    bclinc82_cu_param%d_vtf      => vtf(1, 1, 1, 1)
    bclinc82_cu_param%d_ohbu     => ohbu(1, 1, 1)
    bclinc82_cu_param%d_work     => work(1,1,1)

    bclinc2_cu_param%cuBlock(1) = imt
    bclinc2_cu_param%cuBlock(2) = 1
    bclinc2_cu_param%cuBlock(3) = 1
    bclinc2_cu_param%cuGrid(1)  = 1
    bclinc2_cu_param%cuGrid(2)  = jmt
    bclinc2_cu_param%cuGrid(3)  = km
    bclinc2_cu_param%d_wka        => wka(1,1,1,1)
    bclinc2_cu_param%d_vp         => vp(1,1,1,1)
    bclinc2_cu_param%d_dlv        => dlv(1,1,1,1)
    bclinc2_cu_param%d_dtc2       = dtc2

    bclinc3_cu_param%cuBlock(1) = imt
    bclinc3_cu_param%cuBlock(2) = 1
    bclinc3_cu_param%cuBlock(3) = 1
    bclinc3_cu_param%cuGrid(1)  = (imt + bclinc3_cu_param%cuBlock(1) - 1) / bclinc3_cu_param%cuBlock(1)
    bclinc3_cu_param%cuGrid(2)  = (jmt + bclinc3_cu_param%cuBlock(2) - 1) / bclinc3_cu_param%cuBlock(2)
    bclinc3_cu_param%cuGrid(3)  = 1
    bclinc3_cu_param%d_kmu      => kmu(1, 1, 1)
    bclinc3_cu_param%d_viv      => viv(1, 1, 1, 1)
    bclinc3_cu_param%d_dcb      => akmu(1, 1, 1, 1)
    bclinc3_cu_param%d_dtc2     =  dtc2
    bclinc3_cu_param%d_odzt     => odzt(1)
    bclinc3_cu_param%d_odzp     => odzp(1)


    bclinc_invtriu_param%cuBlock(1) = imt
    bclinc_invtriu_param%cuBlock(2) = 1
    bclinc_invtriu_param%cuBlock(3) = 1
    bclinc_invtriu_param%cuGrid(1)  = (imt + bclinc_invtriu_param%cuBlock(1) - 1) / bclinc_invtriu_param%cuBlock(1)
    bclinc_invtriu_param%cuGrid(2)  = (jmt + bclinc_invtriu_param%cuBlock(2) - 1) / bclinc_invtriu_param%cuBlock(2)
    bclinc_invtriu_param%cuGrid(3)  = 1
    bclinc_invtriu_param%d_kmu      => kmu(1, 1, 1)
    bclinc_invtriu_param%d_viv      => viv(1, 1, 1, 1)
    bclinc_invtriu_param%d_akmu     => akmu(1, 1, 1, 1)
    bclinc_invtriu_param%d_dtc2     =  DTC
    bclinc_invtriu_param%d_bbcx     => bbcx(1, 1, 1)
    bclinc_invtriu_param%d_bbcy     => bbcy(1, 1, 1)
    bclinc_invtriu_param%d_sbcx     => sbcx(1, 1, 1)
    bclinc_invtriu_param%d_sbcy     => sbcy(1, 1, 1)
    bclinc_invtriu_param%d_u        => u(1,1,1,1)
    !bclinc_invtriu_param%d_wk1      => u(1,1,1,1)
    bclinc_invtriu_param%d_v        => v(1,1,1,1)
    !bclinc_invtriu_param%d_wk2      => v(1,1,1,1)
    bclinc_invtriu_param%d_odzt     => odzt(1)
    bclinc_invtriu_param%d_odzp     => odzp(1)


    bclinc4_cu_param%cuBlock(1) = imt
    bclinc4_cu_param%cuBlock(2) = 1
    bclinc4_cu_param%cuBlock(3) = 1
    bclinc4_cu_param%cuGrid(1)  = 1
    bclinc4_cu_param%cuGrid(2)  = jmt
    bclinc4_cu_param%cuGrid(3)  = km
    bclinc4_cu_param%d_ohbu     => ohbu(1, 1, 1)
    bclinc4_cu_param%d_viv      => viv(1, 1, 1, 1)
    bclinc4_cu_param%d_dzp      => dzp(1)
    bclinc4_cu_param%d_odzt     => odzt(1)
    bclinc4_cu_param%d_odzp     => odzp(1)
    bclinc4_cu_param%d_akmu     => akmu(1, 1, 1, 1)
    bclinc4_cu_param%d_afc1     =  AFC1
    bclinc4_cu_param%d_afc2     =  afc2
    bclinc4_cu_param%d_dtc2     =  dtc2
    bclinc4_cu_param%rank       =  my_task
    bclinc4_cu_param%d_v        => v(1, 1, 1, 1)
    bclinc4_cu_param%d_up       => up(1, 1, 1, 1)
    bclinc4_cu_param%d_vp       => vp(1, 1, 1, 1)
    bclinc4_cu_param%d_vb       => vb(1, 1, 1)
    bclinc4_cu_param%d_dlu      => dlu(1,1,1,1)
    bclinc4_cu_param%d_sbcx     => sbcx(1,1,1)
    bclinc4_cu_param%d_bbcx     => bbcx(1,1,1)
    bclinc4_cu_param%d_kmu      => kmu(1,1,1)
    bclinc4_cu_param%d_wka      => wka(1,1,1,1)
!    bclinc4_cu_param%d_a8       => ta8(1,1,1)
!    bclinc4_cu_param%d_b8       => tb8(1,1,1)
!    bclinc4_cu_param%d_c8       => tc8(1,1,1)
!    bclinc4_cu_param%d_d8       => td8(1,1,1)
!    bclinc4_cu_param%d_e8       => te8(1,1,1)
!    bclinc4_cu_param%d_f8       => tf8(1,1,1)
!    bclinc4_cu_param%d_work     => work(1,1,1)


    bclinc4_0_cu_param%cuBlock(1) = imt
    bclinc4_0_cu_param%cuBlock(2) = 1
    bclinc4_0_cu_param%cuBlock(3) = 1
    bclinc4_0_cu_param%cuGrid(1)  = 1
    bclinc4_0_cu_param%cuGrid(2)  = jmt
    bclinc4_0_cu_param%cuGrid(3)  = km
    bclinc4_0_cu_param%d_viv      => viv(1, 1, 1, 1)
    bclinc4_0_cu_param%d_afc1     =  AFC1
    bclinc4_0_cu_param%d_afc2     =  afc2
    bclinc4_0_cu_param%d_dtc2     =  dtc2
    bclinc4_0_cu_param%d_v        => v(1, 1, 1, 1)
    bclinc4_0_cu_param%d_up       => up(1, 1, 1, 1)
    bclinc4_0_cu_param%d_vp       => vp(1, 1, 1, 1)
    bclinc4_0_cu_param%d_vb       => vb(1, 1, 1)
    bclinc4_0_cu_param%d_dlu      => dlu(1,1,1,1)
    bclinc4_0_cu_param%d_wka      => wka(1,1,1,1)
    bclinc4_0_cu_param%d_work     => work(1,1,1)

    bclinc5_cu_param%cuBlock(1) = imt
    bclinc5_cu_param%cuBlock(2) = 1
    bclinc5_cu_param%cuBlock(3) = 1
    bclinc5_cu_param%cuGrid(1)  = 1
    bclinc5_cu_param%cuGrid(2)  = jmt
    bclinc5_cu_param%cuGrid(3)  = km
    bclinc5_cu_param%d_viv      => viv(1, 1, 1, 1)
    bclinc5_cu_param%d_dzp      => dzp(1)
    bclinc5_cu_param%d_afc1     =  AFC1
    bclinc5_cu_param%d_afc2     =  afc2
    bclinc5_cu_param%rank       =  my_task
    bclinc5_cu_param%d_lamda    =  1.0D0/(iau_gamma*86400.D0)
    bclinc5_cu_param%d_iau          = enoi_iau
    bclinc5_cu_param%d_u        => u(1, 1, 1, 1)
    bclinc5_cu_param%d_v        => v(1, 1, 1, 1)
    bclinc5_cu_param%d_up       => up(1, 1, 1, 1)
    bclinc5_cu_param%d_ub       => ub(1, 1, 1)
    bclinc5_cu_param%d_wka      => wka(1,1,1,1)
    bclinc5_cu_param%d_ohbu     => ohbu(1, 1, 1)
    bclinc5_cu_param%d_vtf      => vtf(1, 1, 1, 1)
    bclinc5_cu_param%d_utf      => utf(1, 1, 1, 1)
    bclinc5_cu_param%d_work     => work(1,1,1)
    bclinc5_cu_param%d_ucheck   => ucheck(1,1,1)
    bclinc5_cu_param%d_vcheck   => vcheck(1,1,1)
    bclinc5_cu_param%vv_iau     => vv_iau(1,1,1,1)
    bclinc5_cu_param%uu_iau     => uu_iau(1,1,1,1)

    vinteg_param%cuBlock(1) = imt
    vinteg_param%cuBlock(2) = 1
    vinteg_param%cuBlock(3) = 1
    vinteg_param%cuGrid(1)  = (imt + vinteg_param%cuBlock(1) - 1) / vinteg_param%cuBlock(1)
    vinteg_param%cuGrid(2)  = (jmt + vinteg_param%cuBlock(2) - 1) / vinteg_param%cuBlock(2)
    vinteg_param%cuGrid(3)  = 1
    vinteg_param%rank       =  my_task
    vinteg_param%d_u        => u(1, 1, 1, 1)
    vinteg_param%d_v        => v(1, 1, 1, 1)
    vinteg_param%d_wka      => wka(1,1,1,1)
    vinteg_param%d_work     => work(1,1,1)
    vinteg_param%d_dzp      => dzp(1)
    vinteg_param%d_ohbu     => ohbu(1, 1, 1)
    vinteg_param%d_viv      => viv(1, 1, 1, 1)
end subroutine bclinc_sw_init

subroutine tracer_sw_init
    use pconst_mod
    use tracer_mod
    use dyn_mod
    use work_mod
    use pmix_mod
    use constant_mod
    !    use forc_mod, only: psa,USTAR,BUOYTUR, BUOYSOL,NSWV,SWV
    use forc_mod
    use grid
    use blocks
    use domain
    implicit none

    tracer_upwell_0_param%cuBlock(1)     = imt
    tracer_upwell_0_param%cuBlock(2)     = 1
    tracer_upwell_0_param%cuBlock(3)     = 1
    tracer_upwell_0_param%cuGrid(1)      = 1
    tracer_upwell_0_param%cuGrid(2)      = jmt
    tracer_upwell_0_param%cuGrid(3)      = km
    tracer_upwell_0_param%rank           = my_task
    tracer_upwell_0_param%d_au0            =>AU0(1,1,1)
    tracer_upwell_0_param%d_aus            =>aus(1,1,1)
    tracer_upwell_0_param%d_auw            =>AUW(1,1,1)
    tracer_upwell_0_param%d_ausw           =>AUSW(1,1,1)
    tracer_upwell_0_param%d_stf            =>STF(1,1,1)
    tracer_upwell_0_param%d_ohbu           =>ohbu(1,1,1)
    tracer_upwell_0_param%d_wkd            =>wkd(1,1,1,1)
    tracer_upwell_0_param%d_uk             =>UK(1,1,1,1)
    tracer_upwell_0_param%d_wkb            =>wkb(1,1,1,1)
    tracer_upwell_0_param%d_vk             =>VK(1,1,1,1)
    tracer_upwell_0_param%d_work           =>WORK(1,1,1)

    tracer_upwell_param%cuBlock(1)     = imt
    tracer_upwell_param%cuBlock(2)     = 1
    tracer_upwell_param%cuBlock(3)     = 1
    tracer_upwell_param%cuGrid(1)      = 1
    tracer_upwell_param%cuGrid(2)      = jmt
    tracer_upwell_param%cuGrid(3)      = km
    tracer_upwell_param%rank           = my_task
    tracer_upwell_param%d_dzp            =>DZP(1)
    tracer_upwell_param%d_stf             =>STF(1,1,1)
    tracer_upwell_param%d_hts            =>hts(1,1,1)
    tracer_upwell_param%d_htw            =>HTW(1,1,1)
    tracer_upwell_param%d_ohbt           =>OHBT(1,1,1)
    tracer_upwell_param%d_tarea_r        =>TAREA_R(1,1,1)
    tracer_upwell_param%d_uk             =>uk(1,1,1,1)
    tracer_upwell_param%d_vit            =>VIT(1,1,1,1)
    tracer_upwell_param%d_vk             =>VK(1,1,1,1)
    tracer_upwell_param%d_wka            =>wka(1,1,1,1)
    tracer_upwell_param%d_work           =>WORK(1,1,1)
    tracer_upwell_param%d_ws             =>ws(1,1,1,1)
    tracer_upwell_param%d_kmt            =>KMT(1,1,1)

    tracer_upwell_1_param%cuBlock(1)     = imt
    tracer_upwell_1_param%cuBlock(2)     = 1
    tracer_upwell_1_param%cuBlock(3)     = 1
    tracer_upwell_1_param%cuGrid(1)      = (imt + tracer_upwell_1_param%cuBlock(1) - 1) / tracer_upwell_1_param%cuBlock(1)
    tracer_upwell_1_param%cuGrid(2)      = (jmt + tracer_upwell_1_param%cuBlock(2) - 1) / tracer_upwell_1_param%cuBlock(2)
    tracer_upwell_1_param%cuGrid(3)      = 1
    tracer_upwell_1_param%rank           = my_task
    tracer_upwell_1_param%d_vit            =>VIT(1,1,1,1)
    tracer_upwell_1_param%d_dzp            =>DZP(1)
    tracer_upwell_1_param%d_wka            =>wka(1,1,1,1)
    tracer_upwell_1_param%d_work           =>WORK(1,1,1)

    tracer_upwell_2_param%cuBlock(1)     = imt
    tracer_upwell_2_param%cuBlock(2)     = 1
    tracer_upwell_2_param%cuBlock(3)     = 1
    tracer_upwell_2_param%cuGrid(1)      = (imt + tracer_upwell_2_param%cuBlock(1) - 1) / tracer_upwell_2_param%cuBlock(1)
    tracer_upwell_2_param%cuGrid(2)      = (jmt + tracer_upwell_2_param%cuBlock(2) - 1) / tracer_upwell_2_param%cuBlock(2)
    tracer_upwell_2_param%cuGrid(3)      = (km + tracer_upwell_2_param%cuBlock(3) - 1) / tracer_upwell_2_param%cuBlock(3)
    tracer_upwell_2_param%rank           = my_task
    tracer_upwell_2_param%d_ws             =>ws(1,1,1,1)
    tracer_upwell_2_param%d_stf            =>STF(1,1,1)
    tracer_upwell_2_param%d_dzp            =>DZP(1)
    tracer_upwell_2_param%d_ohbt           =>OHBT(1,1,1)
    tracer_upwell_2_param%d_work           =>WORK(1,1,1)
    tracer_upwell_2_param%d_wsbak          =>wsbak(1,1,1,1)
    tracer_upwell_2_param%d_wka            =>wka(1,1,1,1)
    tracer_upwell_2_param%d_vit            =>VIT(1,1,1,1)

    advection1_0_param%cuBlock(1) = imt
    advection1_0_param%cuBlock(2) = 1
    advection1_0_param%cuBlock(3) = 1
    advection1_0_param%cuGrid(1) = 1
    advection1_0_param%cuGrid(2) = jmt
    advection1_0_param%cuGrid(3) = km
    if ( adv_tracer(1:8) == 'centered') then
        advection1_0_param%Str = 0
    else if (adv_tracer(1:5) == 'tspas') then
        advection1_0_param%Str = 2
    else if (adv_tracer(1:4) == 'flux' ) then
        advection1_0_param%Str = 1
    endif
    advection1_0_param%rank = my_task
    advection1_0_param%d_wkb =>wkb(1,1,1,1)
    advection1_0_param%d_wkd =>wkd(1,1,1,1)
    advection1_0_param%d_u_wface =>u_wface(1,1,1)
    advection1_0_param%d_v_sface =>v_sface(1,1,1)
    advection1_0_param%d_hts =>hts(1,1,1)
    advection1_0_param%d_htw =>htw(1,1,1)
    advection1_0_param%d_dxu =>dxu(1,1,1)
    advection1_0_param%d_dyu =>dyu(1,1,1)
    advection1_0_param%d_adv_tt =>adv_tt(1,1,1)
    advection1_0_param%d_tf =>TF(1,1,1,1)

    advection1_param%cuBlock(1) = imt
    advection1_param%cuBlock(2) = 1
    advection1_param%cuBlock(3) = 1
    advection1_param%cuGrid(1) = 1
    advection1_param%cuGrid(2) = jmt
    advection1_param%cuGrid(3) = km
    advection1_param%Str = advection1_0_param%Str
    advection1_param%d_nss = NSS
    advection1_param%rank = my_task
    advection1_param%d_u_wface =>u_wface(1,1,1)
    advection1_param%d_v_sface =>v_sface(1,1,1)
    advection1_param%d_ws =>ws(1,1,1,1)
    advection1_param%d_at =>at(1,1,1,1,1)
    advection1_param%d_adv_tt =>adv_tt(1,1,1)
    advection1_param%d_tarea_r =>tarea_r(1,1,1)
    advection1_param%d_odzp =>ODZP(1)
    advection1_param%d_ax =>ax(1,1,1,1,1)
    advection1_param%d_ay =>ay(1,1,1,1,1)
    advection1_param%d_az =>az(1,1,1,1,1)

    advection2_param%cuBlock(1) = imt
    advection2_param%cuBlock(2) = 1
    advection2_param%cuBlock(3) = 1
    advection2_param%cuGrid(1) = 1
    advection2_param%cuGrid(2) = jmt
    advection2_param%cuGrid(3) = km
    advection2_param%d_ntra    =ntra
    advection2_param%Str = advection1_0_param%Str
    advection2_param%rank = my_task
    advection2_param%d_wkb                =>wkb(1,1,1,1)
    advection2_param%d_wkd                =>wkd(1,1,1,1)
    advection2_param%d_u_wface            =>u_wface(1,1,1)
    advection2_param%d_v_sface            =>v_sface(1,1,1)
    advection2_param%d_hts                =>hts(1,1,1)
    advection2_param%d_htw                =>htw(1,1,1)
    advection2_param%d_dxu                =>dxu(1,1,1)
    advection2_param%d_dyu                =>dyu(1,1,1)
    advection2_param%d_ws                 =>ws(1,1,1,1)
    advection2_param%d_adv_tt             =>adv_tt(1,1,1)
    advection2_param%d_tarea_r            =>tarea_r(1,1,1)
    advection2_param%d_odzp               =>ODZP(1)
    advection2_param%d_ax                 =>ax(1,1,1,1,1)
    advection2_param%d_ay                 =>ay(1,1,1,1,1)
    advection2_param%d_az                 =>az(1,1,1,1,1)
    advection2_param%d_hun                =>hun(1,1,1)
    advection2_param%d_hue                =>hue(1,1,1)
    advection2_param%d_at00               =>at00(1,1,1)
    advection2_param%d_atmax              =>atmax(1,1,1)
    advection2_param%d_atmin              =>atmin(1,1,1)
    advection2_param%d_odzt               =>ODZT(1)
    advection2_param%d_vit                =>vit(1,1,1,1)
    advection2_param%d_dxt                =>DXT(1,1,1)
    advection2_param%d_dyt                =>DYT(1,1,1)
    advection2_param%d_dxtr               =>DXTR(1,1,1)
    advection2_param%d_dytr               =>DYTR(1,1,1)
    advection2_param%d_adv_x              =>adv_x(1,1,1)
    advection2_param%d_adv_y              =>adv_y(1,1,1)
    advection2_param%d_adv_z              =>adv_z(1,1,1)
!    advection2_param%d_flxqRtmp           =>flxqRtmp(1,1,1)
!    advection2_param%d_flxqLtmp           =>flxqLtmp(1,1,1)
    advection2_param%d_zkp                => ZKP(1)
    advection2_param%d_dts                =dts

#if (defined BIHAR)

    tracer_del4_1_param%cuBlock(1)     = imt
    tracer_del4_1_param%cuBlock(2)     = 1
    tracer_del4_1_param%cuBlock(3)     = 1
    tracer_del4_1_param%cuGrid(1)      = 1
    tracer_del4_1_param%cuGrid(2)      = jmt
    tracer_del4_1_param%cuGrid(3)      = km
    tracer_del4_1_param%d_kmtn         => KMTN(1,1,1)
    tracer_del4_1_param%d_kmts         => KMTS(1,1,1)
    tracer_del4_1_param%d_kmte         => KMTE(1,1,1)
    tracer_del4_1_param%d_kmtw         => KMTW(1,1,1)
    tracer_del4_1_param%d_kmt          => KMT(1,1,1)
    tracer_del4_1_param%d_atb          => ATB(1,1,0,1,1)
    tracer_del4_1_param%d_dtn          => DTN(1,1,1)
    tracer_del4_1_param%d_dts          => DTS_s(1,1,1)
    tracer_del4_1_param%d_dte          => DTE(1,1,1)
    tracer_del4_1_param%d_dtw          => DTW(1,1,1)
    tracer_del4_1_param%d_d2tk         => h_d2tk(1,1,1)
    tracer_del4_1_param%d_ahf          => AHF(1,1,1)
#endif


    tracer_del4_2_param%cuBlock(1)     = imt
    tracer_del4_2_param%cuBlock(2)     = 1
    tracer_del4_2_param%cuBlock(3)     = 1
    tracer_del4_2_param%cuGrid(1)      = 1
    tracer_del4_2_param%cuGrid(2)      = jmt
    tracer_del4_2_param%cuGrid(3)      = km
    tracer_del4_2_param%d_kmtn         => KMTN(1,1,1)
    tracer_del4_2_param%d_kmts         => KMTS(1,1,1)
    tracer_del4_2_param%d_kmte         => KMTE(1,1,1)
    tracer_del4_2_param%d_kmtw         => KMTW(1,1,1)
    tracer_del4_2_param%d_kmt          => KMT(1,1,1)
#ifdef BIHAR
    tracer_del4_2_param%d_ah           = AH_BIHAR
#else
    tracer_del4_2_param%d_ah           = AH_HOR
#endif
    tracer_del4_2_param%d_dwndmix      = DWNDMIX
    tracer_del4_2_param%d_tf           => TF(1,1,1,1)
!    tracer_del4_2_param%d_dx           => dx(1,1,1,1,1)
    tracer_del4_2_param%d_adv_tt       => adv_tt(1,1,1)
    tracer_del4_2_param%d_vit          => VIT(1,1,1,1)
    tracer_del4_2_param%d_akt          => AKT(1,1,1,1,1)
    tracer_del4_2_param%d_wkc          => WKC(1,1,1,1)
    tracer_del4_2_param%d_dtn          => DTN(1,1,1)
    tracer_del4_2_param%d_dts          => DTS_s(1,1,1)
    tracer_del4_2_param%d_dte          => DTE(1,1,1)
    tracer_del4_2_param%d_dtw          => DTW(1,1,1)
    tracer_del4_2_param%d_d2tk         => h_d2tk(1,1,1)
    tracer_del4_2_param%d_atb          => ATB(1,1,0,1,1)
    tracer_del4_2_param%d_factor_vis_tsz0          => factor_vis_tsz0(1,1,1)

    invtrit_param%cuBlock(1)     = imt
    invtrit_param%cuBlock(2)     = 1
    invtrit_param%cuBlock(3)     = 1
    invtrit_param%cuGrid(1)      = 1
    invtrit_param%cuGrid(2)      = jmt
    invtrit_param%cuGrid(3)      = 1
    invtrit_param%rank           = my_task
    invtrit_param%d_kmt          => KMT(1,1,1)
    invtrit_param%d_wkc          => WKC(1,1,1,1)
    invtrit_param%d_vtl          => VTL(1,1,1,1)
    invtrit_param%d_stf          => STF(1,1,1)
    invtrit_param%d_vit          => VIT(1,1,1,1)
    invtrit_param%d_odzt         => ODZT(1)
    invtrit_param%d_odzp         => odzp(1)

    tracer_9_param%cuBlock(1)     = imt
    tracer_9_param%cuBlock(2)     = 1
    tracer_9_param%cuBlock(3)     = 1
    tracer_9_param%cuGrid(1)      = 1
    tracer_9_param%cuGrid(2)      = jmt
    tracer_9_param%cuGrid(3)      = km
    tracer_9_param%rank           = my_task
    tracer_9_param%d_at           => AT(1,1,1,1,1)
    tracer_9_param%d_vtl          => VTL(1,1,1,1)
    tracer_9_param%d_atb          => ATB(1,1,0,1,1)
    tracer_9_param%d_tend         => tend(1,1,1,1,1)
    tracer_9_param%d_vit          => VIT(1,1,1,1)
    !tracer_9_param%d_dt_diff      => dt_diff(1,1,1,1,1)
    tracer_9_param%d_ori         => VTL_ori(1,1,1,1)
    tracer_9_param%d_stf         => STF(1,1,1)
    tracer_9_param%d_odzp         => odzp(1)

    tracer_vtl_parm%cuBlock(1)     = imt
    tracer_vtl_parm%cuBlock(2)     = 1
    tracer_vtl_parm%cuBlock(3)     = 1
    tracer_vtl_parm%cuGrid(1)      = 1
    tracer_vtl_parm%cuGrid(2)      = jmt
    tracer_vtl_parm%cuGrid(3)      = km
    tracer_vtl_parm%d_vtl          => VTL(1,1,1,1)
    tracer_vtl_parm%d_vtl_ori      => vtl_ori(1,1,1,1)

    tracer_wkc_parm%cuBlock(1)     = imt
    tracer_wkc_parm%cuBlock(2)     = 1
    tracer_wkc_parm%cuBlock(3)     = 1
    tracer_wkc_parm%cuGrid(1)      = 1
    tracer_wkc_parm%cuGrid(2)      = jmt
    tracer_wkc_parm%cuGrid(3)      = km
    tracer_wkc_parm%d_ahv          = AHV
    tracer_wkc_parm%d_wkc          => WKC(1,1,1,1)
    tracer_wkc_parm%d_akt          => AKT(1,1,1,1,1)

    tracer_norm2_parm%cuBlock(1)     = imt
    tracer_norm2_parm%cuBlock(2)     = 1
    tracer_norm2_parm%cuBlock(3)     = 1
    tracer_norm2_parm%cuGrid(1)      = 1
    tracer_norm2_parm%cuGrid(2)      = jmt
    tracer_norm2_parm%cuGrid(3)      = 1
    tracer_norm2_parm%d_err_norm2    => d_err_norm2(1)
    tracer_norm2_parm%d_tarea        => tarea(1,1,1)
    tracer_norm2_parm%d_net          => net(1,1,1,1)
    tracer_norm2_parm%d_vit          => vit(1,1,1,1)

    tracer_norm3_parm%cuBlock(1)     = imt
    tracer_norm3_parm%cuBlock(2)     = 1
    tracer_norm3_parm%cuBlock(3)     = 1
    tracer_norm3_parm%cuGrid(1)      = 1
    tracer_norm3_parm%cuGrid(2)      = jmt
    tracer_norm3_parm%cuGrid(3)      = 1
    tracer_norm3_parm%d_tf           => tf(1,1,1,1)
    tracer_norm3_parm%d_vit          => vit(1,1,1,1)

    tracer_norm4_parm%cuBlock(1)     = imt
    tracer_norm4_parm%cuBlock(2)     = 1
    tracer_norm4_parm%cuBlock(3)     = 1
    tracer_norm4_parm%cuGrid(1)      = 1
    tracer_norm4_parm%cuGrid(2)      = jmt
    tracer_norm4_parm%cuGrid(3)      = 1
    tracer_norm4_parm%d_kmt          => kmt(1,1,1)
    tracer_norm4_parm%d_od0cp        = od0cp
    tracer_norm4_parm%d_gamma        = gamma
    tracer_norm4_parm%d_aidif        = 0.5D0
    tracer_norm4_parm%d_stf          => stf(1,1,1)
    tracer_norm4_parm%d_swv          => swv(1,1,1)
    tracer_norm4_parm%d_nswv         => nswv(1,1,1)
    tracer_norm4_parm%d_seaice       => seaice(1,1,1)
    tracer_norm4_parm%d_sst          => sst(1,1,1)
    tracer_norm4_parm%d_atb          => ATB(1,1,1,1,1)
    tracer_norm4_parm%d_odzp         => odzp(1)
    tracer_norm4_parm%d_tf           => tf(1,1,1,1)
    tracer_norm4_parm%d_net          => net(1,1,1,1)
#ifdef SST_restore
    tracer_norm4_parm%d_sst_restore  => sst_restore(1,1,1)
#endif
#ifdef COUP
    tracer_norm4_parm%d_tsf          => tsf(1,1,1)
#endif
end subroutine tracer_sw_init

end module sunway_init
