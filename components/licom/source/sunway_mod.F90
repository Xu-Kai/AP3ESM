!  CVS: $Id: param_mod.F90,v 1.1.1.1 2004/04/29 06:22:39 lhl Exp $
module sunway_mod
#include <def-undef.h>
use precision_mod
!-------------------------------------------------------------------------------
! Author: Maoxue YU  ( 12 July, 2024) !
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! POP_HaloMOD
!-------------------------------------------------------------------------------
type halo_2dr8_send_sw_param
      integer :: numMsgSend
      integer :: bufSizeSend
      integer :: rank
      integer :: i0
      integer ,pointer :: sizeSend
      integer ,pointer :: sendAddr
      real(r8),pointer :: d_array
      real(r8),pointer :: d_bufSendR8
end type

type halo_2dr8_recv_sw_param
      integer :: numMsgRecv
      integer :: bufSizeRecv
      integer :: nxGlobal
      integer :: rank
      integer :: i0
      integer ,pointer :: sizeRecv
      integer ,pointer :: recvAddr
      real(r8),pointer :: bufRecv2DR8
      real(r8),pointer :: bufTripole2DR8
      real(r8),pointer :: array
end type

type halo_2dr8_localcopy_sw_param
      integer :: numLocalCopies
      integer :: nxGlobal
      integer :: rank
      integer :: I0
      integer :: ioffset
      integer :: joffset
      integer :: isign
      integer ,pointer :: srcLocalAddr
      integer ,pointer :: dstLocalAddr
      real(r8)  :: fill
      real(r8),pointer :: bufTripole2D
      real(r8),pointer :: array
end type

type halo_3dr8_send_sw_param
      integer :: numMsgSend
      integer :: bufSizeSend
      integer :: strideArr
      integer :: strideSend
      integer :: I0
      integer :: I2
      integer ,pointer :: sendAddr
      integer ,pointer :: sizeSend
      real(r8),pointer :: bufSend3D
      real(r8),pointer :: array
end type


type halo_3dr8_recv_sw_param
      integer :: numMsgRecv
      integer :: bufSizeRecv
      integer :: strideArr
      integer :: strideBuf
      integer :: strideRecv
      integer :: nxGlobal
      integer :: I0
      integer :: I2
      integer ,pointer :: sizeRecv
      integer ,pointer :: recvAddr
      real(r8),pointer :: bufTripole3D
      real(r8),pointer :: bufRecv3D
      real(r8),pointer :: array
end type


type halo_3dr8_localcopy_sw_param
      integer :: strideArr
      integer :: strideBuf
      integer :: nxGlobal
      integer :: I0
      integer :: I2
      integer :: numLocalCopies
      integer :: ioffset
      integer :: joffset
      integer :: isign
      integer ,pointer :: srcLocalAddr
      integer ,pointer :: dstLocalAddr
      real(r8)  :: fill
      real(r8),pointer :: bufTripole3D
      real(r8),pointer :: array
end type


type(halo_2dr8_send_sw_param)  :: halo_2dr8_send_param
type(halo_2dr8_recv_sw_param)  :: halo_2dr8_recv_param
type(halo_2dr8_localcopy_sw_param)  :: halo_2dr8_localcopy_param
type(halo_3dr8_send_sw_param)  :: halo_3dr8_send_param
type(halo_3dr8_recv_sw_param)  :: halo_3dr8_recv_param
type(halo_3dr8_localcopy_sw_param)  :: halo_3dr8_localcopy_param


external :: slave_halo_2dr8_send_slave
external :: slave_halo_2dr8_recv_slave
external :: slave_halo_update2d_localcopy_slave_1
external :: slave_halo_update2d_localcopy_slave_2
external :: slave_halo_update3d_send_slave
external :: slave_halo_update3d_recv_slave
external :: slave_halo_update3d_localcopy_slave_1
external :: slave_halo_update3d_localcopy_slave_2


!-------------------------------------------------------------------------------
! readyt
!-------------------------------------------------------------------------------

    type readyt_cp
        real(r8),pointer :: d_dlu
        real(r8),pointer :: d_dlv
        real(r8),pointer :: d_rit
        real(r8),pointer :: d_ric
        real(r8),pointer :: d_rict
        real(r8),pointer :: d_ricdtTms
        real(r8),pointer :: d_ricdt
        real(r8),pointer :: d_akt
        integer :: rank
    end type

    type readyt_parm
        integer, dimension(3) :: cuBlock
        integer, dimension(3) :: cuGrid
        integer :: rank
        integer :: pad
         real(r8),pointer :: d_h0l
         real(r8),pointer :: d_h0f
         real(r8),pointer :: d_h0
         real(r8),pointer :: d_utl
         real(r8),pointer :: d_vtl
         real(r8),pointer :: d_utf
         real(r8),pointer :: d_vtf
         real(r8),pointer :: d_u
         real(r8),pointer :: d_v
         real(r8),pointer :: d_dlu
         real(r8),pointer :: d_dlv
         real(r8),pointer :: d_atb
         real(r8),pointer :: d_au0
         real(r8),pointer :: d_aus
         real(r8),pointer :: d_auw
         real(r8),pointer :: d_ausw
    end type

    type dens_parm
        integer, dimension(3) :: cuBlock
        integer, dimension(3) :: cuGrid
        integer :: d_km
        integer :: d_kmm1
        integer :: d_jmt
        integer :: d_imt
        integer :: rank
        integer :: pad
        real(r8) :: d_g
        real(r8) :: d_od0
        real(r8),pointer :: d_dlu
        real(r8),pointer :: d_dlv
        real(r8),pointer :: d_to
        real(r8),pointer :: d_so
        real(r8),pointer :: d_ric
        real(r8),pointer :: d_viv
        real(r8),pointer :: d_odzt
        real(r8),pointer :: d_c
        real(r8),pointer :: d_rict
        real(r8),pointer :: d_vit
        real(r8),pointer :: d_atb
    end type

    type readyt_2_parm
        integer, dimension(3) :: cuBlock
        integer, dimension(3) :: cuGrid
        integer :: d_km
        integer :: d_jmt
        integer :: d_imt
        integer :: rank
        real(r8) :: d_od0
        real(r8) :: d_g
        real(r8),pointer :: d_to
        real(r8),pointer :: d_so
        real(r8),pointer :: d_c
        real(r8),pointer :: d_po
        real(r8),pointer :: d_gg
        real(r8),pointer :: d_pdensity
        real(r8),pointer :: d_vit
        real(r8),pointer :: d_ppb
        real(r8),pointer :: d_ppc
        real(r8),pointer :: d_at
        real(r8),pointer :: d_dzp
    end type

    type readyt_2_2_parm
        integer, dimension(3) :: cuBlock
        integer, dimension(3) :: cuGrid
        integer :: d_km
        integer :: d_jmt
        integer :: d_imt
        integer :: rank
        real(r8) :: d_od0
        real(r8),pointer :: d_psa
        real(r8),pointer :: d_ppb
        real(r8),pointer :: d_ppc
        real(r8),pointer :: d_at
        real(r8),pointer :: d_gg
        real(r8),pointer :: d_vit
        real(r8),pointer :: d_pp
        real(r8),pointer :: d_ppa
        real(r8),pointer :: d_dzp
    end type

    type thermal_cu_parm
        integer, dimension(3) :: cuBlock
        integer, dimension(3) :: cuGrid
        integer :: d_km
        integer :: d_jmt
        integer :: d_imt
        integer :: rank
        real(r8) :: d_od0
        real(r8),pointer :: d_ppb
        real(r8),pointer :: d_ppc
        real(r8),pointer :: d_ppa
        real(r8),pointer :: d_alpha
        real(r8),pointer :: d_beta
        real(r8),pointer :: d_vit
    end type

    type readyt_cu_3_parm
        integer , dimension(3) :: cuBlock
        integer , dimension(3) :: cuGrid
        integer :: d_kmm1
        integer :: d_km
        integer :: d_jmt
        integer :: d_imt
        integer :: rank
        integer ,pointer :: d_kmu
        real(r8) :: d_od0cp
        real(r8) :: d_od0
        real(r8) :: d_g
        real(r8),pointer :: d_ricdttms
        real(r8),pointer :: d_ricdt
        real(r8),pointer :: d_gg
        real(r8),pointer :: d_pdensity
        real(r8),pointer :: d_vit
        real(r8),pointer :: d_pp
        real(r8),pointer :: d_at
        real(r8),pointer :: d_po
        real(r8),pointer :: d_dlu
        real(r8),pointer :: d_dlv
        real(r8),pointer :: d_dxur
        real(r8),pointer :: d_dyur
        real(r8),pointer :: d_buoytur
        real(r8),pointer :: d_buoysol
        real(r8),pointer :: d_nswv
        real(r8),pointer :: d_swv
        real(r8),pointer :: d_alpha
        real(r8),pointer :: d_beta
        real(r8),pointer :: d_odzt
        real(r8),pointer ::  d_akt
#if (defined CANUTOMIXOUT)
        real(r8),pointer :: d_alpha_canuto
        real(r8),pointer :: d_beta_canuto
#endif


    end type

    type readyt_cu_3_2_parm
        integer, dimension(3) :: cuBlock
        integer, dimension(3) :: cuGrid
        integer :: d_km
        integer :: d_jmt
        integer :: d_imt
        integer :: rank
        real(r8),pointer :: d_gg
        real(r8),pointer :: d_dzp
        real(r8),pointer :: d_dlu
        real(r8),pointer :: d_dlv
        real(r8),pointer :: d_ohbu
        real(r8),pointer :: d_ohbt
        real(r8),pointer :: d_zkt
        real(r8),pointer :: d_viv
        real(r8),pointer :: d_pxb
        real(r8),pointer :: d_pyb
    end type

    type readyt_cu_4_parm
        integer, dimension(3) :: cuBlock
        integer, dimension(3) :: cuGrid
        integer :: d_km
        integer :: d_jmt
        integer :: d_imt
        integer :: rank
        real(r8) :: d_g
        real(r8),pointer :: d_dlu
        real(r8),pointer :: d_dlv
        real(r8),pointer :: d_ohbt
    end type

    type readyt_cu_5_parm
        integer, dimension(3) :: cuBlock
        integer, dimension(3) :: cuGrid
        integer :: d_jmt
        integer :: d_imt
        integer :: rank
        integer ,pointer :: d_kmu
        real(r8) :: d_od0
        real(r8),pointer :: d_dlv
        real(r8),pointer :: d_au0
        real(r8),pointer :: d_aus
        real(r8),pointer :: d_auw
        real(r8),pointer :: d_ausw
        real(r8),pointer :: d_wgp
        real(r8),pointer :: d_work
        real(r8),pointer :: d_hbx
        real(r8),pointer :: d_hby
        real(r8),pointer :: d_viv
        real(r8),pointer :: d_whx
        real(r8),pointer :: d_why
        real(r8),pointer :: d_work1
        real(r8),pointer :: d_work2
        real(r8),pointer :: d_psa
        real(r8),pointer :: d_dxur
        real(r8),pointer :: d_dyur
        real(r8),pointer :: d_pax
        real(r8),pointer :: d_pay
    end type
   type pop_halo_parm
        integer :: data_len
        real(r8) :: fill
        real(r8),pointer :: d_data
   end type
   type(readyt_parm) :: mpe_param
   type(readyt_cp) :: readyt_cp_param
   type(pop_halo_parm) :: pop_halo_param
   type(dens_parm) :: dens_para
   type(readyt_2_parm) :: mpe_2_param
   type(readyt_2_2_parm) :: mpe_2_2_param
   type(thermal_cu_parm) :: thermal_param
   type(readyt_cu_3_parm) :: readyt_3_parm
   type(readyt_cu_3_2_parm) :: readyt_3_2_parm
   type(readyt_cu_4_parm) :: readyt_4_parm
   type(readyt_cu_5_parm) :: readyt_5_parm

   external :: slave_readyt_cu_1
   external :: slave_readyt_cp
   external :: slave_readyt_cp_1
   external :: slave_dens_cu
   external :: slave_readyt_cu_2
   external :: slave_readyt_cu_2_2
   external :: slave_thermal_cu
   external :: slave_readyt_cu_3
   external :: slave_readyt_cu_3_2
   external :: slave_readyt_cu_4
   external :: slave_readyt_cu_5
   external :: slave_pop_halo_cp

!-------------------------------------------------------------------------------
! readyc
!-------------------------------------------------------------------------------
 type readyc_vinteg_cu_parm
         integer, dimension(3) :: cuGrid
         integer, dimension(3) :: cuBlock
         real(r8),pointer :: d_dlu
         real(r8),pointer :: d_dlv
         real(r8),pointer :: d_dlub
         real(r8),pointer :: d_dlvb
         real(r8),pointer :: d_dzp
         real(r8),pointer :: d_ohbu
         real(r8),pointer :: d_viv
    end type

    type readyc1_cu_parm
         integer, dimension(3) :: cuGrid
         integer, dimension(3) :: cuBlock
         real(r8),pointer :: d_h0bl
         real(r8),pointer :: d_h0bf
         real(r8),pointer :: d_h0
         real(r8),pointer :: d_up
         real(r8),pointer :: d_vp
         real(r8),pointer :: d_viv
         real(r8),pointer :: d_vit
         real(r8),pointer :: d_wp12
         real(r8),pointer :: d_wp13
         real(r8),pointer :: d_at0
         real(r8),pointer :: d_atn
         real(r8),pointer :: d_ate
         real(r8),pointer :: d_atne
    end type

    type readyc1_2_cu_parm
         integer, dimension(3) :: cuGrid
         integer, dimension(3) :: cuBlock
         real(r8),pointer :: d_wp12
         real(r8),pointer :: d_wp13
         real(r8),pointer :: d_vit
         real(r8),pointer :: d_s2t
         real(r8),pointer :: d_riv1
         real(r8),pointer :: d_riv2
         real(r8),pointer :: d_odzt
         real(r8),pointer :: d_rit
         real(r8),pointer :: d_rict
         real(r8),pointer :: d_up
         real(r8),pointer :: d_vp
         real(r8),pointer :: d_s2u
         real(r8),pointer :: d_viv
         real(r8),pointer :: d_riu
         real(r8),pointer :: d_ric
    end type
#ifdef BCKMEX
    type readyc1_3_cu_parm
         integer, dimension(3) :: cuGrid
         integer, dimension(3) :: cuBlock
         real(r8),pointer :: d_diff_back
         real(r8),pointer :: d_diff_back_sh
         real(r8),pointer :: d_diff_back_nh
         real(r8),pointer :: d_tlat
         real(r8)         :: d_diff_back_coef_max
         real(r8)         :: d_diff_back_eq
         real(r8)         :: d_diff_back_coef
    endtype
#endif
    type readyc_turbo_cu_parm
         integer, dimension(3) :: cuGrid
         integer, dimension(3) :: cuBlock
         integer :: d_km
         integer :: d_imt
         integer :: d_jmt
         integer :: rank
         integer,pointer ::   d_irimax
         integer,pointer ::   d_kmt
         real(r8),pointer ::  d_akmt
         real(r8),pointer ::  d_akt
         real(r8),pointer ::  d_ak_tide
         real(r8),pointer ::  d_amld
         real(r8),pointer ::  d_amtaun2a1
         real(r8),pointer ::  d_and2on2a1
         real(r8),pointer ::  d_at
         real(r8),pointer ::  d_back_ra_r
         real(r8),pointer ::  d_buoysol
         real(r8),pointer ::  d_buoytur
         real(r8),pointer ::  d_dzp
         real(r8),pointer ::  d_fcort
         real(r8),pointer ::  d_fztidal
         real(r8),pointer ::  d_fz_tide
         real(r8),pointer ::  d_pdensity
         real(r8),pointer ::  d_rib
         real(r8),pointer ::  d_ridb
         real(r8),pointer ::  d_ricdttms
         real(r8),pointer ::  d_richardson
         real(r8),pointer ::  d_rict
         real(r8),pointer ::  d_rit
         real(r8),pointer ::  d_s2t
         real(r8),pointer ::  d_sha1
         real(r8),pointer ::  d_shb
         real(r8),pointer ::  d_sh_r1
         real(r8),pointer ::  d_slq2b
         real(r8),pointer ::  d_slq2_r1
         real(r8),pointer ::  d_sma1
         real(r8),pointer ::  d_smb
         real(r8),pointer ::  d_sm_r1
         real(r8),pointer ::  d_ssa1
         real(r8),pointer ::  d_ssb
         real(r8),pointer ::  d_ss_r1
         real(r8),pointer ::  d_turb_param
         real(r8),pointer ::  d_vit
         real(r8),pointer ::  d_wave_dis
         real(r8),pointer ::  d_wp3_tidal
         real(r8),pointer ::  d_zkp
#ifdef BCKMEX
         real(r8),pointer ::  d_diff_back
#endif
         !real(r8),pointer ::  d_wk1
         !real(r8),pointer ::  d_wk2
         !real(r8),pointer ::  d_wk3
    end type

    type readyc_canuto2010_parm
        integer, dimension(3) :: cuGrid
        integer, dimension(3) :: cuBlock
        integer,pointer       :: d_mldlev
        integer,pointer       :: d_kmt
        real(r8)              :: d_pi0_1
        real(r8),pointer      :: d_ak_tide_mixing
        real(r8),pointer      :: d_dzp
        real(r8),pointer      :: d_zkp
        real(r8),pointer      :: d_ulat
        real(r8),pointer      :: d_amld
        real(r8),pointer      :: d_rit
        real(r8),pointer      :: d_rict
        real(r8),pointer      :: d_ricdt
        real(r8),pointer      :: d_s2t
        real(r8),pointer      :: d_vit
        real(r8),pointer      :: d_akt
        real(r8),pointer      :: d_at
        real(r8),pointer      :: d_to
        real(r8),pointer      :: d_so
        real(r8),pointer      :: d_c
        real(r8),pointer      :: d_wave_dis
        real(r8),pointer      :: d_fz_tide
        real(r8),pointer      :: d_akmt
        !real(r8),pointer      :: d_wk1
        !real(r8),pointer      :: d_wk2
        !real(r8),pointer      :: d_wk3
        !real(r8),pointer      :: d_wp3
        integer               :: rank
    end type
    type readyc2_cu_parm
         integer, dimension(3) :: cuGrid
         integer, dimension(3) :: cuBlock
         real(r8),pointer :: d_akmu
         real(r8),pointer :: d_au0
         real(r8),pointer :: d_akmt
         real(r8),pointer :: d_aus
         real(r8),pointer :: d_auw
         real(r8),pointer :: d_ausw
         real(r8),pointer :: d_viv
    end type

    type readyc_upwell_0_parm
         integer, dimension(3) :: cuGrid
         integer, dimension(3) :: cuBlock
         integer          :: rank
         integer          :: pad
         real(r8),pointer :: d_au0
         real(r8),pointer :: d_aus
         real(r8),pointer :: d_auw
         real(r8),pointer :: d_ausw
         real(r8),pointer :: d_h0
         real(r8),pointer :: d_ohbu
         real(r8),pointer :: d_u
         real(r8),pointer :: d_uk
         real(r8),pointer :: d_v
         real(r8),pointer :: d_vk
         real(r8),pointer :: d_work
    end type

    type readyc_upwell_parm
         integer, dimension(3) :: cuGrid
         integer, dimension(3) :: cuBlock
         integer          :: rank
         integer,pointer  :: d_kmt
         real(r8),pointer :: d_dzp
         real(r8),pointer :: d_h0
         real(r8),pointer :: d_hts
         real(r8),pointer :: d_htw
         real(r8),pointer :: d_ohbt
         real(r8),pointer :: d_tarea_r
         real(r8),pointer :: d_uk
         real(r8),pointer :: d_vit
         real(r8),pointer :: d_vk
         real(r8),pointer :: d_wka
         real(r8),pointer :: d_work
         real(r8),pointer :: d_ws
    end type
    type readyc_upwell_1_parm
         integer, dimension(3) :: cuGrid
         integer, dimension(3) :: cuBlock
         integer          :: rank
         integer          :: padd
         real(r8),pointer :: d_vit
         real(r8),pointer :: d_dzp
         real(r8),pointer :: d_wka
         real(r8),pointer :: d_work
    end type
    type readyc_upwell_2_parm
         integer, dimension(3) :: cuGrid
         integer, dimension(3) :: cuBlock
         integer          :: rank
         integer          :: k
         real(r8),pointer :: d_ws
         real(r8),pointer :: d_vit
         real(r8),pointer :: d_h0
         real(r8),pointer :: d_dzp
         real(r8),pointer :: d_ohbt
         real(r8),pointer :: d_wka
         real(r8),pointer :: d_work
         real(r8),pointer :: d_work1
         real(r8),pointer :: d_wsbak
    end type

    type readyc3_parm
         integer, dimension(3) :: cuGrid
         integer, dimension(3) :: cuBlock
         real(r8),pointer :: d_wka
         real(r8),pointer :: d_ws
         real(r8),pointer :: d_au0
         real(r8),pointer :: d_ausw
         real(r8),pointer :: d_aus
         real(r8),pointer :: d_auw
         real(r8),pointer :: d_u_wface
         real(r8),pointer :: d_v_sface
         real(r8),pointer :: d_u
         real(r8),pointer :: d_v
         real(r8),pointer :: d_hue
         real(r8),pointer :: d_hun
    end type

    type advection_momentum_parm
         integer, dimension(3) :: cuGrid
         integer, dimension(3) :: cuBlock
         real(r8),pointer ::d_dlu
         real(r8),pointer ::d_dlv
         real(r8),pointer ::d_odzp
         real(r8),pointer ::d_u
         real(r8),pointer ::d_uarea_r
         real(r8),pointer ::d_u_wface
         real(r8),pointer ::d_v
         real(r8),pointer ::d_v_sface
         real(r8),pointer ::d_wka
!         real(r8),pointer ::d_utest
    end type

    type readyc4_parm
         integer, dimension(3) :: cuGrid
         integer, dimension(3) :: cuBlock
         integer, pointer ::d_kmu
         real(r8),pointer ::d_sbcx
         real(r8),pointer ::d_sbcy
         real(r8),pointer ::d_su
         real(r8),pointer ::d_sv
         real(r8),pointer ::d_bbcx
         real(r8)         ::d_od0
         real(r8)         ::d_c0f
         real(r8),pointer ::d_up
         real(r8),pointer ::d_vp
         real(r8),pointer ::d_snlat
         real(r8)         ::d_cag
         real(r8)         ::d_sag
         real(r8),pointer ::d_dlub
         real(r8),pointer ::d_dlvb
         real(r8),pointer ::d_dzp
         real(r8),pointer ::d_ohbu
         real(r8),pointer ::d_dlu
         real(r8),pointer ::d_viv
         real(r8),pointer ::d_dlv
         real(r8),pointer ::d_bbcy
    end type

    type readyc4_2_parm
         integer, dimension(3) :: cuGrid
         integer, dimension(3) :: cuBlock
         real(r8)         ::d_od0
         real(r8)         ::d_sag
         real(r8)         ::d_cag
         real(r8),pointer ::d_su
         real(r8),pointer ::d_sv
         real(r8),pointer ::d_up
         real(r8),pointer ::d_vp
         real(r8),pointer ::d_akmu
         real(r8),pointer ::d_odzt
         real(r8),pointer ::d_viv
         real(r8),pointer ::d_wka
         real(r8),pointer ::d_snlat
         real(r8),pointer ::d_dlv
         real(r8),pointer ::d_odzp
         real(r8),pointer ::d_dlu
    end type

    type hdel4_div_cu_param
        integer, dimension(3) :: cuBlock
        integer, dimension(3) :: cuGrid
        integer :: rank
        integer ,pointer :: d_kmt
        real(r8),pointer :: d_hts
        real(r8),pointer :: d_htw
        real(r8),pointer :: d_tarea_r
        real(r8),pointer :: d_up
        real(r8),pointer :: d_vp
        real(r8),pointer :: d_wka
    end type

    type hdel4_grad_cu_param
        integer, dimension(3) :: cuBlock
        integer, dimension(3) :: cuGrid
        integer :: rank
        integer ,pointer :: d_kmu
        real(r8),pointer :: d_wka
        real(r8),pointer :: d_dxur
        real(r8),pointer :: d_dyur
        real(r8),pointer :: d_am_factor
        real(r8),pointer :: d_wkz
    end type

    type hdel4_zcurl_cu_param
        integer, dimension(3) :: cuBlock
        integer, dimension(3) :: cuGrid
        integer :: rank
        integer ,pointer :: d_kmt
        real(r8) ,pointer :: d_dxu
        real(r8),pointer :: d_dyu
        real(r8),pointer :: d_up
        real(r8),pointer :: d_vp
        real(r8),pointer :: d_wkz
        real(r8),pointer :: d_tarea_r
    end type

    type hdiffu_del4_0_cu_param
        integer, dimension(3) :: cuBlock
        integer, dimension(3) :: cuGrid
        integer :: d_jb
        integer :: d_je
        integer :: d_ib
        integer :: d_ie
        integer :: rank
        integer ,pointer :: d_kmu
        real(r8) :: d_am
        real(r8),pointer :: d_uarea_r
        real(r8),pointer :: d_up
        real(r8),pointer :: d_vp
        real(r8),pointer :: d_duc
        real(r8),pointer :: d_dum
        real(r8),pointer :: d_dun
        real(r8),pointer :: d_dus
        real(r8),pointer :: d_due
        real(r8),pointer :: d_duw
        real(r8),pointer :: d_dmc
        real(r8),pointer :: d_dmn
        real(r8),pointer :: d_dms
        real(r8),pointer :: d_dme
        real(r8),pointer :: d_dmw
        real(r8),pointer :: d_amf
        real(r8),pointer :: d_am_factor
        real(r8),pointer :: d_d2uk
        real(r8),pointer :: d_d2vk
    end type

   type hdiffu_del4_1_cu_param
        integer, dimension(3) :: cuBlock
        integer, dimension(3) :: cuGrid
        integer :: d_jb
        integer :: d_je
        integer :: d_ib
        integer :: d_ie
        integer :: rank
        integer ,pointer :: d_kmu
        real(r8) :: d_am
        real(r8) :: d_c0f
        real(r8),pointer :: d_wka
        real(r8),pointer :: d_up
        real(r8),pointer :: d_vp
        real(r8),pointer :: d_duc
        real(r8),pointer :: d_dum
        real(r8),pointer :: d_dun
        real(r8),pointer :: d_dus
        real(r8),pointer :: d_due
        real(r8),pointer :: d_duw
        real(r8),pointer :: d_dmc
        real(r8),pointer :: d_dmn
        real(r8),pointer :: d_dms
        real(r8),pointer :: d_dme
        real(r8),pointer :: d_dmw
        real(r8),pointer :: d_dlu
        real(r8),pointer :: d_dlv
        real(r8),pointer :: d_d2uk
        real(r8),pointer :: d_d2vk
        real(r8),pointer :: d_am_factor
        real(r8),pointer :: d_uarea_r
        real(r8),pointer :: d_fac_vis_uv
        !real(r8),pointer :: d_hdukkm
        !real(r8),pointer :: d_hdvkkm
        !real(r8),pointer :: d_am_f_km
        real(r8) :: d_dts
        real(r8) :: d_value
        real(r8) :: d_value2
    end type

   type hdiffu_del2_cu_param
        integer, dimension(3) :: cuBlock
        integer, dimension(3) :: cuGrid
        integer :: d_jb
        integer :: d_je
        integer :: d_ib
        integer :: d_ie
        integer :: rank
        integer ,pointer :: d_kmu
        real(r8) :: d_am
        real(r8) :: d_c0f
        real(r8),pointer :: d_wka
        real(r8),pointer :: d_up
        real(r8),pointer :: d_vp
        real(r8),pointer :: d_duc
        real(r8),pointer :: d_dum
        real(r8),pointer :: d_dun
        real(r8),pointer :: d_dus
        real(r8),pointer :: d_due
        real(r8),pointer :: d_duw
        real(r8),pointer :: d_dmc
        real(r8),pointer :: d_dmn
        real(r8),pointer :: d_dms
        real(r8),pointer :: d_dme
        real(r8),pointer :: d_dmw
        real(r8),pointer :: d_dlu
        real(r8),pointer :: d_dlv
        real(r8),pointer :: d_viv
    end type


    type(readyc1_cu_parm) ::readyc1_param
    type(readyc1_2_cu_parm) ::readyc1_2_param
#ifdef BCKMEX
    type(readyc1_3_cu_parm) ::readyc1_3_param
#endif
    type(readyc_turbo_cu_parm) ::readyc_turbo_parm
    type(readyc_canuto2010_parm) ::canuto2010_parm
    type(readyc2_cu_parm) :: readyc2_parm
    type(readyc_upwell_0_parm) :: readyc_upwell_0_param
    type(readyc_upwell_parm) :: readyc_upwell_param
    type(readyc_upwell_1_parm) :: readyc_upwell_1_param
    type(readyc_upwell_2_parm) :: readyc_upwell_2_param
    type(readyc3_parm) :: readyc3_param
    type(advection_momentum_parm) :: advection_momentum_param
    type(readyc4_parm) :: readyc4_param
    type(readyc4_2_parm) ::readyc4_2_param
    type(readyc_vinteg_cu_parm) ::readyc_vinteg_param
    type(hdel4_div_cu_param) :: hdel4_div_param
    type(hdel4_grad_cu_param) :: hdel4_grad_param
    type(hdel4_zcurl_cu_param) :: hdel4_zcurl_param
    type(hdiffu_del4_0_cu_param) :: hdiffu_del4_0_param
    type(hdiffu_del4_1_cu_param) :: hdiffu_del4_1_param
    type(hdiffu_del2_cu_param) ::hdiffu_del2_param

    external :: slave_readyc1_cu,slave_readyc1_2, slave_readyc1_3, slave_turb_cu, slave_readyc2, slave_upwell_cu_0, slave_upwell_cu, slave_readyc3, slave_advection_momentum_cu, slave_hdel4_div,slave_readyc4,slave_readyc4_2,slave_readyc_vinteg_cu
    external :: slave_upwell_cu_1
    external :: slave_upwell_cu_2
    external :: slave_hdel4_grad
    external :: slave_hdel4_zcurl
    external :: slave_hdiffu_del4_0
    external :: slave_hdiffu_del4_1
    external :: slave_hdiffu_del2
    external :: slave_readyc_canuto2010

!-------------------------------------------------------------------------------
! barotr
!-------------------------------------------------------------------------------

      type barotr1_cu_param
            integer, dimension(3) :: cuBlock
            integer, dimension(3) :: cuGrid
            integer :: d_km
            integer :: d_imt
            integer :: d_jmt
            integer :: rank
            real(r8),pointer :: d_work
            real(r8),pointer :: d_wka
            real(r8),pointer :: d_dzph
            real(r8),pointer :: d_h0
            real(r8),pointer :: d_ub
            real(r8),pointer :: d_vb
            real(r8),pointer :: d_au0
            real(r8),pointer :: d_ausw
            real(r8),pointer :: d_aus
            real(r8),pointer :: d_auw
      end type

      type barotr1_cu1_param
            integer, dimension(3) :: cuBlock
            integer, dimension(3) :: cuGrid
            integer :: d_imt
            integer :: d_jmt
            integer :: rank
            integer ,pointer :: d_kmt
            real(r8),pointer :: d_work
            real(r8),pointer :: d_wka
            real(r8),pointer :: d_vit
            real(r8),pointer :: d_htw
            real(r8),pointer :: d_hts
            real(r8),pointer :: d_tarea_r
      end type

      type barotr2_cu1_param
            integer, dimension(3) :: cuBlock
            integer, dimension(3) :: cuGrid
            integer :: d_imt
            integer :: d_jmt
            real(r8) :: d_dtb
            real(r8),pointer :: d_h0
            real(r8),pointer :: d_h0p
            real(r8),pointer :: d_work
#ifdef TIDEFORCE
            real(r8) :: BETA_TIDE
            real(r8),pointer :: d_h0_tide
            real(r8),pointer :: d_elev            
#endif
      end type

      type bhdel4_div_cu_param
            integer, dimension(3) :: cuBlock
            integer, dimension(3) :: cuGrid
            integer :: rank
            integer ,pointer :: d_kmt
            real(r8),pointer :: d_hts
            real(r8),pointer :: d_htw
            real(r8),pointer :: d_tarea_r
            real(r8),pointer :: d_ubp
            real(r8),pointer :: d_vbp
            real(r8),pointer :: d_wka
            real(r8),pointer :: d_diverg
      end type

      type bhdel4_grad_cu_param
            integer, dimension(3) :: cuBlock
            integer, dimension(3) :: cuGrid
            integer :: rank
            integer ,pointer :: d_kmu
            real(r8),pointer :: d_wka
            real(r8),pointer :: d_dxur
            real(r8),pointer :: d_dyur
            real(r8),pointer :: d_am_factor
            real(r8),pointer :: d_z_curl
            real(r8),pointer :: d_diverg
      end type

      type bhdel4_zcurl_cu_param
            integer, dimension(3) :: cuBlock
            integer, dimension(3) :: cuGrid
            integer :: rank
            integer ,pointer :: d_kmt
            real(r8) ,pointer :: d_dxu
            real(r8),pointer :: d_dyu
            real(r8),pointer :: d_ubp
            real(r8),pointer :: d_vbp
            real(r8),pointer :: d_wka
            real(r8),pointer :: d_tarea_r
            real(r8),pointer :: d_z_curl
      end type


      type barotr2_udel4_0_cu_param
            integer, dimension(3) :: cuBlock
            integer, dimension(3) :: cuGrid
            integer :: d_jb
            integer :: d_je
            integer :: d_ib
            integer :: d_ie
            integer :: d_imt
            integer :: d_jmt
            integer :: rank
            integer ,pointer :: d_kmu
            integer, dimension(2) :: pad
            real(r8) :: d_am
            real(r8),pointer :: d_ubp
            real(r8),pointer :: d_vbp
            real(r8),pointer :: d_uarea_r
            real(r8),pointer :: d_duc
            real(r8),pointer :: d_dum
            real(r8),pointer :: d_dun
            real(r8),pointer :: d_dus
            real(r8),pointer :: d_due
            real(r8),pointer :: d_duw
            real(r8),pointer :: d_dmc
            real(r8),pointer :: d_dmn
            real(r8),pointer :: d_dms
            real(r8),pointer :: d_dme
            real(r8),pointer :: d_dmw
            real(r8),pointer :: d_amf
            real(r8),pointer :: d_am_factor
            real(r8),pointer :: d_d2uk
            real(r8),pointer :: d_d2vk
      end type

      type barotr2_udel4_1_cu_param
            integer, dimension(3) :: cuBlock
            integer, dimension(3) :: cuGrid
            integer :: d_jb
            integer :: d_je
            integer :: d_ib
            integer :: d_ie
            integer :: d_imt
            integer :: d_jmt
            integer :: d_km
            integer :: rank
            integer ,pointer :: d_kmu
            integer :: pad
            real(r8) :: d_am
            real(r8),pointer :: d_viv
            real(r8),pointer :: d_ebea
            real(r8),pointer :: d_ebeb
            real(r8),pointer :: d_fcor
#ifdef TIDEFORCE
            real(r8),pointer :: d_h0_tide
#endif 
            real(r8),pointer :: d_h0
            real(r8),pointer :: d_ubp
            real(r8),pointer :: d_vbp
            real(r8),pointer :: d_dlub
            real(r8),pointer :: d_dlvb
            real(r8),pointer :: d_wka
            real(r8),pointer :: d_work
            real(r8),pointer :: d_wgp
            real(r8),pointer :: d_pax
            real(r8),pointer :: d_pay
            real(r8),pointer :: d_pxb
            real(r8),pointer :: d_pyb
            real(r8),pointer :: d_whx
            real(r8),pointer :: d_why
            real(r8),pointer :: d_duc
            real(r8),pointer :: d_dum
            real(r8),pointer :: d_dun
            real(r8),pointer :: d_dus
            real(r8),pointer :: d_due
            real(r8),pointer :: d_duw
            real(r8),pointer :: d_dmc
            real(r8),pointer :: d_dmn
            real(r8),pointer :: d_dms
            real(r8),pointer :: d_dme
            real(r8),pointer :: d_dmw
            real(r8),pointer :: d_dxur
            real(r8),pointer :: d_dyur
            real(r8),pointer :: d_au0
            real(r8),pointer :: d_ausw
            real(r8),pointer :: d_aus
            real(r8),pointer :: d_auw
            real(r8),pointer :: d_d2uk
            real(r8),pointer :: d_d2vk
            real(r8),pointer :: d_am_factor
            real(r8),pointer :: d_uarea_r
            real(r8),pointer :: d_factor
            real(r8),pointer :: d_fac_vis_uv
            real(r8) :: d_dts
            real(r8) :: d_value
            real(r8) :: d_value2
      end type

      type barotr2_cu_param
            integer, dimension(3) :: cuBlock
            integer, dimension(3) :: cuGrid
            integer :: d_jb
            integer :: d_je
            integer :: d_ib
            integer :: d_ie
            integer :: d_imt
            integer :: d_jmt
            integer :: d_km
            integer :: rank
            integer ,pointer :: d_kmu
            integer :: pad
            real(r8) :: d_am
            real(r8),pointer :: d_viv
            real(r8),pointer :: d_ebea
            real(r8),pointer :: d_ebeb
            real(r8),pointer :: d_fcor
#ifdef TIDEFORCE
            real(r8),pointer :: d_h0_tide
#endif
            real(r8),pointer :: d_h0
            real(r8),pointer :: d_ubp
            real(r8),pointer :: d_vbp
            real(r8),pointer :: d_dlub
            real(r8),pointer :: d_dlvb
            real(r8),pointer :: d_wka
            real(r8),pointer :: d_work
            real(r8),pointer :: d_wgp
            real(r8),pointer :: d_pax
            real(r8),pointer :: d_pay
            real(r8),pointer :: d_pxb
            real(r8),pointer :: d_pyb
            real(r8),pointer :: d_whx
            real(r8),pointer :: d_why
            real(r8),pointer :: d_duc
            real(r8),pointer :: d_dum
            real(r8),pointer :: d_dun
            real(r8),pointer :: d_dus
            real(r8),pointer :: d_due
            real(r8),pointer :: d_duw
            real(r8),pointer :: d_dmc
            real(r8),pointer :: d_dmn
            real(r8),pointer :: d_dms
            real(r8),pointer :: d_dme
            real(r8),pointer :: d_dmw
            real(r8),pointer :: d_dxur
            real(r8),pointer :: d_dyur
            real(r8),pointer :: d_au0
            real(r8),pointer :: d_ausw
            real(r8),pointer :: d_aus
            real(r8),pointer :: d_auw
      end type

      type barotr2_tide_cu_param
            integer, dimension(3) :: cuBlock
            integer, dimension(3) :: cuGrid

            real(r8) :: d_d0
            real(r8),pointer :: d_wka
            real(r8),pointer :: d_ubp
            real(r8),pointer :: d_vbp
            real(r8),pointer :: d_ohbu
            real(r8),pointer :: d_h0
            real(r8),pointer :: d_cd_tide 
            real(r8),pointer :: d_viv 
            integer ,pointer :: d_kmu
            real(r8),pointer :: d_ric
            real(r8),pointer :: d_cd_inwave
            real(r8),pointer :: d_depthstd
      end type

      type barotr3_cu_param
            integer, dimension(3) :: cuBlock
            integer, dimension(3) :: cuGrid
            integer :: d_imt
            integer :: d_jmt
            integer :: rank
            integer :: pad
            real(r8) :: d_dtb
            real(r8) ,pointer :: d_dzph
            real(r8),pointer :: d_ub
            real(r8),pointer :: d_vb
            real(r8),pointer :: d_ubp
            real(r8),pointer :: d_vbp
            real(r8),pointer :: d_wka
            real(r8),pointer :: d_work
      end type

      type barotr3_tdel4_1_cu_param
            integer, dimension(3) :: cuBlock
            integer, dimension(3) :: cuGrid
            integer :: d_jb
            integer :: d_je
            integer :: d_ib
            integer :: d_ie
            integer :: d_nc
            integer :: d_imt
            integer :: d_jmt
            integer :: d_km
            integer :: rank
            integer ,pointer :: d_kmt
            integer ,pointer :: d_kmtn
            integer ,pointer :: d_kmts
            integer ,pointer :: d_kmte
            integer ,pointer :: d_kmtw
            real(r8),pointer :: d_h0p
            real(r8),pointer :: d_dtn
            real(r8),pointer :: d_dts
            real(r8),pointer :: d_dte
            real(r8),pointer :: d_dtw
            real(r8),pointer :: d_d2tk
            real(r8),pointer :: d_ahf
      end type

      type barotr3_cu1_param
            integer, dimension(3) :: cuBlock
            integer, dimension(3) :: cuGrid
            integer :: d_jb
            integer :: d_je
            integer :: d_ib
            integer :: d_ie
            integer :: d_nc
            integer :: d_imt
            integer :: d_jmt
            integer :: d_km
            integer :: rank
            integer :: padd
            real(r8) :: d_ah
            integer ,pointer :: d_kmt
            integer ,pointer :: d_kmtn
            integer ,pointer :: d_kmts
            integer ,pointer :: d_kmte
            integer ,pointer :: d_kmtw
            real(r8),pointer :: d_vit
            real(r8),pointer :: d_h0p
            real(r8),pointer :: d_wka
            real(r8),pointer :: d_work
            real(r8),pointer :: d_htw
            real(r8),pointer :: d_hts
            real(r8),pointer :: d_tarea_r
            real(r8),pointer :: d_dtn
            real(r8),pointer :: d_dts
            real(r8),pointer :: d_dte
            real(r8),pointer :: d_dtw
            real(r8),pointer :: d_d2tk
            real(r8),pointer :: d_factor_vis_tsz0
      end type

      type barotr4_cu_param
            integer, dimension(3) :: cuBlock
            integer, dimension(3) :: cuGrid
            integer :: d_imt
            integer :: d_jmt
            logical :: d_iau
            real(r8) :: d_dtb
            real(r8) :: d_lamda
            real(r8) ,pointer :: d_h0
            real(r8) ,pointer :: d_h0p
            real(r8),pointer :: d_ub
            real(r8),pointer :: d_vb
            real(r8),pointer :: d_ubp
            real(r8),pointer :: d_vbp
            real(r8),pointer :: d_h0f
            real(r8),pointer :: d_h0bf
            real(r8),pointer :: d_work
            real(r8),pointer :: d_h0_iau
            real(r8),pointer :: d_vit
      end type

      type(barotr1_cu_param)  :: barotr1_param
      type(barotr1_cu1_param) :: barotr1_1_param
      type(barotr2_cu1_param) :: barotr2_1_param
      type(bhdel4_div_cu_param) :: bhdel4_div_param
      type(bhdel4_grad_cu_param) :: bhdel4_grad_param
      type(bhdel4_zcurl_cu_param) :: bhdel4_zcurl_param
      type(barotr2_udel4_0_cu_param) :: barotr2_udel4_0_param
      type(barotr2_udel4_1_cu_param) :: barotr2_udel4_1_param
      type(barotr2_cu_param) :: barotr2_param
      type(barotr2_tide_cu_param) :: barotr2_tide_param
      type(barotr3_cu_param) :: barotr3_param
      type(barotr3_tdel4_1_cu_param) :: barotr3_tdel4_1_param
      type(barotr3_cu1_param) :: barotr3_1_param
      type(barotr4_cu_param) :: barotr4_param

      external :: slave_barotr1_cu
      external :: slave_barotr1_cu1
      external :: slave_barotr2_cu1
      external :: slave_bhdel4_div
      external :: slave_bhdel4_grad
      external :: slave_bhdel4_zcurl
      external :: slave_barotr2_udel4_0
      external :: slave_barotr2_udel4_1
      external :: slave_barotr2_cu
      external :: slave_barotr2_tide
      external :: slave_barotr2_tideGW
      external :: slave_barotr3_cu
      external :: slave_barotr3_tdel4_1
      external :: slave_barotr3_cu1
      external :: slave_barotr4_cu

!-------------------------------------------------------------------------------
! bclinc
!-------------------------------------------------------------------------------


      type bclinc12_cu_slave_param
            integer, dimension(3) :: cuBlock
            integer, dimension(3) :: cuGrid
            real(r8):: d_od0
            real(r8):: d_c0f
            real(r8):: d_cag
            real(r8):: d_sag
            real(r8),pointer :: d_snlat
            real(r8),pointer :: d_up
            real(r8),pointer :: d_vp
            real(r8),pointer :: d_bbcy
            real(r8),pointer :: d_sbcy
            real(r8),pointer :: d_sbcx
            real(r8),pointer :: d_bbcx
            real(r8),pointer :: d_su
            real(r8),pointer :: d_sv
            integer,pointer  :: d_kmu
      end type

      type bclinc11_cu_slave_param
            integer, dimension(3) :: cuBlock
            integer, dimension(3) :: cuGrid
            integer :: d_isc
            integer :: rank
            real(r8):: d_od0
            real(r8):: d_onbb
            real(r8),pointer :: d_ohbt
            real(r8),pointer :: d_vit
            real(r8),pointer :: d_dzp
            real(r8),pointer :: d_h0bf
            real(r8),pointer :: d_h0bl
            real(r8),pointer :: d_gg
            real(r8),pointer :: d_psa
            real(r8),pointer :: d_wka
            real(r8),pointer :: d_wka1
            real(r8),pointer :: d_wka2
            real(r8),pointer :: d_zkt
            real(r8),pointer :: d_work
      end type

      type bclinc2_cu_slave_param
            integer, dimension(3) :: cuBlock
            integer, dimension(3) :: cuGrid
            integer :: rank
            integer :: pad
            real(r8) ::d_dtc2
            real(r8),pointer ::d_vp
            real(r8),pointer ::d_dlv
            real(r8),pointer ::d_wka
      end type

      type bclinc3_cu_slave_param
            integer, dimension(3) :: cuBlock
            integer, dimension(3) :: cuGrid
            real(r8) ::d_dtc2
            integer,pointer :: d_kmu
            real(r8),pointer ::d_odzt
            real(r8),pointer ::d_odzp
            real(r8),pointer ::d_viv
            real(r8),pointer ::d_wk
            real(r8),pointer ::d_topbc
            real(r8),pointer ::d_bombc
            real(r8),pointer ::d_dcb
      end type

    type bclinc_invtriu_cu_slave_param
        integer, dimension(3) :: cuBlock
        integer, dimension(3) :: cuGrid
        real(r8) ::d_dtc2
        integer,pointer :: d_kmu
        real(r8),pointer ::d_akmu
        real(r8),pointer ::d_odzt
        real(r8),pointer ::d_odzp
        real(r8),pointer ::d_u
        real(r8),pointer ::d_v
        real(r8),pointer ::d_sbcx
        real(r8),pointer ::d_sbcy
        real(r8),pointer ::d_bbcx
        real(r8),pointer ::d_bbcy
        real(r8),pointer ::d_viv
        !real(r8),pointer ::d_wk1
        !real(r8),pointer ::d_wk2
    end type
      type bclinc4_0_cu_slave_param
            integer, dimension(3) :: cuBlock
            integer, dimension(3) :: cuGrid
            real(r8) ::d_afc1
            real(r8) ::d_afc2
            real(r8) ::d_dtc2
            real(r8),pointer ::d_viv
            real(r8),pointer ::d_v
            real(r8),pointer ::d_up
            real(r8),pointer ::d_vp
            real(r8),pointer ::d_vb
            real(r8),pointer ::d_dlu
            real(r8),pointer ::d_wka
            real(r8),pointer ::d_work
      end type

      type bclinc4_cu_slave_param
            integer, dimension(3) :: cuBlock
            integer, dimension(3) :: cuGrid
            integer :: rank
            integer,pointer :: d_kmu
            real(r8) ::d_afc1
            real(r8) ::d_afc2
            real(r8) ::d_dtc2
            real(r8),pointer ::d_ohbu
            real(r8),pointer ::d_viv
            real(r8),pointer ::d_dzp
            real(r8),pointer ::d_odzt
            real(r8),pointer ::d_odzp
            real(r8),pointer ::d_akmu
            real(r8),pointer ::d_v
            real(r8),pointer ::d_up
            real(r8),pointer ::d_vp
            real(r8),pointer ::d_vb
            real(r8),pointer ::d_dlu
            real(r8),pointer ::d_sbcx
            real(r8),pointer ::d_bbcx
            real(r8),pointer ::d_wka
            real(r8),pointer ::d_a8
            real(r8),pointer ::d_b8
            real(r8),pointer ::d_c8
            real(r8),pointer ::d_d8
            real(r8),pointer ::d_e8
            real(r8),pointer ::d_f8
            real(r8),pointer ::d_work
      end type

      type bclinc5_cu_slave_param
            integer, dimension(3) :: cuBlock
            integer, dimension(3) :: cuGrid
            integer :: rank
            logical :: d_iau
            real(r8) ::d_afc1
            real(r8) ::d_afc2
            real(r8) ::d_lamda
            real(r8),pointer ::d_viv
            real(r8),pointer ::d_dzp
            real(r8),pointer ::d_u
            real(r8),pointer ::d_v
            real(r8),pointer ::d_up
            real(r8),pointer ::d_ub
            real(r8),pointer ::d_wka
            real(r8),pointer ::d_ohbu
            real(r8),pointer ::d_vtf
            real(r8),pointer ::d_utf
            real(r8),pointer ::d_work
            real(r8),pointer ::d_ucheck
            real(r8),pointer ::d_vcheck
            real(r8),pointer ::vv_iau
            real(r8),pointer ::uu_iau
      end type

      type bclinc6_cu_slave_param
            integer, dimension(3) :: cuBlock
            integer, dimension(3) :: cuGrid
            integer :: d_isc
            real(r8),pointer ::d_viv
            real(r8),pointer ::d_epea
            real(r8),pointer ::d_epeb
            real(r8),pointer ::d_epla
            real(r8),pointer ::d_eplb
            real(r8),pointer ::d_up
            real(r8),pointer ::d_vp
            real(r8),pointer ::d_dlu
            real(r8),pointer ::d_dlv
            real(r8),pointer ::d_gg
            real(r8),pointer ::d_fcor
            real(r8),pointer ::d_wka1
            real(r8),pointer ::d_wka2
            real(r8),pointer ::d_dxur
            real(r8),pointer ::d_dyur
            integer,pointer ::d_kmu
      end type

      type bclinc71_cu_slave_param
            integer, dimension(3) :: cuBlock
            integer, dimension(3) :: cuGrid
            integer :: rank
            integer :: padd
            real(r8) ::d_dtc
            integer,pointer  :: d_kmu
            real(r8),pointer ::d_viv
            real(r8),pointer ::d_odzt
            real(r8),pointer ::d_odzp
            real(r8),pointer ::d_akmu
            real(r8),pointer ::d_u
            real(r8),pointer ::d_up
            real(r8),pointer ::d_v
            real(r8),pointer ::d_vp
            real(r8),pointer ::d_dlu
            real(r8),pointer ::d_dlv
            real(r8),pointer ::d_sbcx
            real(r8),pointer ::d_bbcx
      end type

      type bclinc72_cu_slave_param
            integer, dimension(3) :: cuBlock
            integer, dimension(3) :: cuGrid
            real(r8) ::d_dtc
            integer,pointer  :: d_kmu
            real(r8),pointer ::d_viv
            real(r8),pointer ::d_odzt
            real(r8),pointer ::d_odzp
            real(r8),pointer ::d_akmu
            real(r8),pointer ::d_v
            real(r8),pointer ::d_vp
            real(r8),pointer ::d_bbcy
            real(r8),pointer ::d_dlv
            real(r8),pointer ::d_sbcy
            real(r8),pointer ::d_a8
            real(r8),pointer ::d_b8
            real(r8),pointer ::d_c8
            real(r8),pointer ::d_d8
            real(r8),pointer ::d_e8
            real(r8),pointer ::d_f8
            real(r8),pointer ::d_wka2
      end type

      type bclinc81_cu_slave_param
            integer, dimension(3) :: cuBlock
            integer, dimension(3) :: cuGrid
            real(r8),pointer ::d_viv
            real(r8),pointer ::d_dzp
            real(r8),pointer ::d_u1
            real(r8),pointer ::d_u
            real(r8),pointer ::d_ub
            real(r8),pointer ::d_utf
            real(r8),pointer ::d_ohbu
            real(r8),pointer ::d_work
      end type
      type bclinc82_cu_slave_param
            integer, dimension(3) :: cuBlock
            integer, dimension(3) :: cuGrid
            real(r8),pointer ::d_viv
            real(r8),pointer ::d_dzp
            real(r8),pointer ::d_v1
            real(r8),pointer ::d_v
            real(r8),pointer ::d_vb
            real(r8),pointer ::d_vtf
            real(r8),pointer ::d_ohbu
            real(r8),pointer ::d_work
      end type

      type vinteg_cu_slave_param
            integer, dimension(3) :: cuBlock
            integer, dimension(3) :: cuGrid
            integer :: rank
            integer :: vinegtype
            real(r8),pointer ::d_u
            real(r8),pointer ::d_v
            real(r8),pointer ::d_wka
            real(r8),pointer ::d_work
            real(r8),pointer ::d_dzp
            real(r8),pointer ::d_ohbu
            real(r8),pointer ::d_viv
      end type

      type(bclinc12_cu_slave_param) :: bclinc12_cu_param
      type(bclinc11_cu_slave_param) :: bclinc11_cu_param
      type(bclinc6_cu_slave_param) :: bclinc6_cu_param
      type(bclinc71_cu_slave_param) :: bclinc71_cu_param
      type(bclinc72_cu_slave_param) :: bclinc72_cu_param
      type(bclinc81_cu_slave_param) :: bclinc81_cu_param
      type(bclinc82_cu_slave_param) :: bclinc82_cu_param
      type(bclinc2_cu_slave_param) :: bclinc2_cu_param
      type(bclinc3_cu_slave_param) :: bclinc3_cu_param
      type(bclinc_invtriu_cu_slave_param) :: bclinc_invtriu_param
      type(bclinc4_cu_slave_param) :: bclinc4_cu_param
      type(bclinc4_0_cu_slave_param) :: bclinc4_0_cu_param
      type(bclinc5_cu_slave_param) :: bclinc5_cu_param
      type(vinteg_cu_slave_param) :: vinteg_param

      external :: slave_bclinc12_cu_slave_kernel
      external :: slave_bclinc11_cu_slave_kernel
      external :: slave_bclinc6_cu_slave_kernel
      external :: slave_bclinc71_cu_slave_kernel
      external :: slave_bclinc72_cu_slave_kernel
      external :: &
              slave_bclinc81_cu,slave_bclinc82_cu
      external :: slave_bclinc2_cu_slave_kernel
      external :: slave_bclinc3_cu_slave_kernel
      external :: slave_bclinc_invtriu_kernel
      external :: slave_bclinc4_cu_slave_kernel
      external :: slave_bclinc4_cu_0_slave_kernel
      external :: slave_bclinc4_1_cu_slave_kernel
      external :: slave_bclinc5_cu_slave_kernel
      external :: slave_bclinc5_1_cu_slave_kernel
      external :: slave_bclinc5_2_cu_slave_kernel
      external :: slave_vinteg_cu_slave_kernel
      external :: slave_vinteg_cu_slave_kernel_old

!-------------------------------------------------------------------------------
! tracer
!-------------------------------------------------------------------------------

      type tracer_del4_1_cu_param
            integer, dimension(3) :: cuBlock
            integer, dimension(3) :: cuGrid
            integer :: n
            integer ,pointer :: d_kmtn
            integer ,pointer :: d_kmts
            integer ,pointer :: d_kmte
            integer ,pointer :: d_kmtw
            integer ,pointer :: d_kmt
            real(r8),pointer :: d_atb
            real(r8),pointer :: d_dtn
            real(r8),pointer :: d_dts
            real(r8),pointer :: d_dte
            real(r8),pointer :: d_dtw
            real(r8),pointer :: d_d2tk
            real(r8),pointer :: d_ahf
      end type

      type tracer_del4_2_cu_param
            integer, dimension(3) :: cuBlock
            integer, dimension(3) :: cuGrid
            integer :: n
            integer ,pointer :: d_kmtn
            integer ,pointer :: d_kmts
            integer ,pointer :: d_kmte
            integer ,pointer :: d_kmtw
            integer ,pointer :: d_kmt
            real(r8)         :: d_ah
            real(r8)         :: d_dwndmix
            real(r8),pointer :: d_tf
            !real(r8),pointer :: d_dx
            real(r8),pointer :: d_adv_tt
            real(r8),pointer :: d_vit
            real(r8),pointer :: d_akt
            real(r8),pointer :: d_wkc
            real(r8),pointer :: d_dtn
            real(r8),pointer :: d_dts
            real(r8),pointer :: d_dte
            real(r8),pointer :: d_dtw
            real(r8),pointer :: d_d2tk
            real(r8),pointer :: d_atb
            real(r8),pointer :: d_factor_vis_tsz0
#if (defined ISO)
            real(r8),pointer :: d_ahisop
            real(r8),pointer :: d_k3
#endif
      end type

      type invtrit_cu_param
            integer, dimension(3) :: cuBlock
            integer, dimension(3) :: cuGrid
            integer :: rank
            integer ,pointer :: d_kmt
            real(r8) :: d_dts
            real(r8),pointer :: d_wkc
            real(r8),pointer :: d_vtl
            real(r8),pointer :: d_stf
            real(r8),pointer :: d_vit
            real(r8),pointer :: d_odzt
            real(r8),pointer :: d_odzp
      end type

      type tracer_9_cu_param
            integer, dimension(3) :: cuBlock
            integer, dimension(3) :: cuGrid
            integer :: n
            integer :: rank
            real(r8) :: d_dts
            real(r8),pointer :: d_at
            real(r8),pointer :: d_vtl
            real(r8),pointer :: d_atb
            real(r8),pointer :: d_tend
            real(r8),pointer :: d_vit
            !real(r8),pointer :: d_dt_diff
            real(r8),pointer :: d_ori
            real(r8),pointer :: d_stf
            real(r8),pointer :: d_odzp
      end type
    type tracer_vtl_cu_parm
        integer, dimension(3) :: cuGrid
        integer, dimension(3) :: cuBlock
        real(r8),pointer :: d_vtl
        real(r8),pointer :: d_vtl_ori
    end type

    type tracer_wkc_cu_parm
        integer, dimension(3) :: cuGrid
        integer, dimension(3) :: cuBlock
        real(r8)              :: d_ahv
        real(r8),pointer :: d_wkc
        real(r8),pointer :: d_akt
    end type

    type tracer_norm2_cu_parm
        integer, dimension(3) :: cuGrid
        integer, dimension(3) :: cuBlock
        real(r8),pointer      :: d_err_norm2
        real(r8),pointer      :: d_tarea
        real(r8),pointer      :: d_net
        real(r8),pointer      :: d_vit
    end type

    type tracer_norm3_cu_parm
        integer, dimension(3) :: cuGrid
        integer, dimension(3) :: cuBlock
        real(r8)              :: d_err_norm2
        real(r8),pointer      :: d_tf
        real(r8),pointer      :: d_vit
    end type

    type tracer_norm4_cu_parm
        integer, dimension(3) :: cuGrid
        integer, dimension(3) :: cuBlock
        integer, pointer      :: d_kmt
        integer               :: padd
        real(r8)              :: d_od0cp
        real(r8)              :: d_gamma
        real(r8)              :: d_aidif
        real(r8),pointer      :: d_stf
        real(r8),pointer      :: d_swv
        real(r8),pointer      :: d_nswv
        real(r8),pointer      :: d_seaice
        real(r8),pointer      :: d_sst
        real(r8),pointer      :: d_atb
        real(r8),pointer      :: d_odzp
        real(r8),pointer      :: d_tf
        real(r8),pointer      :: d_net
#ifdef SST_restore
        real(r8),pointer      :: d_sst_restore
#endif
#ifdef COUP
        real(r8),pointer      :: d_tsf
#endif
    end type

      type(invtrit_cu_param)  :: invtrit_param
      type(tracer_9_cu_param)  :: tracer_9_param
      type(tracer_del4_1_cu_param)  :: tracer_del4_1_param
      type(tracer_del4_2_cu_param)  :: tracer_del4_2_param
      type(tracer_vtl_cu_parm)  :: tracer_vtl_parm
      type(tracer_wkc_cu_parm)  :: tracer_wkc_parm
      type(tracer_norm2_cu_parm)  :: tracer_norm2_parm
      type(tracer_norm3_cu_parm)  :: tracer_norm3_parm
      type(tracer_norm4_cu_parm)  :: tracer_norm4_parm

      external :: slave_invtrit
      external :: slave_tracer_9
      external :: slave_tracer_del4_1
      external :: slave_tracer_del4_2
      external :: slave_tracer_vtl_kernel
      external :: slave_tracer_wkc_kernel
      external :: slave_tracer_norm2_kernel
      external :: slave_tracer_norm3_kernel
      external :: slave_tracer_norm4_kernel

      type tracer_upwell_0_parm
            integer, dimension(3) :: cuGrid
            integer, dimension(3) :: cuBlock
            integer          :: rank
            integer          :: pad
            real(r8),pointer :: d_au0
            real(r8),pointer :: d_aus
            real(r8),pointer :: d_auw
            real(r8),pointer :: d_ausw
            real(r8),pointer :: d_stf
            real(r8),pointer :: d_ohbu
            real(r8),pointer :: d_wkd
            real(r8),pointer :: d_uk
            real(r8),pointer :: d_wkb
            real(r8),pointer :: d_vk
            real(r8),pointer :: d_work
      end type

      type tracer_upwell_parm
            integer, dimension(3) :: cuGrid
            integer, dimension(3) :: cuBlock
            integer          :: rank
            integer,pointer  :: d_kmt
            real(r8),pointer :: d_dzp
            real(r8),pointer :: d_stf
            real(r8),pointer :: d_hts
            real(r8),pointer :: d_htw
            real(r8),pointer :: d_ohbt
            real(r8),pointer :: d_tarea_r
            real(r8),pointer :: d_uk
            real(r8),pointer :: d_vit
            real(r8),pointer :: d_vk
            real(r8),pointer :: d_wka
            real(r8),pointer :: d_work
            real(r8),pointer :: d_ws
      end type

      type tracer_upwell_1_parm
            integer, dimension(3) :: cuGrid
            integer, dimension(3) :: cuBlock
            integer          :: rank
            integer          :: padd
            real(r8),pointer :: d_vit
            real(r8),pointer :: d_dzp
            real(r8),pointer :: d_wka
            real(r8),pointer :: d_work
      end type

      type tracer_upwell_2_parm
            integer, dimension(3) :: cuGrid
            integer, dimension(3) :: cuBlock
            integer          :: rank
            integer          :: k
            real(r8),pointer :: d_ws
            real(r8),pointer :: d_vit
            real(r8),pointer :: d_stf
            real(r8),pointer :: d_dzp
            real(r8),pointer :: d_ohbt
            real(r8),pointer :: d_wka
            real(r8),pointer :: d_work
            real(r8),pointer :: d_work1
            real(r8),pointer :: d_wsbak
      end type

      type advection1_0_cu_slave_param
            integer, dimension(3) :: cuGrid
            integer, dimension(3) :: cuBlock
            integer :: Str
            integer :: rank
            real(r8),pointer :: d_wkb
            real(r8),pointer :: d_wkd
            real(r8),pointer :: d_u_wface
            real(r8),pointer :: d_v_sface
            real(r8),pointer :: d_hts
            real(r8),pointer :: d_htw
            real(r8),pointer :: d_dxu
            real(r8),pointer :: d_dyu
            real(r8),pointer :: d_adv_tt
            real(r8),pointer :: d_tf
      end type

      type advection1_cu_slave_param
            integer, dimension(3) :: cuGrid
            integer, dimension(3) :: cuBlock
            integer :: d_mtracer
            integer :: Str
            integer :: d_nss
            integer :: rank
            real(r8),pointer :: d_u_wface
            real(r8),pointer :: d_v_sface
            real(r8),pointer :: d_ws
            real(r8),pointer :: d_at
            real(r8),pointer :: d_adv_tt
            real(r8),pointer :: d_tarea_r
            real(r8),pointer :: d_odzp
            real(r8),pointer :: d_ax
            real(r8),pointer :: d_ay
            real(r8),pointer :: d_az
      end type

      type advection2_cu_slave_param
            integer, dimension(3) :: cuGrid
            integer, dimension(3) :: cuBlock
            integer :: d_ntra
            integer :: d_mtracer
            integer :: Str
            integer :: rank
            real(r8),pointer :: d_wkb
            real(r8),pointer :: d_wkd
            real(r8),pointer :: d_u_wface
            real(r8),pointer :: d_v_sface
            real(r8),pointer :: d_hts
            real(r8),pointer :: d_htw
            real(r8),pointer :: d_dxu
            real(r8),pointer :: d_dyu
            real(r8),pointer :: d_ws
            real(r8),pointer :: d_at
            real(r8),pointer :: d_adv_tt
            real(r8),pointer :: d_tarea_r
            real(r8),pointer :: d_odzp
            real(r8),pointer :: d_ax
            real(r8),pointer :: d_ay
            real(r8),pointer :: d_az
            real(r8),pointer :: d_hun
            real(r8),pointer :: d_hue
            real(r8),pointer :: d_at00
            real(r8),pointer :: d_atmax
            real(r8),pointer :: d_atmin
            real(r8),pointer :: d_odzt
            real(r8),pointer :: d_vit
            real(r8),pointer :: d_dxt
            real(r8),pointer :: d_dyt
            real(r8),pointer :: d_dxtr
            real(r8),pointer :: d_dytr
            real(r8),pointer :: d_adv_x
            real(r8),pointer :: d_adv_y
            real(r8),pointer :: d_adv_z
            real(r8),pointer :: d_zkp
            real(r8),pointer :: d_flxqRtmp
            real(r8),pointer :: d_flxqLtmp
            real(r8)         :: d_dts
      end type

      type(tracer_upwell_0_parm) :: tracer_upwell_0_param
      type(tracer_upwell_parm) :: tracer_upwell_param
      type(tracer_upwell_1_parm) :: tracer_upwell_1_param
      type(tracer_upwell_2_parm) :: tracer_upwell_2_param
      type(advection1_0_cu_slave_param) :: advection1_0_param
      type(advection1_cu_slave_param) :: advection1_param
      type(advection2_cu_slave_param) :: advection2_param
      external :: slave_advection_tracer_1_0,slave_advection_tracer_1,slave_advection_tracer_2,slave_advection_tracer_3
      external :: slave_advection_tracer_4
      external :: slave_advection_tracer_5
      external :: slave_advection_tracer_6
      external :: slave_advection_tracer_7
      external :: slave_advection_tracer_8
      external :: slave_advection_tracer_9
      external :: slave_tracer_upwell_cu_0
      external :: slave_tracer_upwell_cu
      external :: slave_tracer_upwell_cu_1
      external :: slave_tracer_upwell_cu_2

end module sunway_mod
