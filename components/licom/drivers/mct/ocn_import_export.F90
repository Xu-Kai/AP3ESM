module ocn_import_export

   use POP_CommMod
   use POP_CplIndices
   use mct_mod
   use blocks
   use domain
   use grid
   use param_mod
   use pconst_mod
   use shr_msg_mod
   use shr_sys_mod
   use control_mod
   use constant_mod
   use buf_mod

   implicit none
   public
   save

   ! accumulated sum of send buffer quantities for averaging before being sent
   real (r8), dimension(:,:,:,:), allocatable ::  SBUFF_SUM 

   real (r8) :: tlast_coupled 

contains

!*************************************************************************
!ROUTINE: ocn_import_mct
!         lihuimin 2012.7.16
!
subroutine ocn_import(x2o_o)
   ! type(mct_aVect)   , intent(inout) :: x2o_o
   real(r8)           , intent(inout) :: x2o_o(:,:)
   ! local
   integer :: j_begin, j_end, iblock
   type (block) :: this_block          ! block information for current block
   integer(8):: tmproff !LPF20150908

    n=0
   do iblock = 1, nblocks_clinic
       this_block = get_block(blocks_clinic(iblock),iblock)
       do j=this_block%jb, this_block%je
       do i=this_block%ib, this_block%ie
          n=n+1
          !--- states ---
          ifrac(i,j,iblock) = x2o_o(index_x2o_Si_ifrac,n) ! ice fraction
          patm (i,j,iblock) = x2o_o(index_x2o_Sa_pslv,n)  ! sea level pressure index_x2o_Sa_pslv 
          !--- fluxes ---
          !taux (i,j,iblock) = x2o_o(index_xao_Faox_taux,n)  ! surface stress, zonal
          !tauy (i,j,iblock) = x2o_o(index_xao_Faox_tauy,n)  ! surface stress, merid
          taux (i,j,iblock) = x2o_o(index_x2o_Foxx_taux,n)  ! surface stress, zonal
          tauy (i,j,iblock) = x2o_o(index_x2o_Foxx_tauy,n)  ! surface stress, merid
          netsw(i,j,iblock) = x2o_o(index_x2o_Foxx_swnet,n) ! net sw rad
          sen  (i,j,iblock) = x2o_o(index_x2o_Foxx_sen,n)   ! sensible
          lwup (i,j,iblock) = x2o_o(index_x2o_Foxx_lwup,n)  ! long-wave up
          lwdn (i,j,iblock) = x2o_o(index_x2o_Faxa_lwdn,n)  ! long-wave down
          melth(i,j,iblock) = x2o_o(index_x2o_Fioi_melth,n) ! melt heat
          salt (i,j,iblock) = x2o_o(index_x2o_Fioi_salt,n)  ! salinity flux
          prec (i,j,iblock) = x2o_o(index_x2o_Faxa_rain,n)+x2o_o(index_x2o_Faxa_snow,n)  !index_x2o_Faxa_prec 
!          prec (i,j,iblock) = x2o_o(index_x2o_Faxa_prec,n)  !index_x2o_Faxa_prec 
          evap (i,j,iblock) = x2o_o(index_x2o_Foxx_evap,n)  ! evaporation
          meltw(i,j,iblock) = x2o_o(index_x2o_Fioi_meltw,n) ! melt water
          roff (i,j,iblock) = x2o_o(index_x2o_Foxx_rofl,n)  ! runoff
          snow1(i,j,iblock) = x2o_o(index_x2o_Faxa_snow,n)  !snow flux due to the melting of snow
          iceoff(i,j,iblock)= x2o_o(index_x2o_Foxx_rofi,n) ! runoff due to ice 
          duu10n(i,j,iblock) = x2o_o(index_x2o_So_duu10n,n)  ! 10m wind speed squared
       end do
       end do
   end do

!LPF20150908
        if(mytid==0) write(16,*) 'test cpl due to runoff inconsistence'
    do iblock = 1, nblocks_clinic
       this_block = get_block(blocks_clinic(iblock),iblock)
        do j=this_block%jb, this_block%je
        do i=this_block%ib, this_block%ie
            tmproff=roff(i,j,iblock)*1.D17
            roff(i,j,iblock)=tmproff*1.D-17
           enddo
        enddo
     enddo
!LPF20150908


       lat1= LATVAP*evap ! latent (derive from evap)
!       write(*,*) 'n= import',n

 end subroutine ocn_import



!*************************************************************************
!ROUTINE: ocn_export_mct
!         lihuimin 2012.7.16
!
 subroutine ocn_export(o2x_o)
   ! type(mct_aVect)   , intent(inout) :: o2x_o
   real(r8)           , intent(inout) :: o2x_o(:,:)
   ! local
   integer :: j_begin,j_end, iblock
   type (block) :: this_block          ! block information for current block
   real (r8) :: ek0

      n=0
! lihuimin, 2012.8.7, consider jst_global
 do iblock =1 , nblocks_clinic
       this_block = get_block(blocks_clinic(iblock),iblock)
       do j=this_block%jb, this_block%je
       do i=this_block%ib, this_block%ie
         n=n+1
         o2x_o(index_o2x_So_t,n)    = T_cpl   (i,j,iblock) ! temperature
         o2x_o(index_o2x_So_s,n)    = S_cpl   (i,j,iblock) ! salinity
         o2x_o(index_o2x_So_u,n)    = U_cpl   (i,j,iblock) ! velocity, zonal
         o2x_o(index_o2x_So_v,n)    = V_cpl   (i,j,iblock) ! velocity, meridional
         o2x_o(index_o2x_So_dhdx,n) = dhdx(i,j,iblock) ! surface slope, zonal
         o2x_o(index_o2x_So_dhdy,n) = dhdy(i,j,iblock) ! surface slope, meridional
         o2x_o(index_o2x_Fioo_q,n)    = qheat   (i,j,iblock) ! heat of fusion xor melt pot
         o2x_o(index_o2x_Sf_ofrac,n) =  vit(i,j,1,iblock)
#ifdef USE_OCN_CARBON
!         buffs(n,cpl_fields_o2c_co2)    = co2_cpl(i,j) ! state: air-sea CO2 of ocean  ~ mol
#endif
      enddo
      enddo
  end do
!       write(*,*) 'n= export',n
 end subroutine ocn_export

end module ocn_import_export
