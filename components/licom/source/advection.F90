!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

 module advection

!BOP
! !MODULE: advection
!
! !DESCRIPTION:
!  This module contains arrays and variables necessary for performing
!  advection of momentum and tracer quantities.  Currently, the
!  module supports leapfrog centered advection of momentum and
!  both leapfrog centered advection and third-order upwinding of
!  tracers.
!
! !REVISION HISTORY:
!  SVN:$Id: advection.F90 28439 2011-05-18 21:40:58Z njn01 $

! !USES:
   use precision_mod
   use param_mod
   use pconst_mod
   use constant_mod
   use grid
   use blocks
   use tracer_mod, only : ax,ay,az
   use dyn_mod, only : h0
   use LICOM_Error_mod

   implicit none
   private
   save

   public :: advection_momentum, advection_tracer
!

!EOP
!BOC

    contains

!-----------------------------------------------------------------------
!EOC
      subroutine advection_momentum(uuu,vvv,www,adv_uu,adv_vv,iblock)
!
      real(r8) ,intent(in) :: uuu(imt,jmt,km), vvv(imt,jmt,km),www(imt,jmt,km)
      real(r8), intent(out):: adv_uu(imt,jmt,km), adv_vv(imt,jmt,km)
      real(r8) :: u_wface(imt,jmt,km), v_sface(imt,jmt,km)
      real(r8) :: adv_z1, adv_z2,adv_z3,adv_z4
      integer, intent(in)  :: iblock
!
     adv_uu = c0
     adv_vv = c0
!
     if ( adv_momentum(1:8) == 'centered' ) then
      do k=1, km
      do j= 1,jmt-1
      do i= 2,imt
            u_wface(i,j,k) = (uuu(i-1,j,k) + uuu(i,j,k))*P25*hue(i-1,j,iblock)
            v_sface(i,j,k) = (vvv(i,j,k) + vvv(i,j+1,k))*P25*hun(i,j+1,iblock)
      end do
      end do
      end do
     else if ( adv_momentum(1:4) == 'flux' ) then
      do k=1, km
      do j= 1,jmt-1
      do i= 2,imt
            u_wface(i,j,k) = (uuu(i-1,j,k)*dyu(i-1,j,iblock) + uuu(i,j,k)*dyu(i,j,iblock))*P25
            v_sface(i,j,k) = (vvv(i,j,k)*dxu(i,j,iblock) + vvv(i,j+1,k)*dxu(i,j+1,iblock))*P25
      end do
      end do
      end do
!     else
!         write(*,*) "adv_momentum222 =", adv_momentum
!        call exit_licom(sigAbort,'The false advection option for tracer')
     end if
!
      if ( adv_momentum(1:8) == 'centered' ) then
        do k = 1,km
        do j= 3,jmt-2
        do i= 3,imt-2
           adv_uu(i,j,k) = (-u_wface(i  ,j,k)*(uuu(i  ,j,k)-uuu(i-1,j,k))                   &
                            -u_wface(i+1,j,k)*(uuu(i+1,j,k)-uuu(i  ,j,k))                   &
                            -v_sface(i,j  ,k)*(uuu(i,j+1,k)-uuu(i,j  ,k))                   &
                            -v_sface(i,j-1,k)*(uuu(i,j  ,k)-uuu(i,j-1,k)))*uarea_r(i,j,iblock)
           adv_vv(i,j,k) = (-u_wface(i  ,j,k)*(vvv(i  ,j,k)-vvv(i-1,j,k))                   &
                            -u_wface(i+1,j,k)*(vvv(i+1,j,k)-vvv(i  ,j,k))                   &
                            -v_sface(i,j  ,k)*(vvv(i,j+1,k)-vvv(i,j  ,k))                   &
                            -v_sface(i,j-1,k)*(vvv(i,j  ,k)-vvv(i,j-1,k)))*uarea_r(i,j,iblock)
!
           if (k==1 )then
               adv_z1=0.0D0
               adv_z3=0.0D0
           else
               adv_z1=www (I,J,K)* (uuu(I,J,K -1) - uuu(I,J,K))
               adv_z3=www (I,J,K)* (vvv(I,J,K -1) - vvv(I,J,K))
           end if
!
           if (k==km )then
                adv_z2=0.0D0
                adv_z4=0.0D0
           else
                adv_z2=www(I,J,K+1)*(uuu(I,J,K ) - uuu(I,J,K+1))
                adv_z4=www(I,J,K+1)*(vvv(I,J,K ) - vvv(I,J,K+1))
           end if
           adv_uu(i,j,k) = adv_uu(i,j,k) - P5*ODZP(K)* (adv_z1+adv_z2)
           adv_vv(i,j,k) = adv_vv(i,j,k) - P5*ODZP(K)* (adv_z3+adv_z4)
        end do
        end do
        end do
      else if ( adv_momentum(1:4) == 'flux' ) then
        do k = 1,km
        do j= 3,jmt-2
        do i= 3,imt-2
           adv_uu(i,j,k) = (-u_wface(i  ,j,k)*(uuu(i  ,j,k)+uuu(i-1,j,k))    &
                            +u_wface(i+1,j,k)*(uuu(i+1,j,k)+uuu(i  ,j,k))    &
                            -v_sface(i,j-1,k)*(uuu(i,  j,k)+uuu(i,j-1,k))    &
                            +v_sface(i,j  ,k)*(uuu(i,j+1,k)+uuu(i,  j,k)))*uarea_r(i,j,iblock)
 
           adv_vv(i,j,k) = (-u_wface(i  ,j,k)*(vvv(i,j  ,k)+vvv(i-1,j,k))    &
                            +u_wface(i+1,j,k)*(vvv(i+1,j,k)+vvv(i  ,j,k))    &
                            -v_sface(i,j-1,k)*(vvv(i,  j,k)+vvv(i,j-1,k))    &
                            +v_sface(i,j  ,k)*(vvv(i,j+1,k)+vvv(i,  j,k)))*uarea_r(i,j,iblock)
 
           if (k ==1 ) then
               adv_z1=0.0D0
               adv_z3=0.0D0
           else
               adv_z1=www (I,J,K)* (uuu(I,J,K -1) + uuu(I,J,K))*P5
               adv_z3=www (I,J,K)* (vvv(I,J,K -1) + vvv(I,J,K))*P5
           end if
!
           if (k == km ) then
                adv_z2=0.0D0
                adv_z4=0.0D0
           else
                adv_z2=www(I,J,K+1)*(uuu(I,J,K ) + uuu(I,J,K+1))*P5
                adv_z4=www(I,J,K+1)*(vvv(I,J,K ) + vvv(I,J,K+1))*P5
           end if
           adv_uu(i,j,k) = adv_uu(i,j,k) - ODZP(K)* (adv_z2-adv_z1)
           adv_vv(i,j,k) = adv_vv(i,j,k) - ODZP(K)* (adv_z4-adv_z3)
        end do
        end do
        end do

      else
        write(6,*) "adv_momentum =", adv_momentum
        call exit_licom(sigAbort,'The false advection option for momentum')
      end if
!
      end subroutine advection_momentum
!
!
!
      !subroutine advection_tracer(uuu,vvv,www,ttt,adv_tt,iblock,mtracer, this_block,ax,ay,az)
      subroutine advection_tracer(uuu,vvv,www,ttt,adv_tt,iblock,mtracer, this_block) !LPF20160811
!
      real(r8) ,intent(in) :: uuu(imt,jmt,km), vvv(imt,jmt,km),www(imt,jmt,km+1),ttt(imt,jmt,km)
      integer, intent(in)  :: iblock,mtracer
      real(r8), intent(out):: adv_tt(imt,jmt,km)
      real(r8) :: u_wface(imt,jmt,km),v_sface(imt,jmt,km)
      real(r8) :: adv_z1, adv_z2, sumv
      !real(r8) :: ucfl, vcfl, wcfl
      real(r8)    :: LAMDA,wt1,wt2,adv_z, ek1,ek0,ek2,sflux(imt,jmt),nflux(imt,jmt)
      real(r8),dimension(:,:,:),allocatable :: adv_xy1,adv_xy2,adv_xy3,adv_xy4, at00
      real(r8),dimension(:,:,:),allocatable :: adv_zz,atz, adv_za,adv_zb1,adv_zb2,atmaxz,atminz
      real(r8),dimension(:,:,:),allocatable :: adv_x0,adv_y0,atmax,atmin,adv_xx,adv_yy
      real(r8),dimension(:,:,:),allocatable :: adv_c1,adv_c2,adv_zc
      type(block) :: this_block
!
     adv_tt=0.0_r8
!
      do k=1, km
      do j=1, jmt
      do i=2, imt-1
         if ( adv_tracer(1:8) == 'centered' .or. adv_tracer(1:5) == 'tspas') then
            v_sface(i,j,k) = (vvv(i,j,k) + vvv(i+1,j,k))*hts(i,j,iblock)*P25
         else if ( adv_tracer(1:4) == 'flux' ) then
            v_sface(i,j,k) = (vvv(i,j,k)*dxu(i,j,iblock) + vvv(i+1,j,k)*dxu(i+1,j,iblock))*P25
         end if
      end do
      end do
      do j= 2,jmt-1
      do i= 1,imt
         if ( adv_tracer(1:8) == 'centered' .or. adv_tracer(1:5) == 'tspas') then
            u_wface(i,j,k) = (uuu(i,j-1,k) + uuu(i,j,k))*htw(i,j,iblock)*P25
         else if ( adv_tracer(1:4) == 'flux' ) then
            u_wface(i,j,k) = (uuu(i,j-1,k)*dyu(i,j-1,iblock) + uuu(i,j,k)*dyu(i,j,iblock))*P25
         end if
      end do
      end do
      end do
!
      if ( adv_tracer(1:8) == 'centered' ) then
        do k = 1,km
        do j= 3,jmt-2
        do i= 3,imt-2
           adv_tt(i,j,k) = (-u_wface(i  ,j,k)*(ttt(i  ,j,k)-ttt(i-1,j,k))                      &
                            -u_wface(i+1,j,k)*(ttt(i+1,j,k)-ttt(i  ,j,k))                      &
                            -v_sface(i,j  ,k)*(ttt(i,j+1,k)-ttt(i,j  ,k))                      &
                            -v_sface(i,j-1,k)*(ttt(i,j  ,k)-ttt(i,j-1,k)))*tarea_r(i,j,iblock)
           ax(i,j,k,mtracer,iblock) =  ax(i,j,k,mtracer,iblock) +  &
                               (-u_wface(i  ,j,k)*(ttt(i  ,j,k)-ttt(i-1,j,k))                      &
                                -u_wface(i+1,j,k)*(ttt(i+1,j,k)-ttt(i  ,j,k)))*tarea_r(i,j,iblock)/dble(NSS)
           ay(i,j,k,mtracer,iblock) =  ay(i,j,k,mtracer,iblock) +  &
                               (-v_sface(i,j  ,k)*(ttt(i,j+1,k)-ttt(i,j  ,k))                      &
                                -v_sface(i,j-1,k)*(ttt(i,j  ,k)-ttt(i,j-1,k)))*tarea_r(i,j,iblock)/dble(NSS)
!
           if (k==1 )then
               adv_z1=0.0D0
           else
               adv_z1=www (I,J,K)* (ttt(I,J,K -1) - ttt(I,J,K))
           end if
!
           if (k==km )then
                adv_z2=0.0D0
           else
                adv_z2=www(I,J,K+1)*(ttt(I,J,K ) - ttt(I,J,K+1))
           end if
!
           adv_tt(i,j,k) = adv_tt(i,j,k) - P5*ODZP(K)* (adv_z1+adv_z2)
           az(i,j,k,mtracer,iblock) =  az(i,j,k,mtracer,iblock) - P5*ODZP(K)* (adv_z1+adv_z2)/dble(NSS)
        end do
        end do
        end do
      else if ( adv_tracer(1:4) == 'flux' ) then
        do k = 1,km
        do j= 3,jmt-2
        do i= 3,imt-2
           adv_tt(i,j,k) = (-u_wface(i  ,j,k)*(ttt(i  ,j,k)+ttt(i-1,j,k))   &
                            +u_wface(i+1,j,k)*(ttt(i+1,j,k)+ttt(i  ,j,k))   &
                            -v_sface(i,j-1,k)*(ttt(i,  j,k)+ttt(i,j-1,k))   &
                            +v_sface(i,j  ,k)*(ttt(i,j+1,k)+ttt(i,  j,k)))*tarea_r(i,j,iblock)
!
           if (k==1 )then
               adv_z1=0.0D0
           else
               adv_z1=www (I,J,K)* (ttt(I,J,K -1) + ttt(I,J,K))*P5
           end if
!
           if (k==km )then
                adv_z2=0.0D0
           else
                adv_z2=www(I,J,K+1)*(ttt(I,J,K ) + ttt(I,J,K+1))*P5
           end if
!
           adv_tt(i,j,k) = adv_tt(i,j,k) - ODZP(K)* (adv_z2-adv_z1)
!
        end do
        end do
        end do
      else if (adv_tracer(1:5) == 'tspas' ) then
!
      allocate ( adv_c1(imt,jmt,km),adv_c2(imt,jmt,km),adv_zc(imt,jmt,km))
      allocate ( adv_x0(imt,jmt,km),adv_y0(imt,jmt,km))
      allocate ( adv_xx(imt,jmt,km), adv_yy(imt,jmt,km))
      allocate ( adv_xy1(imt,jmt,km), adv_xy2(imt,jmt,km), adv_xy3(imt,jmt,km), adv_xy4(imt,jmt,km))
      allocate ( at00(imt,jmt,km), atmax(imt,jmt,km), atmin(imt,jmt,km) )
      allocate (adv_zz(imt,jmt,km), adv_za(imt,jmt,km),adv_zb1(imt,jmt,km) )
      allocate (adv_zb2(imt,jmt,km),atmaxz(imt,jmt,km), atminz(imt,jmt,km), atz(imt,jmt,km))
!
!wr cfl test
!      DO K = 1,km
 !        ucfl = minval(abs(dxt(:,:,iblock)/(uuu(:,:,k)+1d-20)))
  !       vcfl = minval(abs(dyt(:,:,iblock)/(vvv(:,:,k)+1d-20)))
    !     wcfl = minval(abs(dzp(k)/(www(:,:,k)+1d-20)))
   !      write(6,*)"cfl,k,u,v,w",k,ucfl,vcfl,wcfl
     ! ENDDO


      nflux=0.0D0
      sflux=0.0D0
      DO K = 1,km
        DO J = 2, jmt-1
          DO I = 2, imt-1
            adv_x0(i,j,k)=(ttt(i+1,j,k)+ttt(i,j,k))*u_wface(i+1,j,k)*tarea_r(i,j,iblock) &
                       -(ttt(i,j,k)+ttt(i-1,j,k))*u_wface(i,j,k)*tarea_r(i,j,iblock)
            adv_y0(i,j,k)=((ttt(i,j+1,k)+ttt(i,j,k))*v_sface(i,j,k) &
                       -(ttt(i,j,k)+ttt(i,j-1,k))*v_sface(i,j-1,k))*tarea_r(i,j,iblock)
            adv_xy1(i,j,k)=-dts*(ttt(i+1,j,k)-ttt(i,j,k))*2.0_r8*tarea_r(i,j,iblock)* &
                                  u_wface(i+1,j,k)*u_wface(i+1,j,k)/(htw(i+1,j,iblock)*hun(i+1,j,iblock))
            adv_xy2(i,j,k)= dts*(ttt(i,j,k)-ttt(i-1,j,k))*2.0_r8*tarea_r(i,j,iblock)*  &
                                  u_wface(i,j,k)*u_wface(i,j,k)/(htw(i,j,iblock)*hun(i,j,iblock))
            adv_xy3(i,j,k)=-dts*(ttt(i,j+1,k)-ttt(i,j,k))*2.0_r8*tarea_r(i,j,iblock)* &
                                  v_sface(i,j,k)*v_sface(i,j,k)/(hts(i,j,iblock)*hue(i,j,iblock))
            adv_xy4(i,j,k)= dts*(ttt(i,j,k)-ttt(i,j-1,k))*2.0_r8*tarea_r(i,j,iblock)* &
                                  v_sface(i,j-1,k)*v_sface(i,j-1,k)/(hts(i,j-1,iblock)*hue(i,j-1,iblock))
            adv_c1(i,j,k)=-ttt(i,j,k)*(u_wface(i+1,j,k)-u_wface(i,j,k))*tarea_r(i,j,iblock)*2.0_r8
            adv_c2(i,j,k)=-ttt(i,j,k)*(v_sface(i,j,k)-v_sface(i,j-1,k))*tarea_r(i,j,iblock)*2.0_r8


            if (k==1) then
              adv_za (i,j,k)=-0.5*ODZP(1)*WWW(I,J,2)*(ttt(I,J,2)+ttt(I,J,1))
              adv_zb1(i,j,k)=0
              adv_zb2(i,j,k)= 0.5*ODZP(1)*WWW(I,J,2)*WWW(I,J,2)*ODZT(2) &
                                *(ttt(I,J,1)-ttt(I,J,2))*dts
              adv_zc (i,j,k)=     ODZP(1)*ttt(I,J,1)*WWW(I,J,2)
            else if (k==km) then
              adv_za (i,j,k)= 0.5*ODZP(km)*WWW(I,J,km)*(ttt(I,J,km)+ttt(I,J,km-1))
              adv_zb1(i,j,k)=-0.5*ODZP(km)*WWW(I,J,km)*WWW(I,J,km)*ODZT(km  ) &
                                *(ttt(I,J,km-1)-ttt(I,J,km))*dts
              adv_zb2(i,j,k)=0
              adv_zc (i,j,k)=    -ODZP(km)*ttt(I,J,km)*WWW(I,J,km)
            else
              adv_za (i,j,k)= 0.5*ODZP(k)*WWW(I,J,k  )*(ttt(I,J,k)+ttt(I,J,k-1)) &
                         -0.5*ODZP(k)*WWW(I,J,k+1)*(ttt(I,J,k)+ttt(I,J,k+1))
              adv_zb1(i,j,k)=-0.5*ODZP(k)*WWW(I,J,k  )*WWW(I,J,k  )*ODZT(k  ) &
                                *(ttt(I,J,k-1)-ttt(I,J,k))*dts
              adv_zb2(i,j,k)= 0.5*ODZP(k)*WWW(I,J,k+1)*WWW(I,J,k+1)*ODZT(k+1) &
                                *(ttt(I,J,k)-ttt(I,J,k+1))*dts
              adv_zc (i,j,k)=    -ODZP(k)*ttt(I,J,k)*(WWW(I,J,k)-WWW(I,J,k+1))
            endif
!
            adv_xx(i,j,k)=-(adv_x0(i,j,k)+adv_xy1(i,j,k)+adv_xy2(i,j,k)+adv_c1(i,j,k))
            adv_yy(i,j,k)=-(adv_y0(i,j,k)+adv_xy3(i,j,k)+adv_xy4(i,j,k)+adv_c2(i,j,k))
            adv_zz(i,j,k)=-(adv_zb1(i,j,k)+adv_zb2(i,j,k)+adv_za(i,j,k)+adv_zc(i,j,k))
!           adv_xx(i,j,k)=-(adv_x0(i,j,k)+adv_xy1(i,j,k)+adv_xy2(i,j,k))
!           adv_yy(i,j,k)=-(adv_y0(i,j,k)+adv_xy3(i,j,k)+adv_xy4(i,j,k))
!           adv_zz(i,j,k)=-(adv_zb1(i,j,k)+adv_zb2(i,j,k)+adv_za(i,j,k))
            at00(i,j,k)=ttt(i,j,k)+(adv_xx(i,j,k)+adv_yy(i,j,k)+adv_zz(i,j,k))*dts
          ENDDO
        ENDDO
      ENDDO


      !az(:,:,:,mtracer,iblock)=  adv_xy4 !LPF20160812 wrong just for test
!
      WT1 = -1.0D10
      WT2 = +1.0D10
      DO K = 1,km
        DO J = 2, jmt-1
          DO I = 2,imt-1
           if ( k == 1 ) then
             atmax(i,j,k)=max(ttt(i,j,k)*vit(i,j,k,iblock)+(1.0D0-vit(i,j,k,iblock))*WT1, & 
                              ttt(i,j-1,k)*vit(i,j-1,k,iblock)+(1.0D0-vit(i,j-1,k,iblock))*WT1, & 
                              ttt(i,j+1,k)*vit(i,j+1,k,iblock)+(1.0D0-vit(i,j+1,k,iblock))*WT1, & 
                              ttt(i-1,j,k)*vit(i-1,j,k,iblock)+(1.0D0-vit(i-1,j,k,iblock))*WT1, &
                              ttt(i+1,j,k)*vit(i+1,j,k,iblock)+(1.0D0-vit(i+1,j,k,iblock))*WT1, &
                              ttt(i,j,k+1)*vit(i,j,k+1,iblock)+(1.0D0-vit(i,j,k+1,iblock))*WT1)
             atmin(i,j,k)=min(ttt(i,j,k)*vit(i,j,k,iblock)+(1.0D0-vit(i,j,k,iblock))*WT2, & 
                              ttt(i,j-1,k)*vit(i,j-1,k,iblock)+(1.0D0-vit(i,j-1,k,iblock))*WT2, & 
                              ttt(i,j+1,k)*vit(i,j+1,k,iblock)+(1.0D0-vit(i,j+1,k,iblock))*WT2, & 
                              ttt(i-1,j,k)*vit(i-1,j,k,iblock)+(1.0D0-vit(i-1,j,k,iblock))*WT2, &
                              ttt(i+1,j,k)*vit(i+1,j,k,iblock)+(1.0D0-vit(i+1,j,k,iblock))*WT2, &
                              ttt(i,j,k+1)*vit(i,j,k+1,iblock)+(1.0D0-vit(i,j,k+1,iblock))*WT2)
           else if ( k == km) then
             atmax(i,j,k)=max(ttt(i,j,k)*vit(i,j,k,iblock)+(1.0D0-vit(i,j,k,iblock))*WT1, & 
                              ttt(i,j-1,k)*vit(i,j-1,k,iblock)+(1.0D0-vit(i,j-1,k,iblock))*WT1, & 
                              ttt(i,j+1,k)*vit(i,j+1,k,iblock)+(1.0D0-vit(i,j+1,k,iblock))*WT1, & 
                              ttt(i-1,j,k)*vit(i-1,j,k,iblock)+(1.0D0-vit(i-1,j,k,iblock))*WT1, &
                              ttt(i+1,j,k)*vit(i+1,j,k,iblock)+(1.0D0-vit(i+1,j,k,iblock))*WT1, &
                              ttt(i,j,k-1)*vit(i,j,k-1,iblock)+(1.0D0-vit(i,j,k-1,iblock))*WT1)
             atmin(i,j,k)=min(ttt(i,j,k)*vit(i,j,k,iblock)+(1.0D0-vit(i,j,k,iblock))*WT2, & 
                              ttt(i,j-1,k)*vit(i,j-1,k,iblock)+(1.0D0-vit(i,j-1,k,iblock))*WT2, & 
                              ttt(i,j+1,k)*vit(i,j+1,k,iblock)+(1.0D0-vit(i,j+1,k,iblock))*WT2, & 
                              ttt(i-1,j,k)*vit(i-1,j,k,iblock)+(1.0D0-vit(i-1,j,k,iblock))*WT2, &
                              ttt(i+1,j,k)*vit(i+1,j,k,iblock)+(1.0D0-vit(i+1,j,k,iblock))*WT2, &
                              ttt(i,j,k-1)*vit(i,j,k-1,iblock)+(1.0D0-vit(i,j,k-1,iblock))*WT2 )
           else 
             atmax(i,j,k)=max(ttt(i,j,k)*vit(i,j,k,iblock)+(1.0D0-vit(i,j,k,iblock))*WT1, & 
                              ttt(i,j-1,k)*vit(i,j-1,k,iblock)+(1.0D0-vit(i,j-1,k,iblock))*WT1, & 
                              ttt(i,j+1,k)*vit(i,j+1,k,iblock)+(1.0D0-vit(i,j+1,k,iblock))*WT1, & 
                              ttt(i-1,j,k)*vit(i-1,j,k,iblock)+(1.0D0-vit(i-1,j,k,iblock))*WT1, &
                              ttt(i,j,k-1)*vit(i,j,k-1,iblock)+(1.0D0-vit(i,j,k-1,iblock))*WT1, &
                              ttt(i+1,j,k)*vit(i+1,j,k,iblock)+(1.0D0-vit(i+1,j,k,iblock))*WT1, &
                              ttt(i,j,k+1)*vit(i,j,k+1,iblock)+(1.0D0-vit(i,j,k+1,iblock))*WT1)
             atmin(i,j,k)=min(ttt(i,j,k)*vit(i,j,k,iblock)+(1.0D0-vit(i,j,k,iblock))*WT2, & 
                              ttt(i,j-1,k)*vit(i,j-1,k,iblock)+(1.0D0-vit(i,j-1,k,iblock))*WT2, & 
                              ttt(i,j+1,k)*vit(i,j+1,k,iblock)+(1.0D0-vit(i,j+1,k,iblock))*WT2, & 
                              ttt(i-1,j,k)*vit(i-1,j,k,iblock)+(1.0D0-vit(i-1,j,k,iblock))*WT2, &
                              ttt(i+1,j,k)*vit(i+1,j,k,iblock)+(1.0D0-vit(i+1,j,k,iblock))*WT2, &
                              ttt(i,j,k-1)*vit(i,j,k-1,iblock)+(1.0D0-vit(i,j,k-1,iblock))*WT2, &
                              ttt(i,j,k+1)*vit(i,j,k+1,iblock)+(1.0D0-vit(i,j,k+1,iblock))*WT2)
           end if
          ENDDO
        ENDDO
      ENDDO
   
      DO K = 1,km
        DO J = 2, jmt-1
          DO I = 2, imt-1
            if (at00(i,j,k)>atmax(i,j,k).or.at00(i,j,k)<atmin(i,j,k)) then
              adv_xy1(i,j,k)=-(ttt(i+1,j,k)-ttt(i,j,k))*        &
                                 abs(u_wface(i+1,j,k))*tarea_r(i,j,iblock)
              adv_xy2(i,j,k)= (ttt(i,j,k)-ttt(i-1,j,k))*        &
                                 abs(u_wface(i,j,k))*tarea_r(i,j,iblock)
              adv_xy3(i,j,k)=-(ttt(i,j+1,k)-ttt(i,j,k))*        &
                                 abs(v_sface(i,j,k))*tarea_r(i,j,iblock)
              adv_xy4(i,j,k)= (ttt(i,j,k)-ttt(i,j-1,k))*        &
                                 abs(v_sface(i,j-1,k))*tarea_r(i,j,iblock)
              adv_xy2(i+1,j,k)= (ttt(i+1,j,k)-ttt(i,j,k))*        &
                                 abs(u_wface(i+1,j,k))*tarea_r(i+1,j,iblock)
              adv_xy1(i-1,j,k)=-(ttt(i,j,k)-ttt(i-1,j,k))*        &
                                 abs(u_wface(i,j,k))*tarea_r(i-1,j,iblock)
              adv_xy4(i,j+1,k)= (ttt(i,j+1,k)-ttt(i,j,k))*        &
                                 abs(v_sface(i,j,k))*tarea_r(i,j+1,iblock)
              adv_xy3(i,j-1,k)=-(ttt(i,j,k)-ttt(i,j-1,k))*        &
                                 abs(v_sface(i,j-1,k))*tarea_r(i,j-1,iblock)
              if ( k == 1 ) then
                 adv_zb2(i,j,k)= 0.5*abs(WWW(I,J,k+1))*ODZP(k)*(ttt(I,J,k)-ttt(I,J,k+1))
                 adv_zb1(i,j,k+1)=-0.5*abs(WWW(I,J,k+1))*ODZP(k+1)*(ttt(I,J,k)-ttt(I,J,k+1))
              else if ( k==km) then
                 adv_zb2(i,j,k-1)= 0.5*abs(WWW(I,J,k))*ODZP(k-1)*(ttt(I,J,k-1)-ttt(I,J,k))
                 adv_zb1(i,j,k)=-0.5*abs(WWW(I,J,k  ))*ODZP(k)*(ttt(I,J,k-1)-ttt(I,J,k))
              else
                 adv_zb2(i,j,k-1)= 0.5*abs(WWW(I,J,k))*ODZP(k-1)*(ttt(I,J,k-1)-ttt(I,J,k))
                 adv_zb1(i,j,k)=-0.5*abs(WWW(I,J,k  ))*ODZP(k)*(ttt(I,J,k-1)-ttt(I,J,k))
                 adv_zb2(i,j,k)= 0.5*abs(WWW(I,J,k+1))*ODZP(k)*(ttt(I,J,k)-ttt(I,J,k+1))
                 adv_zb1(i,j,k+1)=-0.5*abs(WWW(I,J,k+1))*ODZP(k+1)*(ttt(I,J,k)-ttt(I,J,k+1))
              end if
            endif
          END DO
        END DO
      END DO
   
      adv_xx=-(adv_x0+adv_xy1+adv_xy2+adv_c1)
      adv_yy=-(adv_y0+adv_xy3+adv_xy4+adv_c2)
      adv_zz=-(adv_za+adv_zb1+adv_zb2+adv_zc)
      adv_tt= adv_xx+adv_yy+adv_zz
#if (defined LOWRES)
!LPF20160829
      ax(:,:,:,mtracer,iblock)= adv_xx(:,:,:) !/dble(NSS)
      ay(:,:,:,mtracer,iblock)= adv_yy(:,:,:) !/dble(NSS)
      az(:,:,:,mtracer,iblock)= adv_zz(:,:,:) !/dble(NSS)
!      ax(:,:,:,mtracer,iblock)= ax(:,:,:,mtracer,iblock) + adv_xx(:,:,:) !/dble(NSS)
!      ay(:,:,:,mtracer,iblock)= ay(:,:,:,mtracer,iblock) + adv_yy(:,:,:) !/dble(NSS)
!      az(:,:,:,mtracer,iblock)= az(:,:,:,mtracer,iblock) + adv_zz(:,:,:) !/dble(NSS)
!      ax(:,:,:,mtracer,iblock)= ax(:,:,:,mtracer,iblock) + adv_xx(:,:,:)/dble(NSS)
!      ay(:,:,:,mtracer,iblock)= ay(:,:,:,mtracer,iblock) + adv_yy(:,:,:)/dble(NSS)
!      az(:,:,:,mtracer,iblock)= az(:,:,:,mtracer,iblock) + adv_zz(:,:,:)/dble(NSS)
!LPF20160829
#endif
      deallocate ( adv_c1,adv_c2,adv_zc)
      deallocate ( adv_x0, adv_y0, adv_xx, adv_yy)
      deallocate ( adv_xy1, adv_xy2, adv_xy3, adv_xy4)
      deallocate ( at00, atmax, atmin )
      deallocate (adv_zz, adv_za, adv_zb1 )
      deallocate (adv_zb2, atmaxz, atminz, atz)
!
!WR20241015

      else if ( adv_tracer(1:3) == 'prm' ) then

        allocate ( adv_xx(imt,jmt,km), adv_yy(imt,jmt,km), adv_zz(imt,jmt,km) )

        call ADVEC_PRM_TS( uuu, vvv, www, ttt, adv_tt, iblock, adv_xx, adv_yy, adv_zz )

        ax(:,:,:,mtracer,iblock)= adv_xx(:,:,:)
        ay(:,:,:,mtracer,iblock)= adv_yy(:,:,:)
        az(:,:,:,mtracer,iblock)= adv_zz(:,:,:)

        deallocate ( adv_xx, adv_yy, adv_zz )
!        
      else
          write(*,*) "adv_tracer666 =", adv_tracer
        call exit_licom(sigAbort,'The false advection option for tracer')
      end if
!

      end subroutine advection_tracer

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      SUBROUTINE ADVEC_PRM_TS( UPRM, VPRM, WPRM, VARTMP, DLTMP, iblock, AXTMP, AYTMP, AZTMP )
#include <def-undef.h>
use   precision_mod
use   param_mod,  only: IMT, JMT, KM, KMP1, I, J, K
use   pconst_mod, only: ZKT, ZKP, DTS
use   grid, only: DXT, DYT, DXTR, DYTR, TLAT, TLON
use   constant_mod, only : DEGtoRAD

      IMPLICIT NONE

      REAL(r8) :: UPRM(IMT,JMT,KM)
      REAL(r8) :: VPRM(IMT,JMT,KM)
      REAL(r8) :: WPRM(IMT,JMT,KMP1)
      REAL(r8) :: VARTMP(IMT,JMT,KM)
      REAL(r8) :: DLTMP(IMT,JMT,KM)

      integer, intent(in) :: iblock

      REAL(r8) :: AXTMP(IMT,JMT,KM)
      REAL(r8) :: AYTMP(IMT,JMT,KM)
      REAL(r8) :: AZTMP(IMT,JMT,KM)

      REAL(r8), ALLOCATABLE :: VQC(:)
      REAL(r8), ALLOCATABLE :: DXYZ(:)
      INTEGER :: INDXYZ

      REAL(r8), ALLOCATABLE :: UVWhaf(:)
      REAL(r8), ALLOCATABLE :: flxqR(:), flxqL(:)

      REAL(r8) :: CMIU
      REAL(r8) :: uhafR, uhafL

      DO K = 1, KM
         DO J = 1, JMT
            DO I = 1, IMT
               DLTMP(I,J,K) = 0.0D0
               AXTMP(I,J,K) = 0.0D0
               AYTMP(I,J,K) = 0.0D0
               AZTMP(I,J,K) = 0.0D0
            ENDDO
         ENDDO
      ENDDO

      DO J = 2, JMT

         DO K = 1, KM

            IF ( ALLOCATED( VQC ) ) DEALLOCATE( VQC )
            ALLOCATE( VQC(IMT) )

            IF ( ALLOCATED( DXYZ ) ) DEALLOCATE( DXYZ )
            ALLOCATE( DXYZ(IMT) )

            IF ( ALLOCATED( UVWhaf ) ) DEALLOCATE( UVWhaf )
            ALLOCATE( UVWhaf(IMT) )

            IF ( ALLOCATED( flxqR ) ) DEALLOCATE( flxqR )
            ALLOCATE( flxqR(IMT) )

            IF ( ALLOCATED( flxqL ) ) DEALLOCATE( flxqL )
            ALLOCATE( flxqL(IMT) )

            DO I = 1, IMT
               VQC(I) = VARTMP(I,J,K)
               DXYZ(I) = DXT(I,J,iblock)
            ENDDO

            DO I = 1, IMT-1
               UVWhaf(I) = 0.50D0*(UPRM(I+1,J-1,K)+UPRM(I+1,J,K))
            ENDDO
            UVWhaf(IMT) = 0.0D0

            INDXYZ = 1

            CALL QRLFLX_PRM_NQGRID( INDXYZ, IMT, VQC, DXYZ, UVWhaf, DTS, flxqR, flxqL)

            DO I = 1, IMT

               IF ( I.EQ.1 ) THEN

                       uhafR = UVWhaf(I)
                       uhafL = 0.0D0

               ELSE IF ( I.GE.2 .AND. I.LE.(IMT-1) ) THEN

                       uhafR = UVWhaf(I)
                       uhafL = UVWhaf(I-1)

               ELSE IF ( I.EQ.IMT ) THEN

                       uhafR = 0.0D0
                       uhafL = UVWhaf(I-1)

               ENDIF

                       CMIU = DXTR(I,J,iblock)
                       AXTMP(I,J,K) = ( CMIU*( flxqL(I)-flxqR(I) ) + VARTMP(I,J,K)*(uhafR-uhafL)*DTS*CMIU )/DTS

            ENDDO

         ENDDO

      ENDDO

      DO I = 1, IMT-1

         DO K = 1, KM

            IF ( ALLOCATED( VQC ) ) DEALLOCATE( VQC )
            ALLOCATE( VQC(JMT) )

            IF ( ALLOCATED( DXYZ ) ) DEALLOCATE( DXYZ )
            ALLOCATE( DXYZ(JMT) )

            IF ( ALLOCATED( UVWhaf ) ) DEALLOCATE( UVWhaf )
            ALLOCATE( UVWhaf(JMT) )

            IF ( ALLOCATED( flxqR ) ) DEALLOCATE( flxqR )
            ALLOCATE( flxqR(JMT) )

            IF ( ALLOCATED( flxqL ) ) DEALLOCATE( flxqL )
            ALLOCATE( flxqL(JMT) )

            DO J = 1, JMT
               VQC(J) =  VARTMP(I,J,K)
               DXYZ(J) = DYT(I,J,iblock)
               UVWhaf(J) = 0.50D0*(VPRM(I,J,K)+VPRM(I+1,J,K))
            ENDDO

            INDXYZ = 2

            CALL QRLFLX_PRM_NQGRID( INDXYZ, JMT, VQC, DXYZ, UVWhaf, DTS, flxqR, flxqL)

            DO J = 1, JMT

               IF ( J.EQ.1 ) THEN

                       uhafR = UVWhaf(J)
                       uhafL = 0.0D0

               ELSE IF ( J.GE.2 .AND. J.LE.(JMT-1) ) THEN

                       uhafR = UVWhaf(J)
                       uhafL = UVWhaf(J-1)

               ELSE IF ( J.EQ.JMT ) THEN

                       uhafR = 0.0D0
                       uhafL = UVWhaf(J-1)

               ENDIF

                       CMIU = DYTR(I,J,iblock)
                       AYTMP(I,J,K) = ( CMIU*( flxqL(J)-flxqR(J) ) + VARTMP(I,J,K)*(uhafR-uhafL)*DTS*CMIU )/DTS

            ENDDO

         ENDDO

      ENDDO

      DO I = 1, IMT

         DO J = 1, JMT

            IF ( ALLOCATED( VQC ) ) DEALLOCATE( VQC )
            ALLOCATE( VQC(KM) )

            IF ( ALLOCATED( DXYZ ) ) DEALLOCATE( DXYZ )
            ALLOCATE( DXYZ(KM) )

            IF ( ALLOCATED( UVWhaf ) ) DEALLOCATE( UVWhaf )
            ALLOCATE( UVWhaf(KM) )

            IF ( ALLOCATED( flxqR ) ) DEALLOCATE( flxqR )
            ALLOCATE( flxqR(KM) )

            IF ( ALLOCATED( flxqL ) ) DEALLOCATE( flxqL )
            ALLOCATE( flxqL(KM) )

            DO K = 1, KM
               VQC(K) = VARTMP(I,J,K)
               DXYZ(K) = ZKP(K)-ZKP(K+1)
            ENDDO

            DO K = 1, KM-1
               UVWhaf(K) = -WPRM(I,J,K+1)
            ENDDO
            UVWhaf(KM) = 0.0D0

            INDXYZ = 3

            CALL QRLFLX_PRM_NQGRID( INDXYZ, KM, VQC, DXYZ, UVWhaf, DTS, flxqR, flxqL)

            DO K = 1, KM

               CMIU = 1.0D0/DXYZ(K)

               IF (K .EQ. 1) THEN
                       uhafR = UVWhaf(K)
                       uhafL = 0.0D0
               ELSE
                       uhafR = UVWhaf(K)
                       uhafL = UVWhaf(K-1)
               ENDIF

               AZTMP(I,J,K) = (CMIU*( flxqL(K)-flxqR(K)) + VARTMP(I,J,K)*(uhafR-uhafL)*DTS*CMIU)/DTS

            ENDDO

         ENDDO

      ENDDO


      DO K = 1, KM
         DO J = 1, JMT
            DO I = 1, IMT

               DLTMP(I,J,K) = AZTMP(I,J,K)+AXTMP(I,J,K)+AYTMP(I,J,K)

            ENDDO
         ENDDO
      ENDDO

      DEALLOCATE( VQC )
      DEALLOCATE( DXYZ )
      DEALLOCATE( UVWhaf )

      DEALLOCATE( flxqL )
      DEALLOCATE( flxqR )

      END SUBROUTINE ADVEC_PRM_TS
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
      SUBROUTINE QRLFLX_PRM_NQGRID( KXYZ, IXYZ, VQ, DltXYZ, Suvw, Dtsc, flxR, flxL)

       use   precision_mod

       IMPLICIT NONE

       INTEGER :: KXYZ
       INTEGER :: IXYZ
       REAL(r8) :: VQ(IXYZ), DltXYZ(IXYZ)
       REAL(r8) :: Suvw(IXYZ)
       REAL(r8) :: Dtsc
       REAL(r8) :: flxR(IXYZ), flxL(IXYZ)

       REAL(r8), ALLOCATABLE :: dltQ(:)
       REAL(r8), ALLOCATABLE :: qR(:), qL(:)

       REAL(r8) :: COF

       REAL(r8) :: deltX1, deltX2, deltX3, deltX4
       REAL(r8) :: deltX5, deltX6, deltX7

       REAL(r8) :: dltc, dltm, dltp
       REAL(r8) :: dltca, dltma, dltpa

       REAL(r8) :: VQmin, VQmax

       INTEGER ::  N

       REAL(r8) :: uhafR, uhafL

       REAL(r8) :: FMIU
       REAL(r8) :: Atmp, Btmp
       REAL(r8) :: Ctmp, Dtmp
       REAL(r8) :: Gama, Beta

       REAL(r8) :: Fmin, Fmax, Fave

       REAL(r8) :: EPStmp

       IF ( ALLOCATED( dltQ ) ) DEALLOCATE( dltQ )
       ALLOCATE( dltQ(IXYZ) )

       IF ( ALLOCATED( qR ) ) DEALLOCATE( qR )
       ALLOCATE( qR(IXYZ) )

       IF ( ALLOCATED( qL ) ) DEALLOCATE( qL )
       ALLOCATE( qL(IXYZ) )

       COF = 3.0D0
       EPStmp = 1.0D-20

       DO N = 2, IXYZ-1

          deltX1 = DltXYZ(N)/( DltXYZ(N-1)+DltXYZ(N)+DltXYZ(N+1) )
          deltX2 = ( 2.0D0*DltXYZ(N-1)+DltXYZ(N) )/( DltXYZ(N+1)+DltXYZ(N) )
          deltX3 = ( DltXYZ(N)+2.0D0*DltXYZ(N+1) )/( DltXYZ(N-1)+DltXYZ(N) )
          dltc = deltX1*(deltX2*(VQ(N+1)-VQ(N))+deltX3*(VQ(N)-VQ(N-1)))
          dltm = VQ(N)-VQ(N-1)
          dltp = VQ(N+1)-VQ(N)

          IF ( (dltm*dltp) .LE. 0.0 ) THEN
               dltQ(N) = 0.0
          ELSE
               dltma = COF*DABS(dltm)
               dltpa = COF*DABS(dltp)
               dltca = DABS(dltc)
               dltQ(N) = DMIN1( dltca, dltma, dltpa ) * DSIGN(1.0D0,dltc)
          ENDIF

       ENDDO
       dltQ(1) = 0.0D0
       dltQ(IXYZ) = 0.0D0

!             Compute the values at the left and right edges of cells

       DO N = 2, IXYZ-2
          deltX1 = DltXYZ(N)/( DltXYZ(N)+DltXYZ(N+1) )
          deltX2 = 1.0D0/( DltXYZ(N-1)+DltXYZ(N)+DltXYZ(N+1)+DltXYZ(N+2) )
          deltX3 = ( 2.0D0*DltXYZ(N+1)*DltXYZ(N) )/( DltXYZ(N)+DltXYZ(N+1) )
          deltX4 = ( DltXYZ(N-1)+DltXYZ(N) )/( 2.0D0*DltXYZ(N)+DltXYZ(N+1) )
          deltX5 = ( DltXYZ(N+2)+DltXYZ(N+1) )/( 2.0D0*DltXYZ(N+1)+DltXYZ(N) )
          deltX6 = DltXYZ(N)*( DltXYZ(N-1)+DltXYZ(N) )/(2.0D0*DltXYZ(N)+DltXYZ(N+1))
          deltX7 = DltXYZ(N+1)*( DltXYZ(N+1)+DltXYZ(N+2) )/(DltXYZ(N)+2.0D0*DltXYZ(N+1) )

          qR(N) = VQ(N)+deltX1*(VQ(N+1)-VQ(N))+deltX2*(deltX3*(deltX4-deltX5)*(VQ(N+1)-VQ(N))-&
                  deltX6*dltQ(N+1)+deltX7*dltQ(N) )

          VQmin = DMIN1( VQ(N-1), VQ(N), VQ(N+1) )
          VQmax = DMAX1( VQ(N-1), VQ(N), VQ(N+1) )
          qR(N) = DMAX1( VQmin, qR(N) )
          qR(N) = DMIN1( VQmax, qR(N) )
       ENDDO
       qR(1) = 0.5D0*(VQ(1)+VQ(2))
       qR(IXYZ) = VQ(IXYZ)
       qR(IXYZ-1) = 0.5D0*(VQ(IXYZ-1)+VQ(IXYZ))

       DO N = 2, IXYZ
          qL(N) = qR(N-1)
       ENDDO
       qL(1) = VQ(1)

!            To Compute the Flux Through the Cell

       DO N = 1, IXYZ

          IF ( N.EQ.1 ) THEN

                          uhafR = Suvw(N)
                          uhafL = 0.0D0

                  IF( uhafR .GT. 0.0D0 ) THEN
                          

                          FMIU = DABS(uhafR*Dtsc)

                          Atmp = qR(N)

                          Fmin = DMIN1( qR(N), qL(N) )
                          Fmax = DMAX1( qR(N), qL(N) )
                          Fave = DMAX1( DMIN1( VQ(N), Fmax ), Fmin )
                          Ctmp = DABS( qR(N)-Fave )
                          Dtmp = DABS( Fave-qL(N) )
                          Gama = (Ctmp+EPStmp)/(Dtmp+EPStmp)
                          Beta = ( Gama-1.0D0 )/(-DltXYZ(N))
                          Btmp = ( Gama*Fave-qR(N) )/(-DltXYZ(N))

                          flxR(N) = ( Atmp*FMIU - Btmp*FMIU*FMIU )/( 1.0D0-Beta*FMIU )

                  ELSE IF ( uhafR.EQ.0.0D0 ) THEN

                          flxR(N) = 0.0D0

                  ELSE

                          FMIU = DABS(-uhafR*Dtsc)

                          Atmp = qR(N)

                          Fmin = DMIN1( qR(N+1), qL(N+1) )
                          Fmax = DMAX1( qR(N+1), qL(N+1) )
                          Fave = DMAX1( DMIN1( VQ(N+1), Fmax ), Fmin )
                          Ctmp = DABS( qR(N)-Fave )
                          Dtmp = DABS( Fave-qR(N+1) )
                          Gama = (Ctmp+EPStmp)/(Dtmp+EPStmp)
                          Beta = ( Gama-1.0D0 )/DltXYZ(N+1)
                          Btmp = ( Gama*Fave-qR(N) )/DltXYZ(N+1)

                          flxR(N) = -1.0D0*( Atmp*FMIU + Btmp*FMIU*FMIU )/ ( 1.0D0+Beta*FMIU )
                  ENDIF


                          flxL(N) = 0.0D0


          ELSE IF ( N.GE.2 .AND. N.LE.(IXYZ-1) ) THEN

                          uhafR = Suvw(N)
                          uhafL = Suvw(N-1)
                       
                  IF( uhafR .GT. 0.0D0 ) THEN

                          FMIU = DABS(uhafR*Dtsc) 

                          Atmp = qR(N)

                          Fmin = DMIN1( qR(N), qL(N) )
                          Fmax = DMAX1( qR(N), qL(N) )
                          Fave = DMAX1( DMIN1( VQ(N), Fmax ), Fmin )
                          Ctmp = DABS( qR(N)-Fave )
                          Dtmp = DABS( Fave-qL(N) )
                          Gama = (Ctmp+EPStmp)/(Dtmp+EPStmp)
                          Beta = ( Gama-1.0D0)/(-DltXYZ(N))
                          Btmp = ( Gama*Fave-qR(N) )/(-DltXYZ(N))

                          flxR(N) = ( Atmp*FMIU - Btmp*FMIU*FMIU )/( 1.0D0-Beta*FMIU)

                  ELSE IF ( uhafR.EQ.0.0D0 ) THEN

                          flxR(N) = 0.0D0

                  ELSE

                          FMIU = DABS(-uhafR*Dtsc)

                          Atmp = qR(N)

                          Fmin = DMIN1( qR(N+1), qL(N+1) )
                          Fmax = DMAX1( qR(N+1), qL(N+1) )
                          Fave = DMAX1( DMIN1( VQ(N+1), Fmax ), Fmin )
                          Ctmp = DABS( qR(N)-Fave )
                          Dtmp = DABS( Fave-qR(N+1) )
                          Gama = (Ctmp+EPStmp)/(Dtmp+EPStmp)
                          Beta = ( Gama-1.0D0 )/DltXYZ(N+1)
                          Btmp = ( Gama*Fave-qR(N) )/DltXYZ(N+1)

                          flxR(N) = -1.0D0*( Atmp*FMIU + Btmp*FMIU*FMIU )/ ( 1.0D0+Beta*FMIU )

                  ENDIF

                  IF( uhafL .GT. 0.0D0 ) THEN

                          FMIU = DABS(uhafL*Dtsc)

                          Atmp = qL(N)

                          Fmin = DMIN1( qL(N-1), qR(N-1) )
                          Fmax = DMAX1( qL(N-1), qR(N-1) )
                          Fave = DMAX1( DMIN1( VQ(N-1), Fmax ), Fmin )
                          Ctmp = DABS( qL(N)-Fave )
                          Dtmp = DABS( Fave-qL(N-1) )
                          Gama = (Ctmp+EPStmp)/(Dtmp+EPStmp)
                          Beta = ( Gama-1.0D0 )/(-DltXYZ(N-1))
                          Btmp = ( Gama*Fave-qL(N) )/(-DltXYZ(N-1))

                          flxL(N) = ( Atmp*FMIU-Btmp*FMIU*FMIU )/( 1.0D0-Beta*FMIU )

                  ELSE IF ( uhafL.EQ.0.0D0 ) THEN

                          flxL(N) = 0.0D0

                  ELSE

                          FMIU = DABS(-uhafL*Dtsc) 

                          Atmp = qL(N)

                          Fmin = DMIN1( qL(N), qR(N) )
                          Fmax = DMAX1( qL(N), qR(N) )
                          Fave = DMAX1( DMIN1( VQ(N), Fmax ), Fmin )
                          Ctmp = DABS( qL(N)-Fave )
                          Dtmp = DABS( Fave-qR(N) )
                          Gama = (Ctmp+EPStmp)/(Dtmp+EPStmp)
                          Beta = ( Gama-1.0D0 )/DltXYZ(N)
                          Btmp = ( Gama*Fave-qL(N) )/DltXYZ(N)

                          flxL(N) = -1.0D0*( Atmp*FMIU+Btmp*FMIU*FMIU )/( 1.0D0+Beta*FMIU )

                  ENDIF

          ELSE IF ( N.EQ.IXYZ ) THEN

                          uhafR = 0.0D0
                          uhafL = Suvw(N-1)

                          flxR(N) = 0.0D0

                  IF( uhafL .GT. 0.0D0 ) THEN

                          FMIU = DABS(uhafL*Dtsc)

                          Atmp = qL(N)

                          Fmin = DMIN1( qL(N-1), qR(N-1) )
                          Fmax = DMAX1( qL(N-1), qR(N-1) )
                          Fave = DMAX1( DMIN1( VQ(N-1), Fmax ), Fmin )
                          Ctmp = DABS( qL(N)-Fave )
                          Dtmp = DABS( Fave-qL(N-1) )
                          Gama = (Ctmp+EPStmp)/(Dtmp+EPStmp)
                          Beta = ( Gama-1.0D0 )/(-DltXYZ(N-1))
                          Btmp = ( Gama*Fave-qL(N) )/(-DltXYZ(N-1))

                          flxL(N) = ( Atmp*FMIU-Btmp*FMIU*FMIU )/( 1.0D0-Beta*FMIU )

                  ELSE IF ( uhafL.EQ.0.0D0 ) THEN

                          flxL(N) = 0.0D0

                  ELSE

                          FMIU = DABS(-uhafL*Dtsc)

                          Atmp = qL(N)

                          Fmin = DMIN1( qL(N), qR(N) )
                          Fmax = DMAX1( qL(N), qR(N) )
                          Fave = DMAX1( DMIN1( VQ(N), Fmax ), Fmin )
                          Ctmp = DABS( qL(N)-Fave )
                          Dtmp = DABS( Fave-qR(N) )
                          Gama = (Ctmp+EPStmp)/(Dtmp+EPStmp)
                          Beta = ( Gama-1.0D0 )/DltXYZ(N)
                          Btmp = ( Gama*Fave-qL(N) )/DltXYZ(N)

                          flxL(N) = -1.0D0*( Atmp*FMIU+Btmp*FMIU*FMIU )/( 1.0D0+Beta*FMIU )

                  ENDIF

          ENDIF

       ENDDO

       DEALLOCATE( dltQ )
       DEALLOCATE( qR )
       DEALLOCATE( qL )

       RETURN

      END SUBROUTINE QRLFLX_PRM_NQGRID


 end module advection
!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
