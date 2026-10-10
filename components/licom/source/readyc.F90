!     =================
      SUBROUTINE READYC
!     =================
!     PREPARATION OF BAROTROPIC AND BAROCLINIC INTEGRATION
 
!     ADVECTION + DIFFUSION
 
#include <def-undef.h>
#ifdef SW_ATHREAD
use sunway_mod
use sunway_init
#endif
use precision_mod
use param_mod
use pconst_mod
use dyn_mod
use work_mod
use tracer_mod
use pmix_mod
use POP_HaloMod
use POP_GridHorzMod
#if ( defined TIDEMIX )
#if ( defined TIDEFORCEGW )
use forc_mod, only: su,sv,USTAR,BUOYTUR, BUOYSOL,wave_dis,depthstd
#else
use forc_mod, only: su,sv,USTAR,BUOYTUR, BUOYSOL,wave_dis
#endif
#else
use forc_mod, only: su,sv,USTAR,BUOYTUR, BUOYSOL
#endif
use domain
use grid
use blocks
use advection
use operators
use LICOM_Error_mod
#ifdef BIHAR
use hmix_del4
#else
use hmix_del2
#endif
use msg_mod 
use gather_scatter
use distribution
use constant_mod
use canuto_2010_mod
!use canuto_mod
use var_check
      IMPLICIT NONE
!      REAL(r8)  :: WKP (KMP1)
      INTEGER   :: IWK,n2, iblock,kmb
      REAL(r8)  :: WK1 (KM) ,WK2 (KM), WK3 (KM),WK4(KM)
      REAL(r8)  :: WP1 (KM) ,WP2 (KM), WP3 (KM)
      REAL(r8)  :: WP4 (KMM1) ,WP5 (KMM1), WP6 (KMM1)
      REAL(r8)  :: WP7 (KMM1) ,WP8(KM)
      REAL(r8)  :: WP9 ,WP10, WP11
     ! REAL(r8),dimension(IMT,JMT,KM,MAX_BLOCKS_CLINIC),target :: WP12,WP13,riv1,riv2,wk1_tmp,wk2_tmp,wk3_tmp
      REAL(r8)  :: epsln,RKV,RKV1,ek0
      REAL(r8)  :: adv_x1,adv_x2,adv_y1,adv_y2,diff_u1,diff_u2,diff_v1,diff_v2
      REAL(r8)  :: dlux,dlvx,dluy,dlvy,dluz,dlvz,adv_z1,adv_z2,adv_z3,adv_z4
      real(r8)  :: hdvk(imt,jmt), hduk(imt,jmt), adv_uu(imt,jmt,km), adv_vv(imt,jmt,km)

!LPF20160505
      real(r8)  :: tq,sq,ak_tide_mixing(km)
#ifdef BCKMEX
      real(r8),target :: diff_back(imt,jmt,max_blocks_clinic),&
                  diff_back_sh(imt,jmt,max_blocks_clinic),&
                  diff_back_nh(imt,jmt,max_blocks_clinic) 
#endif
      REAL(r8)    :: DENS
      EXTERNAL DENS

!
!     real (r8) :: ttt(imt_global, jmt_global)
      type (block) :: this_block          ! block information for current block
      integer ::  ErrorCode

    
#ifdef CANUTO2010
!#if (defined CANUTO2010)
      REAL(r8)  :: AIDIF
#endif

        real (r8):: akt_back(km-1),aks_back(km-1),akm_back(km-1)


!YU 
      !real (r8):: akt_back(2*km),aks_back(2*km),akm_back(2*km) !LPF20160728
!Yu   if (ierr /= 0) then
!Yu      write(6,*)'allocation error---riu'
!Yu      stop
!Yu   end if
      epsln = 1.0D-25
!lhl0711
#if (defined CANUTO2010)
       AIDIF=0.5
!      if (ISC/=0)  AIDIF=0.5  
#endif
!lhl0711
!      s2t  = c0
!      ridt = c0
!      riu  = c0
!      wp12 = c0
!      wp13 = c0
!      wk1_tmp = c0
!      wk2_tmp = c0
!      wk3_tmp = c0
#ifdef SW_ATHREAD
#ifdef ENABLE_GPLT
    call t_startf("readyc1"//CHAR(0))
#endif
    call athread_spawn(slave_readyc1_cu,readyc1_param)
    call athread_join()
#ifdef ENABLE_GPLT
    call t_stopf("readyc1"//CHAR(0))
#endif
#else
!        ts1 = clock0f() 
!!$OMP PARALLEL DO PRIVATE (IBLOCK,J,I) COLLAPSE(3)
   DO IBLOCK = 1, NBLOCKS_CLINIC
      DO J = JST,JET
         DO I = 1,IMT
            H0BL (I,J,IBLOCK)= H0BF (I,J,IBLOCK)
            H0BF (I,J,IBLOCK)= H0 (I,J,IBLOCK)
         END DO
      END DO
   END DO
 

 
!---------------------------------------------------------------------
!     Calculating Richardson number riu (at U/V-point);
!---------------------------------------------------------------------

!!$OMP PARALLEL DO PRIVATE (IBLOCK,K) COLLAPSE(2)
   DO IBLOCK = 1, NBLOCKS_CLINIC
      DO K = 1,KM
         call ugrid_to_tgrid(wp12(:,:,k,iblock),up(:,:,k,iblock),iblock,k)
         call ugrid_to_tgrid(wp13(:,:,k,iblock),vp(:,:,k,iblock),iblock,k)
      END DO
   END DO
!    te1 = clock0f()
!    if (my_task == master_task) then
!        write(*,*) "readyc1 mpe spend (te1-ts1)=",(te1-ts1),' s'
!    endif

!
#endif 
#ifdef SW_ATHREAD
#ifdef ENABLE_GPLT
    call t_startf("readyc1_2"//CHAR(0))
#endif
    call athread_spawn(slave_readyc1_2,readyc1_2_param)
    call athread_join()
#ifdef ENABLE_GPLT
    call t_stopf("readyc1_2"//CHAR(0))
#endif
#else
!!$OMP PARALLEL DO PRIVATE (IBLOCK) COLLAPSE(4)
   DO IBLOCK = 1, NBLOCKS_CLINIC
      DO K = 1,KMM1
         DO J = 2, JMT
            DO I = 1,IMT-1
               riv1(i,j,k,iblock) = wp12 (I,J,K,iblock)*vit(i,j,k,iblock) - wp12 (I,J,K +1,iblock)*vit(i,j,k+1,iblock)
               riv2(i,j,k,iblock) = wp13 (I,J,K,iblock)*vit(i,j,k,iblock) - wp13 (I,J,K +1,iblock)*vit(i,j,k+1,iblock)
               s2t (i,j,k,iblock) =vit(i,j,k+1,iblock)*(riv1(i,j,k,iblock)*riv1(i,j,k,iblock)+riv2(i,j,k,iblock)*riv2(i,j,k,iblock))*ODZT(K+1)*ODZT(K+1)
#ifdef CANUTO2010  
               rit (i,j,k,iblock)= VIT (I,J,K +1,iblock)*rict(i,j,k,iblock)/(s2t(i,j,k,iblock)+epsln)
#else          
               rit (i,j,k,iblock)= rit (i,j,k,iblock,iblock) +VIT (I,J,K +1,iblock,iblock)* &
                                   rict(i,j,k,iblock,iblock)/(s2t(i,j,k,iblock,iblock)+epsln)
#endif
            END DO
         END DO
      END DO
   END DO
!!$OMP PARALLEL DO PRIVATE (IBLOCK) COLLAPSE(4)
   DO IBLOCK = 1, NBLOCKS_CLINIC
      DO K = 1,KMM1
         DO J = 1, JMT
            DO I = 1, IMT
               riv1(i,j,k,iblock) = UP (I,J,K,iblock) - UP (I,J,K +1,iblock)
               riv2(i,j,k,iblock) = VP (I,J,K,iblock) - VP (I,J,K +1,iblock)
               s2u (i,j,k,iblock) =viv(i,j,k+1,iblock)*(riv1(i,j,k,iblock)*riv1(i,j,k,iblock)+riv2(i,j,k,iblock)*riv2(i,j,k,iblock))*ODZT(K+1)*ODZT(K+1)
! calculate the shear square and the T component minus S component of Richardson Number
!lhl1204
               riu (i,j,k,iblock) = VIV (I,J,K +1,iblock)*ric (i,j,k,iblock)/(s2u(i,j,k,iblock)+epsln)
!lhl1204
            END DO
         END DO
      END DO
   END DO
!
!
#endif
!LPF20160505
#ifdef BCKMEX
#ifdef SW_ATHREAD
!add the vertical diffusivity due to internal wave mixing based on Jochum(2009)
!!$OMP PARALLEL
!!$OMP WORKSHARE
    call athread_spawn(slave_readyc1_3,readyc1_3_param)
    call athread_join()

#else
             diff_back=0.0d0
             diff_back_sh=0.0d0
             diff_back_nh=0.0d0
!!$OMP END WORKSHARE
!!$OMP END PARALLEL
!!!$OMP PARALLEL DO PRIVATE(iblock,j,i)
   DO iblock = 1, nblocks_clinic
     DO J = 1, JMT
         DO I = 1, IMT
          if(tlat(i,j,iblock)<0.0) then
           diff_back_sh(i,j,iblock) =diff_back_coef_max*exp(-(0.4d0*(tlat(i,j,iblock)/DEGtoRAD+28.9))**2)
          else
           diff_back_nh(i,j,iblock) = diff_back_coef_max*exp(-(0.4d0*(tlat(i,j,iblock)/DEGtoRAD-28.9))**2)
          endif
          diff_back(i,j,iblock)=diff_back_eq+diff_back_sh(i,j,iblock)+diff_back_nh(i,j,iblock)
        if ( tlat(i,j,iblock) .lt. -10.0*degtorad ) then
          diff_back(i,j,iblock) = diff_back(i,j,iblock) + diff_back_coef
        elseif  ( abs(tlat(i,j,iblock)) .le. 10.0*degtorad ) then
          diff_back(i,j,iblock) = diff_back(i,j,iblock) + diff_back_coef*(abs(tlat(i,j,iblock)/DEGtoRAD)/10.0)**2
        else
          diff_back(i,j,iblock) = diff_back(i,j,iblock) + diff_back_coef
        endif
       ENDDO
     ENDDO
   ENDDO
!OMP END PARALLEL DO
#endif
#endif 
#ifdef CANUTO2010

     amld = c0 
     akmt = c0 
     akmu = c0 
     wk1=c0
     wk2=c0
     wk3=c0

#ifdef SW_ATHREAD

#ifdef ENABLE_GPLT
    call t_startf("canuto"//CHAR(0))
#endif
    call athread_spawn(slave_readyc_canuto2010,canuto2010_parm)
    call athread_join()

#ifdef ENABLE_GPLT
    call t_stopf("canuto"//CHAR(0))
#endif
!    if (my_task == 71) then
!      do iblock = 1, nblocks_clinic
!      DO J = 2, JMT-1
!         DO I = 2, IMT-1
!
!         if (VIT(I,J,1,iblock).gt.0.5) then
!            DO k = 1, km
!                !write(*,*) "canuto i,j,k,wk1,wk2,wk3",i,j,k,wk1tmp(i,j,k,iblock),WK2tmp(i,j,k,iblock),WK3tmp(i,j,k,iblock)
!                !write(*,*) "canuto i,j,k,wp3",i,j,k,wp3tmp(i,j,k,iblock)
!                if(akmt(i,j,k,iblock) .ne. 0.0) then
!                    write(*,*) "canuto i,j,k,akt1,akt2,akmt",i,j,k,akt(i,j,k,1,iblock),akt(i,j,k,2,iblock),akmt(i,j,k,iblock)
!                endif
!            end do
!           !write(*,*) "canuto j,i,mldlev",j,i,mldlev(i,j,iblock)
!           ! write(*,*) "canuto i,j,amld,mldlev",i,j,amld(i,j,iblock),mldlev(i,j,iblock)
!         end if
!        end do
!      end do
!    end do
!    end if
#else


!!$OMP PARALLEL DO PRIVATE (iblock,tq,sq,wp1,wp2,wp3,wp4,wp5,wp6,wp7,wp8,wk1,wk2,wk3,wk4)
   do iblock = 1, nblocks_clinic
      DO J = 2, JMT-1
         DO I = 2, IMT-1

        if (VIT(I,J,1,iblock).gt.0.5) then
!
         wp1=0.D0
         wp2=0.D0
         wp3=0.D0
         wp4=0.D0
         wp5=0.D0
         wp6=0.D0
         wp7=0.D0
         wp8=0.D0
         wk1=0.D0
         wk2=0.D0
         wk3=0.D0
!
         do k=1,km - 1
            wp8(k)=-vit(i,j,k+1,iblock)*ZKP(k+1)
         end do
!
         DO K = 1,KMT(i,j,iblock)-1
               wp1(K)= VIT (I,J,K+1,iblock)* (AT (I,J,K,1,iblock)-(AT (I,J,K,1,iblock)-AT (I,J,K+1,1,iblock))* &
                               DZP (K)/(DZP(K)+DZP(K+1)))
               wp2(K)= VIT (I,J,K+1,iblock)* (AT (I,J,K,2,iblock)-(AT (I,J,K,2,iblock)-AT (I,J,K+1,2,iblock))* &
                               DZP (K)/(DZP(K)+DZP(K+1)))*1.0D3+35.
               wp4(k)=vit(i,j,k+1,iblock)*RIT(i,j,k,iblock)
               wp5(k)=vit(i,j,k+1,iblock)*RICDT(i,j,k,iblock)
               wp6(k)=vit(i,j,k+1,iblock)*S2T(i,j,k,iblock)
               wp7(k)=vit(i,j,k+1,iblock)*RICT(i,j,k,iblock)
         END DO
!
         DO K = 1,KMT(i,j,iblock)
               tq = at(i,j,k,1,iblock) - to(1)
               sq = at(i,j,k,2,iblock) - so(1)
               wp3(k) = dens(tq,sq,1)+1.0D3
         end do
!
!YU  May 10th, 2017
         where ( wp7 < 0.0 )  wp7 = 0.0_r8
!YU  May 10th, 2017
         call  canuto_2010_interface(wk1,    wk2,     wk3,     wk4,     amld(i,j,iblock) ,&
                                     wp1,    wp2,     wp3,     wp4,     wp5,              &
                                     wp7,    wp6,     ulat(i,j,iblock)/degtorad ,     wp8,& 
                                     kmt(i,j,iblock) ,i,j, iblock, isc,ustar(i,j,iblock))

!        if (my_task == 71) then
!            DO k = 1, km
!                write(*,*) "canuto i,j,k,wk1,wk2,wk3",i,j,k,wk1(k),wk2(k),wk3(k)
!                !write(*,*) "canuto i,j,k,wp3",i,j,k,wp3(k)
!            end do
!            !write(*,*) "canuto i,j,amld,mldlev",i,j,amld(i,j,iblock),mldlev(i,j,iblock)
!        end if
!        if (mytid == 6 .and. isc ==41 .and. iday == 9 .and. i == 32 .and. j == 17) then
!           write(160,*) i,j,kmt(i,j,1)
!           write(160,*) wk1, wk2, wk3
!           write(160,*) wp1
!           write(160,*) wp2
!           write(160,*) wp3
!           write(160,*) wp4
!           write(160,*) wp5
!           write(160,*) wp6
!           write(160,*) wp7
!           write(160,*) wp8
!           close(160)
!        end if
!        if (isc ==41 .and. iday == 9 .and. i == 32 .and. j == 17) stop
!        CALL  TURB_2(wp8,wp1,wp2,wp3,&
!input RIT and RIDT no unit, S2T in 1/s^2, DFRICMX and DWNDMIX in cm^2/s
!               wp4,wp5,wp6, DFRICMX*1.0d+4, DWNDMIX*1.0d+4,&
!output in cm^2/s, so 1d-4 should be multipled
!              AKM_BACK,AKT_BACK,AKS_BACK,&
!input  RICT in 1/s^2 USTAR in cm/s, BUOYTUR,BUOYSOL in cm^2/s^3,FF in 1/s
!              wp7,wp9,wp10,wp11,FCORT(i,j,iblock) ,& !OK
!output amld in cm, akmt,akh, and aks in cm^2/s
!               AMLD(I,J,iblock),AKMT(I,J,1,iblock),AKT(I,J,1,1,iblock),AKT(I,J,1,2,iblock),&
!              AMLD(I,J,iblock),WK1,WK2,WK3,&
!input int
!              IWK,kmt(I,J,iblock)-1,KM,1,0,0,i,j) !OK!

!output amld in cm, akmt,akh, and aks in cm^2/s
!               AMLD(I,J,iblock),AKMT(I,J,1,iblock),AKT(I,J,1,1,iblock),AKT(I,J,1,2,iblock),&
!              AMLD(I,J,iblock),WK1,WK2,WK3,&
!input int
!              IWK,kmt(I,J,iblock)-1,KM,1,0,0,i,j) !OK!
#if ( defined TIDEMIX )
          
          ak_tide_mixing(:)=0.0
          !LPF2016Nov12 IF (1./OHBT(I,J,IBLOCK).GT.SHELF_CUTOFF) THEN
          !IF (1./OHBT(I,J,IBLOCK).GT.SHELF_CUTOFF) THEN
          DO K = 1,INT(KMT(I,J,IBLOCK))-1
              ak_tide_mixing(K)=back_tidalmixing+mixing_ef*local_mixing_fraction* & 
                                    WAVE_DIS(I,J,IBLOCK)*FZ_TIDE(I,J,K,IBLOCK) &
                                   /(dmax1(rict(I,J,K,iblock),1.0d-8)*WP3(K))     !yuzp-2016/11/23--12/2
              ak_tide_mixing(K)=DMIN1(ak_tide_mixing(K),max_tidalmixing)   !LPF2016Nov12
              !richardson(I,J,K,iblock)=rict(I,J,K,IBLOCK)     !yuzp-2016/11/13--12/2
              !fztidal(I,J,K,iblock)=FZ_TIDE(I,J,K,IBLOCK)     !yuzp-2016/11/13
              !wp3_tidal(I,J,K,iblock)=WP3(K)     !yuzp-2016/11/13
          END DO
          
#if ( defined CANUTOMIXOUT )
          DO K = 1,KMT(i,j,iblock)-1
          wp1_canuto(I,J,K,iblock)=WP1(K)     !yuzp-2016/12/4
          wp2_canuto(I,J,K,iblock)=WP2(K)     !yuzp-2016/12/4
          wp3_canuto(I,J,K,iblock)=WP3(K)     !yuzp-2016/12/4
          wp4_canuto(I,J,K,iblock)=WP4(K)     !yuzp-2016/12/4
          wp5_canuto(I,J,K,iblock)=WP5(K)     !yuzp-2016/12/4
          wp6_canuto(I,J,K,iblock)=WP6(K)     !yuzp-2016/12/4
          wp7_canuto(I,J,K,iblock)=WP7(K)     !yuzp-2016/12/4
          wp8_canuto(I,J,K,iblock)=WP8(K)     !yuzp-2016/12/4
          wp12_canuto(I,J,K,iblock)=wp12(I,J,K,iblock)     !yuzp-2016/12/4
          wp13_canuto(I,J,K,iblock)=wp13(I,J,K,iblock)     !yuzp-2016/12/4
          wk4_canuto(I,J,K,iblock)=WK4(K)     !yuzp-2016/12/4
          END DO
#endif
          DO K = INT(KMT(I,J,IBLOCK))-2,1,-1     !yuzp-2016/11/23
           ak_tide_mixing(K)=DMIN1(ak_tide_mixing(K),ak_tide_mixing(K+1))     !yuzp-2016/11/23
          END DO     !yuzp-2016/11/23
#endif
         DO K = 1,KM
           AKMT(I,J,K,iblock)=WK1(K)
           AKT(I,J,K,1,iblock)=AKT(I,J,K,1,iblock)+ WK2(K)
           AKT(I,J,K,2,iblock)=AKT(I,J,K,2,iblock)+ WK3(K)
#ifdef BCKMEX
           AKMT(I,J,K,iblock) =AKMT(I,J,K,iblock)+diff_back(i,j,iblock)*10.0*1.d-4 !10--Pr_number
           AKT(I,J,K,1,iblock)=AKT(I,J,K,1,iblock)+diff_back(i,j,iblock)*1.d-4
           AKT(I,J,K,2,iblock)=AKT(I,J,K,2,iblock)+diff_back(i,j,iblock)*1.d-4
#endif
#if ( defined TIDEMIX )
           AKMT(I,J,K,iblock)=AKMT(I,J,K,iblock)+  ak_tide_mixing(K)*5.0
           AKT(I,J,K,1,iblock)=AKT(I,J,K,1,iblock)+ak_tide_mixing(K)
           AKT(I,J,K,2,iblock)=AKT(I,J,K,2,iblock)+ak_tide_mixing(K)
#endif 
           akmt(i,j,k,iblock) = min(akmt(i,j,k,iblock), 8.0D-3*DZP(K)*DZP(K))
           akt(i,j,k,1,iblock)= min(akt(i,j,k,1,iblock),8.0D-3*DZP(K)*DZP(K))
           akt(i,j,k,2,iblock)= min(akt(i,j,k,2,iblock),8.0D-3*DZP(K)*DZP(K))
!           if (my_task == 71) then
!                if(akmt(i,j,k,iblock) .ne. 0.0) then
!                    write(*,*) "canuto i,j,k,akt1,akt2,akmt",i,j,k,akt(i,j,k,1,iblock),akt(i,j,k,2,iblock),akmt(i,j,k,iblock)
!                end if
!           end if
         END DO
!       do k =1 ,km
!          ak_tide(i,j,k,iblock) = ak_tide_mixing(k)
!       end do
        endif
         END DO
      END DO
   END DO
#endif 
#endif

#ifdef SW_ATHREAD

#ifdef ENABLE_GPLT
    call t_startf("readyc2"//CHAR(0))
#endif
    call athread_spawn(slave_readyc2,readyc2_parm)
    call athread_join()
#ifdef ENABLE_GPLT
    call t_stopf("readyc2"//CHAR(0))
#endif
#else 
!turb_tu end
!  stop
!!
! calculate the vertical mixing on U-grid

!!!$OMP PARALLEL DO PRIVATE (IBLOCK) 
   do iblock = 1, nblocks_clinic
      DO K = 1,KMM1
         call tgrid_to_ugrid(akmu(:,:,k,iblock), akmt(:,:,k,iblock),iblock)
         do j= 1,jmt
         do i= 1,imt
            akmu(i,j,k,iblock) = akmu(i,j,k,iblock)* viv(i,j,k+1,iblock)
         end do
         end do
      END DO
   end do
!lhl241204
#endif

#if (defined LOWRES)
!LPF20160515
   call POP_HaloUpdate(amld , POP_haloClinic, POP_gridHorzLocCenter,&
                         POP_fieldKindScalar, errorCode, fillValue = 0.0_r8)
!LPF20160515

        !   rict_ref = rict(:,:,14,:) !why 14layer
   DO IBLOCK = 1,NBLOCKS_CLINIC
       DO J = 1,JMT
          DO I = 1, IMT
             rict_ref(i,j,iblock) = rict(i,j,14,iblock) !why 14layer
            DO K = 1,KMM1
            IF (abs(ZKP(K+1)).GT.AMLD(I,J,iblock)) THEN
            rict_ref(I,J,IBLOCK) = rict(I,J,K,IBLOCK)
            exit
            ENDIF
            END DO
         END DO
      END DO
   END DO
#endif
!---------------------------------------------------------------------
!     COMPUTE THE ADVECTIVE TERM: ZONAL COMPONENT
!---------------------------------------------------------------------
#ifdef SW_ATHREAD

#ifdef ENABLE_GPLT
    call t_startf("readyc_upwell"//CHAR(0))
#endif
     work = c0

    call athread_spawn(slave_upwell_cu_0,readyc_upwell_0_param)
    call athread_join()

    call athread_spawn(slave_upwell_cu,readyc_upwell_param)
    call athread_join()

    call athread_spawn(slave_upwell_cu_1,readyc_upwell_1_param)
    call athread_join()

    call athread_spawn(slave_upwell_cu_2,readyc_upwell_2_param)
    call athread_join()

#ifdef ENABLE_GPLT
    call t_stopf("readyc_upwell"//CHAR(0))
#endif
#else
      CALL UPWELL (U,V,H0)
      
#endif      
!---------------------------------------------------------------------
!     INITIALIZE WORK ARRAYS
!---------------------------------------------------------------------
 
     dlu = c0
     dlv = c0
     wka = c0
#ifdef SW_ATHREAD

#ifdef ENABLE_GPLT
    call t_startf("readyc3"//CHAR(0))
#endif
    call athread_spawn(slave_readyc3,readyc3_param)
    call athread_join()
    call athread_spawn(slave_advection_momentum_cu,advection_momentum_param)
    call athread_join()
#ifdef ENABLE_GPLT
    call t_stopf("readyc3"//CHAR(0))
#endif
#else
!!$OMP PARALLEL DO PRIVATE (IBLOCK,K) COLLAPSE(2)
   do iblock = 1, nblocks_clinic
       DO K = 1,KM
          call tgrid_to_ugrid(wka(:,:,k,iblock),ws(:,:,k,iblock), iblock)
       END DO
   end do

!---------------------------------------------------------------------
!     COMPUTE THE ADVECTIVE TERMS
!---------------------------------------------------------------------
 
!!!$OMP PARALLEL DO PRIVATE (IBLOCK)
    DO IBLOCK = 1, NBLOCKS_CLINIC
        call advection_momentum(u(:,:,:,iblock),v(:,:,:,iblock),wka(:,:,:,iblock), &
                                dlu(:,:,:,iblock),dlv(:,:,:,iblock),iblock)
!       dlu(:,:,:,IBLOCK)=adv_uu
!       dlv(:,:,:,IBLOCK)=adv_vv
   END DO
#endif
!
!---------------------------------------------------------------------
!     COMPUTE THE HORIZONTAL VISCOSITY
!---------------------------------------------------------------------
#ifdef SW_ATHREAD
#ifdef ENABLE_GPLT
    call t_startf("hdel4"//CHAR(0))
#endif
   DO IBLOCK = 1, NBLOCKS_CLINIC
      this_block = get_block(blocks_clinic(iblock),iblock)
    am_factor_m = c0

#ifdef BIHAR
    hdiffu_del4_0_param%d_jb           = this_block%JB
    hdiffu_del4_0_param%d_je           = this_block%JE
    hdiffu_del4_0_param%d_ib           = this_block%IB
    hdiffu_del4_0_param%d_ie           = this_block%IE

    hdiffu_del4_1_param%d_jb           = this_block%JB
    hdiffu_del4_1_param%d_je           = this_block%JE
    hdiffu_del4_1_param%d_ib           = this_block%IB
    hdiffu_del4_1_param%d_ie           = this_block%IE

    call athread_spawn(slave_hdel4_div,hdel4_div_param)
    call athread_join()

    call athread_spawn(slave_hdel4_zcurl,hdel4_zcurl_param)
    call athread_join()

    call athread_spawn(slave_hdel4_grad,hdel4_grad_param)
    call athread_join()

    call athread_spawn(slave_hdiffu_del4_0,hdiffu_del4_0_param)
    call athread_join()

    call athread_spawn(slave_hdiffu_del4_1,hdiffu_del4_1_param)
    call athread_join()

!      DO K = 1,KM
!    !      do j = 3, jmt-2
!    !      do i = 3, imt-2
!           do j=this_block%jb-1,this_block%je+1
!            do i=this_block%ib-1,this_block%ie+1
!             if(mytid == 10) then
!               ! if(dlv(i,j,k,iblock) .ne. 0.0) then
!               !     write(*,*) "i,j,k,dlv",i,j,k,dlv (i,j,k,iblock),dlu (i,j,k,iblock)
!                !if(am_factor_m(i,j,k) .ne. 0.0) then
!                 !   write(*,*) "i,j,k,am_factor cpe",i,j,k,am_factor_m (i,j,k)
!                if(d2ukm(i,j,k) .ne. 0.0) then
!                    write(*,*) "i,j,k,D2UK,D2VK",i,j,k,d2ukm (i,j,k),d2vkm (i,j,k)
!                !if(hdukkm(i,j,k) .ne. 0.0) then
!                !    write(*,*) "i,j,k,HDUK,HDVK",i,j,k,hdukkm (i,j,k),hdvkkm (i,j,k),am_f_km(i,j,k)
!                endif
!             endif
!          end do
!          end do
!      END DO

#else
    hdiffu_del2_param%d_jb           = this_block%JB
    hdiffu_del2_param%d_je           = this_block%JE
    hdiffu_del2_param%d_ib           = this_block%IB
    hdiffu_del2_param%d_ie           = this_block%IE
    call athread_spawn(slave_hdiffu_del2,hdiffu_del2_param)
    call athread_join()

#endif

    END DO

#ifdef ENABLE_GPLT
    call t_stopf("hdel4"//CHAR(0))
#endif
#else

#if ( defined SMAG)
!
         CALL SMAG2 (K)
!
#if (defined SMAG_FZ )
#else
#endif
#else
#if (defined BIHAR)

!!!$OMP PARALLEL DO PRIVATE (IBLOCK,this_block,K,HDUK,HDVK)
   DO IBLOCK = 1, NBLOCKS_CLINIC
      this_block = get_block(blocks_clinic(iblock),iblock)
      DO K = 1,KM
          call hdiffu_del4(k, HDUK, HDVK, up(:,:,k,iblock), vp(:,:,k,iblock), this_block)
          do j = 3, jmt-2
          do i = 3, imt-2
             dlv (i,j,k,iblock) = dlv(i,j,k,iblock) + hdvk(i,j)
             dlu (i,j,k,iblock) = dlu(i,j,k,iblock) + hduk(i,j)
          end do
          end do
      END DO
   END DO

#else
!!!$OMP PARALLEL DO PRIVATE (IBLOCK,this_block,hduk,hdvk,K)
   DO IBLOCK = 1, NBLOCKS_CLINIC
      this_block = get_block(blocks_clinic(iblock),iblock)
      DO K = 1,KM
         call hdiffu_del2(k, HDUK, HDVK, up(:,:,k,iblock), vp(:,:,k,iblock), this_block)
         DO J = 3, JMT-2
            DO I = 3,IMT-2
               dlv (i,j,k,iblock) = dlv(i,j,k,iblock) + hdvk(i,j)
               dlu (i,j,k,iblock) = dlu(i,j,k,iblock) + hduk(i,j)
            END DO
         END DO
      END DO
   END DO
 
#endif
#endif

    !!$OMP PARALLEL DO PRIVATE (IBLOCK,K,J,I) COLLAPSE(4)
    DO IBLOCK = 1, NBLOCKS_CLINIC
        DO K = 1,KM
            DO J = 1, JMT
                DO I = 1, IMT
                    WKA (I,J,K,IBLOCK)= C0F * SQRT (UP (I,J,K,IBLOCK)* UP (I,J,K,IBLOCK) + VP (I, &
                            J,K,IBLOCK)* VP (I,J,K,IBLOCK))
                END DO
            END DO
        END DO
    END DO
#endif

#ifdef SW_ATHREAD
#ifdef ENABLE_GPLT
    call t_startf("readyc_vinteg"//CHAR(0))
#endif
    dlub = c0
    dlvb = c0
    call athread_spawn(slave_readyc_vinteg_cu,readyc_vinteg_param)
    call athread_join()
    call athread_spawn(slave_readyc4,readyc4_param)
    call athread_join()
#ifdef ENABLE_GPLT
    call t_stopf("readyc_vinteg"//CHAR(0))

#endif
#else

      if (.not. allocated(dlub)) allocate(dlub(imt,jmt,max_blocks_clinic),dlvb(imt,jmt,max_blocks_clinic),stat=ierr)
!---------------------------------------------------------------------
!     VERTICAL INTEGRATION
!---------------------------------------------------------------------
      CALL VINTEG (DLU,DLUB)
      CALL VINTEG (DLV,DLVB)
!!!$OMP PARALLEL DO PRIVATE (iblock,J,I,kmb)
      do iblock = 1, nblocks_clinic
         DO J = 2, jmt-1
            DO I = 2,imt-1
               kmb = kmu(i,j,iblock)
              if (kmb>=1) then !LPF20170621
               SBCX(I,J,IBLOCK) = SU (I,J,IBLOCK)* OD0
               SBCY(I,J,IBLOCK) = SV (I,J,IBLOCK)* OD0
               BBCX(I,J,IBLOCK)= C0F*SQRT(UP(I,J,kmb,IBLOCK)*UP(I,J,kmb,IBLOCK)+  &
                                          VP(I,J,kmb,IBLOCK)*VP(I,J,kmb,IBLOCK))   &
                          *(UP(I,J,kmb,IBLOCK)*CAG+SNLAT(I,J,IBLOCK)*VP (I,J,kmb,IBLOCK)*SAG)
               BBCY(I,J,IBLOCK)= C0F*SQRT(UP(I,J,kmb,IBLOCK)*UP(I,J,kmb,IBLOCK)+  &  
                                          VP(I,J,kmb,IBLOCK)*VP(I,J,kmb,IBLOCK))&
                          *(-SNLAT(I,J,IBLOCK)*UP(I,J,kmb,IBLOCK)*SAG+VP(I,J,kmb,IBLOCK)*CAG)
              else
                SBCX(I,J,IBLOCK)=0.0D0
                SBCY(I,J,IBLOCK)=0.0D0
                BBCX(I,J,IBLOCK)=0.0D0
                BBCY(I,J,IBLOCK)=0.0D0
              endif
              dlub(i,j,iblock) = dlub(i,j,iblock) +(SBCX(i,j,iblock)-BBCX(i,j,iblock))*ohbu(i,j,iblock)
              dlvb(i,j,iblock) = dlvb(i,j,iblock) +(SBCY(i,j,iblock)-BBCY(i,j,iblock))*ohbu(i,j,iblock)
            ENDDO
         ENDDO
     ENDDO

#endif

 
#ifdef SW_ATHREAD
#ifdef ENABLE_GPLT
    call t_startf("readyc4_2"//CHAR(0))
#endif
    call athread_spawn(slave_readyc4_2,readyc4_2_param)
    call athread_join()
#ifdef ENABLE_GPLT
    call t_stopf("readyc4_2"//CHAR(0))

#endif
#else
!!$OMP PARALLEL DO PRIVATE (IBLOCK,diff_u1,diff_v1,diff_u2,diff_v2) COLLAPSE(4)
    DO IBLOCK = 1, NBLOCKS_CLINIC
      DO K = 1,KM
      DO J = 2, JMT-1
         DO I = 2,IMT-1
!lhl1204
#ifdef CANUTO2010
!       print*,akmU(i,j,k)
            if (k==1) then
               diff_v1 = SV (I,J,IBLOCK)* OD0*(1-AIDIF)
               diff_u1 = SU (I,J,IBLOCK)* OD0*(1-AIDIF)
            else
               diff_v1= AKMU(I,J,K-1,IBLOCK)*(1-AIDIF)*(VP(I,J,K-1,IBLOCK)- VP(I,J,K,IBLOCK))* &
                        ODZT(K)*VIV(I,J,K,IBLOCK)+ &
                        (1.0D0- VIV (I,J,K,IBLOCK))* WKA (I,J,K -1,IBLOCK)*(1-AIDIF) &
                      * (-SNLAT(I,J,IBLOCK)*UP(I,J,K-1,IBLOCK)*SAG+VP(I,J,K-1,IBLOCK)*CAG)
               diff_u1= AKMU(I,J,K-1,IBLOCK)*(1-AIDIF)* (UP(I,J,K-1,IBLOCK)-UP(I,J,K,IBLOCK))*  &
                        ODZT (K)*VIV (I,J,K,IBLOCK) + &
                       (1.0D0- VIV (I,J,K,IBLOCK))* WKA (I,J,K -1,IBLOCK)*(1-AIDIF) &
                      *(UP(I,J,K-1,IBLOCK)*CAG+SNLAT(I,J,IBLOCK)*VP(I,J,K-1,IBLOCK)* SAG)
            end if
            if (k==km) then
               diff_v2= WKA (I,J,KM,IBLOCK)* ( - SNLAT (I,J,IBLOCK)* UP (I,J,KM,IBLOCK)        &
                        * SAG + VP (I,J,KM,IBLOCK)* CAG)*(1-AIDIF)
               diff_u2= WKA (I,J,KM,IBLOCK)* ( UP (I,J,KM,IBLOCK)* CAG + SNLAT (I,J,IBLOCK)    &
                        * VP (I,J,KM,IBLOCK)* SAG)*(1-AIDIF)
            else
               diff_v2= AKMU(I,J,K,IBLOCK)*(1-AIDIF)*(VP(I,J,K,IBLOCK)- VP(I,J,K+1,IBLOCK))* &
                        ODZT(K+1)*VIV(I,J,K+1,IBLOCK)+ &
                        (1.0D0- VIV (I,J,K+1,IBLOCK))* WKA (I,J,K,IBLOCK)*(1-AIDIF) &
                      * (-SNLAT(I,J,IBLOCK)*UP(I,J,K,IBLOCK)*SAG+VP(I,J,K,IBLOCK)*CAG)
               diff_u2= AKMU(I,J,K,IBLOCK)*(1-AIDIF)*(UP(I,J,K,IBLOCK)-UP(I,J,K+1,IBLOCK))*  & 
                        ODZT(K+1)*VIV (I,J,K+1,IBLOCK) + &
                       (1.0D0- VIV (I,J,K+1,IBLOCK))* WKA (I,J,K,IBLOCK)*(1-AIDIF) &
                      *(UP(I,J,K,IBLOCK)*CAG+SNLAT(I,J,IBLOCK)*VP(I,J,K,IBLOCK)* SAG)
            end if
#else
!
            call exit_licom(sigAbort,'The false mixing option')
!
#endif
!lhl1204
            DLV (I,J,K,IBLOCK) = DLV (I,J,K,IBLOCK) + ODZP (K)* (diff_v1-diff_v2)
            DLU (I,J,K,IBLOCK) = DLU (I,J,K,IBLOCK) + ODZP (K)* (diff_u1-diff_u2)
        END DO
      END DO
      END DO
   END DO
#endif
      RETURN
      END SUBROUTINE READYC
 
 
