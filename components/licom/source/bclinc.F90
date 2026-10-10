!  CVS: $Id: bclinc.F90,v 1.1.1.1 2004/04/29 06:22:39 lhl Exp $
!     =================
      SUBROUTINE BCLINC
!     =================
!     INTEGRATION OF MOMENTUM EQUATION
 
#include <def-undef.h>
#ifdef SW_ATHREAD
use sunway_mod
#endif
use precision_mod
use param_mod
use pconst_mod
use dyn_mod
use work_mod
use tracer_mod
use forc_mod, only: psa,su,sv
use domain
use grid
use blocks
use smuvh
use operators

      IMPLICIT NONE
      REAL(r8)    :: AA,GGU,WK1,WK2,fil_lat1,fil_lat2,eee
      integer     :: iblock,kt,irec,kmb
      REAL(r8)    :: gradx(imt,jmt), grady(imt,jmt)
      REAL(r8)    :: LAMDA1(KM)
      integer     :: zwploc(3)
      type (block) :: this_block
      REAL(r8)   :: ts1,te1,clock0f
#if (defined CANUTO) ||(defined CANUTO2010)
      REAL(r8)    :: AIDIF
#endif

!    real(r8), target:: tA8 (imt, jmt, KM),tB8(imt, jmt, KM),tC8(imt,&
!            jmt, KM),tD8(imt, jmt, KM),tE8(imt, jmt, KM+1),tF8(imt, jmt, KM+1)

!lhl0711
#if (defined CANUTO) ||(defined CANUTO2010)
!M
#ifdef SW_ATHREAD
      call athread_spawn(slave_bclinc12_cu_slave_kernel,bclinc12_cu_param)
      call athread_join()
#else 
!!$OMP PARALLEL DO PRIVATE (iblock,J,I)  collapse(3)
      do iblock = 1, nblocks_clinic
         DO J = 2, jmt-1
            DO I = 2,imt-1
              kmb = kmu(i,j,iblock)
              if( kmb>=1 ) then !LPF20170621
               SBCX(I,J,IBLOCK) = SU (I,J,IBLOCK)* OD0
               SBCY(I,J,IBLOCK) = SV (I,J,IBLOCK)* OD0
               BBCX(I,J,IBLOCK)= C0F*SQRT(UP(I,J,kmb,IBLOCK)*UP(I,J,kmb,IBLOCK)+VP(I,J,kmb,IBLOCK)*VP(I,J,kmb,IBLOCK))&
                          *(UP(I,J,kmb,IBLOCK)*CAG+SNLAT(I,J,IBLOCK)*VP (I,J,kmb,IBLOCK)*SAG)
               BBCY(I,J,IBLOCK)= C0F*SQRT(UP(I,J,kmb,IBLOCK)*UP(I,J,kmb,IBLOCK)+VP(I,J,kmb,IBLOCK)*VP(I,J,kmb,IBLOCK))&
                          *(-SNLAT(I,J,IBLOCK)*UP(I,J,kmb,IBLOCK)*SAG+VP(I,J,kmb,IBLOCK)*CAG)
              end if
            ENDDO
         ENDDO
     ENDDO
#endif
      AIDIF=0.5  
!      if (ISC/=0)  AIDIF=0.5
 
#endif
!
!---------------------------------------------------------------------
!      Define the threthold latitute for zonal smoother
       fil_lat1=63.0D0
       fil_lat2=63.0D0
!lhl0711

!---------------------------------------------------------------------
!     ADVECTION + DIFFUSION + CORIOLIS
!---------------------------------------------------------------------
 
!    if (mytid ==68 ) then
!       irec = 16*isc +2
!       write(25,rec=irec) ((buffer_real4(i,j),i=3,imt-2),j=3,jmt-2)
!    end if
!
!---------------------------------------------------------------------
!     PRESSURE GRADIENT FORCES
!---------------------------------------------------------------------
!  if (isc < 10 .and. mytid == 5) then
!    write(100+mytid,*) "OK--2, ISC=", ISC, dlu(14,16,1,1), dlv(14,16,1,1)
!  end if
 
!!@@@@ DP'/DX
#ifdef SW_ATHREAD
      bclinc11_cu_param%d_isc      = isc
      call athread_spawn(slave_bclinc11_cu_slave_kernel,bclinc11_cu_param)
      call athread_join()
      bclinc6_cu_param%d_isc      = isc
      call athread_spawn(slave_bclinc6_cu_slave_kernel,bclinc6_cu_param)
      call athread_join()
#else
      AA = 0.0
      IF (ISC /= 0) AA = 0.5
!!$OMP PARALLEL DO PRIVATE (IBLOCK,J,I) collapse(3)
   DO IBLOCK = 1, NBLOCKS_CLINIC
      DO J = 1,jmt
         DO I = 1,IMT
            H0BF (I,J,IBLOCK)= H0BF (I,J,IBLOCK)* ONBB
         END DO
      END DO
  END DO
 
!!$OMP PARALLEL DO PRIVATE (IBLOCK,J,I) collapse(3)
   DO IBLOCK = 1, NBLOCKS_CLINIC
      DO J = 1, jmt
         DO I = 1,IMT
            WORK (I,J,IBLOCK)= AA * H0BF (I,J,IBLOCK) + (1.0- AA)* H0BL (I,J,IBLOCK)
         END DO
      END DO
   END DO
 
!!$OMP PARALLEL DO PRIVATE (IBLOCK,J,I,WKK) collapse(3)
#if 1
   DO IBLOCK = 1, NBLOCKS_CLINIC
      DO J = 1,jmt
         DO I = 1,IMT
            WKK (1)= (PSA (I,J,IBLOCK)* OD0+ WORK (I,J,IBLOCK)* G)* VIT (I,J,1,IBLOCK)
            !wktmp(i,j,1,IBLOCK) = WKK (1)
            DO K = 1,KM
               WKK (K +1)= WKK (K) - GG (I,J,K,IBLOCK)* DZP (K)* VIT (I,J,K,IBLOCK)
               !wktmp(i,j,k+1,IBLOCK) = WKK (K +1)
            END DO
 
            DO K = 1,KM
               WKA (I,J,K,IBLOCK)= P5* (WKK (K) + WKK (K +1))
            END DO
         END DO
      END DO
  END DO
#else
! Maoxue Yu extract k loop to the outer 20230930
    DO IBLOCK = 1, NBLOCKS_CLINIC
      DO K = 2,KM+1
        DO J = 1,jmt
            DO I = 1,IMT
                if(K .eq. 2) then
                  wktmp(i,j,1,IBLOCK)= (PSA (I,J,IBLOCK)* OD0+ WORK (I,J,IBLOCK)* G)* VIT (I,J,1,IBLOCK)
                endif
                    wktmp(i,j,k,IBLOCK)= wktmp(i,j,k-1,IBLOCK) - GG (I,J,K-1,IBLOCK)* DZP (K-1)* VIT (I,J,K-1,IBLOCK)
                WKA (I,J,K-1,IBLOCK)= P5* (wktmp(i,j,k-1,IBLOCK) + wktmp(i,j,k,IBLOCK))
            END DO
        END DO
      END DO
    END DO

#endif

!!$OMP PARALLEL DO PRIVATE (IBLOCK,K,J,I) collapse(4)
   DO IBLOCK = 1, NBLOCKS_CLINIC
      DO K = 1,KM 
         DO J = 2,jmt-1
            DO I = 2,imt-1
               DLU (I,J,K,IBLOCK) = DLU (I,J,K,IBLOCK) - FCOR(I,J,IBLOCK) * VP (I,J,K,IBLOCK)
               DLV (I,J,K,IBLOCK) = DLV (I,J,K,IBLOCK) + FCOR(I,J,IBLOCK) * UP (I,J,K,IBLOCK)
!              buffer_real4(i,j)= - FCOR(I,J,IBLOCK) * VP (I,J,1,IBLOCK)
            END DO
         END DO
      END DO
  END DO
 
 
!!!$OMP PARALLEL DO PRIVATE (iblock, this_block,gradx,grady,K)
   DO IBLOCK = 1, NBLOCKS_CLINIC
      this_block = get_block(blocks_clinic(iblock),iblock)
      DO K = 1,KM
         call grad(k, GRADX, GRADY, wka(:,:,k,iblock), this_block)
         DO J = 2,jmt-1
            DO I = 2,imt-1
               DLV (I,J,K,IBLOCK) = DLV (I,J,K,IBLOCK) - grady(i,j)
               DLU (I,J,K,IBLOCK) = DLU (I,J,K,IBLOCK) - gradx(i,j)
!              if (k ==1 ) buffer_real4(i,j) =  -gradx(i,j)
            END DO
         END DO
      END DO
  END DO
!
!!@@@@ G'DH/DX
 
!!$OMP PARALLEL DO PRIVATE (IBLOCK,K) collapse(4)
   DO IBLOCK = 1, NBLOCKS_CLINIC
      DO K = 1,KM
         DO J = 1, JMT
            DO I = 1,IMT
               WKA (I,J,K,IBLOCK)= (1.0+ OHBT (I,J,IBLOCK)* ZKT (K))* WORK (I,J,IBLOCK)      &
                            * VIT (I,J,K,IBLOCK)
            END DO
         END DO
      END DO
   END DO
 
 
 
 
!!!$OMP PARALLEL DO PRIVATE (IBLOCK,this_block, K,GGU, gradx,grady)
   DO IBLOCK = 1, NBLOCKS_CLINIC
      this_block = get_block(blocks_clinic(iblock),iblock)
      DO K = 1,KM
         call grad(k, GRADX, GRADY, wka(:,:,k,iblock), this_block)
         DO J = 2, JMT-1
            DO I = 2,IMT-1
               GGU = P25* (GG (I,J,K,IBLOCK) + GG (I -1,J,K,IBLOCK) + GG (I,J +1,K,IBLOCK) &
                     + GG (I -1,J +1,K,IBLOCK))
               DLV (I,J,K,IBLOCK) = DLV (I,J,K,IBLOCK) + GGU * grady(i,j)
               DLU (I,J,K,IBLOCK) = DLU (I,J,K,IBLOCK) + GGU * gradx(i,j)
!              if (k==1) buffer_real4(i,j)= GGU*gradx(i,j)
            END DO
         END DO
      END DO
   END DO

!---------------------------------------------------------------------
!     CORIOLIS ADJUSTMENT
!---------------------------------------------------------------------
 
      IF (ISC == 0) THEN
!!$OMP PARALLEL DO PRIVATE (IBLOCK,K,J,I,WK1,WK2) collapse(4)
   DO IBLOCK = 1, NBLOCKS_CLINIC
         DO K = 1,KM
            DO J = 2, JMT-1
               DO I = 2,IMT-1
                  WK1 = EPEA(I,J,IBLOCK)* DLV(I,J,K,IBLOCK) + EPEB(I,J,IBLOCK)* DLU (I,J,K,IBLOCK)
                  WK2 = EPEA(I,J,IBLOCK)* DLU(I,J,K,IBLOCK) - EPEB(I,J,IBLOCK)* DLV (I,J,K,IBLOCK)
                  DLV (I,J,K,IBLOCK)= WK1* VIV (I,J,K,IBLOCK)
                  DLU (I,J,K,IBLOCK)= WK2* VIV (I,J,K,IBLOCK)

               END DO
            END DO
         END DO
  END DO
      ELSE
!!$OMP PARALLEL DO PRIVATE (IBLOCK,K,J,I,WK1,WK2) collapse(4)
   DO IBLOCK = 1, NBLOCKS_CLINIC
         DO K = 1,KM
            DO J = 2, JMT-1
               DO I = 2,IMT-1
                  WK1 = EPLA(I,J,IBLOCK)* DLV(I,J,K,IBLOCK) + EPLB(I,J,IBLOCK)*DLU(I,J,K,IBLOCK)
                  WK2 = EPLA(I,J,IBLOCK)* DLU(I,J,K,IBLOCK) - EPLB(I,J,IBLOCK)*DLV(I,J,K,IBLOCK)
                  DLV (I,J,K,IBLOCK)= WK1* VIV(I,J,K,IBLOCK)
                  DLU (I,J,K,IBLOCK)= WK2* VIV(I,J,K,IBLOCK)

               END DO
            END DO
         END DO
  END DO
      END IF

#endif
    if (trim(horiz_grid_opt) == 'lat_lon') then
         call POP_HaloUpdate(DLU, POP_haloClinic, POP_gridHorzLocSWcorner , &
                       POP_fieldKindVector, errorCode, fillValue = 0.0_r8)
!
         call POP_HaloUpdate(DLV, POP_haloClinic, POP_gridHorzLocSWcorner , &
                       POP_fieldKindVector, errorCode, fillValue = 0.0_r8)
!
      CALL SMUV_3D (DLV,VIV,fil_lat1)
      CALL SMUV_3D (DLV,VIV,fil_lat1)
    end if
!YU 

!---------------------------------------------------------------------
!     PREDICTING VC & UC
!---------------------------------------------------------------------
 
      IF (ISC < 1) THEN
#ifdef SW_ATHREAD
      call athread_spawn(slave_bclinc71_cu_slave_kernel,bclinc71_cu_param)
      call athread_join()
#else 
!!$OMP PARALLEL DO PRIVATE (IBLOCK,K,J,I) collapse(4)
!ts1 = clock0f()
   DO IBLOCK = 1, NBLOCKS_CLINIC
      DO K = 1,KM
         DO J = 3,JMT-2
            DO I = 3,IMT-2
               V (I,J,K,IBLOCK)= VP (I,J,K,IBLOCK) + DLV (I,J,K,IBLOCK)* DTC
               U (I,J,K,IBLOCK)= UP (I,J,K,IBLOCK) + DLU (I,J,K,IBLOCK)* DTC
            END DO
         END DO
      END DO
  END DO
#endif

!
!    if (mytid ==68 ) then
!       buffer_real4(:,:)= u(:,:,1,1)
!       irec = 16*isc +12
!       write(25,rec=irec) ((buffer_real4(i,j),i=3,imt-2),j=3,jmt-2)
!    end if
!    write(mytid+120,*) ((u(i,j,30,2),i=3,imt-2),j=3,jmt-2)
!    write(mytid+120,*) 
!    write(mytid+120,*) ((sbcx(i,j,2),i=3,imt-2),j=3,jmt-2)
!    write(mytid+120,*) 
!    write(mytid+120,*) ((bbcx(i,j,2),i=3,imt-2),j=3,jmt-2)
!    write(mytid+120,*) 
!    write(mytid+120,*) ((akmu(i,j,29,2),i=3,imt-2),j=3,jmt-2)
!    write(mytid+120,*) 
!    write(mytid+120,*) ((akmu(i,j,30,2),i=3,imt-2),j=3,jmt-2)
!    write(mytid+120,*) 
!    write(mytid+120,*) ((kmu(i,j,2),i=3,imt-2),j=3,jmt-2)
#ifdef SW_ATHREAD
      call athread_spawn(slave_bclinc_invtriu_kernel,bclinc_invtriu_param)
      call athread_join()
#else
     CALL INVTRIU(U,SBCX,BBCX,AKMU,AIDIF,DTC)
!    write(mytid+120,*) 
!    write(mytid+120,*) ((u(i,j,30,2),i=3,imt-2),j=3,jmt-2)
!    close(120+mytid)
     CALL INVTRIU(V,SBCY,BBCY,AKMU,AIDIF,DTC)
#endif
!    stop
!
!    if (mytid ==68 ) then
!       buffer_real4(:,:)= u(:,:,1,1)
!       irec = 16*isc +13
!       write(25,rec=irec) ((buffer_real4(i,j),i=3,imt-2),j=3,jmt-2)
!    end if
         call POP_HaloUpdate3DR8_opt(U, POP_haloClinic, POP_gridHorzLocSWcorner , &
                       POP_fieldKindVector, km, errorCode, fillValue = 0.0_r8)
!
         call POP_HaloUpdate3DR8_opt(V, POP_haloClinic, POP_gridHorzLocSWcorner , &
                       POP_fieldKindVector, km, errorCode, fillValue = 0.0_r8)
!
!  if (isc < 10 .and. mytid == 5) then
!    write(100+mytid,*) "OK--6, ISC=", ISC, u(14,16,1,1), v(14,16,1,1)
!  end if
!
!---------------------------------------------------------------------
!@@@  INTERACTION BETWEEN BAROTROPIC AND BAROCLINIC MODES
!---------------------------------------------------------------------
!ts1 = clock0f()
#ifdef SW_ATHREAD

      wka1 = u
      wka2 = v
      work = c0
      vinteg_param%vinegtype       =  1
      call athread_spawn(slave_vinteg_cu_slave_kernel,vinteg_param)
      call athread_join()

      !CALL VINTEG (U,WORK)
      call athread_spawn(slave_bclinc81_cu,bclinc81_cu_param)
      call athread_join()
#else 
!ts1 = clock0f()
      CALL VINTEG (U,WORK)
!!$OMP PARALLEL DO PRIVATE (IBLOCK,K,J,I) collapse(4)
   DO IBLOCK = 1, NBLOCKS_CLINIC
      DO K = 1,KM
         DO J = 1, JMT
            DO I = 1,IMT
               U (I,J,K,IBLOCK)= (U (I,J,K,IBLOCK) - WORK (I,J,IBLOCK) + UB (I,J,IBLOCK))* VIV (I,J,K,IBLOCK)
               !UTF (I,J,K,IBLOCK)= UTF (I,J,K,IBLOCK) + U (I,J,K,IBLOCK)
            END DO
         END DO
      END DO
  END DO
#endif
!if(mytid == 0) write(*,*) "bclinc81 spend",te1 - ts1
!if(my_task == master_task) then
!   DO IBLOCK = 1, NBLOCKS_CLINIC
!         DO J = 1, JMT
!            DO I = 1,IMT
!                if(work(i,j,iblock) .ne. 0.0e0) write(20,*) 'work i,j',i,j,work(i,j,iblock)
!               !UTF (I,J,K,IBLOCK)= UTF (I,J,K,IBLOCK) + U (I,J,K,IBLOCK)
!            END DO
!         END DO
!  END DO
!  endif
!    if (mytid ==68 ) then
!       buffer_real4(:,:)= work(:,:,1)
!       irec = 16*isc +14
!       write(25,rec=irec) ((buffer_real4(i,j),i=3,imt-2),j=3,jmt-2)
!       buffer_real4(:,:)= ub(:,:,1)
!       irec = 16*isc +15
!       write(25,rec=irec) ((buffer_real4(i,j),i=3,imt-2),j=3,jmt-2)
!       buffer_real4(:,:)= gg(:,:,1,1)
!       irec = 16*isc +16
!       write(25,rec=irec) ((buffer_real4(i,j),i=3,imt-2),j=3,jmt-2)
!    end if
!

!!$OMP PARALLEL DO PRIVATE (IBLOCK,K,J,I) collapse(4)
#ifdef SW_ATHREAD
      work = c0
      vinteg_param%vinegtype       =  2
      call athread_spawn(slave_vinteg_cu_slave_kernel,vinteg_param)
      call athread_join()
!      CALL VINTEG (V,WORK)
      call athread_spawn(slave_bclinc82_cu,bclinc82_cu_param)
      call athread_join()
#else
   CALL VINTEG (V,WORK)
   DO IBLOCK = 1, NBLOCKS_CLINIC
      DO K = 1,KM
         DO J = 1, JMT
            DO I = 1,IMT
               V (I,J,K,IBLOCK)= (V (I,J,K,IBLOCK) - WORK (I,J,IBLOCK) + VB (I,J,IBLOCK))* VIV (I,J,K,IBLOCK)
               !VTF (I,J,K,IBLOCK)= VTF (I,J,K,IBLOCK) + V (I,J,K,IBLOCK)
            END DO
         END DO
      END DO
  END DO

#endif

       !CALL VINTEG (V,WORK)
!  if (isc < 10 .and. mytid == 5) then
!    write(100+mytid,*) "OK--6, ISC=", ISC, u(14,16,1,1), v(14,16,1,1)
!  end if
 
      ISC = ISC +1
    
      ELSE 

#ifdef SW_ATHREAD
      call athread_spawn(slave_bclinc2_cu_slave_kernel,bclinc2_cu_param)
      call athread_join()
#else 
!ts1= clock0f()

!!$OMP PARALLEL DO PRIVATE (IBLOCK,K,J,I) collapse(4)
   DO IBLOCK = 1, NBLOCKS_CLINIC
      DO K = 1,KM
         DO J = 3,JMT-2
            DO I = 3,IMT-2
               WKA (I,J,K,IBLOCK)= VP (I,J,K,IBLOCK) + DLV (I,J,K,IBLOCK)* DTC2
            END DO
         END DO
      END DO
   END DO

#endif

#ifdef SW_ATHREAD

      bclinc3_cu_param%d_wk       => wka(1,1,1,1)
      bclinc3_cu_param%d_bombc    => bbcy(1, 1, 1)
      bclinc3_cu_param%d_topbc    => sbcy(1, 1, 1)

      call athread_spawn(slave_bclinc3_cu_slave_kernel,bclinc3_cu_param)
      call athread_join()

#else
     CALL INVTRIU(WKA,SBCY,BBCY,AKMU,AIDIF,DTC2)
#endif
         call POP_HaloUpdate3DR8_opt(WKA, POP_haloClinic, POP_gridHorzLocSWcorner , &
                       POP_fieldKindVector, km, errorCode, fillValue = 0.0_r8)

#ifdef SW_ATHREAD
      work = c0
      vinteg_param%vinegtype       =  3
      call athread_spawn(slave_vinteg_cu_slave_kernel,vinteg_param)
      call athread_join()
!      CALL VINTEG (WKA,WORK)
      call athread_spawn(slave_bclinc4_cu_0_slave_kernel,bclinc4_0_cu_param)
      call athread_join()
      call athread_spawn(slave_bclinc4_cu_slave_kernel,bclinc4_0_cu_param)
      call athread_join()
!      call athread_spawn(slave_bclinc4_1_cu_slave_kernel,bclinc4_cu_param)
!      call athread_join()
      bclinc3_cu_param%d_topbc    => SBCX(1, 1, 1)
      bclinc3_cu_param%d_bombc    => BBCX(1, 1, 1)
      call athread_spawn(slave_bclinc3_cu_slave_kernel,bclinc3_cu_param)
      call athread_join()
     ! CALL INVTRIU(WKA,SBCX,BBCX,AKMU,AIDIF,DTC2)
#else

      CALL VINTEG (WKA,WORK)
!!$OMP PARALLEL DO PRIVATE (IBLOCK,K,J,I) collapse(4)
   DO IBLOCK = 1, NBLOCKS_CLINIC
      DO K = 1,KM
         DO J = 1, JMT
            DO I = 1,IMT
               WKA (I,J,K,IBLOCK)= (WKA (I,J,K,IBLOCK) - WORK (I,J,IBLOCK) + VB (I,J,IBLOCK))* VIV (I,J,K,IBLOCK)
            END DO
         END DO
      END DO
  END DO

!---------------------------------------------------------------------
!     FILTER FORCING AT HIGT LATITUDES
!---------------------------------------------------------------------
 
 
!!$OMP PARALLEL DO PRIVATE (IBLOCK,K,J,I) collapse(4)
   DO IBLOCK = 1, NBLOCKS_CLINIC
      DO K = 1,KM
         DO J = 1,JMT
            DO I = 1,IMT
               VP (I,J,K,IBLOCK) = AFC2* V (I,J,K,IBLOCK) + AFC1* (VP (I,J,K,IBLOCK) + WKA (I,J,K,IBLOCK))
               V (I,J,K,IBLOCK) = WKA (I,J,K,IBLOCK)
            END DO
         END DO
      END DO
  END DO


!!$OMP PARALLEL DO PRIVATE (IBLOCK,K,J,I) collapse(4)
   DO IBLOCK = 1, NBLOCKS_CLINIC
      DO K = 1,KM
         DO J = 3,jmt-2
            DO I = 3,imt-2
               WKA (I,J,K,IBLOCK)= UP (I,J,K,IBLOCK) + DLU (I,J,K,IBLOCK)* DTC2
            END DO
         END DO
      END DO
  END DO
!
     CALL INVTRIU(WKA,SBCX,BBCX,AKMU,AIDIF,DTC2)
!
#endif
         call POP_HaloUpdate3DR8_opt(WKA, POP_haloClinic, POP_gridHorzLocSWcorner , &
                       POP_fieldKindVector, km, errorCode, fillValue = 0.0_r8)
#ifdef SW_ATHREAD
     work = c0
     vinteg_param%vinegtype       =  3
     call athread_spawn(slave_vinteg_cu_slave_kernel,vinteg_param)
     call athread_join()
!    CALL VINTEG (WKA,WORK)
    call athread_spawn(slave_bclinc5_cu_slave_kernel,bclinc5_cu_param)
    call athread_join()
    call athread_spawn(slave_bclinc5_1_cu_slave_kernel,bclinc5_cu_param)
    call athread_join()
#else
      CALL VINTEG (WKA,WORK)
!
!!$OMP PARALLEL DO PRIVATE (IBLOCK,K,J,I) collapse(4)
   DO IBLOCK = 1, NBLOCKS_CLINIC
      DO K = 1,KM
         DO J = 1, JMT
            DO I = 1,IMT
               WKA (I,J,K,IBLOCK)= (WKA (I,J,K,IBLOCK) - WORK (I,J,IBLOCK) + UB (I,J,IBLOCK))* VIV (I,J,K,IBLOCK)
            END DO
         END DO
      END DO
   END DO
 
!---------------------------------------------------------------------
!     FILTER FORCING AT HIGT LATITUDES
!---------------------------------------------------------------------
 
!!$OMP PARALLEL DO PRIVATE (IBLOCK,K,J,I) collapse(4)
   DO IBLOCK = 1, NBLOCKS_CLINIC
      DO K = 1,KM
         DO J = 1, JMT
            DO I = 1,IMT
               UP (I,J,K,IBLOCK) = AFC2* U (I,J,K,IBLOCK) + AFC1* (UP (I,J,K,IBLOCK) + WKA (I,J,K,IBLOCK))
               U (I,J,K,IBLOCK) = WKA (I,J,K,IBLOCK)
            END DO
         END DO
      END DO
   END DO
#endif

!YU  Oct. 24, 1605
!lhl0711     IF (MOD(ISC,60)==0) THEN
    if (trim(horiz_grid_opt) == 'lat_lon') then
     IF (MOD(ISC,160)==1) THEN
        CALL SMUV_3D (U,VIV,fil_lat2)
        CALL SMUV_3D (V,VIV,fil_lat2)
        CALL SMUV_3D (UP,VIV,fil_lat2)
        CALL SMUV_3D (VP,VIV,fil_lat2)
     END IF
    end if
!
      ISC = ISC +1
      END IF
!     deallocate(buffer_real4)
!!$OMP PARALLEL DO PRIVATE (IBLOCK,K,J,I) collapse(4)
#ifdef SW_ATHREAD

    call athread_spawn(slave_bclinc5_2_cu_slave_kernel,bclinc5_cu_param)
    call athread_join()
#else
     if ( enoi_iau ) then
     do k=1,45 !kvt
      LAMDA1(k)=1.0D0/(iau_gamma*86400.D0)
     enddo
     do k=46,km
      LAMDA1(k)=1.0D0/(iau_gamma*86400.D0) !/(k*k/5.0/5.0)
     enddo
!          if(mytid==0) write(*,*) 'into IAU,uu&vv'
    endif

   DO IBLOCK = 1, NBLOCKS_CLINIC
      DO K = 1,KM
         DO J = JST,JET
            DO I = 1,IMT
                if ( enoi_iau ) then
                    UTF (I,J,K,IBLOCK)= UTF (I,J,K,IBLOCK) + U (I,J,K,IBLOCK) &
                            + VIV (I,J,K,iblock)*uu_iau (I,J,K,iblock)*LAMDA1(k) !*DTC
                    VTF (I,J,K,IBLOCK)= VTF (I,J,K,IBLOCK) + V (I,J,K,IBLOCK) &
                            - VIV (I,J,K,iblock)*vv_iau (I,J,K,iblock)*LAMDA1(k) !*DTC
                else
                    UTF (I,J,K,IBLOCK)= UTF (I,J,K,IBLOCK) + U (I,J,K,IBLOCK)
                    VTF (I,J,K,IBLOCK)= VTF (I,J,K,IBLOCK) + V (I,J,K,IBLOCK)
               end if
            END DO
         END DO
      END DO
  END DO
!ZWP20170212
   DO IBLOCK = 1, NBLOCKS_CLINIC
      this_block = get_block(blocks_clinic(iblock),iblock)
!      write(253,*)'isc=',isc,'mytid=',mytid,'i_global=',this_block%i_glob,'j_global=',this_block%j_glob
      DO K = 1,KM
         DO J = JST,JET
            DO I = 1,IMT
      ucheck(i,j,k) = u(i,j,k,iblock) 
      vcheck(i,j,k) = v(i,j,k,iblock) 
            END DO
         END DO
      END DO
   END DO
!ZWP20170212
#endif
    DO IBLOCK = 1, NBLOCKS_CLINIC
           this_block = get_block(blocks_clinic(iblock),iblock)
  if (maxval(abs(ucheck)) > 5.0D0 ) then
         zwploc = maxloc(abs(ucheck))
!        write(6,*) "UMAX", isc, mytid,maxval(abs(u)),maxloc(abs(u))

!huangjian change
!        write(6,'(1x,A5,1x,2I8,1x,F24.20,1x,3I3,1x,2I6)') "UMAX",isc, mytid,&
!             maxval(abs(ucheck)),maxloc(abs(ucheck)),&
!             this_block%i_glob(zwploc(1)),this_block%j_glob(zwploc(2))
!ZWP20170209
  end if
  if (maxval(abs(vcheck)) > 5.0D0 ) then
         zwploc = maxloc(abs(vcheck))
!        write(6,*) "VMAX", isc, mytid,maxval(abs(v)),maxloc(abs(v))
        write(6,'(1x,A5,1x,2I8,1x,F24.20,1x,3I3,1x,2I6)') "VMAX",isc, mytid,&
             maxval(abs(vcheck)),maxloc(abs(vcheck)),&
             this_block%i_glob(zwploc(1)),this_block%j_glob(zwploc(2))
  end if

  END DO
!ZWP2017dd0209

!  call mpi_barrier(mpi_comm_ocn,ierr)

      RETURN
      END SUBROUTINE BCLINC
