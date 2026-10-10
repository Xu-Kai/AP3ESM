!  CVS: $Id: isoadv.F90,v 1.1.1.1 2004/04/29 06:22:39 lhl Exp $
#include <def-undef.h>
 
#if (defined ISO)
!     =================
      SUBROUTINE ISOADV
!     =================
 
!     compute isopycnal transport velocities.
use precision_mod 
use param_mod
use pconst_mod
use isopyc_mod
use msg_mod
use domain
use grid
use pmix_mod
use tracer_mod, only : amld
#ifdef ISOOUT
use dyn_mod,only:vetiso,vntiso,vbtiso
#endif


      IMPLICIT NONE
 
      REAL(r8):: p5,c0,fxa_param
      INTEGER :: jstrt, iblock
       REAL(r8):: fxa(imt,jmt,km,max_blocks_clinic)
      REAL(r8):: athkdfx(imt,jmt,km,max_blocks_clinic),athkdfy(imt,jmt,km,max_blocks_clinic)  
      !REAL(r8):: fxa(imt,jmt,km,nblocks_clinic)
      !REAL(r8):: athkdfx(imt,jmt,km,nblocks_clinic), athkdfy(imt,jmt,km,nblocks_clinic) 
#ifdef LPFDIAG
      real(r8)::outsum,temp_tmp(imt,jmt,km,max_blocks_clinic) 
#endif
!================================================
!        compute athkdf              lyw_20160128
!===============================================

!    athkdf_param = 0.01D0
    athkdf_param = 0.0D0
    athkdfx=0.0D0
    athkdfy=0.0D0

    Do iblock = 1, nblocks_clinic
      DO k = 1,km-1
          DO j = 1,jmt
             DO i = 1,imt
             if (VIT (I,J,K +1,IBLOCK).gt.0.5) then
!                if (k.le.14) then
                IF (abs(ZKP(K+1)).LT.AMLD(I,J,iblock)) THEN
                athkdf_param(i,j,k,iblock) =1.0D0
                else
                athkdf_param(i,j,k,iblock) = rict_replace(i,j,k,iblock)/dmax1(rict_ref(i,j,iblock),0.0001)
                athkdf_param(i,j,k,iblock) = dmin1(athkdf_param(i,j,k,iblock),1.0D0)
                athkdf_param(i,j,k,iblock) = dmax1(athkdf_param(i,j,k,iblock),0.25D0)
                endif
             else
                athkdf_param(i,j,k,iblock) = 0.0D0
             endif
!                if (abs(rict_ref(i,j,iblock)) .lt. 0.0001) then
!                athkdf_param(i,j,k,iblock) = 1.0D0 
!                else
!                athkdf_param(i,j,k,iblock) = rict_replace(i,j,k,iblock)/rict_ref(i,j,iblock)
!                end if
             END DO
          END DO
       END DO
     END DO   

    DO iblock = 1, nblocks_clinic
          DO j = 1,jmt
             DO i = 1,imt
     IF (KMT(i,j,iblock).LT.6) THEN
      DO k = 1,km-1
      athkdf_param(i,j,k,iblock) = 1.0D0
      END DO
     ENDIF
             END DO
          END DO
     END DO   

!    DO iblock = 1, nblocks_clinic
!          DO j = 1,jmt
!             DO i = 1,imt
!     IF (1./OHBT(I,J,IBLOCK).LT.60) THEN
!      DO k = 1,km-1
!      athkdf_param(i,j,k,iblock) = 1.0
!      END DO
!     ENDIF
!             END DO
!          END DO
!     END DO   

!    Do iblock = 1, nblocks_clinic
!      DO k = 2,km-1
!          DO j = 1,jmt
!             DO i = 1,imt
!             if (VIT (I,J,K +1,IBLOCK).gt.0.5) then
!                athkdf_param(i,j,k,iblock) = dmin1(athkdf_param(i,j,k-1,iblock),athkdf_param(i,j,k,iblock))
!             endif
!             END DO
!          END DO
!       END DO
!     END DO   

#if (defined ISO_TYPE_BF)
    DO iblock = 1, nblocks_clinic
      DO k = 1,km-1
          DO j = 1,jmt
             DO i = 1,imt
      athkdf(i,j,k,iblock)=athkdf_param(i,j,k,iblock)*tmask(i,k+1,j,iblock)*2.0D3
             END DO
          END DO
      END DO
     END DO   
#else
      athkdf = 500.0D0
#endif
    
      !deallocate(rict_ref,rict_replace)

    DO iblock = 1, nblocks_clinic
      DO k = 2,km
          DO j = 1,jmt
             DO i = 1,imt-1
      athkdfx(i,j,k,iblock) = 0.25*tmask (i,k,j,iblock)* tmask (i+1,k,j,iblock)&
                            *(athkdf(i,j,k-1,iblock)+athkdf(i+1,j,k-1,iblock)+athkdf(i,j,k,iblock)+athkdf(i+1,j,k,iblock))
             END DO
          END DO
      END DO
     END DO   

    DO iblock = 1, nblocks_clinic
      DO k = 2,km
          DO j = 1,jmt-1
             DO i = 1,imt
      athkdfy(i,j,k,iblock) = 0.25*tmask (i,k,j,iblock)* tmask (i,k,j+1,iblock)&
                            *(athkdf(i,j,k-1,iblock)+athkdf(i,j+1,k-1,iblock)+athkdf(i,j,k,iblock)+athkdf(i,j+1,k,iblock))
             END DO
          END DO
      END DO
     END DO   

    DO iblock = 1, nblocks_clinic
          DO j = 1,jmt
             DO i = 1,imt-1
      athkdfx(i,j,1,iblock) = 0.5*tmask (i,1,j,iblock)* tmask (i+1,1,j,iblock)&
                            * (athkdf(i,j,1,iblock)+athkdf(i+1,j,1,iblock))
             END DO
          END DO
     END DO   

    DO iblock = 1, nblocks_clinic
          DO j = 1,jmt-1 !LPF20160515
             DO i = 1,imt
          !DO j = 1,jmt
             !DO i = 1,imt-1
      athkdfy(i,j,1,iblock) = 0.5*tmask (i,1,j,iblock)* tmask (i,1,j+1,iblock)&
                            * (athkdf(i,j,1,iblock)+athkdf(i,j+1,1,iblock))
             END DO
          END DO
     END DO   

    DO iblock = 1, nblocks_clinic
          DO j = 1,jmt
             DO i = 1,imt-1
            k = min (KMT(i,j,iblock),KMT(i+1,j,iblock))
            IF (k /= 0) THEN
      athkdfx(i,j,k,iblock) = 0.5*tmask (i,k,j,iblock)* tmask (i+1,k,j,iblock)&
                            *(athkdf(i,j,k-1,iblock)+athkdf(i+1,j,k-1,iblock))
            END IF
             END DO
          END DO
     END DO   

    DO iblock = 1, nblocks_clinic
          DO j = 1,jmt-1
             DO i = 1,imt
            k = min (KMT(i,j,iblock),KMT(i,j+1,iblock))
            IF (k /= 0) THEN
      athkdfy(i,j,k,iblock) = 0.5*tmask (i,k,j,iblock)* tmask (i,k,j+1,iblock)&
                            *(athkdf(i,j,k-1,iblock)+athkdf(i,j+1,k-1,iblock))
            END IF
             END DO
          END DO
     END DO   

!      if (athkdf_param(:,:,:,:) > 1 ) then
!          athkdf_param(:,:,:,:) = 1.0D0
!       else
!          athkdf_param(:,:,:,:) = athkdf_param(:,:,:,:)
!       end if
!       if (athkdf_param(:,:,:,:) < 0.01 ) then
!          athkdf_param(:,:,:,:) = 0.01D0
!       else
!          athkdf_param(:,:,:,:) = athkdf_param(:,:,:,:)
!       end if
      

!-----------------------------------------------------------------------
!     set local constants
!-----------------------------------------------------------------------
      p5 = 0.5D0
      c0 = 0.0D0
 
!-----------------------------------------------------------------------
!     compute the meridional component of the isopycnal mixing velocity
!     at the center of the northern face of the "t" cells.
!-----------------------------------------------------------------------
      adv_vntiso=c0 
      adv_vbtiso=c0 
      adv_vetiso=c0 
!!$OMP PARALLEL DO PRIVATE (iblock,fxa)
   do iblock = 1, nblocks_clinic
      DO j = 2,jmt-1
         DO k = 2,km -1
            DO i = 2,imt-1
!               fxa = - p5* dzr(k)*athkdf(i,j,iblock)
!            fxa_param = -p5*dzr(k)    !lyw 20160204  sinu(j) is removed ?
            fxa_param = -1.0/(1.0/dzr(k)+1.0/dzr(k-1)*0.5+1.0/dzr(k+1)*0.5)    !lyw 20160204  sinu(j) is removed ?
!            fxa(i,j,k,iblock) = fxa_param*athkdf(i,j,k,iblock)
            fxa(i,j,k,iblock) = fxa_param*athkdfy(i,j,k,iblock)
!            adv_vntiso (i,k,j,iblock) = fxa(i,j,k,iblock) * tmask (i,k,j,iblock)* tmask (i,k,j +1,iblock)* ( &
!                                     K2 (i,k -1,j,3,iblock) - K2 (i,k +1,j,3,iblock)) &
!                                     +K2 (i,k,j,3,iblock)*fxa_param*(athkdf(i,j,k-1,iblock)-athkdf(i,j,k+1,iblock))
            adv_vntiso (i,k,j,iblock) = fxa(i,j,k,iblock) * tmask (i,k,j,iblock)* tmask (i,k,j+1,iblock)* &
                                     (K2 (i,k -1,j,3,iblock) - K2 (i,k +1,j,3,iblock)) &
                                     +tmask (i,k,j,iblock)* tmask (i,k,j+1,iblock)&
                                     *K2 (i,k,j,3,iblock)*fxa_param*(athkdfy(i,j,k-1,iblock)-athkdfy(i,j,k+1,iblock))
!               adv_vntiso (i,k,j,iblock) = fxa * tmask (i,k,j,iblock)* tmask (i,k,j +1,iblock)* ( &
!                                    K2 (i,k -1,j,3,iblock) - K2 (i,k +1,j,3,iblock))
            END DO
         END DO
      END DO
   END DO
 
 
!     consider the top and bottom levels. "K2" is assumed to be zero
!     at the ocean top and bottom.
 
!     k = 1
!!$OMP PARALLEL DO PRIVATE (iblock,fxa)
   do iblock = 1, nblocks_clinic
      DO j = 2,jmt-1
         DO i = 2,imt-1
!            fxa = - p5* dzr (1)* athkdf(i,j,iblock)
            fxa_param = - p5* dzr (1)
!            fxa_param = - dzr (1)
!            adv_vntiso (i,1,j,iblock) = - fxa_param *athkdf(i,j,1,iblock)* tmask (i,1,j,iblock)* tmask (i,1,j +1,iblock)&
!                                 * (K2 (i,1,j,3,iblock) + K2 (i,2,j,3,iblock))   
            adv_vntiso (i,1,j,iblock) = - fxa_param *athkdfy(i,j,1,iblock)* tmask (i,1,j,iblock)* tmask (i,1,j +1,iblock)&
                                 * (K2 (i,1,j,3,iblock) + K2 (i,2,j,3,iblock))
         END DO
      END DO
   END DO
 
 
 
 
!!$OMP PARALLEL DO PRIVATE (iblock)
   do iblock = 1, nblocks_clinic
      DO j = 2,jmt-1
         DO i = 2, imt-1
            k = min (KMT(i,j,iblock),KMT(i,j +1,iblock))
            IF (k /= 0) THEN
!               adv_vntiso (i,k,j,iblock) = - p5*dzr(k)*athkdf(i,j,k,iblock)*tmask(i,k,j,iblock) &
!                                    * tmask (i,k,j +1,iblock)* (K2 (i,k,j,3,iblock)   &
!                                      + K2 (i,k -1,j,3,iblock))          
!               adv_vntiso (i,k,j,iblock) = -dzr(k)*athkdfy(i,j,k,iblock)*&
               adv_vntiso (i,k,j,iblock) = -p5*dzr(k)*athkdfy(i,j,k,iblock)*&
                                      tmask(i,k,j,iblock)*tmask (i,k,j +1,iblock)* &
                                     (K2 (i,k,j,3,iblock) + K2 (i,k -1,j,3,iblock))
            END IF
         END DO
      END DO
   END DO

 
!-----------------------------------------------------------------------
!     compute the zonal component of the isopycnal mixing velocity
!     at the center of the eastern face of the "t" grid box.
!-----------------------------------------------------------------------
 
      jstrt = 2
 
!!$OMP PARALLEL DO PRIVATE (iblock,fxa)
   do iblock = 1, nblocks_clinic
      DO j = 2, jmt-1
         DO k = 2,km -1
!            fxa_param = -p5*dzr (k)
            fxa_param = -1./(1./dzr(k)+1./dzr(k-1)*0.5+1./dzr(k+1)*0.5)    !lyw 20160204  sinu(j) is removed ?
            DO i = 2,imt-1
!               fxa = - p5* dzr (k)* athkdf(i,j,iblock)
!               fxa(i,j,k,iblock) = fxa_param*athkdf(i,j,k,iblock)
               fxa(i,j,k,iblock) = fxa_param*athkdfx(i,j,k,iblock)
!               adv_vetiso (i,k,j,iblock) = fxa(i,j,k,iblock) * tmask (i,k,j,iblock)* tmask (i +1,k,j,iblock) &
!                                    * (K1 (i,k -1,j,3,iblock) - K1 (i,k +1,j,3,iblock)) &
!                                    +K1(i,k,j,3,iblock)*fxa_param*(athkdf(i,j,k-1,iblock)-athkdf(i,j,k+1,iblock))
               adv_vetiso (i,k,j,iblock) = fxa(i,j,k,iblock) * tmask (i,k,j,iblock)* tmask (i +1,k,j,iblock) &
                                    * (K1 (i,k -1,j,3,iblock) - K1 (i,k +1,j,3,iblock)) &
                                    +tmask (i,k,j,iblock)* tmask (i +1,k,j,iblock)&
                                    * K1(i,k,j,3,iblock)*fxa_param*(athkdfx(i,j,k-1,iblock)-athkdfx(i,j,k+1,iblock))
!               adv_vetiso (i,k,j,iblock) = fxa * tmask (i,k,j,iblock)* tmask (i +1,k,j,iblock) &
!                                    * (K1 (i,k -1,j,3,iblock) - K1 (i,k +1,j,3,iblock))      
            END DO
         END DO
      END DO
   END DO
 
 
!     consider the top and bottom levels. "K1" is assumed to be zero
!     at the ocean top and bottom.
 
!     k = 1
!!$OMP PARALLEL DO PRIVATE (iblock,fxa)
   fxa_param = -p5*dzr(1)
!   fxa_param = -dzr(1)
   do iblock = 1, nblocks_clinic
      DO j = 2, jmt-1
         DO i = 2,imt-1
!            fxa = - p5* dzr (1)* athkdf(i,j,iblock)
!            adv_vetiso (i,1,j,iblock) = - fxa_param *athkdf(i,j,1,iblock) * tmask (i,1,j,iblock)* tmask (i +1,1,j,iblock) &
!                                 * (K1 (i,1,j,3,iblock) + K1 (i,2,j,3,iblock)) 
            adv_vetiso (i,1,j,iblock) = - fxa_param *athkdfx(i,j,1,iblock) * tmask (i,1,j,iblock)* tmask (i +1,1,j,iblock) &
                                 * (K1 (i,1,j,3,iblock) + K1 (i,2,j,3,iblock))
         END DO
      END DO
   END DO
 
 
!!$OMP PARALLEL DO PRIVATE (iblock)
   do iblock = 1, nblocks_clinic
      DO j = 2,jmt-1
         DO i = 2, imt-1
            k = min (kmt (i,j,iblock),kmt(i +1,j,iblock))
            IF (k /= 0) THEN
!               adv_vetiso (i,k,j,iblock) = - p5*dzr(k)*athkdf(i,j,k,iblock)*tmask(i,k,j,iblock) &
!                                    * tmask(i+1,k,j,iblock)* (K1 (i,k,j,3,iblock)   &
!                                      + K1 (i,k -1,j,3,iblock)) 
!               adv_vetiso (i,k,j,iblock) = -dzr(k)*athkdfx(i,j,k,iblock)&
               adv_vetiso (i,k,j,iblock) = -p5*dzr(k)*athkdfx(i,j,k,iblock)&
                                    * tmask(i,k,j,iblock)*tmask(i+1,k,j,iblock)&
                                    * (K1 (i,k,j,3,iblock)+K1 (i,k -1,j,3,iblock))
            END IF
         END DO
      END DO
  END DO
!----------------------------------------------------------------------
!     compute the vertical component of the isopycnal mixing velocity
!     at the center of the bottom face of the "t" cells, using the
!     continuity equation for the isopycnal mixing velocities
!-----------------------------------------------------------------------
 
 
 
!!$OMP PARALLEL DO PRIVATE (iblock)
   do iblock = 1, nblocks_clinic
      DO j = 3, jmt-2
         DO k = 1,km -1
            DO i = 3,imt-2
               adv_vbtiso (i,k,j,iblock) = DZP (k)* ( &
               (adv_vetiso (i,k,j,iblock)*htw(i+1,j,iblock) - adv_vetiso (i-1,k,j,iblock)*htw(i,j,iblock))* tarea_r(i,j,iblock) + &
               (adv_vntiso (i,k,j,iblock)*hts(i,j,iblock) - adv_vntiso (i,k,j -1,iblock)*hts(i,j-1,iblock))*tarea_r(i,j,iblock))     
            END DO
         END DO
      END DO
   END DO
 
 
!!$OMP PARALLEL DO PRIVATE (iblock,j)
   do iblock = 1, nblocks_clinic
      DO j = 3, jmt-2
         DO k = 1,km -1
            DO i = 3,imt-2
               adv_vbtiso (i,k,j,iblock) = (adv_vbtiso (i,k,j,iblock) + adv_vbtiso (i,k -1,j,iblock))*vit(i,j,k,iblock)
               !adv_vbtiso (i,k,j,iblock) = adv_vbtiso (i,k,j,iblock) + adv_vbtiso (i,k -1,j,iblock) !LPF20160804
            END DO
         END DO
      END DO
  END DO
 
 
!!$OMP PARALLEL DO PRIVATE (iblock,j,i)
   do iblock = 1, nblocks_clinic
      DO j = 3, jmt-2
         DO i = 3, imt-2
            adv_vbtiso (i,kmt (i,j,iblock),j,iblock) = c0
         END DO
      END DO
   END DO

#ifdef ISOOUT
!!$OMP PARALLEL DO PRIVATE (iblock,j,k,i)
     do iblock = 1, nblocks_clinic
       DO j = 1,jmt
          DO k = 1,km-1
             DO i = 2,imt-1
                vetiso(i,j,k,iblock) = p5*(adv_vetiso(i-1,k,j,iblock)+adv_vetiso(i,k,j,iblock))*vit(i,j,k,iblock)
             END DO
          END DO
       END DO  
      END DO 
!!$OMP END PARALLEL DO

!!$OMP PARALLEL DO PRIVATE (iblock,j,k,i)
     do iblock = 1, nblocks_clinic
       DO j = 2,jmt-1
          DO k = 1,km
             DO i = 2,imt-1
                vntiso (i,j,k,iblock) = p5*(adv_vntiso(i,k,j-1,iblock) + adv_vntiso(i,k,j,iblock))*vit(i,j,k,iblock)
             END DO
          END DO
        END DO
     END DO
!!$OMP END PARALLEL DO

!!$OMP PARALLEL DO PRIVATE (iblock,j,k,i)
     do iblock = 1, nblocks_clinic
       DO j = 1,jmt
          DO k = 1,km-1
             DO i = 1,imt
                vbtiso (i,j,k,iblock) = adv_vbtiso(i,k,j,iblock)*vit(i,j,k,iblock)
             END DO
          END DO
       END DO
     END DO
!!$OMP END PARALLEL DO

#endif


#ifdef LPFDIAG
!LPF20160515
          call outsumglobal(athkdf,outsum)
      if (mytid ==0) WRITE (16,*)'athkdf',outsum
          call outsumglobal(athkdf_param,outsum)
      if (mytid ==0) WRITE (16,*)'athkdf_param',outsum
     temp_tmp=0.0D0
    do iblock = 1, nblocks_clinic
      DO j = 1, jmt
         DO k = 1,km-1
            DO i = 1,imt
             temp_tmp (i,j,k,iblock)=rict_replace(i,j,k,iblock)
            END DO
          END DO
         END DO
        END DO
          call outsumglobal(temp_tmp,outsum)
      if (mytid ==0) WRITE (16,*)'rict_replace',outsum

     temp_tmp=0.0D0
    do iblock = 1, nblocks_clinic
      DO j = 1, jmt
            DO i = 1,imt
             temp_tmp (i,j,1,iblock)=rict_ref(i,j,iblock)
            END DO
          END DO
         END DO
          call outsumglobal(temp_tmp,outsum)
      if (mytid ==0) WRITE (16,*)'rict_ref',outsum
     

   temp_tmp=0.0D0
    do iblock = 1, nblocks_clinic
      DO j = 1, jmt
            DO i = 1,imt
             temp_tmp (i,j,1,iblock)=amld(i,j,iblock)
            END DO
          END DO
         END DO
          call outsumglobal(temp_tmp,outsum)
      if (mytid ==0) WRITE (16,*)'amld',outsum
!          call outsumglobal(athkdfx,outsum)
!      if (mytid ==0) WRITE (16,*)'athkdfx',outsum
!          call outsumglobal(athkdfy,outsum)
!      if (mytid ==0) WRITE (16,*)'athkdfy',outsum
!LPF20160515
#endif

      !deallocate(rict_ref,rict_replace) !move to tracer done

      RETURN
      END SUBROUTINE ISOADV
 
 
#else
      SUBROUTINE ISOADV
      RETURN
      END SUBROUTINE ISOADV
#endif 
