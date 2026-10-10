!  CVS: $Id: accumm.F90,v 1.1.1.1 2004/04/29 06:22:39 lhl Exp $
!     =================
      SUBROUTINE ACCUMM3HOUR
!     =================
 
#include <def-undef.h>
use param_mod
use dyn_mod
use tracer_mod
use output_mod
use pconst_mod
use isopyc_mod
#if ( defined TIDEMIX )
use forc_mod, only: su,sv,lthf,sshf,lwv,swv,fresh,runoff,seaice,wave_dis,psa
#else
use forc_mod, only: su,sv,lthf,sshf,lwv,swv,fresh,runoff,seaice,psa
#endif
use domain
use buf_mod,only:prec,evap,iceoff,roff,patm
!use grid, only: FCOR, FCORT 
      IMPLICIT NONE
      integer :: iblock
      real :: Num_op2 


      Num_op2= 1 !float(Num_dailyoutput)


#if (defined AVG3HOUR)
   do iblock = 1, nblocks_clinic
      DO J = 1,JMT
         DO I = 1,IMT
           Z03hour (I,J,iblock)= Z03hour (I,J,iblock) + H0(I,J,iblock)/Num_op2
           su3hour (I,J,iblock)= su3hour (I,J,iblock) + su(I,J,iblock)/Num_op2        !U windstress
           sv3hour (I,J,iblock)= sv3hour (I,J,iblock) - sv(I,J,iblock)/Num_op2        !V windstress
           lthf3hour (I,J,iblock)= lthf3hour (I,J,iblock) + lthf (I,J,iblock)/Num_op2 !latent flux
           sshf3hour (I,J,iblock)= sshf3hour (I,J,iblock) + sshf (I,J,iblock)/Num_op2 !sensible flux
           lwv3hour (I,J,iblock)= lwv3hour (I,J,iblock) + lwv (I,J,iblock)/Num_op2    !long wave flux
           swv3hour (I,J,iblock)= swv3hour (I,J,iblock) + swv (I,J,iblock)/Num_op2    !shortwave flux
!           fresh3hour(I,J,iblock)= fresh3hour (I,J,iblock) + fresh (I,J,iblock)/34.7*1.0e3*DZP(1)/OD0/Num_op2 !virtual slat flux
           fresh3hour(I,J,iblock)= fresh3hour (I,J,iblock) + fresh (I,J,iblock)/Num_op2 !virtual slat flux

        DO K = 1,KM
              TS3hour (I,J,K,iblock)= TS3hour (I,J,K,iblock) + AT (I,J,K,1,iblock)*vit(i,j,k,iblock)/Num_op2
              SS3hour (I,J,K,iblock)= SS3hour (I,J,K,iblock) + (AT (I,J,K,2,iblock)*1000.+35.)*vit(i,j,k,iblock)/Num_op2
              US3hour (I,J,K,iblock)= US3hour (I,J,K,iblock) + U (I,J,K,iblock)/Num_op2
              VS3hour (I,J,K,iblock)= VS3hour (I,J,K,iblock) - V(I,J,K,iblock)/Num_op2
              dens3hour (I,J,K,iblock)= dens3hour (I,J,K,iblock) + PDENSITY(I,J,K,iblock)/Num_op2
              WS3hour (I,J,K,iblock)= WS3hour (I,J,K,iblock) + WS(I,J,K,iblock)/Num_op2
        ENDDO

           WS3hour (I,J,KM+1,iblock)= WS3hour (I,J,KM+1,iblock) + WS(I,J,KM+1,iblock)/Num_op2

         END DO
      END DO
     END DO
#endif
#if (defined OBP_AVG3HOUR)
!yuzp20190411
      pbo_save = psa + G*(H0 + DZP(1)) * pdensity(:,:,1,:)

      DO k = 2, KM
      where (vit(:,:,k,:) > 0.5D0)
        pbo_save = pbo_save + G * DZP(k) * pdensity(:,:,k,:)
      elsewhere 
        pbo_save = pbo_save 
      end where
      END DO

      pbo3hour = pbo3hour + pbo_save / Num_op2    
      ssh3hour = ssh3hour + h0 / Num_op2
#endif
      RETURN
      END SUBROUTINE ACCUMM3HOUR
 
 
