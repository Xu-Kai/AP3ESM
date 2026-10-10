!  CVS: $Id: accumm.F90,v 1.1.1.1 2004/04/29 06:22:39 lhl Exp $
!     =================
      SUBROUTINE CALAVG3HOUR
!     =================
 
#include <def-undef.h>
use param_mod
use dyn_mod
use tracer_mod
use output_mod
use pconst_mod
use isopyc_mod

use domain
use buf_mod,only:prec,evap,iceoff,roff
use forc_mod,only:psa
      IMPLICIT NONE
      integer :: iblock
      real :: Num_3hour


      !Num_3hour=float(Num_3houracc)
      Num_3hour=float(fstep)
      if(mytid==0) write(*,*)'Num_3hour,fstep',Num_3houracc,fstep


#if (defined AVG3HOUR)
   do iblock = 1, nblocks_clinic
      DO J = 1,JMT
         DO I = 1,IMT

           Z03hour (I,J,iblock)= Z03hour (I,J,iblock)/Num_3hour
           su3hour (I,J,iblock)= su3hour (I,J,iblock)/Num_3hour        !U windstress
           sv3hour (I,J,iblock)= sv3hour (I,J,iblock)/Num_3hour        !V windstress
           lthf3hour (I,J,iblock)= lthf3hour (I,J,iblock)/Num_3hour !latent flux
           sshf3hour (I,J,iblock)= sshf3hour (I,J,iblock)/Num_3hour !sensible flux
           lwv3hour (I,J,iblock)= lwv3hour (I,J,iblock)/Num_3hour    !long wave flux
           swv3hour (I,J,iblock)= swv3hour (I,J,iblock)/Num_3hour    !shortwave flux
!       fresh3hour(I,J,iblock)= fresh3hour (I,J,iblock) + fresh (I,J,iblock)/34.7*1.0e3*DZP(1)/OD0/Num_3hour !virtual slat flux
           fresh3hour(I,J,iblock)= fresh3hour (I,J,iblock)/Num_3hour !virtual slat flux
        DO K = 1,KM
              TS3hour (I,J,K,iblock)= TS3hour (I,J,K,iblock)/Num_3hour
              SS3hour (I,J,K,iblock)= SS3hour (I,J,K,iblock)/Num_3hour
              US3hour (I,J,K,iblock)= US3hour (I,J,K,iblock)/Num_3hour
              VS3hour (I,J,K,iblock)= VS3hour (I,J,K,iblock)/Num_3hour
              dens3hour (I,J,K,iblock)= dens3hour (I,J,K,iblock)/Num_3hour
              WS3hour (I,J,K,iblock)= WS3hour (I,J,K,iblock)/Num_3hour
        ENDDO

         END DO
      END DO
     END DO
#endif
#if (defined OBP_AVG3HOUR)
   do iblock = 1, nblocks_clinic
      DO J = 1,JMT
         DO I = 1,IMT
           pbo3hour (I,J,iblock)= pbo3hour (I,J,iblock)/Num_3hour
           pso3hour (I,J,iblock)= psa (I,J,iblock)
           ssh3hour (I,J,iblock)= ssh3hour (I,J,iblock)/Num_3hour
         END DO
      END DO
     END DO
#endif         
      RETURN
      END SUBROUTINE CALAVG3HOUR
 
 
