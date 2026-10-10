!  CVS: $Id: icesnow.F90,v 1.1.1.1 2004/04/29 06:22:39 lhl Exp $
!     ==================
      SUBROUTINE ICESNOW
!     ==================
!     Sea Ice Model
#include <def-undef.h> 
use precision_mod
use param_mod
use pconst_mod
use tracer_mod
use domain
use grid, only :  kmt
use shr_const_mod
      IMPLICIT NONE
!
      real(r8) :: sal_ocn, sal_ice, tdiff, heat_ice_fusion
      real(r8) ::  t_mix, s_mix,ek0,ek1
      integer :: iblock
!
      heat_ice_fusion = SHR_CONST_LATICE         !J/Kg
      sal_ocn=34.7D0                     !Reference salinity for ocean
      sal_ice= 4.0D0                     !Reference salinity for sea ice
!
!------------------------------------------------------------
!  if SST exceed -1.8C to restore it to -1.8C
!------------------------------------------------------------
 
      tdiff = 0.0D0

      DO IBLOCK = 1, NBLOCKS_CLINIC
#ifdef SW_ATHREAD
!$OMP TARGET PARALLEL DO PRIVATE (j,k,i,tdiff) collapse(2)
#endif
      JJJ : DO J = 1, JMT
      III : DO I = 1,IMT
#ifdef COUP
      KKK : DO K = 1, 15
#else
      KKK : DO K = 1, 1
#endif
            IF (AT (I,J,K,1,IBLOCK) < TBICE .and. KMT(I,J,IBLOCK) > 0 ) THEN
               tdiff        = TBICE- AT(i,j,k,1,IBLOCK)
               licomqice (I,J  ,IBLOCK) = licomqice(i,j,IBLOCK)+tdiff*dzp(k)/dzp(1)
               at(i,j,k,2,IBLOCK)=at(i,j,k,2,IBLOCK)+tdiff*(sal_ocn-sal_ice)  &
                          *CP/heat_ice_fusion*0.001D0
               AT (I,J,k,1,iblock) = TBICE
            END IF
      END DO KKK
      END DO  III
      END DO JJJ
#ifdef SW_ATHREAD
!$OMP END TARGET PARALLEL DO
#endif      
#ifdef SW_ATHREAD
!$OMP TARGET PARALLEL DO PRIVATE (j,i,tdiff) collapse(2)
#endif
      DO J=1,JMT
      DO I=1,IMT
         IF (licomqice(I,J,IBLOCK) > 0.0 .and. at(i,j,1,1,IBLOCK) > TBICE) THEN
              tdiff=min((at(i,j,1,1,iblock)-tbice),licomqice(i,j,iblock))
              licomqice(i,j,IBLOCK)=licomqice(i,j,IBLOCK)-tdiff
              at(i,j,1,1,IBLOCK)=at(i,j,1,1,IBLOCK)-tdiff
              at(i,j,1,2,IBLOCK)=at(i,j,1,2,IBLOCK)-tdiff*(sal_ocn-sal_ice)   &
                         *CP/heat_ice_fusion*0.001D0
          END IF
      END DO
      END DO 
#ifdef SW_ATHREAD
!$OMP END TARGET PARALLEL DO
#endif 
      DO K=1 ,KM
      DO J=1,JMT
      DO I=1,IMT
         IF (AT(I,J,K,1,IBLOCK) < TBICE) AT(I,J,K,1,IBLOCK)=TBICE
      ENDDO
      ENDDO
      ENDDO
!
   END DO

      RETURN
      END SUBROUTINE ICESNOW
 
 
