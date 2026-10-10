!  CVS: $Id: mm00.F90,v 1.1.1.1 2004/04/29 06:22:39 lhl Exp $
!     ===================
      SUBROUTINE AVG3HOURINI
!     ===================
 
#include <def-undef.h>
use param_mod
use pconst_mod
use constant_mod
use output_mod
      IMPLICIT NONE

#if (defined AVG3HOUR)

           TS3HOUR = 0.0
           SS3HOUR = 0.0
           US3HOUR = 0.0
           VS3HOUR = 0.0
           WS3HOUR = 0.0
           dens3hour = 0.0
           Z03HOUR    = 0.0
           SU3HOUR    = 0.0   !U windstress
           SV3HOUR    = 0.0   !V windstress
           LTHF3HOUR  = 0.0 !latent flux
           SSHF3HOUR  = 0.0 !sensible flux
           LWV3HOUR   = 0.0  !long wave flux
           SWV3HOUR   = 0.0  !shortwave flux
           FRESH3HOUR =0.0
#endif
#if (defined OBP_AVG3HOUR)
           pbo3hour   = 0.0
           pso3hour = 0.0
           ssh3hour = 0.0
#endif
      RETURN
      END SUBROUTINE AVG3HOURINI
 
 
