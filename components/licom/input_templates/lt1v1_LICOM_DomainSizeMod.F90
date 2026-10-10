!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

 module LICOM_DomainSizeMod

!BOP
! !MODULE: LICOM_DomainSizeMod
!
! !DESCRIPTION:
!  This module contains parameters for the global model domain size
!  decomposition block size.  It is used by the domain and block
!  modules for decomposing the model domain across processors.
!
! !REVISION HISTORY:
!  SVN:$Id$
!  2006-08-14: Phil Jones
!              New domain size module following new naming conventions

! !USES:

   use LICOM_KindsMod

   implicit none
   private
   save

! !DEFINED PARAMETERS:

   integer (LICOM_i4), parameter, public ::  &  ! model size parameters
      LICOM_nxGlobal = 360 ,&! extent of horizontal axis in i direction
      LICOM_nyGlobal = 218 ,&! extent of horizontal axis in j direction
      LICOM_km = 30          ,&! number of vertical levels
      LICOM_nt = 2            ! total number of tracers

   integer (LICOM_i4), parameter, public :: &
      LICOM_blockSizeX = BLCKX, &! size of block in first  horizontal dimension
      LICOM_blockSizeY = BLCKY   ! size of block in second horizontal dimension

   !*** The model will inform the user of the correct
   !*** values for the parameters below.  A value higher than
   !*** necessary will not cause the code to fail, but will
   !*** allocate more memory than is necessary.  A value that
   !*** is too low will cause the code to exit.  
   !*** A good initial guess is found using
   !*** max=(nx_global/block_size_x)*(ny_global/block_size_y)/
   !***         num_procs
 
   integer (LICOM_i4), parameter, public :: &
      LICOM_maxBlocksClinic = MXBLCKS,  &! max number of blocks per processor
      LICOM_maxBlocksTropic = MXBLCKS    !   in each distribution

!EOP
!BOC
!EOC
!***********************************************************************

 end module LICOM_DomainSizeMod

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
