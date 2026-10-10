!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

 module domain_size

!BOP
! !MODULE: domain_size
!
! !DESCRIPTION:
!  This module contains parameters for the global model domain size
!  decomposition block size.  It is used by the domain and block
!  modules for decomposing the model domain across processors.
!
! !REVISION HISTORY:
!  SVN:$Id$

! !USES:

   !use kinds_mod

   implicit none
   private
   save

! !DEFINED PARAMETERS:

   integer, parameter, public ::  &  ! model size parameters
      imt_global =  10800 ,&! extent of horizontal axis in i direction
      jmt_global =  6907  ,&! extent of horizontal axis in j direction
      km = 80          ,&! number of vertical levels
      ntra = 2            ! total number of tracers

   integer, parameter, public :: &
      block_size_x = BLCKX, &! size of block in 1st horizontal dimension
      block_size_y = BLCKY   ! size of block in 2nd horizontal dimension

   !*** The model will inform the user of the correct
   !*** values for theparameters below.  A value higher than
   !*** necessary will not cause the code to fail, but will
   !*** allocate more memory than is necessary.  A value that
   !*** is too low will cause the code to exit.  
   !*** A good initial guess is found using
   !*** max=(nx_global/block_size_x)*(ny_global/block_size_y)/
   !***         num_procs
 
   integer, parameter, public :: &
      max_blocks_clinic = MXBLCKS,  &! max number of blocks per processor
      max_blocks_tropic = MXBLCKS    !   in each distribution

!EOP
!BOC
!EOC
!***********************************************************************

 end module domain_size

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
