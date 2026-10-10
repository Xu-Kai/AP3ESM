module ocn_comp_esmf

#ifdef ESMF_INTERFACE
!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!BOP
! !MODULE: ocn_comp_esmf
! !INTERFACE:

! !DESCRIPTION:
!  This is the main driver for the Parallel Ocean Program (POP).
!
! !REVISION HISTORY:
!  SVN:$Id:
!
!
! !PUBLIC MEMBER FUNCTIONS:
  implicit none
  public :: ocn_register_esmf
  !public :: ocn_init_esmf
  !public :: ocn_run_esmf
  !public :: ocn_final_esmf
  SAVE
  private                              ! By default make data private

!
! ! PUBLIC DATA:
!
! !REVISION HISTORY:
! Author: Mariana Vertenstein, Fei Liu
!
!EOP
   

!=======================================================================

contains

!=======================================================================

  subroutine ocn_register_esmf(comp, rc)
    implicit none

    print *, "In ocn register routine"
    ! Register the callback routines.

    !call ESMF_GridCompSetEntryPoint(comp, ESMF_METHOD_INITIALIZE, &
    !     ocn_init_esmf, phase=1, rc=rc)
    !if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)
    !
    !call ESMF_GridCompSetEntryPoint(comp, ESMF_METHOD_RUN, &
    !     ocn_run_esmf, phase=1, rc=rc)
    !if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

    !call ESMF_GridCompSetEntryPoint(comp, ESMF_METHOD_FINALIZE, &
    !     ocn_final_esmf, phase=1, rc=rc)
    !if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

  end subroutine ocn_register_esmf

#endif

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
end module ocn_comp_esmf
