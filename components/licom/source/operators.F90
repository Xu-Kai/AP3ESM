!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

 module operators

!BOP
! !MODULE: operators
!
! !DESCRIPTION:
!  This module contains routines for various common
!  mathematical operators, including gradient, divergence, curl,
!  and a vertical velocity calculation.  Note that these routines
!  do {\em not} update ghost cells so results contain invalid
!  entries in some ghost cells.
!
! !REVISION HISTORY:
!  SVN:$Id: operators.F90 808 2006-04-28 17:06:38Z njn01 $

! !USES:

   use precision_mod
   use blocks
   use param_mod
   use domain
   use constant_mod
   use grid

   implicit none
   private
   save

! !PUBLIC MEMBER FUNCTIONS:
   public :: div,   &
             grad,  &
             zcurl

!EOP
!BOC
!EOC
!***********************************************************************

 contains

!***********************************************************************
!BOP
! !IROUTINE: div
! !INTERFACE:

 subroutine div(k,DIV_OUT,UX,UY,this_block)

! !DESCRIPTION:
!  This routine returns the divergence at T points (times the cell
!  area) of a vector field defined at U points.  The divergence
!  operator is defined as:
!
!  \begin{equation}
!  \nabla\cdot{\bf\rm u} = {1\over{\Delta_y}} \delta_x (\Delta_y u_x) +
!                          {1\over{\Delta_x}} \delta_y (\Delta_x u_y)
!  \end{equation}
!
! !REVISION HISTORY:
!  same as module

! !INPUT PARAMETERS

   integer (int_kind), intent(in) :: k   ! vertical level

   real (r8), dimension(nx_block,ny_block), intent(in) :: & 
      UX,UY              ! vector field defined at U-points

   type (block), intent(in) :: &
      this_block          ! block information for current block

! !OUTPUT PARAMETERS:

   real (r8), dimension(nx_block,ny_block), intent(out) :: & 
      DIV_OUT            ! divergence at T-points times cell area

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

   integer (int_kind) :: &
      i,j,               &! dummy counters
      bid                 ! local block id

!-----------------------------------------------------------------------
!
!  compute divergence using a 4 point stencil
!
!-----------------------------------------------------------------------

   bid = this_block%local_id

   DIV_OUT = c0

!    write(6,*) '+ operators.F90, 173: max(HTW),min(HTW), ', &
!                maxval(HTW(1:nx_block,1:ny_block,bid)),minval(HTW(1:nx_block,1:ny_block,bid))
!    write(6,*) '- operators.F90, 173: max(HTS),min(HTS), ', &
!                maxval(HTS(1:nx_block,1:ny_block,bid)),minval(HTS(1:nx_block,1:ny_block,bid))
!    write(6,*) '- operators.F90, 173: max(TAREA_R),min(TAREA_R), ', &
!                maxval(TAREA_R(1:nx_block,1:ny_block,bid)),minval(TAREA_R(1:nx_block,1:ny_block,bid))

   !!$OMP PARALLEL DO 
   do j=2,ny_block
   do i=1,nx_block-1
      if (k <= KMT(i,j,bid)) then
         DIV_OUT(i,j) = p5*((UX(i+1,j  )+ UX(i+1,j-1))*HTW(i+1,j,bid) -  &
                            (UX(i  ,j  )+ UX(i  ,j-1))*HTW(i  ,j,bid) +  &
                            (UY(i+1,j  )+ UY(i  ,j  ))*HTS(i  ,j,bid) -  &
                            (UY(i+1,j-1)+ UY(i  ,j-1))*HTS(i  ,j-1,bid))*TAREA_R(i,j,bid)
      endif
   end do
   end do

!-----------------------------------------------------------------------
!EOC

 end subroutine div

!***********************************************************************
!BOP
! !IROUTINE: grad
! !INTERFACE:

 subroutine grad(k, GRADX, GRADY, F, this_block)

! !DESCRIPTION:
!  This routine computes the gradient in i,j directions at U points
!  based on field defined at T-points.
!
!  \begin{eqnarray}
!     \nabla_x(F) &=& \delta_x F \\
!     \nabla_y(F) &=& \delta_y F
!  \end{eqnarray}
!
! !REVISION HISTORY:
!  same as module

! !INPUT PARAMETERS:

   integer (int_kind), intent(in) :: k  ! vertical level

   real (r8), dimension(nx_block,ny_block), intent(in) :: & 
      F                  ! field defined at T points

   type (block), intent(in) :: &
      this_block          ! block information for current block

! !OUTPUT PARAMETERS:

   real (r8), dimension(nx_block,ny_block), intent(out) :: & 
      GRADX,GRADY        ! gradient in (i,j) direction at U points

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

   integer (int_kind) :: &
      i,j,               &! dummy counters
      bid                 ! local block id

!-----------------------------------------------------------------------
!
!  compute gradient with a 4 point stencil
!
!-----------------------------------------------------------------------
      
   bid = this_block%local_id

   GRADX = c0
   GRADY = c0

!    write(6,*) '+ operators.F90, 173: max(DYUR),min(DYUR), ', &
!                maxval(DYUR(3:nx_block-2,3:ny_block-2,bid)),minval(DYUR(3:nx_block-2,3:ny_block-2,bid))
!    write(6,*) '- operators.F90, 173: max(DXUR),min(DXUR), ', &
!                maxval(DXUR(3:nx_block-2,3:ny_block-2,bid)),minval(DXUR(3:nx_block-2,3:ny_block-2,bid))

!    write(6,*) '+ operators.F90, 173: max(DYUR),min(DYUR), ', &
!                maxval(DYUR(2:nx_block-1,2:ny_block-1,bid)),minval(DYUR(2:nx_block-1,2:ny_block-1,bid))
!    write(6,*) '- operators.F90, 173: max(DXUR),min(DXUR), ', &
!                maxval(DXUR(2:nx_block-1,2:ny_block-1,bid)),minval(DXUR(2:nx_block-1,2:ny_block-1,bid))

   do j=1,ny_block-1
   do i=2,nx_block  
      if (k <= KMU(i,j,bid)) then
         GRADX(i,j) = DXUR(i,j,bid)*p5*(F(i  ,j+1) - F(i-1,j) - & 
                                        F(i-1,j+1) + F(i  ,j))
         GRADY(i,j) = DYUR(i,j,bid)*p5*(F(i  ,j+1) - F(i-1,j) + &
                                        F(i-1,j+1) - F(i  ,j))
      endif
   end do
   end do

!-----------------------------------------------------------------------
!EOC

 end subroutine grad

!***********************************************************************
!BOP
! !IROUTINE: zcurl
! !INTERFACE:

 subroutine zcurl(k,CURL,UX,UY,this_block)

! !DESCRIPTION:
!  This function returns the z-component of the curl of a vector
!  field defined at U points. 
!
!  \begin{equation}
!     \hat{\bf\rm z}\cdot\nabla\times{\bf\rm u} =
!     {1\over{\Delta_y}} \delta_x(\Delta_y u_y) -
!     {1\over{\Delta_x}} \delta_y(\Delta_x u_x)
!  \end{equation}
!
!  The result is actually multiplied by cell area and returned at 
!  T points. 
!
! !REVISION HISTORY:
!  same as module

! !INPUT PARAMETERS:

   integer (int_kind), intent(in) :: k  ! vertical level

   real (r8), dimension(nx_block,ny_block), intent(in) :: & 
      UX,UY              ! vector field defined at U-points

   type (block), intent(in) :: &
      this_block          ! block information for current block

! !OUTPUT PARAMETERS:

   real (r8), dimension(nx_block,ny_block), intent(out) :: & 
      CURL              ! z.curl(Ux,Uy) at T-points times cell area

   integer :: counter_i  ! for test

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

   integer (int_kind) :: &
      i,j,               &! dummy counters
      bid                 ! local block index

!-----------------------------------------------------------------------
!
!  compute curl using stencil similar to divergence 4 point stencil
!
!-----------------------------------------------------------------------

   bid = this_block%local_id

   CURL = c0

!    if(mytid == 100)then
!      write(6,*) 'operators.F90, 248: tlat, ', &
!                  maxval(tlat(3:nx_block-2,3:ny_block-2,bid))/DEGtoRAD,minval(tlat(3:nx_block-2,3:ny_block-2,bid))/DEGtoRAD
!      write(6,*) 'operators.F90, 248: tlon, ', &
!                  maxval(tlon(3:nx_block-2,3:ny_block-2,bid))/DEGtoRAD,minval(tlon(3:nx_block-2,3:ny_block-2,bid))/DEGtoRAD
!      write(6,*) 'operators.F90, 248: max(DYU),min(DYU), ', &
!                  maxval(DYU(1:nx_block,1:ny_block,bid)),minval(DYU(1:nx_block,1:ny_block,bid))
!      write(6,*) 'operators.F90, 248: max(DXU),min(DXU), ', &
!                  maxval(DXU(1:nx_block,1:ny_block,bid)),minval(DXU(1:nx_block,1:ny_block,bid))
!    end if


   do j=2,ny_block
   do i=1,nx_block-1
      !if(this_block%j_glob(j) .gt. 1200 .and. this_block%j_glob(j) .lt. 1400) then
!       if((this_block%j_glob(j).gt.0).and.(this_block%j_glob(j).lt.61).and.(this_block%i_glob(i).gt.2250).and.(this_block%i_glob(i).lt.2341))then
!           if(dyu(i,j,bid) .eq. c1) then
!               write(*,*) "DYU negetive val:",i,j,this_block%i_glob(i),this_block%j_glob(j),dyu(i,j,bid),my_task
!           endif 
!       endif


      if (k <= KMT(i,j,bid)) then
         CURL(i,j) = p5*(UY(i+1,   j)*DYU(i+1,   j, bid) +  & 
                         UY(i+1, j-1)*DYU(i+1, j-1, bid) -  &
                         UY(i,     j)*DYU(i,     j, bid) -  &
                         UY(i,   j-1)*DYU(i,   j-1, bid) -  &
                         UX(i+1,   j)*DXU(i+1,   j, bid) -  &
                         UX(i,     j)*DXU(i,     j, bid) +  &
                         UX(i+1, j-1)*DXU(i+1, j-1, bid) +  &
                         UX(i,   j-1)*DXU(i,   j-1, bid))*TAREA_R(i,j,bid)
      endif

      ! for test
!       if(CURL(i,j).eq.0.0)then
!           CURL(i,j) = 100.0
!       end if
   end do
   end do

! do counter_i = 574,574
!   if(mytid == counter_i)then
!     write(6,*) '==== mytid,lat_max,lat_min,lon_max,lon_min,max(DYU),min(DYU),max(CURL),min(CURL): ', mytid, &
!                 maxval(tlat(3:nx_block-2,3:ny_block-2,bid))/DEGtoRAD,minval(tlat(3:nx_block-2,3:ny_block-2,bid))/DEGtoRAD, &
!                 maxval(tlon(3:nx_block-2,3:ny_block-2,bid))/DEGtoRAD,minval(tlon(3:nx_block-2,3:ny_block-2,bid))/DEGtoRAD, &
!                 maxval(DYU(1:nx_block,1:ny_block,bid)),minval(DYU(1:nx_block,1:ny_block,bid)), &
!                 maxval(CURL(1:nx_block,1:ny_block)),minval(CURL(1:nx_block,1:ny_block))
!     write(6,*) '-----                            -----' 
!   end if
! end do

   !write(6,*) '==== operators.F90, line 248: max(CURL),min(CURL), ',maxval(abs(CURL(2:nx_block,2:ny_block))),minval(abs(CURL(2:nx_block,2:ny_block)))

!-----------------------------------------------------------------------
!EOC

 end subroutine zcurl


 end module operators

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
