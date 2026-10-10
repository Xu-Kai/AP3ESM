      SUBROUTINE BAROTR
!     =================
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
use msg_mod
use domain
use grid
use blocks
use hmix_del2
use hmix_del4
use operators
use smuvh
use POP_GridHorzMod
use POP_HaloMod
use global_reductions
use gather_scatter
use distribution
use constant_mod
use var_check
#if ( defined TIDEFORCE)
use forc_mod, only: CD_InWave,depthstd
#endif
      IMPLICIT NONE

      INTEGER :: IEB,NC,IEB_LOOP
      real(r8)    :: gstar ,am_viv,fil_lat1,fil_lat2, ek0, maxz0, minz0
!     real(r8)  ::  ttt(imt_global,jmt_global)
      integer :: iblock,ii1,jj1,ii2,jj2, irec,ikm
      real(r8):: hduk(imt,jmt) , hdvk(imt,jmt), gradx(imt,jmt),grady(imt,jmt), div_out(imt,jmt)
      real(r8),target:: hdtk_b(imt,jmt)
      type(block):: this_block
      REAL(r8)   :: ts1,te1,clock0f
      real(r8):: LAMDA11

!      Define the threthold latitute for zonal smoother
!    allocate (buffer(imt_global,jmt_global))
       fil_lat1=65.0D0
       fil_lat2=65.0D0

!---------------------------------------------------------------------
!     INITIALIZE WORK ARRAYS
!---------------------------------------------------------------------
      wka=0
      work=0

      baro_loop : DO NC = 1,NBB
!!!!$OMP PARALLEL DO PRIVATE (iblock)
#ifdef SW_ATHREAD
#ifdef TIDEFORCE
      call tide(NC) 
#endif
    call athread_spawn(slave_barotr1_cu,barotr1_param)
    call athread_join()

#else

   do iblock = 1, nblocks_clinic
        call tgrid_to_ugrid(work(:,:,iblock),h0(:,:,iblock),iblock)
   end do
!---------------------------------------------------------------------
!     COMPUTING DH0
!---------------------------------------------------------------------
!!!$OMP PARALLEL DO PRIVATE (iblock,J) collapse(3)
     DO IBLOCK = 1, NBLOCKS_CLINIC
         DO J = 2, JMT
            DO I = 1,IMT-1
               WKA (I,J,1,IBLOCK)= UB (I,J,IBLOCK)* (DZPH (I,J,IBLOCK) + WORK (I,J,IBLOCK))
               WKA (I,J,2,IBLOCK)= VB (I,J,IBLOCK)* (DZPH (I,J,IBLOCK) + WORK (I,J,IBLOCK))
            END DO 
         END DO
     END DO
 
#endif

#ifdef SW_ATHREAD

    call athread_spawn(slave_barotr1_cu1,barotr1_1_param)
    call athread_join()

#else
    
    DO IBLOCK = 1, NBLOCKS_CLINIC
         this_block = get_block(blocks_clinic(iblock),iblock)
         call div(1,DIV_OUT,wka(:,:,1,iblock),wka(:,:,2,iblock),this_block)
!!!$OMP PARALLEL DO PRIVATE (this_block,J,div_out) collapse(2)
         DO J = 3, jmt-2
            DO I = 3,imt-2
               WORK (I,J,IBLOCK)=VIT(I,J,1,IBLOCK)*(-1)*div_out(i,j)*P25
             END DO 
          ENDDO
    END DO

#endif
!
         call POP_HaloUpdate2DR8_opt(work , POP_haloClinic, POP_gridHorzLocCenter,&
                       POP_fieldKindScalar, errorCode, fillValue = 0.0_r8)
    
#ifdef SW_ATHREAD

    call athread_spawn(slave_barotr2_cu1,barotr2_1_param)
    call athread_join()

#else
    DO IBLOCK = 1, NBLOCKS_CLINIC
         DO J = 1, jmt
            DO I = 1,imt
               H0 (I,J,IBLOCK)= H0P (I,J,IBLOCK) + WORK (I,J,IBLOCK) * DTB
#ifdef TIDEFORCE
               h0_tide(I,J,IBLOCK) = H0(I,J,IBLOCK)*BETA_TIDE-elevat(I,J,IBLOCK)
#endif 
            END DO
         END DO
     END DO

#endif
!
!---------------------------------------------------------------------
!     COMPUTE THE "ARTIFICIAL" HORIZONTAL VISCOSITY
!---------------------------------------------------------------------

#ifdef SW_ATHREAD
    !am_factor = 1.0_r8

#ifdef BIHAR

    call athread_spawn(slave_bhdel4_div,bhdel4_div_param)
    call athread_join()

    call athread_spawn(slave_bhdel4_zcurl,bhdel4_zcurl_param)
    call athread_join()

    call athread_spawn(slave_bhdel4_grad,bhdel4_grad_param)
    call athread_join()


    call athread_spawn(slave_barotr2_udel4_0,barotr2_udel4_0_param)
    call athread_join()

    call athread_spawn(slave_barotr2_udel4_1,barotr2_udel4_1_param)
    call athread_join()

!    DO IBLOCK = 1, NBLOCKS_CLINIC
!         DO J = 3, jmt-2
!            DO I = 3, imt-2
!         if (my_task == 10) then
!             if( WKA (I,J,5,IBLOCK) .ne. 0.0) then
!                !write(*,*) "i,j,D2UK,D2VK", i,j,D2UK(i,j),D2VK(i,j)
!                !write(*,*) "i,j,real am_factor", i,j,am_factor_d(i,j,IBLOCK)
!                write(*,*) "i,j,hduk,hdvk", i,j, WKA (I,J,5,IBLOCK), WKA (I,J,6,IBLOCK)
!             end if
!          endif
!      end do
!      end do
!    end do


#else
    call athread_spawn(slave_barotr2_cu,barotr2_param)
    call athread_join()
#endif

#else
!sw_athread

#if ( defined SMAG1)
#else
#if (defined BIHAR)
     DO IBLOCK = 1, NBLOCKS_CLINIC
         this_block = get_block(blocks_clinic(iblock),iblock)
         call hdiffu_del4(1, HDUK, HDVK, ubp(:,:,iblock), vbp(:,:,iblock), this_block)
!!!$OMP PARALLEL DO PRIVATE (this_block, hduk,hdvk) collapse(2)
         DO J = 3, jmt-2
            DO I = 3, imt-2
               WKA (I,J,5,IBLOCK)= hduk(i,j)
               WKA (I,J,6,IBLOCK)= hdvk(i,j)
            END DO
         END DO
     END DO
!
!!!!!!
#else
     DO IBLOCK = 1, NBLOCKS_CLINIC
         this_block = get_block(blocks_clinic(iblock),iblock)
         call hdiffu_del2(1, HDUK, HDVK, ubp(:,:,iblock), vbp(:,:,iblock), this_block)
!!!$OMP PARALLEL DO PRIVATE ( this_block, hduk, hdvk) collapse(2)
         DO J = 3, jmt-2
            DO I = 3, imt-2
               WKA (I,J,5,IBLOCK)= hduk(i,j)
               WKA (I,J,6,IBLOCK)= hdvk(i,j)
            END DO
         END DO
    END DO

#endif
#endif

!      if (mytid == 97) then
!         buffer_real4(:,:)= hduk(:,:)
!         irec=(NC-1)*8+1
!         write(25,rec=irec) ((buffer_real4(i,j),i=3,imt-2),j=3,jmt-2)
!      end if
!
!     + (g'-1)g*dH/dr
!---------------------------------------------------------------------

     DO IBLOCK = 1, NBLOCKS_CLINIC
         this_block = get_block(blocks_clinic(iblock),iblock)
#ifdef TIDEFORCE
         call grad(1, GRADX, GRADY, h0_tide(:,:,iblock), this_block) 
#else
         call grad(1, GRADX, GRADY, H0(:,:,iblock), this_block)
#endif
!!!$OMP PARALLEL DO PRIVATE (this_block, gradx, grady, gstar)  collapse(2)
         DO J = 3, jmt-2
            DO I = 3, imt-2
               gstar=(WGP (I,J,IBLOCK) -1.0)*G 
               WKA (I,J,1,IBLOCK) = WKA (I,J,5,IBLOCK) + gstar*GRADX(I,J)
               WKA (I,J,2,IBLOCK) = WKA (I,J,6,IBLOCK) + gstar*GRADY(I,J)
            END DO
         END DO
     END DO

!Yu
!      if (mytid == 97) then
!         buffer_real4(:,:)= wka(:,:,1,1)
!         irec=(NC-1)*8+2
!         write(25,rec=irec) ((buffer_real4(i,j),i=3,imt-2),j=3,jmt-2)
!      end if
!!!!$OMP PARALLEL DO PRIVATE (iblock)
   do iblock = 1, nblocks_clinic
        call tgrid_to_ugrid(work(:,:,iblock),h0(:,:,iblock),iblock)
   end do
!
!---------------------------------------------------------------------
!     COMPUTING DU & DV
!---------------------------------------------------------------------
!!!$OMP PARALLEL DO PRIVATE (IBLOCK,J)  collapse(3)
     DO IBLOCK = 1, NBLOCKS_CLINIC
         DO J = 3, jmt-2
            DO I = 3, imt-2
               WKA (I,J,1,IBLOCK)= VIV (I,J,1,IBLOCK)* ( WKA (I,J,1,IBLOCK) + DLUB (I,J,IBLOCK)     &
                              - FCOR(I,J,Iblock)* VBP (I,J,IBLOCK) + &
               PAX (I,J,IBLOCK) + PXB (I,J,IBLOCK) - WORK (I,J,IBLOCK)* WHX (I,J,IBLOCK) )
               WKA (I,J,2,IBLOCK)= VIV (I,J,1,IBLOCK)* ( WKA (I,J,2,IBLOCK) + DLVB (I,J,IBLOCK)     &
                              + FCOR(I,J,Iblock)* UBP (I,J,IBLOCK) + &
               PAY (I,J,IBLOCK) + PYB (I,J,IBLOCK) - WORK (I,J,IBLOCK)* WHY (I,J,IBLOCK) )
            END DO
         END DO
     END DO
!---------------------------------------------------------------------
!     CORIOLIS ADJUSTMENT
!---------------------------------------------------------------------

!!!$OMP PARALLEL DO PRIVATE (IBLOCK,J)  collapse(3)
     DO IBLOCK = 1, NBLOCKS_CLINIC
     DO J = 3, jmt-2
     DO I = 3, imt-2
        WKA (I,J,3,IBLOCK)= EBEA(I,J,IBLOCK)* WKA (I,J,1,IBLOCK) - EBEB(I,J,IBLOCK)* WKA (I,J,2,IBLOCK)
        WKA (I,J,4,IBLOCK)= EBEA(I,J,IBLOCK)* WKA (I,J,2,IBLOCK) + EBEB(I,J,IBLOCK)* WKA (I,J,1,IBLOCK)
     END DO
     END DO
     END DO
#endif

#if defined(TIDEFORCE) || defined(TIDEFORCEGW)
# ifdef SW_ATHREAD
   call athread_spawn(slave_barotr2_tide,barotr2_tide_param)
   call athread_join()   
# else
#  if ( defined TIDEFORCE )
          DO IBLOCK = 1, NBLOCKS_CLINIC 
          DO J = 3, jmt-2 
          DO I = 3, imt-2 
             WKA(I,J,3,IBLOCK)=WKA(I,J,3,IBLOCK)-(UBP(I,J,IBLOCK)**2+VBP(I,J,IBLOCK)**2)**0.5&
                       *UBP(I,J,IBLOCK)/(1/OHBU(I,J,IBLOCK)+H0(I,J,IBLOCK))*CD_TIDE(I,J,IBLOCK) 
             WKA(I,J,4,IBLOCK)=WKA(I,J,4,IBLOCK)-(UBP(I,J,IBLOCK)**2+VBP(I,J,IBLOCK)**2)**0.5&
                       *VBP(I,J,IBLOCK)/(1/OHBU(I,J,IBLOCK)+H0(I,J,IBLOCK))*CD_TIDE(I,J,IBLOCK) 

          END DO 
          END DO 
          END DO 
#  endif
#  if (defined TIDEFORCEGW)
         DO IBLOCK = 1,NBLOCKS_CLINIC 
         DO J = 3,JMT-2 
         DO I = 3,IMT-2 
             CD_InWave(I,J,IBLOCK)=tidetune*(2*pi/10000.0)*depthstd(i,j,iblock)*depthstd(i,j,iblock)*SQRT(abs(ric(i,j,int(kmu(i,j,iblock))-1,iblock)))
     !   Based  on M.C. Buijsman et al 2015, if the depth is less than 1000m,CD_INwave should be set to zero. 
         if((viv(i,j,1,iblock).gt.0.5))then
           if(1/OHBU(I,J,iblock).lt.1000)THEN 
             CD_InWave(I,J,iblock)=0.0D0 
         endif
         endif
         WKA(I,J,3,IBLOCK)=WKA(I,J,3,IBLOCK)-CD_InWave(I,J,IBLOCK)*UBP(I,J,IBLOCK)/((1/OHBU(I,J,IBLOCK)+H0(I,J,IBLOCK))*D0) 
         WKA(I,J,4,IBLOCK)=WKA(I,J,4,IBLOCK)-CD_InWave(I,J,IBLOCK)*VBP(I,J,IBLOCK)/((1/OHBU(I,J,IBLOCK)+H0(I,J,IBLOCK))*D0)

         END DO 
         END DO 
         END DO 
#  endif
# endif
!sw_athread
#endif 

    ikm = 2
    call POP_HaloUpdate3DR8_opt(wka(:,:,3:4,:), POP_haloClinic, POP_gridHorzLocSWcorner , &
                        POP_fieldKindVector, ikm, errorCode, fillValue = 0.0_r8)
!
!    call POP_HaloUpdate(wka(:,:,4,:), POP_haloClinic, POP_gridHorzLocSWcorner , &
!                        POP_fieldKindVector, errorCode, fillValue = 0.0_r8)

#ifdef SW_ATHREAD

    call athread_spawn(slave_barotr3_cu,barotr3_param)
    call athread_join()

#else
!!!$OMP PARALLEL DO PRIVATE (IBLOCK,J,I)  collapse(3)
    DO IBLOCK = 1, NBLOCKS_CLINIC
         DO J = 1, jmt
            DO I = 1,imt
               UB (I,J,IBLOCK)= UBP (I,J,IBLOCK) + WKA (I,J,3,IBLOCK)* DTB
               VB (I,J,IBLOCK)= VBP (I,J,IBLOCK) + WKA (I,J,4,IBLOCK)* DTB
            END DO
         END DO
     END DO
     call data_check(UBP, info_var="UBP")
     call data_check(VBP, info_var="VBP")
!---------------------------------------------------------------------
!     COMPUTING DH0
!---------------------------------------------------------------------
!!!$OMP PARALLEL DO PRIVATE (IBLOCK,J,I)  collapse(3)
     DO IBLOCK = 1, NBLOCKS_CLINIC
         DO J = 2, JMT
            DO I = 1,IMT-1
               WKA (I,J,1,IBLOCK)= UB (I,J,IBLOCK)* (DZPH (I,J,IBLOCK) + WORK (I,J,IBLOCK))
               WKA (I,J,2,IBLOCK)= VB (I,J,IBLOCK)* (DZPH (I,J,IBLOCK) + WORK (I,J,IBLOCK))
            END DO
         END DO
     END DO

#endif

!!!!$OMP PARALLEL DO PRIVATE (IBLOCK,this_block,hdtk_b,div_out) 
    DO IBLOCK = 1, NBLOCKS_CLINIC
         this_block = get_block(blocks_clinic(iblock),iblock)

#ifdef SW_ATHREAD
#ifdef BIHAR
    barotr3_tdel4_1_param%d_nc           = nc
    call athread_spawn(slave_barotr3_tdel4_1,barotr3_tdel4_1_param)
    call athread_join()
#endif
    barotr3_1_param%d_nc           = nc
    call athread_spawn(slave_barotr3_cu1,barotr3_1_param)
    call athread_join()
#else
!
!2020 YYQ
         !if ( mod(nc,2) == 0 ) then
         if ( mod(nc,4) == 0 ) then
!2020 YYQ
#ifdef  BIHAR
!2020 YYQ
             !call hdifft_del4(1,hdtk_b,h0p(:,:,iblock),this_block)
             call hdifft_del4(1,dt2k_b,hdtk_b,h0p(:,:,iblock),this_block)
!2020 YYQ
#else
             call hdifft_del2(1,hdtk_b,h0p(:,:,iblock),this_block)
#endif
         else
             hdtk_b = c0
         end if

         call div(1,DIV_OUT,wka(:,:,1,iblock),wka(:,:,2,iblock),this_block)

         DO J = 2, jmt-2
            DO I = 3,imt-2
               !WORK (I,J,IBLOCK)=VIT(I,J,1,IBLOCK)*(hdtk_b(i,j) - div_out(i,j))
               WORK (I,J,IBLOCK)=VIT(I,J,1,IBLOCK)*(hdtk_b(i,j)*1.0D0 - div_out(i,j))
             END DO
          ENDDO

#endif
    END DO


!
!---------------------------------------------------------------------
!     PREDICTING VB , UB & H0
!---------------------------------------------------------------------
         call POP_HaloUpdate2DR8_opt(work , POP_haloClinic, POP_gridHorzLocCenter,&
                       POP_fieldKindScalar, errorCode, fillValue = 0.0_r8)
!
!!!$OMP PARALLEL DO PRIVATE (IBLOCK,J,I)  collapse(3)
#ifdef SW_ATHREAD
    call athread_spawn(slave_barotr4_cu,barotr4_param)
    call athread_join()
    ISB = ISB +1
#else

    if (enoi_iau ) then
      LAMDA11=1.0D0/(iau_gamma*86400.D0)
    endif
    DO IBLOCK = 1, NBLOCKS_CLINIC
         DO J = 1, jmt
            DO I = 1,imt
                if (enoi_iau ) then
                    H0 (I,J,IBLOCK)= H0P (I,J,IBLOCK) + WORK (I,J,IBLOCK) * DTB&
                            + VIT (I,J,1,iblock)*h0_iau(I,J,iblock)*LAMDA11 !*DTB
                else
                    H0 (I,J,IBLOCK)= H0P (I,J,IBLOCK) + WORK (I,J,IBLOCK) * DTB
                end if

            END DO
         END DO
     END DO

            ISB = ISB +1
            ubp=ub
            vbp=vb
            h0p=h0
!!!$OMP PARALLEL DO PRIVATE (IBLOCK,J,I)  collapse(3)
    DO IBLOCK = 1, NBLOCKS_CLINIC
            DO J = JST,JET ! Dec. 4, 2002, Yongqiang YU
               DO I = 1,IMT
                  H0F (I,J,IBLOCK) = H0F (I,J,IBLOCK) + H0 (I,J,IBLOCK)
                  H0BF (I,J,IBLOCK) = H0BF (I,J,IBLOCK) + H0 (I,J,IBLOCK)
               END DO
            END DO
    END DO
#endif
!    if (maxval(abs(ub)) > 5.0D0 ) then
!         write(6,*) "UBMAX", isc, isb, mytid,maxval(abs(ub)),maxloc(abs(ub))
!         write(6,'(1x,A5,1x,3I8,1x,F24.20,1x,4I4)') "UBMAX", isc, isb, mytid,maxval(abs(ub)),maxloc(abs(ub))
!ZWP20170209
!    end if
!    if (maxval(abs(vb)) > 5.0D0 ) then
!         write(6,*) "VBMAX", isc, isb, mytid,maxval(abs(vb)),maxloc(abs(vb))
!         write(6,'(1x,A5,1x,3I8,1x,F24.20,1x,4I4)') "UBMAX", isc, isb, mytid,maxval(abs(ub)),maxloc(abs(ub))
!ZWP20170209
!    end if

      END DO baro_loop
!    deallocate (buffer)

      RETURN
      END SUBROUTINE BAROTR


