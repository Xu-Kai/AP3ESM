!  CVS: $Id: tracer.F90,v 1.1.1.1 2004/04/29 06:22:39 lhl Exp $
!     =================
      SUBROUTINE TRACER
!     =================
#include <def-undef.h>
#ifdef SW_ATHREAD
use sunway_mod
#endif
use precision_mod 
use param_mod
use pconst_mod
use tracer_mod
use work_mod
use dyn_mod
use isopyc_mod
use forc_mod
use pmix_mod
use msg_mod
use smuvh
use advection
use blocks
use domain
use LICOM_Error_mod
use gather_scatter
use distribution
use buf_mod
!use omp_lib
#ifdef BIHAR
use hmix_del4
#else
use hmix_del2
#endif
use time_parameter
      IMPLICIT NONE
 
#ifdef SW_ATHREAD
interface
   subroutine tracer_sw() bind(c,name='tracer_sw')
   end subroutine tracer_sw

   subroutine tracer5_cu_slave(n) bind(c,name='tracer5_cu_slave')
            use iso_c_binding, only: c_int
            integer(c_int), value :: n
   end subroutine tracer5_cu_slave

   subroutine tracer6_cu_slave(n,e_norm2) bind(c,name='tracer6_cu_slave')
            use iso_c_binding, only: c_int,c_double
            integer(c_int), value :: n
            real(c_double), value :: e_norm2
   end subroutine tracer6_cu_slave

   subroutine tracer7_cu_slave(n) bind(c,name='tracer7_cu_slave')
            use iso_c_binding, only: c_int
            integer(c_int), value :: n
   end subroutine tracer7_cu_slave

!   subroutine upwell_sw() bind(c,name='upwell_sw')
!   end subroutine upwell_sw

!   subroutine invtrit_sw(h_c2dtts) bind(c,name='invtrit_sw')
!            use iso_c_binding, only: c_double
!            real(c_double), value :: h_c2dtts
!   end subroutine invtrit_sw
end interface

    REAL(r8)   :: ts1,te1

#endif
      integer     :: n2, iblock,kt, irec
      REAL(r8)    :: AIDIF,AA,FAW,FIW,ALF,RNCC,ABC,fil_lat1,fil_lat2,ek0 ,ek1
      REAL (r8)   :: DT2K(imt,jmt)
      real(r8)    :: clock0f
#ifdef SSSNORM
!lhl20130913
      REAL(r8)    :: ERR_norm1,ERR_norm2
!lhl20130913
      real (r8),dimension(:,:,:),allocatable::temp11
      real (r8),dimension(:,:),allocatable::temp12 
#endif

!Xiao Chan (Hereinafter XC for short)
      real(r8)    :: LAMDA(imt,jmt,km,max_blocks_clinic),wt1,wt2,adv_x1,adv_x2
      real(r8)    :: LAMDA1(km)
      type (block) :: this_block          ! block information for current block
#ifdef LPFDIAG
      REAL(r8)    :: outsum 
#endif   

      real(r8) :: tmp_sum(600000)
      integer :: rank_size



!XC
!      integer :: t1,t2,t3,t4,t5,t6,t7,t8
!      integer :: b1,b2,b3,b4
 
!---------------------------------------------------------------------
!     SET LOCAL CONSTANT
!---------------------------------------------------------------------
!      deallocate(dlu,dlv,gg,ric,rict)
 
!      allocate(stf(imt,jmt,max_blocks_clinic),tf(imt,jmt,km,max_blocks_clinic))
!      allocate(wkb(imt,jmt,km,max_blocks_clinic),wkc(imt,jmt,km,max_blocks_clinic),wkd(imt,jmt,km,max_blocks_clinic))
!lhl      AIDIF = 0.5*FLOAT(ISOP)
!
!---------------------------------------------------------------------
!      Define the threthold latitute for zonal smoother
       fil_lat1=63.0D0
       fil_lat2=63.0D0
!      GAMMA = 0.0D0

      AIDIF = 0.5D0
 
      LAMDA=1.0D0/(15.D0*86400.D0)
     do k=1,45
      LAMDA1(k)=1.0D0/(5.D0*86400.D0)
     enddo
     do k=46,km
      LAMDA1(k)=1.0D0/(5.D0*86400.D0)/(k*k/5.0/5.0)
     enddo

      RNCC = 1.0D0/ FLOAT (NCC)
!
!
      if ( adv_tracer(1:8) == 'centered' )  then
      IF (IST >= 1)THEN
         C2DTTS = DTS *2.0D0
         AA = 0.5D0
      ELSE
         C2DTTS = DTS
         AA = 0.0D0
      END IF
      else if ( adv_tracer(1:5) == 'tspas' .or. adv_tracer(1:3)=='prm' ) then
         C2DTTS = DTS
         AA = 0.5D0
      else
          WRITE (*,*)'adv_tracer1:',trim(adv_tracer)
        call exit_licom(sigAbort,'The false advection option for tracer')
      end if
 
#ifdef SW_ATHREAD
#ifdef ENABLE_GPLT
    call t_startf("tracer_sw"//CHAR(0))
#endif
      call tracer_sw
#ifdef ENABLE_GPLT
    call t_stopf("tracer_sw"//CHAR(0))
#endif
#else
       ! DO iblock = 1, nblocks_clinic
 !    DO K = 1, 23  
 !    DO J = 1, JMT
 !       DO I = 1,IMT
 !          if ( tlat(i,j,iblock) > 45.0D0*degtorad .and. tlat(i,j,iblock) < &
 !               75.D0*degtorad .and. tlon(i,j,iblock) > -100.0D0*degtorad .and.  &
 !               tlon(i,j,iblock) <  10.0D0*degtorad) then
 !             AKT(I,J,K,1,iblock)= AKT (I,J,K,1,iblock) +8.0D-3
 !             AKT(I,J,K,2,iblock)= AKT (I,J,K,2,iblock) +8.0D-3
 !          end if
 !       END DO
 !    END DO
 !    END DO
 ! END DO

!---------------------------------------------------------------------
!     PREPARATION FOR VERTICAL ADVECTIVE TERM
!---------------------------------------------------------------------

!!$OMP PARALLEL DO PRIVATE (IBLOCK,J,I) collapse(2) 
   DO iblock = 1, nblocks_clinic
      DO J = 1, JMT
         DO I = 1,IMT
            H0F (I,J,iblock)= H0F (I,J,iblock)* ONBC
         END DO
      END DO
   END DO


!!$OMP PARALLEL DO PRIVATE (IBLOCK,J,I) collapse(3)
   DO iblock = 1, nblocks_clinic
      DO J = 1, JMT
         DO I = 1,IMT
            STF (I,J,IBLOCK)= AA * H0F (I,J,IBLOCK) + (1.0D0- AA)* H0L (I,J,IBLOCK)
         END DO
      END DO
   END DO

!!$OMP PARALLEL DO PRIVATE (IBLOCK,K,J,I) collapse(4)
   DO iblock = 1, nblocks_clinic
      DO K = 1,KM
         DO J = 1, JMT
            DO I = 1,IMT
               UTF (I,J,K,IBLOCK)= UTF (I,J,K,IBLOCK)* ONCC
               VTF (I,J,K,IBLOCK)= VTF (I,J,K,IBLOCK)* ONCC
            END DO
         END DO
      END DO
  END DO


!!$OMP PARALLEL DO PRIVATE (IBLOCK,K,J,I) collapse(4)
   DO iblock = 1, nblocks_clinic
      DO K = 1,KM
         DO J = 1, JMT
            DO I = 1,IMT
               WKD (I,J,K,IBLOCK)= AA * UTF (I,J,K,IBLOCK) + (1.0D0- AA)* UTL (I,J,K,IBLOCK)
               WKB (I,J,K,IBLOCK)= AA * VTF (I,J,K,IBLOCK) + (1.0D0- AA)* VTL (I,J,K,IBLOCK)
            END DO
         END DO
      END DO
  END DO
#endif
      !ttt1=clock0f()
      !ttt1=omp_get_wtime()
#ifdef SW_ATHREAD
    !work = c0
 !     call upwell_sw
#ifdef ENABLE_GPLT
    call t_startf("tracer_upwell"//CHAR(0))
#endif
    call athread_spawn(slave_tracer_upwell_cu_0,tracer_upwell_0_param)
    call athread_join()

    call athread_spawn(slave_tracer_upwell_cu,tracer_upwell_param)
    call athread_join()

    call athread_spawn(slave_tracer_upwell_cu_1,tracer_upwell_1_param)
    call athread_join()

    call athread_spawn(slave_tracer_upwell_cu_2,tracer_upwell_2_param)
    call athread_join()
#ifdef ENABLE_GPLT
    call t_stopf("tracer_upwell"//CHAR(0))
#endif
#else 
      CALL UPWELL (WKD,WKB,STF)

#endif


#if (defined NODIAG)
 
!---------------------------------------------------------------------
!     INITIALIZE WORK ARRAYS
!---------------------------------------------------------------------
 
 
!-----------------------------------------------------------------------
!     PREPARATION FOR ISOPYCNAL DIFFUSION & ADVECTION
!-----------------------------------------------------------------------
#if (defined ISO)
      CALL ISOPYC
#endif
 
!@@@  COMPUTING DIFFUSION COEFFICIENT

#ifdef SW_ATHREAD
#ifdef ENABLE_GPLT
    call t_startf("tracer_wkc"//CHAR(0))
#endif
    call athread_spawn(slave_tracer_wkc_kernel,tracer_wkc_parm)
    call athread_join()
#ifdef ENABLE_GPLT
    call t_stopf("tracer_wkc"//CHAR(0))
#endif
#else
!!$OMP PARALLEL DO PRIVATE (IBLOCK,K,J,I) collapse(4)
   DO iblock = 1, nblocks_clinic
      DO K = 1,KM
         DO J = 2, JMT-1
            DO I = 2, IMT-1
#if (defined ISO)
               WKC (I,J,K,IBLOCK) = AHV + AHISOP(i,j,iblock) * K3 (I,K,J,3,IBLOCK)
#else
               WKC (I,J,K,IBLOCK) = AHV
#endif
!
!            AKT(I,J,K,1,IBLOCK)=WKC(I,J,K,IBLOCK)
!            AKT(I,J,K,2,IBLOCK)=WKC(I,J,K,IBLOCK)
!
            END DO
         END DO
      END DO
   END DO

#endif
!
!-----------------------------------------------------------------------
!     SOLVE FOR ONE TRACER AT A TIME
!-----------------------------------------------------------------------
!     NTRA = 1 => TEMPERATURE
!     NTRA = 2 => SALINITY
 
      DO N = 1,NTRA
!---------------------------------------------------------------------
!     COMPUTE THE ADVECTIVE TERM 
!---------------------------------------------------------------------


   do iblock = 1, nblocks_clinic
      !adv_tt = 0.0_r8
      this_block = get_block(blocks_clinic(iblock),iblock)

 !    ttt3=clock0f()
#ifdef SW_ATHREAD
#ifdef ENABLE_GPLT
    call t_startf("tracer_advection"//CHAR(0))
#endif

        if ( adv_tracer(1:8) == 'centered' .OR. adv_tracer(1:4) == 'flux' ) then
            advection1_param%d_mtracer = N-1
            call athread_spawn(slave_advection_tracer_1_0,advection1_0_param)
            call athread_join()
            call athread_spawn(slave_advection_tracer_1,advection1_param)
            call athread_join()
        else if(adv_tracer(1:5) == 'tspas') then
            advection2_param%d_mtracer = N-1
            advection2_param%d_at      =>at(1,1,1,n,iblock)
            call athread_spawn(slave_advection_tracer_1_0,advection1_0_param)
            call athread_join()
            call  athread_spawn(slave_advection_tracer_2,advection2_param)
            call athread_join()
            call athread_spawn(slave_advection_tracer_3,advection2_param)
            call athread_join()
        else if(adv_tracer(1:3)=='prm')   then
            advection2_param%d_mtracer = N-1
            advection2_param%d_at      =>at(1,1,1,n,iblock)
            call  athread_spawn(slave_advection_tracer_4,advection2_param)
            call athread_join()
            call  athread_spawn(slave_advection_tracer_5,advection2_param)
            call athread_join()
            call  athread_spawn(slave_advection_tracer_6,advection2_param)
            call athread_join()
!            call  athread_spawn(slave_advection_tracer_7,advection2_param)
!            call athread_join()
!            call  athread_spawn(slave_advection_tracer_8,advection2_param)
!            call athread_join()
            call  athread_spawn(slave_advection_tracer_9,advection2_param)
            call athread_join()
!            DO K = 1, KM
!              DO J = 1, JMT
!                DO I = 1, IMT
!                    if(mytid == 10) then
!                        if(adv_tt(I,J,K) .ne. 0.0) then
!                            write(*,*) "adv_tt,ax,ay,az,ql,qr:",i,j,k,adv_tt(I,J,K),adv_x(I,J,K),adv_y(I,J,K),adv_z(I,J,K), flxqLtmp(I,J,K), flxqRtmp(I,J,K)
!                        end if
!                   end if
!                ENDDO
!              ENDDO
!            ENDDO

        endif
#ifdef ENABLE_GPLT
    call t_stopf("tracer_advection"//CHAR(0))
#endif
#else

       tf=0.0D0 !LPF20160829
      call advection_tracer(wkd(:,:,:,iblock),wkb(:,:,:,iblock),ws(:,:,:,iblock), &
                          at(:,:,:,n,iblock),adv_tt,iblock,N,this_block)
#endif
!        ttt4=clock0f()

!    call advection_tracer(wkd(:,:,:,iblock),wkb(:,:,:,iblock),ws(:,:,:,iblock), &
!                            at(:,:,:,n,iblock),adv_tt,iblock,N)



#ifndef SW_ATHREAD
      do k=1, km
      do j =3, jmt-2
      do i =3, imt-2
         tf(i,j,k,iblock) = adv_tt(i,j,k)*vit(i,j,k,iblock) 
!        tf(i,j,k,iblock) = adv_tt(i,j,k)*vit(i,j,k,iblock) +  &
!                          0.1*gamma*(restore(i,j,k,n,iblock)- atb(i,j,k,n,iblock))
!ZWP20170210 value before gamma oraiginal = 0.1
      end do
      end do
      end do
#endif

   end do

#ifdef SW_ATHREAD
#ifdef ENABLE_GPLT
    call t_startf("tracer_del4"//CHAR(0))
#endif
#if (defined BIHAR)
    h_d2tk = c0
    tracer_del4_1_param%n              = n-1

    call athread_spawn(slave_tracer_del4_1,tracer_del4_1_param)
    call athread_join()
#endif
    tracer_del4_2_param%n              = n-1
    call athread_spawn(slave_tracer_del4_2,tracer_del4_2_param)
    call athread_join()
#ifdef ENABLE_GPLT
    call t_stopf("tracer_del4"//CHAR(0))
#endif
#else

#if (defined CANUTO || defined CANUTO2010)
!!$OMP PARALLEL DO PRIVATE (IBLOCK,K,J,I) collapse(4)
   DO IBLOCK = 1, NBLOCKS_CLINIC
      DO K = 1,KM
         DO J = 2, JMT-1
            DO I = 2, IMT-1
                 IF (AKT(I,J,1,N,IBLOCK).lt.dwndmix) AKT(I,J,1,N,IBLOCK)=dwndmix !LPF20160409 
#if (defined ISO)       
!#if ( defined TIDEMIX )
!               WKC (I,J,K,IBLOCK) = AKT(I,J,K,N,IBLOCK)+AK_TIDE(I,J,K,IBLOCK) + AHISOP(i,j,iblock) * K3 (I,K,J,3,IBLOCK) 
!#else
               WKC (I,J,K,IBLOCK) = AKT(I,J,K,N,IBLOCK) + AHISOP(i,j,iblock) * K3 (I,K,J,3,IBLOCK) 
!#endif
#else            
!#if ( defined TIDEMIX )
!               WKC (I,J,K,IBLOCK) = AKT(I,J,K,N,IBLOCK)+AK_TIDE(I,J,K,IBLOCK) 
!#else
               WKC (I,J,K,IBLOCK) = AKT(I,J,K,N,IBLOCK)
!#endif
#endif         
            END DO
         END DO
      END DO
   END DO
!
#endif

!-----------------------------------------------------------------------
!     COMPUTE THE ISOPYCNAL/DIPYCNAL MIXING
!-----------------------------------------------------------------------
!     XZ AND YZ ISOPYCNAL DIFFUSIVE FLUX ARE SOLVED EXPLICITLY;
!     WHILE ZZ COMPONENT WILL BE SOLVED IMPLICITLY.
 
#ifdef LPFDIAG
!LPF20160514
          if (mytid ==0) WRITE (16,*)'before ISOFLUX'
          call outsumglobal(TF,outsum)
          if (mytid ==0) WRITE (16,*)'outsum',outsum 
!LPF20160514
#endif
 
#if (defined ISO)
         CALL ISOFLUX (N)

#ifdef LPFDIAG
!LPF20160514
          if (mytid ==0) WRITE (16,*)'after ISOFLUX'
          call outsumglobal(TF,outsum)
          if (mytid ==0) WRITE (16,*)'outsum',outsum
!LPF20160514
#endif

#else 
 
#if ( defined SMAG)
         CALL SMAG3
#else
!-----------------------------------------------------------------------
!     COMPUTE THE HORIZONTAL DIFFUSION TERMS:  ZONAL and Meridional COMPONENTS
!-----------------------------------------------------------------------
 
#if (defined BIHAR)
!
!!!$OMP PARALLEL DO PRIVATE (IBLOCK,this_block,K,hdtk)
   do iblock = 1, nblocks_clinic
   this_block = get_block(blocks_clinic(iblock),iblock)
!!!$OMP PARALLEL DO PRIVATE (this_block,K,hdtk,dt2k)
   do k =1, km
!20200221 LPF-yyq
      ttt5=clock0f()
      !ttt5=omp_get_wtime()
       call hdifft_del4(k,DT2K,HDTK,ATB(:,:,k,n,iblock),this_block)
      ttt6=clock0f()
      !ttt6=omp_get_wtime()
      bb3 =  bb3 + ttt6 - ttt5
!      if(mytid == 0)write(*,*)"huangjian test upwell time is ",bb3,"s" 
       !call hdifft_del4(k,HDTK,ATB(:,:,k,n,iblock),this_block)
 
!xjw20251110
      if(N==1)then
          HDTK = HDTK*1.0  !! BIHAR tt
      else
          HDTK = HDTK*1.5  !! BIHAR ss
      end if

!20200221 LPF-yyq
      do j= 3, jmt-2
      do i= 3, imt-2
           TF (I,J,K,IBLOCK) = TF (I,J,K,IBLOCK) + HDTK(I,J)
#if (defined LOWRES)
           dx(i,j,k,N,iblock)=HDTK(i,j) !LPF20160822
#endif
      end do
      end do
           !dx(:,:,k,N,iblock)=HDTK(:,:) !LPF20160822
   end do
   end do
!
#else
!
!!$OMP PARALLEL DO PRIVATE (IBLOCK,this_block,K,hdtk)
   do iblock = 1, nblocks_clinic
   this_block = get_block(blocks_clinic(iblock),iblock)
   do k =1, km
      ttt7=clock0f()
      !ttt7=omp_get_wtime()
      call hdifft_del2(k,HDTK,ATB(:,:,k,n,iblock),this_block)
      ttt8=clock0f()
      !ttt8=omp_get_wtime()
      bb4 = bb4+ ttt8 - ttt7
!      if(mytid == 0)write(*,*)"huangjian test upwell time is ",bb4,"s" 

!xjw20250418
      if(N==1)then
          HDTK = HDTK*0.05  !! Laplacian tt
      else
          HDTK = HDTK*0.05  !! Laplacian ss
      end if

      do j= 3, jmt-2
      do i= 3, imt-2
           TF (I,J,K,IBLOCK) = TF (I,J,K,IBLOCK) + HDTK(I,J)
#if (defined LOWRES)
           dx(i,j,k,N,iblock)=HDTK(i,j) !LPF20160822
#endif
      end do
      end do
           !dx(:,:,k,N,iblock)=HDTK(:,:) !LPF20160822
   end do
   end do

!
#endif
!BIHAR end
#endif
!SMAG end
#endif
!ISO end

#endif
!SW_ATHREAD end for 394&411


#ifdef LPFDIAG
!LPF20160514
          if (mytid ==0) WRITE (16,*)'after1 ISOFLUX'
          call outsumglobal(TF,outsum)
          if (mytid ==0) WRITE (16,*)'outsum',outsum 
!LPF20160514
#endif

!-----------------------------------------------------------------------
!     VERTICAL COMPONENT
!-----------------------------------------------------------------------

#ifdef SW_ATHREAD
#ifdef ENABLE_GPLT
    call t_startf("tracer5"//CHAR(0))
#endif
  call tracer5_cu_slave(n)
#ifdef ENABLE_GPLT
    call t_stopf("tracer5"//CHAR(0))
#endif
#else
         IF (N == 1)THEN
#if (defined SOLAR)
!     SOLAR SHORTWAVE PENETRATION
!!$OMP PARALLEL DO PRIVATE (IBLOCK,J,I,wt1,wt2) collapse(3)
    DO IBLOCK = 1, NBLOCKS_CLINIC
        DO K=2,KM-1
        DO J=3, JMT-2
        DO I=3,IMT-2
            wt1= SWV (I,J,IBLOCK)*pen(k-1)*VIT(I,J,K,IBLOCK)
            wt2= SWV (I,J,IBLOCK)*pen(k)*VIT(I,J,K+1,IBLOCK)
            TF (I,J,K,IBLOCK)= TF(I,J,K,IBLOCK)+(wt1-wt2)*ODZP(K)
!
!lhl0105           penetrate(i,j,k)=(wt1-wt2)*ODZP(k)
           penetrate(i,j,k,iblock)=(wt1-wt2)*ODZP(k)
!            penetrate(i,j,k,IBLOCK)= wt2*ODZP(k)
!
        END DO
        END DO
        END DO
     END DO
!!$OMP PARALLEL DO PRIVATE (IBLOCK,J,I,wt1,wt2) collapse(3)
     DO IBLOCK = 1, NBLOCKS_CLINIC
        DO J=3, JMT-2
        DO I=3,IMT-2
           wt1= SWV (I,J,IBLOCK)*pen(1)*VIT(I,J,2,IBLOCK)
           wt2= SWV (I,J,IBLOCK)*pen(km-1)*VIT(I,J,km,IBLOCK)
           TF(I,J,1,IBLOCK)=TF(I,J,1,IBLOCK)-ODZP(1)*wt1
           TF(I,J,km,IBLOCK)=TF(I,J,km,IBLOCK)+ODZP(km)*wt2
!
            penetrate(i,j, 1,IBLOCK)= -wt1*ODZP( 1) !LPF20160824
            penetrate(i,j,km,IBLOCK)= wt2*ODZP(km) !LPF20160828
            !penetrate(i,j,km,IBLOCK)= -wt2*ODZP(km) !LPF20160824
!
        END DO
        END DO
     END DO
#endif


!==================================
!From here ie. the code about SOLARCHLORO added by linpf
!==================================

#if (defined SOLARCHLORO)
!!$OMP PARALLEL DO PRIVATE (IBLOCK,J,I,wt1,wt2) collapse(4)
   DO IBLOCK = 1, NBLOCKS_CLINIC
      DO K=2,KM-1
        DO J=3,JMT-2
        DO I=3,IMT-2
            wt1= SWV(I,J,IBLOCK)*pen_chl(I,J,K-1,IBLOCK)*VIT(I,J,K,IBLOCK)
            wt2= SWV(I,J,IBLOCK)*pen_chl(I,J,K,IBLOCK)*VIT(I,J,K+1,IBLOCK)
            TF (I,J,K,IBLOCK)= TF(I,J,K,IBLOCK)+(wt1-wt2)*ODZP(K)
!
            penetrate(I,J,K,IBLOCK)=(wt1-wt2)*ODZP(K)
!
        END DO
        END DO
      END DO
  END DO
!!$OMP PARALLEL DO PRIVATE (IBLOCK,J,I,wt1,wt2) collapse(4)
   DO IBLOCK = 1, NBLOCKS_CLINIC
        DO J=3,JMT-2
        DO I=3,IMT-2
           wt1= SWV(I,J,iblock)*pen_chl(I,J,1,iblock)*VIT(I,J,2,iblock)
           wt2= SWV(I,J,iblock)*pen_chl(I,J,km-1,iblock)*VIT(I,J,km,iblock)
           TF(I,J,1,iblock)=TF(I,J,1,iblock)-ODZP(1)*wt1
           TF(I,J,km,iblock)=TF(I,J,km,iblock)+ODZP(km)*wt2
!
            penetrate(I,J, 1,iblock)=-wt1*ODZP(1)
            penetrate(I,J,km,iblock)= wt2*ODZP(km)
!
        END DO
        END DO
  END DO
#endif

!==================================
!above code about SOLARCHLORO added by linpf
!==================================
         END IF
 
!     EDDY-DIFFUSION
 
 
        wt1=0
!!$OMP PARALLEL DO PRIVATE (IBLOCK,J,I,wt1,wt2) collapse(4)
   DO IBLOCK = 1, NBLOCKS_CLINIC
        DO K=2,KM-1
           DO J=3,JMT-2
           DO I=3,IMT-2
              wt1= WKC(I,J,K-1,IBLOCK)*(ATB(I,J,K-1,N,IBLOCK)-ATB(I,J,K,N,IBLOCK))*ODZT(K)*VIT(I,J,K,IBLOCK)
              wt2= WKC(I,J,K,IBLOCK)*(ATB(I,J,K,N,IBLOCK)-ATB(I,J,K+1,N,IBLOCK))*ODZT(K+1)*VIT(I,J,K+1,IBLOCK)
              TF (I,J,K,IBLOCK)= TF(I,J,K,IBLOCK)+ODZP(K)*(wt1-wt2)*(1.0D0-AIDIF)
#if (defined LOWRES)
              dz(i,j,k,N,IBLOCK)=ODZP(K)*(wt1-wt2)*(1.0D0-AIDIF)
#endif
           END DO
           END DO
        END DO
  END DO
!!$OMP PARALLEL DO PRIVATE (IBLOCK,J,I,wt1,wt2) collapse(3)
   DO IBLOCK = 1, NBLOCKS_CLINIC
       DO J = 3,JMT-2
       DO I = 3,IMT-2
           wt1= WKC(I,J,1,IBLOCK)*(ATB(I,J,1,N,IBLOCK)-ATB(I,J,2,N,IBLOCK))*ODZT(2)*VIT(I,J,2,IBLOCK)
           wt2= WKC(I,J,km-1,IBLOCK)*(ATB(I,J,km-1,N,IBLOCK)-ATB(I,J,km,N,IBLOCK))*ODZT(km)*VIT(I,J,km,IBLOCK)
           TF(I,J,1,IBLOCK)=TF(I,J,1,IBLOCK)-ODZP(1)*wt1*(1.0D0-AIDIF)
           TF(I,J,km,IBLOCK)=TF(I,J,km,IBLOCK)+ODZP(km)*wt2*(1.0D0-AIDIF)
#if (defined LOWRES)
              dz(i,j, 1,N,IBLOCK)=-ODZP( 1)*wt1*(1.0-AIDIF)
              dz(i,j,km,N,IBLOCK)= ODZP(km)*wt2*(1.0-AIDIF)
#endif
       END DO
       END DO
  END DO        
  
          !if(mytid==100) write(16,*) 'TF,0-0',TF(21,10,1,1)
 !needto update the halo !LPF20150914 

!-----------------------------------------------------------------------
!     SET NEWTONIAN SURFACE BOUNDARY CONDITION
!-----------------------------------------------------------------------
 
         IF (N == 2)THEN
 
!!$OMP PARALLEL DO PRIVATE (IBLOCK,J,I)    
        DO IBLOCK = 1, NBLOCKS_CLINIC
            DO J = 1, JMT
               DO I = 1, IMT
                  IF (KMT(I,J,IBLOCK) > 0)THEN
#ifdef COUP
!lhl20130914
! under seaice 50m/30days, all 50m/4yr (GAMMA=1/30days)
!                    STF (I,J,IBLOCK) = SSF(I,J,IBLOCK)/ODZP(1)&
!                      +(GAMMA*10)*(SSS(I,J,IBLOCK) &
!                      -ATB (I,J,1,2,IBLOCK))*ifrac(i,j,IBLOCK)/ODZP(1)&
!                      +(GAMMA)*(1-ifrac(i,j,iblock))*(SSS(I,J,IBLOCK)-ATB(I,J,1,2,IBLOCK))/ODZP(1)

                    STF (I,J,IBLOCK) = SSF(I,J,IBLOCK)/ODZP(1)
!lhl20130914
!            STF (I,J,IBLOCK) = SSF(I,J,IBLOCK)/ODZP(1) + GAMMA * (RESTORE (I,J,1,N,iblock) - ATB (I,J,1,2,IBLOCK))/ODZP(1)

#else
#ifdef FRC_CORE
                    !STF (I,J,IBLOCK) = (fresh(i,j,IBLOCK)*34.7*OD0*1.0D-3 & 
                    !+GAMMA*(SSS(I,J,IBLOCK)-ATB (I,J,1,2,IBLOCK))*seaice(i,j,IBLOCK)/ODZP(1)&
                    !+GAMMA*(SSS(I,J,IBLOCK)-ATB (I,J,1,2,IBLOCK))/ODZP(1)*(1.0-seaice(i,j,IBLOCK)))
                    STF (I,J,IBLOCK) = fresh(i,j,IBLOCK)*34.7*OD0*1.0D-3
#else
!                     STF (I,J,IBLOCK) = GAMMA * (SSS (I,J,IBLOCK) - ATB (I,J,1,2,IBLOCK))
                      STF (I,J,IBLOCK) = GAMMA * (SSS (I,J,IBLOCK) - ATB (I,J,1,2,IBLOCK))/ODZP(1)
#endif
#endif
!                     TF (I,J,1,IBLOCK) = TF (I,J,1,IBLOCK) + STF (I,J,IBLOCK)* (1.0- AIDIF)
                     TF (I,J,1,IBLOCK) = TF (I,J,1,IBLOCK) + STF (I,J,IBLOCK)*(1.0 - AIDIF)*ODZP(1)
!
!                     NET (I,J,2,IBLOCK) = STF (I,J,IBLOCK)*(1.0-AIDIF)
                     NET (I,J,2,IBLOCK) = STF (I,J,IBLOCK)*ODZP(1)
!
                  END IF
               END DO
            END DO
      END DO

#ifdef SSS_restore
      sss_restore = 0.0
      DO IBLOCK = 1, NBLOCKS_CLINIC
         DO J = 1, JMT
            DO I = 1, IMT
               IF (KMT(I,J,IBLOCK) > 0)THEN
#ifdef COUP

                  sss_restore(I,J,IBLOCK) = (GAMMA*10)*(SSS(I,J,IBLOCK) &
                      -ATB (I,J,1,2,IBLOCK))*ifrac(i,j,IBLOCK)/ODZP(1)&
                      +(GAMMA)*(1-ifrac(i,j,iblock))*(SSS(I,J,IBLOCK)-ATB(I,J,1,2,IBLOCK))/ODZP(1)
#else
                  sss_restore(I,J,IBLOCK) = (7.0/3.0)*GAMMA*(SSS(I,J,IBLOCK)-ATB (I,J,1,2,IBLOCK))*seaice(i,j,IBLOCK)&
                                          + GAMMA*(SSS(I,J,IBLOCK)-ATB (I,J,1,2,IBLOCK))*(1.0-seaice(i,j,IBLOCK))
#endif
               end if
               TF (I,J,1,IBLOCK) = TF(I,J,1,IBLOCK) + sss_restore(I,J,IBLOCK)
            end do
         end do
      end do
#endif

          !if(mytid==100) write(16,*) 'TF,0-1',TF(21,10,1,1)

!LPF20160504
!===================================
!lhl20130913
end if

#endif
! for SW_ATHREAD in line 564&572


    IF(N==2) THEN
#ifdef SSSNORM
! normalization for fresh water flux
!       allocate(temp11(imt,jmt,max_blocks_clinic))
      if(.false.) allocate(temp12(imt_global,jmt_global))
      ERR_norm2=0.0
!!$OMP PARALLEL DO PRIVATE (IBLOCK,J,I) collapse(3)

#ifdef SW_ATHREAD1
     d_err_norm2(1)     = ERR_norm2
     call athread_spawn(slave_tracer_norm2_kernel,tracer_norm2_parm)
     call athread_join()
     ERR_norm2 = d_err_norm2(1)
#else
    DO IBLOCK = 1, NBLOCKS_CLINIC
     DO J=3,JMT-2
     DO I=3,IMT-2
    !  DO J = JST,JMT
    !     DO I = 1,IMT
          ERR_norm2=ERR_norm2+tarea(i,j,iblock)*NET(I,J,2,iblock)*VIT(I,J,1,iblock)
     END DO
     END DO
    END DO
#endif
!!$OMP END PARALLEL DO

#ifdef SPMD
      call mpi_comm_size(MPI_COMM_OCN, rank_size, ierr)
      ! call mpi_reduce(ERR_norm2,ERR_norm1,1,MPI_PR,mpi_sum,0,mpi_comm_ocn,ierr)
      call mpi_gather(ERR_norm2, 1, MPI_DOUBLE, tmp_sum, 1, MPI_DOUBLE, 0, MPI_COMM_OCN, ierr)
      ERR_norm1 = 0.0
      do i = 1,  rank_size
         ERR_norm1 = ERR_norm1 + tmp_sum(i)
      enddo


      !call mpi_reduce(temp11,ERR_norm1,1,MPI_PR,mpi_sum,0,mpi_comm_ocn,ierr)

      if(.false.) then
      call gather_global(temp12, temp11, master_task, distrb_clinic) 
      if(mytid==0) then
        ERR_norm1=0.D0
!!$OMP PARALLEL DO PRIVATE (J,I) collapse(2)
       do j=1,jmt_global
          do i=1,imt_global
            ERR_norm1=ERR_norm1+temp12(i,j)
          enddo
       enddo
     endif
    endif !noworkforSUPERHIGH

       CALL MPI_BCAST(ERR_norm1,1,MPI_PR,0,MPI_COMM_OCN,IERR)
        ERR_norm2 = - ERR_norm1/area_t

!      if (mytid.eq.0) then
!      write(*,'(a60,2D25.15)') "Gross of fresh water flux (total psu*m/s,
!      averaged psu/s) =",err1,err2
!      endif


#else
!      ERR_norm1=0.D0
!       do j=jstart+1,jmt_global-1
!          do i=2,imt_global-1
!            ERR_norm1=ERR_norm1+temp11(i,j)
!          enddo
!       enddo
!      ERR_norm2 = - ERR_norm1 / ASEA
!      write(*,'(a60,2D25.15)') "Gross of fresh water flux (total psu*m/s,
!      averaged psu/s) =",err,err2
#endif

       FW_norm2=ERR_norm2 !LPF20160823

#ifdef SW_ATHREAD
#ifdef ENABLE_GPLT
    call t_startf("tracer_norm3"//CHAR(0))
#endif

     tracer_norm3_parm%d_err_norm2     = ERR_norm2
     call athread_spawn(slave_tracer_norm3_kernel,tracer_norm3_parm)
     call athread_join()
#ifdef ENABLE_GPLT
    call t_stopf("tracer_norm3"//CHAR(0))
#endif

#else
! for SW_ATHREAD in line 838
!!$OMP PARALLEL DO PRIVATE (IBLOCK,J,I) collapse(3)
     DO IBLOCK = 1, NBLOCKS_CLINIC
      DO J = JST,JMT
         DO I = 1,IMT
            TF (I,J,1,IBLOCK)=TF(I,J,1,IBLOCK)+ERR_norm2*VIT(I,J,1,IBLOCK)
         END DO
       END DO
      END DO
#endif
! for SW_ATHREAD in line 838&851

!!$OMP END PARALLEL DO
         
         ! !if(mytid==100) then
         !   DO IBLOCK = 1, NBLOCKS_CLINIC
         !    DO J = JST,JMT
         !     DO I = 1,IMT         
         !       if(abs(tf(i,j,k,iblock))>0)write(16,*) 'TF,1-0',TF(i,j,1,iblock),i,j,iblock
         !     ENDDO
         !     ENDDO
         !   ENDDO
         !endif
          !if(mytid==100) write(16,*) 'TF,1-0',TF(21,10,1,1)

!lhl20130913
!===================================
!       deallocate(temp11)
      if(.false.) deallocate(temp12) 
#endif
!for SSSNORM in line 766

    ELSE !for temperature, for line 765, the following is for temperature

#ifdef SW_ATHREAD

      call athread_spawn(slave_tracer_norm4_kernel,tracer_norm4_parm)
      call athread_join()

#else
! !for not SW_ATHREAD

!!$OMP PARALLEL DO PRIVATE (IBLOCK,J,I)
        DO IBLOCK = 1, NBLOCKS_CLINIC
            DO J = 1, JMT
               DO I = 1, IMT
                IF (KMT(I,J,IBLOCK) > 0)THEN
#ifdef COUP
                 STF (I,J,IBLOCK) = TSF(I,J,IBLOCK)
#else
#ifdef FRC_CORE
!                 STF (I,J,IBLOCK) = ((SWV(I,J,IBLOCK)+NSWV(I,J,IBLOCK))*OD0CP+ & 
!                                      SEAICE(I,J,IBLOCK)*GAMMA*(SST(I,J,IBLOCK)-ATB(I,J,1,1,IBLOCK))/ODZP(1))

                 STF(I,J,IBLOCK) = (SWV(I,J,IBLOCK) + NSWV(I,J,IBLOCK))*OD0CP       
#else
                 STF (I,J,IBLOCK) = (SWV(I,J,IBLOCK)+NSWV(I,J,IBLOCK)-DQDT(I,J,IBLOCK)* & 
                                    (SST (I,J,IBLOCK) - ATB (I,J,1,1,IBLOCK)))*OD0CP
#endif
#endif
                 TF(I,J,1,IBLOCK) = TF(I,J,1,IBLOCK) + STF(I,J,IBLOCK)* ODZP(1)*(1.0D0 - AIDIF)
!
!                 NET (I,J,1,IBLOCK) = STF (I,J,IBLOCK)*ODZP(1)*(1.0-AIDIF)
                 NET (I,J,1,IBLOCK) = STF (I,J,IBLOCK)*ODZP(1)
!
                END IF !KMT
               END DO !I
            END DO !J
        END DO !IBLOCK

#ifdef SST_restore
        sst_restore = 0.0
        DO IBLOCK = 1, NBLOCKS_CLINIC
           DO J = 1, JMT
              DO I = 1, IMT
                 IF (KMT(I,J,IBLOCK) > 0)THEN
                    sst_restore(I,J,IBLOCK) = (7.0/3.0)*GAMMA*seaice(i,j,IBLOCK)*(SST(I,J,IBLOCK) - ATB(I,J,1,1,IBLOCK)) + & 
                                              GAMMA*(1.0-seaice(i,j,IBLOCK))*(SST(I,J,IBLOCK) - ATB(I,J,1,1,IBLOCK))  
                 end if
                 TF (I,J,1,IBLOCK) = TF (I,J,1,IBLOCK) + sst_restore(I,J,IBLOCK)
              end do
           end do
        end do
#endif

#endif
! for SW_ATHREAD for line 886&898

    END IF !N==2or1 for line 765&884

          !if(mytid==100) write(16,*) 'TF,1-1',TF(21,10,1,1)

#ifdef SW_ATHREAD
#ifdef ENABLE_GPLT
    call t_startf("tracer7"//CHAR(0))
#endif
  call tracer7_cu_slave(n)
#ifdef ENABLE_GPLT
    call t_stopf("tracer7"//CHAR(0))
#endif

#else

          !LPF20170729
          !for test the boudary net
          call POP_HaloUpdate3DR8_opt(net , POP_haloClinic, POP_gridHorzLocCenter,&
                  POP_fieldKindScalar, km, errorCode, fillValue = 0.0_r8)
          !for test the boudary net
     if ( simple_assm ) then
     !if ( simple_assm==.true.) then
          if(mytid==0) write(16,*) 'into restoring,n=',N
!-----------------------------------------------------------------------
!   restoring to obs temp and salt
!-----------------------------------------------------------------------
!!$OMP PARALLEL DO PRIVATE (IBLOCK,K,J,I) collapse(4)
     DO IBLOCK = 1, NBLOCKS_CLINIC
       DO K = 1,kvt
         DO j=3,jmt-2
            DO I = 3,imt-2
               TF (I,J,K,iblock)= TF (I,J,K,iblock) + VIT (I,J,K,iblock)* (restore_at (I,J,K,&
                          N,iblock) - ATB (I,J,K,N,iblock))*LAMDA1(k)
                !dt_restore(i,j,k,N,iblock)=VIT (I,J,K,iblock)* (restore_at (I,J,K,&
                 !         N,iblock) - ATB (I,J,K,N,iblock))*LAMDA1(k) 
               END DO
            END DO
         END DO
      END DO
   end if !sim_assm=true 
!LPF20170729 
!-----------------------------------------------------------------------
!IAU
!LPF20200805
     if ( enoi_iau ) then
     do k=1,45 !kvt
      LAMDA1(k)=1.0D0/(iau_gamma*86400.D0)
     enddo
     do k=46,km
      LAMDA1(k)=1.0D0/(iau_gamma*86400.D0) !/(k*k/5.0/5.0)
     enddo
!          if(mytid==0) write(*,*) 'into IAU,AT,n=',N
!-----------------------------------------------------------------------
!   IAU to assim temp and salt
!-----------------------------------------------------------------------
!$OMP PARALLEL DO PRIVATE (IBLOCK,K,J,I)
     DO IBLOCK = 1, NBLOCKS_CLINIC
       DO K = 1,km !kvt
         DO j=3,jmt-2
            DO I = 3,imt-2
               TF (I,J,K,iblock)= TF (I,J,K,iblock) + VIT (I,J,K,iblock)* at_iau (I,J,K,&
                          N,iblock)*LAMDA1(k) !*DTS
!                dt_restore(i,j,k,N,iblock)=VIT (I,J,K,iblock)* at_iau (I,J,K,&
!                          N,iblock)*LAMDA1(k) !*DTS
               END DO
            END DO
         END DO
      END DO
   end if !iau=true
!LPF20200805
!-----------------------------------------------------------------------

if ( boundary_restore == 1) then
!-----------------------------------------------------------------------
!    boundary condition
!-----------------------------------------------------------------------
!!$OMP PARALLEL DO PRIVATE (IBLOCK,K,J,I) collapse(4)
     DO IBLOCK = 1, NBLOCKS_CLINIC
         DO K = 2,KM
         do j=3,jmt-2
               DO I = 3,imt-2
               !DO I = 3,jmt-2
                  TF (I,J,K,iblock)= TF (I,J,K,iblock) + VIT (I,J,K,iblock)* (RESTORE (I,J,K,&
                             N,iblock) - ATB (I,J,K,N,iblock))*LAMDA(i,j,k,iblock)
               END DO
            END DO
         END DO
      END DO
 
!!$OMP PARALLEL DO PRIVATE (IBLOCK,K,J,I) collapse(4)
     DO IBLOCK = 1, NBLOCKS_CLINIC
         DO K = 2,KM
         do j=3,jmt-2
               DO I = 3, imt-2
                  TF (I,J,K,iblock)= TF (I,J,K,iblock) + VIT (I,J,K,iblock)* (RESTORE (I,J,K,&
                             N,iblock) - ATB (I,J,K,N,iblock))*LAMDA(i,j,k,iblock)
               END DO
            END DO
         END DO
    END DO
   end if !boundary  

!-----------------------------------------------------------------------
!     SOLVE FOR "TAU+1" TRACER AT CENTER OF "T" CELLS
!-----------------------------------------------------------------------
      
          !if(mytid==100) write(16,*) 'TF,2-0',TF(21,10,1,1)
    if (adv_tracer(1:5) == 'tspas' .or. adv_tracer(1:3)=='prm' ) then
!!$OMP PARALLEL DO PRIVATE (IBLOCK,K,J,I) collapse(4)
      DO IBLOCK = 1, NBLOCKS_CLINIC
         DO K = 1,KM
            DO J = 3, JMT-2
               DO I = 3, IMT-2
                     VTL (I,J,K,iblock) = AT (I,J,K,N,iblock) + DTS * TF (I,J,K,iblock)
               END DO
            END DO
         END DO
     END DO
   else if (adv_tracer(1:8) == 'centered') then
!!$OMP PARALLEL DO PRIVATE (IBLOCK,K,J,I) collapse(4)
      DO IBLOCK = 1, NBLOCKS_CLINIC
         DO K = 1,KM
            DO J = 3, JMT-2
               DO I = 3, IMT-2
                     VTL (I,J,K,iblock) = ATB (I,J,K,N,iblock) + C2DTTS * TF (I,J,K,iblock)
               END DO
            END DO
         END DO
     END DO
   else
        WRITE (*,*)'adv_tracer2:',trim(adv_tracer)
     call exit_licom(sigAbort,'The false advection option for tracer')
   end if
 
#endif 
!tracer7
!-----------------------------------------------------------------------
!     ADD DT/DT COMPONENT DUE TO IMPLICIT VERTICAL DIFFUSION
!-----------------------------------------------------------------------

#ifdef SW_ATHREAD

#ifdef ENABLE_GPLT
    call t_startf("tracer_invtrit"//CHAR(0))
#endif
!    call athread_spawn(slave_tracer_vtl_kernel,tracer_vtl_parm)
!    call athread_join()

     VTL_ori=VTL !for output dt diffusion !LPF20160823

    invtrit_param%d_dts          = C2DTTS
    call athread_spawn(slave_invtrit,invtrit_param)
    call athread_join()

     !call invtrit_sw(c2dtts) !hx invtrit
#ifdef ENABLE_GPLT
    call t_stopf("tracer_invtrit"//CHAR(0))
#endif
#else
     VTL_ori=VTL !for output dt diffusion !LPF20160823
     CALL INVTRIT(VTL,STF,WKC,AIDIF,C2DTTS)

#endif
#ifdef ENABLE_GPLT
    call t_startf("pop_vtl"//CHAR(0))
#endif

     call POP_HaloUpdate3DR8_opt(VTL , POP_haloClinic, POP_gridHorzLocCenter,&
                         POP_fieldKindScalar, km, errorCode, fillValue = 0.0_r8)
#ifdef ENABLE_GPLT
    call t_stopf("pop_vtl"//CHAR(0))
#endif

#ifdef SW_ATHREAD
#ifdef ENABLE_GPLT
    call t_startf("tracer9"//CHAR(0))
#endif
    tracer_9_param%d_dts          = C2DTTS
    tracer_9_param%n              = n-1
    call athread_spawn(slave_tracer_9,tracer_9_param)
    call athread_join()
#ifdef ENABLE_GPLT
    call t_stopf("tracer9"//CHAR(0))
#endif
#else
!     ts1 = clock0f()
          
!     DO K=1,KM
!      if (K==1) then
!       dt_diff(:,:,K,N,:)=(VTL(:,:,K,:)-VTL_ori(:,:,K,:))/C2DTTS*vit(:,:,K,:)-STF (:,:,:)*AIDIF*ODZP(1) !for output dt diffusion !LPF20160823
!      else
!       dt_diff(:,:,K,N,:)=(VTL(:,:,K,:)-VTL_ori(:,:,K,:))/C2DTTS*vit(:,:,K,:) 
!      endif
!     ENDDO
     !dt_diff(:,:,:,N,:)=(VTL(:,:,:,:)-VTL_ori(:,:,:,:))/C2DTTS*vit(:,:,:,:) !for output dt diffusion !LPF20160823
!
!---------------------------------------------------------------------
!     SET CYCLIC CONDITIONS ON EASTERN AND WESTERN BOUNDARY
!---------------------------------------------------------------------
 
!          !if(mytid==100) write(16,*) 'TF,3-0',TF(21,10,1,1)
               
    if (horiz_grid_opt(1:7) == 'lat_lon') then          
   
     if (mod(ist,180) == 1) then
       CALL SMTS (VTL,VIT,fil_lat2)
     else
!
      if (adv_tracer(1:5) == 'tspas' .or. adv_tracer(1:3)=='prm' ) then
!!$OMP PARALLEL DO PRIVATE (IBLOCK,J,I) collapse(3)
          DO IBLOCK = 1, NBLOCKS_CLINIC
          DO J=1, JMT
          DO I=1,IMT  
             VTL (I,J,1,IBLOCK) = VTL(I,J,1,IBLOCK)-AT (I,J,1,N,IBLOCK) - NET(I,J,N,IBLOCK)*DTS
          END DO
          END DO
          END DO
!!$OMP PARALLEL DO PRIVATE (IBLOCK,K,J,I) collapse(4)
          DO  IBLOCK = 1, NBLOCKS_CLINIC
          DO K=2,KM
          DO J=1, JMT
          DO I=1,IMT  
             VTL (I,J,K,IBLOCK) = VTL(I,J,K,IBLOCK)-AT (I,J,K,N,IBLOCK)
          END DO
          END DO
          END DO
          END DO
      else if (adv_tracer(1:8) == 'centered') then
!!$OMP PARALLEL DO PRIVATE (IBLOCK,J,I) collapse(3)
          DO IBLOCK = 1, NBLOCKS_CLINIC
          DO J=1, JMT
          DO I=1,IMT  
             VTL (I,J,1,IBLOCK) = VTL(I,J,1,IBLOCK)- ATB (I,J,1,N,IBLOCK) - NET(I,J,N,IBLOCK)*C2DTTS
          END DO
          END DO
          END DO
!!$OMP PARALLEL DO PRIVATE (IBLOCK,K,J,I) collapse(4)
          DO  IBLOCK = 1, NBLOCKS_CLINIC
          DO K=2,KM
          DO J=1, JMT
          DO I=1,IMT  
             VTL (I,J,K,IBLOCK) = VTL(I,J,K,IBLOCK)- ATB (I,J,K,N,IBLOCK)
          END DO
          END DO
          END DO
          END DO
!
      else
          WRITE (*,*)'adv_tracer3:',trim(adv_tracer)
         call exit_licom(sigAbort,'The false advection option for tracer')
      end if !choose the advection scheme
!
      CALL SMTS (VTL,VIT,fil_lat1)
!
      if (adv_tracer(1:5) == 'tspas' .or. adv_tracer(1:3)=='prm' ) then
!!$OMP PARALLEL DO PRIVATE (IBLOCK,J,I) collapse(3)
          DO IBLOCK = 1, NBLOCKS_CLINIC
          DO J=1, JMT
          DO I=1,IMT  
             VTL (I,J,1,IBLOCK) = AT (I,J,1,N,IBLOCK) + VTL(I,J,1,IBLOCK) + NET(I,J,N,IBLOCK)*DTS
          END DO
          END DO
          END DO
!!$OMP PARALLEL DO PRIVATE (IBLOCK,K,J,I) collapse(4)
          DO  IBLOCK = 1, NBLOCKS_CLINIC
          DO K=2,KM
          DO J=1, JMT
          DO I=1,IMT  
             VTL (I,J,K,IBLOCK) = AT (I,J,K,N,IBLOCK) + VTL(I,J,K,IBLOCK)
          END DO
          END DO
          END DO
          END DO
      else if (adv_tracer(1:8) == 'centered') then
!!$OMP PARALLEL DO PRIVATE (IBLOCK,J,I) collapse(3)
          DO IBLOCK = 1, NBLOCKS_CLINIC
          DO J=1, JMT
          DO I=1,IMT  
             VTL (I,J,1,IBLOCK) = ATB (I,J,1,N,IBLOCK) + VTL(I,J,1,IBLOCK) + NET(I,J,N,IBLOCK)*C2DTTS
          END DO
          END DO
          END DO
!!$OMP PARALLEL DO PRIVATE (IBLOCK,K,J,I) collapse(4)
          DO  IBLOCK = 1, NBLOCKS_CLINIC
          DO K=2,KM
          DO J=1, JMT
          DO I=1,IMT  
             VTL (I,J,K,IBLOCK) = ATB (I,J,K,N,IBLOCK) + VTL(I,J,K,IBLOCK)
          END DO
          END DO
          END DO
          END DO
      else
           WRITE (*,*)'adv_tracer4:',trim(adv_tracer)
         call exit_licom(sigAbort,'The false advection option for tracer')
      end if !advection scheme
   end if  !ist different
   end if !for lon-lat only 

          !if(mytid==100) write(16,*) 'TF,3-1',TF(21,10,1,1),N
!
!-----------------------------------------------------------------------
!     SOLVE FOR "TAU+1" TRACER AT CENTER OF "T" CELLS
!-----------------------------------------------------------------------
!!$OMP PARALLEL DO PRIVATE (K,J,I) collapse(4)
!      DO IBLOCK = 1, NBLOCKS_CLINIC
!         DO K = 1,KM
!            DO J = 1, JMT
!               DO I = 1,IMT
                 !tend(i,j,k,N,iblock)=(VTL(I,J,K,iblock)-AT(I,J,K,N,iblock))/C2DTTS*VIT(I,J,K,iblock)
                 !AT or ATB
!               END DO
!            END DO
!         END DO
!     END DO
          !if(mytid==100) write(16,*) 'TF,3-2,N=',TF(21,10,1,1),N

  if (adv_tracer(1:5) == 'tspas' .or. adv_tracer(1:3)=='prm' ) then
!!$OMP PARALLEL DO PRIVATE (IBLOCK,K,J,I) collapse(4)
      DO IBLOCK = 1, NBLOCKS_CLINIC
         DO K = 1,KM
            DO J = 1, JMT
               DO I = 1,IMT
                  AT (I,J,K,N,IBLOCK) = VTL (I,J,K,IBLOCK)
                  ATB(I,J,K,N,IBLOCK)=AT (I,J,K,N,IBLOCK)
               END DO
            END DO
         END DO
      END DO
  else if (adv_tracer(1:8) == 'centered') then
         IF (IST >= 1)THEN
!!$OMP PARALLEL DO PRIVATE (IBLOCK,K,J,I) collapse(4)
      DO IBLOCK = 1, NBLOCKS_CLINIC
            DO K = 1,KM
               DO J = 1, JMT
                  DO I = 1,IMT
                     ATB (I,J,K,N,IBLOCK) = AFT2* AT (I,J,K,N,IBLOCK) + AFT1* (ATB (I,&
                                    J,K,N,IBLOCK) + VTL (I,J, K,IBLOCK))
                  END DO
               END DO
            END DO
     END DO
         END IF
 
!!$OMP PARALLEL DO PRIVATE (IBLOCK,K,J,I) collapse(4)
      DO IBLOCK = 1, NBLOCKS_CLINIC
         DO K = 1,KM
            DO J = 1, JMT
               DO I = 1,IMT
                  AT (I,J,K,N,IBLOCK) = VTL (I,J,K,IBLOCK)
               END DO
            END DO
         END DO
      END DO
  else
      WRITE (*,*)'adv_tracer5:',trim(adv_tracer)
      call exit_licom(sigAbort,'The false advection option for tracer')
  end if

#endif
  END DO
!================================================================
       
!      CALL CONVADJ !LPF20160829 !move from mainloop to here
!      DO N = 1, NTRA
!      DO K = 1,KM
!       dt_conv(:,:,k,N,:)=(AT(:,:,k,N,:)-ATB(:,:,k,N,:))/C2DTTS*vit(:,:,k,:) !for output dt diffusion !LPF20160823
!      ENDDO
!      ENDDO 
!      if (trim(adv_tracer) == 'tspas') then
!         ATB(:,:,1:km,:,:)=AT(:,:,1:km,:,:)    
!       endif
!       tend=tend+dt_conv !total tendency

!================================================================
#else
         atb(:,:,:,1,:)= 12.0D0
         atb(:,:,:,2,:)= c0
         at (:,:,:,:,:)= atb(:,:,1:km,:,:)
#endif
      IST = IST +1
!for test the boudary net
!       call POP_HaloUpdate(net , POP_haloClinic, POP_gridHorzLocCenter,&
!                       POP_fieldKindScalar, errorCode, fillValue = 0.0_r8)
!for test the boudary net
!
#ifdef ISO
      deallocate(K1,K2,K3,adv_vetiso)
#endif
!      deallocate(rict_replace)
!      deallocate(stf,tf)
!      deallocate(wkb,wkc,wkd)
!      deallocate(rit)

      RETURN
      END SUBROUTINE TRACER
 
