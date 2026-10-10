module gfs_forcing_mod

    use param_mod
    integer,save :: gfs_init_flag = 0

    public :: gfs_forcing_init,&
              GFS_INTERPOLATE,&
              read_hour_gfs_forcing
contains

SUBROUTINE gfs_forcing_init
    use forc_mod
    IMPLICIT NONE
    integer :: irec

    irec = 1
    if (gfs_init_flag.eq.0) then
        gfs_init_flag=1

        allocate(t10(g_imt,g_jmt,4))
        allocate(u10(g_imt,g_jmt,4))
        allocate(v10(g_imt,g_jmt,4))
        allocate(slp(g_imt,g_jmt,4))
        allocate(q10(g_imt,g_jmt,4))
        allocate(swhf(g_imt,g_jmt,4))
        allocate(lwhf(g_imt,g_jmt,4))
        allocate(precr(g_imt,g_jmt,4))
        allocate(precs(g_imt,g_jmt,4))
        allocate(rf(g_imt,g_jmt,4))
        allocate(si(g_imt,g_jmt,4))

        if (mytid.eq.0) write(*,*) "gfs_init_flag, t2m_00.nc"
        call read_gfs(irec,"t2m_00.nc",t10,"t2m")
        if (mytid.eq.0) write(*,*) "gfs_init_flag, u_10_00.nc"
        call read_gfs(irec,"u_10_00.nc",u10,"u_10")
        if (mytid.eq.0) write(*,*) "gfs_init_flag, v_10_00.nc"
        call read_gfs(irec,"v_10_00.nc",v10,"v_10")
        if (mytid.eq.0) write(*,*) "gfs_init_flag, q2m_00.nc"
        call read_gfs(irec,"q2m_00.nc",q10,"q2m")
        if (mytid.eq.0) write(*,*) "gfs_init_flag, slp_00.nc"
        call read_gfs(irec,"slp_00.nc",slp,"slp")
        if (mytid.eq.0) write(*,*) "gfs_init_flag, swdn_00.nc"
        call read_gfs(irec,"swdn_00.nc",swhf,"swdn")
        if (mytid.eq.0) write(*,*) "gfs_init_flag, lwdn_00.nc"
        call read_gfs(irec,"lwdn_00.nc",lwhf,"lwdn")
        if (mytid.eq.0) write(*,*) "gfs_init_flag, prec_00.nc"
        call read_gfs(irec,"prec_00.nc",precr,"prc")
        if (mytid.eq.0) write(*,*) "gfs_init_flag, pres_00.nc"
        call read_gfs_fillvalue(irec,"snow_00.nc",precs,"snow")
        if (mytid.eq.0) write(*,*) "gfs_init_flag, runoff_00.nc"
        call read_gfs(irec,"runoff_00.nc",rf,"runoff")
        if (mytid.eq.0) write(*,*) "gfs_init_flag, ice_00.nc"
        call read_gfs(irec,"ice_00.nc",si,"ice")

    endif
END SUBROUTINE gfs_forcing_init



!     =================
      SUBROUTINE GFS_INTERPOLATE
!     =================

#include <def-undef.h>
 use precision_mod
 use pconst_mod
 use forc_mod
 use constant_mod
 use grid
 use domain
 use gather_scatter
 use POP_GridHorzMod
 use POP_HaloMod
! use dyn_mod, only: u,v,buffer,buffer_s
 use dyn_mod, only: u,v,buffer
 use tracer_mod, only: at
 use msg_mod
 use blocks
 use var_check
      IMPLICIT NONE
#include <netcdf.inc>

      integer :: mon_day,irec,irech,tnum
      integer ::  ErrorCode


      real(r8),dimension(imt,jmt) :: model_sst,es,qs,zz,uu,vv,windx,windy,theta
      real(r8),dimension(imt,jmt) :: core_sensible,core_latent,core_tau
!
      real(r8), parameter :: tok=273.15
      real(r8), parameter :: epsln=1e-25
      real(r8),dimension(imt,jmt) :: tmp1,tmp2
      integer,save:: gfs_init=0
      integer*8 :: data_size, offset, k_n
      integer*4 :: err_status, w_size, rank, &
            local_block_num, block_start, block_end, err, jra_num

       character(len=40) :: a
!      if ( TNUM.gt.1 ) goto 111
      core_sensible = 0.0     
      core_latent = 0.0      
      core_tau = 0.0      
      model_sst = 0.0      
      es = 0.0      
      qs = 0.0      
      zz = 0.0      
      uu = 0.0      
      vv = 0.0      
      windx = 0.0      
      windy = 0.0      
      theta = 0.0     
      tmp1= 0.0      
      tmp2= 0.0      

! decide recode number
!      IYFM,MON0,IDAY
      mon_day=0
      do i=1,mon0-1
      mon_day=mon_day+nmonth(i)
      enddo


      irec=iday
      irech = F6hREC! 6 hour interval
      if (mytid.eq.0) then
            write(*,*) "irec=",irec," irech=",irech
      endif

      if (gfs_init.eq.0) then
          gfs_init=1

            tsa3(:,:,1,1) = 0.0_r8
            wspdu3(:,:,1,1) = 0.0_r8
            wspdv3(:,:,1,1) = 0.0_r8
            psa3(:,:,1,1) = 0.0_r8
            qar3(:,:,1,1) = 0.0_r8
            swv3(:,:,1,1) = 0.0_r8
            lwv3(:,:,1,1) = 0.0_r8
            rain3(:,:,1,1) = 0.0_r8
            snow3(:,:,1,1) = 0.0_r8
            runoff3(:,:,1,1) = 0.0_r8
            seaice3(:,:,1,1) = 0.0_r8
            tmp_u(:,:,1,1) = 0.0_r8
            tmp_v(:,:,1,1) = 0.0_r8

      endif 

      
      call map_interplation(t10(:,:, irech), tsa3(:,:,1,1), block_index(:,:,3:4), bilin_src_coord, bilin_weight, max_n_bilin_weight)
      call POP_HaloUpdate(tsa3(:,:,1,:) , POP_haloClinic, POP_gridHorzLocCenter,&
                       POP_fieldKindScalar, errorCode, fillValue = 0.0_r8)

      call map_interplation(u10(:,:, irech), tmp_u(:,:,1,1), block_index(:,:,3:4), bilin_src_coord, bilin_weight, max_n_bilin_weight)
      call map_interplation(v10(:,:, irech), tmp_v(:,:,1,1), block_index(:,:,3:4), bilin_src_coord, bilin_weight, max_n_bilin_weight)


      wspdu3(:,:,1,1) = tmp_u(:,:,1, 1)*cos(local_rho(:,:)) + tmp_v(:,:,1,1)*sin(local_rho(:,:))
      wspdv3(:,:,1,1) = -tmp_u(:,:,1, 1)*sin(local_rho(:,:)) + tmp_v(:,:,1,1)*cos(local_rho(:,:))



      call POP_HaloUpdate(wspdu3(:,:,1,:) , POP_haloClinic, POP_gridHorzLocSWcorner,&
              POP_fieldKindVector, errorCode, fillValue = 0.0_r8)

      call POP_HaloUpdate(wspdv3(:,:,1,:) , POP_haloClinic, POP_gridHorzLocSWcorner,&
              POP_fieldKindVector, errorCode, fillValue = 0.0_r8)

      call map_interplation(slp(:,:, irech), psa3(:,:,1,1), block_index(:,:,3:4), bilin_src_coord, bilin_weight, max_n_bilin_weight)
      call POP_HaloUpdate(psa3(:,:,1,:) , POP_haloClinic, POP_gridHorzLocCenter,&
                       POP_fieldKindScalar, errorCode, fillValue = 0.0_r8)

      call map_interplation(q10(:,:, irech), qar3(:,:,1,1), block_index(:,:,3:4), bilin_src_coord, bilin_weight, max_n_bilin_weight)
      call POP_HaloUpdate(qar3(:,:,1,:) , POP_haloClinic, POP_gridHorzLocCenter,&
                       POP_fieldKindScalar, errorCode, fillValue = 0.0_r8)

      call map_interplation(swhf(:,:, irech), swv3(:,:,1,1), block_index(:,:,3:4), bilin_src_coord, bilin_weight, max_n_bilin_weight)
      call POP_HaloUpdate(swv3(:,:,1,:) , POP_haloClinic, POP_gridHorzLocCenter,&
                       POP_fieldKindScalar, errorCode, fillValue = 0.0_r8)

      call map_interplation(lwhf(:,:, irech), lwv3(:,:,1,1), block_index(:,:,3:4), bilin_src_coord, bilin_weight, max_n_bilin_weight)
      call POP_HaloUpdate(lwv3(:,:,1,:) , POP_haloClinic, POP_gridHorzLocCenter,&
                       POP_fieldKindScalar, errorCode, fillValue = 0.0_r8)

      call map_interplation(precr(:,:, irech), rain3(:,:,1,1), block_index(:,:,3:4), bilin_src_coord, bilin_weight, max_n_bilin_weight)
      call POP_HaloUpdate(rain3(:,:,1,:) , POP_haloClinic, POP_gridHorzLocCenter,&
                       POP_fieldKindScalar, errorCode, fillValue = 0.0_r8)

      call map_interplation(precs(:,:, irech), snow3(:,:,1,1), block_index(:,:,3:4), bilin_src_coord, bilin_weight, max_n_bilin_weight)
      call POP_HaloUpdate(snow3(:,:,1,:) , POP_haloClinic, POP_gridHorzLocCenter,&
                       POP_fieldKindScalar, errorCode, fillValue = 0.0_r8)

      call map_interplation(rf(:,:, 1), runoff3(:,:,1,1), block_index(:,:,3:4), bilin_src_coord, bilin_weight, max_n_bilin_weight)
      call POP_HaloUpdate(runoff3(:,:,1,:) , POP_haloClinic, POP_gridHorzLocCenter,&
                       POP_fieldKindScalar, errorCode, fillValue = 0.0_r8)

      call map_interplation(si(:,:, 1), seaice3(:,:,1,1), block_index(:,:,3:4), bilin_src_coord, bilin_weight, max_n_bilin_weight)
      call POP_HaloUpdate(seaice3(:,:,1,:) , POP_haloClinic, POP_gridHorzLocCenter,&
                       POP_fieldKindScalar, errorCode, fillValue = 0.0_r8)  


!111    continue

         do j = jsm,jem
            do i = 2,imm
              uu(i,j)=vit(i,j,1,1)*(u(i,j,1,1)+u(i+1,j,1,1)+u(i,j-1,1,1)+u(i+1,j-1,1,1)) &
              /(viv(i,j,1,1)+viv(i+1,j,1,1)+viv(i,j-1,1,1)+viv(i+1,j-1,1,1)+epsln)
              vv(i,j)=vit(i,j,1,1)*(v(i,j,1,1)+v(i+1,j,1,1)+v(i,j-1,1,1)+v(i+1,j-1,1,1)) &
              /(viv(i,j,1,1)+viv(i+1,j,1,1)+viv(i,j-1,1,1)+viv(i+1,j-1,1,1)+epsln)
            end do
         end do

! transfer core data to what the subroutine need
      do j=1,jmt
         do i=1,imt
! relative speed to surface currents
         windx(i,j)=1.0*(wspdu3(i,j,1,1)-uu(i,j))*vit(i,j,1,1)
         windy(i,j)=1.0*(wspdv3(i,j,1,1)+vv(i,j))*vit(i,j,1,1)
!         windy(i,j)=(wspdv3(i,j,1,1)-vv(i,j))*vit(i,j,1,1)
!  1.0 is from mom4
         wspd3(i,j,1,1)=sqrt(windx(i,j)**2+windy(i,j)**2+1.0)*vit(i,j,1,1)
! using a transient temperature, not daily mean
         model_sst(i,j)=(at(i,j,1,1,1)+tok)*vit(i,j,1,1)
         zz(i,j)=10.
         qs(i,j)=0.98*640380*exp(-5107.4/model_sst(i,j))/1.22*vit(i,j,1,1)
! temperature to potential temperature
         theta(i,j)=tsa3(i,j,1,1)*(100000.0/psa3(i,j,1,1))**0.286*vit(i,j,1,1)
         runoff(i,j,1)=runoff3(i,j,1,1)*vit(i,j,1,1)
         if(seaice3(i,j,1,1) > 100.0) seaice3(i,j,1,1) = 100.0
         seaice(i,j,1)=seaice3(i,j,1,1)*vit(i,j,1,1)/100.0
         rain(i,j,1)=rain3(i,j,1,1)*vit(i,j,1,1)
         if(rain(i,j,1) < 0.0) rain(i,j,1) = 0.0
         snow(i,j,1)=snow3(i,j,1,1)*vit(i,j,1,1)
         psa(i,j,1) = psa3(i,j,1,1) * vit(i,j,1,1)
         end do
      end do


! compute heat flux
       call ncar_ocean_fluxes_jra(wspd3(1,1,1,1),theta(1,1),model_sst(1,1),qar3(1,1,1,1),qs(1,1),zz(1,1),vit(1,1,1,1),&
           core_sensible(1,1),core_latent(1,1),core_tau,ustar(1,1,1))

       do j=1,jmt
          do i=1,imt
            sshf(i,j,1)=core_sensible(i,j)*vit(i,j,1,1)*(1.0d0-seaice(i,j,1))
!            lthf(i,j,1)=(core_latent(i,j)-snow3(i,j,1,1)*3.335e+5)*vit(i,j,1,1)*(1.0d0-seaice(i,j,1))
            !mxyu 20251002 gfs snow3 meaning: Snow phase change heat flux, unit: w m-2
            lthf(i,j,1)=(core_latent(i,j)-snow3(i,j,1,1))*vit(i,j,1,1)*(1.0d0-seaice(i,j,1))
            !msxyu 20251002
            lwv(i,j,1)=(0.95*lwv3(i,j,1,1)-0.95*5.67E-8*model_sst(i,j)**4)*vit(i,j,1,1)*(1.0d0-seaice(i,j,1))
            tmp1(i,j)=core_tau(i,j)*windx(i,j)*vit(i,j,1,1)
            tmp2(i,j)=-core_tau(i,j)*windy(i,j)*vit(i,j,1,1)
            nswv(i,j,1)=(lwv(i,j,1)+sshf(i,j,1)+lthf(i,j,1))
            swv(i,j,1)=(1-0.066)*swv3(i,j,1,1)*vit(i,j,1,1)*(1.0d0-seaice(i,j,1))
            if(swv(i,j,1) < 0.0) write(*,*) "swv is negative, i,j swv, swv3, seaice:", i,j,swv(i,j,1),swv3(i,j,1,1),seaice(i,j,1)
            if(runoff(i,j,1) < 0.0) write(*,*) "runoff is negative, i,j, runoff:", i,j,runoff(i,1,1)
            if(rain3(i,j,1,1) < 0.0) then
                write(*,*) "rain3 is negative, set to 0.0 i,j, rain3:", i,j,rain3(i,j,1,1)
                rain3(i,j,1,1) = 0.0
            end if
            !fresh(i,j,1)=-(core_latent(i,j)/(2.5e+6)+rain3(i,j,1,1)+snow3(i,j,1,1)+runoff(i,j,1))*vit(i,j,1,1)*(1.0d0-seaice(i,j,1))
            !mxyu 20251002 gfs rain3: total Precipitation rate
            fresh(i,j,1)=-(core_latent(i,j)/(2.5e+6)+rain3(i,j,1,1)+runoff(i,j,1))*vit(i,j,1,1)*(1.0d0-seaice(i,j,1))
            !mxyu 20251002
            ustar(i,j,1)=ustar(i,j,1)*vit(i,j,1,1)
         end do
      end do

! tau to U/V grid
        do j = jst,jem
           do i = 2,imm
             su(i,j,1)= 0.25*(tmp1(i,j)+tmp1(i-1,j)+tmp1(i,j+1)+tmp1(i-1,j+1))*viv(i,j,1,1)
             sv(i,j,1)= 0.25*(tmp2(i,j)+tmp2(i-1,j)+tmp2(i,j+1)+tmp2(i-1,j+1))*viv(i,j,1,1)
           end do
        end do

!      deallocate(buffer_s)
      return
      end subroutine GFS_INTERPOLATE



   subroutine map_interplation(src_value, local_buffer, block_index, src_coord, weight, max_g_n_weight)
      use precision_mod
      implicit none
      integer(8), intent(in) :: max_g_n_weight
      real(r8), intent(in) :: src_value(g_imt, g_jmt), weight(max_g_n_weight)
      integer(4), intent(in) :: block_index(nx_block, ny_block, 2)
      integer(8), intent(in) :: src_coord(max_g_n_weight)
      real(r8), intent(inout) :: local_buffer(nx_block, ny_block)

      integer :: local_block_num, g_i, g_j, block_start, block_end
      integer :: i, j, m, src_i_id, src_j_id

      do j = 1, ny_block
            do i = 1, nx_block   
                  local_buffer(i, j) = 0.0_r8
                  do m = block_index(i, j, 1), block_index(i, j, 2)
                        src_i_id = mod(src_coord(m) - 1 , g_imt) + 1
                        src_j_id = (src_coord(m) - 1) / g_imt + 1
                        ! if(rank == 0) write(6, *) m, src_i_id, src_j_id
                        local_buffer(i, j) = local_buffer(i, j) &
                        + weight(m) * src_value(src_i_id, src_j_id)
                        
                  enddo
            enddo
      enddo

   end subroutine

!---------------------------------------------



!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
! Over-ocean fluxes following Large and Yeager (used in NCAR models)           !
! Coded by Mike Winton (Michael.Winton@noaa.gov) in 2004
!
! A bug was found by Laurent Brodeau (brodeau@gmail.com) in 2007.
! Stephen.Griffies@noaa.gov updated the code with the bug fix. 
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
!
subroutine ncar_ocean_fluxes_jra (u_del, t, ts, q, qs, z, avail, &
                              sh,lh,tau,ustar)
!                              cd, ch, ce, ustar, bstar       )
use precision_mod
      implicit none
    real(r8)   , intent(in)   , dimension(imt,jmt) :: u_del, t, ts, q, qs, z
    real(r8)   , intent(in)   , dimension(imt,jmt) :: avail
    real(r8)   , intent(inout), dimension(imt,jmt) :: lh,sh,tau
    real(r8)   , dimension(imt,jmt) :: cd, ch, ce, ustar, bstar
!    real   , intent(inout), dimension(imt,jmt) :: cd, ch, ce, ustar, bstar

  real(r8) :: cd_n10, ce_n10, ch_n10, cd_n10_rt    ! neutral 10m drag coefficients
  real(r8) :: cd_rt                                ! full drag coefficients @ z
  real(r8) :: zeta, x2, x, psi_m, psi_h            ! stability parameters
  real(r8) :: u, u10, tv, tstar, qstar, z0, xx, stab

  integer, parameter :: n_itts = 2
  real(r8), parameter :: grav = 9.80, vonkarm = 0.40,  L=2.5e6, cp=1000.5, r0=1.22
  integer               i, j, jj

  cd_n10 = 0.0    
  cd_n10_rt = 0.0    
  ce_n10 = 0.0    
  stab = 0.0    
  cd = 0.0    
  ch = 0.0    
  ce = 0.0

  do j=1,jmt
  do i=1,imt
!  do i=1,size(u_del)
    if (avail(i,j) > 0.5 ) then
      tv = t(i,j)*(1+0.608*q(i,j));
      u = max(u_del(i,j), 0.5);                                 ! 0.5 m/s floor on wind (undocumented NCAR)
      u10 = u;                                                ! first guess 10m wind
    
      cd_n10 = (2.7/u10+0.142+0.0764*u10)/1e3;                ! L-Y eqn. 6a
      cd_n10_rt = sqrt(cd_n10);
      ce_n10 =                     34.6 *cd_n10_rt/1e3;       ! L-Y eqn. 6b
      stab = 0.5 + sign(0.5d0,t(i,j)-ts(i,j))
      ch_n10 = (18.0*stab+32.7*(1-stab))*cd_n10_rt/1e3;       ! L-Y eqn. 6c
  
      cd(i,j) = cd_n10;                                         ! first guess for exchange coeff's at z
      ch(i,j) = ch_n10;
      ce(i,j) = ce_n10;
      do jj=1,n_itts                                           ! Monin-Obukhov iteration
        cd_rt = sqrt(cd(i,j));
        ustar(i,j) = cd_rt*u;                                   ! L-Y eqn. 7a
        tstar    = (ch(i,j)/cd_rt)*(t(i,j)-ts(i,j));                ! L-Y eqn. 7b
        qstar    = (ce(i,j)/cd_rt)*(q(i,j)-qs(i,j));                ! L-Y eqn. 7c
        bstar(i,j) = grav*(tstar/tv+qstar/(q(i,j)+1/0.608));

        zeta     = vonkarm*bstar(i,j)*z(i,j)/(ustar(i,j)*ustar(i,j)); ! L-Y eqn. 8a
        zeta     = sign( min(abs(zeta),10.0), zeta );         ! undocumented NCAR
        x2 = sqrt(abs(1-16*zeta));                            ! L-Y eqn. 8b
        x2 = max(x2, 1.0);                                    ! undocumented NCAR
        x = sqrt(x2);
    
        if (zeta > 0) then
          psi_m = -5*zeta;                                    ! L-Y eqn. 8c
          psi_h = -5*zeta;                                    ! L-Y eqn. 8c
        else
          psi_m = log((1+2*x+x2)*(1+x2)/8)-2*(atan(x)-atan(1.0)); ! L-Y eqn. 8d
          psi_h = 2*log((1+x2)/2);                                ! L-Y eqn. 8e
        end if
    
        u10 = u/(1+cd_n10_rt*(log(z(i,j)/10)-psi_m)/vonkarm);       ! L-Y eqn. 9


        cd_n10 = (2.7/u10+0.142+0.0764*u10)/1e3;                  ! L-Y eqn. 6a again
        cd_n10_rt = sqrt(cd_n10);
        ce_n10 = 34.6*cd_n10_rt/1e3;                              ! L-Y eqn. 6b again
        stab = 0.5 + sign(0.5d0,zeta)
        ch_n10 = (18.0*stab+32.7*(1-stab))*cd_n10_rt/1e3;         ! L-Y eqn. 6c again
        z0 = 10*exp(-vonkarm/cd_n10_rt);                          ! diagnostic
    
        xx = (log(z(i,j)/10)-psi_m)/vonkarm;
        cd(i,j) = cd_n10/(1+cd_n10_rt*xx)**2;                       ! L-Y 10a
        xx = (log(z(i,j)/10)-psi_h)/vonkarm;
!!$        ch(i,j) = ch_n10/(1+ch_n10*xx/cd_n10_rt)**2;                !       b (bug)
!!$        ce(i,j) = ce_n10/(1+ce_n10*xx/cd_n10_rt)**2;                !       c (bug)
        ch(i,j) = ch_n10/(1+ch_n10*xx/cd_n10_rt)*sqrt(cd(i,j)/cd_n10) ! 10b (corrected code aug2007)
        ce(i,j) = ce_n10/(1+ce_n10*xx/cd_n10_rt)*sqrt(cd(i,j)/cd_n10) ! 10c (corrected code aug2007)
      end do
    end if
  end do
  end do

    sh=r0*cp*ch*(t-ts)*u_del
    lh=r0*ce*l*(q-qs)*u_del
    tau=r0*cd*u_del

    return
end subroutine ncar_ocean_fluxes_jra

subroutine read_hour_gfs_forcing(nnn)
    use forc_mod
    implicit none
    integer, intent(in) :: nnn !forcing days
    integer :: irec,couplefreq

    couplefreq = 4  !6hour interval

    irec = (nnn-1) * couplefreq + 1
    if (mytid.eq.0) then
        write(*,*) "read_hour_gfs_forcing, nn,irec:",nnn,irec
    endif

    call read_gfs(irec,"t2m_00.nc",t10,"t2m")
    call read_gfs(irec,"u_10_00.nc",u10,"u_10")
    call read_gfs(irec,"v_10_00.nc",v10,"v_10")
    call read_gfs(irec,"q2m_00.nc",q10,"q2m")
    call read_gfs(irec,"slp_00.nc",slp,"slp")
    call read_gfs(irec,"swdn_00.nc",swhf,"swdn")
    call read_gfs(irec,"lwdn_00.nc",lwhf,"lwdn")
    call read_gfs(irec,"prec_00.nc",precr,"prc")
    if (mytid.eq.0) then
        write(*,*) "read_hour_gfs_forcing,rain minmax:",minval(precr),maxval(precr)
    endif
    call read_gfs_fillvalue(irec,"snow_00.nc",precs,"snow")
    if (mytid.eq.0) then
        write(*,*) "read_hour_gfs_forcing,snow minmax:",minval(precs),maxval(precs)
    endif
    call read_gfs(irec,"runoff_00.nc",rf,"runoff")
    call read_gfs(irec,"ice_00.nc",si,"ice")

end subroutine



subroutine read_gfs(nnn,fname,var,varname)
    use precision_mod
    use pconst_mod, only: s_lon,s_lat
    implicit none
#include <netcdf.inc>
    real (r8), intent(inout), dimension(g_imt,g_jmt,4) :: var
    integer, intent(in) :: nnn
    character (len=*),intent(in) :: fname
    character (len=*),intent(in) :: varname

    integer :: start(3),count(3)
    integer :: ncid,varid,iret,i,j,k
    real(r4),allocatable,dimension(:,:,:) :: ttmp

    allocate (ttmp(g_imt,g_jmt,4))
    start(1)=1;count(1)=g_imt
    start(2)=1;count(2)=g_jmt
    start(3)=nnn;count(3)=4

    iret=nf_open(fname,nf_nowrite,ncid)
    call check_err (iret)
    iret = nf_inq_varid(ncid, trim(varname), varid)
    call check_err(iret)
    iret=nf_get_vara_real(ncid,varid,start,count,ttmp)
    call check_err (iret)

    iret=nf_close(ncid)
    call check_err (iret)

    do k=1,4
    do j=1,g_jmt
       var(:,j,k)=ttmp(:,g_jmt-j+1,k)
    end do
    end do
    !var=ttmp
    deallocate(ttmp)
    return

end subroutine read_gfs

subroutine read_gfs_fillvalue(nnn,fname,var,varname)
        use precision_mod
        use pconst_mod, only: s_lon,s_lat
        implicit none
#include <netcdf.inc>
        real (r8), intent(inout), dimension(g_imt,g_jmt,4) :: var
        integer, intent(in) :: nnn
        character (len=*),intent(in) :: fname
        character (len=*),intent(in) :: varname

        real(r4), parameter :: FILL_VALUE = 1.0D20
        real(r4), parameter :: TOLERANCE = 1.0D15
        integer :: start(3),count(3)
        integer :: ncid,varid,iret,i,j,k
        real(r4),allocatable,dimension(:,:,:) :: ttmp

        allocate (ttmp(g_imt,g_jmt,4))
        start(1)=1;count(1)=g_imt
        start(2)=1;count(2)=g_jmt
        start(3)=nnn;count(3)=4

        iret=nf_open(fname,nf_nowrite,ncid)
        call check_err (iret)
        iret = nf_inq_varid(ncid, trim(varname), varid)
        call check_err(iret)
        iret=nf_get_vara_real(ncid,varid,start,count,ttmp)
        call check_err (iret)

        iret=nf_close(ncid)
        call check_err (iret)

        where (abs(ttmp - FILL_VALUE) < TOLERANCE)
            ttmp = 0.0
        end where

        do k=1,4
        do j=1,g_jmt
           var(:,j,k)=ttmp(:,g_jmt-j+1,k)
        end do
        end do
        !var=ttmp
        deallocate(ttmp)
        return

end subroutine read_gfs_fillvalue

end module gfs_forcing_mod