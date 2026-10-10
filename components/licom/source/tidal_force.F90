!! ===========================================================================================================================
!! Implementation of tides by computation of the Sun's and Moon's full tidal potential
!! This will be used in the pressure gradient calculation
!! 2024.07.20
!! Authors: TianZhiwei 
!! The method of astronomical is from (Duffett-Smith et al. 2011)
!! The method of tidal potential is from (Pugh 1987)
!! ===========================================================================================================================

#include <def-undef.h>
#ifdef TIDEFORCE

subroutine tide(NC)
use precision_mod
use param_mod
use pconst_mod
use dyn_mod
use work_mod
use tracer_mod
use pmix_mod
use POP_HaloMod
use POP_GridHorzMod
use domain
use grid
use blocks
use operators
use LICOM_Error_mod
use msg_mod
use gather_scatter
use distribution
use constant_mod
    IMPLICIT NONE
    integer :: tide_year,tide_month,tide_day,tide_hour,tide_minute,tide_second
    ! integer :: tide_year0,tide_month0,tide_day0
    integer :: steps_day_bt,tide_steps
    real(r8) :: ndays,gst,decl_m,ra_m,h_m,rs,rm
    real(r8) :: decl_s,ra_s,h_s,mpot,spot,beta2
    real(r8) :: smooth
    integer :: i_tide,j_tide,iblock
    integer :: tide_iteration,tides_smooth_iterations
    integer :: NC,tide_dt,steps_perday_tide

!---------------------------------------------------------------------
    !The real time of the tide in model
    !TZW 2024.08.31
    tide_dt=int(DTB) 

    tide_year = iyfm
    tide_month = month
    tide_day = iday
    
    tide_steps = (II-1)*NBB+NC

    tide_hour = (tide_steps-1)*tide_dt / 3600
    tide_minute = (((tide_steps-1)*tide_dt) - (tide_hour*3600)) / 60
    tide_second = (tide_steps-1)*tide_dt - tide_hour*3600 - tide_minute*60
    ! if (mytid.eq.0) write(16,*)'tidal_force,tide_year=',tide_year,'tide_month=',tide_month,'tide_day=',tide_day
    ! if (mytid.eq.0) write(16,*)'tidal_force,tide_steps=',tide_steps,'tide_hour=',tide_hour,'tide_minute',tide_minute
    ! call flush(16)
!---------------------------------------------------------------------

    beta2 = 0.69 !earth tide factor
    ! tide_iteration = 0
    ! tides_smooth_iterations = 200
    ! smooth = 0.0
    ! IF (tide_iteration >= tides_smooth_iterations) THEN 
    !     smooth = 1.0
    ! ELSE
    !     smooth = REAL(tide_iteration,r8) / REAL(tides_smooth_iterations,r8)
    !     tide_iteration = tide_iteration + 1
    ! ENDIF

    call timing(tide_year,tide_month,tide_day,tide_hour,tide_minute,tide_second,ndays)
    call get_gst(tide_year,tide_month,tide_day,tide_hour,tide_minute,tide_second,gst)
    call sun(ndays,gst,decl_s,ra_s,h_s,rs)
    call moon(ndays,gst,decl_m,ra_m,h_m,rm)
    ! if (mytid==500) write(6,*) "mdecl,sdecl",decl_m,decl_s
    ! if (mytid==500) write(6,*) "h_m,h_s",h_m,h_s
    ! if (mytid==500) write(6,*) "rm,rs",rm,rs
    do iblock = 1, nblocks_clinic
        do j_tide = 1,jmt
            do i_tide = 1,imt
                call Sun_pot(decl_s,h_s,tlat(i_tide,j_tide,iblock),tlon(i_tide,j_tide,iblock),rs,spot)
                call Moon_pot(decl_m,h_m,tlat(i_tide,j_tide,iblock),tlon(i_tide,j_tide,iblock),rm,mpot)
                ! elevat(i_tide,j_tide,iblock) = -1.0*(smooth*beta2*(mpot + spot))/9.8
                elevat(i_tide,j_tide,iblock) = -1.0*(beta2*(mpot + spot))/9.8
            enddo
        enddo
    enddo 

end subroutine tide
!---------------------------------------------------------------------
! get Greenwich siderial time
!---------------------------------------------------------------------
subroutine get_gst(year,month,day,hour,minute,second,gst) 
    use precision_mod
    implicit none
    
    integer  :: year,month,day,hour,minute,second
    real(r8) :: gst,d,s,t0,t,ut
    
    call timing(year,month,day,0,0,0,d)
    
    s = d + 3651.5 
    ! number of days between 2000 Jan 1.5 and 2010 Jan 0.0
    ! 3650 +3(leap year) -1(Jan 0.0) -0.5 (Jan 1.5)
    t = s/36525.0 ! Julian century 
    t0 = 6.697374558 + 2400.051336*t + 0.000025862*t**2
    ! Page 23
    call adj_hour(t0)
    
    ut = hour*1.0 + minute*1.0/60.0 + second*1.0/3600.0 

    ut = ut*1.002737909
    ! Universial time
    gst = t0 + ut
    
    call adj_hour(gst)
    
    return
end subroutine get_gst
!-----------------------------------------------------------------------
! find the number of days (d) since 2010 January 0.0
!-----------------------------------------------------------------------
subroutine timing(year,month,day,hour,minute,second,d)
    use precision_mod
    implicit none
    
    real(r8) :: d,d0
    integer  :: year,month,day,hour,minute,second
    integer ::  molen(12),ylen(1904:2500),y,m
    
    data molen /31,28,31,30,31,30,31,31,30,31,30,31/
    
    do y = 1904,2500,1
      ylen(y) = 365
    enddo
    do y = 1904,2500,4
      ylen(y) = 366
      if(year.eq.y) molen(2) = 29
    enddo
    do y = 2100,2500,100
        if (y.eq.2400) then
          ylen(y) = 366 
  !attention the year 2400 is leap year,ICON-O's source code is wrong.
        else
        ylen(y) = 365
        endif
        if(year.eq.y) molen(2) = 28
        if((year.eq.y) .and. (year.eq.2400)) molen(2) = 29
    enddo
    
    
    d0 = 0.
    do y = 1904,2009,1
      d0 = d0 + ylen(y)
    enddo
    
    d = 0.
    do y = 1904,year-1,1
      d = d + ylen(y)
    enddo
    
    
    do m = 1,month-1,1
      d = d + molen(m)
    enddo
    
    d = d + day
    d = d + hour*1.0/24.0
    d = d + minute*1.0/(24.0*60)
    d = d + second*1.0/(24.0*60*60)
    
    d = d - d0
    
    return
    end subroutine timing

!-------------------------------------------------------------------
! Transform to (0,360) 
!--------------------------------------------------------------------
subroutine adj_angle(x)
    use precision_mod
    implicit none
    
    real(r8) :: x
    
    1 continue
    if(x.lt.0.0) then
        x = x + 360.0
        goto 1
    endif
    
    2 continue
    if(x.gt.360.0) then
        x = x - 360.0
        goto 2
    endif
    
    return
    end subroutine  adj_angle
!------------------------------------------------------------------------
! Transform to (0,24)           
!------------------------------------------------------------------------
subroutine adj_hour(x)
    use precision_mod
    implicit none
    
    real(r8) :: x
    
    1 continue
    if(x.lt.0.0) then
        x = x + 24.0
        goto 1
    endif
    
    2 continue
    if(x.gt.24.0) then
        x = x - 24.0
        goto 2
    endif
    
    return
end subroutine adj_hour
!------------------------------------------------------------------------
subroutine conv(at,x,y)
    use precision_mod
    implicit none
    
    real(r8) :: at,x,y
    
    if(at.lt.0.0) at = at + 360.0
    
    if((x.ge.0.0).and.(y.ge.0.0)) then
      if(at.gt.90.0) at = at - 180.0
    endif
    
    if((x.lt.0.0).and.(y.ge.0.0)) then
      if(at.gt.270.0) at = at - 180.0
    endif
    
    if((x.lt.0.0).and.(y.lt.0.0)) then
      if(at.lt.90.0) at = at + 180.0
    endif
    
    if((x.ge.0.0).and.(y.lt.0.0)) then
      if(at.lt.180.0) at = at + 180.0
    endif
    
    return
end subroutine conv



!-----------------------------------------------------------------------
! Calculate the sun's declination (deg) and right ascension (hours).
! Page 107
!-----------------------------------------------------------------------
subroutine sun(d,gst,decl,ra,h,rs)
    use precision_mod
    implicit none

real(r8) :: d,dcl,alp,gst,h,mel,pl,e,ma,ma_rad,eps,rs
real(r8) :: ea,ta,el,rad,pi,t1,t2,ra,decl,bta,y,x,r0
    !page 107
    !t1 is julian century since 1900 January 0.5.
    !t2 is julian century since 2000 January 1.5
    !d is number of days since 2010 January 0.0
    !mel is mean ecliptical longitude of the sun
    !pl is longitude of the sun at perigee
    !e is eccentricity of the earth's orbit
    !ma is mean anomaly of the sun
    !ea is eccentric anomaly of the sun
    !ta is true anomaly of the sun
    !el is ecliptical longitude of the sun
    !eps is obliquity of the ecliptic
    !ra is right ascension of the sun
    !decl is declination of the sun
    !h is hour angle of the sun
    !rs is distance of the sun
pi = 3.14159265358979323846
rad = pi/180.0

t1 = (d + 3651.5)/36525.0+1.0 !d plus days between 1900 Jan 0.5 and 2010 Jan 0.0

mel = 279.69668 + 36000.76892*t1 + 0.0003025*t1**2
call adj_angle(mel)
pl = 281.22084 + 1.719175*t1 + 0.000452778*t1**2
call adj_angle(pl)
e = 0.01675104 - 0.0000418*t1 - 0.000000126*t1**2

ma = mel - pl
call adj_angle(ma)

ma_rad = ma*rad
ea =  ma_rad

do while(abs(ea-e*sin(ea)-ma_rad) .gt. 1.0e-6)
  ea = ma_rad + e*sin(ea)
enddo
ta = 2.0*atan(sqrt((1.0+e)/(1.0-e))*tan(0.5*ea))
ta = ta/rad
call adj_angle(ta)
el = ta + pl
call adj_angle(el)
!need to convert from ecliptical to equatorial coordinates
!page 51
bta = 0.0  !latitude of the sun is 0
t2 = (d + 3651.5)/36525.0 !julia century since 2000 Jan 1.5
eps = 23.43929166663333336 - 0.013004166666666666*t2 - &
1.6666666666666665e-7*t2**2 + 5.027777777777778e-7*t2**3 ! mean obliquity of the ecliptic

y = sin(rad*el)*cos(rad*eps) - tan(rad*bta)*sin(rad*eps)
x = cos(rad*el)

ra = atan2(y,x)*1.0/rad

call conv(ra,x,y)
ra = ra /15.0 !convert to hours
call adj_hour(ra)
decl = asin(sin(rad*bta)*cos(rad*eps) + &
cos(rad*bta)*sin(rad*eps)*sin(rad*el))*1.0/rad ! declination
call adj_angle(decl)
h =  gst - ra !hour angle  page43
call adj_hour(h)
r0 = 1.495985e11   ! sun's semi-major axis (m)
rs = r0*(1.0-e**2)/(1.0+e*cos(ta*rad)) !page 110
return
end subroutine sun
!-----------------------------------------------------------------------
! Calculate the moon's declination (deg) and right ascension (hours).
! page 164
!-----------------------------------------------------------------------
subroutine moon(d,gst,mdecl,mra,h,rm)
    use precision_mod
    implicit none

real(r8) :: d,gst,ra,h,t1,t2,mel,mmel,l,mma,mn,c,ev,ae,a3,mma_c
real(r8) :: pl,e,ma,ma_rad,ea,ta,el,ec,a4,l_c,v,l_cc,n_c,i,y,x
real(r8) :: at,mbta,mra,mdecl,eps,rad,pi,rm  

pi = 3.14159265358979323846
rad = pi/180.0
!-----------------------------------------------------------------------
! Calculate the sun's declination (deg) and right ascension (hours).
t1 = (d + 3651.5)/36525.0+1.0 !d plus days between 1900 Jan 0.5 and 2010 Jan 0.0

mel = 279.69668 + 36000.76892*t1 + 0.0003025*t1**2
call adj_angle(mel)
pl = 281.22084 + 1.719175*t1 + 0.000452778*t1**2
call adj_angle(pl)
e = 0.01675104 - 0.0000418*t1 - 0.000000126*t1**2

ma = mel - pl
call adj_angle(ma)

ma_rad = ma*rad
ea =  ma_rad

do while(abs(ea-e*sin(ea)-ma_rad) .gt. 1.0e-6)
  ea = ma_rad + e*sin(ea)
enddo
ta = 2.0*atan(sqrt((1.0+e)/(1.0-e))*tan(0.5*ea))
ta = ta/rad
call adj_angle(ta)
el = ta + pl
call adj_angle(el)
!-----------------------------------------------------------------------
!page 164
l=13.1763966*d+91.929336 !mean longitude of the moon
call adj_angle(l)
mma = l-0.1114041*d-130.143076 !mean anomaly of the moon
call adj_angle(mma)
mn = 291.682547-0.0529539*d !moon's mean longitude of the ascending node
call adj_angle(mn)
!-----------------------------------------------------------------------
! correction
c = l-el
ev = 1.2739*sin(rad*(2.0*c-mma))!correction for the evection
ae = 0.1858*sin(ma_rad) !correction for the annual equation
a3 = 0.37*sin(ma_rad) !correction for the third inequality
mma_c = mma + ev - ae - a3 !corrected anomaly
call adj_angle(mma_c)
ec = 6.2886*sin(rad*mma_c) !correction for the equation of the center
a4=0.214*sin(2.0*rad*mma_c) !correction for the fourth inequality
l_c = l + ev + ec - ae + a4 !corrected longitude
v = 0.6583*sin(2.0*rad*(l_c-el)) !correction for the variation
l_cc = l_c + v !corrected longitude
!-----------------------------------------------------------------------
n_c = mn-0.16*sin(ma_rad) !corrected longitude of the ascending node
i = 5.145396 !inclination of the moon's orbit

y = sin(rad*(l_cc-n_c))*cos(rad*i)
x = cos(rad*(l_cc-n_c)) + 1.e-8
at = atan2(y,x)*1.0/rad
call conv(at,x,y)
mmel=at + n_c 
call adj_angle(mmel)
mbta = asin(sin(rad*(l_cc-n_c))*sin(rad*i))*1.0/rad

!need to convert from ecliptical to equatorial coordinates
t2 = (d + 3651.5)/36525.0 !julia century since 2000 Jan 1.5
eps = 23.43929166663333336 - 0.013004166666666666*t2 - &
1.6666666666666665e-7*t2**2 + 5.027777777777778e-7*t2**3 ! mean obliquity of the ecliptic

y = sin(rad*mmel)*cos(rad*eps) - tan(rad*mbta)*sin(rad*eps)
x = cos(rad*mmel)

mra = atan2(y,x)*1.0/rad
call conv(mra,x,y)
mra = mra /15.0 !convert to hours
call adj_hour(mra)
mdecl = asin(sin(rad*mbta)*cos(rad*eps) + &
cos(rad*mbta)*sin(rad*eps)*sin(rad*mmel))*1.0/rad ! declination
call adj_angle(mdecl)
h=gst-mra
call adj_hour(h)


! rho = (1.0-0.0549**2)/(1.0+0.0549*cos((mma_c+ec)*rad)) !page 176
! pai = 0.9507/rho
! rm = 6738140.0/sin(pai*rad) !page 176
! r0m = 384400000.0
rm = 384401000.0*(1.0-0.0549**2)/(1.0+0.0549*cos((mma_c+ec)*rad))
! rm = 6378140.0/sin(pai*rad)
return
end subroutine moon

subroutine Moon_pot(decl,h,lat,lon,r,pot) ! computation of the Moon's tidal potential
    use precision_mod
    implicit none
         
    REAL(r8),intent(inout) :: pot
    REAL(r8),intent(in) :: decl,h, lat,lon,r
    
    REAL(r8) :: pi,rad,gam,m,pot0,er,tet,phi,codec,latd,lond
    REAL(r8) :: cosg,wrz
    
    
    pi = 3.141592653589793
    rad = pi/180.0
    
    
    gam = 6.67408e-11 ! gravitational constant [m^3/kg/s]
    m = 7.349e22 ! Moon's mass [kg]
    er = 6371.e3 ! earth radius [m]
    
    pot0 = gam*m/r
    
    latd = lat/rad
    lond = lon/rad
      
    tet = (90.0 - latd)*rad
    codec = (90.0 - decl)*rad
    phi = (h*360.0/24.0 + lond)*rad
    
      
    cosg = dcos(tet)*dcos(codec) + dsin(tet)*dsin(codec)*dcos(phi)  
    wrz = dsqrt(1.0 - 2.0*er/r*cosg + (er/r)**2)
      
    pot = pot0*(1.0 + er/r*cosg - 1.0/wrz)
      
    
    return
    end subroutine Moon_pot
    
    
    ! ==========================================================================
    
    subroutine Sun_pot(decl,h,lat,lon,r,pot) ! computation of the Sun's tidal potential
        use precision_mod
        use param_mod, only:mytid
    implicit none
         
    REAL(r8),intent(inout) :: pot
    REAL(r8),intent(in) :: decl,h, lat,lon,r

    REAL(r8) :: pi,rad,gam,m,pot0,er,tet,phi,codec,latd,lond
    REAL(r8) :: cosg,wrz
    
    
    pi = 3.141592653589793_r8
    rad = pi/180.0_r8
    
    
    gam = 6.67408e-11_r8 ! gravitaional constant [m^3/kg/s]
    m = 1.9884e30_r8 ! Sun's mass [kg]
    er = 6371.e3_r8 ! earth radius [m]
    
    pot0 = gam*m/r
    
    latd = lat/rad
    lond = lon/rad
      
    tet = (90.0_r8 - latd)*rad
    codec = (90.0_r8 - decl)*rad
    phi = (h*360.0/24.0_r8 + lond)*rad
    

    cosg = dcos(tet)*dcos(codec) + dsin(tet)*dsin(codec)*dcos(phi)  
    wrz = dsqrt(1.0_r8 - 2.0_r8*er/r*cosg + (er/r)**2)
      
    pot = pot0*(1.0_r8 + er/r*cosg - 1.0_r8/wrz)
      
    
    return
    end subroutine Sun_pot


#else
subroutine tide
    return
end subroutine tide
#endif
