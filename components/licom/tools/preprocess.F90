!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

program main

   use netcdf

   implicit none

   call preinvert()

end program main

!***********************************************************************

subroutine preinvert()
   
   call invertnc('sss_phc3_monthly_100km.nc', 'salt'            , .true.  , .false., 0)   ! var_id is 7
   call invertnc('tidal_energy_100km.nc'    , 'wave_dissipation', .false. , .false., 0)   ! var_id is 5
   call invertnc('MODEL_CHLFRC_100km'       , 'chl'             , .true.  , .false., 0)   ! var_id is 6
   call invertnc('sss_phc3_monthly_10km.nc' , 'salt'            , .true.  , .false., 1)   ! var_id is 7
   call invertnc('tidal_energy_10km.nc'     , 'wave_dissipation', .false. , .false., 1)   ! var_id is 5
   call invertnc('MODEL_CHLFRC_10km'        , 'chl'             , .true.  , .false., 1)   ! var_id is 6
   call invertnc('TSinitial_100km.nc'       , 'salt'            , .false. , .true. , 0)   ! var_id is 7, from 'TSinitial'
   call invertnc('TSinitial_100km.nc'       , 'temp'            , .false. , .true. , 0)   ! var_id is 8, from 'TSinitial'
   call invertnc('TSinitial_10km.nc'        , 'salt'            , .false. , .true. , 1)   ! var_id is 7, from 'salt_initial'
   call invertnc('TSinitial_10km.nc'        , 'temp'            , .false. , .true. , 1)   ! var_id is 8, from 'temp_initial'

end subroutine preinvert

subroutine invertnc(filenamein, varname, has12month, use_km, res)

   use netcdf

   !-----------------------------------------------------------------------
   !
   !  local variables
   !
   !-----------------------------------------------------------------------
   integer:: var_id
   integer:: imt_global, jmt_global, km
   integer:: i, j, k, mytid
   integer*4   :: ncid, iret
   integer     , intent(in) :: res   ! 0 is 100km, 1 is 10km, 2 is 2km
   logical     , intent(in) :: has12month, use_km
   character(*), intent(in) :: filenamein
   character(*), intent(in) :: varname
   real(kind=8),allocatable,dimension(:,:) :: buffer

   integer(kind=4),  dimension(4) :: start(4)
   integer(kind=4),  dimension(4) :: count(4)
   integer(kind=4),  dimension(3) :: start1
   integer(kind=4),  dimension(3) :: count1
   integer(kind=4),  dimension(2) :: start2
   integer(kind=4),  dimension(2) :: count2
   real(kind=8) xx_r8,yy_r8 

   if(res == 1) then   ! res = 1 for 10km
      imt_global = 3600
      jmt_global = 2302
      km = 55
   else if(res == 2) then ! res = 2 for 2km
      imt_global = 18000
      jmt_global = 11511
      km = 244
   else   ! res = 0 for 100km
      imt_global = 360
      jmt_global = 218
      km = 30
   end if

   write(*,*) "invertnc-->imt_global/jmt_global/km is ",imt_global,jmt_global,km

   allocate(buffer(imt_global,jmt_global))

   iret=nf90_open(trim(filenamein),nf90_write,ncid)
   call check_err (iret)
   iret = nf90_inq_varid (ncid, varname, var_id) 
   call check_err (iret)
   write(*,*) "invertnc-->filenamein/var_id is ",trim(filenamein),var_id
   write(*,*) "invertnc-->has12month/use_km/res is ",has12month,use_km,res

   if(has12month) then
      do k=1, 12    ! add for 'sss_phc3_monthly.nc' and 'MODEL_CHLFRC'
         start1(1)=1 ; count1(1)=imt_global
         start1(2)=1 ; count1(2)=jmt_global
         start1(3)=k ; count1(3)=1

         start(1)=1 ; count(1)=imt_global
         start(2)=1 ; count(2)=jmt_global
         start(3)=1 ; count(3)=1
         start(4)=k ; count(4)=1    

         if(filenamein(1:len_trim('MODEL_CHLFRC')) == 'MODEL_CHLFRC') then
            iret = nf90_get_var(ncid, var_id, buffer, start1, count1)
         else
            iret = nf90_get_var(ncid, var_id, buffer, start , count)
         end if 
         call check_err (iret)
      
         do j=1 ,jmt_global/2
            do i=1, imt_global
               xx_r8 = buffer(i,j)
               yy_r8 = buffer(i,jmt_global+1-j)
               buffer(i,jmt_global+1-j) = xx_r8
               buffer(i,j) = yy_r8
            end do
         end do

         if(filenamein(1:len_trim('MODEL_CHLFRC')) == 'MODEL_CHLFRC') then
            iret = nf90_put_var(ncid, var_id, buffer, start1, count1)
         else
            iret = nf90_put_var(ncid, var_id, buffer, start , count)
         end if     
         CALL check_err (iret)
      end do
   else   ! has12month is false
      if(use_km) then
         do k=1, km    ! add for 'TSinitial'
            start(1)=1 ; count(1)=imt_global
            start(2)=1 ; count(2)=jmt_global
            start(3)=k ; count(3)=1
            start(4)=1 ; count(4)=1

            iret = nf90_get_var(ncid,  var_id, buffer, start ,count)
            call check_err (iret)
         
            do j=1 ,jmt_global/2
               do i=1, imt_global
                  xx_r8 = buffer(i,j)
                  yy_r8 = buffer(i,jmt_global+1-j)
                  buffer(i,jmt_global+1-j) = xx_r8
                  buffer(i,j) = yy_r8
               end do
            end do

            iret = nf90_put_var(ncid, var_id, buffer, start, count)
            CALL check_err (iret)
         end do
      else    ! has12month  and use_km is false
         do k=1,1   ! add for 'tidal_energy.nc'
            start2(1)=1 ; count2(1)=imt_global
            start2(2)=1 ; count2(2)=jmt_global

            iret = nf90_get_var(ncid, var_id, buffer, start2,count2) 
            call check_err (iret)

            do j=1 ,jmt_global/2
               do i=1, imt_global
                  xx_r8 = buffer(i,j)
                  yy_r8 = buffer(i,jmt_global+1-j)
                  buffer(i,jmt_global+1-j) = xx_r8
                  buffer(i,j) = yy_r8
               end do
            end do

            iret = nf90_put_var(ncid, var_id, buffer, start2,count2)
            CALL check_err (iret)
         end do 
      end if
   end if

   iret = nf90_close (ncid)
   call check_err (iret)

   deallocate(buffer)

end subroutine invertnc

subroutine check_err (iret)
   use netcdf

   integer :: iret
   if (iret /= NF90_NOERR) then
      write(*,*)"check_err--> ",nf90_strerror (iret)
      stop
   end if
end subroutine check_err
