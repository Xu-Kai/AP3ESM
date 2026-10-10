
module grist_gcm_io_h0_module

  use grist_constants,                  only: i4
  use grist_hpe_constants,              only: eta_full_a, eta_full_b, eta_face_a, eta_face_b
  use grist_data_types,                 only: scalar_1d_field, scalar_2d_field
  use grist_domain_types,               only: global_domain
  use grist_dycore_vars_module,         only: dycoreVarGeography 

  use grist_nml_module,                 only: outdir, nlev, ntracer,               &
                                              start_ymd, start_tod,                &
                                              working_mode,                        & 
                                              fname_output,                        & 
                                              model_timestep, grid_info
  
  use grist_time_manager,               only: get_curr_cdate

  use grist_gcm_diagnose_h0_module,     only: gcm_h0_accu_physics_variables, &
                                              gcm_h0_dump_physics_variables, &
                                              gcm_h0_rest_physics_variables, &
                                              diag_physics_vars_1d_h0,       &
                                              diag_physics_vars_2d_h0
  use grist_dycore_diagnose_module_2d,  only: dycore_diagnose_variables

  implicit none
  private
  public  :: gcm_output_history_h0

! private
  integer(i4)      :: restart_timestep
  character(len=2) :: old_mon
  character(len=4) :: old_year

  contains

  subroutine gcm_output_history_h0(mesh, itimestep)
! io
  type(global_domain), intent(in)     :: mesh
  integer(i4)        , intent(in)     :: itimestep
! local
  character(len=4)                    :: cyear
  character(len=2)                    :: cmon
  character(len=2)                    :: cday
  character(len=5)                    :: csec
  logical                             :: does_file_exist

  call get_curr_cdate(start_ymd, start_tod, itimestep, model_timestep, &
                      cyear, cmon, cday, csec) 

    inquire(file=trim(outdir)//'step_restart.txt',exist=does_file_exist)
    if(does_file_exist) then
       open(1,file=trim(outdir)//'step_restart.txt',status='old')
       read(1,*) restart_timestep
       close(1)
    else
       restart_timestep = 1
    endif  !if(does_file_exist)

  if(itimestep .eq. 1.or.itimestep.eq.restart_timestep)then
     old_mon = cmon
     old_year= cyear
  end if  

  if(cmon .eq. old_mon)then
    call dycore_diagnose_variables(mesh)
    call gcm_h0_accu_physics_variables
    return
  end if
  
  if(cmon .ne. old_mon)then
     call gcm_h0_dump_physics_variables
     call gcm_mon_output_file(mesh, old_year, old_mon)
     call gcm_h0_rest_physics_variables

     old_mon = cmon
     if(old_year .ne. cyear) old_year= cyear
     
     call dycore_diagnose_variables(mesh)
     call gcm_h0_accu_physics_variables
  end if

  end subroutine gcm_output_history_h0

  subroutine gcm_mon_output_file(mesh, cyear, cmon)
  use grist_util_module,               only: write_string
  use grist_fileio_list_2d_module_par, only: wrap_output_init_2d      ,&
                                             wrap_add_field_2d        ,&
                                             wrap_output_2d_group     ,&
                                             wrap_output_2d_group_sp  ,&
                                             wrap_output_clean_2d
  use grist_fileio_list_1d_module_par, only: wrap_output_init_1d      ,&
                                             wrap_add_field_1d        ,&
                                             wrap_output_1d_group     ,&
                                             wrap_output_1d_group_sp  ,&
                                             wrap_output_clean_1d 


! io
  type(global_domain), intent(in)     :: mesh
  character(len=4),    intent(in)     :: cyear
  character(len=2),    intent(in)     :: cmon
! local
  character(128)                      :: c_glevel
  integer(i4)                         :: ilev
  type(scalar_1d_field)               :: vert_level, vertp_level
  type(scalar_1d_field)               :: hyam, hybm, hyai, hybi
  type(scalar_1d_field)               :: scalar_io_1d_template
  type(scalar_2d_field)               :: scalar_io_2d_template

  if(.not.allocated(vert_level%f))  allocate(vert_level%f(nlev))
  if(.not.allocated(vertp_level%f)) allocate(vertp_level%f(nlev+1))
  if(.not.allocated(hyam%f))        allocate(hyam%f(nlev))
  if(.not.allocated(hybm%f))        allocate(hybm%f(nlev))
  if(.not.allocated(hyai%f))        allocate(hyai%f(nlev+1))
  if(.not.allocated(hybi%f))        allocate(hybi%f(nlev+1))

  vert_level%pos = 8
  vertp_level%pos = 9
  hyam%pos = 8
  hybm%pos = 8
  hyai%pos = 9
  hybi%pos = 9

  do ilev = 1, nlev
    vert_level%f(ilev)  = ilev
    vertp_level%f(ilev) = ilev
  end do
  vertp_level%f(nlev+1) = nlev+1

  hyam%f(:)  = eta_full_a
  hybm%f(:)  = eta_full_b
  hyai%f(:)  = eta_face_a
  hybi%f(:)  = eta_face_b

  call write_string(mesh%glevel, c_glevel)
!
! 2d file
!
#ifndef GB25

  if(trim(grid_info).eq.'null') grid_info = "G"//trim(c_glevel)

  fname_output="GRIST.ATM."//trim(grid_info)//"."//trim(working_mode)//".MonAvg."//trim(cyear)//"-"//trim(cmon)//".2d.h0.nc"

! output
  call wrap_output_init_2d(mesh)

  call wrap_add_field_2d(diag_physics_vars_2d_h0%uwind          ,"uPC",          "zonal wind speed","m/s"   )
  call wrap_add_field_2d(diag_physics_vars_2d_h0%vwind           ,"vPC",          "meridional wind speed","m/s"   )
  call wrap_add_field_2d(diag_physics_vars_2d_h0%omegaFull       ,"omegaFull",    'vertical speed','Pa/s')
  call wrap_add_field_2d(diag_physics_vars_2d_h0%wwwFace         ,"wwwFace" ,     "height-based vertical speed (NDC-only)", "m/s"  )
  call wrap_add_field_2d(diag_physics_vars_2d_h0%phiFace         ,"geopotentialFace", "geopotential at face level", "m^2/s^2")
  call wrap_add_field_2d(diag_physics_vars_2d_h0%mpressureFace   ,"mpressureFace" ,   "moist-mass pressure at face level", "Pa")
  call wrap_add_field_2d(diag_physics_vars_2d_h0%temp            ,"temperature",  "Air temperature","K")
  call wrap_add_field_2d(diag_physics_vars_2d_h0%cloud           ,"cloud3D" ,     "3d cloud fraction","")
  call wrap_add_field_2d(diag_physics_vars_2d_h0%relhum          ,"relhum"  ,     'relative humidity','percent')
  call wrap_add_field_2d(diag_physics_vars_2d_h0%qv              ,"qv"      )
#ifdef AMIPW_PHYSICS
  call wrap_add_field_2d(diag_physics_vars_2d_h0%qc              ,"qc"      )
  call wrap_add_field_2d(diag_physics_vars_2d_h0%qr              ,"qr"      )
  call wrap_add_field_2d(diag_physics_vars_2d_h0%qi              ,"qi"      )
  call wrap_add_field_2d(diag_physics_vars_2d_h0%qs              ,"qs"      )
  call wrap_add_field_2d(diag_physics_vars_2d_h0%qg              ,"qg"      )
#endif

#ifdef SPIO
  call wrap_output_2d_group_sp(mesh, outdir, fname_output)  
#else
  call wrap_output_2d_group(mesh, outdir, fname_output)
#endif

! clean
  call wrap_output_clean_2d()  
#endif

!
! 1d file
!
  if(trim(grid_info).eq.'null') grid_info = "G"//trim(c_glevel)

  fname_output="GRIST.ATM."//trim(grid_info)//"."//trim(working_mode)//".MonAvg."//trim(cyear)//"-"//trim(cmon)//".1d.h0.nc"

! output
  call wrap_output_init_1d(mesh)

  call wrap_add_field_1d(dycoreVarGeography%scalar_lon_at_pc      , "lon_nv")
  call wrap_add_field_1d(dycoreVarGeography%scalar_lat_at_pc      , "lat_nv")
  call wrap_add_field_1d(vert_level                               , "nlev",  'number of full levels','')
  call wrap_add_field_1d(vertp_level                              , "nlevp" ,'number of face levels','')
  call wrap_add_field_1d(hyam                                     , "hyam", 'full-level hybrid-coordinate coefficient, a','')
  call wrap_add_field_1d(hybm                                     , "hybm", 'full-level hybrid-coordinate coefficient, b','')
  call wrap_add_field_1d(hyai                                     , "hyai", 'face-level hybrid-coordinate coefficient, a','')
  call wrap_add_field_1d(hybi                                     , "hybi", 'face-level hybrid-coordinate coefficient, b','')
  call wrap_add_field_1d(diag_physics_vars_1d_h0%ps               , "ps"  , 'surface full-air pressure','Pa')
  call wrap_add_field_1d(diag_physics_vars_1d_h0%ts               , "ts"  , 'surface temperature','K')
  call wrap_add_field_1d(diag_physics_vars_1d_h0%t2m              , "t2m" , '2-meter temperature','K')
  call wrap_add_field_1d(diag_physics_vars_1d_h0%u10m             , "u10m", "10-meter U-wind",'m/s')
  call wrap_add_field_1d(diag_physics_vars_1d_h0%v10m             , "v10m", "10-meter V-wind",'m/s')

  call wrap_add_field_1d(diag_physics_vars_1d_h0%shflx            , "shflx", 'sensible-heat-flux at surface'  ,'W/m^2')
  call wrap_add_field_1d(diag_physics_vars_1d_h0%lhflx            , 'lhflx', 'latent-heat-flux at surface'    ,'W/m^2')
  call wrap_add_field_1d(diag_physics_vars_1d_h0%qflx             , "qflx" , 'water-flux at surface'          ,'kg/m^2/s')

  call wrap_add_field_1d(diag_physics_vars_1d_h0%prect            , "prect", 'total precipitation rate','m/s')

#ifndef GB25
  call wrap_add_field_1d(diag_physics_vars_1d_h0%precc            , "precc" )
  call wrap_add_field_1d(diag_physics_vars_1d_h0%precl            , "precl" )
  call wrap_add_field_1d(diag_physics_vars_1d_h0%snowl            , "precSnowl" )
  call wrap_add_field_1d(diag_physics_vars_1d_h0%grapl            , "precGrapl" )
#endif

  call wrap_add_field_1d(diag_physics_vars_1d_h0%cldtot           , "cldtot" ,'Total cloud fraction' ,'fraction')
  call wrap_add_field_1d(diag_physics_vars_1d_h0%cldlow           , "cldlow" ,'Low-cloud fraction'   ,'fraction')
  call wrap_add_field_1d(diag_physics_vars_1d_h0%cldmed           , "cldmed" ,'Mid-cloud fraction'   ,'fraction')
  call wrap_add_field_1d(diag_physics_vars_1d_h0%cldhgh           , "cldhgh" ,'High-cloud fraction'  ,'fraction')

  call wrap_add_field_1d(diag_physics_vars_1d_h0%swcf             , "swcf"  ,'Short Wave Cloud Forcing:A','W/m^2')
  call wrap_add_field_1d(diag_physics_vars_1d_h0%lwcf             , "lwcf"  ,'Long Wave Cloud Forcing:A','W/m^2')

  call wrap_add_field_1d(diag_physics_vars_1d_h0%flwut            , "flwut" ) 
  call wrap_add_field_1d(diag_physics_vars_1d_h0%flwdt            , "flwdt" )
  call wrap_add_field_1d(diag_physics_vars_1d_h0%flwus            , "flwus" )
  call wrap_add_field_1d(diag_physics_vars_1d_h0%flwds            , "flwds" )
  call wrap_add_field_1d(diag_physics_vars_1d_h0%flwutc           , "flwutc")
  call wrap_add_field_1d(diag_physics_vars_1d_h0%flwdtc           , "flwdtc")
  call wrap_add_field_1d(diag_physics_vars_1d_h0%flwusc           , "flwusc")
  call wrap_add_field_1d(diag_physics_vars_1d_h0%flwdsc           , "flwdsc")

  call wrap_add_field_1d(diag_physics_vars_1d_h0%fswut            , "fswut" ) 
  call wrap_add_field_1d(diag_physics_vars_1d_h0%fswdt            , "fswdt" )
  call wrap_add_field_1d(diag_physics_vars_1d_h0%fswus            , "fswus" )
  call wrap_add_field_1d(diag_physics_vars_1d_h0%fswds            , "fswds" )
  call wrap_add_field_1d(diag_physics_vars_1d_h0%fswutc           , "fswutc")
  call wrap_add_field_1d(diag_physics_vars_1d_h0%fswdtc           , "fswdtc")
  call wrap_add_field_1d(diag_physics_vars_1d_h0%fswusc           , "fswusc")
  call wrap_add_field_1d(diag_physics_vars_1d_h0%fswdsc           , "fswdsc")
#ifdef AMIPW_PHYSICS
  call wrap_add_field_1d(diag_physics_vars_1d_h0%asdir            , "asdir")
  call wrap_add_field_1d(diag_physics_vars_1d_h0%asdif            , "asdif")
  call wrap_add_field_1d(diag_physics_vars_1d_h0%aldir            , "aldir")
  call wrap_add_field_1d(diag_physics_vars_1d_h0%aldif            , "aldif")
  call wrap_add_field_1d(diag_physics_vars_1d_h0%emiss            , "emiss")
  call wrap_add_field_1d(diag_physics_vars_1d_h0%snowhland        , "snowhland")
  call wrap_add_field_1d(diag_physics_vars_1d_h0%qsfc             , "qsfc")

  ! vertical-interpolated
#ifdef UVT_OUT
  call wrap_add_field_1d(diag_physics_vars_1d_h0%zzz200      , "zzz200",'200-hPa geopotential','m^2/s^2')
  call wrap_add_field_1d(diag_physics_vars_1d_h0%zzz500      , "zzz500",'500-hPa geopotential','m^2/s^2')
  call wrap_add_field_1d(diag_physics_vars_1d_h0%zzz700      , "zzz700",'700-hPa geopotential','m^2/s^2')
  call wrap_add_field_1d(diag_physics_vars_1d_h0%zzz850      , "zzz850",'850-hPa geopotential','m^2/s^2')
  call wrap_add_field_1d(diag_physics_vars_1d_h0%relhum200   , "relhum200",'200-hPa relative humidity','percent')
  call wrap_add_field_1d(diag_physics_vars_1d_h0%relhum500   , "relhum500",'500-hPa relative humidity','percent')
  call wrap_add_field_1d(diag_physics_vars_1d_h0%relhum700   , "relhum700",'700-hPa relative humidity','percent')
  call wrap_add_field_1d(diag_physics_vars_1d_h0%relhum850   , "relhum850",'850-hPa relative humidity','percent')
  call wrap_add_field_1d(diag_physics_vars_1d_h0%omega200    , "omega200" ,'200-hPa pressure vertical speed','Pa/s')
  call wrap_add_field_1d(diag_physics_vars_1d_h0%omega500    , "omega500" ,'500-hPa pressure vertical speed','Pa/s')
  call wrap_add_field_1d(diag_physics_vars_1d_h0%omega700    , "omega700" ,'700-hPa pressure vertical speed','Pa/s')
  call wrap_add_field_1d(diag_physics_vars_1d_h0%omega850    , "omega850" ,'850-hPa pressure vertical speed','Pa/s')
  call wrap_add_field_1d(diag_physics_vars_1d_h0%u200        ,"u200","200-hPa u speed", 'm/s')
  call wrap_add_field_1d(diag_physics_vars_1d_h0%u500        ,"u500","500-hPa u speed", 'm/s')
  call wrap_add_field_1d(diag_physics_vars_1d_h0%u700        ,"u700","700-hPa u speed", 'm/s')
  call wrap_add_field_1d(diag_physics_vars_1d_h0%u850        ,"u850","850-hPa u speed", 'm/s')

  call wrap_add_field_1d(diag_physics_vars_1d_h0%v200        ,"v200","200-hPa v speed", 'm/s')
  call wrap_add_field_1d(diag_physics_vars_1d_h0%v500        ,"v500","500-hPa v speed", 'm/s')
  call wrap_add_field_1d(diag_physics_vars_1d_h0%v700        ,"v700","700-hPa v speed", 'm/s')
  call wrap_add_field_1d(diag_physics_vars_1d_h0%v850        ,"v850","850-hPa v speed", 'm/s')

  call wrap_add_field_1d(diag_physics_vars_1d_h0%t200        ,"t200","200-hPa temperature", 'K')
  call wrap_add_field_1d(diag_physics_vars_1d_h0%t500        ,"t500","500-hPa temperature", 'K')
  call wrap_add_field_1d(diag_physics_vars_1d_h0%t700        ,"t700","700-hPa temperature", 'K')
  call wrap_add_field_1d(diag_physics_vars_1d_h0%t850        ,"t850","850-hPa temperature", 'K')
#endif
#endif
!AMIPW_PHYSICS

  !call wrap_add_field_1d(diag_physics_vars_h0%z500                   , "z500"  )
  !call wrap_add_field_1d(diag_physics_vars_h0%u200                   , "u200"  )
  !call wrap_add_field_1d(diag_physics_vars_h0%u850                   , "u850"  )
  !call wrap_add_field_1d(diag_physics_vars_h0%v850                   , "v850"  )
#ifdef SPIO
  call wrap_output_1d_group_sp(mesh, outdir, fname_output)  
#else
  call wrap_output_1d_group(mesh, outdir, fname_output)
#endif

! clean
  call wrap_output_clean_1d()  

  if(allocated(vert_level%f))  deallocate(vert_level%f)
  if(allocated(vertp_level%f)) deallocate(vertp_level%f)
  if(allocated(hyam%f))        deallocate(hyam%f)
  if(allocated(hybm%f))        deallocate(hybm%f)
  if(allocated(hyai%f))        deallocate(hyai%f)
  if(allocated(hybi%f))        deallocate(hybi%f)

  end subroutine gcm_mon_output_file

 end module grist_gcm_io_h0_module
