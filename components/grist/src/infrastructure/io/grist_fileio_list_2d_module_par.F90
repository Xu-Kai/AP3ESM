!-----------------------------------------------------------
! Created on 2018
! Authors/Contact: Liu Zhuang
! Version 1.0
! Description: File Par I/O module tailored for GRIST. 
! It deals with the general conditions, must be called by an 
! explicit suit of statements.
!            How to use:
!             1) wrap_output_init_2d
!             2) wrap_add_field_2d
!             3) wrap_output_2d
!             4) wrap_output_clean_2d
!-----------------------------------------------------------

 module grist_fileio_list_2d_module_par
   use grist_lib
   use grist_constants,      only: i4, i8, r4, r8
   use grist_domain_types  , only: global_domain, block_structure, group_comm
   use grist_data_types    , only: scalar_2d_field, data2d_temp_sp, data2d_temp_dp
   use grist_handle_error,   only: endrun
   use grist_nml_module,     only: nlev
   use grist_wrap_pf,        only: wrap_open,                 &
                                   wrap_create,               &
                                   wrap_enddef,               &
                                   wrap_close,                &
                                   wrap_inq_dim,              &
                                   wrap_inq_dimid,            &
                                   wrap_inq_varid,            &
                                   wrap_inq_vartype,          &
                                   wrap_put_vara_realx,       &
                                   wrap_put_vara_realx_req,   &
                                   wrap_put_vara_real,        &
                                   wrap_put_vara_real_req,    &
                                   wrap_put_att_text,         &
                                   wrap_def_var,              &
                                   wrap_def_dim,              &
                                   wrap_get_vara_realx_req,   &
                                   wrap_get_vara_realx,       &
                                   wrap_get_vara_real_req,    &
                                   wrap_get_vara_real

   use grist_list_array
   
   implicit none

   include 'pnetcdf.inc'

   private

   public :: wrap_output_init_2d      ,&
             wrap_add_field_2d        ,&
             wrap_output_2d_group     ,&
             wrap_output_2d_group_sp  ,&
             wrap_output_clean_2d     ,&
             wrap_read_2d_group       ,& ! inputdata: ndim*ncell, outputdata, ncell*ndim
             wrap_read_2d_group_rst   ,& ! inputdata: ncell*ndim, outputdata, ndim*ncell
             wrap_read_2d             ,&
             wrap_bcast_2d            ,&
             wrap_output_2d_gather_C  ,&
             wrap_output_after_process_2d

   interface
     subroutine wrap_output_2d_gather_sp_c(temp, len_dim1, local_cols, group_id, all_comm_f, outdir, filename, open_mode, &
        my_count, my_start, displs_send, counts_recv, displs_recv, idx_recv, varname) bind(C, name="wrap_output_2d_gather_sp_c")
        use, intrinsic                                                  :: iso_c_binding, only: c_int, c_float, c_char
        real(c_float), intent(in)                                       :: temp(*)
        integer(c_int), value                                           :: len_dim1, local_cols, group_id, all_comm_f, open_mode
        character(kind=c_char), dimension(*), intent(in)                :: outdir
        character(kind=c_char), dimension(*), intent(in)                :: filename
        integer(c_int), value                                           :: my_count, my_start
        integer(c_int), dimension(*), intent(in)                        :: displs_send
        integer(c_int), dimension(*), intent(in)                        :: counts_recv
        integer(c_int), dimension(*), intent(in)                        :: displs_recv
        integer(c_int), dimension(*), intent(in)                        :: idx_recv
        character(kind=c_char), dimension(*), intent(in)                :: varname
     end subroutine wrap_output_2d_gather_sp_c
     subroutine bin_to_nc_group_sp(filepath, comm, color, nvars, varnames) bind(C, name="bin_to_nc_group_sp")
        use, intrinsic                                                  :: iso_c_binding, only: c_int, c_char
        character(kind=c_char), dimension(*), intent(in)                :: filepath
        integer(c_int), value                                           :: comm
        integer(c_int), value                                           :: color
        integer(c_int), value                                           :: nvars
        character(kind=c_char), dimension(*), intent(in)                :: varnames
     end subroutine bin_to_nc_group_sp


     subroutine create_directory(outdir, filename, varname, outpath, outlen) bind(C, name="create_directory")
          use, intrinsic                                                  :: iso_c_binding, only: c_int, c_char
          character(kind=c_char), dimension(*), intent(in)                :: outdir
          character(kind=c_char), dimension(*), intent(in)                :: filename
          character(kind=c_char), dimension(*), intent(in)                :: varname
          character(kind=c_char), dimension(*), intent(out)               :: outpath
          integer(c_int), value                                           :: outlen
     end subroutine create_directory

    ! subroutine C_after_process_2d() bind(C, name="C_after_process_2d")
    !   use, intrinsic                                                  :: iso_c_binding, only: c_int, c_char
    !   character(kind=c_char), dimension(*), intent(in)                :: outdir
    !   character(kind=c_char), dimension(*), intent(in)                :: filename
    !   integer(c_int), value                                           :: nlevi8
    !   integer(c_int), value                                           :: nlevp
    !   integer(c_int), value                                           :: nvars
    !   type(c_ptr), value                                             :: scalar_2d_list
    ! end subroutine C_after_process_2d

   end interface

   type(array_2d_index_list), pointer   :: scalar_2d_list   ! any 2d shape
   integer(i4)                          :: number_of_scalar
   integer(i4)                          :: mesh_maxvnb
   integer(i8)                          :: mesh_nt
   integer(i8)                          :: mesh_ne
   integer(i8)                          :: mesh_nv
   integer(i4), PARAMETER               :: MAXLEN = 268000000 !required by pnetcdf (< 2GB/8B)
   real(8)                              :: default_8 = 1D38
   real(4)                              :: default_4 = 1e38

 contains

   subroutine wrap_output_init_2d(mesh)

      type(global_domain), intent(in) :: mesh

      mesh_maxvnb                 = mesh%maxvnb
      mesh_nt                     = mesh%nt_all
      mesh_ne                     = mesh%ne_all
      mesh_nv                     = mesh%nv_all
      number_of_scalar            = 0

      return
   end subroutine wrap_output_init_2d

   subroutine wrap_add_field_2d(scalar,scalar_name,longname,units)

      type(scalar_2d_field), intent(in)            :: scalar        ! 2d array
      character*(*)        , intent(in)            :: scalar_name
      character*(*)        , intent(in), optional  :: longname
      character*(*)        , intent(in), optional  :: units

      if(present(longname).and.present(units))then
         call insert_var_2d(scalar_2d_list, scalar%f, trim(scalar_name), trim(longname),trim(units),number_of_scalar, scalar%pos)
      else
         call insert_var_2d(scalar_2d_list, scalar%f, trim(scalar_name), '','',number_of_scalar, scalar%pos)
      end if

      return
   end subroutine wrap_add_field_2d

   subroutine wrap_output_2d_group(mesh, outdir, filename)

   implicit none
!io
     type(global_domain), intent(in), target :: mesh
     character*(*),       intent(in)         :: outdir
     character*(*),       intent(in)         :: filename
! local
     type(group_comm),          pointer      :: comm
     type(block_structure),     pointer      :: local_block
     type(array_2d_index_list), pointer      :: ptr
     integer(i4)                             :: cmode, ncid, err
     integer(i4)                             :: fst_dim_id(6)           ! dim id for the 1st dim (location)
     integer(i4)                             :: sec_dim_id(3)           ! dim id for the 2nd dim
     integer(i4)                             :: len_dim1, mpi_newtype
     integer(i4)                             :: i0, i1, i2 ,i, j, lidx
     integer(i8)                             :: nlevp, nlevi8
     real(r8),                  allocatable  :: buff_var(:), temp_var(:,:)
     real(8),                   allocatable  :: buff_put_dp(:,:) 
     real(4),                   allocatable  :: buff_put_sp(:,:) 
     integer(i4),               pointer      :: displs_send(:), counts_send(:)
     integer(i4),               pointer      :: displs_recv(:), counts_recv(:)
     integer(i4),               pointer      :: idx_recv(:), idx_local(:)
     integer(i4),               pointer      :: mystart, mycount
     integer(i8)                             :: start(2)
     integer(i8)                             :: count(2)
     integer(i4)                             :: req(1), st(1), act_put_num
     integer(i4)                             :: max_diml, jres, jstart

    !  integer(i4)                             :: root_rank,err
    !  integer(i4),               allocatable  :: recvcounts(:), displs(:)

     nlevi8 = nlev
     nlevp= nlevi8+1
     comm => mesh%gcomm_write
     local_block => mesh%local_block
!================================================
!                 create file
!================================================
     cmode = IOR(NF_CLOBBER,NF_64BIT_DATA)
     if(comm%color==0) call wrap_create(comm%gcomm,trim(outdir)//trim(filename),cmode, ncid) ! large-file support 
     if(mpi_rank() .eq. 0) print*,trim(outdir)//trim(filename)
!================================================
! 1st: define dim
!================================================
! possible dims for 1st
     if(comm%color==0) then
       call wrap_def_dim(ncid, 'nlev'       , nlevi8      , fst_dim_id(1)) ! gcm
       call wrap_def_dim(ncid, 'nlevp'      , nlevi8+1    , fst_dim_id(2)) ! gcm
       call wrap_def_dim(ncid, 'dim12'      , int8(12)    , fst_dim_id(3)) ! gcm    
       call wrap_def_dim(ncid, 'dim3'       , int8(3)     , fst_dim_id(4)) ! lsm restart, added by yiz 2020dec
       call wrap_def_dim(ncid, 'dim4'       , int8(4)     , fst_dim_id(5)) ! lsm restart, added by yiz 2020dec
       call wrap_def_dim(ncid, 'dim7'       , int8(7)     , fst_dim_id(6)) ! lsm restart, added by yiz 2020dec
! possible dims for 2nd
       call wrap_def_dim(ncid, 'location_nt', mesh_nt     , sec_dim_id(1)) ! gcm
       call wrap_def_dim(ncid, 'location_ne', mesh_ne     , sec_dim_id(2)) ! gcm
       call wrap_def_dim(ncid, 'location_nv', mesh_nv     , sec_dim_id(3)) ! gcm
     endif
!========================================================
! 2nd: define var, need lon_id, lat_id, var_idlist
!========================================================
     ptr => scalar_2d_list
     do while (associated(ptr))
       if(ptr%var_level.eq.nlevi8)then
          i1 = 1
          if(ptr%var_pos.eq.1) i2 = 1
          if(ptr%var_pos.eq.6) i2 = 2
          if(ptr%var_pos.eq.0) i2 = 3
       end if
       if(ptr%var_level.eq.nlevp)then
          i1 = 2
          if(ptr%var_pos.eq.1) i2 = 1
          if(ptr%var_pos.eq.6) i2 = 2
          if(ptr%var_pos.eq.0) i2 = 3
       end if
       if(ptr%var_level.eq.12)then
          i1 = 3
          if(ptr%var_pos.eq.1) i2 = 1
          if(ptr%var_pos.eq.6) i2 = 2
          if(ptr%var_pos.eq.0) i2 = 3
       end if      
! lsm restart, added by yiz 2020dec
       if(ptr%var_level.eq.3)then
          i1 = 4
          if(ptr%var_pos.eq.1) i2 = 1
          if(ptr%var_pos.eq.6) i2 = 2
          if(ptr%var_pos.eq.0) i2 = 3
       end if      
       if(ptr%var_level.eq.4)then
          i1 = 5
          if(ptr%var_pos.eq.1) i2 = 1
          if(ptr%var_pos.eq.6) i2 = 2
          if(ptr%var_pos.eq.0) i2 = 3
       end if      
       if(ptr%var_level.eq.7)then
          i1 = 6
          if(ptr%var_pos.eq.1) i2 = 1
          if(ptr%var_pos.eq.6) i2 = 2
          if(ptr%var_pos.eq.0) i2 = 3
       end if      
! lsm restart, added by yiz 2020dec
       if(comm%color==0) call def_var_2d_pnetcdf(ptr%varname, fst_dim_id(i1),sec_dim_id(i2), ncid, ptr%var_idlist)
       ptr => ptr%next
     end do

!================================================
! 3rd: put att
!================================================    
     ptr => scalar_2d_list
     do while (associated(ptr))
       if(comm%color==0) call put_att_2d_pnetcdf( ncid, ptr%var_idlist ,ptr%varname,ptr%longname,ptr%units,ptr%var_pos)
       ptr => ptr%next
     end do 
     if(comm%color==0) call wrap_enddef(ncid)

!================================================
! 4th, put var
!================================================

     ptr => scalar_2d_list
     do while (associated(ptr))      
       if (ptr%var_pos .eq. 0) then
         displs_send => comm%displs_sendv
         displs_recv => comm%displs_recvv
         counts_send => comm%counts_sendv
         counts_recv => comm%counts_recvv
         idx_recv    => comm%idx_recvv
         idx_local   => comm%idx_localv
         mystart     => comm%mystart_v
         mycount     => comm%mycount_v
       else if (ptr%var_pos .eq. 1) then 
         displs_send => comm%displs_sendt
         displs_recv => comm%displs_recvt
         counts_send => comm%counts_sendt
         counts_recv => comm%counts_recvt
         idx_recv    => comm%idx_recvt
         idx_local   => comm%idx_localt
         mystart     => comm%mystart_t
         mycount     => comm%mycount_t
       else if (ptr%var_pos .eq. 6) then 
         displs_send => comm%displs_sende
         displs_recv => comm%displs_recve
         counts_send => comm%counts_sende
         counts_recv => comm%counts_recve
         idx_recv    => comm%idx_recve
         idx_local   => comm%idx_locale
         mystart     => comm%mystart_e
         mycount     => comm%mycount_e
       end if
 
       len_dim1 = ptr%var_level
       allocate(temp_var(len_dim1, size(idx_local)))
       call find_output_var_2d(ptr%var, ptr%var_level, idx_local, temp_var)

       if (comm%color == 0) then
         allocate(buff_var(1:mycount*len_dim1))
       else
         allocate(buff_var(1))
       endif

       if (r8 .eq. 8) then
         call MPI_Type_contiguous(len_dim1, MPI_REAL8, mpi_newtype, err)
       elseif (r8 .eq. 4) then
         call MPI_Type_contiguous(len_dim1, MPI_REAL, mpi_newtype, err)
       end if
       call MPI_Type_commit(mpi_newtype, err)

       call MPI_Alltoallv(temp_var, counts_send, displs_send, mpi_newtype, buff_var, &
                          counts_recv, displs_recv, mpi_newtype, comm%all_comm, err)
       deallocate(temp_var)

       if (comm%color == 0) then
         ! get the actual number of outputs (due to the limitation of pnetcdf)
         if (r8 .eq. 8) then
           max_diml = MAXLEN/len_dim1
         elseif (r8 .eq. 4) then
           max_diml = 2*MAXLEN/len_dim1
         end if
         act_put_num = ceiling(mycount/real(max_diml, r8))

         if (r8 .eq. 8) then
           allocate(buff_put_dp(len_dim1, mycount), source=default_8)
           do i = 1, size(idx_recv)
              lidx = idx_recv(i)
              buff_put_dp(1:len_dim1,lidx) = buff_var(len_dim1*(i-1)+1:len_dim1*i)
           end do
         elseif (r8 .eq. 4) then
           allocate(buff_put_sp(len_dim1, mycount), source=default_4)
           do i = 1, size(idx_recv)
              lidx = idx_recv(i)
              buff_put_sp(1:len_dim1,lidx) = buff_var(len_dim1*(i-1)+1:len_dim1*i)
           end do
         end if

         start(1) = 1
         count(1) = len_dim1
         ! may need to put multiple times (each time < 2GB)
         jres = mod(mycount, max_diml)
         do j = 1, act_put_num
           if (j .lt. act_put_num) then
               count(2) = max_diml
           else
               count(2) = jres
           end if
           start(2) = mystart + (j-1)*max_diml
           jstart   = 1 + (j-1)*max_diml
           if (r8 .eq. 8) then
             call wrap_put_vara_realx_req(ncid, ptr%var_idlist, start(1:2), count(1:2), &
                                          buff_put_dp(:, jstart:jstart+count(2)-1),req)
           elseif (r8 .eq. 4) then
             call wrap_put_vara_real_req(ncid, ptr%var_idlist, start(1:2), count(1:2), &
                                         buff_put_sp(:, jstart:jstart+count(2)-1),req)
           end if
           err = nfmpi_wait_all(ncid, 1, req, st)
         end do
         if (r8 .eq. 8) then
           deallocate(buff_put_dp)
         elseif (r8 .eq. 4) then
           deallocate(buff_put_sp)
         end if
       endif
       deallocate(buff_var)
       ptr => ptr%next
     end do

     if(comm%color==0) call wrap_close(ncid)
    !  buildtest
    !  write(*,*) '[MARK] wrap_output_2d_group: returning'
     return
   end subroutine wrap_output_2d_group
 
   subroutine wrap_output_2d_group_sp(mesh, outdir, filename)

   implicit none
!io
     type(global_domain), intent(in), target :: mesh
     character*(*),       intent(in)         :: outdir
     character*(*),       intent(in)         :: filename
! local
     type(group_comm),          pointer      :: comm
     type(block_structure),     pointer      :: local_block
     type(array_2d_index_list), pointer      :: ptr
     integer(i4)                             :: cmode, ncid, err
     integer(i4)                             :: fst_dim_id(6)           ! dim id for the 1st dim (location)
     integer(i4)                             :: sec_dim_id(3)           ! dim id for the 2nd dim
     integer(i4)                             :: len_dim1, mpi_newtype
     integer(i4)                             :: i0, i1, i2 ,i, j, lidx
     integer(i8)                             :: nlevp, nlevi8
     real(r8),                  allocatable  :: buff_var(:), temp_var(:,:)
     real(4),                   allocatable  :: buff_put_sp(:,:) 
     integer(i4),               pointer      :: displs_send(:), counts_send(:)
     integer(i4),               pointer      :: displs_recv(:), counts_recv(:)
     integer(i4),               pointer      :: idx_recv(:), idx_local(:)
     integer(i4),               pointer      :: mystart, mycount
     integer(i8)                             :: start(2)
     integer(i8)                             :: count(2)
     integer(i4)                             :: req(1), st(1), act_put_num
     integer(i4)                             :: max_diml, jres, jstart

     nlevi8 = nlev
     nlevp= nlevi8+1
     comm => mesh%gcomm_write
     local_block => mesh%local_block
!================================================
!                 create file
!================================================
     cmode = IOR(NF_CLOBBER,NF_64BIT_DATA)
     if(comm%color==0) call wrap_create(comm%gcomm,trim(outdir)//trim(filename),cmode, ncid) ! large-file support 
     if(mpi_rank() .eq. 0) print*,trim(outdir)//trim(filename)
!================================================
! 1st: define dim
!================================================
! possible dims for 1st
     if(comm%color==0) then
       call wrap_def_dim(ncid, 'nlev'       , nlevi8      , fst_dim_id(1)) ! gcm
       call wrap_def_dim(ncid, 'nlevp'      , nlevi8+1    , fst_dim_id(2)) ! gcm
       call wrap_def_dim(ncid, 'dim12'      , int8(12)    , fst_dim_id(3)) ! gcm    
       call wrap_def_dim(ncid, 'dim3'       , int8(3)     , fst_dim_id(4)) ! lsm restart, added by yiz 2020dec
       call wrap_def_dim(ncid, 'dim4'       , int8(4)     , fst_dim_id(5)) ! lsm restart, added by yiz 2020dec
       call wrap_def_dim(ncid, 'dim7'       , int8(7)     , fst_dim_id(6)) ! lsm restart, added by yiz 2020dec
! possible dims for 2nd
       call wrap_def_dim(ncid, 'location_nt', mesh_nt     , sec_dim_id(1)) ! gcm
       call wrap_def_dim(ncid, 'location_ne', mesh_ne     , sec_dim_id(2)) ! gcm
       call wrap_def_dim(ncid, 'location_nv', mesh_nv     , sec_dim_id(3)) ! gcm
     endif
!========================================================
! 2nd: define var, need lon_id, lat_id, var_idlist
!========================================================
     ptr => scalar_2d_list
     do while (associated(ptr))
       if(ptr%var_level.eq.nlevi8)then
          i1 = 1
          if(ptr%var_pos.eq.1) i2 = 1
          if(ptr%var_pos.eq.6) i2 = 2
          if(ptr%var_pos.eq.0) i2 = 3
       end if
       if(ptr%var_level.eq.nlevp)then
          i1 = 2
          if(ptr%var_pos.eq.1) i2 = 1
          if(ptr%var_pos.eq.6) i2 = 2
          if(ptr%var_pos.eq.0) i2 = 3
       end if
       if(ptr%var_level.eq.12)then
          i1 = 3
          if(ptr%var_pos.eq.1) i2 = 1
          if(ptr%var_pos.eq.6) i2 = 2
          if(ptr%var_pos.eq.0) i2 = 3
       end if      
! lsm restart, added by yiz 2020dec
       if(ptr%var_level.eq.3)then
          i1 = 4
          if(ptr%var_pos.eq.1) i2 = 1
          if(ptr%var_pos.eq.6) i2 = 2
          if(ptr%var_pos.eq.0) i2 = 3
       end if      
       if(ptr%var_level.eq.4)then
          i1 = 5
          if(ptr%var_pos.eq.1) i2 = 1
          if(ptr%var_pos.eq.6) i2 = 2
          if(ptr%var_pos.eq.0) i2 = 3
       end if      
       if(ptr%var_level.eq.7)then
          i1 = 6
          if(ptr%var_pos.eq.1) i2 = 1
          if(ptr%var_pos.eq.6) i2 = 2
          if(ptr%var_pos.eq.0) i2 = 3
       end if      
! lsm restart, added by yiz 2020dec
       if(comm%color==0) call def_var_2d_pnetcdf_sp(ptr%varname, fst_dim_id(i1),sec_dim_id(i2), ncid, ptr%var_idlist)
       ptr => ptr%next
     end do

!================================================
! 3rd: put att
!================================================    
     ptr => scalar_2d_list
     do while (associated(ptr))
       if(comm%color==0) call put_att_2d_pnetcdf( ncid, ptr%var_idlist ,ptr%varname,ptr%longname,ptr%units,ptr%var_pos)
       ptr => ptr%next
     end do 
     if(comm%color==0) call wrap_enddef(ncid)

!================================================
! 4th, put var
!================================================

     ptr => scalar_2d_list
     do while (associated(ptr))      
       if (ptr%var_pos .eq. 0) then
         displs_send => comm%displs_sendv
         displs_recv => comm%displs_recvv
         counts_send => comm%counts_sendv
         counts_recv => comm%counts_recvv
         idx_recv    => comm%idx_recvv
         idx_local   => comm%idx_localv
         mystart     => comm%mystart_v
         mycount     => comm%mycount_v
       else if (ptr%var_pos .eq. 1) then 
         displs_send => comm%displs_sendt
         displs_recv => comm%displs_recvt
         counts_send => comm%counts_sendt
         counts_recv => comm%counts_recvt
         idx_recv    => comm%idx_recvt
         idx_local   => comm%idx_localt
         mystart     => comm%mystart_t
         mycount     => comm%mycount_t
       else if (ptr%var_pos .eq. 6) then 
         displs_send => comm%displs_sende
         displs_recv => comm%displs_recve
         counts_send => comm%counts_sende
         counts_recv => comm%counts_recve
         idx_recv    => comm%idx_recve
         idx_local   => comm%idx_locale
         mystart     => comm%mystart_e
         mycount     => comm%mycount_e
       end if

       len_dim1 = ptr%var_level
       allocate(temp_var(len_dim1, size(idx_local)))
       call find_output_var_2d(ptr%var, ptr%var_level, idx_local, temp_var)

       if (comm%color == 0) then
         allocate(buff_var(1:mycount*len_dim1))
       else
         allocate(buff_var(1))
       endif

       if (r8 .eq. 8) then
         call MPI_Type_contiguous(len_dim1, MPI_REAL8, mpi_newtype, err)
       elseif (r8 .eq. 4) then
         call MPI_Type_contiguous(len_dim1, MPI_REAL, mpi_newtype, err)
       end if
       call MPI_Type_commit(mpi_newtype, err)

        ! if (mpi_rank() == 16) then
        !     do i = 1, size(idx_local)
        !         do j = 1, len_dim1
        !             write(*,*) 'temp_var(', j, ',', i, ') = ', temp_var(j, i)
        !         end do
        !     end do
        ! end if

      ! === Alltoall 各进程数据交换 ===
       call MPI_Alltoallv(temp_var, counts_send, displs_send, mpi_newtype, buff_var, &
                          counts_recv, displs_recv, mpi_newtype, comm%all_comm, err)
       deallocate(temp_var)

      !  call MPI_Barrier(MPI_COMM_WORLD, err)  ! 确保所有进程通信完成

       if (comm%color == 0) then
         ! get the actual number of outputs (due to the limitation of pnetcdf)
         max_diml = 2*MAXLEN/len_dim1
         act_put_num = ceiling(mycount/real(max_diml, r8))

         allocate(buff_put_sp(len_dim1, mycount), source=default_4)
         do i = 1, size(idx_recv)
            lidx = idx_recv(i)
            buff_put_sp(1:len_dim1,lidx) = buff_var(len_dim1*(i-1)+1:len_dim1*i)
         end do

         start(1) = 1
         count(1) = len_dim1
         ! may need to put multiple times (each time < 2GB)
         jres = mod(mycount, max_diml)
         do j = 1, act_put_num
           if (j .lt. act_put_num) then
               count(2) = max_diml
           else
               count(2) = jres
           end if
           start(2) = mystart + (j-1)*max_diml
           jstart   = 1 + (j-1)*max_diml
           call wrap_put_vara_real_req(ncid, ptr%var_idlist, start(1:2), count(1:2), &
                                       buff_put_sp(:, jstart:jstart+count(2)-1),req)
           err = nfmpi_wait_all(ncid, 1, req, st)
         end do
         deallocate(buff_put_sp)
       endif
       deallocate(buff_var)
       ptr => ptr%next
     end do

     if(comm%color==0) call wrap_close(ncid)
     return
   end subroutine wrap_output_2d_group_sp

!================================================
! 5th: gather + write to binary file
!================================================
   subroutine wrap_output_2d_gather_C(mesh, outdir, filename)
      use, intrinsic :: iso_c_binding, only: c_int, c_float, c_char, c_null_char
      use mpi
      use perf_mod
      implicit none
      ! io
      type(global_domain), intent(in), target :: mesh
      character(*),       intent(in)          :: outdir
      character(*),       intent(in)          :: filename
      ! local
      type(group_comm),          pointer      :: comm
      type(array_2d_index_list), pointer      :: ptr
      integer(i4),               pointer      :: idx_local(:), idx_recv(:)
      integer(c_int),            allocatable  :: idx_recv_c(:)
      integer(i4),               pointer      :: displs_send(:), counts_send(:)
      integer(i4),               pointer      :: displs_recv(:), counts_recv(:)
      integer(c_int),            allocatable  :: displs_recv_c(:), counts_recv_c(:)
      integer(c_int),            allocatable  :: displs_send_c(:)
      integer(i4),               pointer      :: mystart, mycount
      integer(c_int)                          :: mystart_c, mycount_c
      integer(i4)                             :: nrecv
      integer(i4)                             :: len_dim1, i, comm_f
      real(r8),      allocatable              :: temp_var(:,:)
      real(c_float), allocatable              :: temp_sp(:,:)
      character(c_char), allocatable          :: out_c(:), file_c(:)
      integer(c_int)                          :: all_comm_f_c, group_id_c, len_dim1_c, local_cols_c, open_mode_c
      logical                                 :: first_open
      character(kind=c_char), allocatable     :: varname_c(:)
      integer                                 :: vlen
      integer(i8)                             :: nlevp, nlevi8
      integer(i8)                             :: nvars
      character(kind=c_char), allocatable     :: ncfile_c(:)
      integer(i4)                             :: cmode, ncid, err
      integer(i4)                             :: i1, i2
      character(len=:), allocatable           :: ncfile_str
      integer(i4)                             :: ncfile_len
      integer(i4)                             :: fst_dim_id(6)           ! dim id for the 1st dim (location)
      integer(i4)                             :: sec_dim_id(3)           ! dim id for the 2nd dim (time)
      
      character(kind=c_char)                  :: path_buf(1024)
      character(len=:), allocatable           :: path_f
      character(len=:), allocatable           :: tmp
      integer(c_int)                          :: buf_len
      integer(i4)                             :: n, ierr
      
      call t_startf('2d_write') 
  ! 初始化
      nlevi8 = nlev
      nlevp = nlevi8 + 1
      comm => mesh%gcomm_write


      ! 构造 C 字符串（NUL 结尾）
      allocate(out_c(len_trim(outdir)+1))
      do i = 1, len_trim(outdir)
          out_c(i) = transfer(outdir(i:i), out_c(i))
      end do
      out_c(len_trim(outdir)+1) = c_null_char

      allocate(file_c(len_trim(filename)+1))
      do i = 1, len_trim(filename)
          file_c(i) = transfer(filename(i:i), file_c(i))
      end do
      file_c(len_trim(filename)+1) = c_null_char

      call t_startf('2d_nc_create')
      ptr => scalar_2d_list
      do while(associated(ptr))
          
          buf_len = int(size(path_buf), c_int)
          ! 先构建对应目录
          vlen = len_trim(ptr%varname)
          allocate(varname_c(vlen+1))
          do i = 1, vlen
              varname_c(i) = ptr%varname(i:i)
          end do
          varname_c(vlen+1) = c_null_char

          if (mpi_rank() == 0) then
            write(*,*) 'variable: ', trim(ptr%varname)
            call create_directory(out_c, file_c, varname_c, path_buf, buf_len)
            ! write(*,*) 'Created directory path: ', trim(transfer(path_buf, 'a'))
          end if 

          call MPI_Bcast(path_buf, buf_len, MPI_CHARACTER, 0, comm%all_comm, ierr)
          call MPI_Barrier(comm%all_comm, ierr)

          n = 0
          do while (n < size(path_buf) .and. path_buf(n+1) /= c_null_char)
            n = n + 1
          end do

          allocate(character(len=n) :: tmp)
          do i = 1, n
            tmp(i:i) = transfer(path_buf(i), 'a')
          end do
          call move_alloc(tmp, path_f)

          ! 输出提示信息
          if (mpi_rank() == 0) then
            write(*,*) 'Creating directory: ', trim(path_f)
          end if

          ! 构造每个变量的独立 nc 文件名
          ! ncfile_str = trim(outdir)//trim(filename)//'_'//trim(ptr%varname)//'.nc'
          ncfile_str = trim(path_f)//trim(ptr%varname)//'.dir/'//trim(filename)//'_'//trim(ptr%varname)//'.nc'
          ncfile_len = len_trim(ncfile_str)
          allocate(ncfile_c(ncfile_len+1))
          do i = 1, ncfile_len
              ncfile_c(i) = ncfile_str(i:i)
          end do
          ncfile_c(ncfile_len+1) = c_null_char
          ! 创建 NetCDF 文件
          cmode = IOR(NF_CLOBBER, NF_64BIT_DATA)
!          if (comm%color == 0) call wrap_create(comm%gcomm, ncfile_str, cmode, ncid)
!          ! 定义维度
!          if (comm%color == 0) then
!              call wrap_def_dim(ncid, 'nlev',    nlevi8,   fst_dim_id(1))
!              call wrap_def_dim(ncid, 'nlevp',   nlevi8+1, fst_dim_id(2))
!              call wrap_def_dim(ncid, 'dim12',   int8(12), fst_dim_id(3))
!              call wrap_def_dim(ncid, 'dim3',    int8(3),  fst_dim_id(4))
!              call wrap_def_dim(ncid, 'dim4',    int8(4),  fst_dim_id(5))
!              call wrap_def_dim(ncid, 'dim7',    int8(7),  fst_dim_id(6))
!              call wrap_def_dim(ncid, 'location_nt', mesh_nt, sec_dim_id(1))
!              call wrap_def_dim(ncid, 'location_ne', mesh_ne, sec_dim_id(2))
!              call wrap_def_dim(ncid, 'location_nv', mesh_nv, sec_dim_id(3))
!          end if
          if(ptr%var_level.eq.nlevi8)then
              i1 = 1
              if(ptr%var_pos.eq.1) i2 = 1
              if(ptr%var_pos.eq.6) i2 = 2
              if(ptr%var_pos.eq.0) i2 = 3
          end if
          if(ptr%var_level.eq.nlevp)then
              i1 = 2
              if(ptr%var_pos.eq.1) i2 = 1
              if(ptr%var_pos.eq.6) i2 = 2
              if(ptr%var_pos.eq.0) i2 = 3
          end if
          if(ptr%var_level.eq.12)then
              i1 = 3
              if(ptr%var_pos.eq.1) i2 = 1
              if(ptr%var_pos.eq.6) i2 = 2
              if(ptr%var_pos.eq.0) i2 = 3
          end if      
          if(ptr%var_level.eq.3)then
              i1 = 4
              if(ptr%var_pos.eq.1) i2 = 1
              if(ptr%var_pos.eq.6) i2 = 2
              if(ptr%var_pos.eq.0) i2 = 3
          end if      
          if(ptr%var_level.eq.4)then
              i1 = 5
              if(ptr%var_pos.eq.1) i2 = 1
              if(ptr%var_pos.eq.6) i2 = 2
              if(ptr%var_pos.eq.0) i2 = 3
          end if      
          if(ptr%var_level.eq.7)then
              i1 = 6
              if(ptr%var_pos.eq.1) i2 = 1
              if(ptr%var_pos.eq.6) i2 = 2
              if(ptr%var_pos.eq.0) i2 = 3
          end if
!          if (comm%color == 0) call def_var_2d_pnetcdf_sp(ptr%varname, fst_dim_id(i1), sec_dim_id(i2), ncid, ptr%var_idlist)
!          if (comm%color == 0) call put_att_2d_pnetcdf(ncid, ptr%var_idlist, ptr%varname, ptr%longname, ptr%units, ptr%var_pos)
!          ! 结束定义并关闭文件
!          if (comm%color == 0) call wrap_enddef(ncid)
!          if (comm%color == 0) call wrap_close(ncid)
          deallocate(ncfile_c)
          deallocate(varname_c)
          ptr => ptr%next
      end do
      call t_stopf('2d_nc_create')
      if (mpi_rank() == 0) then
        open(unit=11, file= trim(path_f)//'scalar_2d_list_info.txt', status='replace')
        ! write(*,*) 'txt Path: ', trim(path_f)
        ! open(unit=11, file= trim(path_f)//'scalar_2d_list_info.txt', status='unknown', position='append')

        write(11, '(A)') 'ncfile='//trim(filename)  ! 新增：写入nc文件名
        ! 统计变量个数
        nvars = 0
        ptr => scalar_2d_list
        do while (associated(ptr))
            nvars = nvars + 1
            ptr => ptr%next
        end do
        write(11, '(A, I0)') 'nvars=', nvars
        write(11, '(A)') 'varnames:'
        ptr => scalar_2d_list
        do i = 1, nvars
            write(11, '(A)') trim(ptr%varname)
            ptr => ptr%next
        end do
        close(11)
      end if

      ! ! 构造 C 字符串（NUL 结尾）
      ! allocate(out_c(len_trim(outdir)+1))
      ! do i = 1, len_trim(outdir)
      !     out_c(i) = transfer(outdir(i:i), out_c(i))
      ! end do
      ! out_c(len_trim(outdir)+1) = c_null_char

      ! allocate(file_c(len_trim(filename)+1))
      ! do i = 1, len_trim(filename)
      !     file_c(i) = transfer(filename(i:i), file_c(i))
      ! end do
      ! file_c(len_trim(filename)+1) = c_null_char

      ! 直接传递 Fortran 的 communicator（MPI_Fint）给 C
      all_comm_f_c = int(comm%all_comm, c_int)
      group_id_c = int(comm%color, c_int)

      first_open = .true.
      ptr => scalar_2d_list
      do while (associated(ptr))
          ! 选择本地列索引
          if (ptr%var_pos .eq. 0) then
            displs_send => comm%displs_sendv
            displs_recv => comm%displs_recvv
            counts_send => comm%counts_sendv
            counts_recv => comm%counts_recvv
            idx_recv    => comm%idx_recvv
            idx_local   => comm%idx_localv
            mystart     => comm%mystart_v
            mycount     => comm%mycount_v
          else if (ptr%var_pos .eq. 1) then 
            displs_send => comm%displs_sendt
            displs_recv => comm%displs_recvt
            counts_send => comm%counts_sendt
            counts_recv => comm%counts_recvt
            idx_recv    => comm%idx_recvt
            idx_local   => comm%idx_localt
            mystart     => comm%mystart_t
            mycount     => comm%mycount_t
          else if (ptr%var_pos .eq. 6) then 
            displs_send => comm%displs_sende
            displs_recv => comm%displs_recve
            counts_send => comm%counts_sende
            counts_recv => comm%counts_recve
            idx_recv    => comm%idx_recve
            idx_local   => comm%idx_locale
            mystart     => comm%mystart_e
            mycount     => comm%mycount_e
          end if

          ! 提取本地变量数据并转为单精度、连续内存
          len_dim1 = ptr%var_level
          allocate(temp_var(len_dim1, size(idx_local)))
          call find_output_var_2d(ptr%var, ptr%var_level, idx_local, temp_var)

          allocate(temp_sp(len_dim1, size(idx_local)))
          temp_sp = real(temp_var, kind=c_float)
          deallocate(temp_var)

          len_dim1_c = int(len_dim1, c_int)
          open_mode_c = merge(1_c_int, 0_c_int, first_open)
          local_cols_c = int(size(idx_local), c_int)

! 这里还需要提取每个进程的displs_send、displs_recv、counts_recv、idx_recv、mystart、mycount

          ! 转换 displs_send -> displs_send_c 
          if (associated(displs_send)) then
             nrecv = size(displs_send)
             allocate(displs_send_c(nrecv))
             do i = 1, nrecv
                displs_send_c(i) = int(displs_send(i), c_int)
             end do
          else
             allocate(displs_send_c(1))
             displs_send_c(1) = -1
             nrecv = 0
          end if

          ! 转换 displs_recv -> displs_recv_c
          if (associated(displs_recv)) then
             nrecv = size(displs_recv)
             allocate(displs_recv_c(nrecv))
             do i = 1, nrecv
                displs_recv_c(i) = int(displs_recv(i), c_int)
             end do
          else
             allocate(displs_recv_c(1))
             displs_recv_c(1) = -1
             nrecv = 0
          end if

          ! 转换 counts_recv -> counts_recv_c
          if (associated(counts_recv)) then
             if (nrecv == 0) nrecv = size(counts_recv)
             allocate(counts_recv_c(nrecv))
             do i = 1, nrecv
                counts_recv_c(i) = int(counts_recv(i), c_int)
             end do
          else
             allocate(counts_recv_c(1))
             counts_recv_c(1) = -1
             nrecv = 0
          end if

          ! 转换 idx_recv -> idx_recv_c（长度为 size(idx_recv)）
          if (associated(idx_recv)) then
             nrecv = size(idx_recv)
             allocate(idx_recv_c(nrecv))
             do i = 1, nrecv
                idx_recv_c(i) = int(idx_recv(i), c_int)
             end do
          else
             allocate(idx_recv_c(1))
             idx_recv_c(1) = -1
             nrecv = 0
          end if

          if (associated(mystart)) then
              mystart_c = int(mystart, c_int)
          else
              mystart_c = 0_c_int
          end if

          if (associated(mycount)) then
              mycount_c = int(mycount, c_int)
          else
              mycount_c = 0_c_int
          end if

          ! 构造 C 字符串
          vlen = len_trim(ptr%varname)
          allocate(varname_c(vlen+1))
          do i = 1, vlen
              varname_c(i) = ptr%varname(i:i)
          end do
          varname_c(vlen+1) = c_null_char

          ! write(*,*) ' idx_recv_size=', size(idx_recv), ' mycount_c=', mycount_c
          ! call wrap_output_2d_gather_sp_c(temp_sp, len_dim1_c, local_cols_c, &
          !                                 group_id_c, all_comm_f_c, &
          !                                 out_c, file_c, &
          !                                 open_mode_c, &
          !                                 mycount_c, mystart_c, &
          !                                 displs_send_c, counts_recv_c, &
          !                                 displs_recv_c, idx_recv_c, varname_c)
          call wrap_output_2d_gather_sp_c(temp_sp, len_dim1_c, local_cols_c, &
                                          group_id_c, all_comm_f_c, &
                                          path_buf, file_c, &
                                          open_mode_c, &
                                          mycount_c, mystart_c, &
                                          displs_send_c, counts_recv_c, &
                                          displs_recv_c, idx_recv_c, varname_c)
          deallocate(varname_c)
          deallocate(temp_sp)
          deallocate(displs_send_c)
          deallocate(displs_recv_c)
          deallocate(counts_recv_c)
          deallocate(idx_recv_c)
          ! first_open = .false.
          ptr => ptr%next
      end do
      if (allocated(path_f)) deallocate(path_f)
      call t_stopf('2d_write') 
   end subroutine wrap_output_2d_gather_C
!================================================
! 6th: bin -> nc
!================================================
   subroutine wrap_output_after_process_2d(mesh, outdir, filename)
    implicit none
! io
      type(global_domain), intent(in), target :: mesh
      character(*),       intent(in)          :: outdir
      character(*),       intent(in)          :: filename
! local
      type(group_comm),          pointer      :: comm
      type(array_2d_index_list), pointer      :: ptr
      integer(i8)                             :: nlevi8, nlevp
      integer(i4)                             :: fst_dim_id(6)           ! dim id for the 1st dim (location)
      integer(i4)                             :: sec_dim_id(3)           ! dim id for the 2nd dim (time)
      integer(i4)                             :: cmode, ncid
      integer(i4)                             :: color, group_count_c, all_comm_f_c
      integer(i4)                             :: i1, i2, err
      character(kind=c_char), allocatable     :: ncfile_c(:)
      character(len=:), allocatable           :: ncfile_str
      integer(i4)                             :: ncfile_len, i
      integer(i4)                             :: nvars, maxlen, curlen
      character(len=:), allocatable           :: varnames(:)
      character(kind=c_char), allocatable     :: varnames_c(:)
      character(len=256), allocatable         :: varnames_tmp(:)
      integer(i4)                             :: total_len, pos, j
      integer                                 :: ios
      character(len=256)                      :: line

! 初始化
    nlevi8 = nlev
    nlevp = nlevi8 + 1
    comm => mesh%gcomm_write

! 获取通信信息
    all_comm_f_c = int(comm%all_comm, c_int)
    color = int(comm%color, c_int)

! ! 调用 C 端程序将二进制文件转换为 NetCDF 文件
! ! 构造 nc 文件名 C 字符串
!     ncfile_str = trim(outdir)//trim(filename)
!     ncfile_len = len_trim(ncfile_str)
!     allocate(ncfile_c(ncfile_len + 1))
!     do i = 1, ncfile_len
!         ncfile_c(i) = ncfile_str(i:i)
!     end do
!     ncfile_c(ncfile_len + 1) = c_null_char

! ! 统计变量个数
!     nvars = 0
!     maxlen = 0
!     ptr => scalar_2d_list
!     do while (associated(ptr))
!         nvars = nvars + 1
!         curlen = len_trim(ptr%varname)
!         if (curlen > maxlen) maxlen = curlen
!         ptr => ptr%next
!     end do
! ! 分配变量名数组
!     allocate(character(len=maxlen) :: varnames(nvars))
! ! 再次遍历链表，收集变量名
!     ptr => scalar_2d_list
!     i = 1
!     do while (associated(ptr))
!         varnames(i) = trim(ptr%varname)
!         i = i + 1
!         ptr => ptr%next
!     end do
    open(unit=20, file=trim(outdir)//'scalar_2d_list_info.txt', status='old', action='read')
    
    ! 读取 nc 文件名
    read(20, '(A)', iostat=ios) line
    ! 离线还原不需要trim(outdir)，直接把脚本和results放在同一目录下即可
    ! 这里是在线还原的测试代码
    ncfile_str = trim(outdir)//adjustl(trim(line(8:)))   ! 跳过'ncfile='，取第8位起
    ncfile_len = len_trim(ncfile_str)
    allocate(ncfile_c(ncfile_len + 1))
    do i = 1, ncfile_len
        ncfile_c(i) = ncfile_str(i:i)
    end do
    ncfile_c(ncfile_len + 1) = c_null_char

    ! 读取变量个数
    read(20, '(A)', iostat=ios) line
    read(line(7:), *) nvars   ! 跳过'nvars='
    ! 跳过'varnames:'这一行
    read(20, '(A)', iostat=ios) line

    allocate(character(len=256) :: varnames_tmp(nvars))
    maxlen = 0
    do i = 1, nvars
        read(20, '(A)', iostat=ios) line
        varnames_tmp(i) = adjustl(trim(line))
        if (len_trim(varnames_tmp(i)) > maxlen) maxlen = len_trim(varnames_tmp(i))
    end do
    close(20)

    ! 分配最终的 varnames，并拷贝内容
    allocate(character(len=maxlen) :: varnames(nvars))
    do i = 1, nvars
        varnames(i) = varnames_tmp(i)(1:maxlen)
    end do
    deallocate(varnames_tmp)

    write(*,*) ' ncfile to process:', trim(ncfile_str)
    write(*,*) ' 2D vars to process nvars=', nvars
    write(*,*) ' Varnames:', varnames 

! 统计所有变量名总长度（每个加1用于c_null_char）
    total_len = 0
    do i = 1, nvars
        total_len = total_len + len_trim(varnames(i)) + 1
    end do
! 分配C字符串数组
    allocate(varnames_c(total_len))
    pos = 1
    do i = 1, nvars
        do j = 1, len_trim(varnames(i))
            varnames_c(pos) = varnames(i)(j:j)
            pos = pos + 1
        end do
        varnames_c(pos) = c_null_char
        pos = pos + 1
    end do

    call bin_to_nc_group_sp(ncfile_c, all_comm_f_c, color, nvars, varnames_c)
    ! 让 C 端完成多进程数据写入
    ! 只传 nc 文件名给 C 端
    deallocate(varnames)
    deallocate(varnames_c)
    deallocate(ncfile_c)

   end subroutine wrap_output_after_process_2d

!
! 清理工作区
!
   subroutine wrap_output_clean_2d()

     type(array_2d_index_list), pointer :: field_list, tmpf
    
     field_list => scalar_2d_list
     do while(associated(field_list))
       tmpf => field_list%next
       deallocate(field_list)
       field_list => tmpf
     end do
     scalar_2d_list => null()

     return
   end subroutine wrap_output_clean_2d
!
! read routines
!

   subroutine wrap_read_2d_group(comm, outdir, filename, varname, dim2_len, option, varout)
! io
     type(group_comm), intent(in), target :: comm
     character*(*),    intent(in)         :: outdir
     character*(*),    intent(in)         :: filename
     character*(*),    intent(in)         :: varname
     integer(i8),      intent(in)         :: dim2_len
     integer                              :: option
     real(r8),         intent(inout)      :: varout(:,:)
! local
     integer(i4),      parameter          :: omode = 0
     integer(i4)                          :: varid, ncid
     integer(i4)                          :: dim2, act_get_num, mpi_newtype
     integer(i4),      pointer            :: displs_send(:), counts_send(:)
     integer(i4),      pointer            :: displs_recv(:), counts_recv(:)
     integer(i4),      pointer            :: idx_send(:), idx_local(:)
     integer(i4),      pointer            :: mystart, mycount
     integer(i8)                          :: start(2)
     integer(i8)                          :: count(2)
     integer(i4)                          :: req(1), st(1)
     real(r8),         allocatable        :: buff_var(:), send_var(:), temp_var(:)
     integer                              :: i, err, j, xtype
     integer                              :: max_diml, jres, jstart
     real(4),          allocatable        :: local_data_r4(:,:)
     real(8),          allocatable        :: local_data_r8(:,:)

     dim2 = dim2_len
     if (option .eq. 0) then
       displs_send => comm%displs_sendv
       displs_recv => comm%displs_recvv
       counts_send => comm%counts_sendv
       counts_recv => comm%counts_recvv
       idx_send    => comm%idx_sendv
       idx_local   => comm%idx_localv
       mystart     => comm%mystart_v
       mycount     => comm%mycount_v
     else if (option .eq. 1) then
       displs_send => comm%displs_sendt
       displs_recv => comm%displs_recvt
       counts_send => comm%counts_sendt
       counts_recv => comm%counts_recvt
       idx_send    => comm%idx_sendt
       idx_local   => comm%idx_localt
       mystart     => comm%mystart_t
       mycount     => comm%mycount_t
     else if (option .eq. 6) then
       displs_send => comm%displs_sende
       displs_recv => comm%displs_recve
       counts_send => comm%counts_sende
       counts_recv => comm%counts_recve
       idx_send    => comm%idx_sende
       idx_local   => comm%idx_locale
       mystart     => comm%mystart_e
       mycount     => comm%mycount_e
     end if

     if (comm%color == 0 ) then
       call wrap_open(MPI_COMM_SELF,trim(outdir)//trim(filename),omode,ncid)

       call wrap_inq_varid(ncid, trim(varname), varid)
       call wrap_inq_vartype(ncid, varid, xtype)
       ! get the actual number of inputs (due to the limitation of pnetcdf)
       if (xtype .eq. NF_DOUBLE) then
         max_diml = MAXLEN/dim2
       elseif (xtype .eq. NF_FLOAT) then
         max_diml = 2*MAXLEN/dim2
       end if
       act_get_num = ceiling(mycount/real(max_diml, 8))

       allocate(buff_var(mycount*dim2))
       start(2)       = 1
       count(2)       = dim2_len
       ! may need to get multiple times (each time < 2GB)
       jres   = mod(mycount, max_diml)
       do j = 1, act_get_num
         if (j .lt. act_get_num) then
             count(1) = max_diml
         else
             count(1) = jres
         end if
         start(1) = mystart + (j-1)*max_diml
         jstart   = 1 + (j-1)*max_diml
         if (xtype .eq. NF_DOUBLE) then
           allocate(local_data_r8(count(1), dim2))
           call wrap_get_vara_realx_req(ncid, varid, start, count, local_data_r8, req)
           err = nfmpi_wait_all(ncid, 1, req, st)
           do i = 1, count(1)
             buff_var((jstart+i-2)*dim2+1:(jstart+i-1)*dim2) = local_data_r8(i, 1:dim2)
           end do
           deallocate(local_data_r8)
         elseif (xtype .eq. NF_FLOAT) then
           allocate(local_data_r4(count(1), dim2))
           call wrap_get_vara_real_req(ncid, varid, start, count, local_data_r4, req)
           err = nfmpi_wait_all(ncid, 1, req, st)
           do i = 1, count(1)
             buff_var((jstart+i-2)*dim2+1:(jstart+i-1)*dim2) = local_data_r4(i, 1:dim2)
           end do
           deallocate(local_data_r4)
         end if
       end do
       call wrap_close(ncid)
       allocate(send_var(size(idx_send)*dim2))
       do i = 1, size(idx_send)
         send_var((i-1)*dim2+1:i*dim2) = buff_var((idx_send(i)-1)*dim2+1:idx_send(i)*dim2)
       end do
       deallocate(buff_var)
     else
       allocate(send_var(1))
     end if

     allocate(temp_var(size(idx_local)*dim2))

     if (r8 .eq. 8) then
       call MPI_Type_contiguous(dim2, MPI_REAL8, mpi_newtype, err)
     elseif (r8 .eq. 4) then
       call MPI_Type_contiguous(dim2, MPI_REAL, mpi_newtype, err)
     end if
     call MPI_Type_commit(mpi_newtype, err)
     call MPI_Alltoallv(send_var, counts_send, displs_send, mpi_newtype, temp_var, &
                        counts_recv, displs_recv, mpi_newtype, comm%all_comm, err)

     do i = 1, size(varout,1)
        varout(i, 1:dim2) = temp_var((idx_local(i)-1)*dim2+1:idx_local(i)*dim2)
     end do

     deallocate(send_var)
     deallocate(temp_var)
     if(mpi_rank().eq.0) write(6,*) "Sucessfully read ",trim(varname)," from ", trim(outdir)//trim(filename)
     return
   end subroutine wrap_read_2d_group

   subroutine wrap_read_2d_group_rst(comm, outdir, filename, varname, dim2_len, option, varout)
! io
     type(group_comm), intent(in), target :: comm
     character*(*),    intent(in)         :: outdir
     character*(*),    intent(in)         :: filename
     character*(*),    intent(in)         :: varname
     integer(i8),      intent(in)         :: dim2_len
     integer                              :: option
     real(r8),         intent(inout)      :: varout(:,:)
! local
     integer(i4),      parameter          :: omode = 0
     integer(i4)                          :: varid, ncid
     integer(i4)                          :: dim2, act_get_num, mpi_newtype
     integer(i4),      pointer            :: displs_send(:), counts_send(:)
     integer(i4),      pointer            :: displs_recv(:), counts_recv(:)
     integer(i4),      pointer            :: idx_send(:), idx_local(:)
     integer(i4),      pointer            :: mystart, mycount
     integer(i8)                          :: start(2)
     integer(i8)                          :: count(2)
     integer(i4)                          :: req(1), st(1)
     real(r8),         allocatable        :: buff_var(:), send_var(:), temp_var(:)
     integer                              :: i, err, j, xtype
     integer                              :: max_diml, jres, jstart
     real(4),          allocatable        :: local_data_r4(:,:)
     real(8),          allocatable        :: local_data_r8(:,:)

     dim2 = dim2_len
     if (option .eq. 0) then
       displs_send => comm%displs_sendv
       displs_recv => comm%displs_recvv
       counts_send => comm%counts_sendv
       counts_recv => comm%counts_recvv
       idx_send    => comm%idx_sendv
       idx_local   => comm%idx_localv
       mystart     => comm%mystart_v
       mycount     => comm%mycount_v
     else if (option .eq. 1) then
       displs_send => comm%displs_sendt
       displs_recv => comm%displs_recvt
       counts_send => comm%counts_sendt
       counts_recv => comm%counts_recvt
       idx_send    => comm%idx_sendt
       idx_local   => comm%idx_localt
       mystart     => comm%mystart_t
       mycount     => comm%mycount_t
     else if (option .eq. 6) then
       displs_send => comm%displs_sende
       displs_recv => comm%displs_recve
       counts_send => comm%counts_sende
       counts_recv => comm%counts_recve
       idx_send    => comm%idx_sende
       idx_local   => comm%idx_locale
       mystart     => comm%mystart_e
       mycount     => comm%mycount_e
     end if

     if (comm%color == 0 ) then
       call wrap_open(MPI_COMM_SELF,trim(outdir)//trim(filename),omode,ncid)

       call wrap_inq_varid(ncid, trim(varname), varid)
       call wrap_inq_vartype(ncid, varid, xtype)
       ! get the actual number of inputs (due to the limitation of pnetcdf)
       if (xtype .eq. NF_DOUBLE) then
         max_diml = MAXLEN/dim2
       elseif (xtype .eq. NF_FLOAT) then
         max_diml = 2*MAXLEN/dim2
       end if
       act_get_num = ceiling(mycount/real(max_diml, 8))

       allocate(buff_var(mycount*dim2))
       start(1) = 1
       count(1) = dim2_len
       ! may need to get multiple times (each time < 2GB)
       jres   = mod(mycount, max_diml)
       do j = 1, act_get_num
         if (j .lt. act_get_num) then
             count(2) = max_diml
         else
             count(2) = jres
         end if
         start(2) = mystart + (j-1)*max_diml
         jstart   = 1 + (j-1)*max_diml
         if (xtype .eq. NF_DOUBLE) then
           allocate(local_data_r8(dim2, count(2)))
           call wrap_get_vara_realx_req(ncid, varid, start, count, local_data_r8, req)
           err = nfmpi_wait_all(ncid, 1, req, st)
           do i = 1, count(2)
             buff_var((jstart+i-2)*dim2+1:(jstart+i-1)*dim2) = local_data_r8(1:dim2, i)
           end do
           deallocate(local_data_r8)
         elseif (xtype .eq. NF_FLOAT) then
           allocate(local_data_r4(dim2, count(2)))
           call wrap_get_vara_real_req(ncid, varid, start, count, local_data_r4, req)
           err = nfmpi_wait_all(ncid, 1, req, st)
           do i = 1, count(2)
             buff_var((jstart+i-2)*dim2+1:(jstart+i-1)*dim2) = local_data_r4(1:dim2, i)
           end do
           deallocate(local_data_r4)
         end if
       end do
       call wrap_close(ncid)
       allocate(send_var(size(idx_send)*dim2))
       do i = 1, size(idx_send)
         send_var((i-1)*dim2+1:i*dim2) = buff_var((idx_send(i)-1)*dim2+1:idx_send(i)*dim2)
       end do
       deallocate(buff_var)
     else
       allocate(send_var(1))
     end if

     allocate(temp_var(size(idx_local)*dim2))

     if (r8 .eq. 8) then
       call MPI_Type_contiguous(dim2, MPI_REAL8, mpi_newtype, err)
     elseif (r8 .eq. 4) then
       call MPI_Type_contiguous(dim2, MPI_REAL, mpi_newtype, err)
     end if
     call MPI_Type_commit(mpi_newtype, err)
     call MPI_Alltoallv(send_var, counts_send, displs_send, mpi_newtype, temp_var, &
                        counts_recv, displs_recv, mpi_newtype, comm%all_comm, err)

     do i = 1, size(varout,2)
        varout(1:dim2, i) = temp_var((idx_local(i)-1)*dim2+1:idx_local(i)*dim2)
     end do

     deallocate(send_var)
     deallocate(temp_var)
     if(mpi_rank().eq.0) write(6,*) "Sucessfully read ",trim(varname)," from ", trim(outdir)//trim(filename)
     return
   end subroutine wrap_read_2d_group_rst

   subroutine wrap_read_2d(comm,outdir,filename,varname,dim1_len,dim2_len,varout)
! io
     integer,               intent(in)    :: comm
     character*(*),         intent(in)    :: outdir
     character*(*),         intent(in)    :: filename
     character*(*),         intent(in)    :: varname
     integer(i8),           intent(in)    :: dim1_len
     integer(i8),           intent(in)    :: dim2_len
     type(scalar_2d_field), intent(inout) :: varout
! local
     integer(i4),           parameter     :: omode = 0
     integer(i4)                          :: varid
     integer(i4)                          :: ncid
     integer(i4)                          :: xtype
     integer(i8)                          :: start(2)
     integer(i8)                          :: count(2)
     integer                              :: max_diml, act_get_num, j, jres
! temp arrays
     real(8),               allocatable   :: temp_dp(:,:)
     real(4),               allocatable   :: temp_sp(:,:)

     call wrap_open(comm,trim(outdir)//trim(filename),omode,ncid)
     call wrap_inq_varid(ncid, trim(varname), varid)
     call wrap_inq_vartype(ncid, varid, xtype)
!================================================
!              READ vars at edge
!================================================
     start(2)       = 1
     count(2)       = dim2_len
     if (xtype .eq. NF_DOUBLE) then
       max_diml = MAXLEN/dim2_len
     elseif (xtype .eq. NF_FLOAT) then
       max_diml = 2*MAXLEN/dim2_len
     end if
     act_get_num = ceiling(dim1_len/real(max_diml, 8))
     jres = mod(dim1_len, max_diml)

     do j = 1, act_get_num
       if (j .lt. act_get_num) then
           count(1) = max_diml
       else
           count(1) = jres
       end if
       start(1) = 1 + (j-1)*max_diml
       if (xtype .eq. NF_DOUBLE) then
         allocate(temp_dp(count(1), dim2_len))
         call wrap_get_vara_realx(ncid, varid, start, count, temp_dp)
         varout%f(start(1):(start(1)+count(1)-1), 1:dim2_len) = temp_dp
         deallocate(temp_dp)
       elseif (xtype .eq. NF_FLOAT) then
         allocate(temp_sp(count(1), dim2_len))
         call wrap_get_vara_real(ncid, varid, start, count, temp_sp)
         varout%f(start(1):(start(1)+count(1)-1), 1:dim2_len) = temp_sp
         deallocate(temp_sp)
       end if
     end do
     call wrap_close(ncid)

     return

   end subroutine wrap_read_2d
   
   subroutine wrap_bcast_2d(comm, irank, dim1_len, dim2_len, varout)
! io
     integer,               intent(in)    :: comm
     integer,               intent(in)    :: irank
     integer(i8),           intent(in)    :: dim1_len
     integer(i8),           intent(in)    :: dim2_len
     type(scalar_2d_field), intent(inout) :: varout

! local
     integer                              :: err
     integer                              :: var_size
     integer                              :: root_rank
     real(r8), allocatable                :: buf(:,:)
     
     allocate(buf(dim1_len, dim2_len))
     
     buf = varout%f
     var_size = dim1_len*dim2_len
     root_rank = irank
     if (r8 .eq. 8) then
       call MPI_Bcast(buf, var_size, MPI_REAL8, root_rank, comm, err)
     elseif (r8 .eq. 4) then
       call MPI_Bcast(buf, var_size, MPI_REAL, root_rank, comm, err)
     end if
     varout%f = buf
     deallocate(buf)
   end subroutine wrap_bcast_2d

!---------------------------------------------------------------------------
!
!                         Private routines below
!
!---------------------------------------------------------------------------

   subroutine def_var_2d_pnetcdf(varname, fst_dim_id, sec_dim_id, ncid, var_idlist)
! io
     character(len=*), intent(in)    :: varname
     integer(i4),      intent(in)    :: fst_dim_id
     integer(i4),      intent(in)    :: sec_dim_id
     integer(i4),      intent(in)    :: ncid
     integer(i4),      intent(inout) :: var_idlist

! local
     integer(i4)                     :: vdims(2)
     integer(i4)                     :: len_two = 2
    
     vdims(1) = fst_dim_id
     vdims(2) = sec_dim_id
     if (r8 .eq. 8) then
       call wrap_def_var(ncid,trim(varname),NF_DOUBLE, len_two, vdims(1:2), var_idlist)
     elseif (r8 .eq. 4) then
       call wrap_def_var(ncid,trim(varname),NF_FLOAT, len_two, vdims(1:2), var_idlist)
     end if

     return
   end subroutine def_var_2d_pnetcdf

   subroutine def_var_2d_pnetcdf_sp(varname, fst_dim_id, sec_dim_id, ncid, var_idlist)
! io
     character(len=*), intent(in)    :: varname
     integer(i4),      intent(in)    :: fst_dim_id
     integer(i4),      intent(in)    :: sec_dim_id
     integer(i4),      intent(in)    :: ncid
     integer(i4),      intent(inout) :: var_idlist

! local
     integer(i4)                     :: vdims(2)
     integer(i4)                     :: len_two = 2
    
     vdims(1) = fst_dim_id
     vdims(2) = sec_dim_id
     call wrap_def_var(ncid,trim(varname),NF_FLOAT, len_two, vdims(1:2), var_idlist)

     return
   end subroutine def_var_2d_pnetcdf_sp

   subroutine put_att_2d_pnetcdf(ncid, var_idlist, varname, longname,units, var_pos)
! io
     integer(i4),              intent(in)   :: ncid
     integer(i4),              intent(in)   :: var_idlist
     character(len=*),         intent(in)   :: varname
     character(len=*),         intent(in)   :: longname
     character(len=*),         intent(in)   :: units
     integer(i4),              intent(in)   :: var_pos

     call wrap_put_att_text(ncid,var_idlist,'long_name', trim(longname))

     if(trim(varname).eq."lon_nt".or.&
        trim(varname).eq."lon_ne".or.&
        trim(varname).eq."lon_nv")then
        call wrap_put_att_text(ncid,var_idlist,'units', "degrees_east")
     else if (trim(varname).eq."lat_nt".or.&
              trim(varname).eq."lat_ne".or.&
              trim(varname).eq."lat_nv")then
              call wrap_put_att_text(ncid,var_idlist,'units', "degrees_north")
     else if(var_pos.eq.1)then
        call wrap_put_att_text(ncid,var_idlist,'coordinates', "lon_nt lat_nt")
        call wrap_put_att_text(ncid,var_idlist,'units', trim(units))
     else if(var_pos.eq.6)then
        call wrap_put_att_text(ncid,var_idlist,'coordinates', "lon_ne lat_ne")
        call wrap_put_att_text(ncid,var_idlist,'units', trim(units))
     else if(var_pos.eq.0)then
        call wrap_put_att_text(ncid,var_idlist,'coordinates', "lon_nv lat_nv")
        call wrap_put_att_text(ncid,var_idlist,'units', trim(units))
     end if 

     return

   end subroutine put_att_2d_pnetcdf

   subroutine find_output_var_2d(local_var,level_var,output_local_index,output_local_var)
! io
     real(r8),    pointer, intent(in)    :: local_var(:,:)
     integer(i8), intent(in)             :: level_var
     integer(i4), intent(in)             :: output_local_index(:)
     real(r8),    intent(out)            :: output_local_var(:,:)
  
! local
     integer(i8) :: idex,idey
     
     do idex = 1,size(output_local_index)
        do idey = 1,level_var 
           output_local_var(idey,idex) = local_var(idey,output_local_index(idex))
        end do
     end do

   end subroutine find_output_var_2d

 end module grist_fileio_list_2d_module_par
