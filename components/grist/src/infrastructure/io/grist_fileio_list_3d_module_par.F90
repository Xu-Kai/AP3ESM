
!-----------------------------------------------------------
! Created on 2018
! Authors/Contact: Liu Zhuang
! Version 1.0
! Description: File Par I/O module tailored for GRIST. It 
! deals with the general conditions, must be called by an 
! explicit suit of statements.
!            How to use:
!             1) wrap_output_init_3d
!             2) wrap_add_field_3d
!             3) wrap_output_3d
!             4) wrap_output_clean_3d
!-----------------------------------------------------------

 module grist_fileio_list_3d_module_par
   use grist_lib
   use grist_constants,      only: i4, i8, r4, r8
   use grist_domain_types  , only: global_domain, block_structure, group_comm
   use grist_data_types    , only: scalar_3d_field, data3d_temp_dp, data3d_temp_sp
   use grist_handle_error,   only: endrun
   use grist_nml_module,     only: nlev, ntracer
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
                                   wrap_put_vara_real,        &
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

   public :: wrap_output_init_3d    ,&
             wrap_add_field_3d      ,&
             wrap_output_3d_group   ,&
             wrap_output_3d_group_sp,&
             wrap_output_3d_gather_C,&
             wrap_output_clean_3d   ,&
             wrap_read_3d_group_rst

   interface
     subroutine wrap_output_3d_gather_sp_c(temp, len_dim1, len_dim2, local_cols, group_id, all_comm_f, outdir, filename, open_mode, &
        my_count, my_start, displs_send, counts_recv, displs_recv, idx_recv, varname) bind(C, name="wrap_output_3d_gather_sp_c")
        use, intrinsic                                                  :: iso_c_binding, only: c_int, c_float, c_char
        real(c_float), intent(in)                                       :: temp(*)
        integer(c_int), value                                           :: len_dim1, len_dim2, local_cols, group_id, all_comm_f, open_mode
        character(kind=c_char), dimension(*), intent(in)                :: outdir
        character(kind=c_char), dimension(*), intent(in)                :: filename
        integer(c_int), value                                           :: my_count, my_start
        integer(c_int), dimension(*), intent(in)                        :: displs_send
        integer(c_int), dimension(*), intent(in)                        :: counts_recv
        integer(c_int), dimension(*), intent(in)                        :: displs_recv
        integer(c_int), dimension(*), intent(in)                        :: idx_recv
        character(kind=c_char), dimension(*), intent(in)                :: varname
     end subroutine wrap_output_3d_gather_sp_c


     subroutine create_directory(outdir, filename, varname, outpath, outlen) bind(C, name="create_directory")
          use, intrinsic                                                  :: iso_c_binding, only: c_int, c_char
          character(kind=c_char), dimension(*), intent(in)                :: outdir
          character(kind=c_char), dimension(*), intent(in)                :: filename
          character(kind=c_char), dimension(*), intent(in)                :: varname
          character(kind=c_char), dimension(*), intent(out)               :: outpath
          integer(c_int), value                                           :: outlen
     end subroutine create_directory

    end interface

   type(array_3d_index_list), pointer   :: scalar_3d_list   ! any 3d shape
   integer(i4)                          :: number_of_scalar
   integer(i4)                          :: mesh_maxvnb
   integer(i8)                          :: mesh_nt
   integer(i8)                          :: mesh_ne
   integer(i8)                          :: mesh_nv
   integer(i4), PARAMETER               :: MAXLEN = 268000000 !required by pnetcdf (< 2GB/8B)
   real(8)                              :: default_8 = 1D38
   real(4)                              :: default_4 = 1e38

 contains

   subroutine wrap_output_init_3d(mesh)

      type(global_domain), intent(in) :: mesh

      mesh_maxvnb                 = mesh%maxvnb
      mesh_nt                     = mesh%nt_all
      mesh_ne                     = mesh%ne_all
      mesh_nv                     = mesh%nv_all
      number_of_scalar            = 0

      return
   end subroutine wrap_output_init_3d

   subroutine wrap_add_field_3d(scalar,scalar_name,longname,units)

      type(scalar_3d_field), intent(in)           :: scalar        ! 3d array
      character*(*)        , intent(in)           :: scalar_name
      character*(*)        , intent(in), optional :: longname
      character*(*)        , intent(in), optional :: units

      if(present(longname).and.present(units))then
         call insert_var_3d(scalar_3d_list, scalar%f, trim(scalar_name),trim(longname),trim(units), number_of_scalar, scalar%pos)
      else
         call insert_var_3d(scalar_3d_list, scalar%f, trim(scalar_name),'','', number_of_scalar, scalar%pos)
      end if

      return
   end subroutine wrap_add_field_3d

   subroutine wrap_output_3d_group(mesh, outdir, filename, dim3_len)

   implicit none
!io
     type(global_domain), intent(in), target   :: mesh
     character*(*),       intent(in)           :: outdir
     character*(*),       intent(in)           :: filename
     integer(i4),         intent(in), optional :: dim3_len 
! local
     type(group_comm),          pointer        :: comm
     type(block_structure),     pointer        :: local_block
     type(array_3d_index_list), pointer        :: ptr
     integer(i4)                               :: cmode, ncid, err
     integer(i4)                               :: fst_dim_id(1)           ! dim id for the 1st dim
     integer(i4)                               :: sec_dim_id(3)           ! dim id for the 2nd dim
     integer(i4)                               :: thr_dim_id(3)           ! dim id for the 3rd dim
     integer(i4)                               :: len_dim1, mpi_newtype
     integer(i4)                               :: i0, i1, i2 ,i, j, k, lidx
     integer(i8)                               :: nlevp, nlevi8, ntraceri8
     real(r8),                  allocatable    :: buff_var(:), temp_var(:,:,:)
     real(8),                   allocatable    :: buff_put_dp(:,:,:)
     real(4),                   allocatable    :: buff_put_sp(:,:,:)
     integer(i4),               pointer        :: displs_send(:), counts_send(:)
     integer(i4),               pointer        :: displs_recv(:), counts_recv(:)
     integer(i4),               pointer        :: idx_recv(:), idx_local(:)
     integer(i4),               pointer        :: mystart, mycount
     integer(i8)                               :: start(3)
     integer(i8)                               :: count(3)
     integer(i4)                               :: req(1), st(1), act_put_num
     integer(i4)                               :: max_diml, jres, jstart, ndim3

     nlevi8 = nlev
     nlevp= nlevi8+1
     ntraceri8 = ntracer ! default is ntracer
     ndim3     = ntracer
     if(present(dim3_len))then
        ntraceri8 = dim3_len
        ndim3     = dim3_len
     end if

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
     if(comm%color==0) then
! possible dims for 1st
       call wrap_def_dim(ncid, 'ntracer'    , ntraceri8   , fst_dim_id(1)) ! gcm
! possible dims for 2nd
       call wrap_def_dim(ncid, 'nlev'       , nlevi8      , sec_dim_id(1)) ! gcm
       call wrap_def_dim(ncid, 'nlevp'      , nlevi8+1    , sec_dim_id(2)) ! gcm
       call wrap_def_dim(ncid, 'dim12'      , int8(12)    , sec_dim_id(3)) ! gcm    
! possible dims for 3rd
       call wrap_def_dim(ncid, 'location_nt', mesh_nt     , thr_dim_id(1)) ! gcm
       call wrap_def_dim(ncid, 'location_ne', mesh_ne     , thr_dim_id(2)) ! gcm
       call wrap_def_dim(ncid, 'location_nv', mesh_nv     , thr_dim_id(3)) ! gcm
     endif
!========================================================
! 2nd: define var, need lon_id, lat_id, var_idlist
!========================================================
     ptr => scalar_3d_list
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
       if(comm%color==0) call def_var_3d_pnetcdf(ptr%varname, fst_dim_id(1), sec_dim_id(i1), thr_dim_id(i2), ncid, ptr%var_idlist)
       ptr => ptr%next
     end do

!================================================
! 3rd: put att
!================================================    
     ptr => scalar_3d_list
     do while (associated(ptr))
       if(comm%color==0) call put_att_3d_pnetcdf( ncid, ptr%var_idlist ,ptr%varname,ptr%longname,ptr%units,ptr%var_pos)
       ptr => ptr%next
     end do 
     if(comm%color==0) call wrap_enddef(ncid)

!================================================
! 4th, put var
!================================================

     ptr => scalar_3d_list
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
       allocate(temp_var(ndim3, len_dim1, size(idx_local)))
       call find_output_var_3d(ptr%var, ptr%var_level, idx_local, temp_var, ndim3)

       if (comm%color == 0) then
         allocate(buff_var(1:mycount*len_dim1*ndim3))
       else
         allocate(buff_var(1))
       endif

       if (r8 .eq. 8) then
         call MPI_Type_contiguous(len_dim1*ndim3, MPI_REAL8, mpi_newtype, err)
       elseif (r8 .eq. 4) then
         call MPI_Type_contiguous(len_dim1*ndim3, MPI_REAL, mpi_newtype, err)
       end if
       call MPI_Type_commit(mpi_newtype, err)
       call MPI_Alltoallv(temp_var, counts_send, displs_send, mpi_newtype, buff_var, &
                          counts_recv, displs_recv, mpi_newtype, comm%all_comm, err)
       deallocate(temp_var)

       if (comm%color == 0) then
         ! get the actual number of outputs (due to the limitation of pnetcdf)
         if (r8 .eq. 8) then
           max_diml = MAXLEN/(ndim3*len_dim1)
         elseif (r8 .eq. 4) then
           max_diml = 2*MAXLEN/(ndim3*len_dim1)
         end if
         act_put_num = ceiling(mycount/real(max_diml, r8))

         if (r8 .eq. 8) then
           allocate(buff_put_dp(ndim3, len_dim1, mycount), source=default_8)
           do i = 1, size(idx_recv)
              lidx = idx_recv(i)
              do j = 1, len_dim1
                do k = 1, ndim3
                  buff_put_dp(k,j,lidx) = buff_var(k+(j-1)*ndim3+ndim3*len_dim1*(i-1))
                enddo
              enddo
           end do
         elseif (r8 .eq. 4) then
           allocate(buff_put_sp(ndim3, len_dim1, mycount), source=default_4)
           do i = 1, size(idx_recv)
              lidx = idx_recv(i)
              do j = 1, len_dim1
                do k = 1, ndim3
                  buff_put_sp(k,j,lidx) = buff_var(k+(j-1)*ndim3+ndim3*len_dim1*(i-1))
                enddo
              enddo
           end do
         end if
         start(1) = 1
         start(2) = 1
         count(1) = ndim3
         count(2) = len_dim1
         ! may need to put multiple times (each time < 2GB)
         jres = mod(mycount, max_diml)
         do j = 1, act_put_num
           if (j .lt. act_put_num) then
               count(3) = max_diml
           else
               count(3) = jres
           end if
           start(3) = mystart + (j-1)*max_diml
           jstart   = 1 + (j-1)*max_diml
           if (r8 .eq. 8) then
             call wrap_put_vara_realx_req(ncid, ptr%var_idlist, start(1:3), count(1:3), &
                                          buff_put_dp(:, :, jstart:jstart+count(3)-1), req)
           elseif (r8 .eq. 4) then
             call wrap_put_vara_real_req(ncid, ptr%var_idlist, start(1:3), count(1:3), &
                                         buff_put_sp(:, :, jstart:jstart+count(3)-1),req)
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

     return
   end subroutine wrap_output_3d_group

   subroutine wrap_output_3d_group_sp(mesh, outdir, filename, dim3_len)

   implicit none
!io
     type(global_domain), intent(in), target   :: mesh
     character*(*),       intent(in)           :: outdir
     character*(*),       intent(in)           :: filename
     integer(i4),         intent(in), optional :: dim3_len 
! local
     type(group_comm),          pointer        :: comm
     type(block_structure),     pointer        :: local_block
     type(array_3d_index_list), pointer        :: ptr
     integer(i4)                               :: cmode, ncid, err
     integer(i4)                               :: fst_dim_id(1)           ! dim id for the 1st dim
     integer(i4)                               :: sec_dim_id(3)           ! dim id for the 2nd dim
     integer(i4)                               :: thr_dim_id(3)           ! dim id for the 3rd dim
     integer(i4)                               :: len_dim1, mpi_newtype
     integer(i4)                               :: i0, i1, i2 ,i, j, k, lidx
     integer(i8)                               :: nlevp, nlevi8, ntraceri8
     real(r8),                  allocatable    :: buff_var(:), temp_var(:,:,:)
     real(4),                   allocatable    :: buff_put_sp(:,:,:)
     integer(i4),               pointer        :: displs_send(:), counts_send(:)
     integer(i4),               pointer        :: displs_recv(:), counts_recv(:)
     integer(i4),               pointer        :: idx_recv(:), idx_local(:)
     integer(i4),               pointer        :: mystart, mycount
     integer(i8)                               :: start(3)
     integer(i8)                               :: count(3)
     integer(i4)                               :: req(1), st(1), act_put_num
     integer(i4)                               :: max_diml, jres, jstart, ndim3

     nlevi8 = nlev
     nlevp= nlevi8+1
     ntraceri8 = ntracer ! default is ntracer
     ndim3     = ntracer
     if(present(dim3_len))then
        ntraceri8 = dim3_len
        ndim3     = dim3_len
     end if

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
     if(comm%color==0) then
! possible dims for 1st
       call wrap_def_dim(ncid, 'ntracer'    , ntraceri8   , fst_dim_id(1)) ! gcm
! possible dims for 2nd
       call wrap_def_dim(ncid, 'nlev'       , nlevi8      , sec_dim_id(1)) ! gcm
       call wrap_def_dim(ncid, 'nlevp'      , nlevi8+1    , sec_dim_id(2)) ! gcm
       call wrap_def_dim(ncid, 'dim12'      , int8(12)    , sec_dim_id(3)) ! gcm    
! possible dims for 3rd
       call wrap_def_dim(ncid, 'location_nt', mesh_nt     , thr_dim_id(1)) ! gcm
       call wrap_def_dim(ncid, 'location_ne', mesh_ne     , thr_dim_id(2)) ! gcm
       call wrap_def_dim(ncid, 'location_nv', mesh_nv     , thr_dim_id(3)) ! gcm
     endif
!========================================================
! 2nd: define var, need lon_id, lat_id, var_idlist
!========================================================
     ptr => scalar_3d_list
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
       if(comm%color==0) call def_var_3d_pnetcdf_sp(ptr%varname, fst_dim_id(1), sec_dim_id(i1), thr_dim_id(i2), ncid, ptr%var_idlist)
       ptr => ptr%next
     end do

!================================================
! 3rd: put att
!================================================    
     ptr => scalar_3d_list
     do while (associated(ptr))
       if(comm%color==0) call put_att_3d_pnetcdf( ncid, ptr%var_idlist ,ptr%varname,ptr%longname,ptr%units,ptr%var_pos)
       ptr => ptr%next
     end do 
     if(comm%color==0) call wrap_enddef(ncid)

!================================================
! 4th, put var
!================================================

     ptr => scalar_3d_list
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
       allocate(temp_var(ndim3, len_dim1, size(idx_local)))
       call find_output_var_3d(ptr%var, ptr%var_level, idx_local, temp_var, ndim3)

       if (comm%color == 0) then
         allocate(buff_var(1:mycount*len_dim1*ndim3))
       else
         allocate(buff_var(1))
       endif

       if (r8 .eq. 8) then
         call MPI_Type_contiguous(len_dim1*ndim3, MPI_REAL8, mpi_newtype, err)
       elseif (r8 .eq. 4) then
         call MPI_Type_contiguous(len_dim1*ndim3, MPI_REAL, mpi_newtype, err)
       end if
       call MPI_Type_commit(mpi_newtype, err)
       call MPI_Alltoallv(temp_var, counts_send, displs_send, mpi_newtype, buff_var, &
                          counts_recv, displs_recv, mpi_newtype, comm%all_comm, err)
       deallocate(temp_var)

       if (comm%color == 0) then
         ! get the actual number of outputs (due to the limitation of pnetcdf)
         max_diml = 2*MAXLEN/(ndim3*len_dim1)
         act_put_num = ceiling(mycount/real(max_diml, r8))

         allocate(buff_put_sp(ndim3, len_dim1, mycount), source=default_4)
         do i = 1, size(idx_recv)
            lidx = idx_recv(i)
            do j = 1, len_dim1
              do k = 1, ndim3
                buff_put_sp(k,j,lidx) = buff_var(k+(j-1)*ndim3+ndim3*len_dim1*(i-1))
              enddo
            enddo
         end do
         start(1) = 1
         start(2) = 1
         count(1) = ndim3
         count(2) = len_dim1
         ! may need to put multiple times (each time < 2GB)
         jres = mod(mycount, max_diml)
         do j = 1, act_put_num
           if (j .lt. act_put_num) then
               count(3) = max_diml
           else
               count(3) = jres
           end if
           start(3) = mystart + (j-1)*max_diml
           jstart   = 1 + (j-1)*max_diml
           call wrap_put_vara_real_req(ncid, ptr%var_idlist, start(1:3), count(1:3), &
                                       buff_put_sp(:, :, jstart:jstart+count(3)-1),req)
           err = nfmpi_wait_all(ncid, 1, req, st)
         end do
         deallocate(buff_put_sp)
       endif
       deallocate(buff_var)
       ptr => ptr%next
     end do

     if(comm%color==0) call wrap_close(ncid)

     return
   end subroutine wrap_output_3d_group_sp

   subroutine wrap_output_3d_gather_C(mesh, outdir, filename, dim3_len)
    use perf_mod
    implicit none
!io
     type(global_domain), intent(in), target   :: mesh
     character*(*),       intent(in)           :: outdir
     character*(*),       intent(in)           :: filename
     integer(i4),         intent(in), optional :: dim3_len
! local
     type(group_comm),          pointer        :: comm
     type(block_structure),     pointer        :: local_block
     type(array_3d_index_list), pointer        :: ptr

     integer(i4)                               :: cmode, ncid, err
     integer(i4)                               :: fst_dim_id(1)
     integer(i4)                               :: sec_dim_id(3)
     integer(i4)                               :: thr_dim_id(3)
     integer(i4)                               :: i1, i2, i, nvars
     integer(i4)                               :: len_dim1
     integer(i8)                               :: nlevp, nlevi8, ntraceri8
     integer(i4)                               :: ndim3
     integer(i4)                               :: n, ierr

     character(len=:), allocatable             :: ncfile_str
     integer(i4)                               :: ncfile_len
     character(kind=c_char), allocatable       :: ncfile_c(:)

     ! 通信索引与计数（Fortran侧）
     integer(i4),               pointer        :: displs_send(:), counts_send(:)
     integer(i4),               pointer        :: displs_recv(:), counts_recv(:)
     integer(i4),               pointer        :: idx_recv(:), idx_local(:)
     integer(i4),               pointer        :: mystart, mycount

     ! C 接口所需的 c_int/c_char 转换缓冲
     integer(c_int)                            :: len_dim1_c, len_dim2_c
     integer(c_int)                            :: all_comm_f_c, group_id_c
     integer(c_int)                            :: open_mode_c, local_cols_c
     integer(c_int)                            :: mystart_c, mycount_c
     integer(i4)                               :: nrecv, vlen
     character(kind=c_char), allocatable       :: out_c(:), file_c(:), varname_c(:)
     integer(c_int),          allocatable      :: displs_send_c(:), counts_recv_c(:), displs_recv_c(:), idx_recv_c(:)
     character(kind=c_char)                    :: path_buf(1024)
     character(len=:), allocatable             :: path_f
     character(len=:), allocatable             :: tmp
     integer(c_int)                            :: buf_len

     ! 数据缓冲
     real(r8),                  allocatable    :: temp_var(:,:,:)
     real(4),                   allocatable    :: temp_sp(:,:,:)

     logical                                   :: first_open

      call t_startf('3d_write')
      nlevi8 = nlev
      nlevp= nlevi8+1
      ntraceri8 = ntracer ! default is ntracer
      ndim3     = ntracer
      if(present(dim3_len))then
          ntraceri8 = dim3_len
          ndim3     = dim3_len
      end if

      comm => mesh%gcomm_write
      local_block => mesh%local_block

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

      ptr => scalar_3d_list
      call t_startf('3d_nc_create')
      do while (associated(ptr))

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
!          cmode = IOR(NF_CLOBBER, NF_64BIT_DATA)
!          if (comm%color == 0) call wrap_create(comm%gcomm, ncfile_str, cmode, ncid)
!          ! 定义维度
!          if(comm%color==0) then
!          ! possible dims for 1st
!              call wrap_def_dim(ncid, 'ntracer'    , ntraceri8   , fst_dim_id(1)) ! gcm
!          ! possible dims for 2nd
!              call wrap_def_dim(ncid, 'nlev'       , nlevi8      , sec_dim_id(1)) ! gcm
!              call wrap_def_dim(ncid, 'nlevp'      , nlevi8+1    , sec_dim_id(2)) ! gcm
!              call wrap_def_dim(ncid, 'dim12'      , int8(12)    , sec_dim_id(3)) ! gcm    
!          ! possible dims for 3rd
!              call wrap_def_dim(ncid, 'location_nt', mesh_nt     , thr_dim_id(1)) ! gcm
!              call wrap_def_dim(ncid, 'location_ne', mesh_ne     , thr_dim_id(2)) ! gcm
!              call wrap_def_dim(ncid, 'location_nv', mesh_nv     , thr_dim_id(3)) ! gcm
!          endif
          ! 定义变量
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
          ! 添加属性
!          if(comm%color==0) call def_var_3d_pnetcdf_sp(ptr%varname, fst_dim_id(1), sec_dim_id(i1), thr_dim_id(i2), ncid, ptr%var_idlist)
!          if(comm%color==0) call put_att_3d_pnetcdf( ncid, ptr%var_idlist ,ptr%varname,ptr%longname,ptr%units,ptr%var_pos)
!          ! 结束定义关闭文件
!          if(comm%color==0) call wrap_enddef(ncid)
!          if(comm%color==0) call wrap_close(ncid)
          deallocate(ncfile_c)
          deallocate(varname_c)
          ptr => ptr%next
      end do
      call t_stopf('3d_nc_create')
      ! write(*,*) 'Output 3D scalar variables to separate files using C interface: '//trim(outdir)//trim(filename)//'_varname.nc'
      ! 输出变量信息到汇总文件
      if (mpi_rank() == 0) then
        open(unit=11, file= trim(path_f)//'scalar_3d_list_info.txt', status='replace')
        ! open(unit=11, file= trim(path_f)//'scalar_3d_list_info.txt', status='unknown', position='append')
        
        write(11, '(A)') 'ncfile='//trim(filename)  ! 新增：写入nc文件名
        ! 统计变量个数
        nvars = 0
        ptr => scalar_3d_list
        do while (associated(ptr))
            nvars = nvars + 1
            ptr => ptr%next
        end do
        write(11, '(A, I0)') 'nvars=', nvars
        write(11, '(A)') 'varnames:'
        ptr => scalar_3d_list
        do i = 1, nvars
            write(11, '(A)') trim(ptr%varname)
            ptr => ptr%next
        end do
        close(11)
      end if

! ! 构造 C 字符串（NUL 结尾）
!       allocate(out_c(len_trim(outdir)+1))
!       do i = 1, len_trim(outdir)
!           out_c(i) = transfer(outdir(i:i), out_c(i))
!       end do
!       out_c(len_trim(outdir)+1) = c_null_char

!       allocate(file_c(len_trim(filename)+1))
!       do i = 1, len_trim(filename)
!           file_c(i) = transfer(filename(i:i), file_c(i))
!       end do
!       file_c(len_trim(filename)+1) = c_null_char

      ! 直接传递 Fortran 的 communicator（MPI_Fint）给 C
      all_comm_f_c = int(comm%all_comm, c_int)
      group_id_c = int(comm%color, c_int)

      first_open = .true.
      ptr => scalar_3d_list
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
        allocate(temp_var(ndim3, len_dim1, size(idx_local)))
        call find_output_var_3d(ptr%var, ptr%var_level, idx_local, temp_var, ndim3)
        

        allocate(temp_sp(ndim3, len_dim1, size(idx_local)))
        temp_sp = real(temp_var, kind=c_float)
        deallocate(temp_var)

        len_dim1_c = int(len_dim1, c_int)
        len_dim2_c = int(ndim3, c_int)
        open_mode_c = merge(1_c_int, 0_c_int, first_open)
        local_cols_c = int(size(idx_local),c_int)

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
        ! 调用 C 接口函数

        ! call wrap_output_3d_gather_sp_c(temp_sp, len_dim1_c, len_dim2_c, &
        !                                local_cols_c, group_id_c, all_comm_f_c, &
        !                                out_c, file_c, &
        !                                open_mode_c, &
        !                                mycount_c, mystart_c, &
        !                                displs_send_c, counts_recv_c, &
        !                                displs_recv_c, idx_recv_c, varname_c)
        call wrap_output_3d_gather_sp_c(temp_sp, len_dim1_c, len_dim2_c, &
                                       local_cols_c, group_id_c, all_comm_f_c, &
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
      call t_stopf('3d_write')
   end subroutine wrap_output_3d_gather_C

   subroutine wrap_output_clean_3d()

     type(array_3d_index_list),pointer :: field_list, tmpf
    
     field_list => scalar_3d_list
     do while(associated(field_list))
       tmpf => field_list%next
       deallocate(field_list)
       field_list => tmpf
     end do
     scalar_3d_list => null()

     return
   end subroutine wrap_output_clean_3d
!
! read routines
!

   subroutine wrap_read_3d_group_rst(comm, outdir, filename, varname, dim2_len, &
                                     option, varout, dim3_len)
! io
     type(group_comm), intent(in),    target   :: comm
     character*(*),    intent(in)              :: outdir
     character*(*),    intent(in)              :: filename
     character*(*),    intent(in)              :: varname
     integer(i8),      intent(in)              :: dim2_len
     integer                                   :: option
     real(r8),         intent(inout)           :: varout(:,:,:)
     integer(i8),      intent(in),    optional :: dim3_len
! local
     integer(i4),      parameter               :: omode = 0
     integer(i4)                               :: varid, ncid
     integer(i4)                               :: dim2, dim3, act_get_num, mpi_newtype
     integer(i4),      pointer                 :: displs_send(:), counts_send(:)
     integer(i4),      pointer                 :: displs_recv(:), counts_recv(:)
     integer(i4),      pointer                 :: idx_send(:), idx_local(:)
     integer(i4),      pointer                 :: mystart, mycount
     integer(i8)                               :: start(3)
     integer(i8)                               :: count(3)
     integer(i4)                               :: req(1), st(1)
     real(r8),         allocatable             :: buff_var(:), send_var(:), temp_var(:)
     integer                                   :: i, err, j, xtype
     integer                                   :: max_diml, jres, jstart
     real(4),          allocatable             :: local_data_r4(:,:)
     real(8),          allocatable             :: local_data_r8(:,:)

     dim2 = dim2_len
     dim3 = ntracer
     if(present(dim3_len)) dim3   = dim3_len

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
         max_diml = MAXLEN/(dim3*dim2)
       elseif (xtype .eq. NF_FLOAT) then
         max_diml = 2*MAXLEN/(dim3*dim2)
       end if
       act_get_num = ceiling(mycount/real(max_diml, 8))

       allocate(buff_var(mycount*dim3*dim2))
       start(1) = 1
       start(2) = 1
       count(1) = dim3
       count(2) = dim2_len
       ! may need to get multiple times (each time < 2GB)
       jstart = 1
       jres   = mod(mycount, max_diml)
       do j = 1, act_get_num
         if (j .lt. act_get_num) then
             count(3) = max_diml
         else
             count(3) = jres
         end if
         start(3) = mystart + (j-1)*max_diml
         jstart   = jstart + (j-1)*max_diml
         if (xtype .eq. NF_DOUBLE) then
           allocate(local_data_r8(dim3*dim2, count(3)))
           call wrap_get_vara_realx_req(ncid, varid, start, count, local_data_r8, req)
           err = nfmpi_wait_all(ncid, 1, req, st)
           do i = 1, count(3)
             buff_var((jstart+i-2)*dim3*dim2+1:(jstart+i-1)*dim3*dim2) = local_data_r8(1:dim3*dim2, i)
           end do
           deallocate(local_data_r8)
         elseif (xtype .eq. NF_FLOAT) then
           allocate(local_data_r4(dim3*dim2, count(3)))
           call wrap_get_vara_real_req(ncid, varid, start, count, local_data_r4, req)
           err = nfmpi_wait_all(ncid, 1, req, st)
           do i = 1, count(3)
             buff_var((jstart+i-2)*dim3*dim2+1:(jstart+i-1)*dim3*dim2) = local_data_r4(1:dim3*dim2, i)
           end do
           deallocate(local_data_r4)
         end if
       enddo
       call wrap_close(ncid)
       allocate(send_var(size(idx_send)*dim3*dim2))
       do i = 1, size(idx_send)
         send_var((i-1)*dim3*dim2+1:i*dim3*dim2) = buff_var((idx_send(i)-1)*dim3*dim2+1:idx_send(i)*dim3*dim2)
       end do
       deallocate(buff_var)
     else
       allocate(send_var(1))
     end if

     allocate(temp_var(size(idx_local)*dim2*dim3))

     if (r8 .eq. 8) then
       call MPI_Type_contiguous(dim2*dim3, MPI_REAL8, mpi_newtype, err)
     elseif (r8 .eq. 4) then
       call MPI_Type_contiguous(dim2*dim3, MPI_REAL, mpi_newtype, err)
     end if
     call MPI_Type_commit(mpi_newtype, err)
     call MPI_Alltoallv(send_var, counts_send, displs_send, mpi_newtype, temp_var, &
                        counts_recv, displs_recv, mpi_newtype, comm%all_comm, err)

     do i = 1, size(varout,3)
        do j = 1, dim2
           varout(1:dim3, j, i) = temp_var(1+(j-1)*dim3+(idx_local(i)-1)*dim2*dim3:j*dim3+(idx_local(i)-1)*dim2*dim3)
        end do
     end do

     deallocate(send_var)
     deallocate(temp_var)
     if(mpi_rank().eq.0) write(6,*) "Sucessfully read ",trim(varname)," from ", trim(outdir)//trim(filename)
     return
   end subroutine wrap_read_3d_group_rst

!---------------------------------------------------------------------------
!
!                         Private routines below
!
!---------------------------------------------------------------------------

   subroutine def_var_3d_pnetcdf(varname, fst_dim_id, sec_dim_id, thr_dim_id, ncid, var_idlist)
! io
     character(len=*), intent(in)    :: varname
     integer(i4),      intent(in)    :: fst_dim_id
     integer(i4),      intent(in)    :: sec_dim_id
     integer(i4),      intent(in)    :: thr_dim_id
     integer(i4),      intent(in)    :: ncid
     integer(i4),      intent(inout) :: var_idlist

! local
     integer(i4)                     :: vdims(3)
     integer(i4)                     :: len_thr = 3
    
     vdims(1) = fst_dim_id
     vdims(2) = sec_dim_id
     vdims(3) = thr_dim_id
     if (r8 .eq. 8) then
       call wrap_def_var(ncid,trim(varname),NF_DOUBLE, len_thr, vdims(1:3), var_idlist)
     elseif (r8 .eq. 4) then
       call wrap_def_var(ncid,trim(varname),NF_FLOAT, len_thr, vdims(1:3), var_idlist)
     end if

     return
   end subroutine def_var_3d_pnetcdf

   subroutine def_var_3d_pnetcdf_sp(varname, fst_dim_id, sec_dim_id, thr_dim_id, ncid, var_idlist)
! io
     character(len=*), intent(in)    :: varname
     integer(i4),      intent(in)    :: fst_dim_id
     integer(i4),      intent(in)    :: sec_dim_id
     integer(i4),      intent(in)    :: thr_dim_id
     integer(i4),      intent(in)    :: ncid
     integer(i4),      intent(inout) :: var_idlist

! local
     integer(i4)                     :: vdims(3)
     integer(i4)                     :: len_thr = 3
    
     vdims(1) = fst_dim_id
     vdims(2) = sec_dim_id
     vdims(3) = thr_dim_id
     call wrap_def_var(ncid,trim(varname),NF_FLOAT, len_thr, vdims(1:3), var_idlist)

     return
   end subroutine def_var_3d_pnetcdf_sp

   subroutine put_att_3d_pnetcdf(ncid, var_idlist, varname,longname,units,var_pos)
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

   end subroutine put_att_3d_pnetcdf

   subroutine find_output_var_3d(local_var,level_var,output_local_index,output_local_var,ndim3)
! io
     real(r8), pointer, intent(in)   :: local_var(:,:,:)
     integer(i8),       intent(in)   :: level_var
     integer(i4),       intent(in)   :: output_local_index(:)
     real(r8),          intent(out)  :: output_local_var(:,:,:)
     integer(i4),       intent(in)   :: ndim3
  
! local
     integer(i8) :: idex,idey
     
     do idex = 1,size(output_local_index)
        do idey = 1,level_var 
           output_local_var(1:ndim3,idey,idex) = local_var(1:ndim3,idey,output_local_index(idex))
        end do
     end do

   end subroutine find_output_var_3d

 end module grist_fileio_list_3d_module_par
