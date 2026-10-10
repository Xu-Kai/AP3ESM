
#define  UN_VAR_CHECK
module var_check
use precision_mod
use msg_mod, only: tag_1d,tag_2d,tag_3d,tag_4d,nproc,status,mpi_comm_ocn
use blocks
use pconst_mod

public data_check, data_check_lz

interface data_check 
    module procedure data_check3d,data_check2d,data_check1d,data_check4d
end interface


contains

subroutine  data_check4d(var_in, info_var) 

 IMPLICIT NONE
    real(r8), dimension(:,:,:,:), intent(in) :: var_in
    character(*), intent(in), optional :: info_var

    real(r8) :: sq_sum, local_sq_sum
    integer :: nnn,iblock !,month1 !LPF20200315

    type (block) :: this_block
    integer :: j1, j2, i1, i2, i_id, j_id, k_id, k1, k2

    ! if(my_task == master_task) write(6,*) "nblocks_clinic is ",  nblocks_clinic

    sq_sum = 0.0
    local_sq_sum = 0.0
#ifdef VAR_CHECK
    ! this_block = get_block(blocks_clinic(1),1)
    iblock = 1
    i1 = lbound(var_in, dim=1)
    i2 = ubound(var_in, dim=1)
    j1 = lbound(var_in, dim=2)
    j2 = ubound(var_in, dim=2)
    k1 = lbound(var_in, dim=3)
    k2 = ubound(var_in, dim=3)
    do k_id=k1, k2
    DO J_id = j1, j2
        DO I_id = i1, i2
            local_sq_sum = local_sq_sum + var_in(i_id,J_id,k_id, iblock)* var_in (i_id,J_id,k_id,iblock)
        END DO
    END DO
    enddo
    call mpi_reduce(local_sq_sum,sq_sum,1,MPI_PR,mpi_sum,0,mpi_comm_ocn,ierr)

    if (mytid==0 )then
        if(present(info_var)) then
            WRITE (16,FMT='(I5,I3,D25.15,A25)') MONTH,IDAY, sq_sum, info_var
        else
            WRITE (16,FMT='(I5,I3,D25.15)') MONTH,IDAY, sq_sum
        endif 
        call flush(16)
    endif

    sq_sum = 0.0
    local_sq_sum = 0.0

    do k_id=k1, k2
    DO J_id = j1, j2
        DO I_id = i1, i2
            if(var_in(i_id,J_id,k_id,iblock) .le. 0.0)  local_sq_sum = local_sq_sum + var_in(i_id,J_id,k_id, iblock)* var_in (i_id,J_id,k_id,iblock)
        END DO
    END DO
    enddo
    call mpi_reduce(local_sq_sum,sq_sum,1,MPI_PR,mpi_sum,0,mpi_comm_ocn,ierr)
    if (mytid==0 )then
        if(present(info_var)) then
            WRITE (16,FMT='(I5,I3,D25.15,A25)') MONTH,IDAY, sq_sum, info_var
        else
            WRITE (16,FMT='(I5,I3,D25.15)') MONTH,IDAY, sq_sum
        endif 
        call flush(16)
    endif
#endif
end subroutine


subroutine  data_check_lz(var_in, info_var) 

      IMPLICIT NONE
    real(r8), dimension(:,:,:), intent(in) :: var_in
    character(*), intent(in), optional :: info_var

    real(r8) :: sq_sum, local_sq_sum
    integer :: nnn,iblock !,month1 !LPF20200315

    type (block) :: this_block
    integer :: j1, j2, i1, i2, i_id, j_id

    ! if(my_task == master_task) write(6,*) "nblocks_clinic is ",  nblocks_clinic

    sq_sum = 0.0
    local_sq_sum = 0.0
#ifdef VAR_CHECK

    ! this_block = get_block(blocks_clinic(1),1)
    iblock = 1
    i1 = lbound(var_in, dim=1)
    i2 = ubound(var_in, dim=1)
    j1 = lbound(var_in, dim=2)
    j2 = ubound(var_in, dim=2)
    DO J_id = j1, j2
        DO I_id = i1, i2
            if(var_in(i_id,J_id,iblock) .le. 0.0)  local_sq_sum = local_sq_sum + var_in(i_id,J_id,iblock)* var_in (i_id,J_id,iblock)
        END DO
    END DO
    
    call mpi_reduce(local_sq_sum,sq_sum,1,MPI_PR,mpi_sum,0,mpi_comm_ocn,ierr)

    if (mytid==0 )then
        if(present(info_var)) then
            WRITE (16,FMT='(I5,I3,D25.15,A25)') MONTH,IDAY, sq_sum, info_var
        else
            WRITE (16,FMT='(I5,I3,D25.15)') MONTH,IDAY, sq_sum
        endif 
        call flush(16)
    endif
#endif
end subroutine


subroutine  data_check3d(var_in, info_var) 
! #include <def-undef.h>
! use precision_mod
! use param_mod
! use pconst_mod
! use tracer_mod
! use dyn_mod
! use work_mod
! use grid
! use constant_mod
! use msg_mod, only: tag_1d,tag_2d,tag_3d,tag_4d,nproc,status,mpi_comm_ocn
! use shr_sys_mod
! use domain
! use blocks
      IMPLICIT NONE
    real(r8), dimension(:,:,:), intent(in) :: var_in
    character(*), intent(in), optional :: info_var

    real(r8) :: sq_sum, local_sq_sum
    integer :: nnn,iblock !,month1 !LPF20200315

    type (block) :: this_block
    integer :: j1, j2, i1, i2, i_id, j_id

    ! if(my_task == master_task) write(6,*) "nblocks_clinic is ",  nblocks_clinic

    sq_sum = 0.0
    local_sq_sum = 0.0
#ifdef VAR_CHECK

    ! this_block = get_block(blocks_clinic(1),1)
    iblock = 1
    i1 = lbound(var_in, dim=1)
    i2 = ubound(var_in, dim=1)
    j1 = lbound(var_in, dim=2)
    j2 = ubound(var_in, dim=2)
    DO J_id = j1, j2
        DO I_id = i1, i2
            local_sq_sum = local_sq_sum + var_in(i_id,J_id,iblock)* var_in (i_id,J_id,iblock)
        END DO
    END DO
    
    call mpi_reduce(local_sq_sum,sq_sum,1,MPI_PR,mpi_sum,0,mpi_comm_ocn,ierr)

    if (mytid==0 )then
        if(present(info_var)) then
            WRITE (16,FMT='(I5,I3,D25.15,A25)') MONTH,IDAY, sq_sum, info_var
        else
            WRITE (16,FMT='(I5,I3,D25.15)') MONTH,IDAY, sq_sum
        endif 
        call flush(16)
    endif
#endif
end subroutine
subroutine  data_check2d(var_in, info_var) 
! #include <def-undef.h>
! use precision_mod
! use param_mod
! use pconst_mod
! use tracer_mod
! use dyn_mod
! use work_mod
! use grid
! use constant_mod
! use msg_mod, only: tag_1d,tag_2d,tag_3d,tag_4d,nproc,status,mpi_comm_ocn
! use shr_sys_mod
! use domain
! use blocks
      IMPLICIT NONE
    real(r8), dimension(:,:), intent(in) :: var_in
    character(*), intent(in), optional :: info_var

    real(r8) :: sq_sum, local_sq_sum
    integer :: nnn,iblock !,month1 !LPF20200315

    type (block) :: this_block
    integer :: j1, j2, i1, i2, i_id, j_id

    ! if(my_task == master_task) write(6,*) "nblocks_clinic is ",  nblocks_clinic

    sq_sum = 0.0
    local_sq_sum = 0.0
#ifdef VAR_CHECK

    ! this_block = get_block(blocks_clinic(1),1)
    iblock = 1
    i1 = lbound(var_in, dim=1)
    i2 = ubound(var_in, dim=1)
    j1 = lbound(var_in, dim=2)
    j2 = ubound(var_in, dim=2)
    DO J_id = j1, j2
        DO i_id = i1, i2
            local_sq_sum = local_sq_sum + var_in(i_id,J_id)* var_in (i_id,J_id)
        END DO
    END DO
    
    call mpi_reduce(local_sq_sum,sq_sum,1,MPI_PR,mpi_sum,0,mpi_comm_ocn,ierr)
    if (mytid==0 )then

        if(present(info_var)) then
            WRITE (16,FMT='(I5,I3,D25.15,A25)') MONTH,IDAY, sq_sum, info_var
        else
            WRITE (16,FMT='(I5,I3,D25.15)') MONTH,IDAY, sq_sum
        endif 
        call flush(16)
    endif
#endif
end subroutine

subroutine  data_check1d(var_in, info_var) 
! #include <def-undef.h>
! use precision_mod
! use param_mod
! use pconst_mod
! use tracer_mod
! use dyn_mod
! use work_mod
! use grid
! use constant_mod
! use msg_mod, only: tag_1d,tag_2d,tag_3d,tag_4d,nproc,status,mpi_comm_ocn
! use shr_sys_mod
! use domain
! use blocks
      IMPLICIT NONE
    real(r8), dimension(:), intent(in) :: var_in
    character(*), intent(in), optional :: info_var

    real(r8) :: sq_sum, local_sq_sum
    integer :: nnn,iblock !,month1 !LPF20200315

    type (block) :: this_block
    integer :: j1, j2, i1, i2, i_id, j_id

    ! if(my_task == master_task) write(6,*) "nblocks_clinic is ",  nblocks_clinic

    sq_sum = 0.0
    local_sq_sum = 0.0
#ifdef VAR_CHECK

    ! this_block = get_block(blocks_clinic(1),1)
    iblock = 1
    i1 = lbound(var_in, dim=1)
    i2 = ubound(var_in, dim=1)

        DO i_id = i1, i2
            local_sq_sum = local_sq_sum + var_in(i_id)* var_in (i_id)
        END DO
    
    call mpi_reduce(local_sq_sum,sq_sum,1,MPI_PR,mpi_sum,0,mpi_comm_ocn,ierr)
    if (mytid==0 )then

        if(present(info_var)) then
            WRITE (16,FMT='(I5,I3,D25.15,A25)') MONTH,IDAY, sq_sum, info_var
        else
            WRITE (16,FMT='(I5,I3,D25.15)') MONTH,IDAY, sq_sum
        endif 
        call flush(16)
    endif
    sq_sum = 0.0
    local_sq_sum = 0.0
    ! this_block = get_block(blocks_clinic(1),1)
    iblock = 1
    i1 = lbound(var_in, dim=1)
    i2 = ubound(var_in, dim=1)

        DO i_id = i1, i2
            if(var_in(i_id) .le. 0.0) local_sq_sum = local_sq_sum + var_in(i_id)* var_in (i_id)
        END DO
    
    call mpi_reduce(local_sq_sum,sq_sum,1,MPI_PR,mpi_sum,0,mpi_comm_ocn,ierr)
    if (mytid==0 )then

        if(present(info_var)) then
            WRITE (16,FMT='(I5,I3,D25.15,A25)') MONTH,IDAY, sq_sum, info_var
        else
            WRITE (16,FMT='(I5,I3,D25.15)') MONTH,IDAY, sq_sum
        endif 
        call flush(16)
    endif
#endif
end subroutine

end module