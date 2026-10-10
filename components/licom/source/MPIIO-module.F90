 !|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

 module MPIIO

   use precision_mod
   use param_mod
   use msg_mod
   use blocks
   use distribution
   use domain
   use constant_mod
   use LICOM_Error_mod
   use module_add
   implicit none
   private
   save

   public ::  mpiio_write,mpiio_read

   interface mpiio_write
     module procedure mpiio_write_dbl
	 module procedure mpiio_write_real
   end interface

   interface mpiio_read
     module procedure mpiio_read_dbl
     module procedure mpiio_read_real	 
   end interface

 contains


 
 !***********************************************************************

subroutine mpiio_write_dbl(fh,ARRAY,this_block,src_dist,subarray_type,date_count,all_offset_km) 
    include 'mpif.h'
    real (r8), dimension(:,:,:), intent(in) :: ARRAY      ! array containing horizontal slab of distributed field
    type (distrb), intent(in) :: src_dist   ! distribution of blocks in the source array
    type (block),intent(in)::this_block  ! block info for current block	
	integer,intent(in):: subarray_type,fh,date_count	
    integer*8, intent(inout)::all_offset_km
    integer(kind=MPI_OFFSET_KIND)::middle_offset
    integer (int_kind), dimension(MPI_STATUS_SIZE) :: status

!-----------------------------------------------------------------------
	  
		middle_offset=all_offset_km
		call MPI_FILE_SET_VIEW(fh,middle_offset,mpi_dbl,subarray_type,"native",mpi_info_null,ierror)
		call mpi_file_write_all(fh,ARRAY(this_block%ib:this_block%ie,this_block%jb:this_block%je, &
			src_dist%local_block(n)),date_count,mpi_dbl,status,ierror)
	  all_offset_km=all_offset_km+imt_global*jmt_global*8

end subroutine mpiio_write_dbl
 
!***********************************************************************

subroutine mpiio_read_dbl(fh,buffer,all_offset)
    include 'mpif.h'
    real (r8), dimension(:,:), intent(inout) :: buffer     ! array containing horizontal slab of distributed field	
    integer*8, intent(inout)::all_offset
	integer,intent(in):: fh
	integer:: date_count	
    integer(kind=MPI_OFFSET_KIND)::middle_offset
	integer whaterror,errr 

		date_count=imt_global*jmt_global
        middle_offset=all_offset
		call MPI_FILE_READ_AT(fh,middle_offset,buffer,date_count,mpi_dbl,status,ierror)
		call MPI_Error_class(ierror,whaterror,errr)
		if(ierror /= 0) write(*,*) "read is wrong",mytid,whaterror

		all_offset=all_offset+imt_global*jmt_global*8
end subroutine mpiio_read_dbl




subroutine mpiio_write_real(fh,ARRAY,this_block,src_dist,subarray_type,date_count,all_offset_km) 
    include 'mpif.h'
    real (r4), dimension(:,:,:), intent(in) :: ARRAY      ! array containing horizontal slab of distributed field
    type (distrb), intent(in) :: src_dist   ! distribution of blocks in the source array
    type (block),intent(in)::this_block  ! block info for current block	
	integer,intent(in):: subarray_type,fh,date_count	
    integer*8, intent(inout)::all_offset_km
    integer(kind=MPI_OFFSET_KIND)::middle_offset
    integer (int_kind), dimension(MPI_STATUS_SIZE) :: status

!-----------------------------------------------------------------------
	  
		middle_offset=all_offset_km
		call MPI_FILE_SET_VIEW(fh,middle_offset,MPI_REAL,subarray_type,"native",mpi_info_null,ierror)
		call mpi_file_write_all(fh,ARRAY(this_block%ib:this_block%ie,this_block%jb:this_block%je, &
			src_dist%local_block(n)),date_count,MPI_REAL,status,ierror)
	  all_offset_km=all_offset_km+imt_global*jmt_global*4

end subroutine mpiio_write_real
 
!***********************************************************************

subroutine mpiio_read_real(fh,buffer,all_offset)
    include 'mpif.h'
    real (r4), dimension(:,:), intent(inout) :: buffer     ! array containing horizontal slab of distributed field	
    integer*8, intent(inout)::all_offset
	integer,intent(in):: fh
	integer:: date_count	
    integer(kind=MPI_OFFSET_KIND)::middle_offset
	integer whaterror,errr 

		date_count=imt_global*jmt_global
        middle_offset=all_offset
		call MPI_FILE_READ_AT(fh,middle_offset,buffer,date_count,MPI_REAL,status,ierror)
		call MPI_Error_class(ierror,whaterror,errr)
		if(ierror /= 0) write(*,*) "read is wrong",mytid,whaterror

		all_offset=all_offset+imt_global*jmt_global*4
end subroutine mpiio_read_real
 
 end module MPIIO

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
