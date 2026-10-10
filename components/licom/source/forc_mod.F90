!  CVS: $Id: forc_mod.F90,v 1.7 2003/08/25 07:47:52 lhl Exp $
module forc_mod
#include <def-undef.h>
!
use precision_mod
use param_mod
!
!
!     ------------------------------------------------------------------
!     Forcing Fields
!     ------------------------------------------------------------------
      real(r8),allocatable,dimension(:,:,:,:) :: su3,sv3,psa3,tsa3,qar3,uva3, &
                   swv3,cld3,sss3,sst3 ,nswv3,dqdt3,chloro3,                  &
                   wspd3,wspdu3,wspdv3,lwv3,seaice3,rain3,snow3,runoff3
!lhl090730
      real(r8),dimension(imt,jmt,max_blocks_clinic),target ::su,sv,psa,tsa,sss,swv,uva,qar,cld,&
                    ddd,qqq,sst, nswv,dqdt,chloro,lwv,seaice,rain,snow,fresh,runoff, &
                    lthf,sshf !only for output
      real(r8),dimension(imt,jmt,max_blocks_clinic),target ::USTAR,BUOYTUR,BUOYSOL
!
      real(r8),allocatable,dimension(:,:,:)::su3_io,sv3_io,psa3_io,tsa3_io,qar3_io,uva3_io, &
                    swv3_io,cld3_io,sss3_io,sst3_io,nswv3_io,dqdt3_io,chloro3_io, &
                    wspdu3_io,wspdv3_io,lwv3_io,seaice3_io,rain3_io,snow3_io,runoff3_io
!
!       real(r8),dimension(imt,jmt,km,ntra,max_blocks_clinic)::restore
       real(r8),dimension(imt,jmt,max_blocks_clinic),target:: tsf,ssf
       real(r8),dimension(:,:,:,:,:),allocatable::restore
!       real(r8),dimension(:,:,:),allocatable:: tsf,ssf
        real(r8),allocatable,dimension(:,:,:) :: buffer2d
#if ( defined TIDEFORCEGW )
        real(r8),allocatable,dimension(:,:,:),target::depthstd 
        real(r8),allocatable,dimension(:,:,:),target::CD_InWave 
#endif  
#if ( defined TIDEMIX )
      real(r8),allocatable,dimension(:,:,:),target:: wave_dis
#endif
! avoid > 2GB common for 5km version
      real(r8),allocatable,dimension(:,:,:) :: t10,u10,v10,slp,q10,swhf,lwhf
      real(r8),allocatable,dimension(:,:,:) :: precr,precs,rf,si
      real(r8),allocatable,dimension(:,:,:,:) :: buffer3d
      real(r8),allocatable,dimension(:,:,:) :: w3d

      integer*8 :: max_n_bilin_weight, max_n_conserve_weight,max_n_runoff_weight
      integer(4), allocatable ::  block_index(:,:,:)
      integer(8), allocatable :: bilin_src_coord(:), conserve_src_coord(:),runoff_src_coord(:)
      real(8), allocatable :: tmp_u(:,:,:,:), tmp_v(:,:,:,:), local_rho(:,:), bilin_weight(:), conserve_weight(:),runoff_weight(:)
      integer,save:: jra_flag=0
      real(r8), allocatable :: local_lon_o(:, :), local_lat_o(:, :)
      real(r8), allocatable :: nw_lat(:, :)

contains
subroutine forc_init()
      implicit none
      integer*8 :: data_size, offset
      integer :: irec
      allocate(restore(imt,jmt,km,ntra,max_blocks_clinic))
!      allocate(tsf(imt,jmt,max_blocks_clinic))
!      allocate(ssf(imt,jmt,max_blocks_clinic))
      allocate(block_index(nx_block,ny_block,6))
      allocate(tmp_u(nx_block,ny_block,12,max_blocks_clinic))
      allocate(tmp_v(nx_block,ny_block,12,max_blocks_clinic))
      allocate(local_rho(nx_block,ny_block))

      allocate(local_lon_o(nx_block,ny_block))
      allocate(local_lat_o(nx_block,ny_block))
      allocate(nw_lat(2, jmt_global))

      data_size = nx_block*ny_block* sizeof(local_lon_o(1,1))
      call c_fopen_read("ocn"//CHAR(0), "lon.dir/lon"//CHAR(0), data_size)
      data_size = nx_block*ny_block* sizeof(local_lon_o(1,1))
      call c_fread(local_lon_o(1,1), data_size)
      call c_fclose()

      data_size = nx_block*ny_block* sizeof(local_lat_o(1,1))
      call c_fopen_read("ocn"//CHAR(0), "lat.dir/lat"//CHAR(0), data_size)
      data_size = nx_block*ny_block* sizeof(local_lat_o(1,1))
      call c_fread(local_lat_o(1,1), data_size)
      call c_fclose()

      offset = 0
      call c_fopen_read("ocn"//CHAR(0), "nw_lat.dir/nw_lat"//CHAR(0), offset)
      data_size = 2 * jmt_global * sizeof(nw_lat(1,1))
      call c_fread(nw_lat(1,1), data_size)
      call c_fclose()

      data_size = nx_block*ny_block* sizeof(local_rho(1,1))
      call c_fopen_read("ocn"//CHAR(0), "rho.dir/rho"//CHAR(0), data_size)
      data_size = nx_block*ny_block* sizeof(local_rho(1,1))
      call c_fread(local_rho(1,1), data_size)
      call c_fclose()

      data_size = nx_block*ny_block* sizeof(block_index(1,1,1)) * 6
      call c_fopen_read("ocn"//CHAR(0), "block_coord.dir/block_coord"//CHAR(0), data_size)
      data_size = nx_block*ny_block* sizeof(block_index(1,1,1)) * 6
      call c_fread(block_index(1,1,1), data_size)
      call c_fclose()


      offset = 0
      call c_fopen_read("ocn"//CHAR(0), "bilin_src_coord.dir/bilin_src_coord"//CHAR(0), offset)
      data_size = 8
      call c_fread(max_n_bilin_weight,  data_size)

       if(mytid .eq. 0) write(6, *) "bilin src coord max_n_bilin_weight is ", max_n_bilin_weight

      allocate(bilin_src_coord(max_n_bilin_weight))
      offset = 8
      data_size = max_n_bilin_weight * sizeof(bilin_src_coord(1))
      call c_seek_fread(bilin_src_coord(1), offset, data_size)
      call c_fclose()

      offset = 0
      call c_fopen_read("ocn"//CHAR(0), "bilin_weight.dir/bilin_weight"//CHAR(0), offset)
      data_size = 8
      call c_fread(max_n_bilin_weight, data_size)
      allocate(bilin_weight(max_n_bilin_weight))

      offset = 8
      data_size = max_n_bilin_weight * sizeof(bilin_weight(1))
      call c_seek_fread(bilin_weight(1), offset, data_size)
      call c_fclose()

!      offset = 0
!      call c_fopen_read("ocn"//CHAR(0), "conserve_src_coord.dir/conserve_src_coord"//CHAR(0), offset)
!      data_size = 8
!      call c_fread(max_n_conserve_weight,  data_size)

!       if(mytid .eq. 0) write(6, *) "conserve src coord max_n_conserve_weight is ", max_n_conserve_weight
!
!      allocate(conserve_src_coord(max_n_conserve_weight))
!      offset = 8
!      data_size = max_n_conserve_weight * sizeof(conserve_src_coord(1))
!      call c_seek_fread(conserve_src_coord(1), offset, data_size)
!      call c_fclose()
!
!      offset = 0
!      call c_fopen_read("ocn"//CHAR(0), "conserve_weight.dir/conserve_weight"//CHAR(0), offset)
!      data_size = 8
!      call c_fread(max_n_conserve_weight, data_size)
!      allocate(conserve_weight(max_n_conserve_weight))
!
!      offset = 8
!      data_size = max_n_conserve_weight * sizeof(conserve_weight(1))
!      call c_seek_fread(conserve_weight(1), offset, data_size)
!      call c_fclose()

      offset = 0
      call c_fopen_read("ocn"//CHAR(0), "runoff_src_coord.dir/runoff_src_coord"//CHAR(0), offset)
      data_size = 8
      call c_fread(max_n_runoff_weight,  data_size)

       if(mytid .eq. 0) write(6, *) "runoff src coord max_n_runoff_weight is ", max_n_runoff_weight

      allocate(runoff_src_coord(max_n_runoff_weight))
      offset = 8
      data_size = max_n_runoff_weight * sizeof(runoff_src_coord(1))
      call c_seek_fread(runoff_src_coord(1), offset, data_size)
      call c_fclose()

      offset = 0
      call c_fopen_read("ocn"//CHAR(0), "runoff_weight.dir/runoff_weight"//CHAR(0), offset)
      data_size = 8
      call c_fread(max_n_runoff_weight, data_size)
      allocate(runoff_weight(max_n_runoff_weight))

      offset = 8
      data_size = max_n_runoff_weight * sizeof(runoff_weight(1))
      call c_seek_fread(runoff_weight(1), offset, data_size)
      call c_fclose()





end subroutine






end module forc_mod
