#include <stdlib.h>
#include <stdio.h>
#include <libgen.h>
#include <mpi.h>

static FILE *f_cice_ptr;
static long c_my_task;
static int file_block_num = 100;
static int my_rank;
static int w_size;
static char file_name[256];
static MPI_Comm comm_ice;
int *remap_ice_blocks = NULL;


void init_ice_read_map_(int *remap)
{
    remap_ice_blocks = remap;
}

void init_ice_com_(MPI_Fint *comm_fortran){
    comm_ice = MPI_Comm_f2c(*comm_fortran);
}


void ice_fopen_restart_read_(char* parent_dir, char *prefix, long* rec_len, int nprocall)
{
    // my_rank = *mytask;
    MPI_Comm_size(comm_ice, &w_size);
    MPI_Comm_rank(comm_ice, &my_rank);
    int block_id = 0;
    block_id = my_rank;

    int file_id = block_id / 100;
    long file_offset = block_id % 100;

    //int dir_no = ((w_size + 99)/ 100) * 100;
    sprintf(file_name, "%s_%06d/%s_%04d", parent_dir, nprocall, prefix, file_id);

    if(my_rank == 0)
        printf("file name is %s, nprocall:%d\n", file_name,nprocall);
    f_cice_ptr = fopen(file_name, "r");
    if(f_cice_ptr == NULL) printf("ERROR: can't open file %s\n", file_name);

    long offset = file_offset * (*rec_len);
    // fseek(f_cice_ptr, 0, SEEK_END);
    // file_len = ftell(f_cice_ptr);
    int s_ret = fseek(f_cice_ptr, offset, SEEK_SET);
    if(s_ret != 0) printf("ERROR: file %s fseek error %ld\n", file_name, offset);

    // cur_block_size = file_len - offset;

    // if(my_rank < 100)
    // {
    //   printf("file name is %s offset %ld %ld %ld\n", file_name, *rec_len,  offset, cur_block_size);
    // }


}


void ice_fopen_read_(char* parent_dir, char *prefix, long* rec_len, long *nprocall)
{
    // my_rank = *mytask;
    MPI_Comm_size(comm_ice, &w_size);
    MPI_Comm_rank(comm_ice, &my_rank);
    int block_id = 0;
    if(remap_ice_blocks == NULL) {
        block_id = my_rank;
        if(my_rank == 0)
            printf("remap_ice_blocks is null\n");
    }
    else {
        block_id = remap_ice_blocks[my_rank];
//        printf("remap_ice_blocks my_rank:%d,block_id:%d\n",my_rank,block_id);
    }

    int file_id = block_id / 100;
    long file_offset = block_id % 100;
    int localnrpoc = (*nprocall);

   //int dir_no = ((w_size + 99)/ 100) * 100; 
    sprintf(file_name, "%s_%06d/%s_%04d", parent_dir, localnrpoc, prefix, file_id);

    if(my_rank == 0)
      printf("ice file name is %s, nprocall:%d, w_size:%d\n", file_name,localnrpoc,w_size);


    f_cice_ptr = fopen(file_name, "r");
    if(f_cice_ptr == NULL) printf("ERROR: can't open file %s\n", file_name);

    long offset = file_offset * (*rec_len);
    // fseek(f_cice_ptr, 0, SEEK_END);
    // file_len = ftell(f_cice_ptr);
    int s_ret = fseek(f_cice_ptr, offset, SEEK_SET);
    if(s_ret != 0) printf("ERROR: ice_fopen_read file %s fseek error %d,%ld, %ld\n", file_name, block_id, (*rec_len),offset);

    // cur_block_size = file_len - offset;

    // if(my_rank < 100)
    // {
    //   printf("file name is %s offset %ld %ld %ld\n", file_name, *rec_len,  offset, cur_block_size);
    // }


}

void ice_fread_( void *data, long* data_size)
{
    size_t data_read = *data_size;
    size_t r_size = fread(data, 1, data_read, f_cice_ptr);


    // if(data_read > cur_block_size) printf("READ ERROR: data read  %ld > block size %d\n", data_read, cur_block_size);


     if(my_rank == 0) printf("ice_fread data size is %ld\n", data_read);
    if(r_size != data_read) printf("read error %s %ld %ld\n", file_name, data_read, r_size);

}


void ice_seek_fread_( void *data, long* rec_offset,  long* data_size)
{
    size_t data_read = *data_size;


    int block_id = 0;
    if(remap_ice_blocks == NULL)
      block_id = my_rank;
    else
      block_id = remap_ice_blocks[my_rank];


    int file_id = block_id / 100;
    long file_offset = block_id % 100;


    long offset = file_offset * (*data_size) + (*rec_offset);

    // long offset = *rec_offset;
    int s_ret = fseek(f_cice_ptr, offset, SEEK_SET);
    if(s_ret != 0) printf("ERROR: file %s fseek error %ld\n", file_name, offset);

    size_t r_size = fread(data, 1, data_read, f_cice_ptr);


    // if(data_read > cur_block_size) printf("READ ERROR: data read  %ld > block size %d\n", data_read, cur_block_size);


    // if(my_rank == 0) printf("data size is %ld\n", data_read);
    if(r_size != data_read) printf("read error %s %ld %ld\n", file_name, data_read, r_size);

}


void ice_fclose_()
{
  fclose(f_cice_ptr);
}
