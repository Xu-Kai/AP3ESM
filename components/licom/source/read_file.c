#include <stdlib.h>
#include <stdio.h>
#include <libgen.h>

// void read_temp_file_(long *mytask, double *data, long* data_size)
// {
//     char file_name[256];
//     sprintf(file_name, "temp.dir/temp_initial_%04d", (*mytask) + 1);
//     FILE * f = fopen(file_name, "r");
//     fread(data, 1, *data_size, f);
//     fclose(f); 
// }



// // void write_temp_file_(long *mytask, double *data, int* data_size)
// // {
// //     char file_name[256];
// //     sprintf(file_name, "temp.dir/temp_initial_%04d", (*mytask) + 1);
// //     FILE * f = fopen(file_name, "w");
// //     fwrite(data, 1, *data_size, f);
// //     fclose(f); 
// // }


// void read_salt_file_(long *mytask, double *data, long* data_size)
// {
//     char file_name[256];
//     sprintf(file_name, "salt.dir/salt_initial_%04d", (*mytask) + 1);
//     FILE * f = fopen(file_name, "r");
//     fread(data, 1, *data_size, f);
//     fclose(f); 
// }





#include <mpi.h>
#include "def-undef.h"

static FILE *f_ptr;
static long c_my_task;
static int file_block_num = 100;
static int my_rank;
static int w_size;
static char file_name[256];
static MPI_Comm comm_ocn=MPI_COMM_WORLD;

// FILE *f_ptr;
// long c_my_task;

// void c_fopen_read_split_(char *prefix,  long* mytask)
// {
//      char file_name[256];
//     sprintf(file_name, "%s_%04d", prefix, *mytask);
//     c_my_task = *mytask;
//     // char *path;
//     // path = strdup(file_name);
//     // dirname(path);
//     // if(NULL==opendir(path))
//     //   mkdir(path, 0775);
//     if(c_my_task == 0)
//       printf("file name is %s\n", file_name);
//     f_ptr = fopen(file_name, "r");
//     if(f_ptr == NULL) printf("ERROR: can't open file %s\n", file_name);
// }


// long file_len;
// long cur_block_size;

int *remap_ocn_blocks = NULL;
void init_read_map_(int *remap)
{
  remap_ocn_blocks = remap;
}

void init_ocn_com_(MPI_Fint *comm_fortran){
    comm_ocn = MPI_Comm_f2c(*comm_fortran);
}

void c_fopen_restart_read_(char* parent_dir, char *prefix, long* rec_len)
{
    // my_rank = *mytask;
    MPI_Comm_size(comm_ocn, &w_size);
    MPI_Comm_rank(comm_ocn, &my_rank);
    int block_id = 0;
    block_id = my_rank;

    int file_id = block_id / 100;
    long file_offset = block_id % 100;

    //int dir_no = ((w_size + 99)/ 100) * 100;
    sprintf(file_name, "%s_%06d/%s_%04d", parent_dir, NPROC_ALL, prefix, file_id);

    if(my_rank == 0)
        printf("file name is %s,w_size:%d\n", file_name,w_size);
    f_ptr = fopen(file_name, "r");
    if(f_ptr == NULL) printf("ERROR: can't open file %s\n", file_name);

    long offset = file_offset * (*rec_len);
    // fseek(f_ptr, 0, SEEK_END);
    // file_len = ftell(f_ptr);
    int s_ret = fseek(f_ptr, offset, SEEK_SET);
    if(s_ret != 0) printf("ERROR: file %s fseek error %ld\n", file_name, offset);

    // cur_block_size = file_len - offset;

    // if(my_rank < 100)
    // {
    //   printf("file name is %s offset %ld %ld %ld\n", file_name, *rec_len,  offset, cur_block_size);
    // }


}


void c_fopen_read_(char* parent_dir, char *prefix, long* rec_len)
{
    // my_rank = *mytask;
    MPI_Comm_size(comm_ocn, &w_size);
    MPI_Comm_rank(comm_ocn, &my_rank);
    int block_id = 0;
    if(remap_ocn_blocks == NULL) 
      block_id = my_rank;
    else
      block_id = remap_ocn_blocks[my_rank];

    int file_id = block_id / 100;
    long file_offset = block_id % 100;

   //int dir_no = ((w_size + 99)/ 100) * 100; 
    sprintf(file_name, "%s_%06d/%s_%04d", parent_dir, NPROC_ALL, prefix, file_id);

    if(my_rank == 0)
      printf("file name is %s,w_size:%d\n", file_name,w_size);
    f_ptr = fopen(file_name, "r");
    if(f_ptr == NULL) printf("ERROR: can't open file %s\n", file_name);

    long offset = file_offset * (*rec_len);
    // fseek(f_ptr, 0, SEEK_END);
    // file_len = ftell(f_ptr);
    int s_ret = fseek(f_ptr, offset, SEEK_SET);
    if(s_ret != 0) printf("ERROR: file %s fseek error %ld\n", file_name, offset);

    // cur_block_size = file_len - offset;

    // if(my_rank < 100)
    // {
    //   printf("file name is %s offset %ld %ld %ld\n", file_name, *rec_len,  offset, cur_block_size);
    // }


}

void c_fread_( void *data, long* data_size)
{
    size_t data_read = *data_size;
    size_t r_size = fread(data, 1, data_read, f_ptr);


    // if(data_read > cur_block_size) printf("READ ERROR: data read  %ld > block size %d\n", data_read, cur_block_size);


    // if(my_rank == 0) printf("data size is %ld\n", data_read);
    if(r_size != data_read) printf("read error %s %ld %ld\n", file_name, data_read, r_size);

}


void c_seek_fread_( void *data, long* rec_offset,  long* data_size)
{
    size_t data_read = *data_size;


    int block_id = 0;
    if(remap_ocn_blocks == NULL) 
      block_id = my_rank;
    else
      block_id = remap_ocn_blocks[my_rank];


    int file_id = block_id / 100;
    long file_offset = block_id % 100;


    long offset = file_offset * (*data_size) + (*rec_offset);

    // long offset = *rec_offset;
    int s_ret = fseek(f_ptr, offset, SEEK_SET);
    if(s_ret != 0) printf("ERROR: file %s fseek error %ld\n", file_name, offset);

    size_t r_size = fread(data, 1, data_read, f_ptr);


    // if(data_read > cur_block_size) printf("READ ERROR: data read  %ld > block size %d\n", data_read, cur_block_size);


    // if(my_rank == 0) printf("data size is %ld\n", data_read);
    if(r_size != data_read) printf("read error %s %ld %ld\n", file_name, data_read, r_size);

}


void c_fclose_()
{
  fclose(f_ptr);
}
