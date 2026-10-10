#include <stdlib.h>
#include <stdio.h>
#include <libgen.h>
#include <mpi.h>
#include "def-undef.h"
static FILE *f_ptr;
static long c_my_task;
static int file_block_num = 100;
static int my_rank;
static int w_size;

static double *data_buffer;
static double *data_buffer_kmp1;
static float *data_buffer_r4;
static float *data_buffer_kmp1_r4;
static MPI_Comm write_ocn_comm;
static MPI_Comm file_group;
static int file_rank;
static int is_init = 0;
static double *other_buffer;
static double *coup_buffer;
static long x_dim, y_dim,  k_dim;
static long local_block_num;

void write_init_(MPI_Fint *comm_fortran,long *nx_block, long *ny_block, long *km)
{
#ifdef COUP
    write_ocn_comm = MPI_Comm_f2c(*comm_fortran);
#else
	write_ocn_comm = MPI_COMM_WORLD;
#endif
  MPI_Comm_rank(write_ocn_comm, &my_rank);
  MPI_Comm_split(write_ocn_comm, my_rank/file_block_num, my_rank%file_block_num, &file_group);
  MPI_Comm_size(write_ocn_comm, &w_size);



  MPI_Comm_rank(file_group, &file_rank);
  x_dim = (*nx_block);
  y_dim = (*ny_block);
  k_dim = (*km);
  long block_start = (my_rank/file_block_num) * file_block_num;
  long block_end = (my_rank/file_block_num + 1) * file_block_num;
  local_block_num = block_end > w_size ? w_size - block_start : file_block_num;
  if(is_init == 0)
  {

      if(file_rank == 0)
      {
        data_buffer = (double*)malloc(x_dim * y_dim  * k_dim * file_block_num * sizeof(double));
        data_buffer_kmp1 = (double*)malloc(x_dim * y_dim  * (k_dim+1) * file_block_num * sizeof(double));
        data_buffer_r4 = (float*)malloc(x_dim * y_dim  * k_dim * file_block_num * sizeof(float));
        data_buffer_kmp1_r4 = (float*)malloc(x_dim * y_dim  * (k_dim+1) * file_block_num * sizeof(float));
        if(data_buffer == NULL)
        {
          printf("ERROR: can't allocate memory data buffer\n");
          fflush(stdout);
          MPI_Abort(write_ocn_comm, MPI_ERR_UNKNOWN);
        }
      }

      other_buffer = (double*)malloc(x_dim * y_dim  * 12 * sizeof(double));
      if(other_buffer == NULL) 
      {
        printf("ERROR: can't allocate memory other buffer\n");
        fflush(stdout);
        MPI_Abort(write_ocn_comm, MPI_ERR_UNKNOWN);
      }

      coup_buffer = (double*)malloc(x_dim * y_dim  * 7 * sizeof(double));
      if(coup_buffer == NULL)
      {
          printf("ERROR: can't allocate memory other buffer\n");
          fflush(stdout);
          MPI_Abort(write_ocn_comm, MPI_ERR_UNKNOWN);
      }


      is_init = 1;
      MPI_Barrier(write_ocn_comm);
  }
};


void ocn_create_directory(char* parten_dir, char* file_prefix)
{
    MPI_Comm_rank(write_ocn_comm, &my_rank);
    // if(my_rank % 100 == 0) printf("in create dir %d\n", my_rank);
    if(my_rank == 0)
    {
        char dir_name[256];
//        printf("parten_dir %s %s\n", parten_dir, file_prefix);
        sprintf(dir_name, "%s/%s", "./", parten_dir);
//        printf("parrent dir_name dir %s\n", dir_name);
      mkdir(dir_name, 0777);
            // sprintf(dir_name,  "%s/%s.dir", parten_dir, file_prefix); 

        // if (mkdir(parten_dir, 0777) != -1) 
        //   printf("create dir %s\n", parten_dir);
        // else
        //   printf("dir %s exists\n", parten_dir);

      sprintf(dir_name,  "%s/%s.dir", dir_name, file_prefix);
      mkdir(dir_name, 0777);
      //   if (mkdir(dir_name, 0777) != -1) 
      //     printf("create dir %s\n", dir_name);
      //   else
      //     printf("dir %s exists\n", dir_name);
      // printf("parent dir %s\n", parten_dir);
      printf("sub dir_name dir %s\n", dir_name);

    }
    MPI_Barrier(write_ocn_comm);
    // if(file_rank == 0)
    // {

    //     // Creating a directory
    //   if (mkdir(parten_dir, 0777) != -1)
    //       printf("create dir %s\n", parten_dir);
    //   else
    //       printf("dir %s exists\n", parten_dir);



    //   if (mkdir(dir_name, 0777) != -1)
    //       printf("create dir %s\n", dir_name);
    //   else
    //       printf("dir %s exists\n", dir_name);
    // }
    // MPI_Barrier(file_group);  
}

void  write_km_data_(double *data, char* dir_name, char* file_name)
{



  ocn_create_directory(dir_name, file_name);


  int send_count = x_dim * y_dim * k_dim;
  int recv_count = x_dim * y_dim * k_dim;

  MPI_Gather(data, send_count, MPI_DOUBLE, data_buffer, recv_count, MPI_DOUBLE, 0, file_group);


   FILE* fptr;
  if(file_rank == 0)
  {
    // printf("my rank is %d %d\n", my_rank, file_rank);
    //       fflush(stdout);
    int gorup_rank = my_rank / file_block_num;
    char data_file_name[256];
    sprintf(data_file_name,  "%s/%s.dir/%s_%04d", dir_name, file_name, file_name, gorup_rank);
    // printf("data_file_name %s\n", data_file_name);

    fptr = fopen(data_file_name, "w");
    if(fptr == NULL) 
    {
      printf("ERROR: can't open file %s\n", data_file_name);
      fflush(stdout);

      MPI_Abort(write_ocn_comm, MPI_ERR_UNKNOWN);
    }
    long data_size = x_dim * y_dim * k_dim * sizeof(double) * local_block_num;
    long write_size = fwrite(data_buffer, 1, data_size, fptr);
    if(write_size != data_size) 
    {
      printf("ERROR: km write size is %ld %ld", data_size, write_size);
      fflush(stdout);
      MPI_Abort(write_ocn_comm, MPI_ERR_UNKNOWN);

    }
    fclose(fptr);

  }
  MPI_Barrier(write_ocn_comm);


}

void  write_kmp1_data_(double *data, char* dir_name, char* file_name)
{



    ocn_create_directory(dir_name, file_name);


    int send_count = x_dim * y_dim * (k_dim+1);
    int recv_count = x_dim * y_dim * (k_dim+1);

    MPI_Gather(data, send_count, MPI_DOUBLE, data_buffer_kmp1, recv_count, MPI_DOUBLE, 0, file_group);


    FILE* fptr;
    if(file_rank == 0)
    {
        // printf("my rank is %d %d\n", my_rank, file_rank);
        //       fflush(stdout);
        int gorup_rank = my_rank / file_block_num;
        char data_file_name[256];
        sprintf(data_file_name,  "%s/%s.dir/%s_%04d", dir_name, file_name, file_name, gorup_rank);
        // printf("data_file_name %s\n", data_file_name);

        fptr = fopen(data_file_name, "w");
        if(fptr == NULL)
        {
            printf("ERROR: can't open file %s\n", data_file_name);
            fflush(stdout);

            MPI_Abort(write_ocn_comm, MPI_ERR_UNKNOWN);
        }
        long data_size = x_dim * y_dim * (k_dim+1) * sizeof(double) * local_block_num;
        long write_size = fwrite(data_buffer_kmp1, 1, data_size, fptr);
        if(write_size != data_size)
        {
            printf("ERROR: km write size is %ld %ld", data_size, write_size);
            fflush(stdout);
            MPI_Abort(write_ocn_comm, MPI_ERR_UNKNOWN);

        }
        fclose(fptr);

    }
    MPI_Barrier(write_ocn_comm);


}

void  write_data_(double *data, char* dir_name, char* file_name, long *k_n)
{



  ocn_create_directory(dir_name, file_name);

  if((*k_n) > k_dim) 
  {
    printf("ERROR: k > NKM k is %ld,  NKM is %ld", *k_n, k_dim);
  }
  int send_count = x_dim * y_dim *(*k_n);
  int recv_count = x_dim * y_dim *(*k_n);

  MPI_Gather(data, send_count, MPI_DOUBLE, data_buffer, recv_count, MPI_DOUBLE, 0, file_group);


   FILE* fptr;
  if(file_rank == 0)
  {
    // printf("my rank is %d %d\n", my_rank, file_rank);
    //       fflush(stdout);
    int gorup_rank = my_rank / file_block_num;
    char data_file_name[256];
    sprintf(data_file_name,  "%s/%s.dir/%s_%04d", dir_name, file_name, file_name, gorup_rank);
    // printf("data_file_name %s\n", data_file_name);

    fptr = fopen(data_file_name, "w");
    if(fptr == NULL) 
    {
      printf("ERROR: can't open file %s\n", data_file_name);
      fflush(stdout);

      MPI_Abort(write_ocn_comm, MPI_ERR_UNKNOWN);
    }
    long data_size = x_dim * y_dim * (*k_n) * sizeof(double) * local_block_num;
    long write_size = fwrite(data_buffer, 1, data_size, fptr);
    if(write_size != data_size) 
    {
      printf("ERROR: km write size is %ld %ld", data_size, write_size);
      fflush(stdout);
      MPI_Abort(write_ocn_comm, MPI_ERR_UNKNOWN);

    }
    fclose(fptr);

  }
  MPI_Barrier(write_ocn_comm);


}

void  write_data_r4_(float *data, char* dir_name, char* file_name, long *k_n)
{



    ocn_create_directory(dir_name, file_name);

    if((*k_n) > k_dim)
    {
        printf("ERROR: k > NKM k is %ld,  NKM is %ld", *k_n, k_dim);
    }
    int send_count = x_dim * y_dim *(*k_n);
    int recv_count = x_dim * y_dim *(*k_n);

    MPI_Gather(data, send_count, MPI_REAL, data_buffer_r4, recv_count, MPI_REAL, 0, file_group);


    FILE* fptr;
    if(file_rank == 0)
    {
        // printf("my rank is %d %d\n", my_rank, file_rank);
        //       fflush(stdout);
        int gorup_rank = my_rank / file_block_num;
        char data_file_name[256];
        sprintf(data_file_name,  "%s/%s.dir/%s_%04d", dir_name, file_name, file_name, gorup_rank);
//        printf("data_file_name %s, k_n : %ld\n", data_file_name, (*k_n));

        fptr = fopen(data_file_name, "w");
        if(fptr == NULL)
        {
            printf("ERROR: can't open file %s\n", data_file_name);
            fflush(stdout);

            MPI_Abort(write_ocn_comm, MPI_ERR_UNKNOWN);
        }
        long data_size = x_dim * y_dim * (*k_n) * sizeof(float) * local_block_num;
        long write_size = fwrite(data_buffer_r4, 1, data_size, fptr);
        if(write_size != data_size)
        {
            printf("ERROR: km write size is %ld %ld", data_size, write_size);
            fflush(stdout);
            MPI_Abort(write_ocn_comm, MPI_ERR_UNKNOWN);

        }
        fclose(fptr);

    }
    MPI_Barrier(write_ocn_comm);


}

void  write_data_kmp1_r4_(float *data, char* dir_name, char* file_name, long *k_n)
{



    ocn_create_directory(dir_name, file_name);

    if((*k_n) > k_dim+1)
    {
        printf("ERROR: k > NKM k is %ld,  NKM is %ld", *k_n, k_dim);
    }
    int send_count = x_dim * y_dim *(*k_n);
    int recv_count = x_dim * y_dim *(*k_n);

    MPI_Gather(data, send_count, MPI_REAL, data_buffer_kmp1_r4, recv_count, MPI_REAL, 0, file_group);


    FILE* fptr;
    if(file_rank == 0)
    {
        // printf("my rank is %d %d\n", my_rank, file_rank);
        //       fflush(stdout);
        int gorup_rank = my_rank / file_block_num;
        char data_file_name[256];
        sprintf(data_file_name,  "%s/%s.dir/%s_%04d", dir_name, file_name, file_name, gorup_rank);
//        printf("data_file_name %s, k_n : %ld\n", data_file_name, (*k_n));

        fptr = fopen(data_file_name, "w");
        if(fptr == NULL)
        {
            printf("ERROR: can't open file %s\n", data_file_name);
            fflush(stdout);

            MPI_Abort(write_ocn_comm, MPI_ERR_UNKNOWN);
        }
        long data_size = x_dim * y_dim * (*k_n) * sizeof(float) * local_block_num;
        long write_size = fwrite(data_buffer_kmp1_r4, 1, data_size, fptr);
        if(write_size != data_size)
        {
            printf("ERROR: km write size is %ld %ld", data_size, write_size);
            fflush(stdout);
            MPI_Abort(write_ocn_comm, MPI_ERR_UNKNOWN);

        }
        fclose(fptr);

    }
    MPI_Barrier(write_ocn_comm);


}

void  write_other_data_(double *d1, double *d2, double *d3,  double *d4, double *d5, double *d6, double *d7,
                        double *d8, double *d9,double *d10,double *d11,double *d12, char* dir_name, char* file_name)
{

  long d_size = x_dim * y_dim * sizeof(double);
  long d_nmemb = x_dim * y_dim;

  memcpy(other_buffer, d1, d_size);

  MPI_Barrier(write_ocn_comm);


  memcpy(&other_buffer[d_nmemb], d2, d_size);
  memcpy(&other_buffer[2 * d_nmemb], d3, d_size);
  memcpy(&other_buffer[3 * d_nmemb], d4, d_size);
  memcpy(&other_buffer[4 * d_nmemb], d5, d_size);
  memcpy(&other_buffer[5 * d_nmemb], d6, d_size);
  memcpy(&other_buffer[6 * d_nmemb], d7, d_size);
  memcpy(&other_buffer[7 * d_nmemb], d8, d_size);
  memcpy(&other_buffer[8 * d_nmemb], d9, d_size);
  memcpy(&other_buffer[9 * d_nmemb], d10, d_size);
  memcpy(&other_buffer[10 * d_nmemb], d11, d_size);
  memcpy(&other_buffer[11 * d_nmemb], d12, d_size);

  int send_count = x_dim * y_dim * 12;
  int recv_count = x_dim * y_dim * 12;
  MPI_Gather(other_buffer, send_count, MPI_DOUBLE, data_buffer, recv_count, MPI_DOUBLE, 0, file_group);
  ocn_create_directory(dir_name, file_name);

  if(file_rank == 0)
  {
    int gorup_rank = my_rank / file_block_num;
    char data_file_name[256];
    sprintf(data_file_name,  "%s/%s.dir/%s_%04d", dir_name, file_name, file_name, gorup_rank);
    long data_size = x_dim * y_dim *12 * sizeof(double) * local_block_num;
    FILE* fptr = fopen(data_file_name, "w");
    if(fptr == NULL) 
    {
      printf("ERROR: can't open file %s\n", data_file_name);
      fflush(stdout);
      MPI_Abort(write_ocn_comm, MPI_ERR_UNKNOWN);
    }
    long write_size = fwrite(data_buffer, 1, data_size, fptr);
    fclose(fptr);
    if(write_size != data_size) 
    {
      printf("ERROR: write size is %ld %ld", data_size, write_size);
      fflush(stdout);
      MPI_Abort(write_ocn_comm, MPI_ERR_UNKNOWN);
    }
  }
  MPI_Barrier(write_ocn_comm);

}


void  write_coup_data_(double *d1, double *d2, double *d3,  double *d4, double *d5, double *d6, double *d7,
                        char* dir_name, char* file_name)
{

    long d_size = x_dim * y_dim * sizeof(double);
    long d_nmemb = x_dim * y_dim;

    memcpy(coup_buffer, d1, d_size);

    MPI_Barrier(write_ocn_comm);


    memcpy(&coup_buffer[d_nmemb], d2, d_size);
    memcpy(&coup_buffer[2 * d_nmemb], d3, d_size);
    memcpy(&coup_buffer[3 * d_nmemb], d4, d_size);
    memcpy(&coup_buffer[4 * d_nmemb], d5, d_size);
    memcpy(&coup_buffer[5 * d_nmemb], d6, d_size);
    memcpy(&coup_buffer[6 * d_nmemb], d7, d_size);


    int send_count = x_dim * y_dim * 7;
    int recv_count = x_dim * y_dim * 7;
    MPI_Gather(coup_buffer, send_count, MPI_DOUBLE, data_buffer, recv_count, MPI_DOUBLE, 0, file_group);
    ocn_create_directory(dir_name, file_name);

    if(file_rank == 0)
    {
        int gorup_rank = my_rank / file_block_num;
        char data_file_name[256];
        sprintf(data_file_name,  "%s/%s.dir/%s_%04d", dir_name, file_name, file_name, gorup_rank);
        long data_size = x_dim * y_dim *7 * sizeof(double) * local_block_num;
        FILE* fptr = fopen(data_file_name, "w");
        if(fptr == NULL)
        {
            printf("ERROR: can't open file %s\n", data_file_name);
            fflush(stdout);
            MPI_Abort(write_ocn_comm, MPI_ERR_UNKNOWN);
        }
        long write_size = fwrite(data_buffer, 1, data_size, fptr);
        fclose(fptr);
        if(write_size != data_size)
        {
            printf("ERROR: write size is %ld %ld", data_size, write_size);
            fflush(stdout);
            MPI_Abort(write_ocn_comm, MPI_ERR_UNKNOWN);
        }
    }
    MPI_Barrier(write_ocn_comm);

}



void  write_date_data_(long *month, long *day, 
                      char* dir_name, char* file_name)
{
    ocn_create_directory(dir_name, file_name);

  if(file_rank == 0)
  {
    long date_buffer[2];
    date_buffer[0] = (*month);
    date_buffer[0] = (*day);

    int gorup_rank = my_rank / file_block_num;
    char data_file_name[256];
    sprintf(data_file_name,  "%s/%s.dir/%s_%04d", dir_name, file_name, file_name, gorup_rank);

    long data_size = 2 * sizeof(long);
    FILE* fptr = fopen(data_file_name, "w");

    if(fptr == NULL) 
    {
      printf("ERROR: can't open file %s\n", data_file_name);
      fflush(stdout);
      MPI_Abort(write_ocn_comm, MPI_ERR_UNKNOWN);
    }
    int k_id = 0; 
    for(k_id = 0; k_id < local_block_num; k_id ++){
      long write_size = fwrite(date_buffer, 1, data_size, fptr);
      if(write_size != data_size) 
      {
        printf("ERROR: write size is %ld %ld", data_size, write_size);
        fflush(stdout);
        MPI_Abort(write_ocn_comm, MPI_ERR_UNKNOWN);
      }
    }
    fclose(fptr);

  }
  MPI_Barrier(write_ocn_comm);

}
