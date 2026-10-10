
#include "def-undef.h"
//#define ENABLE_GPLT
#ifdef ENABLE_GPLT
#include <stdlib.h>
#include <mpi.h>
#include "gptl.h"
extern void gptl_initialize_();
extern void gptl_finalize_(const MPI_Comm comm);
extern void t_startf_(char *name);
extern void t_stopf_(char *name);


int gptl_ret;
void gptl_initialize_(){
    gptl_ret = GPTLsetoption(GPTLpercent, 1);
    gptl_ret = GPTLsetoption(GPTLoverhead, 1);
//    gptl_ret = GPTLsetoption(GPTLverbose, 1);    
    gptl_ret = GPTLinitialize();
    //gptl_ret = GPTLstart ("LICOM3-Kokkos");
}

void gptl_finalize_(const MPI_Comm comm){
    GPTLpr_summary_file(MPI_COMM_WORLD, "./gptl.log");
    gptl_ret = GPTLfinalize();
}

void t_startf_(char *name){
    GPTLstart (name);
}

void t_stopf_(char *name){
    GPTLstop (name);
}

#endif
