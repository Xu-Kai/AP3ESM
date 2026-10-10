
//ATHREAD include file
#include <crts.h>
#include <slave.h>
//#include <simd.h>
#include <math.h>
#include <stdio.h>
#include <float.h>
#include <string.h>
#include <stdbool.h>
#include "swperf.h" 

#include "param_mod.h"	



typedef struct readyt_cp_kernelPara{
    double *d_dlu;
    double *d_dlv;
    double *d_rit;
    double *d_ric;
    double *d_rict;
    double *d_ricdtTms;
    double *d_ricdt;
    double *d_akt;
    int rank;
}readyt_cp_para;


//ATHREAD DMA transfer definition
__thread_local unsigned int D_COUNT = 0;
__thread_local crts_rply_t dma_rply = 0;
//ATHREAD global definition
__thread_local int blockIndex;

void  readyt_cp_(readyt_cp_para *Paras) {

    int i, j, k;
    blockIndex = CRTS_tid;
    int totalnum = km * jmt * imt;
    int numPerCpe = totalnum / 64;
    int remainder = totalnum % 64;

    int start = blockIndex * numPerCpe;
    
    int totalnum1 = kmm1 * jmt * imt;
    int numPerCpe1 = totalnum1 / 64;
    int remainder1 = totalnum1 % 64;
    int start1 = blockIndex * numPerCpe1;
 
    int sum;
    int sum1;
    if (blockIndex < 63){
        sum  = numPerCpe;
        sum1 = numPerCpe1;
    }
    else{
        if(remainder == 0){
            sum  = numPerCpe;
            sum1 = numPerCpe1;
	}
        else{
            sum  = remainder;
            sum1 = remainder1;
	}
    }

    for(int i=0; i<sum; i++) {
        Paras->d_dlu[start + i] = 0.0;
        Paras->d_dlv[start + i] = 0.0;
    }
    for(int i=0; i<sum1; i++) {
        Paras->d_rit[start1 + i] = 0.0;
        Paras->d_ric[start1 + i] = 0.0;
        Paras->d_rict[start1 + i] = 0.0;
        Paras->d_ricdtTms[start1 + i] = 0.0;
        Paras->d_ricdt[start1 + i] = 0.0;
    }
}


void  readyt_cp_1_(readyt_cp_para *Paras){

    int i,j,k;
    blockIndex = CRTS_tid;
    int totalnum = 2* km * jmt * imt;

    int avaiLdm = get_allocatable_size();
    //int numPerCpe = avaiLdm / (64 * sizeof(double));
    int numPerCpe = totalnum / 64;
    int remainder = totalnum % numPerCpe;

    int start = blockIndex * numPerCpe;
    int sum;
    if (blockIndex < 63)
        sum = numPerCpe;
    else
        if(remainder == 0)
            sum = numPerCpe;
        else
            sum = remainder;

//    if(Paras->rank == 0)
//        printf("[slave %d]:numPerCpe:%d sum:%d,start:%d\n",blockIndex,numPerCpe,sum,start);
    for(int i=0; i<sum; i++)
        Paras->d_akt[start+i] = 0.0;
}
