
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



typedef struct pop_halo_cp_kernelPara{
    int data_len;
    double fill;
    double *d_data;
}pop_halo_cp_para;


//ATHREAD DMA transfer definition
__thread_local unsigned int D_COUNT = 0;
__thread_local crts_rply_t dma_rply = 0;
//ATHREAD global definition
__thread_local int blockIndex;

void  pop_halo_cp_(pop_halo_cp_para *Paras) {

    int i, j, k;
    blockIndex = CRTS_tid;
    int data_len = Paras->data_len;
    int numPerCpe = data_len / 64;
    int remainder = data_len % 64;

    int start = blockIndex * numPerCpe;

    int sum;
    double fill=Paras->fill;

    if (blockIndex < 63){
        sum  = numPerCpe;
    }
    else{
        if(remainder == 0){
            sum  = numPerCpe;
	}
        else{
            sum  = remainder;
	}
    }

    for(int i=0; i<sum; i++) {
        Paras->d_data[start + i] = fill;
    }
}
