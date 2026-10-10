//ATHREAD include file
#include <crts.h>
#include <slave.h>
#include <simd.h>
#include <math.h>
#include <stdio.h>
#include <float.h>
#include <string.h>
#include <stdbool.h>
#include "swperf.h"
#include "param_mod.h"

//ATHREAD DMA transfer definition 
__thread_local unsigned int D_COUNT = 0;
__thread_local crts_rply_t dma_rply = 0;
__thread_local crts_rply_t l_rply = 0;
__thread_local crts_rply_t r_rply = 0;

typedef struct kernelPara{
    int    cuBlock[3];
    int    cuGrid[3];
    int    rank;
    int    pad;
    double d_dtc2;
    double *d_vp;
    double *d_dlv;
    double *d_wka;
}Para;

//dynamic atomic definition
#define CP_CORE      64
#define RCV_CORE     0
#define DATAPAKAGE   1

//ATHREAD global definition 
__thread_local int warp_size = 32; 
__thread_local int blockIndex;
__thread_local Para* Para_s;


void bclinc2_cu_slave_kernel_(Para *Paras){
    int blockId_x,blockId_y,blockId_z;
    int i,j,k;
    double *d_wka;
    double *d_vp;
    double *d_dlv;
	//ATHREAD basic operation
    blockIndex = CRTS_tid;
    Para_s = (Para*)ldm_malloc(sizeof(Para));
    CRTS_dma_iget(Para_s, Paras, sizeof(Para), &dma_rply);
    D_COUNT++;
    CRTS_dma_wait_value(&dma_rply, D_COUNT);

    //cuda grid and block affine translation
    int blockSum   = Para_s->cuGrid[0] * Para_s->cuGrid[1] * Para_s->cuGrid[2];
    int threadSum  = Para_s->cuBlock[0] * Para_s->cuBlock[1] * Para_s->cuBlock[2];
    int quotient   = blockSum / 64;
    int remainder  = blockSum % 64;
    int blockStart = quotient * blockIndex;
    if(remainder > 0)
        blockStart += (remainder > blockIndex ? blockIndex : remainder);
    int blockEnd = blockStart + quotient;
    if(remainder > 0)
        blockEnd   += (remainder > blockIndex ? 1 : 0);

    for(int block = blockStart; block < blockEnd; block++){
        blockId_z = block / (Para_s->cuGrid[0] * Para_s->cuGrid[1]);
        blockId_y = (block - blockId_z * (Para_s->cuGrid[0] * Para_s->cuGrid[1])) / Para_s->cuGrid[0];
        blockId_x = block % Para_s->cuGrid[0];
        j = blockId_y * Para_s->cuBlock[1];
        k = blockId_z * Para_s->cuBlock[2];

        d_wka = &Para_s->d_wka[k * imt * jmt + j * imt];
        d_vp  = &Para_s->d_vp[k * imt * jmt + j * imt];
        d_dlv = &Para_s->d_dlv[k * imt * jmt + j * imt];

        for (i = 0; i < imt; i++) {
            if(i >= 2 && i < imt - 2 && j >= 2 && j < jmt - 2 && k < km){
                d_wka[i] = d_vp[i] + d_dlv[i] * Para_s->d_dtc2;
            }
        }

//        CRTS_dma_put(&Para_s->d_wka[k * imt * jmt + j * imt], d_wka, sizeof(double) * imt);

    }

    ldm_free(Para_s, sizeof(Para));
    return;
}
