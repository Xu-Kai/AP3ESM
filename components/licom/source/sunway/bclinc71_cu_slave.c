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
    int cuBlock[3];
    int cuGrid[3];
    int rank;
    int padd;
    double d_dtc;
    int *d_kmu;
    double *d_viv;
    double *d_odzt;
    double *d_odzp;
    double *d_akmu;
    double *d_u;
    double *d_up;
    double *d_v;
    double *d_vp;
    double *d_dlu;
    double *d_dlv;
    double *d_sbcx;
    double *d_bbcx;
}Para;

//dynamic atomic definition
#define CP_CORE      64
#define RCV_CORE     0
#define DATAPAKAGE   1

//ATHREAD global definition 
__thread_local int warp_size = 32; 
__thread_local int blockIndex;
__thread_local Para* Para_s;


void bclinc71_cu_slave_kernel_(Para *Paras){

    int blockId_x,blockId_y,blockId_z;
    int i,j,k;
    double d_u[imt];
    double d_v[imt];
    double *d_dlu;
    double *d_dlv;
    double *d_up;
    double *d_vp;
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

        d_dlu = &Para_s->d_dlu[k * imt * jmt + j * imt];
        d_dlv = &Para_s->d_dlv[k * imt * jmt + j * imt];
        d_up  = &Para_s->d_up[k * imt * jmt + j * imt];
        d_vp  = &Para_s->d_vp[k * imt * jmt + j * imt];

        for (i = 0; i < imt; i++) {
            if (i >= 2 && i < imt - 2 && j >= 2 && j < jmt - 2 && k < km) {
                d_u[i] = d_up[i] + d_dlu[i] * Para_s->d_dtc;
                d_v[i] = d_vp[i] + d_dlv[i] * Para_s->d_dtc;
            }
        }

        CRTS_dma_put(&Para_s->d_u[k * imt * jmt + j * imt], d_u, sizeof(double) * imt);
        CRTS_dma_put(&Para_s->d_v[k * imt * jmt + j * imt], d_v, sizeof(double) * imt);
    }

    ldm_free(Para_s, sizeof(Para));
    return;
}
