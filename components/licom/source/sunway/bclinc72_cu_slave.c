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
    double d_dtc;
    int    *d_kmu;
    double *d_viv;
    double *d_odzt;
    double *d_odzp;
    double *d_akmu;
    double *d_v;
    double *d_vp; 
    double *d_bbcy;
    double *d_dlv;
    double *d_sbcy; 
    double *d_a8;
    double *d_b8;
    double *d_c8;
    double *d_d8;
    double *d_e8;
    double *d_f8;
}Para;

//dynamic atomic definition
#define CP_CORE      64
#define RCV_CORE     0
#define DATAPAKAGE   1

//ATHREAD global definition 
__thread_local int warp_size = 32; 
__thread_local int blockIndex;
__thread_local Para* Para_s;


void bclinc72_cu_slave_kernel_(Para *Paras){
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

    double d_a8[km],d_b8[km],d_c8[km],d_d8[km],d_e8[km+1],d_f8[km+1];

    for(int block = blockStart; block < blockEnd; block++) {
        int blockId_z = block / (Para_s->cuGrid[0] * Para_s->cuGrid[1]);
        int blockId_y = (block - blockId_z * (Para_s->cuGrid[0] * Para_s->cuGrid[1])) / Para_s->cuGrid[0];
        int blockId_x = block % Para_s->cuGrid[0];
        // ATHREAD serial operation replace GPU parallel mutli-thread in single thread-block
        for(int threadId = 0; threadId < threadSum; threadId++) {
            int threadId_x = threadId % Para_s->cuBlock[0];
            int threadId_z = threadId / (Para_s->cuBlock[0] * Para_s->cuBlock[1]);
            int threadId_y = (threadId - threadId_z * (Para_s->cuBlock[0] * Para_s->cuBlock[1])) / Para_s->cuBlock[0];
            

            int i, j, k;

            i = (blockId_x) * Para_s->cuBlock[0] + threadId_x;
            j = (blockId_y) * Para_s->cuBlock[1] + threadId_y;


            if (i >= 2 && i < imt - 2 && j >= 2 && j < jmt - 2) {
                double g0;
                int kz = Para_s->d_kmu[j * imt + i];

                for (k = 0; k < km; k++) {
                    Para_s->d_v[k * jmt * imt + j * imt + i] =
                            Para_s->d_vp[k * jmt * imt + j * imt + i] +
                            Para_s->d_dlv[k * jmt * imt + j * imt + i] * Para_s->d_dtc;
                }

                if (kz > 0) {
                    for (k = 1; k < kz; k++) {
                        d_a8[k] =
                                Para_s->d_akmu[(k - 1) * jmt * imt + j * imt + i] *
                                Para_s->d_odzt[k] * Para_s->d_odzp[k] * Para_s->d_dtc * 0.5;
                        d_d8[k] =
                                Para_s->d_v[k * jmt * imt + j * imt + i];
                    }
                    for (k = 1; k < kz - 1; k++) {
                        d_c8[k] =
                                Para_s->d_akmu[k * jmt * imt + j * imt + i] *
                                Para_s->d_odzt[k + 1] * Para_s->d_odzp[k] * Para_s->d_dtc * 0.5;
                        d_b8[k] =
                                1.0 + d_a8[k] +
                                d_c8[k];
                        d_e8[k] = 0.0;
                        d_f8[k] = 0.0;
                    }

                    //  B. C. AT TOP
                    k = 0;
                    d_a8[k] = Para_s->d_odzp[k] * Para_s->d_dtc * 0.5;
                    d_c8[k] =
                            Para_s->d_akmu[k * jmt * imt + j * imt + i] *
                            Para_s->d_odzt[k + 1] * Para_s->d_odzp[k] * Para_s->d_dtc * 0.5;
                    d_b8[k] = 1.0 + d_c8[k];
                    d_d8[k] = Para_s->d_v[k * jmt * imt + j * imt + i];
                    d_e8[k] = 0.0;
                    d_f8[k] = 0.0;

                    //B. C. AT BOTTOM
                    d_b8[(kz - 1)] = 1.0 + d_a8[(kz - 1)];
                    d_c8[(kz - 1)] = Para_s->d_odzp[kz - 1] * Para_s->d_dtc * 0.5;
                    d_e8[kz] = 0.0;
                    d_f8[kz] = 0.0;
                    d_d8[(kz - 1)] =
                            Para_s->d_v[(kz - 1) * jmt * imt + j * imt + i] -
                            Para_s->d_bbcy[j * imt + i] * Para_s->d_odzp[kz - 1] * Para_s->d_dtc * 0.5;

                    //NOW INVERT
                    for (k = kz - 1; k >= 0; k--) {
                        g0 = 1.0 / (d_b8[k] -
                                    d_c8[k] *
                                    d_e8[(k + 1)]);
                        d_e8[k] =
                                d_a8[k] * g0;
                        d_f8[k] =
                                (d_d8[k] +
                                d_c8[k] *
                                d_f8[(k + 1)]) * g0;
                    }

                    //B.C. AT SURFACE
                    Para_s->d_v[j * imt + i] =
                            (d_e8[0] *
                            Para_s->d_sbcy[j * imt + i] +
                            d_f8[0]) *
                            Para_s->d_viv[j * imt + i];
                    for (k = 1; k < kz; k++) {
                        Para_s->d_v[k * jmt * imt + j * imt + i] =
                                (d_e8[k] *
                                Para_s->d_v[(k - 1) * jmt * imt + j * imt + i] +
                                d_f8[k]) *
                                Para_s->d_viv[k * jmt * imt + j * imt + i];
                    }
                }
            }
        }
        }
        ldm_free(Para_s, sizeof(Para));
        return;
}
