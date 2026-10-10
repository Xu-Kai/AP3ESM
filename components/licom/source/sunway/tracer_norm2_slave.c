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


typedef struct tracer_norm2_Para{
    int     cuGrid[3];
    int     cuBlock[3];
    double  *d_err_norm2;
    double  *d_tarea;
    double  *d_net;
    double  *d_vit;
}norm2Para;

__thread_local int blockIndex;


void tracer_norm2_kernel_(norm2Para *Paras){

    //ATHREAD basic operation
    int blockId_x,blockId_y,blockId_z;
    int i,j,k;
    blockIndex = CRTS_tid;

    //cuda grid and block affine translation
    int blockSum   = Paras->cuGrid[0] * Paras->cuGrid[1] * Paras->cuGrid[2];
    int quotient   = blockSum / 64;
    int remainder  = blockSum % 64;
    int blockStart = quotient * blockIndex;
    if(remainder > 0)
        blockStart += (remainder > blockIndex ? blockIndex : remainder);
    int blockEnd = blockStart + quotient;
    if(remainder > 0)
        blockEnd   += (remainder > blockIndex ? 1 : 0);

    double err_norm2 = 0.0;
    for(int block = blockStart; block < blockEnd; block++){
        blockId_z = block / (Paras->cuGrid[0] * Paras->cuGrid[1]);
        blockId_y = (block - blockId_z * (Paras->cuGrid[0] * Paras->cuGrid[1])) / Paras->cuGrid[0];
        blockId_x = block % Paras->cuGrid[0];
        j = blockId_y * Paras->cuBlock[1];

//        CRTS_dma_get(d_vtl, &Paras->d_vtl[k * imt * jmt + j * imt], sizeof(double) * imt);
//        CRTS_dma_get(d_vtl_ori, &Paras->d_vtl_ori[k * imt * jmt + j * imt], sizeof(double) * imt);
        for (i = 0; i < imt; i++) {
            if(2<=j&&j<(jmt-2)&&2<=i&&i<(imt-2))
                err_norm2 = err_norm2 + Paras->d_tarea[j*imt+i]
                        *Paras->d_net[imt*jmt + j*imt + i] * Paras->d_vit[j*imt+i];
        }

    //    CRTS_dma_put(&Paras->d_vtl_ori[k * imt * jmt + j * imt], d_vtl_ori, sizeof(double) * imt);
    }
    Paras->d_err_norm2[0] = err_norm2;
    return;
}