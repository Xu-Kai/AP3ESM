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
    bool   d_iau;
    double d_afc1;
    double d_afc2;
    double d_lamda;
    double *d_viv;
    double *d_dzp;
    double *d_u;
    double *d_v;
    double *d_up;
    double *d_ub;
    double *d_wka;
    double *d_ohbu;
    double *d_vtf;
    double *d_utf;
    double *d_work;
    double *d_ucheck;
    double *d_vcheck;
    double *vv_iau;
    double *uu_iau;
}Para;

//dynamic atomic definition
#define CP_CORE      64
#define RCV_CORE     0
#define DATAPAKAGE   1
//#define BC5DEBUG

//ATHREAD global definition 
__thread_local int warp_size = 32; 
__thread_local int blockIndex;
__thread_local Para* Para_s;

void bclinc5_cu_slave_kernel_(Para *Paras){

    double d_wka[imt];
    double d_wka1[imt];
    double d_utf[imt];
    double d_vtf[imt];
    double d_up[imt];
    double d_u[imt];
    double d_v[imt];
    double *d_viv;
    double wka_s;
	//ATHREAD basic operation
    int blockId_x,blockId_y,blockId_z;
    int i,j,k;
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

        CRTS_dma_get(d_wka, &Para_s->d_wka[k * imt * jmt + j * imt], sizeof(double) * imt);
        d_viv = &Para_s->d_viv[k * imt * jmt + j * imt];

        for (i = 0; i < imt; i++) {
            d_wka1[i] = (d_wka[i] - Para_s->d_work[j * imt + i] + Para_s->d_ub[j * imt + i]) * d_viv[i];
        }

        CRTS_dma_put(&Para_s->d_wka[k * imt * jmt + j * imt], d_wka1, sizeof(double) * imt);

    }


    ldm_free(Para_s, sizeof(Para));
    return;
}


void bclinc5_1_cu_slave_kernel_(Para *Paras){

    double *d_wka;
    double *d_utf;
    double *d_vtf;
    double *d_up;
    double *d_u;
    double *d_v;
    double *d_viv;
    double wka_s;
    //ATHREAD basic operation
    int blockId_x,blockId_y,blockId_z;
    int i,j,k;
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
        d_up  = &Para_s->d_up[k * imt * jmt + j * imt];
        d_u   = &Para_s->d_u[k * imt * jmt + j * imt];
        d_v   = &Para_s->d_v[k * imt * jmt + j * imt];
//        CRTS_dma_get(d_utf, &Para_s->d_utf[k * imt * jmt + j * imt], sizeof(double) * imt);
//        CRTS_dma_get(d_vtf, &Para_s->d_vtf[k * imt * jmt + j * imt], sizeof(double) * imt);

        for (i = 0; i < imt; i++) {
            d_up[i] = Para_s->d_afc2 * d_u[i] + Para_s->d_afc1 * (d_up[i] + d_wka[i]);
            d_u[i] = d_wka[i];
//            d_utf[i] += d_u[i];
//            d_vtf[i] += d_v[i];
        }

//        CRTS_dma_put(&Para_s->d_up[k * imt * jmt + j * imt], d_up, sizeof(double) * imt);
//        CRTS_dma_put(&Para_s->d_u[k * imt * jmt + j * imt], d_u, sizeof(double) * imt);
//        CRTS_dma_put(&Para_s->d_utf[k * imt * jmt + j * imt], d_utf, sizeof(double) * imt);
//        CRTS_dma_put(&Para_s->d_vtf[k * imt * jmt + j * imt], d_vtf, sizeof(double) * imt);

    }


    ldm_free(Para_s, sizeof(Para));
    return;
}


void bclinc5_2_cu_slave_kernel_(Para *Paras){

    double *d_utf;
    double *d_vtf;
    double *d_u;
    double *d_v;
    double *d_ucheck;
    double *d_vcheck;
    double wka_s;
    //ATHREAD basic operation
    int blockId_x,blockId_y,blockId_z;
    int i,j,k;
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

        d_u      = &Para_s->d_u[k * imt * jmt + j * imt];
        d_v      = &Para_s->d_v[k * imt * jmt + j * imt];
        d_utf    = &Para_s->d_utf[k * imt * jmt + j * imt];
        d_vtf    = &Para_s->d_vtf[k * imt * jmt + j * imt];
        d_ucheck = &Para_s->d_ucheck[k * imt * jmt + j * imt];
        d_vcheck = &Para_s->d_vcheck[k * imt * jmt + j * imt];

        for (i = 0; i < imt; i++) {
            if(Para_s->d_iau){
                d_utf[i] = d_utf[i] + d_u[i] + Para_s->d_viv[k*imt * jmt + j * imt + i] *
                    Para_s->uu_iau[k * imt * jmt + j * imt + i] * Para_s->d_lamda ;
                d_vtf[i] = d_vtf[i] + d_v[i] - Para_s->d_viv[k*imt * jmt + j * imt + i] *
                    Para_s->vv_iau[k * imt * jmt + j * imt + i] * Para_s->d_lamda ;
            }
            else{
                d_utf[i] += d_u[i];
                d_vtf[i] += d_v[i];
            }
            d_ucheck[i] = d_u[i];
            d_vcheck[i] = d_v[i];
        }

//        CRTS_dma_put(&Para_s->d_utf[k * imt * jmt + j * imt], d_utf, sizeof(double) * imt);
//        CRTS_dma_put(&Para_s->d_vtf[k * imt * jmt + j * imt], d_vtf, sizeof(double) * imt);
//        CRTS_dma_put(&Para_s->d_ucheck[k * imt * jmt + j * imt], d_ucheck, sizeof(double) * imt);
//        CRTS_dma_put(&Para_s->d_vcheck[k * imt * jmt + j * imt], d_vcheck, sizeof(double) * imt);
    }


    ldm_free(Para_s, sizeof(Para));
    return;
}