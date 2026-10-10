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

typedef struct bclinc4Para{
    int    cuBlock[3];
    int    cuGrid[3];
    int    rank;
    int    *d_kmu;
    double d_afc1;
    double d_afc2;
    double d_dtc2;
    double *d_ohbu;
    double *d_viv;
    double *d_dzp;
    double *d_odzt;
    double *d_odzp;
    double *d_akmu;
    double *d_v;
    double *d_up;
    double *d_vp;
    double *d_vb;
    double *d_dlu;
    double *d_sbcx;
    double *d_bbcx;
    double *d_wka;
    double *d_a8;
    double *d_b8;
    double *d_c8;
    double *d_d8;
    double *d_e8;
    double *d_f8;
    double *d_work;
}bclinc4_Para;

typedef struct bclinc4_0_Para{
    int    cuBlock[3];
    int    cuGrid[3];
    double d_afc1;
    double d_afc2;
    double d_dtc2;
    double *d_viv;
    double *d_v;
    double *d_up;
    double *d_vp;
    double *d_vb;
    double *d_dlu;
    double *d_wka;
    double *d_work;
}bclinc4_0_Para;

//dynamic atomic definition
#define CP_CORE      64
#define RCV_CORE     0
#define DATAPAKAGE   1
//#define BC4DEBUG

//ATHREAD global definition 
__thread_local int warp_size = 32; 
__thread_local int blockIndex;


void bclinc4_cu_slave_kernel_(bclinc4_0_Para *Paras){

    double d_wka[imt];
    double d_vp[imt];
    double d_v[imt];
    
    //ATHREAD basic operation
    int blockId_x,blockId_y,blockId_z;
    int i,j,k;
    blockIndex = CRTS_tid;
    bclinc4_0_Para *Para_s = (bclinc4_0_Para*)ldm_malloc(sizeof(bclinc4_0_Para));
    CRTS_dma_iget(Para_s, Paras, sizeof(bclinc4_0_Para), &dma_rply);
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
        CRTS_dma_get(d_vp, &Para_s->d_vp[k * imt * jmt + j * imt], sizeof(double) * imt);

        for ( i = 0; i < imt; i++) {
            if (i < imt && j < jmt){
                //VP (I,J,K,IBLOCK) = AFC2* V (I,J,K,IBLOCK) + AFC1* (VP (I,J,K,IBLOCK) + WKA (I,J,K,IBLOCK))
                d_vp[i] = Para_s->d_afc2 * Para_s->d_v[k * jmt * imt + j * imt + i] +
                          Para_s->d_afc1 * (d_vp[i] + d_wka[i]);
                //V (I,J,K,IBLOCK) = WKA (I,J,K,IBLOCK)
                d_v[i] = d_wka[i];
                if (i >= 2 && i < imt - 2 && j >= 2 && j < jmt - 2)
                    d_wka[i] = Para_s->d_up[k * jmt * imt + j * imt + i] +
                            Para_s->d_dlu[k * jmt * imt + j * imt + i] * Para_s->d_dtc2;
            }
        }

        CRTS_dma_put(&Para_s->d_wka[k * imt * jmt + j * imt], d_wka, sizeof(double) * imt);
        CRTS_dma_put(&Para_s->d_vp[k * imt * jmt + j * imt], d_vp, sizeof(double) * imt);
        CRTS_dma_put(&Para_s->d_v[k * imt * jmt + j * imt], d_v, sizeof(double) * imt);

    }

    ldm_free(Para_s, sizeof(bclinc4_0_Para));
    return;    
}


void bclinc4_cu_0_slave_kernel_(bclinc4_0_Para *Paras){

    double d_wka[imt];
    double d_wka1[imt];
    double d_vp[imt];
    double d_v[imt];

    //ATHREAD basic operation
    int i,j,k;
    blockIndex = CRTS_tid;
    bclinc4_0_Para *Para_s = (bclinc4_0_Para*)ldm_malloc(sizeof(bclinc4_0_Para));
    CRTS_dma_iget(Para_s, Paras, sizeof(bclinc4_0_Para), &dma_rply);
    D_COUNT++;
    CRTS_dma_wait_value(&dma_rply, D_COUNT);

    for( k = 0; k < km; k++) {

        for ( j = blockIndex; j < jmt; j += 64) {

            CRTS_dma_get(d_wka, &Para_s->d_wka[k * imt * jmt + j * imt], sizeof(double) * imt);

            for (i = 0; i < imt; i++) {
                // WKA (I,J,K,IBLOCK)= (WKA (I,J,K,IBLOCK) - WORK (I,J,IBLOCK) + VB (I,J,IBLOCK))* VIV (I,J,K,IBLOCK)
                d_wka1[i] = (d_wka[i] - Para_s->d_work[j * imt + i] + Para_s->d_vb[j * imt + i]) *
                            Para_s->d_viv[k * jmt * imt + j * imt + i];

            }

            CRTS_dma_put(&Para_s->d_wka[k * imt * jmt + j * imt], d_wka1, sizeof(double) * imt);
        }

    }

    ldm_free(Para_s, sizeof(bclinc4_0_Para));
    return;
}

void bclinc4_1_cu_slave_kernel_(bclinc4_Para *Paras){

    double d_wka[imt];
    double d_vp[imt];
    double d_v[imt];
    //ATHREAD basic operation
    int i, j;
    blockIndex = CRTS_tid;
    bclinc4_Para *Para_s = (bclinc4_Para*)ldm_malloc(sizeof(bclinc4_Para));
    CRTS_dma_iget(Para_s, Paras, sizeof(bclinc4_Para), &dma_rply);
    D_COUNT++;
    CRTS_dma_wait_value(&dma_rply, D_COUNT);

#ifdef BC4DEBUG
    unsigned long count,counte;
    unsigned long misscount,misscounte;

    penv_slave1_dcache_access_init();
    penv_slave2_dcache_miss_init();

    penv_slave1_dcache_access_count(&count);
    penv_slave2_dcache_miss_count(&misscount);
#endif
    //cuda grid and block affine translation
    int blockSum   = Para_s->cuGrid[0] * Para_s->cuGrid[1];
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
        j = blockId_y * Para_s->cuBlock[1];
        // ATHREAD serial operation replace GPU parallel mutli-thread in single thread-block
        for(i = 0; i < threadSum; i++) {

            if (i >= 2 && i < imt - 2 && j >= 2 && j < jmt - 2) {
                int kmb = Para_s->d_kmu[j * imt + i] - 1, k;

                if (kmb >= 0) {
                    double g0;
                    for (k = 1; k <= kmb; k++) {
                        d_a8[k] =
                                Para_s->d_akmu[(k - 1) * jmt * imt + j * imt + i] *
                                Para_s->d_odzt[k] * Para_s->d_odzp[k] * Para_s->d_dtc2 * 0.5;
                        d_d8[k] =
                                Para_s->d_wka[k * jmt * imt + j * imt + i];
                        d_c8[k] =
                                Para_s->d_akmu[k * jmt * imt + j * imt + i] *
                                Para_s->d_odzt[k + 1] * Para_s->d_odzp[k] * Para_s->d_dtc2 * 0.5;
                        d_b8[k] =
                                1.0 + d_a8[k] +
                                d_c8[k];
                        d_e8[k] = 0.0;
                        d_f8[k] = 0.0;
                    }

                    //  B. C. AT TOP
                    k = 0;
                    d_a8[k] = Para_s->d_odzp[k] * Para_s->d_dtc2 * 0.5;
                    d_c8[k] =
                            Para_s->d_akmu[k * jmt * imt + j * imt + i] *
                            Para_s->d_odzt[k + 1] * Para_s->d_odzp[k] * Para_s->d_dtc2 * 0.5;
                    d_b8[k] = 1.0 + d_c8[k];
                    d_d8[k] = Para_s->d_wka[k * jmt * imt + j * imt + i];

                    //B. C. AT BOTTOM
                    d_b8[kmb] = 1.0 + d_a8[kmb];
                    d_c8[kmb] = Para_s->d_odzp[kmb] * Para_s->d_dtc2 * 0.5;
                    d_e8[(kmb + 1)] = 0.0;
                    d_f8[(kmb + 1)] = 0.0;
                    d_d8[kmb] =
                            Para_s->d_wka[kmb * jmt * imt + j * imt + i] -
                            Para_s->d_bbcx[j * imt + i] * Para_s->d_odzp[kmb] * Para_s->d_dtc2 * 0.5;

                    //NOW INVERT
                    for (k = kmb; k >= 0; k--) {
                        g0 = 1.0 / (d_b8[k] -
                                    d_c8[k] *
                                    d_e8[(k + 1)]);
                        d_e8[k] = d_a8[k] * g0;
                        d_f8[k] =
                                (d_d8[k] +
                                 d_c8[k] *
                                 d_f8[(k + 1)]) * g0;
                    }

                    //B.C. AT SURFACE
                    Para_s->d_wka[j * imt + i] =
                            (d_e8[0] * Para_s->d_sbcx[j * imt + i] +
                             d_f8[0]) * Para_s->d_viv[j * imt + i];
                    for (k = 1; k <= kmb; k++) {
                        Para_s->d_wka[k * jmt * imt + j * imt + i] =
                                (d_e8[k] *
                                 Para_s->d_wka[(k - 1) * jmt * imt + j * imt + i] +
                                 d_f8[k]) *
                                Para_s->d_viv[k * jmt * imt + j * imt + i];
                    }
                }
            }

        }
    }



#ifdef BC4DEBUG
    penv_slave1_dcache_access_count(&counte);
    penv_slave2_dcache_miss_count(&misscounte);
    if(Para_s->rank == 0)
        printf("[slave %d]:bclinc4_1 count:%d,counte:%d,misscount:%d,misscounte:%d,total count:%d, total misscount:%d, misscount rate:%f\n",
               blockIndex,count,counte,misscount,misscounte,counte-count, misscounte-misscount,
               (float)(misscounte-misscount)/(float)(counte-count));
#endif
    ldm_free(Para_s, sizeof(bclinc4_Para));
    return;
}
