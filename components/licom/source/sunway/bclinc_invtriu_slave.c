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


typedef struct bc_invtriu_kernelPara{
    int    cuBlock[3];
    int    cuGrid[3];
    double d_dtc2;
    int   *d_kmu;
    double *d_akmu;
    double *d_odzt;
    double *d_odzp;
    double *d_u;
    double *d_v;
    double *d_sbcx;
    double *d_sbcy;
    double *d_bbcx;
    double *d_bbcy;
    double *d_viv;
    //double *d_wk1;
    //double *d_wk2;
}Invtriu_Para;


__thread_local int blockIndex;

void bclinc_invtriu_kernel_(Invtriu_Para *Paras){
	//ATHREAD basic operation
    blockIndex = CRTS_tid;
    Invtriu_Para *Para_s = (Invtriu_Para*)ldm_malloc(sizeof(Invtriu_Para));
    CRTS_dma_iget(Para_s, Paras, sizeof(Invtriu_Para), &dma_rply);
    D_COUNT++;
    CRTS_dma_wait_value(&dma_rply, D_COUNT);

    int tile  = imt ;
    int avaiLdm = get_allocatable_size();
    int varNum = 4; // variable of dma stride
    int tileBlk = 1; //blocks number in i direction
    int tileLimit = avaiLdm/(varNum*km*sizeof(double));
    int remTile = 0;
    if(tile > tileLimit) {
        remTile = tile % tileLimit;
        if(remTile == 0)
            tileBlk = tile / tileLimit;
        else
            tileBlk = tile / tileLimit + 1;
        tile = tileLimit;
    }

    int start_blk;
    int tile_size  = tile * sizeof(double);
    int tileSum   = tile * km * sizeof(double);
    int size_stride = (imt * jmt - tile) * sizeof(double);
    int tileV;
    int tileSumV;
    int size_strideV;

    double d_a8[km],d_b8[km],d_c8[km],d_e8[km+1];
    double d8_u[km],d8_v[km],f8_u[km+1],f8_v[km+1];

    // target variable for k direction inconsistent access
    double* ldm_akmu = (double*) ldm_malloc (tileSum);
    double* ldm_wk1  = (double*) ldm_malloc (tileSum);
    double* ldm_wk2  = (double*) ldm_malloc (tileSum);
    double* ldm_viv  = (double*) ldm_malloc (tileSum);

    // each cpe is responsible for i*k data
    for(int j = blockIndex; j < jmt; j += 64){
        // blocks in i direction
        for(int idx=0; idx < tileBlk; idx ++) {
            start_blk = j * imt + idx * tile;
            //distinguish tile and remainder for stride dma
            if((remTile != 0) && (idx == tileBlk -1)) {
                tileV = remTile;
                tileSumV   = remTile * km * sizeof(double);
                tile_size  = remTile * sizeof(double);
                size_stride = (imt * jmt - remTile) * sizeof(double);
            }
            else {
                tileV = tile;
                tileSumV   = tile * km * sizeof(double);
                tile_size  = tile * sizeof(double);
                size_stride = (imt * jmt - tile) * sizeof(double);
            }

            athread_dma_iget_stride(ldm_akmu, &Para_s->d_akmu[start_blk],
                                    tileSumV, tile_size, size_stride, &dma_rply);
            D_COUNT++;
            athread_dma_wait_value(&dma_rply, D_COUNT);

            athread_dma_iget_stride(ldm_wk1, &Para_s->d_u[start_blk],
                                    tileSumV, tile_size, size_stride, &dma_rply);
            D_COUNT++;
            athread_dma_wait_value(&dma_rply, D_COUNT);

            athread_dma_iget_stride(ldm_wk2, &Para_s->d_v[start_blk],
                                    tileSumV, tile_size, size_stride, &dma_rply);
            D_COUNT++;
            athread_dma_wait_value(&dma_rply, D_COUNT);

            athread_dma_iget_stride(ldm_viv, &Para_s->d_viv[start_blk],
                                    tileSumV, tile_size, size_stride, &dma_rply);
            D_COUNT++;
            athread_dma_wait_value(&dma_rply, D_COUNT);

            for (int i = 0; i < tileV; ++i) {
                int index = start_blk + i;
                int ii = idx * tile + i;
                if (ii >= 2 && ii < imt - 2 && j >= 2 && j < jmt - 2) {
                    //int kmb = Para_s->d_kmu[j * imt + i] - 1, k;
                    //KZ = kmu (I,J,iblock)
                    int KZ = Para_s->d_kmu[index] - 1;
                    int k;
                    //IF (kmu(I,J,iblock) > 0) then
                    if (KZ >= 0) {
                        double g0;
                        for (k = 1; k <= KZ; k++) {
                            //A8 (K) = DCB (I,J,K -1,iblock)* ODZT (K )* ODZP (K)* C2DTC * AIDIF
                            //D8 (K) = WK (I,J,K,iblock)
                           // d_a8[k] = Para_s->d_akmu[(k - 1) * jmt * imt + j * imt + i] *
                           //         Para_s->d_odzt[k] * Para_s->d_odzp[k] * Para_s->d_dtc2 * 0.5;

                            d_a8[k] = ldm_akmu[(k - 1) * tileV + i] *
                                      Para_s->d_odzt[k] * Para_s->d_odzp[k] * Para_s->d_dtc2 * 0.5;

                            d8_u[k] = ldm_wk1[k * tileV + i];
                            d8_v[k] = ldm_wk2[k * tileV + i];


                            d_c8[k] = ldm_akmu[k * tileV + i] *
                                    Para_s->d_odzt[k + 1] * Para_s->d_odzp[k] * Para_s->d_dtc2 * 0.5;
                            d_b8[k] = 1.0 + d_a8[k] + d_c8[k];
                            d_e8[k] = 0.0;
                            f8_u[k] = 0.0;
                            f8_v[k] = 0.0;
                        }

                        //  B. C. AT TOP
                        k = 0;
                        d_a8[k] = Para_s->d_odzp[k] * Para_s->d_dtc2 * 0.5;
                        d_c8[k] = ldm_akmu[k * tileV + i] *
                                Para_s->d_odzt[k + 1] * Para_s->d_odzp[k] * Para_s->d_dtc2 * 0.5;
                        d_b8[k] = 1.0 + d_c8[k];
                        d8_u[k] = ldm_wk1[k * tileV + i];
                        d8_v[k] = ldm_wk2[k * tileV + i];

                        //B. C. AT BOTTOM
                        d_b8[KZ] = 1.0 + d_a8[KZ];
                        d_c8[KZ] = Para_s->d_odzp[KZ] * Para_s->d_dtc2 * 0.5;
                        d_e8[(KZ + 1)] = 0.0;
                        f8_u[(KZ + 1)] = 0.0;
                        f8_v[(KZ + 1)] = 0.0;
                        d8_u[KZ] = ldm_wk1[KZ * tileV + i] -
                                    Para_s->d_bbcx[index] * Para_s->d_odzp[KZ] * Para_s->d_dtc2 * 0.5;
                        d8_v[KZ] = ldm_wk2[KZ * tileV + i] -
                                    Para_s->d_bbcy[index] * Para_s->d_odzp[KZ] * Para_s->d_dtc2 * 0.5;
                        //NOW INVERT
                        for (k = KZ; k >= 0; k--) {
                            g0 = 1.0 / (d_b8[k] - d_c8[k] * d_e8[(k + 1)]);
                            d_e8[k] = d_a8[k] * g0;
                            f8_u[k] = (d8_u[k] + d_c8[k] * f8_u[(k + 1)]) * g0;
                            f8_v[k] = (d8_v[k] + d_c8[k] * f8_v[(k + 1)]) * g0;
                        }

                        //B.C. AT SURFACE
                       // Para_s->d_u[j * imt + i] = (d_e8[0] * Para_s->d_sbcx[j * imt + i] +
                       //                             f8_u[0]) * Para_s->d_viv[j * imt + i];

                        ldm_wk1[i] = (d_e8[0] * Para_s->d_sbcx[index] + f8_u[0]) * ldm_viv[i];


                        ldm_wk2[i] = (d_e8[0] * Para_s->d_sbcy[index] + f8_v[0]) * ldm_viv[i];
                        for (k = 1; k <= KZ; k++) {
                            ldm_wk1[k * tileV + i] = (d_e8[k] * ldm_wk1[(k - 1) * tileV + i] +
                                     f8_u[k]) * ldm_viv[k * tileV + i];
                            ldm_wk2[k * tileV + i] = (d_e8[k] * ldm_wk2[(k - 1) * tileV + i] +
                                     f8_v[k]) * ldm_viv[k * tileV + i];
                        }
                    }
                }
            }

            athread_dma_iput_stride (&Para_s->d_u[start_blk],
                                     ldm_wk1, tileSumV, tile_size, size_stride,&dma_rply);
            D_COUNT++;
            athread_dma_wait_value(&dma_rply, D_COUNT);

            athread_dma_iput_stride (&Para_s->d_v[start_blk],
                                     ldm_wk2, tileSumV, tile_size, size_stride,&dma_rply);
            D_COUNT++;
            athread_dma_wait_value(&dma_rply, D_COUNT);
        }
    }


    ldm_free (ldm_viv, tileSum);
    ldm_free (ldm_akmu, tileSum);
    ldm_free (ldm_wk1, tileSum);
    ldm_free (ldm_wk2, tileSum);
    ldm_free(Para_s, sizeof(Invtriu_Para));
    return;    
}           
            
