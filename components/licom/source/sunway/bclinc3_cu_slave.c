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


typedef struct bclinc3_kernelPara{
    int    cuBlock[3];
    int    cuGrid[3];
    double d_dtc2;
    int   *d_kmu;
    double *d_odzt;
    double *d_odzp;
    double *d_viv;
    double *d_wk;
    double *d_topbc;
    double *d_bombc;
    double *d_dcb;
}bclinc3_Para;



//ATHREAD global definition 
__thread_local int blockIndex;

void bclinc3_cu_slave_kernel_(bclinc3_Para *Paras){
	//ATHREAD basic operation
    blockIndex = CRTS_tid;
    bclinc3_Para* Para_s = (bclinc3_Para*)ldm_malloc(sizeof(bclinc3_Para));
    CRTS_dma_iget(Para_s, Paras, sizeof(bclinc3_Para), &dma_rply);
    D_COUNT++;
    CRTS_dma_wait_value(&dma_rply, D_COUNT);

    int tile  = imt ;
    int avaiLdm = get_allocatable_size();
    int varNum = 3; // variable of dma stride
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

    double d_a8[km],d_b8[km],d_c8[km],d_d8[km],d_e8[km+1],d_f8[km+1];

    // target variable for k direction inconsistent access
    double* ldm_wk  = (double*) ldm_malloc (tileSum);
    double* ldm_dcb  = (double*) ldm_malloc (tileSum);
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

            athread_dma_iget_stride(ldm_wk, &Para_s->d_wk[start_blk],
                                    tileSumV, tile_size, size_stride, &dma_rply);
            D_COUNT++;
            athread_dma_wait_value(&dma_rply, D_COUNT);

            athread_dma_iget_stride(ldm_viv, &Para_s->d_viv[start_blk],
                                    tileSumV, tile_size, size_stride, &dma_rply);
            D_COUNT++;
            athread_dma_wait_value(&dma_rply, D_COUNT);

            athread_dma_iget_stride(ldm_dcb, &Para_s->d_dcb[start_blk],
                                    tileSumV, tile_size, size_stride, &dma_rply);
            D_COUNT++;
            athread_dma_wait_value(&dma_rply, D_COUNT);
            for (int i = 0; i < tileV; ++i) {
                int index = start_blk + i;
                int ii = idx * tile + i;
                if (ii >= 2 && ii < imt - 2 && j >= 2 && j < jmt - 2) {
                    int kmb = Para_s->d_kmu[index] - 1, k;

                    if (kmb >= 0) {
                        double g0;
                        for (k = 1; k <= kmb; k++) {
                            d_a8[k] = ldm_dcb[(k - 1) * tileV + i] *
                                      Para_s->d_odzt[k] * Para_s->d_odzp[k] * Para_s->d_dtc2 * 0.5;
                            d_d8[k] = ldm_wk[k * tileV + i];
                            d_c8[k] = ldm_dcb[k * tileV + i] *
                                    Para_s->d_odzt[k + 1] * Para_s->d_odzp[k] * Para_s->d_dtc2 * 0.5;
                            d_b8[k] = 1.0 + d_a8[k] + d_c8[k];
                            d_e8[k] = 0.0;
                            d_f8[k] = 0.0;
                        }

                        //  B. C. AT TOP
                        k = 0;
                        d_a8[k] = Para_s->d_odzp[k] * Para_s->d_dtc2 * 0.5;
                        d_c8[k] = ldm_dcb[k * tileV + i] *
                                Para_s->d_odzt[k + 1] * Para_s->d_odzp[k] * Para_s->d_dtc2 * 0.5;
                        d_b8[k] = 1.0 + d_c8[k];
                        d_d8[k] = ldm_wk[k * tileV + i];

                        //B. C. AT BOTTOM
                        d_b8[kmb] = 1.0 + d_a8[kmb];
                        d_c8[kmb] = Para_s->d_odzp[kmb] * Para_s->d_dtc2 * 0.5;
                        d_e8[(kmb + 1)] = 0.0;
                        d_f8[(kmb + 1)] = 0.0;
                        d_d8[kmb] = ldm_wk[kmb * tileV + i] -
                                Para_s->d_bombc[index] * Para_s->d_odzp[kmb] * Para_s->d_dtc2 * 0.5;
                        //Para_s->d_bbcy[j * imt + i] * Para_s->d_odzp[kmb] * Para_s->d_dtc2 * 0.5;

                        //NOW INVERT
                        for (k = kmb; k >= 0; k--) {
                            g0 = 1.0 / (d_b8[k] - d_c8[k] * d_e8[(k + 1)]);
                            d_e8[k] = d_a8[k] * g0;
                            d_f8[k] = (d_d8[k] + d_c8[k] * d_f8[(k + 1)]) * g0;
                        }

                        //B.C. AT SURFACE
                        ldm_wk[i] = (d_e8[0] * Para_s->d_topbc[index] + d_f8[0]) * ldm_viv[i];

                        for (k = 1; k <= kmb; k++) {
                            ldm_wk[k * tileV + i] = (d_e8[k] *
                                     ldm_wk[(k - 1) * tileV + i] + d_f8[k]) * ldm_viv[k * tileV + i];
                        }
                    }
                }
            }

            athread_dma_iput_stride (&Para_s->d_wk[start_blk],
                                     ldm_wk, tileSumV, tile_size, size_stride,&dma_rply);
            D_COUNT++;
            athread_dma_wait_value(&dma_rply, D_COUNT);

        }
    }



    ldm_free (ldm_wk, tileSum);
    ldm_free (ldm_dcb, tileSum);
    ldm_free (ldm_viv, tileSum);
    ldm_free(Para_s, sizeof(bclinc3_Para));
    return;    
}           
            
