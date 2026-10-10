#include <crts.h>
#include <slave.h>
#include <simd.h>
#include <math.h>
#include <stdio.h>
#include <float.h>
#include <string.h>
#include <stdbool.h>
#include "param_mod.h"
#include "swperf.h"

__thread_local unsigned int D_COUNT = 0;
__thread_local crts_rply_t dma_rply = 0;


typedef struct bclinc_vinteg_kernelPara{
    int cuGrid[3];
    int cuBlock[3];
    int rank;
    int vinegtype;
    double *d_u;
    double *d_v;
    double *d_wka;
    double *d_work;
    double *d_dzp;
    double *d_ohbu;
    double *d_viv;
}bclinc_vinteg_Para;

void  vinteg_cu_slave_kernel_(bclinc_vinteg_Para *Paras){

    //ATHREAD basic operation
    int blockIndex = CRTS_tid;
    bclinc_vinteg_Para* Para_s = (bclinc_vinteg_Para*)ldm_malloc(sizeof(bclinc_vinteg_Para));
    CRTS_dma_get(Para_s, Paras, sizeof(bclinc_vinteg_Para));

    double *wk3;
    if(Para_s->vinegtype == 1)
        wk3 = Para_s->d_u;
    else if (Para_s->vinegtype == 2)
        wk3 = Para_s->d_v;
    else if (Para_s->vinegtype == 3)
        wk3 = Para_s->d_wka;

    double *viv = Para_s->d_viv;
    double *wk2 = Para_s->d_work;
    double *dzp = Para_s->d_dzp;
    double *ohbu = Para_s->d_ohbu;

    int tile  = imt ;
    int avaiLdm = get_allocatable_size();
    int varNum = 2; // variable of dma stride
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

    // target variable for k direction inconsistent access
    double* ldm_wk3 = (double*) ldm_malloc (tileSum);
    double* ldm_viv = (double*) ldm_malloc (tileSum);

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

            athread_dma_iget_stride(ldm_wk3, &wk3[start_blk],
                                    tileSumV, tile_size, size_stride, &dma_rply);
            D_COUNT++;
            athread_dma_wait_value(&dma_rply, D_COUNT);

            athread_dma_iget_stride(ldm_viv, &viv[start_blk],
                                    tileSumV, tile_size, size_stride, &dma_rply);
            D_COUNT++;
            athread_dma_wait_value(&dma_rply, D_COUNT);

            for (int i = 0; i < tileV; ++i) {
                int index = start_blk + i;
                for (int k = 0; k < km; ++k) {
                    wk2[index] += dzp[k] * ohbu[index] *
                                  ldm_wk3[k * tileV + i] * ldm_viv[k * tileV + i];
                }
            }
        }
    }

//                    WK2 (I,J,IBLOCK)= WK2 (I,J,IBLOCK) + DZP (K)* OHBU(I,J,IBLOCK)*  &
//                            WK3(I,J,K,IBLOCK) *VIV (I,J,K,IBLOCK)

    ldm_free (ldm_viv, tileSum);
    ldm_free (ldm_wk3, tileSum);

    ldm_free(Para_s, sizeof(bclinc_vinteg_Para));
}


void  vinteg_cu_slave_kernel_old_(bclinc_vinteg_Para *Paras){

    double d_work[imt];
    double d_wk3[imt];

    //ATHREAD basic operation
    int blockIndex = CRTS_tid;
    bclinc_vinteg_Para* Para_s = (bclinc_vinteg_Para*)ldm_malloc(sizeof(bclinc_vinteg_Para));
    CRTS_dma_get(Para_s, Paras, sizeof(bclinc_vinteg_Para));

    for(int k = 0; k < km; k++) {
        for (int j = blockIndex; j < jmt; j += 64) {

            CRTS_dma_get(d_work, &Para_s->d_work[j * imt], sizeof(double) * imt);
            if(Para_s->vinegtype == 1)
                CRTS_dma_get(d_wk3, &Para_s->d_u[k * imt * jmt + j * imt], sizeof(double) * imt);
            else if(Para_s->vinegtype == 2)
                CRTS_dma_get(d_wk3, &Para_s->d_v[k * imt * jmt + j * imt], sizeof(double) * imt);
            else if(Para_s->vinegtype == 3)
                CRTS_dma_get(d_wk3, &Para_s->d_wka[k * imt * jmt + j * imt], sizeof(double) * imt);

            for (int i = 0; i < imt; i++) {
//                    WK2 (I,J,IBLOCK)= WK2 (I,J,IBLOCK) + DZP (K)* OHBU(I,J,IBLOCK)*  &
//                            WK3(I,J,K,IBLOCK) *VIV (I,J,K,IBLOCK)

                d_work[i] = d_work[i] + Para_s->d_dzp[k] *Para_s->d_ohbu[j * imt + i] *
                            d_wk3[i] *Para_s->d_viv[k *jmt *imt+j * imt + i];

            }

            CRTS_dma_put(&Para_s->d_work[j*imt], d_work, sizeof(double) * imt);
        }
    }


    ldm_free(Para_s, sizeof(bclinc_vinteg_Para));
}

void  vinteg_cu_slave_kernel_gb_(bclinc_vinteg_Para *Paras){

//    double d_work[imt];
//    double d_wk3[imt];
    double* wk3;
    double *viv;
    //ATHREAD basic operation
    int blockIndex = CRTS_tid;
    bclinc_vinteg_Para* Para_s = (bclinc_vinteg_Para*)ldm_malloc(sizeof(bclinc_vinteg_Para));
    CRTS_dma_get(Para_s, Paras, sizeof(bclinc_vinteg_Para));

    if(Para_s->vinegtype == 1)
        wk3 = Para_s->d_u;
    else if (Para_s->vinegtype == 2)
        wk3 = Para_s->d_v;
    else if (Para_s->vinegtype == 3)
        wk3 = Para_s->d_wka;

    viv = Para_s->d_viv;
    double* wk2 = Para_s->d_work;
    double* dzp = Para_s->d_dzp;
    double* ohbu = Para_s->d_ohbu;

    int total_len = imt * jmt;
    int blk_per_cpe = (total_len + 63) / 64;
    int varNum = 2;
    int avail_ldm = 96 * 1024;
    int allocatable_num_blk = avail_ldm / (varNum * km * sizeof(double));

    int num_blk = (blk_per_cpe + allocatable_num_blk - 1)
                  / allocatable_num_blk;
    int local_blk = blk_per_cpe;
    if (blk_per_cpe > allocatable_num_blk) {
        local_blk = allocatable_num_blk;
    }
    int remainder_blk = blk_per_cpe % local_blk;

    size_t size_tile = local_blk * sizeof(double);
    size_t size_blk = size_tile * km;
    size_t size_stride_layer = (total_len - local_blk) * sizeof(double);

    double* ldm_wk3 = (double*) ldm_malloc (size_blk);
    double* ldm_viv = (double*) ldm_malloc (size_blk);

    int start_cpe = CRTS_tid * blk_per_cpe;

    for (int idx_blk = 0; idx_blk < num_blk; ++idx_blk) {
        int start_blk = idx_blk * local_blk;
        if ((idx_blk != num_blk - 1) || (idx_blk == 0)
            || (idx_blk == num_blk - 1 && remainder_blk == 0)) {
            // No blocked
            athread_dma_iget_stride (ldm_wk3, &wk3[start_cpe + start_blk],
                                     size_blk, size_tile, size_stride_layer, &dma_rply);
            D_COUNT++;
            athread_dma_iget_stride (ldm_viv, &viv[start_cpe + start_blk],
                                     size_blk, size_tile, size_stride_layer, &dma_rply);
            D_COUNT++;

            for (int idx = 0; idx < local_blk; ++idx) {
                int ii = start_blk + idx;
                if (ii >= blk_per_cpe) {break;}
                int index = start_cpe + ii;
                if (index >= total_len) {break;}
                wk2[index] = 0.0;
            }

            athread_dma_wait_value(&dma_rply, D_COUNT);
            athread_dma_wait_value(&dma_rply, D_COUNT);

            for (int idx = 0; idx < local_blk; ++idx) {
                int ii = start_blk + idx;
                if (ii >= blk_per_cpe) {break;}
                int index = start_cpe + ii;
                if (index >= total_len) {break;}
                for (int k = 0; k < km; ++k) {
                    wk2[index] += dzp[k] * ohbu[index] *
                                  ldm_wk3[k*local_blk + idx] * ldm_viv[k*local_blk + idx];
                }
            }
        } else {
            // No blocked
            athread_dma_iget_stride (ldm_wk3, &wk3[start_cpe + start_blk],
                                     remainder_blk * sizeof(double) * km,
                                     remainder_blk * sizeof(double),
                                     (total_len - remainder_blk) * sizeof(double), &dma_rply);
            D_COUNT++;
            athread_dma_iget_stride (ldm_viv, &viv[start_cpe + start_blk],
                                     remainder_blk * sizeof(double) * km,
                                     remainder_blk * sizeof(double),
                                     (total_len - remainder_blk) * sizeof(double), &dma_rply);
            D_COUNT++;

            for (int idx = 0; idx < local_blk; ++idx) {
                int ii = start_blk + idx;
                if (ii >= blk_per_cpe) {break;}
                int index = start_cpe + ii;
                if (index >= total_len) {break;}
                wk2[index] = 0.0;
            }

            athread_dma_wait_value(&dma_rply, D_COUNT);
            athread_dma_wait_value(&dma_rply, D_COUNT);

            for (int idx = 0; idx < local_blk; ++idx) {
                int ii = start_blk + idx;
                if (ii >= blk_per_cpe) {break;}
                int index = start_cpe + ii;
                if (index >= total_len) {break;}
                for (int k = 0; k < km; ++k) {
                    wk2[index] += dzp[k] * ohbu[index] *
                                  ldm_wk3[k*remainder_blk + idx] * ldm_viv[k*remainder_blk + idx];
                }
            }
        }
    }
    ldm_free (ldm_viv, size_blk);
    ldm_free (ldm_wk3, size_blk);



    ldm_free(Para_s, sizeof(bclinc_vinteg_Para));
}