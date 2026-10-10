
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

//__userDef__ paraVarAttr(int,d_nx,sizeof(int),readonly,no,no,no);	
//__userDef__ paraVarAttr(int,d_jmt,sizeof(int),readonly,no,no,no);	
//__userDef__ paraVarAttr(int,d_imt,sizeof(int),readonly,no,no,no);	
//__userDef__ paraVarAttr(int,d_kmt,sizeof(int)*ny_block*nx_block,readonly,no,no,no);	
//__userDef__ paraVarAttr(double,d_wkc,sizeof(double)*km*jmt*imt,readonly,no,no,no);	
//__userDef__ paraVarAttr(double,d_vtl,sizeof(double)*km*jmt*imt,outonly,yes,no,no);	
//__userDef__ paraVarAttr(double,d_stf,sizeof(double)*jmt*imt,readonly,no,no,no);	
//__userDef__ paraVarAttr(double,d_vit,sizeof(double)*jmt*imt*km,readonly,no,no,no);	
//__userDef__ paraVarAttr(double,d_odzt,sizeof(double)*km,readonly,no,yes,no);	
//__userDef__ paraVarAttr(double,d_odzp,sizeof(double)*km,readonly,no,yes,no);	
//__userDef__ paraVarAttr(double,d_dts,sizeof(double),readonly,no,no,no);	
typedef struct invtrit_kernelPara{
	int cuBlock[3];
	int cuGrid[3];
	int rank;
	int *d_kmt;
	double d_dts;
	double *d_wkc;
	double *d_vtl;
	double *d_stf;
	double *d_vit;
	double *d_odzt;
	double *d_odzp;
}invtrit_Para;


__thread_local int blockIndex;

void invtrit_(invtrit_Para *Paras){


  //ATHREAD basic operation
  int i,j,k;
  blockIndex = CRTS_tid;
  invtrit_Para *Para_s = (invtrit_Para*)ldm_malloc(sizeof(invtrit_Para));
  CRTS_dma_iget(Para_s, Paras, sizeof(invtrit_Para), &dma_rply);
  D_COUNT++;
  CRTS_dma_wait_value(&dma_rply, D_COUNT);

  double* d_odzt_s = (double*)ldm_malloc(sizeof(double) *km);
  CRTS_dma_iget(d_odzt_s, Para_s->d_odzt, sizeof(double) *km, &dma_rply);
  D_COUNT++;
  CRTS_dma_wait_value(&dma_rply, D_COUNT);

  double* d_odzp_s = (double*)ldm_malloc(sizeof(double) *km);
  CRTS_dma_iget(d_odzp_s, Para_s->d_odzp, sizeof(double) *km, &dma_rply);
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

    double a_8[km],b_8[km],c_8[km],d_8[km];
    double e_8[km+1],f_8[km+1];

    // target variable for k direction inconsistent access
    double* ldm_wkc  = (double*) ldm_malloc (tileSum);
    double* ldm_vtl  = (double*) ldm_malloc (tileSum);
    double* ldm_vit  = (double*) ldm_malloc (tileSum);


    // each cpe is responsible for i*k data
    for(int j = blockIndex; j < jmt; j += 64){



        // blocks in i direction
        for(int idx=0; idx < tileBlk; idx ++) {
      
            double g0;
            int kz;
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

            athread_dma_iget_stride(ldm_wkc, &Para_s->d_wkc[start_blk],
                                    tileSumV, tile_size, size_stride, &dma_rply);
            D_COUNT++;
            athread_dma_wait_value(&dma_rply, D_COUNT);

            athread_dma_iget_stride(ldm_vit, &Para_s->d_vit[start_blk],
                                    tileSumV, tile_size, size_stride, &dma_rply);
            D_COUNT++;
            athread_dma_wait_value(&dma_rply, D_COUNT);

            athread_dma_iget_stride(ldm_vtl, &Para_s->d_vtl[start_blk],
                                    tileSumV, tile_size, size_stride, &dma_rply);
            D_COUNT++;
            athread_dma_wait_value(&dma_rply, D_COUNT);

            for (int i = 0; i < tileV; ++i) {
                int index = start_blk + i;
                int ii = idx * tile + i;
                if (j >= 2 && j < jmt - 2 && ii >= 2 && ii < imt - 2) {
                    if (Para_s->d_kmt[index] > 0) {
                        kz = Para_s->d_kmt[index];
                        for (k = 1; k < kz; k++) {
                            a_8[k] = ldm_wkc[(k - 1) * tileV + i] * d_odzt_s[k] * d_odzp_s[k] *
                                     Para_s->d_dts * 0.5e0;
                            d_8[k] = ldm_vtl[k * tileV + i];
                        }
                        for (k = 1; k < kz - 1; k++) {
                            c_8[k] = ldm_wkc[k * tileV + i] * d_odzt_s[k + 1] * d_odzp_s[k] *
                                     Para_s->d_dts * 0.5e0;
                            b_8[k] = 1.0 + a_8[k] + c_8[k];
                            e_8[k] = 0.0;
                            f_8[k] = 0.0;
                        }

                        k = 0;
                        a_8[k] = d_odzp_s[k] * Para_s->d_dts * 0.5e0;
                        c_8[k] = ldm_wkc[k * tileV + i] * d_odzt_s[k + 1] * d_odzp_s[k] *
                                 Para_s->d_dts * 0.5e0;
                        b_8[k] = 1.0 + c_8[k];
                        d_8[k] = ldm_vtl[k * tileV + i];
                        e_8[k] = 0.0;
                        f_8[k] = 0.0;
                        b_8[kz - 1] = 1.0 + a_8[kz - 1];
                        c_8[kz - 1] = d_odzp_s[kz - 1] * Para_s->d_dts * 0.5e0;
                        e_8[kz] = 0.0;
                        f_8[kz] = 0.0;

                        for (k = kz - 1; k >= 0; k--) {
                            g0 = 1.0 / (b_8[k] - c_8[k] * e_8[k + 1]);
                            e_8[k] = a_8[k] * g0;
                            f_8[k] = (d_8[k] + c_8[k] * f_8[k + 1]) * g0;
                        }
                        ldm_vtl[i] = (e_8[0] * Para_s->d_stf[index] + f_8[0]) * ldm_vit[i];

                        for (k = 1; k < kz; k++) {
                            ldm_vtl[k * tileV + i] = (e_8[k] * ldm_vtl[(k - 1) * tileV + i] + f_8[k]) *
                                            ldm_vit[k * tileV + i];
                        }
                    }
                }
            }
            athread_dma_iput_stride (&Para_s->d_vtl[start_blk],
                                     ldm_vtl, tileSumV, tile_size, size_stride,&dma_rply);
            D_COUNT++;
            athread_dma_wait_value(&dma_rply, D_COUNT);
    }

  }
    ldm_free (ldm_wkc, tileSum);
    ldm_free (ldm_vtl, tileSum);
    ldm_free (ldm_vit, tileSum);
    ldm_free(d_odzt_s, sizeof(double) *km);

    ldm_free(d_odzp_s, sizeof(double) *km);

    ldm_free(Para_s, sizeof(invtrit_Para));
}

void invtrit_old_(invtrit_Para *Paras){

    //ATHREAD local variable definition
    int swGrid[3];
    int swBlock[3];
    int stride;
    int blockId_x,blockId_y,blockId_z;
    int threadId_x,threadId_y,threadId_z;


    //ATHREAD basic operation
    int i,j,k;
    blockIndex = CRTS_tid;
    invtrit_Para *Para_s = (invtrit_Para*)ldm_malloc(sizeof(invtrit_Para));
    CRTS_dma_iget(Para_s, Paras, sizeof(invtrit_Para), &dma_rply);
    D_COUNT++;
    CRTS_dma_wait_value(&dma_rply, D_COUNT);

    double* d_odzt_s = (double*)ldm_malloc(sizeof(double) *km);
    CRTS_dma_iget(d_odzt_s, Para_s->d_odzt, sizeof(double) *km, &dma_rply);
    D_COUNT++;
    CRTS_dma_wait_value(&dma_rply, D_COUNT);

    double* d_odzp_s = (double*)ldm_malloc(sizeof(double) *km);
    CRTS_dma_iget(d_odzp_s, Para_s->d_odzp, sizeof(double) *km, &dma_rply);
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

    double *d_vtl_s = Para_s->d_vtl;

    for(int block = blockStart; block < blockEnd; block++){
        blockId_z = block / (Para_s->cuGrid[0] * Para_s->cuGrid[1]);
        blockId_y = (block - blockId_z * (Para_s->cuGrid[0] * Para_s->cuGrid[1])) / Para_s->cuGrid[0];
        blockId_x = block % Para_s->cuGrid[0];
        j = blockId_y * Para_s->cuBlock[1];


        //ATHREAD serial operation replace GPU parallel mutli-thread in single thread-block
        for(i = 0; i < threadSum; i++) {

            double a_8[km],b_8[km],c_8[km],d_8[km];
            double e_8[km+1],f_8[km+1];
            double g0;
            int kz;

            if(j>=2&&j<jmt-2&&i>=2&&i<imt-2)
            {
                if(Para_s->d_kmt[j *imt+i]> 0)
                {
                    kz=Para_s->d_kmt[j *imt+i];
                    for(k=1;k<kz; k++)
                    {
                        a_8[k]=Para_s->d_wkc[(k-1) *jmt *imt+j *imt+i] *d_odzt_s[k] *d_odzp_s[k] *Para_s->d_dts *0.5e0;
                        d_8[k]=d_vtl_s[k *jmt *imt+j *imt+i];
                    }
                    for(k=1;k<kz-1; k++)
                    {
                        c_8[k]=Para_s->d_wkc[k *jmt *imt+j *imt+i] *d_odzt_s[k+1] *d_odzp_s[k] *Para_s->d_dts *0.5e0;
                        b_8[k]=1.0+a_8[k]+c_8[k];
                        e_8[k]=0.0;
                        f_8[k]=0.0;
                    }

                    k=0;
                    a_8[k]=d_odzp_s[k] *Para_s->d_dts *0.5e0;
                    c_8[k]=Para_s->d_wkc[k *jmt *imt+j *imt+i] *d_odzt_s[k+1] *d_odzp_s[k] *Para_s->d_dts *0.5e0;
                    b_8[k]=1.0+c_8[k];
                    d_8[k]=d_vtl_s[k *jmt *imt+j *imt+i];
                    e_8[k]=0.0;
                    f_8[k]=0.0;
                    b_8[kz-1]=1.0+a_8[kz-1];
                    c_8[kz-1]=d_odzp_s[kz-1] *Para_s->d_dts *0.5e0;
                    e_8[kz]=0.0;
                    f_8[kz]=0.0;

                    for(k=kz-1;k>=0; k--)
                    {
                        g0=1.0 /(b_8[k]-c_8[k] *e_8[k+1]);
                        e_8[k]=a_8[k] *g0;
                        f_8[k]=(d_8[k]+c_8[k] *f_8[k+1]) *g0;
                    }
                    d_vtl_s[j *imt+i]=(e_8[0] *Para_s->d_stf[j *imt+i]+f_8[0]) *Para_s->d_vit[j *imt+i];

                    for(k=1;k<kz; k++)
                    {
                        d_vtl_s[k *jmt *imt+j *imt+i]=
                                (e_8[k] *d_vtl_s[(k-1) *jmt *imt+j *imt+i]+f_8[k]) *
                                Para_s->d_vit[k *jmt *imt+j *imt+i];
                    }
                }
            }

        }

    }

    ldm_free(d_odzt_s, sizeof(double) *km);

    ldm_free(d_odzp_s, sizeof(double) *km);

    ldm_free(Para_s, sizeof(invtrit_Para));
}
