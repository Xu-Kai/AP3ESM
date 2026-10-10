
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

//__userDef__ paraVarAttr(double,d_work,sizeof(double)*jmt*imt,outonly,no,yes,no);	
//__userDef__ paraVarAttr(double,d_wka,sizeof(double)*jmt*imt*km,readonly,no,no,no);	
//__userDef__ paraVarAttr(double,d_vit,sizeof(double)*jmt*imt*km,readonly,no,no,no);	
//__userDef__ paraVarAttr(double,d_htw,sizeof(double)*ny_block*nx_block,readonly,no,no,no);	
//__userDef__ paraVarAttr(double,d_hts,sizeof(double)*ny_block*nx_block,readonly,no,no,no);	
//__userDef__ paraVarAttr(int,d_kmt,sizeof(int)*ny_block*nx_block,readonly,no,no,no);	
//__userDef__ paraVarAttr(double,d_tarea_r,sizeof(double)*ny_block*nx_block,readonly,no,no,no);	
//__userDef__ paraVarAttr(int,d_jmt,sizeof(int),readonly,no,no,no);	
//__userDef__ paraVarAttr(int,d_imt,sizeof(int),readonly,no,no,no);	
typedef struct barotr1_cu1_kernelPara{
	int cuBlock[3];
	int cuGrid[3];
	int d_imt;
	int d_jmt;
	int rank;
	int *d_kmt;
	double *d_work;
	double *d_wka;
	double *d_vit;
	double *d_htw;
	double *d_hts;
	double *d_tarea_r;
}barotr1_cu1_Para;


//dynamic atomic definition
#define CP_CORE      64
#define RCV_CORE     0
#define DATAPAKAGE   1

//ATHREAD global definition 
__thread_local int warp_size = 32; 
__thread_local int blockIndex;

void barotr1_cu1_(barotr1_cu1_Para *Paras){	

  //ATHREAD local variable definition
  int swGrid[3];
  int swBlock[3];
  int stride;
  double d_work[imt];
  double d_wka[2*imt];
  double d_wka1[2*imt];
  double d_hts[2*imt];
  double d_htw[imt];
  double d_tarea_r[imt];
  double d_vit[imt];

  //ATHREAD basic operation
    int i,j,k;
    int blockId_x,blockId_y,blockId_z;
  blockIndex = CRTS_tid;
  barotr1_cu1_Para *Para_s = (barotr1_cu1_Para*)ldm_malloc(sizeof(barotr1_cu1_Para));
  CRTS_dma_iget(Para_s, Paras, sizeof(barotr1_cu1_Para), &dma_rply);
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

        if(j>0 && j < jmt) {
            CRTS_dma_get(d_wka, &Para_s->d_wka[(j-1) * imt], sizeof(double) * imt * 2);
            CRTS_dma_get(d_wka1, &Para_s->d_wka[imt * jmt + (j-1) * imt], sizeof(double) * imt * 2);
            CRTS_dma_get(d_work, &Para_s->d_work[j * imt], sizeof(double) * imt);
            CRTS_dma_get(d_hts, &Para_s->d_hts[(j-1) * imt], sizeof(double) * imt * 2);
            CRTS_dma_get(d_htw, &Para_s->d_htw[j * imt], sizeof(double) * imt);
            CRTS_dma_get(d_tarea_r, &Para_s->d_tarea_r[j * imt], sizeof(double) * imt);
            CRTS_dma_get(d_vit, &Para_s->d_vit[j * imt], sizeof(double) * imt);
        }

        for (i = 0; i < imt; i++) {

            if(i<imt-2&&i> 1&&j> 0&&j<jmt)
            {
                if(0<=Para_s->d_kmt[j *imt+i]-1)
                {
                    d_work[i] = 0.5 *((d_wka[imt+ i+1]+d_wka[i+1]) * d_htw[i+1]-(d_wka[imt+ i]+d_wka[i]) *
                            d_htw[i]+(d_wka1[imt+i+1]+ d_wka1[imt+i]) * d_hts[imt+i]-(d_wka1[i+1]+d_wka1[i]) *
                            d_hts[i]) * d_tarea_r[i] *d_vit[i] *-0.250;
                }
                else{
                    d_work[i]=0;
                }
            }

        }
        if(j>0 && j < jmt)
            CRTS_dma_put(&Para_s->d_work[j * imt], d_work, sizeof(double) * imt);
    }

    ldm_free(Para_s, sizeof(barotr1_cu1_Para));
}
