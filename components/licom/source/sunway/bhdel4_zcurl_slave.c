
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

//__userDef__ paraVarAttr(double,d_dxu,sizeof(double)*ny_block*nx_block,readonly,no,no,no);	
//__userDef__ paraVarAttr(double,d_dyu,sizeof(double)*ny_block*nx_block,readonly,no,no,no);	
//__userDef__ paraVarAttr(double,d_ubp,sizeof(double)*jmt*imt,readonly,no,no,no);	
//__userDef__ paraVarAttr(double,d_vbp,sizeof(double)*jmt*imt,readonly,no,no,no);	
//__userDef__ paraVarAttr(double,d_wka,sizeof(double)*jmt*imt*km,outonly,yes,no,no);	
//__userDef__ paraVarAttr(int,d_kmt,sizeof(int)*ny_block*nx_block,readonly,no,no,no);	
typedef struct bhdel4_zcurl_kernelPara{
	int cuBlock[3];
	int cuGrid[3];
	int rank;
	int *d_kmt;
	double *d_dxu;
	double *d_dyu;
	double *d_ubp;
	double *d_vbp;
	double *d_wka;
	double *d_tarea_r;
	double *d_z_curl;
}bhdel4_zcurl_Para;


//dynamic atomic definition
#define CP_CORE      64
#define RCV_CORE     0
#define DATAPAKAGE   1

//ATHREAD global definition 
__thread_local int warp_size = 32; 
__thread_local int blockIndex;

void bhdel4_zcurl_(bhdel4_zcurl_Para *Paras){	

  //ATHREAD local variable definition
  int swGrid[3];
  int swBlock[3];
  int stride;
  int blockId_x,blockId_y,blockId_z;
  int threadId_x,threadId_y,threadId_z;
  double d_wka[imt];
  double d_ubp[imt];
  double d_ubpmj1[imt];
  double d_vbp[imt];
  double d_vbpmj1[imt];
  double d_dxu[imt];
  double d_dxumj1[imt];
  double d_dyu[imt];
  double d_dyumj1[imt];

  //ATHREAD basic operation
    int i,j,k;
  blockIndex = CRTS_tid;
  bhdel4_zcurl_Para *Para_s = (bhdel4_zcurl_Para*)ldm_malloc(sizeof(bhdel4_zcurl_Para));
  CRTS_dma_iget(Para_s, Paras, sizeof(bhdel4_zcurl_Para), &dma_rply);
  D_COUNT++;
  CRTS_dma_wait_value(&dma_rply, D_COUNT);

  double *d_wka_s = Para_s->d_wka;

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

        //CRTS_dma_get(d_wka, &Para_s->d_wka[j * imt], sizeof(double) * imt);
        CRTS_dma_get(d_ubp, &Para_s->d_ubp[j * imt], sizeof(double) * imt);
        CRTS_dma_get(d_ubpmj1, &Para_s->d_ubp[(j - 1) * imt], sizeof(double) * imt);
        CRTS_dma_get(d_vbp, &Para_s->d_vbp[j * imt], sizeof(double) * imt);
        CRTS_dma_get(d_vbpmj1, &Para_s->d_vbp[(j - 1) * imt], sizeof(double) * imt);
        CRTS_dma_get(d_dxu, &Para_s->d_dxu[j * imt], sizeof(double) * imt);
        CRTS_dma_get(d_dxumj1, &Para_s->d_dxu[(j - 1) * imt], sizeof(double) * imt);
        CRTS_dma_get(d_dyu, &Para_s->d_dyu[j * imt], sizeof(double) * imt);
        CRTS_dma_get(d_dyumj1, &Para_s->d_dyu[(j - 1) * imt], sizeof(double) * imt);

        for (i = 0; i < imt-1; i++) {
         //   d_wka[i]=0.0e0;
            if(j>=1){
                if(0<=Para_s->d_kmt[j * imt + i]-1)
                {
                    Para_s->d_z_curl[j*imt+i] = 0.5 * ( d_vbp[i+1] * d_dyu[i+1] +
							d_vbpmj1[i+1] * d_dyumj1[i+1] -
                            d_vbp[i] * d_dyu[i]-
							d_vbpmj1[i] * d_dyumj1[i] -
                            d_ubp[i+1] * d_dxu[i+1] -
							d_ubp[i] * d_dxu[i] +
                            d_ubpmj1[i+1] * d_dxumj1[i+1] +
							d_ubpmj1[i] * d_dxumj1[i]) *
                            Para_s->d_tarea_r[j * imt + i];
                }
            }
        }

        //CRTS_dma_put(&Para_s->d_wka[j * imt], d_wka, sizeof(double) * imt);
    }


    ldm_free(Para_s, sizeof(bhdel4_zcurl_Para));
}
