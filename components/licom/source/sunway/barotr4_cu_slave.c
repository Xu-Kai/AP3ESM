
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

//__userDef__ paraVarAttr(double,d_dtb,sizeof(double),readonly,no,no,no);	
//__userDef__ paraVarAttr(double,d_h0,sizeof(double)*jmt*imt,readonly,no,no,no);	
//__userDef__ paraVarAttr(double,d_h0p,sizeof(double)*jmt*imt,readonly,no,no,no);	
//__userDef__ paraVarAttr(double,d_ub,sizeof(double)*jmt*imt,readonly,no,no,no);	
//__userDef__ paraVarAttr(double,d_vb,sizeof(double)*jmt*imt,readonly,no,no,no);	
//__userDef__ paraVarAttr(double,d_ubp,sizeof(double)*jmt*imt,outonly,yes,no,no);	
//__userDef__ paraVarAttr(double,d_vbp,sizeof(double)*jmt*imt,outonly,yes,no,no);	
//__userDef__ paraVarAttr(double,d_h0f,sizeof(double)*jmt*imt,outonly,yes,no,no);	
//__userDef__ paraVarAttr(double,d_h0bf,sizeof(double)*jmt*imt,outonly,yes,no,no);	
//__userDef__ paraVarAttr(double,d_work,sizeof(double)*jmt*imt,readonly,no,no,no);	
//__userDef__ paraVarAttr(int,d_jmt,sizeof(int),readonly,no,no,no);	
//__userDef__ paraVarAttr(int,d_imt,sizeof(int),readonly,no,no,no);	
typedef struct barotr4_cu_kernelPara{
    int cuBlock[3];
    int cuGrid[3];
	int d_imt;
	int d_jmt;
	bool d_iau;
	double d_dtb;
	double d_lamda;
	double *d_h0;
	double *d_h0p;
	double *d_ub;
	double *d_vb;
	double *d_ubp;
	double *d_vbp;
	double *d_h0f;
	double *d_h0bf;
	double *d_work;
	double *d_h0_iau;
	double *d_vit;
}barotr4_cu_Para;


//dynamic atomic definition
#define CP_CORE      64
#define RCV_CORE     0
#define DATAPAKAGE   1

//ATHREAD global definition 
__thread_local int warp_size = 32; 
__thread_local int blockIndex;

void barotr4_cu_(barotr4_cu_Para *Paras){

//    double d_h0p[imt];

//    double d_h0f[imt];
//    double d_h0bf[imt];
//    double d_work[imt];

	double d_ubp[imt];
	double d_vbp[imt];
    double d_h0[imt];

	double h0;
	double *d_h0p;
	double *d_h0f;
	double *d_h0bf;
	double *d_work;
  //ATHREAD basic operation
    int blockId_x,blockId_y,blockId_z;
    int i,j,k;
  blockIndex = CRTS_tid;
  barotr4_cu_Para *Para_s = (barotr4_cu_Para*)ldm_malloc(sizeof(barotr4_cu_Para));
  CRTS_dma_iget(Para_s, Paras, sizeof(barotr4_cu_Para), &dma_rply);
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

//        CRTS_dma_get(d_h0p, &Para_s->d_h0p[j * imt], sizeof(double) * imt);
//        CRTS_dma_get(d_h0f, &Para_s->d_h0f[j * imt], sizeof(double) * imt);
//        CRTS_dma_get(d_h0bf, &Para_s->d_h0bf[j * imt], sizeof(double) * imt);
//        CRTS_dma_get(d_work, &Para_s->d_work[j * imt], sizeof(double) * imt);

    	d_h0p  = &Para_s->d_h0p[j * imt];
    	d_h0f  = &Para_s->d_h0f[j * imt];
    	d_h0bf = &Para_s->d_h0bf[j * imt];
    	d_work = &Para_s->d_work[j * imt];

        for (i = 0; i < imt; i++) {

			if(Para_s->d_iau)
                h0 = d_h0p[i] + (d_work[i] * Para_s->d_dtb) +
						Para_s->d_vit[j*imt+i] * Para_s->d_h0_iau[j*imt+i] * Para_s->d_lamda;
			else
				h0 = d_h0p[i] + (d_work[i] * Para_s->d_dtb);

            d_h0p[i] = h0;
            d_ubp[i] = Para_s->d_ub[j * imt + i];
            d_vbp[i] = Para_s->d_vb[j * imt + i];
            d_h0f[i] += h0;
            d_h0bf[i] += h0;
            d_h0[i] = h0;
        }

//        CRTS_dma_put(&Para_s->d_h0p[j * imt], d_h0p, sizeof(double) * imt);
        CRTS_dma_put(&Para_s->d_ubp[j * imt], d_ubp, sizeof(double) * imt);
        CRTS_dma_put(&Para_s->d_vbp[j * imt], d_vbp, sizeof(double) * imt);
//        CRTS_dma_put(&Para_s->d_h0f[j * imt], d_h0f, sizeof(double) * imt);
//        CRTS_dma_put(&Para_s->d_h0bf[j * imt], d_h0bf, sizeof(double) * imt);
        CRTS_dma_put(&Para_s->d_h0[j * imt], d_h0, sizeof(double) * imt);
    }
      

  

    ldm_free(Para_s, sizeof(barotr4_cu_Para));
}
