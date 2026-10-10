
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
//__userDef__ paraVarAttr(double,d_dzph,sizeof(double)*jmt*imt,readonly,no,no,no);	
//__userDef__ paraVarAttr(double,d_ub,sizeof(double)*jmt*imt,outonly,yes,no,no);	
//__userDef__ paraVarAttr(double,d_vb,sizeof(double)*jmt*imt,outonly,yes,no,no);	
//__userDef__ paraVarAttr(double,d_ubp,sizeof(double)*jmt*imt,readonly,no,no,no);	
//__userDef__ paraVarAttr(double,d_vbp,sizeof(double)*jmt*imt,readonly,no,no,no);	
//__userDef__ paraVarAttr(double,d_wka,sizeof(double)*km*jmt*imt,outonly,yes,no,no);	
//__userDef__ paraVarAttr(double,d_work,sizeof(double)*jmt*imt,readonly,no,no,no);	
//__userDef__ paraVarAttr(int,d_jmt,sizeof(int),readonly,no,no,no);	
//__userDef__ paraVarAttr(int,d_imt,sizeof(int),readonly,no,no,no);	
typedef struct barotr3_cu_kernelPara{
	int cuBlock[3];
	int cuGrid[3];
	int d_imt;
	int d_jmt;
	int rank;
	int pad;
	double d_dtb;
	double *d_dzph;
	double *d_ub;
	double *d_vb;
	double *d_ubp;
	double *d_vbp;
	double *d_wka;
	double *d_work;
}barotr3_cu_Para;


//dynamic atomic definition
#define CP_CORE      64
#define RCV_CORE     0
#define DATAPAKAGE   1

//ATHREAD global definition 
__thread_local int warp_size = 32; 
__thread_local int blockIndex;

void barotr3_cu_(barotr3_cu_Para *Paras){

  double d_ub[imt];
  double d_vb[imt];
  double d_wka0[imt];
  double d_wka1[imt];
  double *d_ubp;
  double *d_vbp;
  double *d_wka2;
  double *d_wka3;
  double *d_dzph;
  double *d_work;

  //ATHREAD basic operation
    int blockId_x,blockId_y,blockId_z;
    int i,j,k;
  blockIndex = CRTS_tid;
  barotr3_cu_Para *Para_s = (barotr3_cu_Para*)ldm_malloc(sizeof(barotr3_cu_Para));
  CRTS_dma_iget(Para_s, Paras, sizeof(barotr3_cu_Para), &dma_rply);
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

        d_ubp  = &Para_s->d_ubp[j * imt];
        d_vbp  = &Para_s->d_vbp[j * imt];
        d_wka2 = &Para_s->d_wka[2 *jmt *imt+j *imt];
        d_wka3 = &Para_s->d_wka[3 *jmt *imt+j *imt];
        d_dzph = &Para_s->d_dzph[j * imt];
        d_work = &Para_s->d_work[j * imt];


        for (i = 0; i < imt; i++) {
            d_ub[i]=d_ubp[i]+ d_wka2[i] *Para_s->d_dtb;
            d_vb[i]=d_vbp[i]+d_wka3[i] *Para_s->d_dtb;

            if(j>=1&&i<imt-1){
                d_wka0[i]=d_ub[i] * (d_dzph[i]+d_work[i]);
                d_wka1[i]=d_vb[i] * (d_dzph[i]+d_work[i]);
            }
        }

        CRTS_dma_put(&Para_s->d_ub[j * imt], d_ub, sizeof(double) * imt);
        CRTS_dma_put(&Para_s->d_vb[j * imt], d_vb, sizeof(double) * imt);
        CRTS_dma_put(&Para_s->d_wka[j * imt], d_wka0, sizeof(double) * imt);
        CRTS_dma_put(&Para_s->d_wka[1 *jmt *imt+j *imt], d_wka1, sizeof(double) * imt);
    }

    ldm_free(Para_s, sizeof(barotr3_cu_Para));
}
