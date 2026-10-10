
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

//__userDef__ paraVarAttr(int,n,sizeof(int),readonly,no,no,no);	
//__userDef__ paraVarAttr(double,d_at,sizeof(double)*ntra*km*jmt*imt,outonly,yes,no,no);	
//__userDef__ paraVarAttr(double,d_vtl,sizeof(double)*km*jmt*imt,readonly,no,no,no);	
//__userDef__ paraVarAttr(double,d_atb,sizeof(double)*ntra*(km+1)*jmt*imt,outonly,yes,no,no);	
//__userDef__ paraVarAttr(double,d_tend,sizeof(double)*ntra*km*jmt*imt,outonly,yes,no,no);	
//__userDef__ paraVarAttr(double,d_vit,sizeof(double)*jmt*imt*km,readonly,no,no,no);	
//__userDef__ paraVarAttr(double,d_dt_diff,sizeof(double)*ntra*km*jmt*imt,outonly,yes,no,no);	
//__userDef__ paraVarAttr(double,d_ori,sizeof(double)*km*jmt*imt,readonly,no,no,no);	
//__userDef__ paraVarAttr(double,d_stf,sizeof(double)*jmt*imt,readonly,no,no,no);	
//__userDef__ paraVarAttr(double,d_odzp,sizeof(double)*km,readonly,no,yes,no);	
//__userDef__ paraVarAttr(double,d_dts,sizeof(double),readonly,no,no,no);	
typedef struct tracer_9_kernelPara{
	int cuBlock[3];
	int cuGrid[3];
	int n;
	int rank;
	double d_dts;
	double *d_at;
	double *d_vtl;
	double *d_atb;
	double *d_tend;
	double *d_vit;
	//double *d_dt_diff;
	double *d_ori;
	double *d_stf;
	double *d_odzp;
}tracer_9_Para;


//dynamic atomic definition
#define CP_CORE      64
#define RCV_CORE     0
#define DATAPAKAGE   1

//ATHREAD global definition 
__thread_local int blockIndex;

void tracer_9_(tracer_9_Para *Paras){	

  //ATHREAD basic operation
    int blockId_x,blockId_y,blockId_z;
    int i,j,k;
  blockIndex = CRTS_tid;
  tracer_9_Para *Para_s = (tracer_9_Para*)ldm_malloc(sizeof(tracer_9_Para));
  CRTS_dma_iget(Para_s, Paras, sizeof(tracer_9_Para), &dma_rply);
  D_COUNT++;
  CRTS_dma_wait_value(&dma_rply, D_COUNT);

  double* d_odzp_s = (double*)ldm_malloc(sizeof(double) *km);
  CRTS_dma_iget(d_odzp_s, Para_s->d_odzp, sizeof(double) *km, &dma_rply);
  D_COUNT++;
  CRTS_dma_wait_value(&dma_rply, D_COUNT);

  int n = Para_s->n;
  double d_dt_diff_s[imt],d_vtl[imt],d_ori[imt],d_vit[imt],d_tend_s[imt],d_at_s[imt];

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

 	    CRTS_dma_get(d_vtl, &Para_s->d_vtl[k *jmt *imt+j *imt], sizeof(double) * imt);
	    CRTS_dma_get(d_ori, &Para_s->d_ori[k *jmt *imt+j *imt], sizeof(double) * imt);
	    CRTS_dma_get(d_vit, &Para_s->d_vit[k *jmt *imt+j *imt], sizeof(double) * imt);
	    CRTS_dma_get(d_at_s, &Para_s->d_at[n *km *jmt *imt+k *jmt *imt+j *imt], sizeof(double) * imt);
/*
	    for(i=0;i<imt;i++)
	    {
              //d_dt_diff_s[i]= (d_vtl[i]-d_ori[i]) /Para_s->d_dts *d_vit[i];
  
              if(k==0)	
              {	
              //	d_dt_diff_s[i]=(d_vtl[i]-d_ori[i]) /
               // 		Para_s->d_dts *d_vit[i]-Para_s->d_stf[j *imt+i] *0.5e0 *d_odzp_s[0];	
              }	
      	      d_tend_s[i]= (d_vtl[i]-d_at_s[i]) /Para_s->d_dts *d_vit[i];
        }

	    //CRTS_dma_put(&Para_s->d_dt_diff[n *km *jmt *imt+k *jmt *imt+j *imt], d_dt_diff_s, sizeof(double) * imt);
	    CRTS_dma_put(&Para_s->d_tend[n *km *jmt *imt+k *jmt *imt+j *imt], d_tend_s, sizeof(double) * imt);
*/
	    CRTS_dma_put(&Para_s->d_at[n *km *jmt *imt+k *jmt *imt+j *imt], d_vtl, sizeof(double) * imt);
	    CRTS_dma_put(&Para_s->d_atb[n *(km+1) *jmt *imt+(k+1) *jmt *imt+j *imt], d_vtl, sizeof(double) * imt);

  }
    //todo: user need modify data size according to the practical situation

    ldm_free(d_odzp_s, sizeof(double) *km);

    ldm_free(Para_s, sizeof(tracer_9_Para));
}
