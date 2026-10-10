
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
//__userDef__ paraVarAttr(double,d_atb,sizeof(double)*ntra*(km+1)*jmt*imt,readonly,no,no,no);	
//__userDef__ paraVarAttr(int,d_kmtn,sizeof(int)*ny_block*nx_block,readonly,no,no,no);	
//__userDef__ paraVarAttr(int,d_kmts,sizeof(int)*ny_block*nx_block,readonly,no,no,no);	
//__userDef__ paraVarAttr(int,d_kmte,sizeof(int)*ny_block*nx_block,readonly,no,no,no);	
//__userDef__ paraVarAttr(int,d_kmtw,sizeof(int)*ny_block*nx_block,readonly,no,no,no);	
//__userDef__ paraVarAttr(double,d_dtn,sizeof(double)*ny_block*nx_block,readonly,no,no,no);	
//__userDef__ paraVarAttr(double,d_dts,sizeof(double)*ny_block*nx_block,readonly,no,no,no);	
//__userDef__ paraVarAttr(double,d_dte,sizeof(double)*ny_block*nx_block,readonly,no,no,no);	
//__userDef__ paraVarAttr(double,d_dtw,sizeof(double)*ny_block*nx_block,readonly,no,no,no);	
//__userDef__ paraVarAttr(int,d_kmt,sizeof(int)*ny_block*nx_block,readonly,no,no,no);	
//__userDef__ paraVarAttr(double,d_d2tk,sizeof(double)*km*jmt*imt,outonly,yes,no,no);	
//__userDef__ paraVarAttr(double,d_ahf,sizeof(double)*ny_block*nx_block,readonly,no,yes,no);	
typedef struct tracer_del4_1_kernelPara{
	int cuBlock[3];
	int cuGrid[3];
	int n;
	int *d_kmtn;
	int *d_kmts;
	int *d_kmte;
	int *d_kmtw;
	int *d_kmt;
	double *d_atb;
	double *d_dtn;
	double *d_dts;
	double *d_dte;
	double *d_dtw;
	double *d_d2tk;
	double *d_ahf;
}tracer_del4_1_Para;


//dynamic atomic definition
#define CP_CORE      64
#define RCV_CORE     0
#define DATAPAKAGE   1

//ATHREAD global definition
__thread_local int blockIndex;

void tracer_del4_1_(tracer_del4_1_Para *Paras){	

  //ATHREAD local variable definition
  int swGrid[3];
  int swBlock[3];
  int stride;
  int blockId_x,blockId_y,blockId_z;
  int threadId_x,threadId_y,threadId_z;
  
  //ATHREAD basic operation
    int i,j,k;
  blockIndex = CRTS_tid;
  tracer_del4_1_Para *Para_s = (tracer_del4_1_Para*)ldm_malloc(sizeof(tracer_del4_1_Para));
  CRTS_dma_iget(Para_s, Paras, sizeof(tracer_del4_1_Para), &dma_rply);
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

  double *d_d2tk_s = Para_s->d_d2tk;

  for(int block = blockStart; block < blockEnd; block++){
    blockId_z = block / (Para_s->cuGrid[0] * Para_s->cuGrid[1]);
    blockId_y = (block - blockId_z * (Para_s->cuGrid[0] * Para_s->cuGrid[1])) / Para_s->cuGrid[0];
    blockId_x = block % Para_s->cuGrid[0];
      j = blockId_y * Para_s->cuBlock[1];
      k = blockId_z * Para_s->cuBlock[2];

    //ATHREAD serial operation replace GPU parallel mutli-thread in single thread-block
    for(i = 0; i < threadSum; i++) {
        
      if(k>=0&&k<km)	
      {	
        if(j>=1&&j<jmt-1&&i>=1&&i<imt-1)	
        {	
          double hdtk=0.0e0;	
          double cc=0.0e0;
          double cn=0.0e0;
          double cs=0.0e0;
          double ce=0.0e0;
          double cw=0.0e0;
          if(k<=Para_s->d_kmt[j *imt+i]-1)	
          {	
            double c;

              if(k<=Para_s->d_kmtn[j *imt+i]-1)
            {	
      	      cn=Para_s->d_dtn[j *imt+i];
      	      //hdtk+=c *Para_s->d_atb[Para_s->n *(km+1) *jmt *imt+(k+1) *jmt *imt+(j-1) *imt+i];
      	      //cc-=c;
            }	
            if(k<=Para_s->d_kmts[j *imt+i]-1)	
            {	
      	      cs=Para_s->d_dts[j *imt+i];
      	      //hdtk+=c *Para_s->d_atb[Para_s->n *(km+1) *jmt *imt+(k+1) *jmt *imt+(j+1) *imt+i];
      	      //cc-=c;
            }	
            if(k<=Para_s->d_kmte[j *imt+i]-1)	
            {	
      	      ce=Para_s->d_dte[j *imt+i];
      	      //hdtk+=c *Para_s->d_atb[Para_s->n *(km+1) *jmt *imt+(k+1) *jmt *imt+j *imt+i+1];
      	      //cc-=c;
            }	
            if(k<=Para_s->d_kmtw[j *imt+i]-1)	
            {	
      	      cw=Para_s->d_dtw[j *imt+i];
      	      //hdtk+=c *Para_s->d_atb[Para_s->n *(km+1) *jmt *imt+(k+1) *jmt *imt+j *imt+i-1];
      	      //cc-=c;
            }
            cc = -(cn + cs + ce + cw);
          }	
      	  hdtk=Para_s->d_ahf[j *imt+i] *(cc *Para_s->d_atb[Para_s->n *(km+1) *jmt *imt+(k+1) *jmt *imt+j *imt+i]+
                  cn *Para_s->d_atb[Para_s->n *(km+1) *jmt *imt+(k+1) *jmt *imt+(j-1) *imt+i] +
                  cs *Para_s->d_atb[Para_s->n *(km+1) *jmt *imt+(k+1) *jmt *imt+(j+1) *imt+i] +
                  ce *Para_s->d_atb[Para_s->n *(km+1) *jmt *imt+(k+1) *jmt *imt+j *imt+i+1] +
                  cw *Para_s->d_atb[Para_s->n *(km+1) *jmt *imt+(k+1) *jmt *imt+j *imt+i-1]);

      	  d_d2tk_s[k *jmt *imt+j *imt+i]=hdtk;	
        }	
      }	

    }

  }

    ldm_free(Para_s, sizeof(tracer_del4_1_Para));
}
