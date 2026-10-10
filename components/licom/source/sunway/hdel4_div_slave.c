
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

//__userDef__ paraVarAttr(double,d_hts,sizeof(double)*ny_block*nx_block,readonly,no,no,no);	
//__userDef__ paraVarAttr(double,d_htw,sizeof(double)*ny_block*nx_block,readonly,no,no,no);	
//__userDef__ paraVarAttr(double,d_tarea_r,sizeof(double)*ny_block*nx_block,readonly,no,no,no);	
//__userDef__ paraVarAttr(double,d_up,sizeof(double)*km*jmt*imt,readonly,no,no,no);	
//__userDef__ paraVarAttr(double,d_vp,sizeof(double)*km*jmt*imt,readonly,no,no,no);	
//__userDef__ paraVarAttr(double,d_wka,sizeof(double)*jmt*imt*km,outonly,yes,no,no);	
//__userDef__ paraVarAttr(int,d_kmt,sizeof(int)*ny_block*nx_block,readonly,no,no,no);	
typedef struct hdel4_div_kernelPara{
	int cuBlock[3];
	int cuGrid[3];
	int rank;
	int *d_kmt;
	double *d_hts;
	double *d_htw;
	double *d_tarea_r;
	double *d_up;
	double *d_vp;
	double *d_wka;
}hdel4_div_Para;

//dynamic atomic definition
#define CP_CORE      64
#define RCV_CORE     0
#define DATAPAKAGE   1

//ATHREAD global definition 
__thread_local int warp_size = 32; 
__thread_local int blockIndex;

void hdel4_div_(hdel4_div_Para *Paras){	

  //ATHREAD local variable definition
  int swGrid[3];
  int swBlock[3];
  int stride;
  int blockId_x,blockId_y,blockId_z;
  int threadId_x,threadId_y,threadId_z;


  //ATHREAD basic operation
    int i,j,k;
  blockIndex = CRTS_tid;
  hdel4_div_Para *Para_s = (hdel4_div_Para*)ldm_malloc(sizeof(hdel4_div_Para));
  CRTS_dma_iget(Para_s, Paras, sizeof(hdel4_div_Para), &dma_rply);
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

  double *d_wka_s = Para_s->d_wka;

  for(int block = blockStart; block < blockEnd; block++){
    blockId_z = block / (Para_s->cuGrid[0] * Para_s->cuGrid[1]);
    blockId_y = (block - blockId_z * (Para_s->cuGrid[0] * Para_s->cuGrid[1])) / Para_s->cuGrid[0];
    blockId_x = block % Para_s->cuGrid[0];
      j = blockId_y * Para_s->cuBlock[1];
      k = blockId_z * Para_s->cuBlock[2];

    //ATHREAD serial operation replace GPU parallel mutli-thread in single thread-block
    for(i = 0; i < threadSum; i++) {
      unsigned int tid=j *imt+i;	
	
      if(i<imt&&j<jmt&&k<km)	
      {	
      	d_wka_s[k *jmt *imt+tid]=0.0e0;	
        if(i<imt-1&&j>=1)	
        {	
          if(k<=Para_s->d_kmt[tid]-1)	
          {	
      	    d_wka_s[k *jmt *imt+tid]=0.5 *((Para_s->d_up[k *jmt *imt+tid+1]+
					Para_s->d_up[k *jmt *imt+tid-imt+1]) *Para_s->d_htw[tid+1]-
					(Para_s->d_up[k *jmt *imt+tid]+Para_s->d_up[k *jmt *imt+tid-imt]) *Para_s->d_htw[tid]+
					(Para_s->d_vp[k *jmt *imt+tid+1]+Para_s->d_vp[k *jmt *imt+tid]) *Para_s->d_hts[tid]-
					(Para_s->d_vp[k *jmt *imt+tid-imt+1]+Para_s->d_vp[k *jmt *imt+tid-imt]) *Para_s->d_hts[tid-imt]) *Para_s->d_tarea_r[tid];	
          }	
        }	
      }	
    }

  }

    ldm_free(Para_s, sizeof(hdel4_div_Para));
}
