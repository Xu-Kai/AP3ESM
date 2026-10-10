
//ATHREAD include file
#include <crts.h>
#include <slave.h>
//#include <simd.h>
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

//__userDef__ paraVarAttr(double,d_dlu,sizeof(double)*jmt*imt*km,readonly,no,no,no);	
//__userDef__ paraVarAttr(double,d_dlv,sizeof(double)*jmt*imt*km,outonly,yes,no,no);	
//__userDef__ paraVarAttr(double,d_ohbt,sizeof(double)*jmt*imt,readonly,no,no,no);	
//__userDef__ paraVarAttr(double,d_g,sizeof(double),readonly,no,no,no);	
//__userDef__ paraVarAttr(int,d_jmt,sizeof(int),readonly,no,no,no);	
//__userDef__ paraVarAttr(int,d_imt,sizeof(int),readonly,no,no,no);	
typedef struct readyt_4_kernelPara{
	int cuBlock[3];
	int cuGrid[3];
	int d_km;
	int d_jmt;
	int d_imt;
	int rank;
	double d_g;
	double *d_dlu;
	double *d_dlv;
	double *d_ohbt;
}readyt_4_para;


//dynamic atomic definition
#define CP_CORE      64
#define RCV_CORE     0
#define DATAPAKAGE   1

//ATHREAD global definition
__thread_local int blockIndex;
//__thread_local Para* Para_s;

void  readyt_cu_4_(readyt_4_para *Paras){

  //ATHREAD local variable definition
  int swGrid[3];
  int swBlock[3];
  int stride;
  int blockId_x,blockId_y,blockId_z;
  int threadId_x,threadId_y,threadId_z;


  //ATHREAD basic operation
  int i,j;
  blockIndex = CRTS_tid;
  readyt_4_para* Para_s = (readyt_4_para*)ldm_malloc(sizeof(readyt_4_para));
  CRTS_dma_iget(Para_s, Paras, sizeof(readyt_4_para), &dma_rply);
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

  double *d_dlv_s = Para_s->d_dlv;

  for(int block = blockStart; block < blockEnd; block++){
    blockId_z = block / (Para_s->cuGrid[0] * Para_s->cuGrid[1]);
    blockId_y = (block - blockId_z * (Para_s->cuGrid[0] * Para_s->cuGrid[1])) / Para_s->cuGrid[0];
    blockId_x = block % Para_s->cuGrid[0];
      j = blockId_y * Para_s->cuBlock[1];

    //ATHREAD serial operation replace GPU parallel mutli-thread in single thread-block
    for(i = 0; i < threadSum; i++) {

      if(i<Para_s->d_imt&&j<Para_s->d_jmt)	
      {	
      	d_dlv_s[j *Para_s->d_imt+i]=(Para_s->d_dlu[j *Para_s->d_imt+i]+Para_s->d_dlu[1 *Para_s->d_jmt *Para_s->d_imt+j *Para_s->d_imt+i] *Para_s->d_ohbt[j *Para_s->d_imt+i]) /Para_s->d_g;	
      	d_dlv_s[1 *Para_s->d_jmt *Para_s->d_imt+j *Para_s->d_imt+i]=Para_s->d_dlu[1 *Para_s->d_jmt *Para_s->d_imt+j *Para_s->d_imt+i] *Para_s->d_ohbt[j *Para_s->d_imt+i] *Para_s->d_ohbt[j *Para_s->d_imt+i];	
      }	

    }

  }

    ldm_free(Para_s, sizeof(readyt_4_para));
}
