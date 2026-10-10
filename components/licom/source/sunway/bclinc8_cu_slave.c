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

typedef struct bclinc81Para{
    int cuBlock[3];
	int cuGrid[3];
	double *d_viv;
	double *d_dzp;
	double *d_u1;
	double *d_u;
	double *d_ub;
	double *d_utf;
	double *d_ohbu;
	double *d_work;
}bclinc81_Para;

typedef struct bclinc82Para{
	int cuBlock[3];
	int cuGrid[3];
	double *d_viv;
	double *d_dzp;
	double *d_v1;
	double *d_v;
	double *d_vb;
	double *d_vtf;
	double *d_ohbu;
	double *d_work;
}bclinc82_Para;

//dynamic atomic definition
#define CP_CORE      64
#define RCV_CORE     0
#define DATAPAKAGE   1

//ATHREAD global definition 
__thread_local int warp_size = 32; 
__thread_local int blockIndex;
//__thread_local Para* Para_s;

 void  bclinc81_cu_(bclinc81_Para *Paras){

  //ATHREAD local variable definition
     int i,j,k;
  int swGrid[3];
  int swBlock[3];
  int stride;
  int blockId_x,blockId_y,blockId_z;
  int threadId_x,threadId_y,threadId_z;


  //ATHREAD basic operation
  blockIndex = CRTS_tid;
  bclinc81_Para *Para_s = (bclinc81_Para*)ldm_malloc(sizeof(bclinc81_Para));
  CRTS_dma_get(Para_s, Paras, sizeof(bclinc81_Para));

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

  double* d_utf = Para_s->d_utf;
  double* d_ub = Para_s->d_ub;
  double* d_u = Para_s->d_u;

  double work=0.0e0,u=0.0e0;

  for(int block = blockStart; block < blockEnd; block++){
    blockId_z = block / (Para_s->cuGrid[0] * Para_s->cuGrid[1]);
    blockId_y = (block - blockId_z * (Para_s->cuGrid[0] * Para_s->cuGrid[1])) / Para_s->cuGrid[0];
    blockId_x = block % Para_s->cuGrid[0];
      j = blockId_y * Para_s->cuBlock[1];
      k = blockId_z * Para_s->cuBlock[2];

    //ATHREAD serial operation replace GPU parallel mutli-thread in single thread-block
    for(i = 0; i < threadSum; i++) {

      if(i<imt&&j<jmt&&k<km)
      {
	u = (d_u[k *jmt *imt+j *imt+i]-Para_s->d_work[j * imt + i]+d_ub[j *imt+i]) *Para_s->d_viv[k *jmt *imt+j *imt+i];
//	d_utf[k *jmt *imt+j *imt+i]+=u;
	d_u[k *jmt *imt+j *imt+i]= u;
      }	

    }

  }

    ldm_free(Para_s, sizeof(bclinc81_Para));
}

 void  bclinc82_cu_(bclinc82_Para *Paras){	

  //ATHREAD local variable definition
  int swGrid[3];
  int swBlock[3];
  int stride;
  int blockId_x,blockId_y,blockId_z;
  int threadId_x,threadId_y,threadId_z;


  //ATHREAD basic operation
     int i,j,k;
     blockIndex = CRTS_tid;
  bclinc82_Para* Para_s = (bclinc82_Para*)ldm_malloc(sizeof(bclinc82_Para));
  CRTS_dma_get(Para_s, Paras, sizeof(bclinc82_Para));
  //D_COUNT++;
  //CRTS_dma_wait_value(&dma_rply, D_COUNT);

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

  //double *d_wp12_s = Para_s->d_wp12;
  //double *d_wp13_s = Para_s->d_wp13;
  double* d_vtf = Para_s->d_vtf;
  double* d_v1 = Para_s->d_v1;
  double* d_v =  Para_s->d_v;
  double* d_vb = Para_s->d_vb;

  for(int block = blockStart; block < blockEnd; block++){
    blockId_z = block / (Para_s->cuGrid[0] * Para_s->cuGrid[1]);
    blockId_y = (block - blockId_z * (Para_s->cuGrid[0] * Para_s->cuGrid[1])) / Para_s->cuGrid[0];
    blockId_x = block % Para_s->cuGrid[0];
      j = blockId_y * Para_s->cuBlock[1];
      k = blockId_z * Para_s->cuBlock[2];

    //ATHREAD serial operation replace GPU parallel mutli-thread in single thread-block
    for(i = 0; i < threadSum; i++) {

      if(i<imt&&j<jmt&&k<km)
      {	
        double v = (d_v[k *jmt *imt+j *imt+i]-Para_s->d_work[j * imt + i] +d_vb[j *imt+i]) *Para_s->d_viv[k *jmt *imt+j *imt+i];
//        d_vtf[k *jmt *imt+j *imt+i]+=v;
	Para_s->d_v[k *jmt *imt+j *imt+i]=v;
      }	

    }

  }

    ldm_free(Para_s, sizeof(bclinc82_Para));
}

