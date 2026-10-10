
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

//__userDef__ paraVarAttr(double,d_ppb,sizeof(double)*km*jmt*imt,readonly,no,no,no);	
//__userDef__ paraVarAttr(double,d_ppc,sizeof(double)*km*jmt*imt,readonly,no,no,no);	
//__userDef__ paraVarAttr(double,d_ppa,sizeof(double)*km*jmt*imt,readonly,no,no,no);	
//__userDef__ paraVarAttr(double,d_alpha,sizeof(double)*km*jmt*imt,outonly,yes,no,no);	
//__userDef__ paraVarAttr(double,d_beta,sizeof(double)*km*jmt*imt,outonly,yes,no,no);	
//__userDef__ paraVarAttr(double,d_vit,sizeof(double)*jmt*imt*km,readonly,no,no,no);	
//__userDef__ paraVarAttr(double,d_od0,sizeof(double),readonly,no,no,no);	
//__userDef__ paraVarAttr(int,d_imt,sizeof(int),readonly,no,no,no);	
//__userDef__ paraVarAttr(int,d_jmt,sizeof(int),readonly,no,no,no);	
//__userDef__ paraVarAttr(int,d_km,sizeof(int),readonly,no,no,no);	
typedef struct thermal_cu_para{
	int cuBlock[3];
	int cuGrid[3];
	int d_km;
	int d_jmt;
	int d_imt;
	int rank;
	double d_od0;
	double *d_ppb;
	double *d_ppc;
	double *d_ppa;
	double *d_alpha;
	double *d_beta;
	double *d_vit;
}thermal_para;


//dynamic atomic definition
#define CP_CORE      64
#define RCV_CORE     0
#define DATAPAKAGE   1

//ATHREAD global definition 
__thread_local int warp_size = 32; 
__thread_local int blockIndex;
//__thread_local Para* Para_s;

void  thermal_cu_(thermal_para *Paras){

  //ATHREAD local variable definition
    int i,j,k;
  int swGrid[3];
  int swBlock[3];
  int stride;
  int blockId_x,blockId_y,blockId_z;
  int threadId_x,threadId_y,threadId_z;
      
  double tt1,tt2,tt3,ss1,ss2,tm1,tm2,tm3;	
  double immediate;	

  //ATHREAD basic operation
  blockIndex = CRTS_tid;
  thermal_para* Para_s = (thermal_para*)ldm_malloc(sizeof(thermal_para));
  CRTS_dma_iget(Para_s, Paras, sizeof(thermal_para), &dma_rply);
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

  double *d_alpha_s = Para_s->d_alpha;
  double *d_beta_s = Para_s->d_beta;

  for(int block = blockStart; block < blockEnd; block++){
    blockId_z = block / (Para_s->cuGrid[0] * Para_s->cuGrid[1]);
    blockId_y = (block - blockId_z * (Para_s->cuGrid[0] * Para_s->cuGrid[1])) / Para_s->cuGrid[0];
    blockId_x = block % Para_s->cuGrid[0];

      j = blockId_y * Para_s->cuBlock[1];
      k = blockId_z * Para_s->cuBlock[2];

    //ATHREAD serial operation replace GPU parallel mutli-thread in single thread-block
    for(i = 0; i < threadSum; i++) {

      if(k<Para_s->d_km&&i<Para_s->d_imt-1&&i> 0&&j<Para_s->d_jmt)	
      {	
      	tt1=Para_s->d_ppb[k *Para_s->d_imt *Para_s->d_jmt+j *Para_s->d_imt+i];	
      	tt2=tt1 *tt1;	
      	tt3=tt2 *tt1;	
      	ss1=Para_s->d_ppc[k *Para_s->d_imt *Para_s->d_jmt+j *Para_s->d_imt+i] *1000.0e0;	
      	ss2=ss1 *ss1;	
      	tm1=-Para_s->d_ppa[k *Para_s->d_imt *Para_s->d_jmt+j *Para_s->d_imt+i] /Para_s->d_od0 /10000.0e0 *Para_s->d_vit[k *Para_s->d_imt *Para_s->d_jmt+j *Para_s->d_imt+i];	
      	tm2=tm1 *tm1;	
      	tm3=tm2 *tm1;	
      	immediate = (0.785567e-3-0.301985e-5 *tt1+0.555579e-7 *tt2-0.415613e-9 *tt3+
			ss1 *(-0.356603e-6+0.788212e-8 *tt1+0.408195e-10 *tm1-0.602281e-15 *tm2)+ss2 *(0.515032e-8)+
			tm1 *(-0.121555e-7+0.192867e-9 *tt1-0.2131127e-11 *tt2)+tm2 *(0.176621e-12-0.175379e-14 *tt1)+
			tm3 *(0.12155e-17)) *Para_s->d_vit[k *Para_s->d_imt *Para_s->d_jmt+j *Para_s->d_imt+i];

	d_beta_s[k *Para_s->d_imt *Para_s->d_jmt+j *Para_s->d_imt+i] = immediate;
      	//d_alpha_s[k *Para_s->d_imt *Para_s->d_jmt+j *Para_s->d_imt+i]=(0.665157e-1+0.170907e-1 *tt1-0.203814e-3 *tt2+0.298357e-5 *tt3-0.255019e-7 *tt3 *tt1+ss1 *(0.378110e-2-0.846960e-4 *tt1-0.164759e-6 *tm1-0.251520e-11 *tm2)+ss2 *(-0.678662e-5)+tm1 *(0.380374e-4-0.933746e-6 *tt1+0.791325e-8 *tt2)+0.512857e-12 *tm2 *tt2-0.302285e-13 *tm3) *d_beta_s[k *Para_s->d_imt *Para_s->d_jmt+j *Para_s->d_imt+i] *Para_s->d_vit[k *Para_s->d_imt *Para_s->d_jmt+j *Para_s->d_imt+i];	
      	d_alpha_s[k *Para_s->d_imt *Para_s->d_jmt+j *Para_s->d_imt+i]=(0.665157e-1+0.170907e-1 *tt1-0.203814e-3 *tt2+
			0.298357e-5 *tt3-0.255019e-7 *(tt2*tt2)+ss1 *(0.378110e-2-0.846960e-4 *tt1-0.164759e-6 *tm1-0.251520e-11 *tm2)+
			ss2 *(-0.678662e-5)+tm1 *(0.380374e-4-0.933746e-6 *tt1+0.791325e-8 *tt2)+
            0.512857e-12 *tm2 *tt2-0.302285e-13 *tm3) *
			immediate * Para_s->d_vit[k *Para_s->d_imt *Para_s->d_jmt+j *Para_s->d_imt+i];	
      }	

    }

  }

    ldm_free(Para_s, sizeof(thermal_para));
}
