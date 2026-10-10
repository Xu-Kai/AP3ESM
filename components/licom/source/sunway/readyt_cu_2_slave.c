
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

//__userDef__ paraVarAttr(double,d_to,sizeof(double)*km,readonly,no,yes,no);	
//__userDef__ paraVarAttr(double,d_so,sizeof(double)*km,readonly,no,yes,no);	
//__userDef__ paraVarAttr(double,d_c,sizeof(double)*km*9,readonly,no,yes,no);	
//__userDef__ paraVarAttr(double,d_po,sizeof(double)*km,readonly,no,yes,no);	
//__userDef__ paraVarAttr(double,d_gg,sizeof(double)*km*jmt*imt,outonly,yes,no,no);	
//__userDef__ paraVarAttr(double,d_pdensity,sizeof(double)*(km+1)*jmt*imt,outonly,yes,no,no);	
//__userDef__ paraVarAttr(double,d_vit,sizeof(double)*jmt*imt*km,readonly,no,no,no);	
//__userDef__ paraVarAttr(double,d_ppb,sizeof(double)*km*jmt*imt,outonly,yes,no,no);	
//__userDef__ paraVarAttr(double,d_ppc,sizeof(double)*km*jmt*imt,outonly,yes,no,no);	
//__userDef__ paraVarAttr(double,d_at,sizeof(double)*ntra*km*jmt*imt,readonly,no,no,no);	
//__userDef__ paraVarAttr(double,d_dzp,sizeof(double)*km,readonly,no,yes,no);	
//__userDef__ paraVarAttr(double,d_od0,sizeof(double),readonly,no,no,no);	
//__userDef__ paraVarAttr(double,d_g,sizeof(double),readonly,no,no,no);	
//__userDef__ paraVarAttr(int,d_km,sizeof(int),readonly,no,no,no);	
//__userDef__ paraVarAttr(int,d_jmt,sizeof(int),readonly,no,no,no);	
//__userDef__ paraVarAttr(int,d_imt,sizeof(int),readonly,no,no,no);	
typedef struct readyt_2_kernelPara{
	int cuBlock[3];
	int cuGrid[3];
	int d_km;
	int d_jmt;
	int d_imt;
	int rank;
	double d_od0;
	double d_g;
	double *d_to;
	double *d_so;
	double *d_c;
	double *d_po;
	double *d_gg;
	double *d_pdensity;
	double *d_vit;
	double *d_ppb;
	double *d_ppc;
	double *d_at;
	double *d_dzp;
}readyt_2_para;


//dynamic atomic definition
#define CP_CORE      64
#define RCV_CORE     0
#define DATAPAKAGE   1

//ATHREAD global definition
__thread_local int blockIndex;
//__thread_local Para* Paras;

 void  readyt_cu_2_(readyt_2_para *Paras){

  //ATHREAD local variable definition
  int swGrid[3];
  int swBlock[3];
  int stride;
  int blockId_x,blockId_y,blockId_z;
  int threadId_x,threadId_y,threadId_z;


  //ATHREAD basic operation
     int i,j,k;
  blockIndex = CRTS_tid;
//  readyt_2_para *Paras = (readyt_2_para*)ldm_malloc(sizeof(readyt_2_para));
//  CRTS_dma_iget(Paras, Paras, sizeof(readyt_2_para), &dma_rply);
//  D_COUNT++;
//  CRTS_dma_wait_value(&dma_rply, D_COUNT);


 	double* d_to_s = Paras->d_to;
	double*	d_so_s = Paras->d_so;
	double*	d_c_s = Paras->d_c;
	double*	d_po_s = Paras->d_po;
	double*	d_dzp_s = Paras->d_dzp;

  //cuda grid and block affine translation
  int blockSum   = Paras->cuGrid[0] * Paras->cuGrid[1] * Paras->cuGrid[2];
  int threadSum  = Paras->cuBlock[0] * Paras->cuBlock[1] * Paras->cuBlock[2];
  int quotient   = blockSum / 64; 
  int remainder  = blockSum % 64; 
  int blockStart = quotient * blockIndex;
  if(remainder > 0)
    blockStart += (remainder > blockIndex ? blockIndex : remainder);
  int blockEnd = blockStart + quotient;
  if(remainder > 0)
    blockEnd   += (remainder > blockIndex ? 1 : 0);

  double *d_gg_s = Paras->d_gg;
  double *d_pdensity_s = Paras->d_pdensity;
  double *d_ppb_s = Paras->d_ppb;
  double *d_ppc_s = Paras->d_ppc;

  for(int block = blockStart; block < blockEnd; block++){
    blockId_z = block / (Paras->cuGrid[0] * Paras->cuGrid[1]);
    blockId_y = (block - blockId_z * (Paras->cuGrid[0] * Paras->cuGrid[1])) / Paras->cuGrid[0];
    blockId_x = block % Paras->cuGrid[0];

      j = blockId_y * Paras->cuBlock[1];
      k = blockId_z * Paras->cuBlock[2];

    //ATHREAD serial operation replace GPU parallel mutli-thread in single thread-block
    for(i = 0; i < threadSum; i++) {

      if(k<km&&j<jmt&&i<imt)	
      {	
      if(Paras->d_vit[k *jmt *imt+j *imt+i]> 0.0)	
      {	
      double tq,sq;	
	
      	tq=Paras->d_at[k *jmt *imt+j *imt+i]-d_to_s[k];	
      	sq=Paras->d_at[1 *km *jmt *imt+k *jmt *imt+j *imt+i]-d_so_s[k];	
      	d_pdensity_s[k *jmt *imt+j *imt+i]=1.0e+3+d_po_s[k]+(d_c_s[0 *km+k]+(d_c_s[3 *km+k]+d_c_s[6 *km+k] *sq) *sq+(d_c_s[2 *km+k]+d_c_s[7 *km+k] *sq+d_c_s[5 *km+k] *tq) *tq) *tq+(d_c_s[1 *km+k]+(d_c_s[4 *km+k]+d_c_s[8 *km+k] *sq) *sq) *sq;	
      }	
      else	{
	      	d_pdensity_s[k *jmt *imt+j *imt+i]=0.0;	
      }	
      	d_gg_s[k *jmt *imt+j *imt+i]=-Paras->d_od0 *Paras->d_g *d_pdensity_s[k *jmt *imt+j *imt+i] *Paras->d_vit[k *jmt *imt+j *imt+i];	
      if(k>=1)	
      {	
      	d_ppb_s[k *jmt *imt+j *imt+i]=Paras->d_vit[k *jmt *imt+j *imt+i] *(Paras->d_at[(k-1) *jmt *imt+j *imt+i]-(Paras->d_at[(k-1) *jmt *imt+j *imt+i]-Paras->d_at[k *jmt *imt+j *imt+i]) *d_dzp_s[k-1] /(d_dzp_s[k-1]+d_dzp_s[k]));	
      	d_ppc_s[k *jmt *imt+j *imt+i]=Paras->d_vit[k *jmt *imt+j *imt+i] *(Paras->d_at[1 *km *jmt *imt+(k-1) *jmt *imt+j *imt+i]-(Paras->d_at[1 *km *jmt *imt+(k-1) *jmt *imt+j *imt+i]-Paras->d_at[1 *km *jmt *imt+k *jmt *imt+j *imt+i]) *d_dzp_s[k-1] /(d_dzp_s[k-1]+d_dzp_s[k]));	
      }	
      }	

    }

  }


//    ldm_free(Paras, sizeof(readyt_2_para));
}
