
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

//__userDef__ paraVarAttr(double,d_psa,sizeof(double)*jmt*imt,readonly,no,no,no);	
//__userDef__ paraVarAttr(double,d_ppb,sizeof(double)*km*jmt*imt,outonly,yes,no,no);	
//__userDef__ paraVarAttr(double,d_ppc,sizeof(double)*km*jmt*imt,outonly,yes,no,no);	
//__userDef__ paraVarAttr(double,d_at,sizeof(double)*ntra*km*jmt*imt,readonly,no,no,no);	
//__userDef__ paraVarAttr(double,d_gg,sizeof(double)*km*jmt*imt,readonly,no,no,no);	
//__userDef__ paraVarAttr(double,d_vit,sizeof(double)*jmt*imt*km,readonly,no,no,no);	
//__userDef__ paraVarAttr(double,d_pp,sizeof(double)*km*jmt*imt,outonly,yes,no,no);	
//__userDef__ paraVarAttr(double,d_ppa,sizeof(double)*km*jmt*imt,outonly,yes,no,no);	
//__userDef__ paraVarAttr(double,d_dzp,sizeof(double)*km,readonly,no,yes,no);	
//__userDef__ paraVarAttr(int,d_km,sizeof(int),readonly,no,no,no);	
//__userDef__ paraVarAttr(int,d_jmt,sizeof(int),readonly,no,no,no);	
//__userDef__ paraVarAttr(int,d_imt,sizeof(int),readonly,no,no,no);	
typedef struct readyt_2_2_kernelPara{
	int cuBlock[3];
	int cuGrid[3];
	int d_km;
	int d_jmt;
	int d_imt;
	int rank;
	double  d_od0;
	double *d_psa;
	double *d_ppb;
	double *d_ppc;
	double *d_at;
	double *d_gg;
	double *d_vit;
	double *d_pp;
	double *d_ppa;
	double *d_dzp;
}readyt_2_2_para;


//dynamic atomic definition
#define CP_CORE      64
#define RCV_CORE     0
#define DATAPAKAGE   1

//ATHREAD global definition
__thread_local int blockIndex;
//__thread_local Para* Paras;

 void  readyt_cu_2_2_(readyt_2_2_para *Paras){

  //ATHREAD local variable definition
  int swGrid[3];
  int swBlock[3];
  int stride;
  int blockId_x,blockId_y,blockId_z;
  int threadId_x,threadId_y,threadId_z;


  //ATHREAD basic operation
  int i,j,k;
  blockIndex = CRTS_tid;
//  readyt_2_2_para *Paras = (readyt_2_2_para*)ldm_malloc(sizeof(readyt_2_2_para));
//  CRTS_dma_iget(Paras, Paras, sizeof(readyt_2_2_para), &dma_rply);
//  D_COUNT++;
//  CRTS_dma_wait_value(&dma_rply, D_COUNT);

  double* d_dzp_s = (double*)ldm_malloc(sizeof(double) *km);
  CRTS_dma_iget(d_dzp_s, Paras->d_dzp, sizeof(double) *km, &dma_rply);
  D_COUNT++;
  CRTS_dma_wait_value(&dma_rply, D_COUNT);


  double *d_ppb_s = Paras->d_ppb;
  double *d_ppc_s = Paras->d_ppc;
  double *d_pp_s = Paras->d_pp;
  double *d_ppa_s = Paras->d_ppa;

  for(int k = 1; k < km; k++){
      for (j = blockIndex; j < jmt; j += 64) {
        for(i = 0; i < imt; i++) {
          if(j<jmt&&i<imt)
          {
            if(k == 1) {
                d_pp_s[j * imt + i] = 0.5e0 * d_dzp_s[0] * Paras->d_gg[j * imt + i] * Paras->d_vit[j * imt + i];
                d_ppa_s[j * imt + i] = -Paras->d_od0 * Paras->d_psa[j * imt + i] * Paras->d_vit[j * imt + i];
                d_ppb_s[j * imt + i] = Paras->d_at[j * imt + i] * Paras->d_vit[j * imt + i];
                d_ppc_s[j * imt + i] = Paras->d_at[1 * km * jmt * imt + j * imt + i] * Paras->d_vit[j * imt + i];
            }
      	    d_pp_s[k *jmt *imt+j *imt+i]=Paras->d_vit[k *jmt *imt+j *imt+i] *
                    (d_pp_s[(k-1) *jmt *imt+j *imt+i]+0.5e0 *(Paras->d_gg[k *jmt *imt+j *imt+i] *
                    d_dzp_s[k]+Paras->d_gg[(k-1) *jmt *imt+j *imt+i] *d_dzp_s[k-1]));
      	    d_ppa_s[k *jmt *imt+j *imt+i]=Paras->d_vit[k *jmt *imt+j *imt+i] *
                    (d_ppa_s[(k-1) *jmt *imt+j *imt+i]+Paras->d_gg[(k-1) *jmt *imt+j *imt+i] *d_dzp_s[k-1]);
          }
        }

      }
  }

    ldm_free(d_dzp_s, sizeof(double) *km);

//    ldm_free(Paras, sizeof(readyt_2_2_para));
}
