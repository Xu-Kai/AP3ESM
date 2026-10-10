
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

//__userDef__ paraVarAttr(double,d_gg,sizeof(double)*km*jmt*imt,readonly,no,no,no);	
//__userDef__ paraVarAttr(double,d_dzp,sizeof(double)*km,readonly,no,yes,no);	
//__userDef__ paraVarAttr(double,d_dlu,sizeof(double)*jmt*imt*km,outonly,yes,no,no);	
//__userDef__ paraVarAttr(double,d_dlv,sizeof(double)*jmt*imt*km,readonly,no,no,no);	
//__userDef__ paraVarAttr(double,d_ohbu,sizeof(double)*jmt*imt,readonly,no,no,no);	
//__userDef__ paraVarAttr(double,d_ohbt,sizeof(double)*jmt*imt,readonly,no,no,no);	
//__userDef__ paraVarAttr(double,d_zkt,sizeof(double)*km,readonly,no,yes,no);	
//__userDef__ paraVarAttr(double,d_viv,sizeof(double)*jmt*imt*km,readonly,no,no,no);	
//__userDef__ paraVarAttr(double,d_pxb,sizeof(double)*jmt*imt,outonly,yes,no,no);	
//__userDef__ paraVarAttr(double,d_pyb,sizeof(double)*jmt*imt,outonly,yes,no,no);	
//__userDef__ paraVarAttr(int,d_km,sizeof(int),readonly,no,no,no);	
//__userDef__ paraVarAttr(int,d_jmt,sizeof(int),readonly,no,no,no);	
//__userDef__ paraVarAttr(int,d_imt,sizeof(int),readonly,no,no,no);	
typedef struct readyt_3_2_kernelPara{
	int cuBlock[3];
	int cuGrid[3];
	int d_km;
	int d_jmt;
	int d_imt;
	int rank;
	double *d_gg;
	double *d_dzp;
	double *d_dlu;
	double *d_dlv;
	double *d_ohbu;
	double *d_ohbt;
	double *d_zkt;
	double *d_viv;
	double *d_pxb;
	double *d_pyb;
}readyt_3_2_Para;


//dynamic atomic definition
#define CP_CORE      64
#define RCV_CORE     0
#define DATAPAKAGE   1

//ATHREAD global definition
__thread_local int blockIndex;
//__thread_local Para* Paras;

 void  readyt_cu_3_2_(readyt_3_2_Para *Paras){

  //ATHREAD local variable definition
  int swGrid[3];
  int swBlock[3];
  int stride;
  int blockId_x,blockId_y,blockId_z;
  int threadId_x,threadId_y,threadId_z;


  //ATHREAD basic operation
  int i,j,k;
  blockIndex = CRTS_tid;
//  readyt_3_2_Para* Paras = (readyt_3_2_Para*)ldm_malloc(sizeof(readyt_3_2_Para));
//  CRTS_dma_iget(Paras, Paras, sizeof(readyt_3_2_Para), &dma_rply);
//  D_COUNT++;
//  CRTS_dma_wait_value(&dma_rply, D_COUNT);

  double* d_dzp_s = (double*)ldm_malloc(sizeof(double) *km);
  CRTS_dma_iget(d_dzp_s, Paras->d_dzp, sizeof(double) *km, &dma_rply);
  D_COUNT++;
  CRTS_dma_wait_value(&dma_rply, D_COUNT);

  double* d_zkt_s = (double*)ldm_malloc(sizeof(double) *km);
  CRTS_dma_iget(d_zkt_s, Paras->d_zkt, sizeof(double) *km, &dma_rply);
  D_COUNT++;
  CRTS_dma_wait_value(&dma_rply, D_COUNT);


  double *d_dlu_s = Paras->d_dlu;
  double *d_pxb_s = Paras->d_pxb;
  double *d_pyb_s = Paras->d_pyb;

  for(k = 0; k < 2; k++) {
      for (j = blockIndex; j < jmt; j += 64) {
          for (i = 0; i < imt; i++) {
              if (i < imt && j < jmt) {
                  if(k == 0) {
                      d_pxb_s[j * imt + i] = 0.0;
                      d_pyb_s[j * imt + i] = 0.0;
                  }
                  d_pxb_s[j * imt + i] +=
                              d_dzp_s[k] * Paras->d_ohbu[j * imt + i] * d_dlu_s[k * jmt * imt + j * imt + i] *
                              Paras->d_viv[k * jmt * imt + j * imt + i];
                  d_pyb_s[j * imt + i] +=
                              d_dzp_s[k] * Paras->d_ohbu[j * imt + i] * Paras->d_dlv[k * jmt * imt + j * imt + i] *
                              Paras->d_viv[k * jmt * imt + j * imt + i];

              }
          }
      }
  }

     for(k = 0; k < 2; k++) {
         for (j = blockIndex; j < jmt; j += 64) {
             for (i = 0; i < imt; i++) {
                 double abcd;
                 if (i < imt && j < jmt) {
                     if(k == 0){
                         d_dlu_s[j * imt + i] = 0.0e0;
                         d_dlu_s[1 * jmt * imt + j * imt + i] = 0.0e0;
                     }
                     abcd = Paras->d_gg[k * jmt * imt + j * imt + i] * Paras->d_ohbt[j * imt + i] * d_dzp_s[k];
                     d_dlu_s[j * imt + i] += abcd;
                     d_dlu_s[1 * jmt * imt + j * imt + i] += abcd * d_zkt_s[k];
                 }
             }
         }
     }


     for(k = 2; k < km; k++) {
         for (j = blockIndex; j < jmt; j += 64) {
             for (i = 0; i < imt; i++) {
                 double abcd;
                 if (i < imt && j < jmt) {
                     d_pxb_s[j * imt + i] +=
                             d_dzp_s[k] * Paras->d_ohbu[j * imt + i] * d_dlu_s[k * jmt * imt + j * imt + i] *
                             Paras->d_viv[k * jmt * imt + j * imt + i];
                     d_pyb_s[j * imt + i] +=
                             d_dzp_s[k] * Paras->d_ohbu[j * imt + i] * Paras->d_dlv[k * jmt * imt + j * imt + i] *
                             Paras->d_viv[k * jmt * imt + j * imt + i];

                     abcd = Paras->d_gg[k * jmt * imt + j * imt + i] * Paras->d_ohbt[j * imt + i] * d_dzp_s[k];
                     d_dlu_s[j * imt + i] += abcd;
                     d_dlu_s[1 * jmt * imt + j * imt + i] += abcd * d_zkt_s[k];
                 }

             }

         }
     }

    ldm_free(d_dzp_s, sizeof(double) *km);

    ldm_free(d_zkt_s, sizeof(double) *km);

//    ldm_free(Paras, sizeof(readyt_3_2_Para));
}
