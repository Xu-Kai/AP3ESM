
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
#include "constant_mod.h"
#include "param_mod.h"	

//ATHREAD DMA transfer definition 
__thread_local unsigned int D_COUNT = 0;
__thread_local crts_rply_t dma_rply = 0;
__thread_local crts_rply_t l_rply = 0;
__thread_local crts_rply_t r_rply = 0;

//__userDef__ paraVarAttr(double,d_work,sizeof(double)*jmt*imt,outonly,yes,no,no);	
//__userDef__ paraVarAttr(double,d_wka,sizeof(double)*jmt*imt*km,outonly,yes,no,no);	
//__userDef__ paraVarAttr(double,d_dzph,sizeof(double)*jmt*imt,readonly,no,no,no);	
//__userDef__ paraVarAttr(double,d_h0,sizeof(double)*jmt*imt,readonly,no,no,no);	
//__userDef__ paraVarAttr(double,d_ub,sizeof(double)*jmt*imt,readonly,no,no,no);	
//__userDef__ paraVarAttr(double,d_vb,sizeof(double)*jmt*imt,readonly,no,no,no);	
//__userDef__ paraVarAttr(double,d_au0,sizeof(double)*ny_block*nx_block,readonly,no,no,no);	
//__userDef__ paraVarAttr(double,d_ausw,sizeof(double)*ny_block*nx_block,readonly,no,no,no);	
//__userDef__ paraVarAttr(double,d_aus,sizeof(double)*ny_block*nx_block,readonly,no,no,no);	
//__userDef__ paraVarAttr(double,d_auw,sizeof(double)*ny_block*nx_block,readonly,no,no,no);	
//__userDef__ paraVarAttr(int,d_jmt,sizeof(int),readonly,no,no,no);	
//__userDef__ paraVarAttr(int,d_imt,sizeof(int),readonly,no,no,no);	
typedef struct barotri_cu_kernelPara{
	int cuBlock[3];	
	int cuGrid[3];
	int d_km;
	int d_imt;
	int d_jmt;
	int rank;
	double *d_work;
	double *d_wka;
	double *d_dzph;
	double *d_h0;
	double *d_ub;
	double *d_vb;
	double *d_au0;
	double *d_ausw;
	double *d_aus;
	double *d_auw;
}barotri_cu_Para;


//dynamic atomic definition
#define CP_CORE      64
#define RCV_CORE     0
#define DATAPAKAGE   1

//ATHREAD global definition
__thread_local int blockIndex;
//__thread_local Para* Para_s;

 void  barotr1_cu_(barotri_cu_Para *Paras){

  //ATHREAD local variable definition
  int swGrid[3];
  int swBlock[3];
  int stride;
  double work_s = 0.0e0;
  double d_wka[imt];
  double d_wka1[imt];
  double d_ub[imt];
  double d_vb[imt];
  double d_dzph[imt];

  //ATHREAD basic operation
     int i,j,k;
     int blockId_x,blockId_y,blockId_z;
  blockIndex = CRTS_tid;
  barotri_cu_Para *Para_s = (barotri_cu_Para*)ldm_malloc(sizeof(barotri_cu_Para));
  CRTS_dma_iget(Para_s, Paras, sizeof(barotri_cu_Para), &dma_rply);
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
             
         CRTS_dma_get(d_wka, &Para_s->d_wka[j * imt], sizeof(double) * imt);
         CRTS_dma_get(d_wka1, &Para_s->d_wka[imt * jmt + j * imt], sizeof(double) * imt);
//         CRTS_dma_get(d_au0, &Para_s->d_au0[j * imt], sizeof(double) * imt);
//         CRTS_dma_get(d_aus, &Para_s->d_aus[j * imt], sizeof(double) * imt);
//         CRTS_dma_get(d_auw, &Para_s->d_auw[j * imt], sizeof(double) * imt);
         CRTS_dma_get(d_ub, &Para_s->d_ub[j * imt], sizeof(double) * imt);
//         CRTS_dma_get(d_ausw, &Para_s->d_ausw[j * imt], sizeof(double) * imt);
         CRTS_dma_get(d_vb, &Para_s->d_vb[j * imt], sizeof(double) * imt);
         CRTS_dma_get(d_dzph, &Para_s->d_dzph[j * imt], sizeof(double) * imt);

         for (i = 0; i < imt; i++) {
             if(j==jmt-1||i==0)
                 work_s = 0.0e0;
             else{
                 work_s =p25 *Para_s->d_h0[j *imt+i]+
                         p25 *Para_s->d_h0[(j+1) *imt+i]+
                         p25 *Para_s->d_h0[j *imt+i-1]+
                         p25 *Para_s->d_h0[(j+1) *imt+i-1];
             }

             if(i<imt-1&&j> 0)
             {
                 d_wka[i]=d_ub[i] * (d_dzph[i]+work_s);
                 d_wka1[i]=d_vb[i] * (d_dzph[i]+work_s);
             }
         }
         CRTS_dma_put(&Para_s->d_wka[j * imt], d_wka, sizeof(double) * imt);
         CRTS_dma_put(&Para_s->d_wka[imt * jmt + j * imt], d_wka1, sizeof(double) * imt);
     }

    ldm_free(Para_s, sizeof(barotri_cu_Para));
}
