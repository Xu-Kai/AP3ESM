
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
//__userDef__ paraVarAttr(double,d_dlv,sizeof(double)*jmt*imt*km,readonly,no,no,no);	
//__userDef__ paraVarAttr(double,d_to,sizeof(double)*km,readonly,no,yes,no);	
//__userDef__ paraVarAttr(double,d_so,sizeof(double)*km,readonly,no,yes,no);	
//__userDef__ paraVarAttr(double,d_ric,sizeof(double)*jmt*imt*kmm1,outonly,yes,no,no);	
//__userDef__ paraVarAttr(double,d_viv,sizeof(double)*jmt*imt*km,readonly,no,no,no);	
//__userDef__ paraVarAttr(double,d_odzt,sizeof(double)*km,readonly,no,yes,no);	
//__userDef__ paraVarAttr(double,d_c,sizeof(double)*km*9,readonly,no,yes,no);	
//__userDef__ paraVarAttr(double,d_rict,sizeof(double)*jmt*imt*kmm1,outonly,yes,no,no);	
//__userDef__ paraVarAttr(double,d_vit,sizeof(double)*jmt*imt*km,readonly,no,no,no);	
//__userDef__ paraVarAttr(double,d_atb,sizeof(double)*ntra*(km+1)*jmt*imt,readonly,no,no,no);	
//__userDef__ paraVarAttr(double,d_od0,sizeof(double),readonly,no,no,no);	
//__userDef__ paraVarAttr(double,d_g,sizeof(double),readonly,no,no,no);	
//__userDef__ paraVarAttr(int,d_km,sizeof(int),readonly,no,no,no);	
//__userDef__ paraVarAttr(int,d_kmm1,sizeof(int),readonly,no,no,no);	
//__userDef__ paraVarAttr(int,d_jmt,sizeof(int),readonly,no,no,no);	
//__userDef__ paraVarAttr(int,d_imt,sizeof(int),readonly,no,no,no);	
 double  dens_c_cu(double *tq, double *sq, int *kk, double *d_c){	
      double outcome;	
	
      int k_tmpt=*kk;	
	
      	outcome=(d_c[k_tmpt]+(d_c[3 *km+k_tmpt]+d_c[6 *km+k_tmpt] *(*sq)) *(*sq)+(d_c[2 *km+k_tmpt]+d_c[7 *km+k_tmpt] *(*sq)+d_c[5 *km+k_tmpt] *(*tq)) *(*tq)) *(*tq)+(d_c[1 *km+k_tmpt]+(d_c[4 *km+k_tmpt]+d_c[8 *km+k_tmpt] *(*sq)) *(*sq)) *(*sq);	
      return outcome;	

}
typedef struct dens_kernelPara{
	int cuBlock[3];
    int cuGrid[3];
	int d_km;
	int d_kmm1;
	int d_jmt;
	int d_imt;
	int rank;
	int pad;
	double d_g;
	double d_od0;
	double *d_dlu;
	double *d_dlv;
	double *d_to;
	double *d_so;
	double *d_ric;
	double *d_viv;
	double *d_odzt;
	double *d_c;
	double *d_rict;
	double *d_vit;
	double *d_atb;
}dens_para;


//dynamic atomic definition
#define CP_CORE      64
#define RCV_CORE     0
#define DATAPAKAGE   1

//ATHREAD global definition 
__thread_local int warp_size = 32; 
__thread_local int blockIndex;
//__thread_local Para* Para_s;

 void  dens_cu_(dens_para *Paras){

  //ATHREAD local variable definition
  int swGrid[3];
  int swBlock[3];
  int stride;
  int blockId_x,blockId_y,blockId_z;
  int threadId_x,threadId_y,threadId_z;

  double tup,sup,tlo,slo;
  double rholo,rhoup;
  int kp1;

     //ATHREAD basic operation
     int i,j,k;
     blockIndex = CRTS_tid;
  dens_para* Para_s = (dens_para*)ldm_malloc(sizeof(dens_para));
  CRTS_dma_iget(Para_s, Paras, sizeof(dens_para), &dma_rply);
  D_COUNT++;
  CRTS_dma_wait_value(&dma_rply, D_COUNT);

  double* d_to_s = (double*)ldm_malloc(sizeof(double) *km);
  CRTS_dma_iget(d_to_s, Para_s->d_to, sizeof(double) *km, &dma_rply);
  D_COUNT++;
  CRTS_dma_wait_value(&dma_rply, D_COUNT);

  double* d_so_s = (double*)ldm_malloc(sizeof(double) *km);
  CRTS_dma_iget(d_so_s, Para_s->d_so, sizeof(double) *km, &dma_rply);
  D_COUNT++;
  CRTS_dma_wait_value(&dma_rply, D_COUNT);

  double* d_odzt_s = (double*)ldm_malloc(sizeof(double) *km);
  CRTS_dma_iget(d_odzt_s, Para_s->d_odzt, sizeof(double) *km, &dma_rply);
  D_COUNT++;
  CRTS_dma_wait_value(&dma_rply, D_COUNT);

  double* d_c_s = (double*)ldm_malloc(sizeof(double) *km *9);
  CRTS_dma_iget(d_c_s, Para_s->d_c, sizeof(double) *km *9, &dma_rply);
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

  double *d_ric_s = Para_s->d_ric;
  double *d_rict_s = Para_s->d_rict;

  for(int block = blockStart; block < blockEnd; block++){
    blockId_z = block / (Para_s->cuGrid[0] * Para_s->cuGrid[1]);
    blockId_y = (block - blockId_z * (Para_s->cuGrid[0] * Para_s->cuGrid[1])) / Para_s->cuGrid[0];
    blockId_x = block % Para_s->cuGrid[0];
      j = blockId_y * Para_s->cuBlock[1];
      k = blockId_z * Para_s->cuBlock[2];

    //ATHREAD serial operation replace GPU parallel mutli-thread in single thread-block
    for(i = 0; i < threadSum; i++) {

      if(k<Para_s->d_kmm1&&j<Para_s->d_jmt-1&&1<=i&&i<Para_s->d_imt)	
      {
	
      	tup=Para_s->d_dlu[k *Para_s->d_jmt *Para_s->d_imt+j *Para_s->d_imt+i]-d_to_s[k+1];	
      	sup=Para_s->d_dlv[k *Para_s->d_jmt *Para_s->d_imt+j *Para_s->d_imt+i]-d_so_s[k+1];	
      	tlo=Para_s->d_dlu[(k+1) *Para_s->d_jmt *Para_s->d_imt+j *Para_s->d_imt+i]-d_to_s[k+1];	
      	slo=Para_s->d_dlv[(k+1) *Para_s->d_jmt *Para_s->d_imt+j *Para_s->d_imt+i]-d_so_s[k+1];	
      	kp1=k+1;	
      	rhoup=dens_c_cu(&tup,&sup,&kp1,d_c_s);	
      	rholo=dens_c_cu(&tlo,&slo,&kp1,d_c_s);	
      	d_ric_s[k *Para_s->d_jmt *Para_s->d_imt+j *Para_s->d_imt+i]=Para_s->d_viv[(k+1) *Para_s->d_jmt *Para_s->d_imt+j *Para_s->d_imt+i] *Para_s->d_od0 *Para_s->d_g *(rholo-rhoup) *d_odzt_s[k+1];
        //  if(Para_s->rank == 0) {
        //      double temp = Para_s->d_viv[(k+1) *Para_s->d_jmt *Para_s->d_imt+j *Para_s->d_imt+i] *Para_s->d_od0 *Para_s->d_g *(rholo-rhoup) *d_odzt_s[k+1];
        //      if(temp > 1e-8)
        //      printf("[slave %d]: cuGrid:%d,%d,%d, i,j,k:%d,%d,%d,od0:%lf,d_g:%lf,rhoup:%.16e,rholo:%.16e,d_ric_s:%.16e\n",
        //             blockIndex, Para_s->cuGrid[0], Para_s->cuGrid[1], Para_s->cuGrid[2],
        //             i,j,k,Para_s->d_od0,Para_s->d_g,rhoup,rholo,temp);
        //  }

      	tup=Para_s->d_atb[(k+1) *Para_s->d_jmt *Para_s->d_imt+j *Para_s->d_imt+i]-d_to_s[k+1];
      	sup=Para_s->d_atb[1 *(Para_s->d_km+1) *Para_s->d_jmt *Para_s->d_imt+(k+1) *Para_s->d_jmt *Para_s->d_imt+j *Para_s->d_imt+i]-d_so_s[k+1];	
      	tlo=Para_s->d_atb[(k+2) *Para_s->d_jmt *Para_s->d_imt+j *Para_s->d_imt+i]-d_to_s[k+1];	
      	slo=Para_s->d_atb[1 *(Para_s->d_km+1) *Para_s->d_jmt *Para_s->d_imt
                            +(k+2) *Para_s->d_jmt *Para_s->d_imt+j *Para_s->d_imt+i]-d_so_s[k+1];
      	kp1=k+1;	
      	rhoup=dens_c_cu(&tup,&sup,&kp1,d_c_s);	
      	rholo=dens_c_cu(&tlo,&slo,&kp1,d_c_s);	
      	d_rict_s[k *Para_s->d_jmt *Para_s->d_imt+j *Para_s->d_imt+i]=Para_s->d_vit[(k+1) *Para_s->d_jmt *Para_s->d_imt+j *Para_s->d_imt+i] *Para_s->d_od0 *Para_s->d_g *(rholo-rhoup) *d_odzt_s[k+1];
         // if(Para_s->rank == 0) {
         //     double temp = Para_s->d_vit[(k+1) *Para_s->d_jmt *Para_s->d_imt+j *Para_s->d_imt+i] *Para_s->d_od0 *Para_s->d_g *(rholo-rhoup) *d_odzt_s[k+1];
         //     if(temp > 1e-8)
         //     printf("[slave %d]: cuGrid:%d,%d,%d,i,j,k:%d,%d,%d,od0:%lf,d_g:%lf,rhoup:%.16e,rholo:%.16e,d_rict_s:%.16e\n",
         //            blockIndex, Para_s->cuGrid[0], Para_s->cuGrid[1], Para_s->cuGrid[2],
         //            i,j,k,Para_s->d_od0,Para_s->d_g,rhoup,rholo,temp);
         // }
      }	

    }

  }
    //todo: user need modify data size according to the practical situation

    ldm_free(d_to_s, sizeof(double) *km);
    //todo: user need modify data size according to the practical situation

    ldm_free(d_so_s, sizeof(double) *km);
    //todo: user need modify data size according to the practical situation

    ldm_free(d_odzt_s, sizeof(double) *km);
    //todo: user need modify data size according to the practical situation

    ldm_free(d_c_s, sizeof(double) *km *9);

    ldm_free(Para_s, sizeof(dens_para));
}
