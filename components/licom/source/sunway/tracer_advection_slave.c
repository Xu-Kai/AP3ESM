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
#include "constant_mod.h"

//#define ADV_DEBUG
typedef struct kernelPara{
	int cuGrid[3];
	int cuBlock[3];
        int Str;
	int rank;
	double *d_wkb;
	double *d_wkd;
	double *d_u_wface;
	double *d_v_sface;
	double *d_hts;
	double *d_htw;
	double *d_dxu;
	double *d_dyu;
  double *d_adv_tt;
  double *d_tf;

}Para;

typedef struct advection1Para{
	int cuGrid[3];
	int cuBlock[3];
	int d_mtracer;
	int Str;
	int d_nss;
	int rank;	
	double *d_u_wface;
	double *d_v_sface;
	double *d_ws;
	double *d_at;
	double *d_adv_tt;
	double *d_tarea_r;
	double *d_odzp;
	double *d_ax;
	double *d_ay;
	double *d_az;

}advection1_Para;

typedef struct advection2Para{
	int cuGrid[3];
	int cuBlock[3];
        int d_ntra;
        int d_mtracer;
        int Str;
	int rank;
	double *d_wkb;
	double *d_wkd;
	double *d_u_wface;
	double *d_v_sface;
	double *d_hts;
	double *d_htw;
	double *d_dxu;
	double *d_dyu;
	double *d_ws;
	double *d_at;
	double *d_adv_tt;
	double *d_tarea_r;
	double *d_odzp;
	double *d_ax;
	double *d_ay;
	double *d_az;
	double *d_hun;
	double *d_hue;
	double *d_at00;
	double *d_atmax;
	double *d_atmin;
	double *d_odzt;
	double *d_vit;
  double *d_dxt;
  double *d_dyt;
  double *d_dxtr;
  double *d_dytr;
  double *d_adv_x;
  double *d_adv_y;
  double *d_adv_z;
  double *d_zkp;
  double *d_flxqRtmp;
  double *d_flxqLtmp;
	double d_dts;

}advection2_Para;

//dynamic atomic definition
#define CP_CORE      64
#define RCV_CORE     0
#define DATAPAKAGE   1

//ATHREAD global definition 
__thread_local int warp_size = 32; 
__thread_local int blockIndex;
__thread_local Para* Para_s;

__inline__ double min1(double a, double b, double c, double d, double e, double f) { 
    double i = a; 
    if (i > b) i = b; 
    if (i > c) i = c; 
    if (i > d) i = d; 
    if (i > e) i = e; 
    if (i > f) i = f; 

    return i;
}

__inline__ double min2(double a, double b) { 
    double i = a; 
    if (i > b) i = b; 

    return i;
}

__inline__ double max1(double a, double b, double c, double d, double e, double f) { 
    double i = a; 
    if (i < b) i = b; 
    if (i < c) i = c; 
    if (i < d) i = d; 
    if (i < e) i = e; 
    if (i < f) i = f; 

    return i;
}

__inline__ double max2(double a, double b) { 
    double i = a; 
    if (i < b) i = b; 

    return i;
}

__inline__ double findMin(double a, double b, double c) {
    return fmin(fmin(a, b), c);
}

__inline__ double findMax(double a, double b, double c) {
    return fmax(fmax(a, b), c);
}

__inline__  double  sign(double x, double y){
    return (y) >= 0.0e0 ? fabs(x) : -fabs(x);

}
 void  advection_tracer_1_0_(Para *Paras){

     int blockId_x,blockId_y,blockId_z;
#ifdef ADV_DEBUG
     unsigned long count,counte;
  unsigned long misscount,misscounte;

  penv_slave1_dcache_access_init();
  penv_slave2_dcache_miss_init();

  penv_slave1_dcache_access_count(&count);
  penv_slave2_dcache_miss_count(&misscount);
#endif
  //ATHREAD basic operation
  blockIndex = CRTS_tid;
  Para_s = (Para*)ldm_malloc(sizeof(Para));
  CRTS_dma_get(Para_s, Paras, sizeof(Para));

  int Str = Para_s->Str;

  //int ret = CRTS_memcpy_sldm(array_2d,Para_s->d_dyu,sizeof(double) * imt * jmt, MEM_TO_LDM);

  //double d_wkd[imt],d_wkdjm1[imt];
  double d_v_sface[imt],d_wkb[imt],d_u_wface[imt];
  double dd_wkd[imt*2];

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

     int i,j,k;

     for(int block = blockStart; block < blockEnd; block++) {
         blockId_z = block / (Para_s->cuGrid[0] * Para_s->cuGrid[1]);
         blockId_y = (block - blockId_z * (Para_s->cuGrid[0] * Para_s->cuGrid[1])) / Para_s->cuGrid[0];
         blockId_x = block % Para_s->cuGrid[0];
         j = blockId_y * Para_s->cuBlock[1];
         k = blockId_z * Para_s->cuBlock[2];

         CRTS_dma_get(d_wkb, &Para_s->d_wkb[k * imt * jmt + j * imt], sizeof(double) * Para_s->cuBlock[0]);
         CRTS_dma_get(dd_wkd, &Para_s->d_wkd[k *jmt *imt+(j-1) *imt], sizeof(double) * Para_s->cuBlock[0] * 2);

         //ATHREAD serial operation replace GPU parallel mutli-thread in single thread-block
         for (i = 0; i < threadSum; i++) {
             Para_s->d_adv_tt[k * imt * jmt + j * imt + i] = 0.0;
             Para_s->d_tf[k * imt * jmt + j * imt + i] = 0.0;
             if(1<=i&&i<(imt-1)&& k< km)
             {
                 d_v_sface[i]=(d_wkb[i]+d_wkb[i+1]) *Para_s->d_hts[j *imt+i] *0.25;
             }
             if(1<=j&&j<jmt-1&& k< km)
             {
                 //d_u_wface[i]=(d_wkdjm1[i]+d_wkd[i]) *Para_s->d_htw[j *imt+i] *0.25;
                 d_u_wface[i]=(dd_wkd[i]+dd_wkd[imt+i]) *Para_s->d_htw[j *imt+i] *0.25;
             }

         }
         CRTS_dma_put(&Para_s->d_v_sface[k * imt * jmt + j * imt], d_v_sface, sizeof(double) * Para_s->cuBlock[0]);
         CRTS_dma_put(&Para_s->d_u_wface[k * imt * jmt + j * imt], d_u_wface, sizeof(double) * Para_s->cuBlock[0]);
     }

#ifdef ADV_DEBUG
     penv_slave1_dcache_access_count(&counte);
    penv_slave2_dcache_miss_count(&misscounte);
    if(Para_s->rank == 0)
        printf("[slave %d]:advection_0 count:%d,counte:%d,misscount:%d,misscounte:%d,total count:%d, total misscount:%d, misscount rate:%f\n",
		blockIndex,count,counte,misscount,misscounte,counte-count, misscounte-misscount,
		 (float)(misscounte-misscount)/(float)(counte-count));
#endif
    ldm_free(Para_s, sizeof(Para));
}

 void  advection_tracer_1_(advection1_Para *Paras){	

  //ATHREAD local variable definition
  int swGrid[3];
  int swBlock[3];
  int stride;
  int blockId_x,blockId_y,blockId_z;
  int threadId_x,threadId_y,threadId_z;
#ifdef ADV_DEBUG
     unsigned long count,counte;
  unsigned long misscount,misscounte;

  penv_slave1_dcache_access_init();
  penv_slave2_dcache_miss_init();

  penv_slave1_dcache_access_count(&count);
  penv_slave2_dcache_miss_count(&misscount);
#endif
  //ATHREAD basic operation
  blockIndex = CRTS_tid;
  advection1_Para* Para_s = (advection1_Para*)ldm_malloc(sizeof(advection1_Para));
  CRTS_dma_get(Para_s, Paras, sizeof(advection1_Para));
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
   int Str = Para_s->Str;
   int NSS = Para_s->d_nss;
   double* d_adv_tt = Para_s->d_adv_tt;
   double* d_ax = Para_s->d_ax;
   double* d_ay = Para_s->d_ay;
   double* d_az = Para_s->d_az;
   int d_mtracer = Para_s->d_mtracer;

  for(int block = blockStart; block < blockEnd; block++){
    blockId_z = block / (Para_s->cuGrid[0] * Para_s->cuGrid[1]);
    blockId_y = (block - blockId_z * (Para_s->cuGrid[0] * Para_s->cuGrid[1])) / Para_s->cuGrid[0];
    blockId_x = block % Para_s->cuGrid[0];

    bool dataLoad = true;

    //ATHREAD serial operation replace GPU parallel mutli-thread in single thread-block
    for(int threadId = 0; threadId < threadSum; threadId++) {
      threadId_x = threadId % Para_s->cuBlock[0];
      threadId_z = threadId / (Para_s->cuBlock[0] * Para_s->cuBlock[1]);
      threadId_y = (threadId - threadId_z * (Para_s->cuBlock[0] * Para_s->cuBlock[1])) / Para_s->cuBlock[0];

      int i=blockId_x *Para_s->cuBlock[0]+threadId_x;	
	
      int j=blockId_y *Para_s->cuBlock[1]+threadId_y;	
	
      int k=blockId_z *Para_s->cuBlock[2]+threadId_z;	
	
      double adv_z1,adv_z2;	
	
      if(Str==0)	
      {	
      if(2<=j&&j<(jmt-2)&&2<=i&&i<(imt-2))	
      {	
      	d_adv_tt[j *imt+i]=(-Para_s->d_u_wface[j *imt+i] *(Para_s->d_at[j *imt+i]-Para_s->d_at[j *imt+i-1])-Para_s->d_u_wface[j *imt+i+1] *(Para_s->d_at[j *imt+i+1]-Para_s->d_at[j *imt+i])-Para_s->d_v_sface[j *imt+i] *(Para_s->d_at[(j+1) *imt+i]-Para_s->d_at[j *imt+i])-Para_s->d_v_sface[(j-1) *imt+i] *(Para_s->d_at[j *imt+i]-Para_s->d_at[(j-1) *imt+i])) *Para_s->d_tarea_r[j *nx_block+i];	
//      	d_ax[d_mtracer *km *jmt *imt+j *imt+i]=d_ax[d_mtracer *km *jmt *imt+j *imt+i]+(-Para_s->d_u_wface[j *imt+i] *(Para_s->d_at[j *imt+i]-Para_s->d_at[j *imt+i-1])-Para_s->d_u_wface[j *imt+i+1] *(Para_s->d_at[j *imt+i+1]-Para_s->d_at[j *imt+i])) *Para_s->d_tarea_r[j *nx_block+i] /(double)(NSS);
//      	d_ay[d_mtracer *km *jmt *imt+j *imt+i]=d_ay[d_mtracer *km *jmt *imt+j *imt+i]+(-Para_s->d_v_sface[j *imt+i] *(Para_s->d_at[(j+1) *imt+i]-Para_s->d_at[j *imt+i])-Para_s->d_v_sface[(j-1) *imt+i] *(Para_s->d_at[j *imt+i]-Para_s->d_at[(j-1) *imt+i])) *Para_s->d_tarea_r[j *nx_block+i] /(double)(NSS);
      	adv_z2=Para_s->d_ws[1 *jmt *imt+j *imt+i] *(Para_s->d_at[j *imt+i]-Para_s->d_at[1 *jmt *imt+j *imt+i]);	
      	d_adv_tt[j *imt+i]=d_adv_tt[j *imt+i]-p5 *Para_s->d_odzp[0] *adv_z2;	
//      	d_az[d_mtracer *km *jmt *imt+j *imt+i]=d_az[d_mtracer *km *jmt *imt+j *imt+i]-p5 *Para_s->d_odzp[0] *adv_z2 /(double)(NSS);
      }	
      if(1<=k&&k<km-1&&2<=j&&j<(jmt-2)&&2<=i&&i<(imt-2))	
      {	
      	d_adv_tt[k *jmt *imt+j *imt+i]=(-Para_s->d_u_wface[k *jmt *imt+j *imt+i] *(Para_s->d_at[k *jmt *imt+j *imt+i]-Para_s->d_at[k *jmt *imt+j *imt+i-1])-Para_s->d_u_wface[k *jmt *imt+j *imt+i+1] *(Para_s->d_at[k *jmt *imt+j *imt+i+1]-Para_s->d_at[k *jmt *imt+j *imt+i])-Para_s->d_v_sface[k *jmt *imt+j *imt+i] *(Para_s->d_at[k *jmt *imt+(j+1) *imt+i]-Para_s->d_at[k *jmt *imt+j *imt+i])-Para_s->d_v_sface[k *jmt *imt+(j-1) *imt+i] *(Para_s->d_at[k *jmt *imt+j *imt+i]-Para_s->d_at[k *jmt *imt+(j-1) *imt+i])) *Para_s->d_tarea_r[j *nx_block+i];	
//      	d_ax[d_mtracer *km *jmt *imt+k *jmt *imt+j *imt+i]=d_ax[d_mtracer *km *jmt *imt+k *jmt *imt+j *imt+i]+(-Para_s->d_u_wface[k *jmt *imt+j *imt+i] *(Para_s->d_at[k *jmt *imt+j *imt+i]-Para_s->d_at[k *jmt *imt+j *imt+i-1])-Para_s->d_u_wface[k *jmt *imt+j *imt+i+1] *(Para_s->d_at[k *jmt *imt+j *imt+i+1]-Para_s->d_at[k *jmt *imt+j *imt+i])) *Para_s->d_tarea_r[j *nx_block+i] /(double)(NSS);
//      	d_ay[d_mtracer *km *jmt *imt+k *jmt *imt+j *imt+i]=d_ay[d_mtracer *km *jmt *imt+k *jmt *imt+j *imt+i]+(-Para_s->d_v_sface[k *jmt *imt+j *imt+i] *(Para_s->d_at[k *jmt *imt+(j+1) *imt+i]-Para_s->d_at[k *jmt *imt+j *imt+i])-Para_s->d_v_sface[k *jmt *imt+(j-1) *imt+i] *(Para_s->d_at[k *jmt *imt+j *imt+i]-Para_s->d_at[k *jmt *imt+(j-1) *imt+i])) *Para_s->d_tarea_r[j *nx_block+i] /(double)(NSS);
      	adv_z1=Para_s->d_ws[k *jmt *imt+j *imt+i] *(Para_s->d_at[(k-1) *jmt *imt+j *imt+i]-Para_s->d_at[k *jmt *imt+j *imt+i]);	
      	adv_z2=Para_s->d_ws[(k+1) *jmt *imt+j *imt+i] *(Para_s->d_at[k *jmt *imt+j *imt+i]-Para_s->d_at[(k+1) *jmt *imt+j *imt+i]);	
      	d_adv_tt[k *jmt *imt+j *imt+i]=d_adv_tt[k *jmt *imt+j *imt+i]-p5 *Para_s->d_odzp[k] *(adv_z1+adv_z2);	
//      	d_az[d_mtracer *km *jmt *imt+k *jmt *imt+j *imt+i]=d_az[d_mtracer *km *jmt *imt+k *jmt *imt+j *imt+i]-p5 *Para_s->d_odzp[k] *(adv_z1+adv_z2) /(double)(NSS);
      }	
      if(2<=j&&j<(jmt-2)&&2<=i&&i<(imt-2))	
      {	
      	d_adv_tt[(km-1) *jmt *imt+j *imt+i]=(-Para_s->d_u_wface[(km-1) *jmt *imt+j *imt+i] *(Para_s->d_at[(km-1) *jmt *imt+j *imt+i]-Para_s->d_at[(km-1) *jmt *imt+j *imt+i-1])-Para_s->d_u_wface[(km-1) *jmt *imt+j *imt+i+1] *(Para_s->d_at[(km-1) *jmt *imt+j *imt+i+1]-Para_s->d_at[(km-1) *jmt *imt+j *imt+i])-Para_s->d_v_sface[(km-1) *jmt *imt+j *imt+i] *(Para_s->d_at[(km-1) *jmt *imt+(j+1) *imt+i]-Para_s->d_at[(km-1) *jmt *imt+j *imt+i])-Para_s->d_v_sface[(km-1) *jmt *imt+(j-1) *imt+i] *(Para_s->d_at[(km-1) *jmt *imt+j *imt+i]-Para_s->d_at[(km-1) *jmt *imt+(j-1) *imt+i])) *Para_s->d_tarea_r[j *nx_block+i];	
//      	d_ax[d_mtracer *km *jmt *imt+(km-1) *jmt *imt+j *imt+i]=d_ax[d_mtracer *km *jmt *imt+(km-1) *jmt *imt+j *imt+i]+(-Para_s->d_u_wface[(km-1) *jmt *imt+j *imt+i] *(Para_s->d_at[(km-1) *jmt *imt+j *imt+i]-Para_s->d_at[(km-1) *jmt *imt+j *imt+i-1])-Para_s->d_u_wface[(km-1) *jmt *imt+j *imt+i+1] *(Para_s->d_at[(km-1) *jmt *imt+j *imt+i+1]-Para_s->d_at[(km-1) *jmt *imt+j *imt+i])) *Para_s->d_tarea_r[j *nx_block+i] /(double)(NSS);
//      	d_ay[d_mtracer *km *jmt *imt+(km-1) *jmt *imt+j *imt+i]=d_ay[d_mtracer *km *jmt *imt+(km-1) *jmt *imt+j *imt+i]+(-Para_s->d_v_sface[(km-1) *jmt *imt+j *imt+i] *(Para_s->d_at[(km-1) *jmt *imt+(j+1) *imt+i]-Para_s->d_at[(km-1) *jmt *imt+j *imt+i])-Para_s->d_v_sface[(km-1) *jmt *imt+(j-1) *imt+i] *(Para_s->d_at[(km-1) *jmt *imt+j *imt+i]-Para_s->d_at[(km-1) *jmt *imt+(j-1) *imt+i])) *Para_s->d_tarea_r[j *nx_block+i] /(double)(NSS);
      	adv_z1=Para_s->d_ws[(km-1) *jmt *imt+j *imt+i] *(Para_s->d_at[(km-2) *jmt *imt+j *imt+i]-Para_s->d_at[(km-1) *jmt *imt+j *imt+i]);	
      	d_adv_tt[(km-1) *jmt *imt+j *imt+i]=d_adv_tt[(km-1) *jmt *imt+j *imt+i]-p5 *Para_s->d_odzp[km-1] *adv_z1;	
//      	d_az[d_mtracer *km *jmt *imt+(km-1) *jmt *imt+j *imt+i]=d_az[d_mtracer *km *jmt *imt+(km-1) *jmt *imt+j *imt+i]-p5 *Para_s->d_odzp[km-1] *adv_z1 /(double)(NSS);
      }	
      }	
      else	
      if(Str==1)	
      {	
      if(k<km&&2<=j&&j<(jmt-2)&&2<=i&&i<(imt-2))	
      {	
      	d_adv_tt[k *jmt *imt+j *imt+i]=(-Para_s->d_u_wface[k *jmt *imt+j *imt+i] *(Para_s->d_at[k *jmt *imt+j *imt+i]+Para_s->d_at[k *jmt *imt+j *imt+i-1])+Para_s->d_u_wface[k *jmt *imt+j *imt+i+1] *(Para_s->d_at[k *jmt *imt+j *imt+i+1]+Para_s->d_at[k *jmt *imt+j *imt+i])-Para_s->d_v_sface[k *jmt *imt+(j-1) *imt+i] *(Para_s->d_at[k *jmt *imt+j *imt+i]+Para_s->d_at[k *jmt *imt+(j-1) *imt+i])+Para_s->d_v_sface[k *jmt *imt+j *imt+i] *(Para_s->d_at[k *jmt *imt+(j+1) *imt+i]+Para_s->d_at[k *jmt *imt+j *imt+i])) *Para_s->d_tarea_r[j *nx_block+i];	
      if(k==0)	
	adv_z1=0.0e0;	
	
      else	
      	adv_z1=Para_s->d_ws[k *jmt *imt+j *imt+i] *(Para_s->d_at[(k-1) *jmt *imt+j *imt+i]+Para_s->d_at[k *jmt *imt+j *imt+i]) *p5;	
	
      if(k==km-1)	
	adv_z2=0.0e0;	
	
      else	
      	adv_z2=Para_s->d_ws[(k+1) *jmt *imt+j *imt+i] *(Para_s->d_at[k *jmt *imt+j *imt+i]+Para_s->d_at[(k+1) *jmt *imt+j *imt+i]) *p5;	
	
      	d_adv_tt[k *jmt *imt+j *imt+i]=d_adv_tt[k *jmt *imt+j *imt+i]-Para_s->d_odzp[k] *(adv_z2-adv_z1);	
      }	
      }	
	

    }

  }
#ifdef ADV_DEBUG
     penv_slave1_dcache_access_count(&counte);
    penv_slave2_dcache_miss_count(&misscounte);
    if(Para_s->rank == 0)
        printf("[slave %d]:advection_1 count:%d,counte:%d,misscount:%d,misscounte:%d,total count:%d, total misscount:%d, misscount rate:%f\n",
		blockIndex,count,counte,misscount,misscounte,counte-count, misscounte-misscount,
		 (float)(misscounte-misscount)/(float)(counte-count));
#endif
    ldm_free(Para_s, sizeof(advection1_Para));
}

 void  advection_tracer_2_(advection2_Para *Paras){

     int blockId_x,blockId_y,blockId_z;

#ifdef ADV_DEBUG
     unsigned long count,counte;
  unsigned long misscount,misscounte;

  penv_slave1_dcache_access_init();
  penv_slave2_dcache_miss_init();

  penv_slave1_dcache_access_count(&count);
  penv_slave2_dcache_miss_count(&misscount);
#endif

  //ATHREAD basic operation
  blockIndex = CRTS_tid;
  advection2_Para* Para_s = (advection2_Para*)ldm_malloc(sizeof(advection2_Para));
  CRTS_dma_get(Para_s, Paras, sizeof(advection2_Para));
  //D_COUNT++;
  //CRTS_dma_wait_value(&dma_rply, D_COUNT);

 // int ret = CRTS_memcpy_sldm(array_2d,Para_s->d_tarea_r,sizeof(double) * imt * jmt, MEM_TO_LDM);
//  ret = CRTS_memcpy_sldm(array1_2d,Para_s->d_hue,sizeof(double) * imt * jmt, MEM_TO_LDM);
 // ret = CRTS_memcpy_sldm(array2_2d,Para_s->d_hts,sizeof(double) * imt * jmt, MEM_TO_LDM);

  double d_dts = Para_s->d_dts;
  double *d_odzt = Para_s->d_odzt;

  double adv_xy1,adv_xy2,adv_xy3,adv_xy4;   
   
  double wt1,wt2;   
  
  double adv_zz,adv_za,adv_zb1,adv_zb2;    
  
  double adv_x0,adv_y0,adv_xx,adv_yy;    
  
  double adv_c1,adv_c2,adv_zc; 


//  double d_u_wface[imt],d_ws[imt],d_wskp1[imt],d_atmax[imt],d_atmin[imt],d_at00[imt],d_atkm1[imt],d_atkp1[imt];
  double d_atmax[imt],d_atmin[imt],d_at00[imt];
  double *d_u_wface, *d_ws,*d_wskp1,*d_atkm1,*d_atkp1;
  //double dd_at[imt*3];
  double *dd_at;
 // double dd_v_sface[2*imt];
  double *dd_v_sface;
  //double dd_vit[imt*3];
  double *dd_vit;

  wt1=-1.0e10;    
  wt2=+1.0e10;

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

     int i,j,k;

     for(int block = blockStart; block < blockEnd; block++) {
         blockId_z = block / (Para_s->cuGrid[0] * Para_s->cuGrid[1]);
         blockId_y = (block - blockId_z * (Para_s->cuGrid[0] * Para_s->cuGrid[1])) / Para_s->cuGrid[0];
         blockId_x = block % Para_s->cuGrid[0];
         j = blockId_y * Para_s->cuBlock[1];
         k = blockId_z * Para_s->cuBlock[2];

//         CRTS_dma_get(d_u_wface, &Para_s->d_u_wface[k * imt * jmt + j * imt], sizeof(double) * Para_s->cuBlock[0]);
//         CRTS_dma_get(d_ws, &Para_s->d_ws[k * imt * jmt + j * imt], sizeof(double) * Para_s->cuBlock[0]);
//         CRTS_dma_get(d_wskp1, &Para_s->d_ws[(k+1) * imt * jmt + j * imt], sizeof(double) * Para_s->cuBlock[0]);
//         CRTS_dma_get(d_atkm1, &Para_s->d_at[(k-1) * imt * jmt + j * imt], sizeof(double) * Para_s->cuBlock[0]);
//         CRTS_dma_get(d_atkp1, &Para_s->d_at[(k+1) * imt * jmt + j * imt], sizeof(double) * Para_s->cuBlock[0]);
//
//         CRTS_dma_get(dd_at, &Para_s->d_at[k * imt * jmt + (j-1) * imt], sizeof(double) * Para_s->cuBlock[0] * 3);
//         CRTS_dma_get(dd_v_sface, &Para_s->d_v_sface[k * imt * jmt + (j-1) * imt], sizeof(double) * Para_s->cuBlock[0] * 2);
//         CRTS_dma_get(dd_vit, &Para_s->d_vit[k * imt * jmt + (j-1) * imt], sizeof(double) * Para_s->cuBlock[0] * 3);

         d_u_wface = &Para_s->d_u_wface[k * imt * jmt + j * imt];
         d_ws = &Para_s->d_ws[k * imt * jmt + j * imt];
         d_wskp1 =  &Para_s->d_ws[(k+1) * imt * jmt + j * imt];
         d_atkm1 = &Para_s->d_at[(k-1) * imt * jmt + j * imt];
         d_atkp1 =&Para_s->d_at[(k+1) * imt * jmt + j * imt];
         dd_at = &Para_s->d_at[k * imt * jmt + (j-1) * imt];
         dd_v_sface = &Para_s->d_v_sface[k * imt * jmt + (j-1) * imt];
         dd_vit = &Para_s->d_vit[k * imt * jmt + (j-1) * imt];

         for (i = 0; i < threadSum; i++) {
             if(1<=i&&i<(imt-1)&&1<=j&&j<(jmt-1)&& k< km){

                 adv_x0=(dd_at[imt+i+1]+dd_at[imt+i]) *d_u_wface[i+1] *Para_s->d_tarea_r[j *imt+i]-
                         (dd_at[imt+i]+dd_at[imt+i-1]) *d_u_wface[i] *Para_s->d_tarea_r[j *imt+i];
                 adv_y0=((dd_at[2*imt+i]+dd_at[imt+i]) *dd_v_sface[imt+i]-(dd_at[imt+i]+dd_at[i]) *dd_v_sface[i]) *
                         Para_s->d_tarea_r[j *imt+i];
                 adv_xy1=-d_dts *(dd_at[imt+i+1]-dd_at[imt+i]) *2.0 *Para_s->d_tarea_r[j *imt+i] *
                         d_u_wface[i+1] *d_u_wface[i+1] /(Para_s->d_htw[j *imt+i+1] *
                         Para_s->d_hun[j *imt+i+1]);
                 adv_xy2=d_dts *(dd_at[imt+i]-dd_at[imt+i-1]) *2.0 *Para_s->d_tarea_r[j *imt+i] *
                         d_u_wface[i] *d_u_wface[i] /(Para_s->d_htw[j *imt+i] *
                         Para_s->d_hun[j *imt+i]);
                 adv_xy3=-d_dts *(dd_at[2*imt+i]-dd_at[imt+i]) *2.0 *Para_s->d_tarea_r[j *imt+i] *
                         dd_v_sface[imt+i] *dd_v_sface[imt+i] /(Para_s->d_hts[j *imt+i] *
                         Para_s->d_hue[j *imt+i]);
                 adv_xy4=d_dts *(dd_at[imt+i]-dd_at[i]) *2.0 *Para_s->d_tarea_r[j *imt+i] *
                         dd_v_sface[i] *dd_v_sface[i] /(Para_s->d_hts[(j-1) *imt+i] *
                         Para_s->d_hue[(j-1) *imt+i]);
                 adv_c1=-dd_at[imt+i] *(d_u_wface[i+1]-d_u_wface[i]) *Para_s->d_tarea_r[j *imt+i] *2.0;
                 adv_c2=-dd_at[imt+i] *(dd_v_sface[imt+i]-dd_v_sface[i]) *Para_s->d_tarea_r[j *imt+i] *2.0;

                 if(k>0 && k < (km-1))
                 {
                     adv_za=0.5 *Para_s->d_odzp[k] *d_ws[i] *(dd_at[imt+i]+d_atkm1[i])-
                             0.5 *Para_s->d_odzp[k] *d_wskp1[i] *(dd_at[imt+i]+d_atkp1[i]);
                     adv_zb1=-0.5 *Para_s->d_odzp[k] *d_ws[i] *d_ws[i] *d_odzt[k] *
                             (d_atkm1[i]-dd_at[imt+i]) *d_dts;
                     adv_zb2=0.5 *Para_s->d_odzp[k] *d_wskp1[i] *d_wskp1[i] *d_odzt[k+1] *
                             (dd_at[imt+i]-d_atkp1[i]) *d_dts;
                     adv_zc=-Para_s->d_odzp[k] *dd_at[imt+i] *(d_ws[i]-d_wskp1[i]);

                     d_atmax[i]=max2(max1(dd_at[imt+i] *dd_vit[imt+i]+(1.0e0-dd_vit[imt+i]) *wt1,
                             dd_at[i] *dd_vit[i]+(1.0e0-dd_vit[i]) *wt1,
                             dd_at[2*imt+i] *dd_vit[2*imt+i]+(1.0e0-dd_vit[2*imt+i]) *wt1,
                             dd_at[imt+i-1] *dd_vit[imt+i-1]+(1.0e0-dd_vit[imt+i-1]) *wt1,
                             dd_at[imt+i+1] *dd_vit[imt+i+1]+(1.0e0-dd_vit[imt+i+1]) *wt1,
                             d_atkm1[i] *Para_s->d_vit[(k-1) *jmt *imt+j *imt+i]+
                             (1.0e0-Para_s->d_vit[(k-1) *jmt *imt+j *imt+i]) *wt1),
                              d_atkp1[i] *Para_s->d_vit[(k+1) *jmt *imt+j *imt+i]+
                              (1.0e0-Para_s->d_vit[(k+1) *jmt *imt+j *imt+i]) *wt1);

                     d_atmin[i]=min2(min1(dd_at[imt+i] *dd_vit[imt+i]+(1.0e0-dd_vit[imt+i]) *wt2,
                             dd_at[i] *dd_vit[i]+(1.0e0-dd_vit[i]) *wt2,
                             dd_at[2*imt+i] *dd_vit[2*imt+i]+(1.0e0-dd_vit[2*imt+i]) *wt2,
                             dd_at[imt+i-1] *dd_vit[imt+i-1]+(1.0e0-dd_vit[imt+i-1]) *wt2,
                             dd_at[imt+i+1] *dd_vit[imt+i+1]+(1.0e0-dd_vit[imt+i+1]) *wt2,
                             d_atkm1[i] *Para_s->d_vit[(k-1) *jmt *imt+j *imt+i]+
                             (1.0e0-Para_s->d_vit[(k-1) *jmt *imt+j *imt+i]) *wt2),
                             d_atkp1[i] *Para_s->d_vit[(k+1) *jmt *imt+j *imt+i]+
                             (1.0e0-Para_s->d_vit[(k+1) *jmt *imt+j *imt+i]) *wt2);

                 }
                 else if(k==0)
                 {
                     adv_za=-0.5 *Para_s->d_odzp[0] *d_wskp1[i] *(d_atkp1[i]+dd_at[imt+i]);
                     adv_zb1=0.0;
                     adv_zb2=0.5 *Para_s->d_odzp[0] *d_wskp1[i] *d_wskp1[i] *d_odzt[1] *
                             (dd_at[imt+i]-d_atkp1[i]) *d_dts;
                     adv_zc=Para_s->d_odzp[0] *dd_at[imt+i] *d_wskp1[i];

                     d_atmax[i]=max1(dd_at[imt+i] *dd_vit[imt+i]+(1.0e0-dd_vit[imt+i]) *wt1,
                             dd_at[i] *dd_vit[i]+(1.0e0-dd_vit[i]) *wt1,
                             dd_at[2*imt+i] *dd_vit[2*imt+i]+(1.0e0-dd_vit[2*imt+i]) *wt1,
                             dd_at[imt+i-1] *dd_vit[imt+i-1]+(1.0e0-dd_vit[imt+i-1]) *wt1,
                             dd_at[imt+i+1] *dd_vit[imt+i+1]+(1.0e0-dd_vit[imt+i+1]) *wt1,
                             d_atkp1[i] *Para_s->d_vit[(k+1) *jmt *imt+j *imt+i]+
                             (1.0e0-Para_s->d_vit[(k+1) *jmt *imt+j *imt+i]) *wt1);

                     d_atmin[i]=min1(dd_at[imt+i] *dd_vit[imt+i]+(1.0e0-dd_vit[imt+i]) *wt2,
                             dd_at[i] *Para_s->d_vit[(j-1) *imt+i]+(1.0e0-dd_vit[i]) *wt2,
                             dd_at[2*imt+i] *dd_vit[2*imt+i]+(1.0e0-dd_vit[2*imt+i]) *wt2,
                             dd_at[imt+i-1] *dd_vit[imt+i-1]+(1.0e0-dd_vit[imt+i-1]) *wt2,
                             dd_at[imt+i+1] *dd_vit[imt+i+1]+(1.0e0-dd_vit[imt+i+1]) *wt2,
                             d_atkp1[i] *Para_s->d_vit[(k+1) *jmt *imt+j *imt+i]+
                             (1.0e0-Para_s->d_vit[(k+1) *jmt *imt+j *imt+i]) *wt2);

                 }
                 else
                 {
                     adv_za=0.5 *Para_s->d_odzp[km-1] *d_ws[i] *(dd_at[imt+i]+d_atkm1[i]);
                     adv_zb1=-0.5 *Para_s->d_odzp[km-1] *d_ws[i] *d_ws[i] *d_odzt[km-1] *
                             (d_atkm1[i]-dd_at[imt+i]) *d_dts;
                     adv_zb2=0.0;
                     adv_zc=-Para_s->d_odzp[km-1] *dd_at[imt+i] *d_ws[i];

                     d_atmax[i]=max1(dd_at[imt+i] *dd_vit[imt+i]+(1.0e0-dd_vit[imt+i]) *wt1,
                            dd_at[i] *dd_vit[i]+(1.0e0-dd_vit[i]) *wt1,
                            dd_at[2*imt+i] *dd_vit[2*imt+i]+(1.0e0-dd_vit[2*imt+i]) *wt1,
                            dd_at[imt+i-1] *dd_vit[imt+i-1]+(1.0e0-dd_vit[imt+i-1]) *wt1,
                            dd_at[imt+i+1] *dd_vit[imt+i+1]+(1.0e0-dd_vit[imt+i+1]) *wt1,
                            d_atkm1[i] *Para_s->d_vit[(k-1) *jmt *imt+j *imt+i]+
                            (1.0e0-Para_s->d_vit[(k-1) *jmt *imt+j *imt+i]) *wt1);

                     d_atmin[i]=min1(dd_at[imt+i] *dd_vit[imt+i]+(1.0e0-dd_vit[imt+i]) *wt2,
                            dd_at[i] *dd_vit[i]+(1.0e0-dd_vit[i]) *wt2,
                            dd_at[2*imt+i] *dd_vit[2*imt+i]+(1.0e0-dd_vit[2*imt+i]) *wt2,
                            dd_at[imt+i-1] *dd_vit[imt+i-1]+(1.0e0-dd_vit[imt+i-1]) *wt2,
                            dd_at[imt+i+1] *dd_vit[imt+i+1]+(1.0e0-dd_vit[imt+i+1]) *wt2,
                            d_atkm1[i] *Para_s->d_vit[(k-1) *jmt *imt+j *imt+i]+
                            (1.0e0-Para_s->d_vit[(k-1) *jmt *imt+j *imt+i]) *wt2);

                 }

                 adv_xx=-(adv_x0+adv_xy1+adv_xy2+adv_c1);
                 adv_yy=-(adv_y0+adv_xy3+adv_xy4+adv_c2);
                 adv_zz=-(adv_zb1+adv_zb2+adv_za+adv_zc);
                 d_at00[i]=dd_at[imt+i]+(adv_xx+adv_yy+adv_zz) *d_dts;

             }

         }

         CRTS_dma_put(&Para_s->d_atmax[k * imt * jmt + j * imt], d_atmax, sizeof(double) * Para_s->cuBlock[0]);
         CRTS_dma_put(&Para_s->d_atmin[k * imt * jmt + j * imt], d_atmin, sizeof(double) * Para_s->cuBlock[0]);
         CRTS_dma_put(&Para_s->d_at00[k * imt * jmt + j * imt], d_at00, sizeof(double) * Para_s->cuBlock[0]);

     }

#ifdef ADV_DEBUG
    penv_slave1_dcache_access_count(&counte);
    penv_slave2_dcache_miss_count(&misscounte);
    if(Para_s->rank == 0)
        printf("[slave %d]:advection_2 count:%d,counte:%d,misscount:%d,misscounte:%d,total count:%d, total misscount:%d, misscount rate:%f\n",
		blockIndex,count,counte,misscount,misscounte,counte-count, misscounte-misscount,
		 (float)(misscounte-misscount)/(float)(counte-count));
#endif
    ldm_free(Para_s, sizeof(advection2_Para));
}

 void  advection_tracer_3_(advection2_Para *Paras){	

  //ATHREAD local variable definition
  int swGrid[3];
  int swBlock[3];
  int stride;
  int blockId_x,blockId_y,blockId_z;
#ifdef ADV_DEBUG
     unsigned long count,counte;
  unsigned long misscount,misscounte;

  penv_slave1_dcache_access_init();
  penv_slave2_dcache_miss_init();

  penv_slave1_dcache_access_count(&count);
  penv_slave2_dcache_miss_count(&misscount);
#endif
  //ATHREAD basic operation
  blockIndex = CRTS_tid;
  advection2_Para* Para_s = (advection2_Para*)ldm_malloc(sizeof(advection2_Para));
  CRTS_dma_get(Para_s, Paras, sizeof(advection2_Para));
  //D_COUNT++;
  //CRTS_dma_wait_value(&dma_rply, D_COUNT);

 // int ret = CRTS_memcpy_sldm(array_2d,Para_s->d_tarea_r,sizeof(double) * imt * jmt, MEM_TO_LDM);
  //ret = CRTS_memcpy_sldm(array1_2d,Para_s->d_hue,sizeof(double) * imt * jmt, MEM_TO_LDM);
  //ret = CRTS_memcpy_sldm(array2_2d,Para_s->d_hts,sizeof(double) * imt * jmt, MEM_TO_LDM);

//  double dd_at00[3*imt];
//
//  double d_at00pk1[imt];
//  double d_at00mk1[imt];
//
//  double dd_atmax[3*imt];
//  double d_atmaxpk1[imt];
//  double d_atmaxmk1[imt];
//
//  double dd_atmin[3*imt];
//  double d_atminpk1[imt];
//  double d_atminmk1[imt];
//
//  double dd_at[3*imt];
//  double d_atpk1[imt];
//  double d_atmk1[imt];
//  double d_u_wface[imt];
//
//  double dd_v_wface[2*imt];
//  double d_ws[imt];
//  double d_wspk1[imt];
//  double d_tarea_r[imt];
//
//  double dd_hue[2*imt];
//
//  double dd_hts[2*imt];

     double *dd_at00;

     double *d_at00pk1;
     double *d_at00mk1;

     double *dd_atmax;
     double *d_atmaxpk1;
     double *d_atmaxmk1;

     double *dd_atmin;
     double *d_atminpk1;
     double *d_atminmk1;

     double *dd_at;
     double *d_atpk1;
     double *d_atmk1;
     double *d_u_wface;

     double *dd_v_wface;
     double *d_ws;
     double *d_wspk1;
     double *d_tarea_r;

     double *dd_hue;
     double *dd_hts;
  double d_dts = Para_s->d_dts;
  double* d_odzt = Para_s->d_odzt;
  //double* d_adv_tt = Para_s->d_adv_tt;
  double d_adv_tt[imt];
  double adv_xy1,adv_xy2,adv_xy3,adv_xy4;
  double adv_zz,adv_za,adv_zb1,adv_zb2;
  double adv_x0,adv_y0,adv_xx,adv_yy;
  double adv_c1,adv_c2,adv_zc;

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

     int i,j,k;

     for(int block = blockStart; block < blockEnd; block++) {
         blockId_z = block / (Para_s->cuGrid[0] * Para_s->cuGrid[1]);
         blockId_y = (block - blockId_z * (Para_s->cuGrid[0] * Para_s->cuGrid[1])) / Para_s->cuGrid[0];
         blockId_x = block % Para_s->cuGrid[0];
         j = blockId_y * Para_s->cuBlock[1];
         k = blockId_z * Para_s->cuBlock[2];

//         CRTS_dma_get(d_at00pk1, &Para_s->d_at00[(k+1) * imt * jmt + j*imt], sizeof(double) * Para_s->cuBlock[0]);
//         CRTS_dma_get(d_at00mk1, &Para_s->d_at00[(k-1) * imt * jmt + j*imt], sizeof(double) * Para_s->cuBlock[0]);
//         CRTS_dma_get(d_atmaxpk1, &Para_s->d_atmax[(k+1) * imt * jmt + j*imt], sizeof(double) * Para_s->cuBlock[0]);
//         CRTS_dma_get(d_atmaxmk1, &Para_s->d_atmax[(k-1) * imt * jmt + j*imt], sizeof(double) * Para_s->cuBlock[0]);
//         CRTS_dma_get(d_atminpk1, &Para_s->d_atmin[(k+1) * imt * jmt + j*imt], sizeof(double) * Para_s->cuBlock[0]);
//         CRTS_dma_get(d_atminmk1, &Para_s->d_atmin[(k-1) * imt * jmt + j*imt], sizeof(double) * Para_s->cuBlock[0]);
//         CRTS_dma_get(d_atpk1, &Para_s->d_at[(k+1) * imt * jmt + j*imt], sizeof(double) * Para_s->cuBlock[0]);
//         CRTS_dma_get(d_atmk1, &Para_s->d_at[(k-1) * imt * jmt + j*imt], sizeof(double) * Para_s->cuBlock[0]);
//         CRTS_dma_get(d_u_wface, &Para_s->d_u_wface[k * imt * jmt + j*imt], sizeof(double) * Para_s->cuBlock[0]);
//         CRTS_dma_get(d_ws, &Para_s->d_ws[k * jmt * imt + j * imt], sizeof(double) * Para_s->cuBlock[0]);
//         CRTS_dma_get(d_wspk1, &Para_s->d_ws[(k + 1) * jmt * imt + j * imt], sizeof(double) * Para_s->cuBlock[0]);
//         CRTS_dma_get(d_tarea_r, &Para_s->d_tarea_r[j*imt], sizeof(double) * Para_s->cuBlock[0]);
//
//         CRTS_dma_get(dd_at00, &Para_s->d_at00[k * imt * jmt + (j-1)*imt], sizeof(double) * Para_s->cuBlock[0] * 3);
//         CRTS_dma_get(dd_atmax, &Para_s->d_atmax[k * imt * jmt + (j-1)*imt], sizeof(double) * Para_s->cuBlock[0] * 3);
//         CRTS_dma_get(dd_atmin, &Para_s->d_atmin[k * imt * jmt + (j-1)*imt], sizeof(double) * Para_s->cuBlock[0] * 3);
//         CRTS_dma_get(dd_at, &Para_s->d_at[k * imt * jmt + (j-1)*imt], sizeof(double) * Para_s->cuBlock[0] * 3);
//         CRTS_dma_get(dd_v_wface, &Para_s->d_v_sface[k * imt * jmt + (j-1)*imt], sizeof(double) * Para_s->cuBlock[0] * 2);
//         CRTS_dma_get(dd_hue, &Para_s->d_hue[(j-1)*imt], sizeof(double) * Para_s->cuBlock[0] * 2);
//         CRTS_dma_get(dd_hts, &Para_s->d_hts[(j-1)*imt], sizeof(double) * Para_s->cuBlock[0] * 2);

         d_at00pk1 = &Para_s->d_at00[(k+1) * imt * jmt + j*imt];
         d_at00mk1 = &Para_s->d_at00[(k-1) * imt * jmt + j*imt];
         d_atmaxpk1 =  &Para_s->d_atmax[(k+1) * imt * jmt + j*imt];
         d_atmaxmk1 = &Para_s->d_atmax[(k-1) * imt * jmt + j*imt];
         d_atminpk1 = &Para_s->d_atmin[(k+1) * imt * jmt + j*imt];
         d_atminmk1 = &Para_s->d_atmin[(k-1) * imt * jmt + j*imt];
         d_atpk1 = &Para_s->d_at[(k+1) * imt * jmt + j*imt];
         d_atmk1 = &Para_s->d_at[(k-1) * imt * jmt + j*imt];
         d_u_wface = &Para_s->d_u_wface[k * imt * jmt + j*imt];
         d_ws = &Para_s->d_ws[k * jmt * imt + j * imt];
         d_wspk1 = &Para_s->d_ws[(k + 1) * jmt * imt + j * imt];
         d_tarea_r = &Para_s->d_tarea_r[j*imt];
         dd_at00 = &Para_s->d_at00[k * imt * jmt + (j-1)*imt];
         dd_atmax = &Para_s->d_atmax[k * imt * jmt + (j-1)*imt];
         dd_atmin =  &Para_s->d_atmin[k * imt * jmt + (j-1)*imt];
         dd_at = &Para_s->d_at[k * imt * jmt + (j-1)*imt];
         dd_v_wface = &Para_s->d_v_sface[k * imt * jmt + (j-1)*imt];
         dd_hue = &Para_s->d_hue[(j-1)*imt];
         dd_hts = &Para_s->d_hts[(j-1)*imt];

         for (i = 0; i < threadSum; i++) {
             d_adv_tt[i] = 0.0;
             if (1 <= i && i < (imt - 1) && 1 <= j && j < (jmt - 1) && k==0) {

                 if (dd_at00[imt+i] > dd_atmax[imt+i] || dd_at00[imt+i] < dd_atmin[imt+i] ||
                     ((dd_at00[imt+i + 1] > dd_atmax[imt+i+1] || dd_at00[imt+i + 1] < dd_atmin[imt+i+1]) && i <= imt - 3)) {
                     adv_xy1 = -(dd_at[imt+i + 1] -dd_at[imt+i]) * fabs(d_u_wface[i + 1]) * d_tarea_r[i];
                 } else {
                     adv_xy1 = -d_dts * (dd_at[imt+i + 1] - dd_at[imt+i]) * 2.0 *
                               d_tarea_r[i] * d_u_wface[i + 1] * d_u_wface[i + 1] /
                               (Para_s->d_htw[j * imt + i + 1] * Para_s->d_hun[j * imt + i + 1]);
                 }
                 if (dd_at00[imt+i] > dd_atmax[imt+i] || dd_at00[imt+i] < dd_atmin[imt+i] ||
                     ((dd_at00[imt+i - 1] > dd_atmax[imt+i - 1] || dd_at00[imt+i - 1] < dd_atmin[imt+i - 1]) &&i >= 2)) {
                     adv_xy2 = (dd_at[imt+i] - dd_at[imt+i - 1]) * fabs(d_u_wface[i]) * d_tarea_r[i];
                 } else {
                     adv_xy2 = d_dts * (dd_at[imt+i] - dd_at[imt+i - 1]) * 2.0 *
                               d_tarea_r[i] * d_u_wface[i] * d_u_wface[i] /
                               (Para_s->d_htw[j * imt + i] * Para_s->d_hun[j * imt + i]);
                 }
                 if (dd_at00[imt+i] > dd_atmax[imt+i] || dd_at00[imt+i] < dd_atmin[imt+i] ||
                     ((dd_at00[2*imt+i] > dd_atmax[2*imt+i] || dd_at00[2*imt+i] < dd_atmin[2*imt+i]) && j <= jmt - 3)) {
                     adv_xy3 = -(dd_at[2*imt+i] - dd_at[imt+i]) * fabs(dd_v_wface[imt+i]) * d_tarea_r[i];
                 } else {
                     adv_xy3 = -d_dts * (dd_at[2*imt+i] - dd_at[imt+i]) * 2.0 *
                               d_tarea_r[i] * dd_v_wface[imt+i] * dd_v_wface[imt+i] /
                               (dd_hts[imt+i] * dd_hue[imt+i]);
                 }
                 if (dd_at00[imt+i] > dd_atmax[imt+i] || dd_at00[imt+i] < dd_atmin[imt+i] ||
                     ((dd_at00[i] > dd_atmax[i] || dd_at00[i] < dd_atmin[i]) && j >= 2)) {
                     adv_xy4 = (dd_at[imt+i] - dd_at[i]) * fabs(dd_v_wface[i]) * d_tarea_r[i];
                 } else {
                     adv_xy4 = d_dts * (dd_at[imt+i] - dd_at[i]) * 2.0 *
                               d_tarea_r[i] * dd_v_wface[i] * dd_v_wface[i] /
                               (dd_hts[i] * dd_hue[i]);
                 }
                 adv_zb1 = 0.0;
                 if (dd_at00[imt+i] > dd_atmax[imt+i] || dd_at00[imt+i] < dd_atmin[imt+i] ||
                     ((d_at00pk1[i] > d_atmaxpk1[i] ||
                       d_at00pk1[i] < d_atminpk1[i]) && k <= km - 2)) {
                     adv_zb2 = 0.5 * fabs(d_wspk1[i]) * Para_s->d_odzp[k] *
                               (dd_at[imt+i] - d_atpk1[i]);
                 } else {
                     adv_zb2 = 0.5 * Para_s->d_odzp[0] * Para_s->d_ws[1 * jmt * imt + j * imt + i] *
                               Para_s->d_ws[1 * jmt * imt + j * imt + i] * d_odzt[1] *
                               (Para_s->d_at[j * imt + i] - Para_s->d_at[1 * jmt * imt + j * imt + i]) * d_dts;
                 }
                 adv_za = -0.5 * Para_s->d_odzp[0] * Para_s->d_ws[1 * jmt * imt + j * imt + i] *
                          (Para_s->d_at[1 * jmt * imt + j * imt + i] + Para_s->d_at[j * imt + i]);
                 adv_zc = Para_s->d_odzp[0] * Para_s->d_at[j * imt + i] * Para_s->d_ws[1 * jmt * imt + j * imt + i];
                 adv_c1 = -dd_at[imt+i] * (d_u_wface[i + 1] - d_u_wface[i]) * d_tarea_r[i] * 2.0;
                 adv_c2 = -dd_at[imt+i] * (dd_v_wface[imt+i] - dd_v_wface[i]) * d_tarea_r[i] * 2.0;
                 adv_x0 = (dd_at[imt+i + 1] +dd_at[imt+i]) * d_u_wface[i + 1] * d_tarea_r[i] -
                          (dd_at[imt+i] + dd_at[imt+i - 1]) * d_u_wface[i] * d_tarea_r[i];
                 adv_y0 = ((dd_at[2*imt+i] + dd_at[imt+i]) * dd_v_wface[imt+i] -
                           (dd_at[imt+i] + dd_at[i]) * dd_v_wface[i]) * d_tarea_r[i];
                 adv_xx = -(adv_x0 + adv_xy1 + adv_xy2 + adv_c1);
                 adv_yy = -(adv_y0 + adv_xy3 + adv_xy4 + adv_c2);
                 adv_zz = -(adv_za + adv_zb1 + adv_zb2 + adv_zc);
                 //d_adv_tt[k * jmt * imt + j * imt + i] = adv_xx + adv_yy + adv_zz;
                 d_adv_tt[i] = adv_xx + adv_yy + adv_zz;


             }
             if (k >= 1 && k < km - 1) {
                 if (1 <= j && j < (jmt - 1) && 1 <= i && i < (imt - 1)) {
                     if ((dd_at00[imt+i] > dd_atmax[imt+i] || dd_at00[imt+i] < dd_atmin[imt+i] ||
                          dd_at00[imt+i + 1] > dd_atmax[imt+i+1] || dd_at00[imt+i + 1] < dd_atmin[imt+i+1]) && i <= imt - 3) {
                         adv_xy1 = -(dd_at[imt+i + 1] -dd_at[imt+i]) *fabs(d_u_wface[i + 1]) * d_tarea_r[i];
                     } else {
                         adv_xy1 = -d_dts * (dd_at[imt+i + 1] - dd_at[imt+i]) * 2.0 *
                                   d_tarea_r[i] * d_u_wface[i + 1] * d_u_wface[i + 1] /
                                   (Para_s->d_htw[j * imt + i + 1] * Para_s->d_hun[j * imt + i + 1]);
                     }
                     if ((dd_at00[imt+i] > dd_atmax[imt+i] || dd_at00[imt+i] < dd_atmin[imt+i] ||
                          dd_at00[imt+i - 1] > dd_atmax[imt+i - 1] || dd_at00[imt+i - 1] < dd_atmin[imt+i - 1]) && i >= 2) {
                         adv_xy2 = (dd_at[imt+i] - dd_at[imt+i - 1]) * fabs(d_u_wface[i]) * d_tarea_r[i];
                     } else {
                         adv_xy2 = d_dts * (dd_at[imt+i] - dd_at[imt+i - 1]) * 2.0 *
                                   d_tarea_r[i] * d_u_wface[i] * d_u_wface[i] /
                                   (Para_s->d_htw[j * imt + i] * Para_s->d_hun[j * imt + i]);
                     }
                     if ((dd_at00[imt+i] > dd_atmax[imt+i] || dd_at00[imt+i] < dd_atmin[imt+i] ||
                          dd_at00[2*imt+i] > dd_atmax[2*imt+i] || dd_at00[2*imt+i] < dd_atmin[2*imt+i]) && j <= jmt - 3) {
                         adv_xy3 = -(dd_at[2*imt+i] - dd_at[imt+i]) *
                                   fabs(dd_v_wface[imt+i]) * d_tarea_r[i];
                     } else {
                         adv_xy3 = -d_dts * (dd_at[2*imt+i] - dd_at[imt+i]) * 2.0 *
                                   d_tarea_r[i] * dd_v_wface[imt+i] * dd_v_wface[imt+i] /
                                   (dd_hts[imt+i] * dd_hue[imt+i]);
                     }
                     if ((dd_at00[imt+i] > dd_atmax[imt+i] || dd_at00[imt+i] < dd_atmin[imt+i] ||
                          dd_at00[i] > dd_atmax[i] || dd_at00[i] < dd_atmin[i]) && j >= 2) {
                         adv_xy4 = (dd_at[imt+i] - dd_at[i]) * fabs(dd_v_wface[i]) * d_tarea_r[i];
                     } else {
                         adv_xy4 = d_dts * (dd_at[imt+i] - dd_at[i]) * 2.0 *
                                   d_tarea_r[i] * dd_v_wface[i] * dd_v_wface[i] /
                                   (dd_hts[i] * dd_hue[i]);
                     }
                     if ((dd_at00[imt+i] > dd_atmax[imt+i] || dd_at00[imt+i] < dd_atmin[imt+i] ||
                          d_at00mk1[i] > d_atmaxmk1[i] || d_at00mk1[i] < d_atminmk1[i]) && k >= 1) {
                         adv_zb1 = -0.5 * fabs(d_ws[i]) * Para_s->d_odzp[k] * (d_atmk1[i] - dd_at[imt+i]);
                     } else {
                         adv_zb1 = -0.5 * Para_s->d_odzp[k] * d_ws[i] *
                                   d_ws[i] * d_odzt[k] * (d_atmk1[i] - dd_at[imt+i]) * d_dts;
                     }
                     if ((dd_at00[imt+i] > dd_atmax[imt+i] || dd_at00[imt+i] < dd_atmin[imt+i] ||
                          d_at00pk1[i] > d_atmaxpk1[i] || d_at00pk1[i] < d_atminpk1[i]) && k <= km - 2) {
                         adv_zb2 = 0.5 * fabs(d_wspk1[i]) * Para_s->d_odzp[k] * (dd_at[imt+i] - d_atpk1[i]);
                     } else {
                         adv_zb2 = 0.5 * Para_s->d_odzp[k] * d_wspk1[i] * d_wspk1[i] * d_odzt[k + 1] *
                                   (dd_at[imt+i] - d_atpk1[i]) * d_dts;
                     }
                     adv_c1 = -dd_at[imt+i] * (d_u_wface[i + 1] - d_u_wface[i]) * d_tarea_r[i] * 2.0;
                     adv_c2 = -dd_at[imt+i] * (dd_v_wface[imt+i] - dd_v_wface[i]) * d_tarea_r[i] * 2.0;
                     adv_za = 0.5 * Para_s->d_odzp[k] * d_ws[i] *
                              (dd_at[imt+i] + d_atmk1[i]) - 0.5 * Para_s->d_odzp[k] * d_wspk1[i] *
                                                       (dd_at[imt+i] + d_atpk1[i]);
                     adv_zc = -Para_s->d_odzp[k] *dd_at[imt+i] * (d_ws[i] - d_wspk1[i]);
                     adv_x0 = (dd_at[imt+i + 1] + dd_at[imt+i]) * d_u_wface[i + 1] * d_tarea_r[i] -
                              (dd_at[imt+i] + dd_at[imt+i - 1]) * d_u_wface[i] * d_tarea_r[i];
                     adv_y0 = ((dd_at[2*imt+i] + dd_at[imt+i]) * dd_v_wface[imt+i] -
                               (dd_at[imt+i] + dd_at[i]) * dd_v_wface[i]) * d_tarea_r[i];
                     adv_xx = -(adv_x0 + adv_xy1 + adv_xy2 + adv_c1);
                     adv_yy = -(adv_y0 + adv_xy3 + adv_xy4 + adv_c2);
                     adv_zz = -(adv_za + adv_zb1 + adv_zb2 + adv_zc);
                     // d_adv_tt[k * jmt * imt + j * imt + i] = adv_xx + adv_yy + adv_zz;
                     d_adv_tt[i] = adv_xx + adv_yy + adv_zz;
                 }
             }
             if (k == km - 1) {
                 if (dd_at00[imt+i] > dd_atmax[imt+i] || dd_at00[imt+i] < dd_atmin[imt+i] ||
                     ((dd_at00[imt+i + 1] > dd_atmax[imt+i+1] || dd_at00[imt+i + 1] < dd_atmin[imt+i+1]) && i <= imt - 3)) {
                     adv_xy1 = -(dd_at[imt+i + 1] - dd_at[imt+i]) *
                               fabs(d_u_wface[i + 1]) * d_tarea_r[i];
                 } else {
                     adv_xy1 = -d_dts * (dd_at[imt+i + 1] - dd_at[imt+i]) * 2.0 *
                               d_tarea_r[i] * d_u_wface[i + 1] *
                               d_u_wface[i + 1] / (Para_s->d_htw[j * imt + i + 1] * Para_s->d_hun[j * imt + i + 1]);
                 }
                 if (dd_at00[imt+i] > dd_atmax[imt+i] || dd_at00[imt+i] < dd_atmin[imt+i] ||
                     ((dd_at00[imt+i - 1] > dd_atmax[imt+i - 1] || dd_at00[imt+i - 1] < dd_atmin[imt+i - 1]) && i >= 2)) {
                     adv_xy2 = (dd_at[imt+i] - dd_at[imt+i - 1]) * fabs(d_u_wface[i]) * d_tarea_r[i];
                 } else {
                     adv_xy2 = d_dts * (dd_at[imt+i] - dd_at[imt+i - 1]) * 2.0 *
                               d_tarea_r[i] * d_u_wface[i] *
                               d_u_wface[i] / (Para_s->d_htw[j * imt + i] * Para_s->d_hun[j * imt + i]);
                 }
                 if (dd_at00[imt+i] > dd_atmax[imt+i] || dd_at00[imt+i] < dd_atmin[imt+i] ||
                     ((dd_at00[2*imt+i] > dd_atmax[2*imt+i] || dd_at00[2*imt+i] < dd_atmin[2*imt+i]) && j <= jmt - 3)) {
                     adv_xy3 = -(dd_at[2*imt+i] - dd_at[imt+i]) *
                               fabs(dd_v_wface[imt+i]) * d_tarea_r[i];
                 } else {
                     adv_xy3 = -d_dts * (dd_at[2*imt+i] - dd_at[imt+i]) * 2.0 *
                               d_tarea_r[i] * dd_v_wface[imt+i] * dd_v_wface[imt+i] /
                               (dd_hts[imt+i] * dd_hue[imt+i]);
                 }
                 if (dd_at00[imt+i] > dd_atmax[imt+i] || dd_at00[imt+i] < dd_atmin[imt+i] ||
                     ((dd_at00[i] > dd_atmax[i] || dd_at00[i] < dd_atmin[i]) && j >= 2)) {
                     adv_xy4 = (dd_at[imt+i] - dd_at[i]) *
                               fabs(dd_v_wface[i]) * d_tarea_r[i];
                 } else {
                     adv_xy4 = d_dts * (dd_at[imt+i] - dd_at[i]) * 2.0 *
                               d_tarea_r[i] * dd_v_wface[i] * dd_v_wface[i] /
                               (dd_hts[i] * dd_hue[i]);
                 }
                 if (dd_at00[imt+i] > dd_atmax[imt+i] || dd_at00[imt+i] < dd_atmin[imt+i] ||
                     d_at00mk1[i] > d_atmaxmk1[i] || d_at00mk1[i] < d_atminmk1[i]) {
                     adv_zb1 = -0.5 * fabs(d_ws[i]) * Para_s->d_odzp[k] * (d_atmk1[i] - dd_at[imt+i]);
                 } else {
                     adv_zb1 = -0.5 * Para_s->d_odzp[k] * d_ws[i] * d_ws[i] * d_odzt[k] *
                               (d_atmk1[i] - dd_at[imt+i]) * d_dts;
                 }
                 adv_zb2 = 0.0;
                 adv_c1 = -dd_at[imt+i] * (d_u_wface[i + 1] - d_u_wface[i]) * d_tarea_r[i] * 2.0;
                 adv_c2 = -dd_at[imt+i] * (dd_v_wface[imt+i] - dd_v_wface[i]) * d_tarea_r[i] * 2.0;
                 adv_za = 0.5 * Para_s->d_odzp[km - 1] * Para_s->d_ws[(km - 1) * jmt * imt + j * imt + i] *
                          (Para_s->d_at[(km - 1) * jmt * imt + j * imt + i] +
                           Para_s->d_at[(km - 2) * jmt * imt + j * imt + i]);
                 adv_zc = -Para_s->d_odzp[km - 1] * Para_s->d_at[(km - 1) * jmt * imt + j * imt + i] *
                          Para_s->d_ws[(km - 1) * jmt * imt + j * imt + i];
                 adv_x0 = (dd_at[imt+i + 1] +dd_at[imt+i]) * d_u_wface[i + 1] * d_tarea_r[i] -
                          (dd_at[imt+i] + dd_at[imt+i - 1]) * d_u_wface[i] * d_tarea_r[i];
                 adv_y0 = ((dd_at[2*imt+i] + dd_at[imt+i]) * dd_v_wface[imt+i] - (dd_at[imt+i] + dd_at[i]) *
                                                                   dd_v_wface[i]) * d_tarea_r[i];
                 adv_xx = -(adv_x0 + adv_xy1 + adv_xy2 + adv_c1);
                 adv_yy = -(adv_y0 + adv_xy3 + adv_xy4 + adv_c2);
                 adv_zz = -(adv_za + adv_zb1 + adv_zb2 + adv_zc);
                 // d_adv_tt[k * jmt * imt + j * imt + i] = adv_xx + adv_yy + adv_zz;
                     d_adv_tt[i] = adv_xx + adv_yy + adv_zz;
             }
         }

         CRTS_dma_put(&Para_s->d_adv_tt[k * imt * jmt + j * imt], d_adv_tt, sizeof(double) * Para_s->cuBlock[0]);
     }


#ifdef ADV_DEBUG 
     penv_slave1_dcache_access_count(&counte);
    penv_slave2_dcache_miss_count(&misscounte);
    if(Para_s->rank == 0)
        printf("[slave %d]:advection_3 count:%d,counte:%d,misscount:%d,misscounte:%d,total count:%d, total misscount:%d, misscount rate:%f\n",
		blockIndex,count,counte,misscount,misscounte,counte-count, misscounte-misscount,
		 (float)(misscounte-misscount)/(float)(counte-count));
#endif
    ldm_free(Para_s, sizeof(advection2_Para));
}



void QRLFLX_PRM_NQGRID( int KXYZ, int IXYZ, double *VQ, double *DltXYZ, double *Suvw, double Dtsc,  double *flxR, double *flxL){
    double COF;
    double  deltX1, deltX2, deltX3, deltX4;
    double  deltX5, deltX6, deltX7;

    double  dltc, dltm, dltp;
    double  dltca, dltma, dltpa;

    double  VQmin, VQmax;
    int  N;

    double uhafR, uhafL;

    double FMIU;
    double Atmp, Btmp;
    double Ctmp, Dtmp;
    double Gama, Beta;
    double Fmin, Fmax, Fave;

    double EPStmp;
    double dltQ[IXYZ];
    double qR[IXYZ];
    double qL[IXYZ];
    double minval;
    double signval;
    COF = 3.0e0;
    EPStmp = 1.0e-20;

    //DO N = 2, IXYZ-1
    for( N = 1; N<IXYZ-1; N++) {

        deltX1 = DltXYZ[N] / (DltXYZ[N - 1] + DltXYZ[N] + DltXYZ[N + 1]);
        deltX2 = (2.0e0 * DltXYZ[N - 1] + DltXYZ[N]) / (DltXYZ[N + 1] + DltXYZ[N]);
        deltX3 = (DltXYZ[N] + 2.0e0 * DltXYZ[N + 1]) / (DltXYZ[N - 1] + DltXYZ[N]);
        dltc = deltX1 * (deltX2 * (VQ[N + 1] - VQ[N]) + deltX3 * (VQ[N] - VQ[N - 1]));
        dltm = VQ[N] - VQ[N - 1];
        dltp = VQ[N + 1] - VQ[N];

        if((dltm * dltp) <= 0.0e0 )
            dltQ[N] = 0.0e0;
        else {
            dltma = COF * fabs(dltm);
            dltpa = COF * fabs(dltp);
            dltca = fabs(dltc);
            //dltQ(N) = DMIN1(dltca, dltma, dltpa) * DSIGN(1.0D0, dltc)
            minval = findMin(dltca, dltma, dltpa);
            if(dltc >= 0.0)
                signval = 1.0e0;
            else
                signval = -1.0e0;
            dltQ[N] = minval * signval;
        }
    }

    dltQ[0] = 0.0e0;
    dltQ[IXYZ-1] = 0.0e0;

    // DO N = 2, IXYZ-2
    for(N = 1; N < IXYZ-2; N++) {
        deltX1 = DltXYZ[N] / (DltXYZ[N] + DltXYZ[N + 1]);
        deltX2 = 1.0e0 / (DltXYZ[N - 1] + DltXYZ[N] + DltXYZ[N + 1] + DltXYZ[N + 2]);
        deltX3 = (2.0e0 * DltXYZ[N + 1] * DltXYZ[N]) / (DltXYZ[N] + DltXYZ[N + 1]);
        deltX4 = (DltXYZ[N - 1] + DltXYZ[N]) / (2.0e0 * DltXYZ[N] + DltXYZ[N + 1]);
        deltX5 = (DltXYZ[N + 2] + DltXYZ[N + 1]) / (2.0e0 * DltXYZ[N + 1] + DltXYZ[N]);
        deltX6 = DltXYZ[N] * (DltXYZ[N - 1] + DltXYZ[N]) / (2.0e0 * DltXYZ[N] + DltXYZ[N + 1]);
        deltX7 = DltXYZ[N + 1] * (DltXYZ[N + 1] + DltXYZ[N + 2]) / (DltXYZ[N] + 2.0e0 * DltXYZ[N + 1]);

        qR[N] = VQ[N] + deltX1 * (VQ[N + 1] - VQ[N]) + deltX2 * (deltX3 * (deltX4 - deltX5) * (VQ[N + 1] - VQ[N]) -
                                                                 deltX6 *dltQ[N + 1] + deltX7 * dltQ[N]);

        VQmin = findMin(VQ[N - 1], VQ[N], VQ[N + 1]);
        VQmax = findMax(VQ[N - 1], VQ[N], VQ[N + 1]);
        qR[N] = fmax(VQmin, qR[N]);
        //qR[N] = DMIN1(VQmax, qR(N))
        qR[N] = fmin(VQmax, qR[N]);
    }

    qR[0] = 0.5e0*(VQ[0]+VQ[1]);
    qR[IXYZ-1] = VQ[IXYZ-1];
    qR[IXYZ-2] = 0.5e0*(VQ[IXYZ-2]+VQ[IXYZ-1]);

    //DO N = 2, IXYZ
    for( N = 1; N < IXYZ; N++)
        qL[N] = qR[N-1];

    qL[0] = VQ[0];


    //To Compute the Flux Through the Cell

    //DO N = 1, IXYZ
    for( N = 0; N<IXYZ; N++) {

        if(N == 0 ){

            uhafR = Suvw[N];
            uhafL = 0.0e0;

            if(uhafR > 0.0e0 ) {
                FMIU = fmin(1.0e0 * DltXYZ[N], fabs(uhafR * Dtsc));
                Atmp = qR[N];

                Fmin = fmin(qR[N], qL[N]);
                Fmax = fmax(qR[N], qL[N]);
                Fave = fmax(fmin(VQ[N], Fmax), Fmin);
                Ctmp = fabs(qR[N] - Fave);
                Dtmp = fabs(Fave - qL[N]);
                Gama = (Ctmp + EPStmp) / (Dtmp + EPStmp);
                Beta = (Gama - 1.0e0) / (-DltXYZ[N]);
                Btmp = (Gama * Fave - qR[N]) / (-DltXYZ[N]);

                flxR[N] = (Atmp * FMIU - Btmp * FMIU * FMIU) / (1.0e0 - Beta * FMIU);
            }
            else if (uhafR == 0.0e0)
                flxR[N] = 0.0e0;
            else{

                FMIU = fabs(-uhafR * Dtsc);

                Atmp = qR[N];

                Fmin = fmin(qR[N + 1], qL[N + 1]);
                Fmax = fmax(qR[N + 1], qL[N + 1]);
                Fave = fmax(fmin(VQ[N + 1], Fmax), Fmin);
                Ctmp = fabs(qR[N] - Fave);
                Dtmp = fabs(Fave - qR[N + 1]);
                Gama = (Ctmp + EPStmp) / (Dtmp + EPStmp);
                Beta = (Gama - 1.0e0) / DltXYZ[N + 1];
                Btmp = (Gama * Fave - qR[N]) / DltXYZ[N + 1];

                flxR[N] = -1.0e0 * (Atmp * FMIU + Btmp * FMIU * FMIU) / (1.0e0 + Beta * FMIU);
            }

                flxL[N] = 0.0e0;
        }
            //N.GE.2 .AND. N.LE.(IXYZ-1)
        else if(N >= 1 && (N < (IXYZ - 1))) {

            uhafR = Suvw[N];
            uhafL = Suvw[N - 1];

            if (uhafR > 0.0e0) {

                FMIU = fmin(1.0e0 * DltXYZ[N], fabs(uhafR * Dtsc));

                Atmp = qR[N];

                Fmin = fmin(qR[N], qL[N]);
                Fmax = fmax(qR[N], qL[N]);
                Fave = fmax(fmin(VQ[N], Fmax), Fmin);
                Ctmp = fabs(qR[N] - Fave);
                Dtmp = fabs(Fave - qL[N]);
                Gama = (Ctmp + EPStmp) / (Dtmp + EPStmp);
                Beta = (Gama - 1.0e0) / (-DltXYZ[N]);
                Btmp = (Gama * Fave - qR[N]) / (-DltXYZ[N]);

                flxR[N] = (Atmp * FMIU - Btmp * FMIU * FMIU) / (1.0e0 - Beta * FMIU);
            } else if (uhafR == 0.0e0)
                flxR[N] = 0.0e0;

            else {
                FMIU =  fabs(-uhafR * Dtsc);

                Atmp = qR[N];

                Fmin = fmin(qR[N + 1], qL[N + 1]);
                Fmax = fmax(qR[N + 1], qL[N + 1]);
                Fave = fmax(fmin(VQ[N + 1], Fmax), Fmin);
                Ctmp = fabs(qR[N] - Fave);
                Dtmp = fabs(Fave - qR[N + 1]);
                Gama = (Ctmp + EPStmp) / (Dtmp + EPStmp);
                Beta = (Gama - 1.0e0) / DltXYZ[N + 1];
                Btmp = (Gama * Fave - qR[N]) / DltXYZ[N + 1];

                flxR[N] = -1.0e0 * (Atmp * FMIU + Btmp * FMIU * FMIU) / (1.0e0 + Beta * FMIU);

            }

            if (uhafL > 0.0e0) {

                FMIU =  fabs(uhafL * Dtsc);

                Atmp = qL[N];

                Fmin = fmin(qL[N - 1], qR[N - 1]);
                Fmax = fmax(qL[N - 1], qR[N - 1]);
                Fave = fmax(fmin(VQ[N - 1], Fmax), Fmin);
                Ctmp = fabs(qL[N] - Fave);
                Dtmp = fabs(Fave - qL[N - 1]);
                Gama = (Ctmp + EPStmp) / (Dtmp + EPStmp);
                Beta = (Gama - 1.0e0) / (-DltXYZ[N - 1]);
                Btmp = (Gama * Fave - qL[N]) / (-DltXYZ[N - 1]);

                flxL[N] = (Atmp * FMIU - Btmp * FMIU * FMIU) / (1.0e0 - Beta * FMIU);
            } else if (uhafL == 0.0e0)
                flxL[N] = 0.0e0;

            else {
                FMIU =  fabs(-uhafL * Dtsc);
                Atmp = qL[N];


                Fmin = fmin(qL[N], qR[N]);
                Fmax = fmax(qL[N], qR[N]);
                Fave = fmax(fmin(VQ[N], Fmax), Fmin);
                Ctmp = fabs(qL[N] - Fave);
                Dtmp = fabs(Fave - qR[N]);
                Gama = (Ctmp + EPStmp) / (Dtmp + EPStmp);
                Beta = (Gama - 1.0e0) / DltXYZ[N];
                Btmp = (Gama * Fave - qL[N]) / DltXYZ[N];

                flxL[N] = -1.0e0 * (Atmp * FMIU + Btmp * FMIU * FMIU) / (1.0e0 + Beta * FMIU);

            }
        }
        else if(N == IXYZ-1) {

            uhafR = 0.0e0;
            uhafL = Suvw[N - 1];

            flxR[N] = 0.0e0;

            if (uhafL > 0.0e0) {

                FMIU = fabs(uhafL * Dtsc);

                Atmp = qL[N];

                Fmin = fmin(qL[N - 1], qR[N - 1]);
                Fmax = fmax(qL[N - 1], qR[N - 1]);
                Fave = fmax(fmin(VQ[N - 1], Fmax), Fmin);
                Ctmp = fabs(qL[N] - Fave);
                Dtmp = fabs(Fave - qL[N - 1]);
                Gama = (Ctmp + EPStmp) / (Dtmp + EPStmp);
                Beta = (Gama - 1.0e0) / (-DltXYZ[N - 1]);
                Btmp = (Gama * Fave - qL[N]) / (-DltXYZ[N - 1]);

                flxL[N] = (Atmp * FMIU - Btmp * FMIU * FMIU) / (1.0e0 - Beta * FMIU);
            } else if (uhafL == 0.0e0)
                flxL[N] = 0.0e0;

            else {

                FMIU =  fabs(-uhafL * Dtsc);

                Atmp = qL[N];

                Fmin = fmin(qL[N], qR[N]);
                Fmax = fmax(qL[N], qR[N]);
                Fave = fmax(fmin(VQ[N], Fmax), Fmin);
                Ctmp = fabs(qL[N] - Fave);
                Dtmp = fabs(Fave - qR[N]);
                Gama = (Ctmp + EPStmp) / (Dtmp + EPStmp);
                Beta = (Gama - 1.0e0) / DltXYZ[N];
                Btmp = (Gama * Fave - qL[N]) / DltXYZ[N];

                flxL[N] = -1.0e0 * (Atmp * FMIU + Btmp * FMIU * FMIU) / (1.0e0 + Beta * FMIU);

            }

        }
    }
}
void  advection_tracer_4_(advection2_Para *Paras) {

    int i,j,k;

    double VQC[imt];
    double DXYZ[imt];
    double UVWhaf[imt];
    double flxqR[imt];
    double flxqL[imt];
    int INDXYZ = 1;
    double FPAVE,FCAVE,CMIU;
    double uhafR, uhafL;
    //ATHREAD basic operation
    blockIndex = CRTS_tid;
    advection2_Para *Para_s = (advection2_Para *) ldm_malloc(sizeof(advection2_Para));
    CRTS_dma_get(Para_s, Paras, sizeof(advection2_Para));

    double d_dts = Para_s->d_dts;

    int d_mtracer = Para_s->d_mtracer;

    double* d_adv_tt = Para_s->d_adv_tt;
    double* d_adv_x = Para_s->d_adv_x;

    double *dd_at = Para_s->d_at;
    double *d_dxtr = Para_s->d_dxtr;
    for(j = 1; j < jmt; j++){
       for(k =blockIndex; k < km; k+=64 ){

          for(i = 0; i < imt; i++){

              //VQC(I) = VARTMP(I,J,K)
              //DXYZ(I) = DXT(I,J,iblock)
              VQC[i] = dd_at[k * imt * jmt + j * imt + i];
              DXYZ[i] = Para_s->d_dxt[j*imt + i];
              //UVWhaf(I) = 0.50D0*(UPRM(I+1,J-1,K)+UPRM(I+1,J,K))
              if(i<(imt-1))
                  UVWhaf[i] = 0.5e0 * (Para_s->d_wkd[k*imt*jmt + (j-1)*imt + i + 1] +
                          Para_s->d_wkd[k*imt*jmt + j*imt + i + 1] );
              else
                  UVWhaf[i] = 0.0e0;

          }
           INDXYZ = 1;
           QRLFLX_PRM_NQGRID( INDXYZ, imt, VQC, DXYZ, UVWhaf, d_dts, flxqR, flxqL);

           for(i = 0; i < imt; i++){

               if ( i == 0 ) {

                   uhafR = UVWhaf[i];
                   uhafL = 0.0e0;

               }
               else if ( (i >= 1) && (i < (imt-1)) ) {

                   uhafR = UVWhaf[i];
                   uhafL = UVWhaf[i-1];

               }
               else if ( i == (imt-1) ) {

                   uhafR = 0.0e0;
                   uhafL = UVWhaf[i-1];

               }
//               Para_s->d_flxqLtmp[k*imt*jmt + j*imt + i] = flxqL[i];
//               Para_s->d_flxqRtmp[k*imt*jmt + j*imt + i] = flxqR[i];

               CMIU = d_dxtr[j*imt+i];

               d_adv_x[k*imt*jmt + j*imt + i] = (CMIU * (flxqL[i] - flxqR[i]) + dd_at[k * imt * jmt + j * imt + i] *
                      (uhafR - uhafL) * d_dts *CMIU) / d_dts;
           }
       }
    }
    ldm_free(Para_s, sizeof(advection2_Para));
}


void  advection_tracer_5_(advection2_Para *Paras) {

    int i,j,k;

    double VQC[jmt];
    double DXYZ[jmt];
    double UVWhaf[jmt];
    double flxqR[jmt];
    double flxqL[jmt];
    int INDXYZ = 2;
    double FPAVE,FCAVE,CMIU;
    double uhafR, uhafL;
    //ATHREAD basic operation
    blockIndex = CRTS_tid;
    advection2_Para *Para_s = (advection2_Para *) ldm_malloc(sizeof(advection2_Para));
    CRTS_dma_get(Para_s, Paras, sizeof(advection2_Para));

    double d_dts = Para_s->d_dts;
    double *dd_at = Para_s->d_at;
    int d_mtracer = Para_s->d_mtracer;

    double* d_adv_tt = Para_s->d_adv_tt;
    double* d_adv_y = Para_s->d_adv_y;
    double *d_dytr = Para_s->d_dytr;
    for(i = 0; i < imt-1; i++){
        for(k =blockIndex; k < km; k += 64 ){
            for(j = 0; j < jmt; j++){
                //VQC(I) = VARTMP(I,J,K)
                // DXYZ(J) = DYT(I,J,iblock)
                VQC[j] = dd_at[k*imt*jmt + j*imt + i];
                DXYZ[j] = Para_s->d_dyt[j*imt + i];
                //UVWhaf(J) = 0.50D0*(VPRM(I,J,K)+VPRM(I+1,J,K))
                UVWhaf[j] = 0.5e0 * (Para_s->d_wkb[k*imt*jmt + j*imt + i] +
                        Para_s->d_wkb[k*imt*jmt + j*imt + i + 1] );
            }

            QRLFLX_PRM_NQGRID( INDXYZ, jmt, VQC, DXYZ, UVWhaf, d_dts, flxqR, flxqL);

            for(j = 0; j < jmt; j++){

                if ( j == 0 ) {

                    uhafR = UVWhaf[j];
                    uhafL = 0.0e0;

                }
                else if ( (j >= 1) && (j < (jmt-1)) ) {

                    uhafR = UVWhaf[j];
                    uhafL = UVWhaf[j-1];

                }
                else if ( j == (jmt-1) ) {

                    uhafR = 0.0e0;
                    uhafL = UVWhaf[j-1];
                }

                CMIU = d_dytr[j*imt+i];
                d_adv_y[k*imt*jmt + j*imt + i] = (CMIU * (flxqL[j] - flxqR[j]) + dd_at[k * imt * jmt + j * imt + i] *
                           (uhafR - uhafL) * d_dts *CMIU) / d_dts;
            }
        }
    }
    ldm_free(Para_s, sizeof(advection2_Para));
}

void  advection_tracer_6_(advection2_Para *Paras) {

    int i,j,k;

    double VQC[km];
    double DXYZ[km];
    double UVWhaf[km];
    double flxqR[km];
    double flxqL[km];
    int INDXYZ = 3;
    double FPAVE,FCAVE,CMIU;
    double uhafR, uhafL;
    //ATHREAD basic operation
    blockIndex = CRTS_tid;
    advection2_Para *Para_s = (advection2_Para *) ldm_malloc(sizeof(advection2_Para));
    CRTS_dma_get(Para_s, Paras, sizeof(advection2_Para));

    double d_dts = Para_s->d_dts;
    double *dd_at = Para_s->d_at;
    int d_mtracer = Para_s->d_mtracer;

    double* d_adv_tt = Para_s->d_adv_tt;
    double* d_adv_z = Para_s->d_adv_z;
    double *d_dytr = Para_s->d_dytr;
    double *d_ws = Para_s->d_ws;
    for(i = 0; i < imt; i++){
        for(j =blockIndex; j < jmt; j += 64 ){
            for(k = 0; k < km; k++){
                //VQC(I) = VARTMP(I,J,K)
                //  DXYZ(K) = ZKP(K)-ZKP(K+1)
                VQC[k] = dd_at[k*imt*jmt + j*imt + i];
                DXYZ[k] = Para_s->d_zkp[k] - Para_s->d_zkp[k+1];
                //UVWhaf(K) = -WPRM(I,J,K+1)
                if(k<km-1)
                    UVWhaf[k] = -d_ws[(k+1)*imt*jmt + j*imt + i];
                else
                    UVWhaf[k] = 0.0e0;
            }

            QRLFLX_PRM_NQGRID( INDXYZ, km, VQC, DXYZ, UVWhaf, d_dts, flxqR, flxqL);

            for(k = 0; k < km; k++){
                CMIU = 1.0e0/DXYZ[k];
                if ( k == 0 ) {

                    uhafR = UVWhaf[k];
                    uhafL = 0.0e0;
                }
                else {
                    uhafR = UVWhaf[k];
                    uhafL = UVWhaf[k-1];
                }
                //AZTMP(I,J,K) = (CMIU*( flxqL(K)-flxqR(K)) + VARTMP(I,J,K)*(uhafR-uhafL)*DTS/DXYZ(K))/DTS
                d_adv_z[k*imt*jmt + j*imt + i] = (CMIU*( flxqL[k]-flxqR[k]) +
                        dd_at[k * imt * jmt + j * imt + i]*(uhafR-uhafL)*d_dts * CMIU)/d_dts;

            }
        }
    }
    ldm_free(Para_s, sizeof(advection2_Para));
}

void  advection_tracer_7_(advection2_Para *Paras) {

    int i,j,k;
    //ATHREAD basic operation
    blockIndex = CRTS_tid;
    advection2_Para *Para_s = (advection2_Para *) ldm_malloc(sizeof(advection2_Para));
    CRTS_dma_get(Para_s, Paras, sizeof(advection2_Para));

    double localX[imt];
    double* d_adv_x = Para_s->d_adv_x;
    for(k = blockIndex; k < km; k += 64){
        CRTS_dma_get(localX, &Para_s->d_adv_x[k * imt * jmt +  imt], sizeof(double) * imt);
            for(i = 0; i < imt; i++){
                //AXTMP(I,1,K) = AXTMP(I,2,K)
                //d_adv_x[k*imt*jmt + i] = d_adv_x[k*imt*jmt + imt + i];
                d_adv_x[k*imt*jmt + i] = localX[i];
            }

    }
    ldm_free(Para_s, sizeof(advection2_Para));
}

void  advection_tracer_8_(advection2_Para *Paras) {

    int i,j,k;
    //ATHREAD basic operation
    blockIndex = CRTS_tid;
    advection2_Para *Para_s = (advection2_Para *) ldm_malloc(sizeof(advection2_Para));
    CRTS_dma_get(Para_s, Paras, sizeof(advection2_Para));

//    double localY[imt];

    double* d_adv_y = Para_s->d_adv_y;
    for(k = blockIndex; k < km; k += 64){
        for(j =0; j < jmt; j ++ ){
            //AYTMP(IMT,J,K) = AYTMP(IMT-1,J,K)
            d_adv_y[k*imt*jmt + j*imt + imt -1] = d_adv_y[k*imt*jmt + j*imt + imt -2];
        }
    }
    ldm_free(Para_s, sizeof(advection2_Para));
}

void  advection_tracer_9_(advection2_Para *Paras) {

    int i,j,k;
    //ATHREAD basic operation
    blockIndex = CRTS_tid;
    advection2_Para *Para_s = (advection2_Para *) ldm_malloc(sizeof(advection2_Para));
    CRTS_dma_get(Para_s, Paras, sizeof(advection2_Para));


    double* d_adv_tt = Para_s->d_adv_tt;
    double* d_adv_x = Para_s->d_adv_x;
    double* d_adv_y = Para_s->d_adv_y;
    double* d_adv_z = Para_s->d_adv_z;
    for(k = 0; k < km; k++){
        for(j =blockIndex; j < jmt; j += 64 ){
            for(i = 0; i < imt; i++){
                //DLTMP(I,J,K) = AZTMP(I,J,K)+AXTMP(I,J,K)+AYTMP(I,J,K)
                d_adv_tt[k*imt*jmt + j*imt + i] = d_adv_z[k*imt*jmt + j*imt + i] +
                                                  d_adv_x[k*imt*jmt + j*imt + i] + d_adv_y[k*imt*jmt + j*imt + i];
            }

        }
    }
    ldm_free(Para_s, sizeof(advection2_Para));
}