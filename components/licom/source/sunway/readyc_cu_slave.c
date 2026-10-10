
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
#include "constant_mod.h"

//#define DEGtoRAD M_PI/180.0
typedef struct vinteg_kernelPara{
    int cuGrid[3];
    int cuBlock[3];
    double *d_dlu;
    double *d_dlv;
    double *d_dlub;
    double *d_dlvb;
    double *d_dzp;
    double *d_ohbu;
    double *d_viv;
}vinteg_Para;

typedef struct kernelPara{
    int cuGrid[3];
    int cuBlock[3];
	double *d_h0bl;
	double *d_h0bf;
	double *d_h0;
	double *d_up;
	double *d_vp;
	double *d_viv;
	double *d_vit;
	double *d_wp12;
	double *d_wp13;
	double *d_at0;
	double *d_atn;
	double *d_ate;
	double *d_atne;

}Para;

typedef struct readyc1_2Para{
	int cuGrid[3];
	int cuBlock[3];
	double *d_wp12;
	double *d_wp13;
	double *d_vit;
	double *d_s2t;
	double *d_riv1;
	double *d_riv2;
	double *d_odzt;
	double *d_rit;
	double *d_rict;
	double *d_up;
	double *d_vp;
	double *d_s2u;
	double *d_viv;
	double *d_riu;
	double *d_ric;
}readyc1_2_Para;

typedef struct readyc1_3Para{
	int cuGrid[3];
	int cuBlock[3];
	double *d_diff_back;
	double *d_diff_back_sh;
	double *d_diff_back_nh;
	double *d_tlat;
	double d_diff_back_coef_max;
	double d_diff_back_eq;
	double d_diff_back_coef;
}readyc1_3_Para;

typedef struct readyc2_Para{
	int cuGrid[3];
	int cuBlock[3];
	double *d_akmu;
	double *d_au0;
	double *d_akmt;
	double *d_aus;
	double *d_auw;
	double *d_ausw;
	double *d_viv;

}readyc2_Para;

typedef struct readyc3Para{
	int cuGrid[3];
	int cuBlock[3];
	double *d_wka;
	double *d_ws;
	double *d_au0;
	double *d_ausw;
	double *d_aus;
	double *d_auw;
	double *d_u_wface;
	double *d_v_sface;
	double *d_u;
	double *d_v;
	double *d_hue;
	double *d_hun;

}readyc3_Para;

typedef struct advection_momentumPara{
	int cuGrid[3];
	int cuBlock[3];
	double *d_dlu;
	double *d_dlv;
	double *d_odzp;
	double *d_u;
	double *d_uarea_r;
	double *d_u_wface;
	double *d_v;
	double *d_v_sface;
	double *d_wka;
	//double *d_utest;

}advection_momentum_Para;

typedef struct readyc4Para{
	int cuGrid[3];
	int cuBlock[3];
	int *d_kmu;
	double *d_sbcx;
	double *d_sbcy;
	double *d_su;
	double *d_sv;
	double *d_bbcx;
	double  d_od0;
	double  d_c0f;
	double *d_up;
	double *d_vp;
	double *d_snlat;
	double  d_cag;
	double  d_sag;
	double *d_dlub;
	double *d_dlvb;
	double *d_dzp;
	double *d_ohbu;
	double *d_dlu;
	double *d_viv;
	double *d_dlv;
	double *d_bbcy;

}readyc4_Para;

typedef struct readyc4_2Para{
	int cuGrid[3];
	int cuBlock[3];
        double d_od0;
	double d_sag;
	double d_cag;
	double *d_su;
	double *d_sv;
	double *d_up;
	double *d_vp;
	double *d_akmu;
	double *d_odzt;
	double *d_viv;
	double *d_wka;
	double *d_snlat;
	double *d_dlv;
	double *d_odzp;
	double *d_dlu;

}readyc4_2_Para;

typedef struct hdiffu_del2_kernelPara{
        int cuBlock[3];
        int cuGrid[3];
        int d_jb;
        int d_je;
        int d_ib;
        int d_ie;
        int rank;
        int *d_kmu;
        double d_am;
        double d_c0f;
        double *d_wka;
        double *d_up;
        double *d_vp;
        double *d_duc;
        double *d_dum;
        double *d_dun;
        double *d_dus;
        double *d_due;
        double *d_duw;
        double *d_dmc;
        double *d_dmn;
        double *d_dms;
        double *d_dme;
        double *d_dmw;
        double *d_dlu;
        double *d_dlv;
	double *d_viv;
}hdiffu_del2_Para;


//dynamic atomic definition
#define CP_CORE      64
#define RCV_CORE     0
#define DATAPAKAGE   1

//ATHREAD global definition 
__thread_local int warp_size = 32; 
//__thread_local int blockIndex;
//__thread_local Para* Para_s;

 void  readyc1_cu_(Para *Paras){    

  Para* Para_s;

  //ATHREAD basic operation
  int blockId_x,blockId_y,blockId_z;
  int i,j,k;
  int blockIndex = CRTS_tid;
  Para_s = (Para*)ldm_malloc(sizeof(Para));
  CRTS_dma_get(Para_s, Paras, sizeof(Para));
  //CRTS_dma_iget(Para_s, Paras, sizeof(Para), &dma_rply);
  //D_COUNT++;
  //CRTS_dma_wait_value(&dma_rply, D_COUNT);

  double *d_h0bl_s = Para_s->d_h0bl;
  double *d_h0bf_s = Para_s->d_h0bf;
  //double *d_wp12_s = Para_s->d_wp12;
  //double *d_wp13_s = Para_s->d_wp13;

  double *d_viv,*d_vivjm1,*d_vit,*d_vp,*d_up,*d_upjm1,*d_vpjm1;
  double d_wp12_s[imt],d_wp13_s[imt];

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
         k = blockId_z * Para_s->cuBlock[2];

        d_viv    = &Para_s->d_viv[k * imt * jmt + j * imt];
        d_vivjm1 = &Para_s->d_viv[k * imt * jmt + (j-1) * imt];
        d_vit    = &Para_s->d_vit[k * imt * jmt + j * imt];
        d_vp     = &Para_s->d_vp[k * imt * jmt + j * imt];
        d_up     = &Para_s->d_up[k * imt * jmt + j * imt];
        d_upjm1  = &Para_s->d_up[k * imt * jmt + (j-1) * imt];
        d_vpjm1  = &Para_s->d_vp[k * imt * jmt + (j-1) * imt];

        for(i=0;i<imt;i++) {

            if(jst-1<=j&&j<jet&&k==0)
            {
                d_h0bl_s[j *imt+i]=d_h0bf_s[j *imt+i];
                d_h0bf_s[j *imt+i]=Para_s->d_h0[j *imt+i];
            }
            if(1<=j&&i<(imt-1))
            {
                d_wp12_s[i]=d_vit[i] *(Para_s->d_at0[j *imt+i] *d_up[i] *d_viv[i]+Para_s->d_atn[j *imt+i] *d_upjm1[i] *d_vivjm1[i]+Para_s->d_ate[j *imt+i] *d_up[i+1] *d_viv[i+1]+Para_s->d_atne[j *imt+i] *d_upjm1[i+1] *d_vivjm1[i+1]) /(d_viv[i] *Para_s->d_at0[j *imt+i]+d_vivjm1[i] *Para_s->d_atn[j *imt+i]+d_viv[i+1] *Para_s->d_ate[j *imt+i]+d_vivjm1[i+1] *Para_s->d_atne[j *imt+i]+1.0e-25);
      	        d_wp13_s[i]=d_vit[i] *(Para_s->d_at0[j *imt+i] *d_vp[i] *d_viv[i]+Para_s->d_atn[j *imt+i] *d_vpjm1[i] *d_vivjm1[i]+Para_s->d_ate[j *imt+i] *d_vp[i+1] *d_viv[i+1]+Para_s->d_atne[j *imt+i] *d_vpjm1[i+1] *d_vivjm1[i+1]) /(d_viv[i] *Para_s->d_at0[j *imt+i]+d_vivjm1[i] *Para_s->d_atn[j *imt+i]+d_viv[i+1] *Para_s->d_ate[j *imt+i]+d_vivjm1[i+1] *Para_s->d_atne[j *imt+i]+1.0e-25);
            }

        }
        CRTS_dma_put(&Para_s->d_wp12[k * imt * jmt + j * imt], d_wp12_s, sizeof(double) * imt);
        CRTS_dma_put(&Para_s->d_wp13[k * imt * jmt + j * imt], d_wp13_s, sizeof(double) * imt);

  }


    ldm_free(Para_s, sizeof(Para));
}


void  readyc1_2_(readyc1_2_Para *Paras){	

    readyc1_2_Para* Para_s;

    //ATHREAD basic operation
    int blockId_x,blockId_y,blockId_z;
    int i,j,k;
    int blockIndex = CRTS_tid;
    Para_s = (readyc1_2_Para*)ldm_malloc(sizeof(readyc1_2_Para));
    CRTS_dma_get(Para_s, Paras, sizeof(readyc1_2_Para));
    //D_COUNT++;
    //CRTS_dma_wait_value(&dma_rply, D_COUNT);


    double *d_wp12_s = Para_s->d_wp12;
    double *d_wp13_s = Para_s->d_wp13;

    double riv1 = 0.0;
    double riv2 = 0.0;

    double epsln=1.0e-25;
    double *d_wp12kp1,*d_vitkp1,*d_wp13kp1,*d_vivkp1,*d_upkp1,*d_vpkp1,*d_vit,*d_rict,*d_ric;
    double *d_up,*d_vp;
    double d_s2t[imt],d_s2u[imt],d_rit[imt],d_riu[imt];

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
        k = blockId_z * Para_s->cuBlock[2];

       d_wp12kp1 = &Para_s->d_wp12[(k+1) * imt * jmt + j * imt];
       d_wp13kp1 = &Para_s->d_wp13[(k+1) * imt * jmt + j * imt];
       d_vitkp1  = &Para_s->d_vit[(k+1) * imt * jmt + j * imt];
       d_vivkp1  = &Para_s->d_viv[(k+1) * imt * jmt + j * imt];
       d_upkp1   = &Para_s->d_up[(k+1) * imt * jmt + j * imt];
       d_vpkp1   = &Para_s->d_vp[(k+1) * imt * jmt + j * imt];
       d_up      = &Para_s->d_up[k * imt * jmt + j * imt];
       d_vp      = &Para_s->d_vp[k * imt * jmt + j * imt];
       d_ric     = &Para_s->d_ric[k * imt * jmt + j * imt];
       d_rict    = &Para_s->d_rict[k * imt * jmt + j * imt];

       for( i=0;i<imt;i++) {

            if(1<=j&&i<imt-1)
            {
              riv1 = d_wp12_s[k *jmt *imt+j *imt+i] *
                      Para_s->d_vit[k *jmt *imt+j *imt+i]-d_wp12kp1[i] *d_vitkp1[i];

              riv2 = d_wp13_s[k *jmt *imt+j *imt+i] *
                      Para_s->d_vit[k *jmt *imt+j *imt+i]-d_wp13kp1[i] *d_vitkp1[i];

              d_s2t[i]=d_vitkp1[i] *
                      (riv1 * riv1 +
                       riv2 * riv2) *Para_s->d_odzt[k+1] * Para_s->d_odzt[k+1];
           // #ifdef CANUTO
              d_rit[i]=d_vitkp1[i] *
                      d_rict[i] /(d_s2t[i]+epsln);
            //#else
            //printf("need logic to handle\n");
            //#endif
            }
            riv1 = d_up[i]-d_upkp1[i];
            riv2 = d_vp[i]-d_vpkp1[i];
            d_s2u[i]=d_vivkp1[i] *
                    (riv1 * riv1 +
                     riv2 * riv2) *Para_s->d_odzt[k+1] *Para_s->d_odzt[k+1];
            d_riu[i]=d_vivkp1[i] *d_ric[i] /(d_s2u[i]+epsln);
        }
       CRTS_dma_put(&Para_s->d_riu[(k+1) * imt * jmt + j * imt],d_riu, sizeof(double) * imt);
       CRTS_dma_put(&Para_s->d_s2t[k * imt * jmt + j * imt],d_s2t, sizeof(double) * imt);
       CRTS_dma_put(&Para_s->d_s2u[k * imt * jmt + j * imt],d_s2u, sizeof(double) * imt);
//#ifdef CANUTO
       CRTS_dma_put(&Para_s->d_rit[k * imt * jmt + j * imt],d_rit, sizeof(double) * imt);
//#endif


    }

    ldm_free(Para_s, sizeof(readyc1_2_Para));
}

 void  readyc1_3_(readyc1_3_Para *Paras){	

  //ATHREAD local variable definition
  int swGrid[3];
  int swBlock[3];
  int stride;
  int blockId_x,blockId_y,blockId_z;
  int threadId_x,threadId_y,threadId_z;
  readyc1_3_Para* Para_s;


  //ATHREAD basic operation
     int i,j,k;
  int blockIndex = CRTS_tid;
  Para_s = (readyc1_3_Para*)ldm_malloc(sizeof(readyc1_3_Para));
  CRTS_dma_get(Para_s, Paras, sizeof(readyc1_3_Para));
  

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

  double *d_diff_back = Para_s->d_diff_back;
  double *d_diff_back_sh = Para_s->d_diff_back_sh;
  double *d_diff_back_nh = Para_s->d_diff_back_nh;
  double *d_tlat = Para_s->d_tlat;

  for(int block = blockStart; block < blockEnd; block++){
    blockId_z = block / (Para_s->cuGrid[0] * Para_s->cuGrid[1]);
    blockId_y = (block - blockId_z * (Para_s->cuGrid[0] * Para_s->cuGrid[1])) / Para_s->cuGrid[0];
    blockId_x = block % Para_s->cuGrid[0];
    j = blockId_y * Para_s->cuBlock[1];


    //ATHREAD serial operation replace GPU parallel mutli-thread in single thread-block
    for(i = 0; i < threadSum; i++) {

      if(j<jmt&&i<imt)	
      {	
      	d_diff_back[j *imt+i]=0.0e0;	
      	d_diff_back_sh[j *imt+i]=0.0e0;	
      	d_diff_back_nh[j *imt+i]=0.0e0;	
      if(d_tlat[j *imt+i]<0.0)	
      {	
      	d_diff_back_sh[j *imt+i]=Para_s->d_diff_back_coef_max *exp(-(0.4e0 *(d_tlat[j *imt+i] /DEGtoRAD+28.9)) *(0.4e0 *(d_tlat[j *imt+i] /DEGtoRAD+28.9)));	
      }	
      else	{
	      	d_diff_back_nh[j *imt+i]=Para_s->d_diff_back_coef_max *exp(-(0.4e0 *(d_tlat[j *imt+i] /DEGtoRAD-28.9)) *(0.4e0 *(d_tlat[j *imt+i] /DEGtoRAD-28.9)));	
      }	
      	d_diff_back[j *imt+i]=Para_s->d_diff_back_eq+d_diff_back_sh[j *imt+i]+d_diff_back_nh[j *imt+i];	
      if(d_tlat[j *imt+i]<-10.0 *DEGtoRAD)	
      {	
      	d_diff_back[j *imt+i]=d_diff_back[j *imt+i]+Para_s->d_diff_back_coef;	
      }	
      else	
      if(fabs(d_tlat[j *imt+i])<=10.0 *DEGtoRAD)	
      {	
      	d_diff_back[j *imt+i]=d_diff_back[j *imt+i]+Para_s->d_diff_back_coef *((fabs(d_tlat[j *imt+i]) /DEGtoRAD /10.0) *(fabs(d_tlat[j *imt+i]) /DEGtoRAD /10.0));	
      }	
      else	{
	      	d_diff_back[j *imt+i]=d_diff_back[j *imt+i]+Para_s->d_diff_back_coef;	
      }	
	
      }	

    }

  }

    ldm_free(Para_s, sizeof(readyc1_3_Para));
}

void  readyc2_(readyc2_Para *Paras){

  //ATHREAD basic operation
    int blockId_x,blockId_y,blockId_z;
    int i,j,k;
  int blockIndex = CRTS_tid;
  readyc2_Para* Para_s = (readyc2_Para*)ldm_malloc(sizeof(readyc2_Para));
  CRTS_dma_get(Para_s, Paras, sizeof(readyc2_Para));

  double *d_aus = Para_s->d_aus;
  double *d_auw = Para_s->d_auw;
  double *d_au0 = Para_s->d_au0;
  double *d_ausw = Para_s->d_ausw;

  double  d_akmu[imt];
  double  *d_akmt, *d_akmtjp1, *d_vivkp1;
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
        k = blockId_z * Para_s->cuBlock[2];

        if(j<jmt-1) {
            d_akmt    = &Para_s->d_akmt[k * imt * jmt + j * imt];
            d_vivkp1  = &Para_s->d_viv[(k + 1) * imt * jmt + j * imt];
            d_akmtjp1 = &Para_s->d_akmt[k * imt * jmt + (j + 1) * imt];

            for (i = 1; i < imt; i++) {
                d_akmu[i] = (p25 * d_akmt[i] + p25 * d_akmtjp1[i] +
                        p25 * d_akmt[i - 1] + p25 * d_akmtjp1[i - 1]) * d_vivkp1[i];
            }
            CRTS_dma_put(&Para_s->d_akmu[k * imt * jmt + j * imt], d_akmu, sizeof(double) * imt);
        }

    }

    ldm_free(Para_s, sizeof(readyc2_Para));
}


void  readyc3_(readyc3_Para *Paras){

  //ATHREAD basic operation
    int blockId_x,blockId_y,blockId_z;
    int i,j,k;
  int blockIndex = CRTS_tid;
  readyc3_Para* Para_s = (readyc3_Para*)ldm_malloc(sizeof(readyc3_Para));
  CRTS_dma_get(Para_s, Paras, sizeof(readyc3_Para));

  double d_wka[imt],d_u_wface[imt],d_v_sface[imt];
  double *d_ws, *d_u;
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
        k = blockId_z * Para_s->cuBlock[2];

        if(j<jmt-1) {
            d_ws = &Para_s->d_ws[k * imt * jmt + j * imt];
            d_u  = &Para_s->d_u[k * imt * jmt + j * imt];
            for( i=1;i<imt;i++) {
                unsigned int tid=j *imt+i;
                unsigned int ttid=k *jmt *imt+j *imt+i;

                d_wka[i]=p25 *(d_ws[i])+p25 *(Para_s->d_ws[ttid+imt])+p25 *(d_ws[i-1])+p25 *(Para_s->d_ws[ttid+imt-1]);
                d_u_wface[i]=(d_u[i-1]+d_u[i]) *p25 *Para_s->d_hue[tid-1];
                d_v_sface[i]=(Para_s->d_v[ttid]+Para_s->d_v[ttid+imt]) *p25 *Para_s->d_hun[(j+1) *imt+i];
            }
          CRTS_dma_put(&Para_s->d_wka[k * imt * jmt + j * imt],d_wka, sizeof(double) * imt);
          CRTS_dma_put(&Para_s->d_u_wface[k * imt * jmt + j * imt],d_u_wface, sizeof(double) * imt);
          CRTS_dma_put(&Para_s->d_v_sface[k * imt * jmt + j * imt],d_v_sface, sizeof(double) * imt);

        }

    }

    ldm_free(Para_s, sizeof(readyc3_Para));
}

void  advection_momentum_cu_(advection_momentum_Para *Paras){

  //ATHREAD basic operation
    int blockId_x,blockId_y,blockId_z;
    int i,j,k;
  int blockIndex = CRTS_tid;
  advection_momentum_Para* Para_s = (advection_momentum_Para*)ldm_malloc(sizeof(advection_momentum_Para));
  CRTS_dma_get(Para_s, Paras, sizeof(advection_momentum_Para));

  double d_dlu[imt],d_dlv[imt];
  double *d_u_wface, *d_u, *d_v;
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
        k = blockId_z * Para_s->cuBlock[2];

        if (j>1 && j<jmt-2) {
            d_u_wface = &Para_s->d_u_wface[k * imt * jmt + j * imt];
            d_u       = &Para_s->d_u[k * imt * jmt + j * imt];
            d_v       = &Para_s->d_v[k * imt * jmt + j * imt];
            for(i=2;i<imt-2;i++) {
                double adv_z1=0.0e0,adv_z2=0.0e0,adv_z3=0.0e0,adv_z4=0.0e0;
                unsigned int tid=j *imt+i;
                unsigned int ttid=k *jmt *imt+j *imt+i;

                d_dlu[i]=(-d_u_wface[i] *(d_u[i]-d_u[i-1])-d_u_wface[i+1] *(d_u[i+1]-d_u[i])-Para_s->d_v_sface[ttid] *(Para_s->d_u[k *jmt *imt+(j+1) *imt+i]-d_u[i])-Para_s->d_v_sface[ttid-imt] *(d_u[i]-Para_s->d_u[ttid-imt])) *Para_s->d_uarea_r[tid];
                d_dlv[i]=(-d_u_wface[i] *(d_v[i]-d_v[i-1])-d_u_wface[i+1] *(d_v[i+1]-d_v[i])-Para_s->d_v_sface[ttid] *(Para_s->d_v[k *jmt *imt+(j+1) *imt+i]-d_v[i])-Para_s->d_v_sface[ttid-imt] *(d_v[i]-Para_s->d_v[ttid-imt])) *Para_s->d_uarea_r[tid];
      //if(k==0)

//                      Para_s->d_utest[tid]=Para_s->d_u[ttid]-Para_s->d_u[ttid-1];

              if(k> 0)
              {
                adv_z1=Para_s->d_wka[ttid] *(Para_s->d_u[(k-1) *jmt *imt+tid]-d_u[i]);
                adv_z3=Para_s->d_wka[ttid] *(Para_s->d_v[(k-1) *jmt *imt+tid]-d_v[i]);
              }
              if(k<km-1)
              {
                adv_z2=Para_s->d_wka[jmt *imt+ttid] *(d_u[i]-Para_s->d_u[jmt *imt+ttid]);
                adv_z4=Para_s->d_wka[jmt *imt+ttid] *(d_v[i]-Para_s->d_v[jmt *imt+ttid]);
              }
                d_dlu[i]+=-p5 *Para_s->d_odzp[k] *(adv_z1+adv_z2);
                d_dlv[i]+=-p5 *Para_s->d_odzp[k] *(adv_z3+adv_z4);
            }
            CRTS_dma_put(&Para_s->d_dlu[k * imt * jmt + j * imt],d_dlu, sizeof(double) * imt);
            CRTS_dma_put(&Para_s->d_dlv[k * imt * jmt + j * imt],d_dlv, sizeof(double) * imt);

        }
    }

    ldm_free(Para_s, sizeof(advection_momentum_Para));
}


 void  readyc4_(readyc4_Para *Paras){

  //ATHREAD basic operation
    int blockId_x,blockId_y,blockId_z;
    int i,j,k;
    int blockIndex = CRTS_tid;
    readyc4_Para* Para_s = (readyc4_Para*)ldm_malloc(sizeof(readyc4_Para));
    CRTS_dma_get(Para_s, Paras, sizeof(readyc4_Para));
    int kmb;

    double d_sbcx[imt],d_sbcy[imt],d_bbcx[imt],d_bbcy[imt],d_dlub[imt],d_dlvb[imt];

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

         if(j>0 && j<jmt-1) {
             CRTS_dma_get(d_dlub, &Para_s->d_dlub[j * imt], sizeof(double) * imt);
             CRTS_dma_get(d_dlvb, &Para_s->d_dlvb[j * imt], sizeof(double) * imt);

             for (i = 1; i < imt - 1; i++) {

                 kmb = Para_s->d_kmu[j * imt + i] - 1;
                 if (kmb >= 0) {
                     d_sbcx[i] = Para_s->d_su[j * imt + i] * Para_s->d_od0;
                     d_sbcy[i] = Para_s->d_sv[j * imt + i] * Para_s->d_od0;
                     d_bbcx[i] = Para_s->d_c0f * sqrt(Para_s->d_up[kmb * jmt * imt + j * imt + i] *
                                                      Para_s->d_up[kmb * jmt * imt + j * imt + i] +
                                                      Para_s->d_vp[kmb * jmt * imt + j * imt + i] *
                                                      Para_s->d_vp[kmb * jmt * imt + j * imt + i]) *
                                 (Para_s->d_up[kmb * jmt * imt + j * imt + i] * Para_s->d_cag +
                                  Para_s->d_snlat[j * imt + i] * Para_s->d_vp[kmb * jmt * imt + j * imt + i] *
                                  Para_s->d_sag);
                     d_bbcy[i] = Para_s->d_c0f * sqrt(Para_s->d_up[kmb * jmt * imt + j * imt + i] *
                                                      Para_s->d_up[kmb * jmt * imt + j * imt + i] +
                                                      Para_s->d_vp[kmb * jmt * imt + j * imt + i] *
                                                      Para_s->d_vp[kmb * jmt * imt + j * imt + i]) *
                                 (-Para_s->d_snlat[j * imt + i] * Para_s->d_up[kmb * jmt * imt + j * imt + i] *
                                  Para_s->d_sag + Para_s->d_vp[kmb * jmt * imt + j * imt + i] * Para_s->d_cag);
                 } else {
                     d_sbcx[i] = 0.0e0;
                     d_sbcy[i] = 0.0e0;
                     d_bbcx[i] = 0.0e0;
                     d_bbcy[i] = 0.0e0;
                 }
                 d_dlub[i] = d_dlub[i] + (d_sbcx[i] - d_bbcx[i]) * Para_s->d_ohbu[j * imt + i];
                 d_dlvb[i] = d_dlvb[i] + (d_sbcy[i] - d_bbcy[i]) * Para_s->d_ohbu[j * imt + i];

             }
             CRTS_dma_put(&Para_s->d_dlub[j * imt], d_dlub, sizeof(double) * imt);
             CRTS_dma_put(&Para_s->d_dlvb[j * imt], d_dlvb, sizeof(double) * imt);
         }
    }

    ldm_free(Para_s, sizeof(readyc4_Para));
}

void  readyc4_2_(readyc4_2_Para *Paras){

    //ATHREAD basic operation
    int blockId_x,blockId_y,blockId_z;
    int i,j,k;
    int blockIndex = CRTS_tid;
    readyc4_2_Para* Para_s = (readyc4_2_Para*)ldm_malloc(sizeof(readyc4_2_Para));
    CRTS_dma_get(Para_s, Paras, sizeof(readyc4_2_Para));

    double *d_akmukm1, *d_wkakm1, *d_vpkm1, *d_vp, *d_wka, *d_vivkp1, *d_viv, *d_upkm1, *d_up, *d_upkp1, *d_akmu;

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
        k = blockId_z * Para_s->cuBlock[2];
        if(j>0 && j < jmt-1) {
            d_akmukm1 = &Para_s->d_akmu[(k - 1) * imt * jmt + j * imt];
            d_wkakm1  = &Para_s->d_wka[(k - 1) * imt * jmt + j * imt];
            d_vpkm1   = &Para_s->d_vp[(k - 1) * imt * jmt + j * imt];
            d_vp      = &Para_s->d_vp[k * imt * jmt + j * imt];
            d_wka     = &Para_s->d_wka[k * imt * jmt + j * imt];
            d_vivkp1  = &Para_s->d_viv[(k + 1) * imt * jmt + j * imt];
            d_viv     = &Para_s->d_viv[k * imt * jmt + j * imt];
            d_upkm1   = &Para_s->d_up[(k - 1) * imt * jmt + j * imt];
            d_up      = &Para_s->d_up[k * imt * jmt + j * imt];
            d_upkp1   = &Para_s->d_up[(k + 1) * imt * jmt + j * imt];
            d_akmu    = &Para_s->d_akmu[k * imt * jmt + j * imt];

            for (i = 1; i < imt - 1; i++) {
                double diff_u1, diff_u2, diff_v1, diff_v2;
                double aidif = 0.5;

                if (k == 0) {
                    diff_v1 = Para_s->d_sv[j * imt + i] * Para_s->d_od0 * (1 - aidif);
                    diff_u1 = Para_s->d_su[j * imt + i] * Para_s->d_od0 * (1 - aidif);
                } else {
                    diff_v1 = d_akmukm1[i] * (1 - aidif) * (d_vpkm1[i] - d_vp[i]) * Para_s->d_odzt[k] *
                            d_viv[i] + (1.0e0 - d_viv[i]) * d_wkakm1[i] * (1 - aidif) *
                            (-Para_s->d_snlat[j * imt + i] * d_upkm1[i] * Para_s->d_sag +
                            d_vpkm1[i] * Para_s->d_cag);
                    diff_u1 = d_akmukm1[i] * (1 - aidif) * (d_upkm1[i] - d_up[i]) * Para_s->d_odzt[k] *
                            d_viv[i] + (1.0e0 - d_viv[i]) * d_wkakm1[i] * (1 - aidif) *
                              (d_upkm1[i] * Para_s->d_cag + Para_s->d_snlat[j * imt + i] *
                              d_vpkm1[i] * Para_s->d_sag);
                }
                if (k == km - 1) {
                    diff_v2 = d_wka[i] * (-Para_s->d_snlat[j * imt + i] * d_up[i] * Para_s->d_sag +
                            d_vp[i] * Para_s->d_cag) *(1 - aidif);
                    diff_u2 = d_wka[i] * (d_up[i] * Para_s->d_cag + Para_s->d_snlat[j * imt + i] *
                            d_vp[i] * Para_s->d_sag) * (1 - aidif);
                } else {
                    diff_v2 = d_akmu[i] * (1 - aidif) * (d_vp[i] - Para_s->d_vp[(k+1)*jmt*imt+j*imt+i]) *
                              Para_s->d_odzt[k + 1] * d_vivkp1[i] + (1.0e0 - d_vivkp1[i]) *
                              d_wka[i] * (1 - aidif) * (-Para_s->d_snlat[j * imt + i] *
                              d_up[i] * Para_s->d_sag + d_vp[i] * Para_s->d_cag);
                    diff_u2 = d_akmu[i] * (1 - aidif) * (d_up[i] - d_upkp1[i]) *
                            Para_s->d_odzt[k + 1] * d_vivkp1[i] +
                            (1.0e0 - d_vivkp1[i]) * d_wka[i] * (1 - aidif) *
                              (d_up[i] * Para_s->d_cag + Para_s->d_snlat[j * imt + i] *
                              d_vp[i] * Para_s->d_sag);
                }
                Para_s->d_dlv[k * jmt * imt + j * imt + i] =
                        Para_s->d_dlv[k * jmt * imt + j * imt + i] + Para_s->d_odzp[k] * (diff_v1 - diff_v2);
                Para_s->d_dlu[k * jmt * imt + j * imt + i] =
                        Para_s->d_dlu[k * jmt * imt + j * imt + i] + Para_s->d_odzp[k] * (diff_u1 - diff_u2);
            }
        }

    }

    ldm_free(Para_s, sizeof(readyc4_2_Para));
}

void  readyc_vinteg_cu_(vinteg_Para *Paras){

    double d_dlub[imt];
    double d_dlvb[imt];

    //ATHREAD basic operation
    int blockIndex = CRTS_tid;
    vinteg_Para* Para_s = (vinteg_Para*)ldm_malloc(sizeof(vinteg_Para));
    CRTS_dma_get(Para_s, Paras, sizeof(vinteg_Para));
    double *d_dlu, *d_dlv, *d_viv;

    for(int k = 0; k < km; k++) {
        for (int j = blockIndex; j < jmt; j += 64) {

            d_dlu = &Para_s->d_dlu[k * imt * jmt + j * imt];
            d_dlv = &Para_s->d_dlv[k * imt * jmt + j * imt];
	          d_viv = &Para_s->d_viv[k * imt * jmt + j * imt];
	          CRTS_dma_get(d_dlub, &Para_s->d_dlub[j * imt], sizeof(double) * imt);
	          CRTS_dma_get(d_dlvb, &Para_s->d_dlvb[j * imt], sizeof(double) * imt);

            for (int i = 0; i < imt; i++) {
		            d_dlub[i]+=Para_s->d_dzp[k] *Para_s->d_ohbu[j *imt+i] *d_dlu[i] *d_viv[i];
		            d_dlvb[i]+=Para_s->d_dzp[k] *Para_s->d_ohbu[j *imt+i] *d_dlv[i] *d_viv[i];
            }

            CRTS_dma_put(&Para_s->d_dlub[j*imt], d_dlub, sizeof(double) * imt);
            CRTS_dma_put(&Para_s->d_dlvb[j*imt], d_dlvb, sizeof(double) * imt);
        }
    }

    ldm_free(Para_s, sizeof(vinteg_Para));
}

void hdiffu_del2_(hdiffu_del2_Para *Paras){ 
  int blockIndex = CRTS_tid;
  hdiffu_del2_Para *Para_s = (hdiffu_del2_Para*)ldm_malloc(sizeof(hdiffu_del2_Para));
  CRTS_dma_get(Para_s, Paras, sizeof(hdiffu_del2_Para));

  double hduk = 0.0, hdvk = 0.0;

  for(int k=0;k<km;k++)
  {
     for(int j=blockIndex;j<jmt;j+=64)
     {
	for(int i=0;i<imt;i++)
	{
	   hduk = 0.0; 
           hdvk = 0.0; 

           if (j >= Para_s->d_jb - 1 && j < Para_s->d_je && i >= Para_s->d_ib - 1 && i < Para_s->d_ie) {
                 hduk = Para_s->d_am * (((Para_s->d_duc[j * imt + i] + 
                                          Para_s->d_dum[j * imt + i]) * Para_s->d_up[k * jmt * imt + j * imt + i] + 
                                          Para_s->d_dun[j * imt + i] * Para_s->d_up[k * jmt * imt + (j - 1) * imt + i] + 
                                          Para_s->d_dus[j * imt + i] * Para_s->d_up[k * jmt * imt + (j + 1) * imt + i] + 
                                          Para_s->d_due[j * imt + i] * Para_s->d_up[k * jmt * imt + j * imt + i + 1] + 
                                          Para_s->d_duw[j * imt + i] * Para_s->d_up[k * jmt * imt + j * imt + i - 1]) +
                                          (Para_s->d_dmc[j * imt + i] * Para_s->d_vp[k * jmt * imt + j * imt + i] + 
                                           Para_s->d_dmn[j * imt + i] * Para_s->d_vp[k * jmt * imt + (j - 1) * imt + i] + 
                                           Para_s->d_dms[j * imt + i] * Para_s->d_vp[k * jmt * imt + (j + 1) * imt + i] + 
                                           Para_s->d_dme[j * imt + i] * Para_s->d_vp[k * jmt * imt + j * imt + i + 1] + 
                                           Para_s->d_dmw[j * imt + i] * Para_s->d_vp[k * jmt * imt + j * imt + i - 1])) *
                                         Para_s->d_viv[k * jmt * imt + j * imt + i];

                  hdvk = Para_s->d_am * (((Para_s->d_duc[j * imt + i] + 
                                           Para_s->d_dum[j * imt + i]) * Para_s->d_vp[k * jmt * imt + j * imt + i] + 
                                           Para_s->d_dun[j * imt + i] * Para_s->d_vp[k * jmt * imt + (j - 1) * imt + i] + 
                                           Para_s->d_dus[j * imt + i] * Para_s->d_vp[k * jmt * imt + (j + 1) * imt + i] + 
                                           Para_s->d_due[j * imt + i] * Para_s->d_vp[k * jmt * imt + j * imt + i + 1] + 
                                           Para_s->d_duw[j * imt + i] * Para_s->d_vp[k * jmt * imt + j * imt + i - 1]) -
                                          (Para_s->d_dmc[j * imt + i] * Para_s->d_up[k * jmt * imt + j * imt + i] + 
                                           Para_s->d_dmn[j * imt + i] * Para_s->d_up[k * jmt * imt + (j - 1) * imt + i] + 
                                           Para_s->d_dms[j * imt + i] * Para_s->d_up[k * jmt * imt + (j + 1) * imt + i] + 
                                           Para_s->d_dme[j * imt + i] * Para_s->d_up[k * jmt * imt + j * imt + i + 1] + 
                                           Para_s->d_dmw[j * imt + i] * Para_s->d_up[k * jmt * imt + j * imt + i - 1])) *
                                         Para_s->d_viv[k * jmt * imt + j * imt + i];
                }    

                if ((k > Para_s->d_kmu[j * imt + i] - 1)) {
                        hduk = 0.0; 
                        hdvk = 0.0; 
                }    

                if (2 <= j && j < jmt - 2 && 2 <= i && i < imt - 2) { 
                        Para_s->d_dlv[k * jmt * imt + j * imt + i] = 
                                Para_s->d_dlv[k * jmt * imt + j * imt + i] + 
                                hdvk;

                        Para_s->d_dlu[k * jmt * imt + j * imt + i] = 
                                Para_s->d_dlu[k * jmt * imt + j * imt + i] + 
                                hduk;
                }    

                Para_s->d_wka[k * jmt * imt + j * imt + i] = 
                        Para_s->d_c0f * sqrt(Para_s->d_up[k * jmt * imt + j * imt + i] * 
                                        Para_s->d_up[k * jmt * imt + j * imt + i] + 
                                        Para_s->d_vp[k * jmt * imt + j * imt + i] * 
                                        Para_s->d_vp[k * jmt * imt + j * imt + i]); 

	}
     }
  }
  ldm_free(Para_s, sizeof(hdiffu_del2_Para));
}
