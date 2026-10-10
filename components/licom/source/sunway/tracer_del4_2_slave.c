
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

//ATHREAD DMA transfer definition 
__thread_local unsigned int D_COUNT = 0;
__thread_local crts_rply_t dma_rply = 0;
__thread_local crts_rply_t l_rply = 0;
__thread_local crts_rply_t r_rply = 0;

//__userDef__ paraVarAttr(int,n,sizeof(int),readonly,no,no,no);	
//__userDef__ paraVarAttr(double,d_tf,sizeof(double)*km*jmt*imt,outonly,yes,no,no);	
//__userDef__ paraVarAttr(double,d_adv_tt,sizeof(double)*km*jmt*imt,readonly,no,no,no);	
//__userDef__ paraVarAttr(double,d_vit,sizeof(double)*km*jmt*imt,readonly,no,no,no);	
//__userDef__ paraVarAttr(double,d_akt,sizeof(double)*ntra*km*jmt*imt,readonly,no,no,no);	
//__userDef__ paraVarAttr(double,d_wkc,sizeof(double)*km*jmt*imt,outonly,yes,no,no);	
//__userDef__ paraVarAttr(double,d_dwndmix,sizeof(double),readonly,no,no,no);	
//__userDef__ paraVarAttr(double,d_ah,sizeof(double),readonly,no,no,no);	
//__userDef__ paraVarAttr(int,d_kmtn,sizeof(int)*ny_block*nx_block,readonly,no,no,no);	
//__userDef__ paraVarAttr(int,d_kmts,sizeof(int)*ny_block*nx_block,readonly,no,no,no);	
//__userDef__ paraVarAttr(int,d_kmte,sizeof(int)*ny_block*nx_block,readonly,no,no,no);	
//__userDef__ paraVarAttr(int,d_kmtw,sizeof(int)*ny_block*nx_block,readonly,no,no,no);	
//__userDef__ paraVarAttr(double,d_dtn,sizeof(double)*ny_block*nx_block,readonly,no,no,no);	
//__userDef__ paraVarAttr(double,d_dts,sizeof(double)*ny_block*nx_block,readonly,no,no,no);	
//__userDef__ paraVarAttr(double,d_dte,sizeof(double)*ny_block*nx_block,readonly,no,no,no);	
//__userDef__ paraVarAttr(double,d_dtw,sizeof(double)*ny_block*nx_block,readonly,no,no,no);	
//__userDef__ paraVarAttr(int,d_kmt,sizeof(int)*ny_block*nx_block,readonly,no,no,no);	
//__userDef__ paraVarAttr(double,d_d2tk,sizeof(double)*km*jmt*imt,outonly,yes,no,no);	
typedef struct tracer_del4_2_kernelPara{
	int cuBlock[3];
	int cuGrid[3];
	int n;
	int *d_kmtn;
	int *d_kmts;
	int *d_kmte;
	int *d_kmtw;
	int *d_kmt;
	double d_ah;
	double d_dwndmix;
	double *d_tf;
  //double *d_dx;
	double *d_adv_tt;
	double *d_vit;
	double *d_akt;
	double *d_wkc;
	double *d_dtn;
	double *d_dts;
	double *d_dte;
	double *d_dtw;
	double *d_d2tk;
    double *d_atb;
    double *d_factor_vis_tsz0;
#if (defined ISO)
    double *d_ahisop;
    double *d_k3;
#endif
}tracer_del4_2_Para;


//dynamic atomic definition
#define CP_CORE      64
#define RCV_CORE     0
#define DATAPAKAGE   1

//ATHREAD global definition 
__thread_local int blockIndex;

void tracer_del4_2_(tracer_del4_2_Para *Paras){	

  //ATHREAD local variable definition
  int swGrid[3];
  int swBlock[3];
  int stride;
  int blockId_x,blockId_y,blockId_z;
  int threadId_x,threadId_y,threadId_z;


  //ATHREAD basic operation
    int i,j,k;
  blockIndex = CRTS_tid;
  tracer_del4_2_Para *Para_s = (tracer_del4_2_Para*)ldm_malloc(sizeof(tracer_del4_2_Para));
  CRTS_dma_iget(Para_s, Paras, sizeof(tracer_del4_2_Para), &dma_rply);
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

  double *d_tf_s = Para_s->d_tf;
  double *d_wkc_s = Para_s->d_wkc;
  double *d_d2tk_s = Para_s->d_d2tk;

  for(int block = blockStart; block < blockEnd; block++){
    blockId_z = block / (Para_s->cuGrid[0] * Para_s->cuGrid[1]);
    blockId_y = (block - blockId_z * (Para_s->cuGrid[0] * Para_s->cuGrid[1])) / Para_s->cuGrid[0];
    blockId_x = block % Para_s->cuGrid[0];
      j = blockId_y * Para_s->cuBlock[1];
      k = blockId_z * Para_s->cuBlock[2];

    //ATHREAD serial operation replace GPU parallel mutli-thread in single thread-block
    for(i = 0; i < threadSum; i++) {

      if(k>=0&&k<km)	
      {	
        if(j>=2&&j<jmt-2&&i>=2&&i<imt-2)	
        {	
      	  d_tf_s[k *jmt *imt+j *imt+i]=Para_s->d_adv_tt[k *jmt *imt+j *imt+i] *Para_s->d_vit[k *jmt *imt+j *imt+i];	
        }	
      }	
      if(k>=0&&k<km)	
      {	
        if(j>=1&&j<jmt-1&&i>=1&&i<imt-1)	
        {	
          if((Para_s->d_akt[Para_s->n *km *jmt *imt+j *imt+i]<Para_s->d_dwndmix)&&(k==0))
          {
      	      Para_s->d_akt[Para_s->n *km *jmt *imt+j *imt+i]=Para_s->d_dwndmix;
              d_wkc_s[k *jmt *imt+j *imt+i] = Para_s->d_dwndmix;
          }
          else
      	      d_wkc_s[k *jmt *imt+j *imt+i]=Para_s->d_akt[Para_s->n *km *jmt *imt+k *jmt *imt+j *imt+i];
          double hdtk=0.0e0;	
	
          double cc=0.0e0;
          double cn=0.0e0;
          double cs=0.0e0;
          double ce=0.0e0;
          double cw=0.0e0;

          if(k<=Para_s->d_kmt[j *imt+i]-1)	
          {	
            double c;
	
            if(k<=Para_s->d_kmtn[j *imt+i]-1)	
            {	
      	      cn=Para_s->d_dtn[j *imt+i];
			  //hdtk+=c *d_d2tk_s[k *jmt *imt+(j-1) *imt+i];
			  //cc-=c;
			}	
            if(k<=Para_s->d_kmts[j *imt+i]-1)	
		    {	
			  cs=Para_s->d_dts[j *imt+i];
			  //hdtk+=c *d_d2tk_s[k *jmt *imt+(j+1) *imt+i];
			  //cc-=c;
            }	
		    if(k<=Para_s->d_kmte[j *imt+i]-1)	
		    {	
			  ce=Para_s->d_dte[j *imt+i];
			  //hdtk+=c *d_d2tk_s[k *jmt *imt+j *imt+i+1];
			  //cc-=c;
			}	
			if(k<=Para_s->d_kmtw[j *imt+i]-1)	
            {	
			  cw=Para_s->d_dtw[j *imt+i];
			  //hdtk+=c *d_d2tk_s[k *jmt *imt+j *imt+i-1];
			  //cc-=c;
		    }
            cc = -(cn + cs + ce + cw);
		  }	
          if(j>=2&&j<jmt-2&&i>=2&&i<imt-2)	
          {
#ifdef  BIHAR
      	    hdtk=Para_s->d_ah * Para_s->d_factor_vis_tsz0[j *imt+i] *
				(cc *d_d2tk_s[k *jmt *imt+j *imt+i]+cn *d_d2tk_s[k *jmt *imt+(j-1) *imt+i] +
                    cs *d_d2tk_s[k *jmt *imt+(j+1) *imt+i] + ce *d_d2tk_s[k *jmt *imt+j *imt+i+1] +
                    cw *d_d2tk_s[k *jmt *imt+j *imt+i-1]);

          	if(Para_s->n == 0)
          		hdtk = hdtk * 1.0;
          	else
          		hdtk = hdtk * 1.5;
#else
              hdtk=Para_s->d_ah *(cc *Para_s->d_atb[Para_s->n *(km+1) *jmt *imt+(k+1) *jmt *imt+j *imt+i]+
                      cn *Para_s->d_atb[Para_s->n *(km+1) *jmt *imt+(k+1) *jmt *imt+(j-1) *imt+i] +
                      cs *Para_s->d_atb[Para_s->n *(km+1) *jmt *imt+(k+1) *jmt *imt+(j+1) *imt+i] +
                      ce *Para_s->d_atb[Para_s->n *(km+1) *jmt *imt+(k+1) *jmt *imt+j *imt+i+1] +
                      cw *Para_s->d_atb[Para_s->n *(km+1) *jmt *imt+(k+1) *jmt *imt+j *imt+i-1]);

          	if(Para_s->n == 0)
          		hdtk = hdtk * 0.05;
          	else
          		hdtk = hdtk * 0.05;
#endif
          }	
	  else{
	      	hdtk=0.0e0;
	  }


      	  d_tf_s[k *jmt *imt+j *imt+i]+=hdtk;
          //Para_s->d_dx[Para_s->n * jmt * imt + k * jmt * imt + j * imt + i] = hdtk;
        }	
      }	

    }

  }

    ldm_free(Para_s, sizeof(tracer_del4_2_Para));
}
