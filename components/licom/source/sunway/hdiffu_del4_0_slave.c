
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

//__userDef__ paraVarAttr(double,d_hduk2,sizeof(double)*km*jmt*imt,outonly,yes,no,no);	
//__userDef__ paraVarAttr(double,d_hdvk2,sizeof(double)*km*jmt*imt,outonly,yes,no,no);	
//__userDef__ paraVarAttr(int,d_jb,sizeof(int),readonly,no,no,no);	
//__userDef__ paraVarAttr(int,d_je,sizeof(int),readonly,no,no,no);	
//__userDef__ paraVarAttr(int,d_ib,sizeof(int),readonly,no,no,no);	
//__userDef__ paraVarAttr(int,d_ie,sizeof(int),readonly,no,no,no);	
//__userDef__ paraVarAttr(double,d_uarea_r,sizeof(double)*ny_block*nx_block,readonly,no,no,no);	
//__userDef__ paraVarAttr(double,d_up,sizeof(double)*km*jmt*imt,readonly,no,no,no);	
//__userDef__ paraVarAttr(double,d_vp,sizeof(double)*km*jmt*imt,readonly,no,no,no);	
//__userDef__ paraVarAttr(double,d_duc,sizeof(double)*ny_block*nx_block,readonly,no,no,no);	
//__userDef__ paraVarAttr(double,d_dum,sizeof(double)*ny_block*nx_block,readonly,no,no,no);	
//__userDef__ paraVarAttr(double,d_dun,sizeof(double)*ny_block*nx_block,readonly,no,no,no);	
//__userDef__ paraVarAttr(double,d_dus,sizeof(double)*ny_block*nx_block,readonly,no,no,no);	
//__userDef__ paraVarAttr(double,d_due,sizeof(double)*ny_block*nx_block,readonly,no,no,no);	
//__userDef__ paraVarAttr(double,d_duw,sizeof(double)*ny_block*nx_block,readonly,no,no,no);	
//__userDef__ paraVarAttr(double,d_dmc,sizeof(double)*ny_block*nx_block,readonly,no,no,no);	
//__userDef__ paraVarAttr(double,d_dmn,sizeof(double)*ny_block*nx_block,readonly,no,no,no);	
//__userDef__ paraVarAttr(double,d_dms,sizeof(double)*ny_block*nx_block,readonly,no,no,no);	
//__userDef__ paraVarAttr(double,d_dme,sizeof(double)*ny_block*nx_block,readonly,no,no,no);	
//__userDef__ paraVarAttr(double,d_dmw,sizeof(double)*ny_block*nx_block,readonly,no,no,no);	
//__userDef__ paraVarAttr(double,d_am,sizeof(double),readonly,no,no,no);	
//__userDef__ paraVarAttr(int,d_kmu,sizeof(int)*ny_block*nx_block,readonly,no,no,no);	
//__userDef__ paraVarAttr(double,d_amf,sizeof(double)*ny_block*nx_block,readonly,no,no,no);	
//__userDef__ paraVarAttr(double,d_am_factor,sizeof(double)*km*jmt*imt,readonly,no,no,no);	
//__userDef__ paraVarAttr(double,d_d2uk,sizeof(double)*km*jmt*imt,outonly,yes,no,no);	
//__userDef__ paraVarAttr(double,d_d2vk,sizeof(double)*km*jmt*imt,outonly,yes,no,no);	
typedef struct hdiffu_del4_0_kernelPara{
	int cuBlock[3];
	int cuGrid[3];
	int d_jb;
	int d_je;
	int d_ib;
	int d_ie;
	int rank;
	int *d_kmu;
	double d_am;
	double *d_uarea_r;
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
	double *d_amf;
	double *d_am_factor;
	double *d_d2uk;
	double *d_d2vk;
}hdiffu_del4_0_Para;


//dynamic atomic definition
#define CP_CORE      64
#define RCV_CORE     0
#define DATAPAKAGE   1

//ATHREAD global definition 
__thread_local int blockIndex;

void hdiffu_del4_0_(hdiffu_del4_0_Para *Paras){	

  //ATHREAD local variable definition
  int swGrid[3];
  int swBlock[3];
  int stride;
  int blockId_x,blockId_y,blockId_z;
  int threadId_x,threadId_y,threadId_z;


  //ATHREAD basic operation
    int i,j,k;
    blockIndex = CRTS_tid;
  hdiffu_del4_0_Para *Para_s = (hdiffu_del4_0_Para*)ldm_malloc(sizeof(hdiffu_del4_0_Para));
  CRTS_dma_iget(Para_s, Paras, sizeof(hdiffu_del4_0_Para), &dma_rply);
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

  double *d_d2uk_s = Para_s->d_d2uk;
  double *d_d2vk_s = Para_s->d_d2vk;

  for(int block = blockStart; block < blockEnd; block++){
    blockId_z = block / (Para_s->cuGrid[0] * Para_s->cuGrid[1]);
    blockId_y = (block - blockId_z * (Para_s->cuGrid[0] * Para_s->cuGrid[1])) / Para_s->cuGrid[0];
    blockId_x = block % Para_s->cuGrid[0];
      j = blockId_y * Para_s->cuBlock[1];
      k = blockId_z * Para_s->cuBlock[2];

    //ATHREAD serial operation replace GPU parallel mutli-thread in single thread-block
    for(i = 0; i < threadSum; i++) {

	  if(j<jmt&&i<imt&&k<km)	
      {	
        double am_f,dxdy,t_area;	
		double hduk = 0.0, hdvk = 0.0;	
      	am_f=1.0e0;	
     // 	d_hduk2_s[k *jmt *imt+j *imt+i]=0.0;	
      //	d_hdvk2_s[k *jmt *imt+j *imt+i]=0.0;
	    hduk = 0.0;
		hdvk = 0.0;
        if(j>=(Para_s->d_jb-2)&&j<=Para_s->d_je&&i>=(Para_s->d_ib-2)&&i<=Para_s->d_ie)	
        {	
        	if((k<=Para_s->d_kmu[j *imt+i]-1)){

			//d_hduk2_s[k *jmt *imt+j *imt+i]=(((Para_s->d_duc[j *imt+i]+
			hduk = (((Para_s->d_duc[j *imt+i]+
						Para_s->d_dum[j *imt+i]) *Para_s->d_up[k *jmt *imt+j *imt+i]+
						Para_s->d_dun[j *imt+i] *Para_s->d_up[k *jmt *imt+(j-1) *imt+i]+
						Para_s->d_dus[j *imt+i] *Para_s->d_up[k *jmt *imt+(j+1) *imt+i]+
						Para_s->d_due[j *imt+i] *Para_s->d_up[k *jmt *imt+j *imt+i+1]+
						Para_s->d_duw[j *imt+i] *Para_s->d_up[k *jmt *imt+j *imt+i-1])+
						(Para_s->d_dmc[j *imt+i] *Para_s->d_vp[k *jmt *imt+j *imt+i]+
						Para_s->d_dmn[j *imt+i] *Para_s->d_vp[k *jmt *imt+(j-1) *imt+i]+
						Para_s->d_dms[j *imt+i] *Para_s->d_vp[k *jmt *imt+(j+1) *imt+i]+
						Para_s->d_dme[j *imt+i] *Para_s->d_vp[k *jmt *imt+j *imt+i+1]+
						Para_s->d_dmw[j *imt+i] *Para_s->d_vp[k *jmt *imt+j *imt+i-1]));	
      	
			//d_hdvk2_s[k *jmt *imt+j *imt+i]=(((Para_s->d_duc[j *imt+i]+
			hdvk=(((Para_s->d_duc[j *imt+i]+
						Para_s->d_dum[j *imt+i]) *Para_s->d_vp[k *jmt *imt+j *imt+i]+
						Para_s->d_dun[j *imt+i] *Para_s->d_vp[k *jmt *imt+(j-1) *imt+i]+
						Para_s->d_dus[j *imt+i] *Para_s->d_vp[k *jmt *imt+(j+1) *imt+i]+
						Para_s->d_due[j *imt+i] *Para_s->d_vp[k *jmt *imt+j *imt+i+1]+
						Para_s->d_duw[j *imt+i] *Para_s->d_vp[k *jmt *imt+j *imt+i-1])-
						(Para_s->d_dmc[j *imt+i] *Para_s->d_up[k *jmt *imt+j *imt+i]+
						Para_s->d_dmn[j *imt+i] *Para_s->d_up[k *jmt *imt+(j-1) *imt+i]+
						Para_s->d_dms[j *imt+i] *Para_s->d_up[k *jmt *imt+(j+1) *imt+i]+
						Para_s->d_dme[j *imt+i] *Para_s->d_up[k *jmt *imt+j *imt+i+1]+
						Para_s->d_dmw[j *imt+i] *Para_s->d_up[k *jmt *imt+j *imt+i-1]));
			}
        }
      	
    		//d_d2uk_s[k *jmt *imt+j *imt+i]=am_f *Para_s->d_amf[j *imt+i] *d_hduk2_s[k *jmt *imt+j *imt+i];
		d_d2uk_s[k *jmt *imt+j *imt+i]=hduk;
      	//d_d2vk_s[k *jmt *imt+j *imt+i]=am_f *Para_s->d_amf[j *imt+i] *d_hdvk2_s[k *jmt *imt+j *imt+i];	
      	d_d2vk_s[k *jmt *imt+j *imt+i]=hdvk;
      }	

    }

  }

    ldm_free(Para_s, sizeof(hdiffu_del4_0_Para));
}
