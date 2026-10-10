
//ATHREAD include file
#include <crts.h>
#include <slave.h>
#include <simd.h>
#include <math.h>
#include <stdio.h>
#include <float.h>
#include <string.h>
#include <stdbool.h>
//#include "gmp.h"
//#include "mpfr.h"
#include "swperf.h"

#include "param_mod.h"	
#include "constant_mod.h"
//ATHREAD DMA transfer definition 
__thread_local unsigned int D_COUNT = 0;
__thread_local crts_rply_t dma_rply = 0;
__thread_local crts_rply_t l_rply = 0;
__thread_local crts_rply_t r_rply = 0;


typedef struct hdiffu_del4_1_kernelPara{
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
	double *d_d2uk;
	double *d_d2vk;
	double *d_am_factor;
	double *d_uarea_r;
	double *d_fac_vis_uv;
	//double *d_hdukkm;
	//double *d_hdvkkm;
	//double *d_am_f_km;
	double d_dts;
	double d_value;
	double d_value2;
}hdiffu_del4_1_Para;


//dynamic atomic definition
#define CP_CORE      64
#define RCV_CORE     0
#define DATAPAKAGE   1

//ATHREAD global definition 
__thread_local int warp_size = 32; 
__thread_local int blockIndex;

void hdiffu_del4_1_(hdiffu_del4_1_Para *Paras){

  //ATHREAD local variable definition
  int swGrid[3];
  int swBlock[3];
  int stride;
  int blockId_x,blockId_y,blockId_z;
  int threadId_x,threadId_y,threadId_z;
  double am_f,t_area,t_area2,t_area4,scaling,scalbase,u_area2,d_am,t2,t4;
  //ATHREAD basic operation
    int i,j,k;
  blockIndex = CRTS_tid;
  hdiffu_del4_1_Para *Para_s = (hdiffu_del4_1_Para*)ldm_malloc(sizeof(hdiffu_del4_1_Para));
  CRTS_dma_iget(Para_s, Paras, sizeof(hdiffu_del4_1_Para), &dma_rply);
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

  double *d_wka_s = Para_s->d_wka;
  double *d_dlu_s = Para_s->d_dlu;
  double *d_dlv_s = Para_s->d_dlv;
	d_am = Para_s->d_am;


  for(int block = blockStart; block < blockEnd; block++){
    blockId_z = block / (Para_s->cuGrid[0] * Para_s->cuGrid[1]);
    blockId_y = (block - blockId_z * (Para_s->cuGrid[0] * Para_s->cuGrid[1])) / Para_s->cuGrid[0];
    blockId_x = block % Para_s->cuGrid[0];
      j = blockId_y * Para_s->cuBlock[1];
      k = blockId_z * Para_s->cuBlock[2];

    //ATHREAD serial operation replace GPU parallel mutli-thread in single thread-block
    for(i = 0; i < threadSum; i++) {

	  double hduk = 0.0, hdvk = 0.0;	
    	am_f = 1.0;
	  if(j<jmt&&i<imt&&k<km)	
      {	
      	//d_hduk2_s[k *jmt *imt+j *imt+i]=0.0;	
      	hduk=0.0;	
      	//d_hdvk2_s[k *jmt *imt+j *imt+i]=0.0;	
      	hdvk=0.0;
	  	if(j>=(Para_s->d_jb-2)&&j<(Para_s->d_je+1)&&i>=(Para_s->d_ib-2)&&i<(Para_s->d_ie+1))
	  	{
	  		scalbase = sqrt(Para_s->d_uarea_r[j *imt+i]);
	  		t_area = d_am * scalbase / (pi);
	  	//	t_area2 = t_area * t_area;
	  	//	t_area4 = t_area2 * t_area2;
	  		u_area2 = scalbase * scalbase;
	  		scaling = t_area * t_area * t_area * t_area * t_area;
	  		am_f = scaling * sqrt(Para_s->d_am_factor[k*jmt*imt+j*imt+i]);

	  		am_f = fmax(1.3e-5 * (u_area2 * u_area2)*Para_s->d_value2, am_f);
	  		am_f = fmin(1.0/80.0 * (u_area2 * u_area2 ) / Para_s->d_dts, am_f);
//			if(Para_s->rank == 10){
//				if(k == 2 && j == 26 && (i>=16) && (i<24)){
//					printf("i:%d,j:%d,k:%d,am:%f,scalbase:%.15lf,pi:%.15lf,t_area:%.16lf, scaling:%.16E, am_factor:%.16E, amf:%.16lf\n",
//						i,j,k,d_am,scalbase,pi,t_area,scaling,sqrt(Para_s->d_am_factor[k*jmt*imt+j*imt+i]),am_f);
//				}
//			}
	  	}
	  	am_f = fmax(Para_s->d_value, am_f);

        if(j>=Para_s->d_jb-1&&j<Para_s->d_je&&i>=Para_s->d_ib-1&&i<Para_s->d_ie)
        {	
      	  //d_hduk2_s[k *jmt *imt+j *imt+i]=Para_s->d_am *(((Para_s->d_duc[j *imt+i]+
      	  hduk=Para_s->d_fac_vis_uv[j *imt+i]*(((Para_s->d_duc[j *imt+i]+ Para_s->d_dum[j *imt+i]) *Para_s->d_d2uk[k *jmt *imt+j *imt+i]+
					Para_s->d_dun[j *imt+i] *Para_s->d_d2uk[k *jmt *imt+(j-1) *imt+i]+
					Para_s->d_dus[j *imt+i] *Para_s->d_d2uk[k *jmt *imt+(j+1) *imt+i]+
					Para_s->d_due[j *imt+i] *Para_s->d_d2uk[k *jmt *imt+j *imt+i+1]+
					Para_s->d_duw[j *imt+i] *Para_s->d_d2uk[k *jmt *imt+j *imt+i-1])+
				    (Para_s->d_dmc[j *imt+i] *Para_s->d_d2vk[k *jmt *imt+j *imt+i]+
					Para_s->d_dmn[j *imt+i] *Para_s->d_d2vk[k *jmt *imt+(j-1) *imt+i]+
					Para_s->d_dms[j *imt+i] *Para_s->d_d2vk[k *jmt *imt+(j+1) *imt+i]+
					Para_s->d_dme[j *imt+i] *Para_s->d_d2vk[k *jmt *imt+j *imt+i+1]+
					Para_s->d_dmw[j *imt+i] *Para_s->d_d2vk[k *jmt *imt+j *imt+i-1]));
      	
		  //d_hdvk2_s[k *jmt *imt+j *imt+i]=Para_s->d_am *(((Para_s->d_duc[j *imt+i]+
		  hdvk=Para_s->d_fac_vis_uv[j *imt+i]*(((Para_s->d_duc[j *imt+i]+
					Para_s->d_dum[j *imt+i]) *Para_s->d_d2vk[k *jmt *imt+j *imt+i]+
					Para_s->d_dun[j *imt+i] *Para_s->d_d2vk[k *jmt *imt+(j-1) *imt+i]+
					Para_s->d_dus[j *imt+i] *Para_s->d_d2vk[k *jmt *imt+(j+1) *imt+i]+
					Para_s->d_due[j *imt+i] *Para_s->d_d2vk[k *jmt *imt+j *imt+i+1]+
					Para_s->d_duw[j *imt+i] *Para_s->d_d2vk[k *jmt *imt+j *imt+i-1])-
				    (Para_s->d_dmc[j *imt+i] *Para_s->d_d2uk[k *jmt *imt+j *imt+i]+
					Para_s->d_dmn[j *imt+i] *Para_s->d_d2uk[k *jmt *imt+(j-1) *imt+i]+
					Para_s->d_dms[j *imt+i] *Para_s->d_d2uk[k *jmt *imt+(j+1) *imt+i]+
					Para_s->d_dme[j *imt+i] *Para_s->d_d2uk[k *jmt *imt+j *imt+i+1]+
					Para_s->d_dmw[j *imt+i] *Para_s->d_d2uk[k *jmt *imt+j *imt+i-1]));
        }

	  	hduk = -1.0 * am_f * hduk;
	  	hdvk = -1.0 * am_f * hdvk;

	  	if((k> Para_s->d_kmu[j *imt+i]-1))
        {	
      	  //d_hduk2_s[k *jmt *imt+j *imt+i]=0.0;	
      	  hduk=0.0;	
      	  //d_hdvk2_s[k *jmt *imt+j *imt+i]=0.0;	
      	  hdvk=0.0;	
        }
//		Para_s->d_hdukkm[k*jmt*imt+j*imt+i]=hduk;
//		Para_s->d_hdvkkm[k*jmt*imt+j*imt+i]=hdvk;
//		Para_s->d_am_f_km[k*jmt*imt+j*imt+i]=am_f;
        if(2<=j&&j<jmt-2&&2<=i&&i<imt-2)
        {	
      	  //d_dlv_s[k *jmt *imt+j *imt+i]=d_dlv_s[k *jmt *imt+j *imt+i]+d_hdvk2_s[k *jmt *imt+j *imt+i];	
      	  d_dlv_s[k *jmt *imt+j *imt+i]=d_dlv_s[k *jmt *imt+j *imt+i]+hdvk;	
      	  //d_dlu_s[k *jmt *imt+j *imt+i]=d_dlu_s[k *jmt *imt+j *imt+i]+d_hduk2_s[k *jmt *imt+j *imt+i];	
      	  d_dlu_s[k *jmt *imt+j *imt+i]=d_dlu_s[k *jmt *imt+j *imt+i]+hduk;	
        }	
      	d_wka_s[k *jmt *imt+j *imt+i]=Para_s->d_c0f *sqrt(Para_s->d_up[k *jmt *imt+j *imt+i] *
				Para_s->d_up[k *jmt *imt+j *imt+i]+Para_s->d_vp[k *jmt *imt+j *imt+i] *Para_s->d_vp[k *jmt *imt+j *imt+i]);	
      }	

    }

  }

    ldm_free(Para_s, sizeof(hdiffu_del4_1_Para));
}
