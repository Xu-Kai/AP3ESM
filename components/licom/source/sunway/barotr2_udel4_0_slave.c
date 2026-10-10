
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


typedef struct barotr2_udel4_0_kernelPara{
	int cuBlock[3];
	int cuGrid[3];
	int d_jb;
	int d_je;
	int d_ib;
	int d_ie;
	int d_imt;
	int d_jmt;
	int rank;
	int *d_kmu;
    int pad[2];
	double d_am;
	double *d_ubp;
	double *d_vbp;
	double *d_uarea_r;
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
}barotr2_udel4_0_Para;


//dynamic atomic definition
#define CP_CORE      64
#define RCV_CORE     0
#define DATAPAKAGE   1

//ATHREAD global definition 
__thread_local int warp_size = 32; 
__thread_local int blockIndex;

void barotr2_udel4_0_(barotr2_udel4_0_Para *Paras){	

  //ATHREAD local variable definition
  int swGrid[3];
  int swBlock[3];
  int stride;
  int blockId_x,blockId_y,blockId_z;
  int threadId_x,threadId_y,threadId_z;
//  double d_ubp[3 * imt];
//  double d_vbp[3 * imt];
  double *d_ubp;
  double *d_vbp;
  double d_d2uk[imt];
  double d_d2vk[imt];
  double *d_amf;
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
  double hduk=0.0e0;
  double hdvk=0.0e0;
  double am_f,t_area,scaling,u_area2,temp;


  //ATHREAD basic operation
    int i,j,k;
  blockIndex = CRTS_tid;
  barotr2_udel4_0_Para *Para_s = (barotr2_udel4_0_Para*)ldm_malloc(sizeof(barotr2_udel4_0_Para));
  CRTS_dma_iget(Para_s, Paras, sizeof(barotr2_udel4_0_Para), &dma_rply);
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

        d_ubp = &Para_s->d_ubp[(j-1) * imt];
        d_vbp = &Para_s->d_vbp[(j-1) * imt];
        d_duc = &Para_s->d_duc[j * imt];
        d_dum = &Para_s->d_dum[j * imt];
        d_dun = &Para_s->d_dun[j * imt];
        d_dus = &Para_s->d_dus[j * imt];
        d_due = &Para_s->d_due[j * imt];
        d_duw = &Para_s->d_duw[j * imt];
        d_dmc = &Para_s->d_dmc[j * imt];
        d_dmn = &Para_s->d_dmn[j * imt];
        d_dms = &Para_s->d_dms[j * imt];
        d_dme = &Para_s->d_dme[j * imt];
        d_dmw = &Para_s->d_dmw[j * imt];


        for (int i = 0; i < imt; i++) {
//            am_f=1.0e0;
	    	hduk=0.0e0;
  	    	hdvk=0.0e0;

            if(j>=Para_s->d_jb-2&&j<=Para_s->d_je&&i>=Para_s->d_ib-2&&i<=Para_s->d_ie)
            {
                if(Para_s->d_kmu[j *imt+i]-1>=0)
                {
                    hduk=((d_duc[i] + d_dum[i]) * d_ubp[imt+i]+d_dun[i] * d_ubp[i]+
                           d_dus[i] * d_ubp[2*imt+i] + d_due[i] * d_ubp[imt+i+1] + d_duw[i] * d_ubp[imt+i-1])+
                          (d_dmc[i] * d_vbp[imt+i] + d_dmn[i] * d_vbp[i] + d_dms[i] * d_vbp[2*imt+i]+
                           d_dme[i] * d_vbp[imt+i+1] + d_dmw[i] * d_vbp[imt+i-1]);

                    hdvk=((d_duc[i] + d_dum[i]) * d_vbp[imt+i] + d_dun[i] * d_vbp[i]+
                           d_dus[i] * d_vbp[2*imt+i] + d_due[i] * d_vbp[imt+i+1] + d_duw[i] * d_vbp[imt+i-1])-
                          (d_dmc[i] * d_ubp[imt+i] + d_dmn[i] * d_ubp[i] + d_dms[i] * d_ubp[2*imt+i]+
                           d_dme[i] * d_ubp[imt+i+1] + d_dmw[i] * d_ubp[imt+i-1]);
                }
            }
            //d_d2uk[i]=am_f *d_amf[i] *hduk;
            //d_d2vk[i]=am_f *d_amf[i] *hdvk;
        	d_d2uk[i]=hduk;
        	d_d2vk[i]=hdvk;

        }

        CRTS_dma_put(&Para_s->d_d2uk[j * imt], d_d2uk, sizeof(double) * imt);
        CRTS_dma_put(&Para_s->d_d2vk[j * imt], d_d2vk, sizeof(double) * imt);
    }


    ldm_free(Para_s, sizeof(barotr2_udel4_0_Para));
}
