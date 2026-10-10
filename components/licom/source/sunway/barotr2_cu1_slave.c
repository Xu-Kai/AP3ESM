
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
//ATHREAD DMA transfer definition 
__thread_local unsigned int D_COUNT = 0;
__thread_local crts_rply_t dma_rply = 0;
__thread_local crts_rply_t l_rply = 0;
__thread_local crts_rply_t r_rply = 0;

//__userDef__ paraVarAttr(double,d_dtb,sizeof(double),readonly,no,no,no);	
//__userDef__ paraVarAttr(double,d_h0,sizeof(double)*jmt*imt,outonly,yes,no,no);	
//__userDef__ paraVarAttr(double,d_h0p,sizeof(double)*jmt*imt,readonly,no,no,no);	
//__userDef__ paraVarAttr(double,d_work,sizeof(double)*jmt*imt,readonly,no,no,no);	
//__userDef__ paraVarAttr(int,d_jmt,sizeof(int),readonly,no,no,no);	
//__userDef__ paraVarAttr(int,d_imt,sizeof(int),readonly,no,no,no);	
typedef struct barotr2_cu1_kernelPara{
	int cuBlock[3];
	int cuGrid[3];
	int d_imt;
	int d_jmt;
	double d_dtb;
	double *d_h0;
	double *d_h0p;
	double *d_work;
#ifdef TIDEFORCE
  double BETA_TIDE;
  double *d_h0_tide;
  double *d_elev;
#endif
}barotr2_cu1_Para;

typedef struct barotr2_kernelPara{
        int cuBlock[3];
        int cuGrid[3];
        int d_jb;
        int d_je;
        int d_ib;
        int d_ie;
        int d_imt;
        int d_jmt;
        int d_km;
        int rank;
        int *d_kmu;
        int pad;
        double d_am;
        double *d_viv;
        double *d_ebea;
        double *d_ebeb;
        double *d_fcor;
#ifdef TIDEFORCE
        double *d_h0_tide;
#endif 
        double *d_h0;
        double *d_ubp;
        double *d_vbp;
        double *d_dlub;
        double *d_dlvb;
        double *d_wka;
        double *d_work;
        double *d_wgp;
        double *d_pax;
        double *d_pay;
        double *d_pxb;
        double *d_pyb;
        double *d_whx;
        double *d_why;
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
        double *d_dxur;
        double *d_dyur;
        double *d_au0;
        double *d_ausw;
        double *d_aus;
        double *d_auw;
}barotr2_Para;




//dynamic atomic definition
#define CP_CORE      64
#define RCV_CORE     0
#define DATAPAKAGE   1

//ATHREAD global definition 
__thread_local int warp_size = 32; 
__thread_local int blockIndex;

void barotr2_cu1_(barotr2_cu1_Para *Paras){	

  //ATHREAD local variable definition
  int swGrid[3];
  int swBlock[3];
  int stride;
  int blockId_x,blockId_y,blockId_z;
  int threadId_x,threadId_y,threadId_z;

  //ATHREAD basic operation
    int i,j,k;
  blockIndex = CRTS_tid;
  barotr2_cu1_Para *Para_s = (barotr2_cu1_Para*)ldm_malloc(sizeof(barotr2_cu1_Para));
  CRTS_dma_iget(Para_s, Paras, sizeof(barotr2_cu1_Para), &dma_rply);
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

  double *d_h0_s = Para_s->d_h0;
#ifdef TIDEFORCE
  double *d_h0_tide = Para_s->d_h0_tide;
  double *d_elev = Para_s->d_elev;
#endif
  for(int block = blockStart; block < blockEnd; block++){
    blockId_z = block / (Para_s->cuGrid[0] * Para_s->cuGrid[1]);
    blockId_y = (block - blockId_z * (Para_s->cuGrid[0] * Para_s->cuGrid[1])) / Para_s->cuGrid[0];
    blockId_x = block % Para_s->cuGrid[0];
    j = blockId_y * Para_s->cuBlock[1];

    //ATHREAD serial operation replace GPU parallel mutli-thread in single thread-block
    for(i = 0; i < threadSum; i++) {
      
	  if(i<Para_s->d_imt&&j<Para_s->d_jmt)	
      {	 
      	d_h0_s[j *Para_s->d_imt+i]=Para_s->d_h0p[j *Para_s->d_imt+i]+Para_s->d_work[j *Para_s->d_imt+i] *Para_s->d_dtb;	
#ifdef TIDEFORCE
      	d_h0_tide[j *Para_s->d_imt+i]=d_h0_s[j *Para_s->d_imt+i]*Para_s->BETA_TIDE-Para_s->d_elev[j *Para_s->d_imt+i];
#endif
      }	

    }
  }

    ldm_free(Para_s, sizeof(barotr2_cu1_Para));
}

void barotr2_cu_(barotr2_Para *Paras){
  blockIndex = CRTS_tid;
  barotr2_Para *Para_s = (barotr2_Para*)ldm_malloc(sizeof(barotr2_Para));
  CRTS_dma_iget(Para_s, Paras, sizeof(barotr2_Para), &dma_rply);
  D_COUNT++;
  CRTS_dma_wait_value(&dma_rply, D_COUNT);
  
  double hduk = 0.0e0;
  double hdvk = 0.0e0;
  double gradx=0.0e0,grady=0.0e0;
  //double d_wka0[imt],d_wka1[imt];

  for(int j=blockIndex;j < jmt;j+=64)
  {
	  /*
	CRTS_dma_get(d_wka0, &Para_s->d_wka[j * imt], sizeof(double) * imt);
	CRTS_dma_get(d_wka1, &Para_s->d_wka[1 *jmt *imt+j *imt], sizeof(double) * imt);
	CRTS_dma_get(d_wka2, &Para_s->d_wka[2 *jmt *imt+j *imt], sizeof(double) * imt);
	CRTS_dma_get(d_wka3, &Para_s->d_wka[3 *jmt *imt+j *imt], sizeof(double) * imt);
	CRTS_dma_get(d_wka4, &Para_s->d_wka[4 *jmt *imt+j *imt], sizeof(double) * imt);
	CRTS_dma_get(d_wka5, &Para_s->d_wka[5 *jmt *imt+j *imt], sizeof(double) * imt);
*/
	for(int i=0;i < imt;i++)
	{
	   hduk = 0.0e0, hdvk = 0.0e0; 
            if (Para_s->d_kmu[j * imt + i] - 1 >= 0) { 
                if (j >= Para_s->d_jb - 1 && j < Para_s->d_je && i >= Para_s->d_ib - 1 && i < Para_s->d_ie) {
                    hduk = Para_s->d_am *
                           (((Para_s->d_duc[j * imt + i] + Para_s->d_dum[j * imt + i]) * Para_s->d_ubp[j * imt + i] + 
                             Para_s->d_dun[j * imt + i] * Para_s->d_ubp[(j - 1) * imt + i] + 
                             Para_s->d_dus[j * imt + i] * Para_s->d_ubp[(j + 1) * imt + i] + 
                             Para_s->d_due[j * imt + i] * Para_s->d_ubp[j * imt + i + 1] + 
                             Para_s->d_duw[j * imt + i] * Para_s->d_ubp[j * imt + i - 1]) +
                            (Para_s->d_dmc[j * imt + i] * Para_s->d_vbp[j * imt + i] + 
                             Para_s->d_dmn[j * imt + i] * Para_s->d_vbp[(j - 1) * imt + i] + 
                             Para_s->d_dms[j * imt + i] * Para_s->d_vbp[(j + 1) * imt + i] + 
                             Para_s->d_dme[j * imt + i] * Para_s->d_vbp[j * imt + i + 1] + 
                             Para_s->d_dmw[j * imt + i] * Para_s->d_vbp[j * imt + i - 1])) * Para_s->d_viv[j * imt + i];

                    hdvk = Para_s->d_am *
                           (((Para_s->d_duc[j * imt + i] + Para_s->d_dum[j * imt + i]) * Para_s->d_vbp[j * imt + i] + 
                             Para_s->d_dun[j * imt + i] * Para_s->d_vbp[(j - 1) * imt + i] + 
                             Para_s->d_dus[j * imt + i] * Para_s->d_vbp[(j + 1) * imt + i] + 
                             Para_s->d_due[j * imt + i] * Para_s->d_vbp[j * imt + i + 1] + 
                             Para_s->d_duw[j * imt + i] * Para_s->d_vbp[j * imt + i - 1]) -
                            (Para_s->d_dmc[j * imt + i] * Para_s->d_ubp[j * imt + i] + 
                             Para_s->d_dmn[j * imt + i] * Para_s->d_ubp[(j - 1) * imt + i] + 
                             Para_s->d_dms[j * imt + i] * Para_s->d_ubp[(j + 1) * imt + i] + 
                             Para_s->d_dme[j * imt + i] * Para_s->d_ubp[j * imt + i + 1] + 
                             Para_s->d_dmw[j * imt + i] * Para_s->d_ubp[j * imt + i - 1])) * Para_s->d_viv[j * imt + i];  
              }

            }
	    gradx = 0.0e0, grady = 0.0e0;
	    double wka0,wka0_0,wka1,wka1_1;
        if (i > 1 && i < imt-2 && j > 1 && j < jmt-2) { //jjr bug
#ifdef TIDEFORCE
            if (Para_s->d_kmu[j * imt + i] - 1 >= 0) { 
                gradx = Para_s->d_dxur[j * imt + i] * 0.5 *
                        (Para_s->d_h0_tide[(j + 1) * imt + i] -
                         Para_s->d_h0_tide[j * imt + i - 1] - 
                         Para_s->d_h0_tide[(j + 1) * imt + i - 1] +
                         Para_s->d_h0_tide[j * imt + i]);

                grady = Para_s->d_dyur[j * imt + i] * 0.5 *
                        (Para_s->d_h0_tide[(j + 1) * imt + i] - 
                         Para_s->d_h0_tide[j * imt + i - 1] + 
                         Para_s->d_h0_tide[(j + 1) * imt + i - 1] -
                         Para_s->d_h0_tide[j * imt + i]);
            }
#else
            if (Para_s->d_kmu[j * imt + i] - 1 >= 0) { 
                gradx = Para_s->d_dxur[j * imt + i] * 0.5 *
                        (Para_s->d_h0[(j + 1) * imt + i] -
                         Para_s->d_h0[j * imt + i - 1] - 
                         Para_s->d_h0[(j + 1) * imt + i - 1] +
                         Para_s->d_h0[j * imt + i]);

                grady = Para_s->d_dyur[j * imt + i] * 0.5 *
                        (Para_s->d_h0[(j + 1) * imt + i] - 
                         Para_s->d_h0[j * imt + i - 1] + 
                         Para_s->d_h0[(j + 1) * imt + i - 1] -
                          Para_s->d_h0[j * imt + i]);
            }
#endif       
            double gstar = (Para_s->d_wgp[j * imt + i] - 1.0) * g;
	    wka0 = hduk + gstar * gradx;
            Para_s->d_wka[j * imt + i] = wka0;
	    wka1 = hdvk+gstar *grady;
            Para_s->d_wka[1 * jmt * imt + j * imt + i] = wka1;

            Para_s->d_wka[4 * jmt * imt + j * imt + i] = hduk;
            Para_s->d_wka[5 * jmt * imt + j * imt + i] = hdvk;

                }

//        Para_s->d_h0[j * imt + i] = Para_s->d_h0p[j * imt + i]+ Para_s->d_work[j * imt + i] * Para_s->d_dtb;
        if (j == jmt - 1 || i == 0) {
            Para_s->d_work[j * imt + i] = 0.0e0;
        } else {
            Para_s->d_work[j * imt + i] =
                    p25 * Para_s->d_h0[j * imt + i] +
                    p25 * Para_s->d_h0[(j + 1) * imt + i] +
                    p25 * Para_s->d_h0[j * imt + i - 1] +
                    p25 * Para_s->d_h0[(j + 1) * imt + i - 1];
        }

        if (i > 1 && i < imt - 2 && j > 1 && j < jmt - 2) {
            wka0_0  = Para_s->d_viv[j * imt + i] *
                                   (wka0 + Para_s->d_dlub[j * imt + i] -
                                    Para_s->d_fcor[j * imt + i] * Para_s->d_vbp[j * imt + i] +
                                    Para_s->d_pax[j * imt + i] + Para_s->d_pxb[j * imt + i] -
                                    Para_s->d_work[j * imt + i] * Para_s->d_whx[j * imt + i]);
	    Para_s->d_wka[j * imt + i] = wka0_0;

            wka1_1 = 
                    Para_s->d_viv[j * imt + i] * 
                    (wka1 + Para_s->d_dlvb[j * imt + i] + 
                     Para_s->d_fcor[j * imt + i] * Para_s->d_ubp[j * imt + i] + 
                     Para_s->d_pay[j * imt + i] + Para_s->d_pyb[j * imt + i] - 
                     Para_s->d_work[j * imt + i] * Para_s->d_why[j * imt + i]); 
	    Para_s->d_wka[1 * jmt * imt + j * imt + i] = wka1_1;

            Para_s->d_wka[2 * jmt * imt + j * imt + i] = 
                    Para_s->d_ebea[j * imt + i] *  wka0_0- 
                    Para_s->d_ebeb[j * imt + i] * wka1_1;

            Para_s->d_wka[3 * jmt * imt + j * imt + i] = 
                    Para_s->d_ebea[j * imt + i] * wka1_1 + 
                    Para_s->d_ebeb[j * imt + i] * wka0_0;
        }    
//
    }
	}
  
  ldm_free(Para_s, sizeof(barotr2_Para));
}



