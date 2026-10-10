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

typedef struct barotr2_tide_kernelPara{
  int cuBlock[3];
	int cuGrid[3];

  double d_d0;
	double *d_wka;
  double *d_ubp;
	double *d_vbp;
  double *d_ohbu;
  double *d_h0;
  double *d_cd_tide;
  double *d_viv;
  int    *d_kmu;
  double *d_ric;
  double *d_cd_inwave;
  double *d_depthstd;
}barotr2_tide_Para;

//dynamic atomic definition
#define CP_CORE      64
#define RCV_CORE     0
#define DATAPAKAGE   1

//ATHREAD global definition 
__thread_local int warp_size = 32; 
__thread_local int blockIndex;

void barotr2_tide_(barotr2_tide_Para *Paras){

    double d_ric;
    double d_d0;
    double *d_wka2;
    double *d_wka3;
    double *d_ubp;
    double *d_vbp;
    double *d_ohbu;
    double *d_h0;
    double *d_cd_tide;
    double *d_viv;
    int    *d_kmu;
    double *d_cd_inwave;
    double *d_depthstd;


    //ATHREAD basic operation
    int i,j,k;
    int blockId_x,blockId_y,blockId_z;
    blockIndex = CRTS_tid;
    barotr2_tide_Para *Para_s = (barotr2_tide_Para*)ldm_malloc(sizeof(barotr2_tide_Para));
    CRTS_dma_iget(Para_s, Paras, sizeof(barotr2_tide_Para), &dma_rply);
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

        d_wka2      = &Para_s->d_wka[2 *jmt *imt+j *imt];
        d_wka3      = &Para_s->d_wka[3 *jmt *imt+j *imt];
        d_vbp       = &Para_s->d_vbp[j * imt];
        d_ubp       = &Para_s->d_ubp[j * imt];
        d_ohbu      = &Para_s->d_ohbu[j*imt];
        d_h0        = &Para_s->d_h0[j*imt];
        d_cd_tide   = &Para_s->d_cd_tide[j*imt];
        d_viv       = &Para_s->d_viv[j*imt];
        d_cd_inwave = &Para_s->d_cd_inwave[j*imt];
        d_depthstd  = &Para_s->d_depthstd[j*imt];

        for (int i = 0; i < imt; i++) {


#ifdef TIDEFORCE
            if(i> 1&&i<imt-2&&j> 1&&j<jmt-2)
            { 
                d_wka2[i]=d_wka2[i]-pow((pow(d_ubp[i],2)+pow(d_vbp[i],2)),0.5)*
                          d_ubp[i]/(1/d_ohbu[i]+d_h0[i])*d_cd_tide[i];

                d_wka3[i]=d_wka3[i]-pow((pow(d_ubp[i],2)+pow(d_vbp[i],2)),0.5)*
                          d_vbp[i]/(1/d_ohbu[i]+d_h0[i])*d_cd_tide[i];
#if 0  
//使用sqrt(d_ubp[i]*d_ubp[i]+d_vbp[i]*d_vbp[i])，数据不一致                        
                d_wka2[i]=d_wka2[i]-sqrt(d_ubp[i]*d_ubp[i]+d_vbp[i]*d_vbp[i])*
                          d_ubp[i]/(1/d_ohbu[i]+d_h0[i])*d_cd_tide[i];

                d_wka3[i]=d_wka3[i]-sqrt(d_ubp[i]*d_ubp[i]+d_vbp[i]*d_vbp[i])*
                          d_vbp[i]/(1/d_ohbu[i]+d_h0[i])*d_cd_tide[i];
#endif

            }
#endif

#ifdef TIDEFORCEGW
            if(i> 1&&i<imt-2&&j> 1&&j<jmt-2)
            {
//CD_InWave(I,J,IBLOCK)=tidetune*(2*pi/10000.0)*depthstd(i,j,iblock)*depthstd(i,j,iblock)*SQRT(abs(ric(i,j,int(kmu(i,j,iblock))-1,iblock)))
                double tidetune=0.5;
                d_ric = Para_s->d_ric[((int)(Para_s->d_kmu[j *imt+i])-2)* imt * jmt + j * imt + i];
                d_cd_inwave[i]=tidetune*(2*pi/10000.0)*d_depthstd[i]*d_depthstd[i]*sqrt(fabs(d_ric));

                if (d_viv[i] > 0.5)
                {
                    if (1/d_ohbu[i] < 1000) 
                        d_cd_inwave[i] = 0.0;
                }
        //WKA(I,J,3,IBLOCK)=WKA(I,J,3,IBLOCK)-CD_InWave(I,J,IBLOCK)*UBP(I,J,IBLOCK)/((1/OHBU(I,J,IBLOCK)+H0(I,J,IBLOCK))*D0) 
        // WKA(I,J,4,IBLOCK)=WKA(I,J,4,IBLOCK)-CD_InWave(I,J,IBLOCK)*VBP(I,J,IBLOCK)/((1/OHBU(I,J,IBLOCK)+H0(I,J,IBLOCK))*D0)

                d_wka2[i] = d_wka2[i]-d_cd_inwave[i]*d_ubp[i]/((1/d_ohbu[i]+d_h0[i])*Para_s->d_d0);
                d_wka3[i] = d_wka3[i]-d_cd_inwave[i]*d_vbp[i]/((1/d_ohbu[i]+d_h0[i])*Para_s->d_d0);

            }
#endif
          }
    }
    ldm_free(Para_s, sizeof(barotr2_tide_Para));
}