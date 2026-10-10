
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

//__userDef__ paraVarAttr(double,d_wka,sizeof(double)*jmt*imt*km,readonly,no,no,no);	
//__userDef__ paraVarAttr(double,d_dxur,sizeof(double)*ny_block*nx_block,readonly,no,no,no);	
//__userDef__ paraVarAttr(double,d_dyur,sizeof(double)*ny_block*nx_block,readonly,no,no,no);	
//__userDef__ paraVarAttr(int,d_kmu,sizeof(int)*ny_block*nx_block,readonly,no,no,no);	
//__userDef__ paraVarAttr(double,d_am_factor,sizeof(double)*km*jmt*imt,outonly,yes,no,no);	
typedef struct bhdel4_grad_kernelPara{
	int cuBlock[3];
	int cuGrid[3];
	int rank;
	int *d_kmu;
	double *d_wka;
	double *d_dxur;
	double *d_dyur;
	double *d_am_factor;
	double *d_z_curl;
	double *d_diverg;
}bhdel4_grad_Para;


//dynamic atomic definition
#define CP_CORE      64
#define RCV_CORE     0
#define DATAPAKAGE   1

//ATHREAD global definition 
__thread_local int warp_size = 32; 
__thread_local int blockIndex;

void bhdel4_grad_(bhdel4_grad_Para *Paras){

  //ATHREAD local variable definition
  int swGrid[3];
  int swBlock[3];
  int blockId_x,blockId_y,blockId_z;
  int threadId_x,threadId_y,threadId_z;
  double gradx1 = 0.0e0;
  double gradx2 = 0.0e0;
  double grady1 = 0.0e0;
  double grady2 = 0.0e0;
  double d_wka[imt];
  double d_wkapj1[imt];
  double d_am_factor[imt];
  double d_dxur[imt];
  double d_dyur[imt];

  //ATHREAD basic operation
    int i,j,k;
  blockIndex = CRTS_tid;
  bhdel4_grad_Para *Para_s = (bhdel4_grad_Para*)ldm_malloc(sizeof(bhdel4_grad_Para));
  CRTS_dma_iget(Para_s, Paras, sizeof(bhdel4_grad_Para), &dma_rply);
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

//        CRTS_dma_get(d_am_factor, &Para_s->d_am_factor[j * imt], sizeof(double) * imt);
//        CRTS_dma_get(d_wka, &Para_s->d_wka[j * imt], sizeof(double) * imt);
//        CRTS_dma_get(d_wkapj1, &Para_s->d_wka[(j + 1) * imt], sizeof(double) * imt);
//        CRTS_dma_get(d_dxur, &Para_s->d_dxur[j * imt], sizeof(double) * imt);
//        CRTS_dma_get(d_dyur, &Para_s->d_dyur[j * imt], sizeof(double) * imt);

        for (i = 0; i < imt; i++) {
            gradx1=0.0;
            gradx2=0.0;
            grady1=0.0;
            grady2=0.0;
            d_am_factor[i] = 1.0;
            if((1<=i)&&(j<(jmt-1)))
            {
                if(0<=Para_s->d_kmu[j *imt+i]-1)
                {
                    gradx1=Para_s->d_dxur[j*imt+i] *p5 *(Para_s->d_z_curl[(j+1)*imt+i]- Para_s->d_z_curl[j*imt+i-1]-
                            Para_s->d_z_curl[(j+1)*imt+i-1]+Para_s->d_z_curl[j*imt+i]);

                    gradx2=Para_s->d_dxur[j*imt+i] *p5 *(Para_s->d_diverg[(j+1)*imt+i]- Para_s->d_diverg[j*imt+i-1]-
                            Para_s->d_diverg[(j+1)*imt+i-1]+Para_s->d_diverg[j*imt+i]);

                    grady1=Para_s->d_dyur[j*imt+i] *p5 *(Para_s->d_z_curl[(j+1)*imt+i]-Para_s->d_z_curl[j*imt+i-1]+
                            Para_s->d_z_curl[(j+1)*imt+i-1]-Para_s->d_z_curl[j*imt+i]);

                    grady2=Para_s->d_dyur[j*imt+i] *p5 *(Para_s->d_diverg[(j+1)*imt+i]-Para_s->d_diverg[j*imt+i-1]+
                            Para_s->d_diverg[(j+1)*imt+i-1]-Para_s->d_diverg[j*imt+i]);


                    d_am_factor[i] = (gradx2 - grady1)*(gradx2 - grady1) +
                                (grady2 + gradx1) * (grady2 + gradx1);
                }
            }

        }

        CRTS_dma_put(&Para_s->d_am_factor[j * imt], d_am_factor, sizeof(double) * imt);
    }


    ldm_free(Para_s, sizeof(bhdel4_grad_Para));
}
