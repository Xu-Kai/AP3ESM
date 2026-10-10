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

typedef struct kernelPara{
    int    cuBlock[3];
    int    cuGrid[3];
    double d_od0;
    double d_c0f;
    double d_cag;
    double d_sag;
    double *d_snlat;
    double *d_up;
    double *d_vp;
    double *d_bbcy;
    double *d_sbcy;
    double *d_sbcx;
    double *d_bbcx;
    double *d_su;
    double *d_sv;
    int *d_kmu;
}Para;

//dynamic atomic definition
#define CP_CORE      64
#define RCV_CORE     0
#define DATAPAKAGE   1

//ATHREAD global definition 
__thread_local int warp_size = 32; 
__thread_local int blockIndex;
__thread_local Para* Para_s;

void bclinc12_cu_slave_kernel_(Para *Paras){
	//ATHREAD basic operation

    int i,j,k;
    blockIndex = CRTS_tid;
    Para_s = (Para*)ldm_malloc(sizeof(Para));
    CRTS_dma_iget(Para_s, Paras, sizeof(Para), &dma_rply);
    D_COUNT++;
    CRTS_dma_wait_value(&dma_rply, D_COUNT);

    double d_od0 = Para_s->d_od0;
    double d_cag = Para_s->d_cag;
    double d_c0f = Para_s->d_c0f;
    double d_sag = Para_s->d_sag;

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

    for(int block = blockStart; block < blockEnd; block++) {
        int blockId_z = block / (Para_s->cuGrid[0] * Para_s->cuGrid[1]);
        int blockId_y = (block - blockId_z * (Para_s->cuGrid[0] * Para_s->cuGrid[1])) / Para_s->cuGrid[0];
        int blockId_x = block % Para_s->cuGrid[0];
        j = blockId_y * Para_s->cuBlock[1];
        // ATHREAD serial operation replace GPU parallel mutli-thread in single thread-block
        for(i = 0; i < threadSum; i++) {

        if (i >= 1 && i < imt - 1 && j >= 1 && j < jmt - 1) {
            int kmb = Para_s->d_kmu[j * imt + i] - 1;

            if (kmb >= 0) {
                double d_up_v = Para_s->d_up[kmb * jmt * imt + j * imt + i], d_vp_v = Para_s->d_vp[kmb * jmt * imt + j * imt + i];
		        double temp = d_c0f * sqrt(d_up_v * d_up_v + d_vp_v * d_vp_v);
		
                Para_s->d_sbcx[j * imt + i] = Para_s->d_su[j * imt + i] * d_od0;
                Para_s->d_sbcy[j * imt + i] = Para_s->d_sv[j * imt + i] * d_od0;

                Para_s->d_bbcx[j * imt + i] = temp *
                        (d_up_v * d_cag + Para_s->d_snlat[j * imt + i] * d_vp_v * d_sag);

                Para_s->d_bbcy[j * imt + i] = temp *
                        (d_vp_v * d_cag - Para_s->d_snlat[j * imt + i] * d_up_v * d_sag);
            }
        }
        }
    }
    ldm_free(Para_s, sizeof(Para));
    return;
}
