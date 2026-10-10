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


typedef struct tracer_norm4_Para{
    int     cuGrid[3];
    int     cuBlock[3];
    int     *d_kmt;
    int     padd;
    double  d_od0cp;
    double  d_gamma;
    double  d_aidif;
    double  *d_stf;
    double  *d_swv;
    double  *d_nswv;
    double  *d_seaice;
    double  *d_sst;
    double  *d_atb;
    double  *d_odzp;
    double  *d_tf;
    double  *d_net;
#ifdef SST_restore
    double *d_sst_restore;
#endif
#ifdef COUP
    double *d_tsf;
#endif
}norm4Para;

__thread_local int blockIndex;


void tracer_norm4_kernel_(norm4Para *Paras){

    //ATHREAD basic operation
    int blockId_x,blockId_y,blockId_z;
    int i,j,k;
    blockIndex = CRTS_tid;
//    double d_tf[imt];
//    double d_stf[imt];
//    double d_net[imt];
    double *d_tf;
    double *d_stf;
    double *d_net;
    double *d_tsf;
    double *d_sst_restore;
    double  d_od0cp = Paras->d_od0cp;
    double  d_gamma = Paras->d_gamma;
    double  d_aidif = Paras->d_aidif;
    double  d_odzp = Paras->d_odzp[0];
    //cuda grid and block affine translation
    int blockSum   = Paras->cuGrid[0] * Paras->cuGrid[1] * Paras->cuGrid[2];
    int quotient   = blockSum / 64;
    int remainder  = blockSum % 64;
    int blockStart = quotient * blockIndex;
    if(remainder > 0)
        blockStart += (remainder > blockIndex ? blockIndex : remainder);
    int blockEnd = blockStart + quotient;
    if(remainder > 0)
        blockEnd   += (remainder > blockIndex ? 1 : 0);

    for(int block = blockStart; block < blockEnd; block++){
        blockId_z = block / (Paras->cuGrid[0] * Paras->cuGrid[1]);
        blockId_y = (block - blockId_z * (Paras->cuGrid[0] * Paras->cuGrid[1])) / Paras->cuGrid[0];
        blockId_x = block % Paras->cuGrid[0];
        j = blockId_y * Paras->cuBlock[1];

//        CRTS_dma_get(d_stf, &Paras->d_stf[ j * imt], sizeof(double) * imt);
//        CRTS_dma_get(d_tf, &Paras->d_tf[j * imt], sizeof(double) * imt);
    //    CRTS_dma_get(d_net, &Paras->d_net[j * imt], sizeof(double) * imt);
        d_stf = &Paras->d_stf[ j * imt];
        d_tf  = &Paras->d_tf[j * imt];
        d_net = &Paras->d_net[j * imt];
#ifdef COUP
        d_tsf = &Paras->d_tsf[j * imt];
#endif
#ifdef SST_restore
        d_sst_restore = &Paras->d_sst_restore[j * imt];
#endif
        for (i = 0; i < imt; i++) {
            int tid = j*imt+i;
            d_net[i] = 0.0;
            if(Paras->d_kmt[tid] > 0){
#ifdef COUP
                d_stf[i] = d_tsf[i];
#else
               // STF(I,J,IBLOCK) = (SWV(I,J,IBLOCK) + NSWV(I,J,IBLOCK))*OD0CP
                d_stf[i] = (Paras->d_swv[tid] + Paras->d_nswv[tid])*d_od0cp;
#endif
                //+Paras->d_seaice[tid]*d_gamma*(Paras->d_sst[tid]-Paras->d_atb[tid])/d_odzp;

                //sst_restore(I,J,IBLOCK) = (7.0/3.0)*GAMMA*seaice(i,j,IBLOCK)*(SST(I,J,IBLOCK) - ATB(I,J,1,1,IBLOCK)) + &
                //                          GAMMA*(1.0-seaice(i,j,IBLOCK))*(SST(I,J,IBLOCK) - ATB(I,J,1,1,IBLOCK))
#ifdef SST_restore
                d_sst_restore[i]  = (7.0/3.0) * d_gamma * Paras->d_seaice[tid] *
                            (Paras->d_sst[tid] - Paras->d_atb[tid]) +
                                d_gamma * (1.0 - Paras->d_seaice[tid]) *
                                    (Paras->d_sst[tid] - Paras->d_atb[tid]);
#endif
                d_tf[i] = d_tf[i] + d_stf[i] * d_odzp * (1.0 - d_aidif);
#ifdef SST_restore
                d_tf[i] = d_tf[i] + d_sst_restore[i];
#endif
                d_net[i] = d_stf[i] * d_odzp;
            }
        }

//        CRTS_dma_put(&Paras->d_tf[j * imt], d_tf, sizeof(double) * imt);
//        CRTS_dma_put(&Paras->d_stf[j * imt], d_stf, sizeof(double) * imt);
//        CRTS_dma_put(&Paras->d_net[j * imt], d_net, sizeof(double) * imt);
    }
    return;
}