
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


typedef struct halo_2dr8_kernelPara{
    int numMsgSend;
    int bufSizeSend;
    int rank;
    int i0;
    int *sizeSend;
    int *sendAddr;
    double *d_array;
    double *d_bufSendR8;
}halo_2dr8_Para;


typedef struct halo_2dr8_recv_kernelPara{
    int numMsgRecv;
    int bufSizeRecv;
    int nxGlobal;
    int rank;
    int i0;
    int *sizeRecv;
    int *recvAddr;
    double *bufRecv2DR8;
    double *bufTripole2DR8;
    double *array;
}halo_2dr8_recv_Para;

typedef struct halo_2dr8_localcopy_1_kernelPara{
    int numLocalCopies;
    int nxGlobal;
    int rank;
    int I0;
    int ioffset;
    int joffset;
    int isign;
    int *srcLocalAddr;
    int *dstLocalAddr;
    double fill;
    double *bufTripole2D;
    double *array;
}halo_2dr8_localcopy_1_Para;


typedef struct halo_3dr8_send_kernelPara{
    int numMsgSend;
    int bufSizeSend;
    int strideArr;
    int strideSend;
    int I0;
    int I2;
    int *sendAddr;
    int *sizeSend;
    double *bufSend3D;
    double *array;
}halo_3dr8_send_Para;


typedef struct halo_3dr8_recv_kernelPara{
    int numMsgRecv;
    int bufSizeRecv;
    int strideArr;
    int strideBuf;
    int strideRecv;
    int nxGlobal;
    int I0;
    int I2;
    int *sizeRecv;
    int *recvAddr;
    double *bufTripole3D;
    double *bufRecv3D;
    double *array;
}halo_ddr8_recv_Para;


typedef struct halo_3dr8_localcopy_kernelPara{
    int strideArr;
    int strideBuf;
    int nxGlobal;
    int I0;
    int I2;
    int numLocalCopies;
    int ioffset;
    int joffset;
    int isign;
    int *srcLocalAddr;
    int *dstLocalAddr;
    double fill;
    double *bufTripole3D;
    double *array;
}halo_3dr8_localcopy_Para;

double  sign(double x, double y){
    return (y) >= 0.0e0 ? fabs(x) : -fabs(x);
}

void halo_2dr8_send_slave_(halo_2dr8_Para *Paras){

    halo_2dr8_Para *Para_s = (halo_2dr8_Para*)ldm_malloc(sizeof(halo_2dr8_Para));
    CRTS_dma_iget(Para_s, Paras, sizeof(halo_2dr8_Para), &dma_rply);
    D_COUNT++;
    CRTS_dma_wait_value(&dma_rply, D_COUNT);
    int bufSizeSend = Para_s->bufSizeSend;
    int numMsgSend = Para_s->numMsgSend;
    int i0 = Para_s->i0;
    const int (*sendAddr)[3] = Para_s->sendAddr;

    for (int nmsg = CRTS_tid; nmsg < numMsgSend; nmsg += 64) {
        int msgSize = Para_s->sizeSend[nmsg];
        for (int n = 0; n < msgSize; ++n) {
            //         iSrc     = halo%sendAddr(1,n,nmsg)
            //         jSrc     = halo%sendAddr(2,n,nmsg)
            // halo%sendAddr(3,bufSizeSend,numMsgSend)
            //int iSrc     = Para_s->sendAddr[nmsg * bufSizeSend * 3 + n * 3];
            //int jSrc     = Para_s->sendAddr[nmsg * bufSizeSend * 3 + n * 3 + 1];
            int iSrc     = sendAddr[nmsg * bufSizeSend + n][0] - 1;
            int jSrc     = sendAddr[nmsg * bufSizeSend + n][1] - 1;
            // srcBlock = sendAddr[nmsg * bufSizeSend + n][2];
            // bufSendR8(bufSizeSend, numMsgSend)
            Para_s->d_bufSendR8[nmsg * bufSizeSend + n] = Para_s->d_array[jSrc * i0 + iSrc];
        }
    }

    ldm_free(Para_s, sizeof(halo_2dr8_Para));
}


void halo_2dr8_recv_slave_(halo_2dr8_recv_Para *Paras){

    halo_2dr8_recv_Para *Para_s = (halo_2dr8_recv_Para*)ldm_malloc(sizeof(halo_2dr8_recv_Para));
    CRTS_dma_iget(Para_s, Paras, sizeof(halo_2dr8_recv_Para), &dma_rply);
    D_COUNT++;
    CRTS_dma_wait_value(&dma_rply, D_COUNT);
    int bufSizeRecv = Para_s->bufSizeRecv;
    int numMsgRecv = Para_s->numMsgRecv;
    int I0 = Para_s->i0;
    const int (*recvAddr)[3] = Para_s->recvAddr;
    int *sizeRecv = Para_s->sizeRecv;
    int nxGlobal = Para_s->nxGlobal;
    double *bufTripole2D = Para_s->bufTripole2DR8;
    double *bufRecv2D = Para_s->bufRecv2DR8;
    double *array = Para_s->array;

    for (int nmsg = CRTS_tid; nmsg < numMsgRecv; nmsg += 64) {
        const int sizeRecv_tmp = sizeRecv[nmsg];
        for (int n = 0; n < sizeRecv_tmp; ++n) {
            int iDst     = recvAddr[nmsg * bufSizeRecv + n][0] - 1;
            int jDst     = recvAddr[nmsg * bufSizeRecv + n][1] - 1;
            int dstBlock = recvAddr[nmsg * bufSizeRecv + n][2];
            if (dstBlock > 0) {
                // if(numLocalCopies == 0) printf("dst block >0 %d\n", numLocalCopies);
                array[jDst * I0 + iDst] = bufRecv2D[nmsg * bufSizeRecv + n];
            } else if (dstBlock < 0) {
                // if(numLocalCopies == 0)
                //  printf("dst block <0 %d %d\n", numLocalCopies, rank);

                // tripole
                bufTripole2D[jDst * nxGlobal + iDst] =  bufRecv2D[nmsg * bufSizeRecv + n];
            }
        }
    }

    ldm_free(Para_s, sizeof(halo_2dr8_recv_Para));
}

void halo_update2d_localcopy_slave_1_(halo_2dr8_localcopy_1_Para *args)
{

    double *array = args->array;
    double *bufTripole2D = args->bufTripole2D;
    const int (*srcLocalAddr)[3] = args->srcLocalAddr;
    const int (*dstLocalAddr)[3] = args->dstLocalAddr;
    int numLocalCopies = args->numLocalCopies;
    int nxGlobal = args->nxGlobal;
    int I0 = args->I0;
    double fill = args->fill;

    int local_block = ((numLocalCopies  + 63) / 64) > (128 / sizeof(array[0])) ? (numLocalCopies  + 63) / 64 : (128 / sizeof(array[0]));
    int local_start = CRTS_tid * local_block;
    int local_end = (CRTS_tid + 1) * local_block < numLocalCopies?(CRTS_tid + 1) * local_block:numLocalCopies;


    for (int nmsg = local_start; nmsg < local_end; ++nmsg) {
        int iSrc     = srcLocalAddr[nmsg][0] - 1;
        int jSrc     = srcLocalAddr[nmsg][1] - 1;
        int srcBlock = srcLocalAddr[nmsg][2];
        int iDst     = dstLocalAddr[nmsg][0] - 1;
        int jDst     = dstLocalAddr[nmsg][1] - 1;
        int dstBlock = dstLocalAddr[nmsg][2];

        if (srcBlock > 0) {
            if (dstBlock > 0) {
                array[jDst * I0 + iDst] = array[jSrc * I0 + iSrc];
            } else if (dstBlock < 0) {
                // tripole copy into buffer
                bufTripole2D[jDst * nxGlobal + iDst] = array[jSrc * I0 + iSrc];
            }
        } else if (srcBlock == 0) {
            array[jDst * I0 + iDst] = fill;
        }
    }
}

void halo_update2d_localcopy_slave_2_( halo_2dr8_localcopy_1_Para *args)
{


    double *array = args->array;
    double *bufTripole2D = args->bufTripole2D;
    const int (*srcLocalAddr)[3] = args->srcLocalAddr;
    const int (*dstLocalAddr)[3] = args->dstLocalAddr;
    int numLocalCopies = args->numLocalCopies;
    int nxGlobal = args->nxGlobal;
    int I0 = args->I0;
    double fill = args->fill;
    int ioffset = args->ioffset;
    int joffset = args->joffset;
    int isign = args->isign;


    int local_block = ((numLocalCopies  + 63) / 64) > (128 / sizeof(array[0])) ? (numLocalCopies  + 63) / 64 : (128 / sizeof(array[0]));
    int local_start = CRTS_tid * local_block;
    int local_end = (CRTS_tid + 1) * local_block < numLocalCopies?(CRTS_tid + 1) * local_block:numLocalCopies;



    for (int nmsg = local_start; nmsg < local_end; ++nmsg) {
        int srcBlock = srcLocalAddr[nmsg][2];
        if (srcBlock < 0) {
            int iSrc = srcLocalAddr[nmsg][0] - 1;
            int jSrc = srcLocalAddr[nmsg][1] - 1;

            int iDst = dstLocalAddr[nmsg][0] - 1;
            int jDst = dstLocalAddr[nmsg][1] - 1;
            // dstBlock = dstLocalAddr[nmsg][2];

            iSrc -= ioffset;
            jSrc -= joffset;
            if (iSrc == -1) {
                iSrc = nxGlobal - 1;
            }
            if (jSrc > -1) {

                if(isign == 1)
                {
                    if(jSrc != 0){
                        array[jDst * I0 + iDst] = isign * bufTripole2D[jSrc * nxGlobal + iSrc];
                    }
                    else{
                        if(iSrc < nxGlobal - 1)
                        {
                            double x1 = bufTripole2D[iSrc];
                            double x2 =  bufTripole2D[nxGlobal - iSrc - 2];

                            if(isign == 1)
                            {
                                array[jDst * I0 + iDst] = isign * 0.5 * (x1 + x2);
                            }
                        }
                        else
                        {
                            array[jDst * I0 + iDst] = bufTripole2D[iSrc];
                        }
                    }
                }
                else if(isign == -1)
                {
                    array[jDst * I0 + iDst] = isign * bufTripole2D[jSrc * nxGlobal + iSrc];
                }
            }
        }
    }

}

void halo_update3d_send_slave_(halo_3dr8_send_Para *args)
{
    double *array = args->array;
    double *bufSend3D = args->bufSend3D;
    const int (*sendAddr)[3] = args->sendAddr;
    int *sizeSend = args->sizeSend;
    int numMsgSend = args->numMsgSend;
    int bufSizeSend = args->bufSizeSend;
    int strideArr = args->strideArr;
    int strideSend = args->strideSend;
    int I0 = args->I0;
    int I2 = args->I2;

    int total_send = I2 * numMsgSend;
    int local_block = (total_send + 63) / 64;
    int local_start = CRTS_tid * local_block;
    int local_end = (CRTS_tid + 1) * local_block < total_send?(CRTS_tid + 1) * local_block:total_send;

    for (int t = local_start; t < local_end; t += 1) {

        int nmsg = t / I2;
        int k = t % I2;
        const int size_Send = sizeSend[nmsg];
        for (int n = 0; n < size_Send; ++n) {
            int i =  k*size_Send + n;
            int iSrc     = sendAddr[nmsg * bufSizeSend + n][0] - 1;
            int jSrc     = sendAddr[nmsg * bufSizeSend + n][1] - 1;
            // srcBlock = sendAddr[nmsg * bufSizeSend + n][2];

            bufSend3D[nmsg * strideSend + i] = array[k*strideArr + jSrc*I0 + iSrc];

        }
    }

}

void halo_update3d_recv_slave_(halo_ddr8_recv_Para *args)
{


    double *array = args->array;
    double *bufTripole3D = args->bufTripole3D;
    double *bufRecv3D = args->bufRecv3D;
    const int (*recvAddr)[3] = args->recvAddr;
    int *sizeRecv = args->sizeRecv;
    int numMsgRecv = args->numMsgRecv;
    int bufSizeRecv = args->bufSizeRecv;
    int strideArr = args->strideArr;
    int strideBuf = args->strideBuf;
    int strideRecv = args->strideRecv;
    int nxGlobal = args->nxGlobal;
    int I0 = args->I0;
    int I2 = args->I2;

    int total_recv = I2 * numMsgRecv;
    int local_block = (total_recv + 63) / 64;
    int local_start = CRTS_tid * local_block;
    int local_end = (CRTS_tid + 1) * local_block < total_recv?(CRTS_tid + 1) * local_block:total_recv;


    for (int t = local_start; t < local_end; t += 1) {

        int nmsg = t / I2;
        int k = t % I2;
        int sizeRecv_tmp = sizeRecv[nmsg];
        for (int n = 0; n < sizeRecv_tmp; ++n) {
            int iDst     = recvAddr[nmsg * bufSizeRecv + n][0] - 1;
            int jDst     = recvAddr[nmsg * bufSizeRecv + n][1] - 1;
            int dstBlock = recvAddr[nmsg * bufSizeRecv + n][2];
            int i = k * sizeRecv_tmp + n;
            if (dstBlock > 0) {
                // for (int k = 0; k < I2; ++k) {
                array[k*strideArr + jDst*I0 + iDst] = bufRecv3D[nmsg * strideRecv + i];
                // i += 1;
                // }
            } else if (dstBlock < 0) {
                // tripole
                // for (int k = 0; k < I2; ++k) {

                bufTripole3D[k * strideBuf + jDst * nxGlobal + iDst] =  bufRecv3D[nmsg * strideRecv + i];
                // i += 1;
                // }
            }
        }

    }

}

void halo_update3d_localcopy_slave_1_(halo_3dr8_localcopy_Para *args)
{


    double *array = args->array;
    double *bufTripole3D = args->bufTripole3D;
    const int (*srcLocalAddr)[3] = args->srcLocalAddr;
    const int (*dstLocalAddr)[3] = args->dstLocalAddr;
    int strideArr = args->strideArr;
    int strideBuf = args->strideBuf;
    int nxGlobal = args->nxGlobal;
    int I0 = args->I0;
    int I2 = args->I2;
    int numLocalCopies = args->numLocalCopies;
    double fill = args->fill;

    int total_copy = I2 * numLocalCopies;
    int local_block = (total_copy + 63) / 64;
    int local_start = CRTS_tid * local_block;
    int local_end = (CRTS_tid + 1) * local_block < total_copy?(CRTS_tid + 1) * local_block:total_copy;


    for(int t = local_start; t < local_end; t++)
    {
        int nmsg = t / I2;
        int k = t % I2;
        int iSrc     = srcLocalAddr[nmsg][0] - 1;
        int jSrc     = srcLocalAddr[nmsg][1] - 1;
        int srcBlock = srcLocalAddr[nmsg][2];
        int iDst     = dstLocalAddr[nmsg][0] - 1;
        int jDst     = dstLocalAddr[nmsg][1] - 1;
        int dstBlock = dstLocalAddr[nmsg][2];
        if (srcBlock > 0) {
            if (dstBlock > 0) {
                array[k*strideArr + jDst*I0 + iDst] = array[k*strideArr + jSrc*I0 + iSrc];
            } else if (dstBlock < 0) {
                // tripole copy into buffer
                bufTripole3D[k * strideBuf + jDst * nxGlobal + iDst] = array[k*strideArr + jSrc*I0 + iSrc];
            }
        } else if (srcBlock == 0) {
            array[k*strideArr + jDst*I0 + iDst] = fill;
        }

    }

}

void halo_update3d_localcopy_slave_2_(halo_3dr8_localcopy_Para *args)
{

    double *array = args->array;
    double *bufTripole3D = args->bufTripole3D;
    const int (*srcLocalAddr)[3] = args->srcLocalAddr;
    const int (*dstLocalAddr)[3] = args->dstLocalAddr;
    int strideArr = args->strideArr;
    int strideBuf = args->strideBuf;
    int nxGlobal = args->nxGlobal;
    int I0 = args->I0;
    int I2 = args->I2;
    int numLocalCopies = args->numLocalCopies;
    double fill = args->fill;
    int ioffset = args->ioffset;
    int joffset = args->joffset;
    int isign = args->isign;


    int total_copy = I2 * numLocalCopies;
    int local_block = (total_copy + 63) / 64;
    int local_start = CRTS_tid * local_block;
    int local_end = (CRTS_tid + 1) * local_block < total_copy?(CRTS_tid + 1) * local_block:total_copy;


    for(int t = local_start; t < local_end; t++)
    {
        int nmsg = t / I2;
        int k = t % I2;
        int srcBlock = srcLocalAddr[nmsg][2];
        if (srcBlock < 0) {
            int iSrc = srcLocalAddr[nmsg][0] - 1;
            int jSrc = srcLocalAddr[nmsg][1] - 1;

            int iDst = dstLocalAddr[nmsg][0] - 1;
            int jDst = dstLocalAddr[nmsg][1] - 1;
            // dstBlock = dstLocalAddr[nmsg][2];

            iSrc -= ioffset;
            jSrc -= joffset;
            if (iSrc == -1) {
                iSrc = nxGlobal - 1;
            }
            if (jSrc > -1) {
                if(isign == 1)
                {
                    if(jSrc == 0){
                        if(iSrc < nxGlobal - 1)
                        {
                            double x1 = bufTripole3D[k * strideBuf + iSrc];
                            double x2 =  bufTripole3D[k * strideBuf + nxGlobal - iSrc - 2];
                            if(isign == 1)
                            {
                                array[k*strideArr + jDst*I0 + iDst] = isign * 0.5 * (x1 + x2);
                            }
                            else
                            {
                                double xavg = (abs(x1) + abs(x2));
                                array[k*strideArr + jDst*I0 + iDst] = isign * 0.5 *  sign(xavg, x2);
                            }
                        }
                        else if(iSrc == nxGlobal - 1)
                        {
                            array[k*strideArr + jDst*I0 + iDst] = isign
                                                                  * bufTripole3D[k * strideBuf + jSrc * nxGlobal + iSrc];
                        }
                    }
                    else
                    {
                        array[k*strideArr + jDst*I0 + iDst] = isign
                                                              * bufTripole3D[k * strideBuf + jSrc * nxGlobal + iSrc];
                    }
                }
                else if(isign == -1)
                {
                    array[k*strideArr + jDst*I0 + iDst] = isign
                                                          * bufTripole3D[k * strideBuf + jSrc * nxGlobal + iSrc];
                }
            }
        }
    }

}
