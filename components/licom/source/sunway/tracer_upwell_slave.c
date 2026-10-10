#include <crts.h>
#include <slave.h>
#include <simd.h>
#include <math.h>
#include <stdio.h>
#include <float.h>
#include <string.h>
#include <stdbool.h>
#include "param_mod.h"
#include "constant_mod.h"
#include "swperf.h"

//#define UPWELLDEBUG 1
typedef struct tracer_upwell_cu_0_kernelPara{
	int cuGrid[3];
	int cuBlock[3];
        int rank;
        int pad;
	double *d_au0;
	double *d_aus;
	double *d_auw;
	double *d_ausw;
	double *d_stf;
	double *d_ohbu;
	double *d_wkd;
	double *d_uk;
	double *d_wkb;
	double *d_vk;
	double *d_work;
}tracer_upwell_cu_0_Para;

typedef struct tracer_upwellPara{
	int cuGrid[3];
	int cuBlock[3];
        int rank;
	int *d_kmt;
	double *d_dzp;
	double *d_stf;
	double *d_hts;
	double *d_htw;
	double *d_ohbt;
	double *d_tarea_r;
	double *d_uk;
	double *d_vit;
	double *d_vk;
	double *d_wka;
	double *d_work;
	double *d_ws;
}tracer_upwell_Para;

typedef struct tracer_upwell_cu_1_kernelPara{
    int cuGrid[3];
    int cuBlock[3];
    int rank;
    int padd;
    double *d_vit;
    double *d_dzp;
    double *d_wka;
    double *d_work;
}tracer_upwell_cu_1_Para;

typedef struct tracer_upwell_cu_2_kernelPara{
    int cuGrid[3];
    int cuBlock[3];
    int rank;
    int k;
    double *d_ws;
    double *d_vit;
    double *d_stf;
    double *d_dzp;
    double *d_ohbt;
    double *d_wka;
    double *d_work;
    double *d_work1;
    double *d_wsbak;
}tracer_upwell_cu_2_Para;

__thread_local unsigned int D_COUNT = 0;
__thread_local crts_rply_t dma_rply = 0;

void  tracer_upwell_cu_0_(tracer_upwell_cu_0_Para *Paras){

  //ATHREAD local variable definition
  int blockId_x,blockId_y,blockId_z;
  int threadId_x,threadId_y,threadId_z;
  double work_s = 0.0e0;

#ifdef UPWELLDEBUG
  unsigned long count,counte;
  unsigned long misscount,misscounte;

  penv_slave1_dcache_access_init();
  penv_slave2_dcache_miss_init();

  penv_slave1_dcache_access_count(&count);
  penv_slave2_dcache_miss_count(&misscount);
#endif

  //ATHREAD basic operation
    int i,j,k;
  int blockIndex = CRTS_tid;
  tracer_upwell_cu_0_Para* Para_s = (tracer_upwell_cu_0_Para*)ldm_malloc(sizeof(tracer_upwell_cu_0_Para));
  CRTS_dma_get(Para_s, Paras, sizeof(tracer_upwell_cu_0_Para));

//    int ret = CRTS_memcpy_sldm(array_2d,Para_s->d_stf,sizeof(double) * imt * jmt, MEM_TO_LDM);

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
    k = blockId_z * Para_s->cuBlock[2];

    //ATHREAD serial operation replace GPU parallel mutli-thread in single thread-block
    for(i = 0; i < threadSum; i++) {

      unsigned int tid=j *imt+i;
      work_s = 0.0e0;

      if(i<imt && i>=1 && j<jmt-1 && k< km) {

          //Para_s->d_work[tid] = Para_s->d_au0[tid] * Para_s->d_stf[tid] +
          work_s = p25 * Para_s->d_stf[tid] +
                  p25 * Para_s->d_stf[(j + 1) * imt + i] +
                  p25 * Para_s->d_stf[tid - 1] +
                  p25 * Para_s->d_stf[(j + 1) * imt + i - 1];

          Para_s->d_uk[k *jmt *imt+tid]=(1.0e0+work_s *Para_s->d_ohbu[tid]) *
                  Para_s->d_wkd[k *jmt *imt+tid];
          Para_s->d_vk[k *jmt *imt+tid]=(1.0e0+work_s *Para_s->d_ohbu[tid]) *
                  Para_s->d_wkb[k *jmt *imt+tid];
      }

    }

  }
#ifdef UPWELLDEBUG
    penv_slave1_dcache_access_count(&counte);
    penv_slave2_dcache_miss_count(&misscounte);
    if(Para_s->rank == 0)
        printf("[slave %d]:upwell_0 count:%d,counte:%d,misscount:%d,misscounte:%d,total count:%d, total misscount:%d, misscount rate:%f\n",
		blockIndex,count,counte,misscount,misscounte,counte-count, misscounte-misscount,
		 (float)(misscounte-misscount)/(float)(counte-count));
#endif
    ldm_free(Para_s, sizeof(tracer_upwell_cu_0_Para));
}

 void  tracer_upwell_cu_(tracer_upwell_Para *Paras){

  //ATHREAD local variable definition
  int blockId_x,blockId_y,blockId_z;
  int threadId_x,threadId_y,threadId_z;
#ifdef UPWELLDEBUG
  unsigned long count,counte;
  unsigned long misscount,misscounte;

  penv_slave1_dcache_access_init();
  penv_slave2_dcache_miss_init();

  penv_slave1_dcache_access_count(&count);
  penv_slave2_dcache_miss_count(&misscount);
#endif

  //ATHREAD basic operation
  int i,j,k;
  int blockIndex = CRTS_tid;
  tracer_upwell_Para* Para_s = (tracer_upwell_Para*)ldm_malloc(sizeof(tracer_upwell_Para));
  CRTS_dma_get(Para_s, Paras, sizeof(tracer_upwell_Para));

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
      k = blockId_z * Para_s->cuBlock[2];

    //ATHREAD serial operation replace GPU parallel mutli-thread in single thread-block
      for (i = 0; i < threadSum; i++) {
      unsigned int tid=j *imt+i;	
	
      Para_s->d_work[tid]=0.0e0;
      if(i<imt&&j<jmt&&k<km)
      {	
        if(i<imt-1&&j>=1)
        {
      	  Para_s->d_wka[k *jmt *imt+tid]=0.0e0;
          if(k<=Para_s->d_kmt[tid]-1)
          {
      	    Para_s->d_wka[k *jmt *imt+tid]=0.5 *((Para_s->d_uk[k *jmt *imt+tid+1]+
                      Para_s->d_uk[k *jmt *imt+tid-imt+1]) *Para_s->d_htw[tid+1]-
                      (Para_s->d_uk[k *jmt *imt+tid]+Para_s->d_uk[k *jmt *imt+tid-imt]) *
                      Para_s->d_htw[tid]+(Para_s->d_vk[k *jmt *imt+tid+1]+Para_s->d_vk[k *jmt *imt+tid]) *
                      Para_s->d_hts[tid]-(Para_s->d_vk[k *jmt *imt+tid-imt+1]+Para_s->d_vk[k *jmt *imt+tid-imt]) *
                      Para_s->d_hts[tid-imt]) *Para_s->d_tarea_r[tid];
          }

        }
      }

    }

  }
#ifdef UPWELLDEBUG
    penv_slave1_dcache_access_count(&counte);
    penv_slave2_dcache_miss_count(&misscounte);
    if(Para_s->rank == 0)
        printf("[slave %d]:upwell count:%d,counte:%d,misscount:%d,misscounte:%d,total count:%d, total misscount:%d, misscount rate:%f\n",
		blockIndex,count,counte,misscount,misscounte,counte-count, misscounte-misscount,
		 (float)(misscounte-misscount)/(float)(counte-count));
#endif
    ldm_free(Para_s, sizeof(tracer_upwell_Para));
}

void  tracer_upwell_cu_1_(tracer_upwell_cu_1_Para *Paras){

    //ATHREAD local variable definition
    int swGrid[3];
    int swBlock[3];
    int blockId_x,blockId_y,blockId_z;
    int threadId_x,threadId_y,threadId_z;
    double d_work[imt];
    double d_wka[imt];

#ifdef UPWELLDEBUG
    unsigned long count,counte;
    unsigned long misscount,misscounte;

    penv_slave1_dcache_access_init();
    penv_slave2_dcache_miss_init();

    penv_slave1_dcache_access_count(&count);
    penv_slave2_dcache_miss_count(&misscount);
#endif

    //ATHREAD basic operation
    int blockIndex = CRTS_tid;
    tracer_upwell_cu_1_Para* Para_s = (tracer_upwell_cu_1_Para*)ldm_malloc(sizeof(tracer_upwell_cu_1_Para));
    CRTS_dma_get(Para_s, Paras, sizeof(tracer_upwell_cu_1_Para));

    for(int k = 0; k < km; k++) {
        for (int j = blockIndex; j < jmt - 1; j += 64) {

            CRTS_dma_get(d_work, &Para_s->d_work[j * imt], sizeof(double) * imt);
            CRTS_dma_get(d_wka, &Para_s->d_wka[k * imt * jmt + j * imt], sizeof(double) * imt);

            for (int i = 1; i < imt - 1; i++) {
                if (1 <= j && j < (jmt - 1) && 1 <= i && i < (imt - 1)){
                    //WORK (I,J,IBLOCK)= WORK(I,J,IBLOCK) - DZP (K)* WKA (I,J,K,IBLOCK)* VIT (I,J,K,IBLOCK)
                    d_work[i] = d_work[i] - Para_s->d_dzp[k] *d_wka[i] *Para_s->d_vit[k *jmt *imt+j * imt + i];
                }
            }

            CRTS_dma_put(&Para_s->d_work[j*imt], d_work, sizeof(double) * imt);
        }
    }

#ifdef UPWELLDEBUG
    penv_slave1_dcache_access_count(&counte);
    penv_slave2_dcache_miss_count(&misscounte);
    if(Para_s->rank == 0)
        printf("[slave %d]:upwell_1 count:%d,counte:%d,misscount:%d,misscounte:%d,total count:%d, total misscount:%d, misscount rate:%f\n",
		blockIndex,count,counte,misscount,misscounte,counte-count, misscounte-misscount,
		 (float)(misscounte-misscount)/(float)(counte-count));
#endif
    ldm_free(Para_s, sizeof(tracer_upwell_cu_1_Para));
}


void  tracer_upwell_cu_2_new_(tracer_upwell_cu_2_Para *Paras){

    double work1_s = 0.0e0;
    double d_wsbak = 0.0e0;
    //ATHREAD basic operation
    int blockIndex = CRTS_tid;
    tracer_upwell_cu_2_Para* Para_s = (tracer_upwell_cu_2_Para*)ldm_malloc(sizeof(tracer_upwell_cu_2_Para));
    CRTS_dma_get(Para_s, Paras, sizeof(tracer_upwell_cu_2_Para));

    int tile  = imt ;
    int avaiLdm = get_allocatable_size();
    int varNum = 4; // variable of dma stride
    int tileBlk = 1; //blocks number in i direction
    int tileLimit = avaiLdm/(varNum*km*sizeof(double));
    int remTile = 0;
    if(tile > tileLimit) {
        remTile = tile % tileLimit;
        if(remTile == 0)
            tileBlk = tile / tileLimit;
        else
            tileBlk = tile / tileLimit + 1;
        tile = tileLimit;
    }

    int start_blk;
    int tile_size  = tile * sizeof(double);
    int tileSum   = tile * km * sizeof(double);
    int size_stride = (imt * jmt - tile) * sizeof(double);
    int tileV;
    int tileSumV;
    int size_strideV;

    // target variable for k direction inconsistent access
    double* ldm_wsbak  = (double*) ldm_malloc (tileSum);
    double* ldm_wka  = (double*) ldm_malloc (tileSum);
    double* ldm_vit  = (double*) ldm_malloc (tileSum);
    double* ldm_ws   = (double*) ldm_malloc (tileSum);

    double work_s[tile];
    // each cpe is responsible for i*k data
    for(int j = blockIndex; j < jmt; j += 64){


        // blocks in i direction
        for(int idx=0; idx < tileBlk; idx ++) {
            start_blk = j * imt + idx * tile;
            //distinguish tile and remainder for stride dma
            if((remTile != 0) && (idx == tileBlk -1)) {
                tileV = remTile;
                tileSumV   = remTile * km * sizeof(double);
                tile_size  = remTile * sizeof(double);
                size_stride = (imt * jmt - remTile) * sizeof(double);
            }
            else {
                tileV = tile;
                tileSumV   = tile * km * sizeof(double);
                tile_size  = tile * sizeof(double);
                size_stride = (imt * jmt - tile) * sizeof(double);
            }

            athread_dma_iget_stride(ldm_wsbak, &Para_s->d_wsbak[start_blk],
                                    tileSumV, tile_size, size_stride, &dma_rply);
            D_COUNT++;
            athread_dma_wait_value(&dma_rply, D_COUNT);

            athread_dma_iget_stride(ldm_vit, &Para_s->d_vit[start_blk],
                                    tileSumV, tile_size, size_stride, &dma_rply);
            D_COUNT++;
            athread_dma_wait_value(&dma_rply, D_COUNT);

            athread_dma_iget_stride(ldm_wka, &Para_s->d_wka[start_blk],
                                    tileSumV, tile_size, size_stride, &dma_rply);
            D_COUNT++;
            athread_dma_wait_value(&dma_rply, D_COUNT);


            for (int i = 0; i < tileV; ++i) {

                int index = start_blk + i;
                int ii = idx * tile + i;


                if (1 <= j && j < (jmt - 1) && 1 <= ii && ii < (imt - 1)) {

                    work_s[i] = 0.0;
                    //WORK (I,J,IBLOCK)= WORK(I,J,IBLOCK) - DZP (K)* WKA (I,J,K,IBLOCK)* VIT (I,J,K,IBLOCK)
                    for(int k = 0; k < km; k++)
                        work_s[i] = work_s[i] - Para_s->d_dzp[k] * ldm_wka[k * tileV + i] * ldm_vit[k * tileV + i];

//                    WS (I,J,K,IBLOCK)= VIT (I,J,K,IBLOCK)* (WS (I,J,K -1,IBLOCK) + &
//                            DZP (K -1)* (WORK (I,J,IBLOCK)* OHBT (I,J,IBLOCK) + WKA (I,J,K -1,IBLOCK)))

                    work1_s = 1.0e0 / (1.0e0 + Para_s->d_stf[index] * Para_s->d_ohbt[index]);
                    for(int k=1; k<km; k++) {

                        ldm_wsbak[k * tileV + i] = ldm_vit[k * tileV + i] *
                                     (ldm_wsbak[(k-1) * tileV + i] + Para_s->d_dzp[k - 1] *
                                     (work_s[i] * Para_s->d_ohbt[index] + ldm_wka[(k-1) * tileV + i]));

                        ldm_ws[k * tileV + i] = ldm_wsbak[k * tileV + i] * work1_s;
                    }
                    //WS (I,J,K,IBLOCK)= WS (I,J,K,IBLOCK)* WORK (I,J,IBLOCK)
                }
            }

            athread_dma_iput_stride (&Para_s->d_ws[start_blk],
                                     ldm_ws, tileSumV, tile_size, size_stride,&dma_rply);
            D_COUNT++;
            athread_dma_wait_value(&dma_rply, D_COUNT);

            athread_dma_iput_stride (&Para_s->d_wsbak[start_blk],
                                     ldm_wsbak, tileSumV, tile_size, size_stride,&dma_rply);
            D_COUNT++;
            athread_dma_wait_value(&dma_rply, D_COUNT);

        }
    }

    ldm_free (ldm_wsbak, tileSum);
    ldm_free (ldm_wka, tileSum);
    ldm_free (ldm_vit, tileSum);
    ldm_free (ldm_ws, tileSum);
    ldm_free(Para_s, sizeof(tracer_upwell_cu_2_Para));
}


void  tracer_upwell_cu_2_(tracer_upwell_cu_2_Para *Paras){

    //ATHREAD local variable definition
    int swGrid[3];
    int swBlock[3];
    int blockId_x,blockId_y,blockId_z;
    int threadId_x,threadId_y,threadId_z;
    double work_s = 0.0e0;
    double work1_s = 0.0e0;
    double wsk = 0.0e0;
    double wsk_1 = 0.0e0;
    double d_ws[imt];
    double d_wsmk1[imt];
    double d_wsbak[imt];
    double d_wka[imt];

#ifdef UPWELLDEBUG
    unsigned long count,counte;
    unsigned long misscount,misscounte;

    penv_slave1_dcache_access_init();
    penv_slave2_dcache_miss_init();

    penv_slave1_dcache_access_count(&count);
    penv_slave2_dcache_miss_count(&misscount);
#endif

    //ATHREAD basic operation
    int blockIndex = CRTS_tid;
    tracer_upwell_cu_2_Para* Para_s = (tracer_upwell_cu_2_Para*)ldm_malloc(sizeof(tracer_upwell_cu_2_Para));
    CRTS_dma_get(Para_s, Paras, sizeof(tracer_upwell_cu_2_Para));

    for(int k = 1; k < km; k++) {
        for (int j = blockIndex; j < jmt - 1; j += 64) {

            //CRTS_dma_get(d_ws, &Para_s->d_ws[k * imt * jmt + j*imt], sizeof(double) * imt);
            CRTS_dma_get(d_wsmk1, &Para_s->d_wsbak[(k-1) * imt * jmt + j*imt], sizeof(double) * imt);
            CRTS_dma_get(d_wka, &Para_s->d_wka[(k-1) * imt * jmt + j*imt], sizeof(double) * imt);

            for (int i = 1; i < imt - 1; i++) {
                if (1 <= j && j < (jmt - 1) && 1 <= i && i < (imt - 1)) {
//                    WS (I,J,K,IBLOCK)= VIT (I,J,K,IBLOCK)* (WS (I,J,K -1,IBLOCK) + &
//                            DZP (K -1)* (WORK (I,J,IBLOCK)* OHBT (I,J,IBLOCK) + WKA (I,J,K -1,IBLOCK)))
                    d_wsbak[i] = Para_s->d_vit[k * imt * jmt + j * imt + i] *
                              (d_wsmk1[i] + Para_s->d_dzp[k - 1] *
                              (Para_s->d_work[j * imt + i] * Para_s->d_ohbt[j * imt + i] + d_wka[i]));
                    work1_s = 1.0e0 / (1.0e0 + Para_s->d_stf[j * imt + i] * Para_s->d_ohbt[j * imt + i]);
                    //WS (I,J,K,IBLOCK)= WS (I,J,K,IBLOCK)* WORK (I,J,IBLOCK)
                    d_ws[i] = d_wsbak[i] * work1_s;
                }
            }
            CRTS_dma_put(&Para_s->d_ws[k * imt * jmt + j*imt], d_ws, sizeof(double) * imt);
            CRTS_dma_put(&Para_s->d_wsbak[k * imt * jmt + j*imt], d_wsbak, sizeof(double) * imt);

        }
    }

#ifdef UPWELLDEBUG
    penv_slave1_dcache_access_count(&counte);
    penv_slave2_dcache_miss_count(&misscounte);
    if(Para_s->rank == 0)
        printf("[slave %d]:upwell_2 count:%d,counte:%d,misscount:%d,misscounte:%d,total count:%d, total misscount:%d, misscount rate:%f\n",
		blockIndex,count,counte,misscount,misscounte,counte-count, misscounte-misscount,
		 (float)(misscounte-misscount)/(float)(counte-count));
#endif
    ldm_free(Para_s, sizeof(tracer_upwell_cu_2_Para));
}
