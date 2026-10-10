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

typedef struct bclinc11Para{
	int cuBlock[3];
	int cuGrid[3];
	int d_isc;
	int rank;
	double d_od0;
	double d_onbb;
	double *d_ohbt;
	double *d_vit;
	double *d_dzp;
	double *d_h0bf;
	double *d_h0bl;
	double *d_gg;
	double *d_psa;
	double *d_wka;
	double *d_wka1;
	double *d_wka2;
	double *d_zkt;
    double *d_work;
}bclinc11_Para;

typedef struct bclinc6Para{
	int cuBlock[3];
	int cuGrid[3];	
	int d_isc;
	double *d_viv;
	double *d_epea;
	double *d_epeb;
	double *d_epla;
	double *d_eplb;
	double *d_up;
	double *d_vp;
	double *d_dlu;
	double *d_dlv;
	double *d_gg;
	double *d_fcor;
	double *d_wka1;
	double *d_wka2;
	double *d_dxur;
	double *d_dyur;
	int    *d_kmu;
}bclinc6_Para;

//dynamic atomic definition
#define CP_CORE      64
#define RCV_CORE     0
#define DATAPAKAGE   1

//ATHREAD global definition 
__thread_local int warp_size = 32; 
__thread_local int blockIndex;
//__thread_local Para* Para_s;

void  bclinc11_cu_slave_kernel_(bclinc11_Para *Paras){

    double work = 0.0;
//    double wkk_km1[imt * jmt];
    double wkk_k = 0.0;
    double *gg_km1;
    double *vit_km1;
    double d_wka1[imt];
    double d_wka[imt];
    double aa = 0.0;
    double d_h0bf = 0.0;

    //ATHREAD basic operation
    blockIndex = CRTS_tid;
    bclinc11_Para *Para_s = (bclinc11_Para*)ldm_malloc(sizeof(bclinc11_Para));
    CRTS_dma_iget(Para_s, Paras, sizeof(bclinc11_Para), &dma_rply);
    D_COUNT++;
    CRTS_dma_wait_value(&dma_rply, D_COUNT);

    double d_od0 = Para_s->d_od0;
    double d_onbb = Para_s->d_onbb;
    int d_isc = Para_s->d_isc;

    int quotient   = jmt / 64;
    double wkk_km1[imt * (quotient+1)];
    int jj;

    if(d_isc != 0)
        aa = 0.5;


    for(int k = 1; k < km+1; k++) {
        jj = 0;
        for (int j = blockIndex; j < jmt; j += 64) {

            gg_km1 = &Para_s->d_gg[(k-1) * imt * jmt + j * imt];
            CRTS_dma_get(d_wka1, &Para_s->d_wka1[(k-1) * imt * jmt + j * imt], sizeof(double) * imt);
            CRTS_dma_get(d_wka, &Para_s->d_wka[(k-1) * imt * jmt + j * imt], sizeof(double) * imt);
            vit_km1 = &Para_s->d_vit[(k-1) * imt * jmt + j * imt];

            for (int i = 0; i < imt; i++) {
                if( k == 1){
                    //H0BF (I,J,IBLOCK)= H0BF (I,J,IBLOCK)* ONBB
                    d_h0bf = Para_s->d_h0bf[j *imt+i] * d_onbb;
                    // WORK (I,J,IBLOCK)= AA * H0BF (I,J,IBLOCK) + (1.0- AA)* H0BL (I,J,IBLOCK)
                    work = aa *d_h0bf+(1.0-aa) *Para_s->d_h0bl[j *imt+i];
                    Para_s->d_work[j *imt+i] = work;
                    Para_s->d_h0bf[j *imt+i] = d_h0bf;
                    //wktmp(i,j,1,IBLOCK)= (PSA (I,J,IBLOCK)* OD0+ WORK (I,J,IBLOCK)* G)* VIT (I,J,1,IBLOCK)
                    wkk_km1[jj*imt + i] = (Para_s->d_psa[j *imt+i] * d_od0 + work * 9.806e0) *
                            vit_km1[i];
                }

                //wktmp(i,j,k,IBLOCK)= wktmp(i,j,k-1,IBLOCK) - GG (I,J,K-1,IBLOCK)* DZP (K-1)* VIT (I,J,K-1,IBLOCK)
                wkk_k = wkk_km1[jj*imt + i] - gg_km1[i] * Para_s->d_dzp[k - 1] * vit_km1[i];
                //WKA (I,J,K-1,IBLOCK)= P5* (wktmp(i,j,k-1,IBLOCK) + wktmp(i,j,k,IBLOCK))
                d_wka1[i] = 0.5 * (wkk_km1[jj*imt + i] + wkk_k);
                wkk_km1[jj*imt + i] = wkk_k;
                //WKA (I,J,K,IBLOCK)= (1.0+ OHBT (I,J,IBLOCK)* ZKT (K))* WORK (I,J,IBLOCK)      &
                //                            * VIT (I,J,K,IBLOCK)
                d_wka[i] = (1.0 + Para_s->d_ohbt[j*imt + i] * Para_s->d_zkt[k-1]) *
                        Para_s->d_work[j *imt+i] * vit_km1[i];
            }
            CRTS_dma_put(&Para_s->d_wka1[(k-1) * imt * jmt + j * imt], d_wka1, sizeof(double) * imt);
            CRTS_dma_put(&Para_s->d_wka[(k-1) * imt * jmt + j * imt], d_wka, sizeof(double) * imt);
            jj++;
        }

    }

    ldm_free(Para_s, sizeof(bclinc11_Para));
}

void  bclinc6_cu_slave_kernel_(bclinc6_Para *Paras){

    int blockId_x,blockId_y,blockId_z;
    int i,j,k;
  blockIndex = CRTS_tid;
  bclinc6_Para *Para_s = (bclinc6_Para*)ldm_malloc(sizeof(bclinc6_Para));
  CRTS_dma_iget(Para_s, Paras, sizeof(bclinc6_Para), &dma_rply);
  D_COUNT++;
  CRTS_dma_wait_value(&dma_rply, D_COUNT);

  int d_isc = Para_s->d_isc;
  //double* d_wka2 = Para_s->d_wka2;
  double* d_epla = Para_s->d_epla;
  double* d_eplb = Para_s->d_eplb;
  double* d_epea = Para_s->d_epea;
  double* d_epeb = Para_s->d_epeb;
//  double d_gg[imt],d_ggjp1[imt],d_wka1[imt],d_wka1jp1[imt],d_wka2[imt],d_wka2jp1[imt],d_vp[imt],d_up[imt],d_viv[imt],d_dlu[imt],d_dlv[imt];
  double *d_gg,*d_ggjp1,*d_wka1,*d_wka1jp1,*d_wka2,*d_wka2jp1,*d_vp,*d_up,*d_viv,*d_dlu,*d_dlv;
  double dlu_arr[imt],dlv_arr[imt];

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

     if(j<jmt-1)
     {
	      d_gg      = &Para_s->d_gg[k * imt * jmt + j*imt];
	      d_ggjp1   = &Para_s->d_gg[k * imt * jmt + (j+1)*imt];
	      d_wka1    = &Para_s->d_wka1[k * imt * jmt + j*imt];
	      d_wka1jp1 = &Para_s->d_wka1[k * imt * jmt + (j+1)*imt];
	      d_wka2    = &Para_s->d_wka2[k * imt * jmt + j*imt];
        d_wka2jp1 = &Para_s->d_wka2[k * imt * jmt + (j+1)*imt];
	      d_vp      = &Para_s->d_vp[k * imt * jmt + j*imt];
	      d_up      = &Para_s->d_up[k * imt * jmt + j*imt];
	      d_viv     = &Para_s->d_viv[k * imt * jmt + j*imt];
	      d_dlu     = &Para_s->d_dlu[k * imt * jmt + j*imt];
	      d_dlv     = &Para_s->d_dlv[k * imt * jmt + j*imt];

	//CRTS_dma_get(d_hue, &Para_s->d_hue[j*imt], sizeof(double) * imt);

	for(i=1;i<imt;i++)
        {

      int kmb=Para_s->d_kmu[j *imt+i]-1;	
	
      double gradx,grady,ggu,wk1,wk2,dlu,dlv;	
      dlu = d_dlu[i] - Para_s->d_fcor[j *imt+i] *d_vp[i];
      dlv = d_dlv[i] + Para_s->d_fcor[j *imt+i] *d_up[i];
	
      if(kmb>=k)	
      {	
      	ggu=0.25 *(d_gg[i]+d_gg[i-1]+d_ggjp1[i]+d_ggjp1[i-1]);	

      	gradx=Para_s->d_dxur[j *imt+i] *0.5 *(d_wka1jp1[i]-d_wka1[i-1]-d_wka1jp1[i-1]+d_wka1[i]);	

      	grady=Para_s->d_dyur[j *imt+i] *0.5 *(d_wka1jp1[i]-d_wka1[i-1]+d_wka1jp1[i-1]-d_wka1[i]);	
	}
      dlu -= gradx;
      dlv -= grady;
      if(kmb>=k)
      {
	   gradx=Para_s->d_dxur[j *imt+i] *0.5 *(d_wka2jp1[i]-d_wka2[i-1]-d_wka2jp1[i-1]+d_wka2[i]);  
       	   grady=Para_s->d_dyur[j *imt+i] *0.5 *(d_wka2jp1[i]-d_wka2[i-1]+d_wka2jp1[i-1]-d_wka2[i]);
      }	
	dlu += ggu*gradx;
	dlv += ggu*grady;
      if(d_isc==0)	
      {	
      	wk1=(d_epea[j *imt+i] *dlv +d_epeb[j *imt+i] *dlu);	
      	wk2=(d_epea[j *imt+i] *dlu -d_epeb[j *imt+i] *dlv);	
      }	
      else	{
	      	wk1=(d_epla[j *imt+i] *dlv +d_eplb[j *imt+i] *dlu);	
      	wk2=(d_epla[j *imt+i] *dlu-d_eplb[j *imt+i] *dlv);	
      }	
      	dlv_arr[i]=wk1 *d_viv[i];	
      	dlu_arr[i]=wk2 *d_viv[i];	
      }
      CRTS_dma_put(&Para_s->d_dlv[k * imt * jmt + j*imt],dlv_arr, sizeof(double) * imt);
      CRTS_dma_put(&Para_s->d_dlu[k * imt * jmt + j*imt],dlu_arr, sizeof(double) * imt);
}
}	


    ldm_free(Para_s, sizeof(bclinc6_Para));
}


