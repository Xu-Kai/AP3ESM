
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
#include "def-undef.h"	

//ATHREAD DMA transfer definition 
__thread_local unsigned int D_COUNT = 0;
__thread_local crts_rply_t dma_rply = 0;
__thread_local crts_rply_t l_rply = 0;
__thread_local crts_rply_t r_rply = 0;

//__userDef__ paraVarAttr(double,d_vit,sizeof(double)*jmt*imt*km,readonly,no,no,no);	
//__userDef__ paraVarAttr(double,d_h0p,sizeof(double)*jmt*imt,readonly,no,no,no);	
//__userDef__ paraVarAttr(double,d_wka,sizeof(double)*km*jmt*imt,readonly,no,no,no);	
//__userDef__ paraVarAttr(double,d_work,sizeof(double)*jmt*imt,outonly,yes,no,no);	
//__userDef__ paraVarAttr(double,d_htw,sizeof(double)*ny_block*nx_block,readonly,no,no,no);	
//__userDef__ paraVarAttr(double,d_hts,sizeof(double)*ny_block*nx_block,readonly,no,no,no);	
//__userDef__ paraVarAttr(int,d_kmt,sizeof(int)*ny_block*nx_block,readonly,no,no,no);	
//__userDef__ paraVarAttr(double,d_tarea_r,sizeof(double)*ny_block*nx_block,readonly,no,no,no);	
//__userDef__ paraVarAttr(double,d_ah,sizeof(double),readonly,no,no,no);	
//__userDef__ paraVarAttr(int,d_kmtn,sizeof(int)*ny_block*nx_block,readonly,no,no,no);	
//__userDef__ paraVarAttr(int,d_kmts,sizeof(int)*ny_block*nx_block,readonly,no,no,no);	
//__userDef__ paraVarAttr(int,d_kmte,sizeof(int)*ny_block*nx_block,readonly,no,no,no);	
//__userDef__ paraVarAttr(int,d_kmtw,sizeof(int)*ny_block*nx_block,readonly,no,no,no);	
//__userDef__ paraVarAttr(double,d_dtn,sizeof(double)*ny_block*nx_block,readonly,no,no,no);	
//__userDef__ paraVarAttr(double,d_dts,sizeof(double)*ny_block*nx_block,readonly,no,no,no);	
//__userDef__ paraVarAttr(double,d_dte,sizeof(double)*ny_block*nx_block,readonly,no,no,no);	
//__userDef__ paraVarAttr(double,d_dtw,sizeof(double)*ny_block*nx_block,readonly,no,no,no);	
//__userDef__ paraVarAttr(int,d_jb,sizeof(int),readonly,no,no,no);	
//__userDef__ paraVarAttr(int,d_je,sizeof(int),readonly,no,no,no);	
//__userDef__ paraVarAttr(int,d_ib,sizeof(int),readonly,no,no,no);	
//__userDef__ paraVarAttr(int,d_ie,sizeof(int),readonly,no,no,no);	
//__userDef__ paraVarAttr(int,d_nc,sizeof(int),readonly,no,no,no);	
//__userDef__ paraVarAttr(int,d_jmt,sizeof(int),readonly,no,no,no);	
//__userDef__ paraVarAttr(int,d_imt,sizeof(int),readonly,no,no,no);	
//__userDef__ paraVarAttr(int,d_km,sizeof(int),readonly,no,no,no);	
//__userDef__ paraVarAttr(double,d_d2tk,sizeof(double)*jmt*imt*km,readonly,no,no,no);	
typedef struct barotr3_cu1_kernelPara{
	int cuBlock[3];
	int cuGrid[3];
	int d_jb;
	int d_je;
	int d_ib;
	int d_ie;
	int d_nc;
	int d_imt;
	int d_jmt;
	int d_km;
	int rank;
	int padd;
	double d_ah;
	int *d_kmt;
	int *d_kmtn;
	int *d_kmts;
	int *d_kmte;
	int *d_kmtw;
	double *d_vit;
	double *d_h0p;
	double *d_wka;
	double *d_work;
	double *d_htw;
	double *d_hts;
	double *d_tarea_r;
	double *d_dtn;
	double *d_dts;
	double *d_dte;
	double *d_dtw;
	double *d_d2tk;
	double *d_factor_vis_tsz0;
}barotr3_cu1_Para;


//dynamic atomic definition
#define CP_CORE      64
#define RCV_CORE     0
#define DATAPAKAGE   1

//ATHREAD global definition 
__thread_local int warp_size = 32; 
__thread_local int blockIndex;

void  barotr3_cu1_(barotr3_cu1_Para *Paras){

    double div_out=0.0e0,hdtk=0.0e0;
    double cc=0.0e0;
    double cn=0.0e0;
    double ce=0.0e0;
    double cs=0.0e0;
    double cw=0.0e0;
//    double d_wka0[imt];
//    double d_wka0mj1[imt];
//    double d_htw[imt];
//    double d_wka1[imt];
//    double d_wka1mj1[imt];
//    double d_hts[imt];
//    double d_htsmj1[imt];
//    double d_tarea_r[imt];
//    double d_d2tk[imt];
//    double d_d2tkpj1[imt];
//    double d_d2tkmj1[imt];
//    double d_work[imt];
    double *d_wka0,*d_wka0mj1, *d_htw, *d_wka1, *d_wka1mj1;
    double *d_hts, *d_htsmj1, *d_tarea_r, *d_d2tk, *d_d2tkpj1, *d_d2tkmj1, *d_work;
  //ATHREAD basic operation
    int blockId_x,blockId_y,blockId_z;
    int i,j,k;
  blockIndex = CRTS_tid;
  barotr3_cu1_Para *Para_s = (barotr3_cu1_Para*)ldm_malloc(sizeof(barotr3_cu1_Para));
  CRTS_dma_iget(Para_s, Paras, sizeof(barotr3_cu1_Para), &dma_rply);
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

//        CRTS_dma_get(d_wka0, &Para_s->d_wka[j * imt], sizeof(double) * imt);
//        CRTS_dma_get(d_wka0mj1, &Para_s->d_wka[(j - 1) * imt], sizeof(double) * imt);
//        CRTS_dma_get(d_htw, &Para_s->d_htw[j * imt], sizeof(double) * imt);
//        CRTS_dma_get(d_wka1, &Para_s->d_wka[1 *imt *jmt+j * imt], sizeof(double) * imt);
//        CRTS_dma_get(d_wka1mj1, &Para_s->d_wka[1 *imt *jmt+(j-1) * imt], sizeof(double) * imt);
//        CRTS_dma_get(d_hts, &Para_s->d_hts[j * imt], sizeof(double) * imt);
//        CRTS_dma_get(d_htsmj1, &Para_s->d_hts[(j - 1) * imt], sizeof(double) * imt);
//        CRTS_dma_get(d_tarea_r, &Para_s->d_tarea_r[j * imt], sizeof(double) * imt);
//        CRTS_dma_get(d_d2tk, &Para_s->d_d2tk[j * imt], sizeof(double) * imt);
//        CRTS_dma_get(d_d2tkpj1, &Para_s->d_d2tk[(j+1) * imt], sizeof(double) * imt);
//        CRTS_dma_get(d_d2tkmj1, &Para_s->d_d2tk[(j-1) * imt], sizeof(double) * imt);
//        CRTS_dma_get(d_work, &Para_s->d_work[j * imt], sizeof(double) * imt);

        d_wka0    = &Para_s->d_wka[j * imt];
        d_wka0mj1 = &Para_s->d_wka[(j - 1) * imt];
        d_htw     = &Para_s->d_htw[j * imt];
        d_wka1    = &Para_s->d_wka[1 *imt *jmt+j * imt];
        d_wka1mj1 = &Para_s->d_wka[1 *imt *jmt+(j-1) * imt];
        d_hts     = &Para_s->d_hts[j * imt];
        d_htsmj1  = &Para_s->d_hts[(j - 1) * imt];
        d_tarea_r = &Para_s->d_tarea_r[j * imt];
        d_d2tk    = &Para_s->d_d2tk[j * imt];
        d_d2tkpj1 = &Para_s->d_d2tk[(j+1) * imt];
        d_d2tkmj1 = &Para_s->d_d2tk[(j-1) * imt];
        d_work    = &Para_s->d_work[j * imt];

        for (int i = 0; i < imt; i++) {
            if(i<imt-1&&j>=1&&j<jmt) {
                div_out=0.0e0;
                hdtk=0.0e0;

                if(0<=Para_s->d_kmt[j *imt+i]-1)
                {
                    div_out=0.5 *((d_wka0[i+1]+d_wka0mj1[i+1]) * d_htw[i+1]-(d_wka0[i]+
                            d_wka0mj1[i]) *d_htw[i]+(d_wka1[i+1]+d_wka1[i]) *d_hts[i]-
                                  (d_wka1mj1[i+1]+d_wka1mj1[i]) *d_htsmj1[i]) *d_tarea_r[i];
                }
                if(Para_s->d_nc %4==0)
                {
                    cc=0.0e0;
                    cn=0.0e0;
                    cs=0.0e0;
                    ce=0.0e0;
                    cw=0.0e0;
#ifdef  BIHAR
                    if(j>=Para_s->d_jb-1&&j<Para_s->d_je&&i>=Para_s->d_ib-1&&i<Para_s->d_ie)
                    {
                        if(0<=Para_s->d_kmt[j *imt+i]-1)
                        {
                            double c;

                            if(0<=Para_s->d_kmtn[j *imt+i]-1)
                            {
                                cn=Para_s->d_dtn[j *imt+i];
                                //hdtk+=c *Para_s->d_d2tk[(j-1) *imt+i];
                                //cc-=c;
                            }
                            if(0<=Para_s->d_kmts[j *imt+i]-1)
                            {
                                cs=Para_s->d_dts[j *imt+i];
                                //hdtk+=c *Para_s->d_d2tk[(j+1) *imt+i];
                                //cc-=c;
                            }
                            if(0<=Para_s->d_kmte[j *imt+i]-1)
                            {
                                ce=Para_s->d_dte[j *imt+i];
                                //hdtk+=c *Para_s->d_d2tk[j *imt+i+1];
                                //cc-=c;
                            }
                            if(0<=Para_s->d_kmtw[j *imt+i]-1)
                            {
                                cw=Para_s->d_dtw[j *imt+i];
                                //hdtk+=c *Para_s->d_d2tk[j *imt+i-1];
                                //cc-=c;
                            }
                            cc = -(cn + cs + ce +cw);
                        }
                        hdtk=Para_s->d_ah * Para_s->d_factor_vis_tsz0[j*imt+i] *
                                    (cc *d_d2tk[i]+cn *d_d2tkmj1[i] +
                                     cs *d_d2tkpj1[i] + ce *d_d2tk[i+1] +
                                     cw *d_d2tk[i-1]);
                    }
                    else{
                        hdtk=0.0e0;
                    }
#else
                if (0 <= Para_s->d_kmt[j * imt + i] - 1) { 
                    double c;

                    if (0 <= Para_s->d_kmtn[j * imt + i] - 1) { 
                        cn = Para_s->d_dtn[j * imt + i];
                    }    


                    if (0 <= Para_s->d_kmts[j * imt + i] - 1) { 
                        cs = Para_s->d_dts[j * imt + i];
                    }    

                    if (0 <= Para_s->d_kmte[j * imt + i] - 1) { 
                        ce = Para_s->d_dte[j * imt + i];
                    }    

                    if (0 <= Para_s->d_kmtw[j * imt + i] - 1) { 
                        cw = Para_s->d_dtw[j * imt + i];
                    }
		    cc = -(cn + cs + ce +cw);    
                }    


                if (j >= Para_s->d_jb - 1 && j < Para_s->d_je && i >= Para_s->d_ib - 1 && i < Para_s->d_ie) {
                    hdtk = Para_s->d_ah * (cc * Para_s->d_h0p[j * imt + i] + cn*Para_s->d_h0p[(j - 1) * imt + i] +
				    cs*Para_s->d_h0p[(j + 1) * imt + i] +ce*Para_s->d_h0p[j * imt + i + 1] +
				    cw*Para_s->d_h0p[j * imt + i - 1]);
                } else {
                    hdtk = 0.0e0;
                }    
#endif
                }
                else
                    hdtk=0.0e0;
                //d_work_s[j *imt+i]=(hdtk *1.0e0-div_out) *Para_s->d_vit[j *imt+i];
                d_work[i] = Para_s->d_vit[j *imt+i] * (hdtk * 1.0e0 - div_out);

            }
        }

//        CRTS_dma_put(&Para_s->d_work[j * imt], d_work, sizeof(double) * imt);
    }

    ldm_free(Para_s, sizeof(barotr3_cu1_Para));
}
