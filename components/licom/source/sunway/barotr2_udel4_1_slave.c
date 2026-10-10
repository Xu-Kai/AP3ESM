
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

typedef struct barotr2_udel4_1_kernelPara{
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
	double *d_d2uk;
	double *d_d2vk; 
	double *d_am_factor;
	double *d_uarea_r;
	double *d_factor;
    double *d_fac_vis_uv;
	double d_dts;
	double d_value;
    double d_value2;
}barotr2_udel4_1_Para;


//dynamic atomic definition
#define CP_CORE      64
#define RCV_CORE     0
#define DATAPAKAGE   1

//ATHREAD global definition 
__thread_local int warp_size = 32; 
__thread_local int blockIndex;

void barotr2_udel4_1_(barotr2_udel4_1_Para *Paras){

//    double d_d2uk[3 * imt];
//    double d_d2vk[3 * imt];
    double *d_d2uk;
    double *d_d2vk;
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
#ifdef TIDEFORCE
    double *d_h0_tide;
    double *d_h0_tidepj1;
#endif
    double *d_h0;
    double *d_h0pj1;
    double *d_wka0;
    double *d_wka1;
    double *d_wka2;
    double *d_wka3;
    double *d_wka4;
    double *d_wka5;
    double d_work[imt];
    double *d_dxur;
    double *d_dyur;
    double *d_ebea;
    double *d_ebeb;
    double *d_fcor;
    double *d_dlub;
    double *d_dlvb;
    double *d_vbp;
    double *d_ubp;
    double *d_viv;
    double *d_pax;
    double *d_pay;
    double *d_pxb;
    double *d_pyb;
    double *d_whx;
    double *d_why;
    double hduk = 0.0e0;
    double hdvk = 0.0e0;
    double gradx=0.0e0,grady=0.0e0;
	double am_f,t_area,scaling,u_area2,temp;


    //ATHREAD basic operation
    int i,j,k;
    int blockId_x,blockId_y,blockId_z;
    blockIndex = CRTS_tid;
    barotr2_udel4_1_Para *Para_s = (barotr2_udel4_1_Para*)ldm_malloc(sizeof(barotr2_udel4_1_Para));
    CRTS_dma_iget(Para_s, Paras, sizeof(barotr2_udel4_1_Para), &dma_rply);
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
        k = blockId_z * Para_s->cuBlock[2];

        d_d2uk = &Para_s->d_d2uk[(j-1) * imt];
        d_d2vk = &Para_s->d_d2vk[(j-1) * imt];
        d_duc  = &Para_s->d_duc[j * imt];
        d_dum  = &Para_s->d_dum[j * imt];
        d_dun  = &Para_s->d_dun[j * imt];
        d_dus  = &Para_s->d_dus[j * imt];
        d_due  = &Para_s->d_due[j * imt];
        d_duw  = &Para_s->d_duw[j * imt];
        d_dmc  = &Para_s->d_dmc[j * imt];
        d_dmn  = &Para_s->d_dmn[j * imt];
        d_dms  = &Para_s->d_dms[j * imt];
        d_dme  = &Para_s->d_dme[j * imt];
        d_dmw  = &Para_s->d_dmw[j * imt];
#ifdef TIDEFORCE
        d_h0_tide    = &Para_s->d_h0_tide[j*imt];
        d_h0_tidepj1 = &Para_s->d_h0_tide[(j+1)*imt];
#endif 
        d_h0   = &Para_s->d_h0[j * imt];
        d_h0pj1= &Para_s->d_h0[(j + 1) * imt];
        d_dxur = &Para_s->d_dxur[j * imt];
        d_dyur = &Para_s->d_dyur[j * imt];
        d_ebea = &Para_s->d_ebea[j * imt];
        d_ebeb = &Para_s->d_ebeb[j * imt];
        d_fcor = &Para_s->d_fcor[j * imt];

        d_dlub = &Para_s->d_dlub[j * imt];
        d_dlvb = &Para_s->d_dlvb[j * imt];
        d_vbp  = &Para_s->d_vbp[j * imt];
        d_ubp  = &Para_s->d_ubp[j * imt];
        d_pax  = &Para_s->d_pax[j * imt];
        d_pay  = &Para_s->d_pay[j * imt];
        d_pxb  = &Para_s->d_pxb[j * imt];
        d_pyb  = &Para_s->d_pyb[j * imt];
        d_viv  = &Para_s->d_viv[j * imt];
        d_whx  = &Para_s->d_whx[j * imt];
        d_why  = &Para_s->d_why[j * imt];
        d_wka0  = &Para_s->d_wka[j * imt];
        d_wka1  = &Para_s->d_wka[1 *jmt *imt+j *imt];
        d_wka2  = &Para_s->d_wka[2 *jmt *imt+j *imt];
        d_wka3  = &Para_s->d_wka[3 *jmt *imt+j *imt];
        d_wka4  = &Para_s->d_wka[4 *jmt *imt+j *imt];
        d_wka5  = &Para_s->d_wka[5 *jmt *imt+j *imt];
 

        for (int i = 0; i < imt; i++) {

            am_f=1.0e0;
            hduk = 0.0e0;
            hdvk = 0.0e0;

        	if(j>=Para_s->d_jb-2&&j<=Para_s->d_je&&i>=Para_s->d_ib-2&&i<=Para_s->d_ie)
        	{
        			temp = sqrt(Para_s->d_uarea_r[j *imt+i]);
        			u_area2 = temp * temp;
        			t_area = Para_s->d_am * temp / (pi);
        			scaling = t_area * t_area * t_area * t_area * t_area;
        			am_f = scaling*sqrt(Para_s->d_am_factor[j *imt+i]);
        			am_f = fmax(1.3e-5 * (u_area2 * u_area2)*Para_s->d_value2, am_f);
        			am_f = fmin(1.0/80.0 * (u_area2 * u_area2) / Para_s->d_dts, am_f);
        	}
        	am_f = fmax(Para_s->d_value, am_f);

            if(j>=Para_s->d_jb-1&&j<Para_s->d_je&&i>=Para_s->d_ib-1&&i<Para_s->d_ie)
            {
                if(Para_s->d_kmu[j *imt+i]-1>=0)
                {
                    hduk=Para_s->d_fac_vis_uv[j *imt+i]*((((d_duc[i]+
                            d_dum[i]) *d_d2uk[imt+i]+
                            d_dun[i] *d_d2uk[i]+
                            d_dus[i] *d_d2uk[2*imt+i]+
                            d_due[i] *d_d2uk[imt+i+1]+
                            d_duw[i] *d_d2uk[imt+i-1])+
                           (d_dmc[i] *d_d2vk[imt+i]+
                            d_dmn[i] *d_d2vk[i]+
                            d_dms[i] *d_d2vk[2*imt+i]+
                            d_dme[i] *d_d2vk[imt+i+1]+
                            d_dmw[i] *d_d2vk[imt+i-1])));

                    hdvk=Para_s->d_fac_vis_uv[j *imt+i]*((((d_duc[i]+
                              d_dum[i]) *d_d2vk[imt+i]+
                              d_dun[i] *d_d2vk[i]+
                              d_dus[i] *d_d2vk[2*imt+i]+
                              d_due[i] *d_d2vk[imt+i+1]+
                              d_duw[i] *d_d2vk[imt+i-1])-
                             (d_dmc[i] *d_d2uk[imt+i]+
                              d_dmn[i] *d_d2uk[i]+
                              d_dms[i] *d_d2uk[2*imt+i]+
                              d_dme[i] *d_d2uk[imt+i+1]+
                              d_dmw[i] *d_d2uk[imt+i-1])));
                }
            }
        	//HDUK = -1.0*am_factor(:,:,bid)*HDUK
			//HDVK = -1.0*am_factor(:,:,bid)*HDVK
			hduk = -1.0 * am_f * hduk;
			hdvk = -1.0 * am_f * hdvk;
			//Para_s->d_factor[j*imt+i] = am_f;

            double gradx=0.0e0,grady=0.0e0;
            double wka0,wka0_0,wka1,wka1_1;

            if(i> 1&&i<imt-2&&j> 1&&j<jmt-2)
            {

#ifdef TIDEFORCE
                if(Para_s->d_kmu[j *imt+i]-1>=0)
                {
// barotr.F90 <TIDEFORCE> "call grad" d_h0>h0(*,j)  d_h0pj1>h0(*,j+1)
                    gradx=d_dxur[i] *0.5 *(d_h0_tidepj1[i]-d_h0_tide[i-1]-d_h0_tidepj1[i-1]+d_h0_tide[i]);

                    grady=d_dyur[i] *0.5 *(d_h0_tidepj1[i]-d_h0_tide[i-1]+d_h0_tidepj1[i-1]-d_h0_tide[i]);
                }
#else 
                if(Para_s->d_kmu[j *imt+i]-1>=0)
                {
// barotr.F90 "call grad" d_h0>h0(*,j)  d_h0pj1>h0(*,j+1)
                    gradx=d_dxur[i] *0.5 *(d_h0pj1[i]-d_h0[i-1]-d_h0pj1[i-1]+d_h0[i]);

                    grady=d_dyur[i] *0.5 *(d_h0pj1[i]-d_h0[i-1]+d_h0pj1[i-1]-d_h0[i]);
                }
#endif 
// barote.F90 Line210
                double gstar=(Para_s->d_wgp[j *imt+i]-1.0) *g;
                wka0 = hduk+gstar *gradx;
                d_wka0[i]=hduk+gstar *gradx;
                wka1 = hdvk+gstar *grady;
                d_wka1[i]=hdvk+gstar *grady;
                d_wka4[i]=hduk;
                d_wka5[i]=hdvk;
            }
            if(j==jmt-1||i==0)
            {
                d_work[i]=0.0e0;
            }
            else{
// barotr.F90  call tgrid_to_ugrid
                d_work[i]=p25 *d_h0[i]+ p25 *d_h0pj1[i]+
                        p25 *d_h0[i-1]+ p25 *d_h0pj1[i-1];
            }
            if(i> 1&&i<imt-2&&j> 1&&j<jmt-2)
            {

// barotr.F90 "COMPUTING DU & DV"

                //d_wka_s[j *imt+i]=Para_s->d_viv[j *imt+i] *(d_wka_s[j *imt+i]+
                wka0_0=d_viv[i] *(wka0+d_dlub[i]-d_fcor[i] *d_vbp[i]+
                        d_pax[i]+d_pxb[i]-d_work[i] *d_whx[i]);

                d_wka0[i] = wka0_0;

                //d_wka_s[1 *jmt *imt+j *imt+i]=Para_s->d_viv[j *imt+i] *
                wka1_1=d_viv[i] *
                       //(d_wka_s[1 *jmt *imt+j *imt+i]+Para_s->d_dlvb[j *imt+i]+
                       (wka1+d_dlvb[i]+d_fcor[i] *d_ubp[i]+
                        d_pay[i]+d_pyb[i]-d_work[i] *d_why[i]);

                d_wka1[i] = wka1_1;

// barotr.F90 "CORIOLIS ADJUSTMENT"

                d_wka2[i]=d_ebea[i] *
                                              //d_wka_s[j *imt+i]-Para_s->d_ebeb[j *imt+i] *
                                              wka0_0-d_ebeb[i] *
                                                     //d_wka_s[1 *jmt *imt+j *imt+i];
                                                     wka1_1;

                d_wka3[i]=d_ebea[i] *
                                              //d_wka_s[1 *jmt *imt+j *imt+i]+
                                              wka1_1+
                                              //Para_s->d_ebeb[j *imt+i] *d_wka_s[j *imt+i];
                                              d_ebeb[i] *wka0_0;
            }

        }//for i

        CRTS_dma_put(&Para_s->d_work[j * imt], d_work, sizeof(double) * imt);
//        CRTS_dma_put(&Para_s->d_wka[j * imt], d_wka0, sizeof(double) * imt);
//        CRTS_dma_put(&Para_s->d_wka[1 *jmt *imt+j *imt], d_wka1, sizeof(double) * imt);
//        CRTS_dma_put(&Para_s->d_wka[2 *jmt *imt+j *imt], d_wka2, sizeof(double) * imt);
//        CRTS_dma_put(&Para_s->d_wka[3 *jmt *imt+j *imt], d_wka3, sizeof(double) * imt);
//        CRTS_dma_put(&Para_s->d_wka[4 *jmt *imt+j *imt], d_wka4, sizeof(double) * imt);
//        CRTS_dma_put(&Para_s->d_wka[5 *jmt *imt+j *imt], d_wka5, sizeof(double) * imt);

    }//for block

    ldm_free(Para_s, sizeof(barotr2_udel4_1_Para));
}
