//ATHREAD include file
#include <crts.h>
#include <slave.h>
//#include <simd.h>
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


//__userDef__ paraVarAttr(double,d_h0l,sizeof(double)*jmt*imt,inout,yes,no,no);	
//__userDef__ paraVarAttr(double,d_h0f,sizeof(double)*jmt*imt,inout,yes,no,no);	
//__userDef__ paraVarAttr(double,d_h0,sizeof(double)*jmt*imt,readonly,no,no,no);	
//__userDef__ paraVarAttr(double,d_utl,sizeof(double)*jmt*imt*km,inout,yes,no,no);	
//__userDef__ paraVarAttr(double,d_vtl,sizeof(double)*jmt*imt*km,inout,yes,no,no);	
//__userDef__ paraVarAttr(double,d_utf,sizeof(double)*jmt*imt*km,inout,yes,no,no);	
//__userDef__ paraVarAttr(double,d_vtf,sizeof(double)*jmt*imt*km,inout,yes,no,no);	
//__userDef__ paraVarAttr(double,d_u,sizeof(double)*jmt*imt*km,readonly,no,no,no);	
//__userDef__ paraVarAttr(double,d_v,sizeof(double)*jmt*imt*km,readonly,no,no,no);	
//__userDef__ paraVarAttr(double,d_dlu,sizeof(double)*jmt*imt*km,inout,yes,no,no);	
//__userDef__ paraVarAttr(double,d_dlv,sizeof(double)*jmt*imt*km,inout,yes,no,no);	
//__userDef__ paraVarAttr(double,d_atb,sizeof(double)*ntra*(km+1)*jmt*imt,readonly,no,no,no);	
//__userDef__ paraVarAttr(double,d_au0,sizeof(double)*ny_block*nx_block,readonly,no,no,no);	
//__userDef__ paraVarAttr(double,d_aus,sizeof(double)*ny_block*nx_block,readonly,no,no,no);	
//__userDef__ paraVarAttr(double,d_auw,sizeof(double)*ny_block*nx_block,readonly,no,no,no);	
//__userDef__ paraVarAttr(double,d_ausw,sizeof(double)*ny_block*nx_block,readonly,no,no,no);	
typedef struct readyt_1_kernelPara{
    int    cuBlock[3];
    int    cuGrid[3];
    int    rank;
    int    pad;
	double *d_h0l;
	double *d_h0f;
	double *d_h0;
	double *d_utl;
	double *d_vtl;
	double *d_utf;
	double *d_vtf;
	double *d_u;
	double *d_v;
	double *d_dlu;
	double *d_dlv;
	double *d_atb;
	double *d_au0;
	double *d_aus;
	double *d_auw;
	double *d_ausw;
}readyt_1_para;


//dynamic atomic definition
#define CP_CORE      64
#define RCV_CORE     0
#define DATAPAKAGE   1

//ATHREAD global definition
__thread_local int blockIndex;

void  readyt_cu_1_(readyt_1_para *Paras){

  //ATHREAD local variable definition
  int swGrid[3];
  int swBlock[3];
  int stride;
  int blockId_x,blockId_y,blockId_z;
  int threadId_x,threadId_y,threadId_z;

  //ATHREAD basic operation
  int i,j,k;
  blockIndex = CRTS_tid;
//  readyt_1_para* Paras = (readyt_1_para*)ldm_malloc(sizeof(readyt_1_para));
//  CRTS_dma_iget(Paras, Paras, sizeof(readyt_1_para), &dma_rply);
//  D_COUNT++;
//  CRTS_dma_wait_value(&dma_rply, D_COUNT);

  //cuda grid and block affine translation
  int blockSum   = Paras->cuGrid[0] * Paras->cuGrid[1] * Paras->cuGrid[2];
  int threadSum  = Paras->cuBlock[0] * Paras->cuBlock[1] * Paras->cuBlock[2];
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
      k = blockId_z * Paras->cuBlock[2];
             
      for(i = 0; i < threadSum; i++) {

			int n;


            //if(0&&jst-1<=j&&j<jet&&0<=i&&i<imt)
            if(k==0 && jst-1<=j && j<jet && 0<= i && i<imt )
            {
                //d_h0l_s[j *imt+i]=d_h0f_s[j *imt+i];
                //d_h0l_s[threadId_x]=d_h0f_s[threadId_x];
                Paras->d_h0l[j *imt+i]=Paras->d_h0f[j *imt+i];
                //d_h0f_s[j *imt+i]=Paras->d_h0[j *imt+i];
                //d_h0f_s[threadId_x]=Paras->d_h0[j *imt+i];
                Paras->d_h0f[j *imt+i]=Paras->d_h0[j *imt+i];
            }
            if(0<=k&&k<km&&jst-1<=j&&j<jet&&0<=i&&i<imt)
            {
                //d_utl_s[k *jmt *imt+j *imt+i]=d_utf_s[k *jmt *imt+j *imt+i];
                //d_vtl_s[k *jmt *imt+j *imt+i]=d_vtf_s[k *jmt *imt+j *imt+i];
                //d_utf_s[k *jmt *imt+j *imt+i]=Paras->d_u[k *jmt *imt+j *imt+i];
                //d_vtf_s[k *jmt *imt+j *imt+i]=Paras->d_v[k *jmt *imt+j *imt+i];
							
//                d_utl_s[threadId_x]=d_utf_s[threadId_x];
//                d_vtl_s[threadId_x]=d_vtf_s[threadId_x];
//                d_utf_s[threadId_x]=Paras->d_u[k *jmt *imt+j *imt+i];
//                d_vtf_s[threadId_x]=Paras->d_v[k *jmt *imt+j *imt+i];

                Paras->d_utl[k *jmt *imt+j *imt+i]=Paras->d_utf[k *jmt *imt+j *imt+i];
                Paras->d_vtl[k *jmt *imt+j *imt+i]=Paras->d_vtf[k *jmt *imt+j *imt+i];
                Paras->d_utf[k *jmt *imt+j *imt+i]=Paras->d_u[k *jmt *imt+j *imt+i];
                Paras->d_vtf[k *jmt *imt+j *imt+i]=Paras->d_v[k *jmt *imt+j *imt+i];
            }
            if(k<km&&j<jmt&&i<imt)
            {
                if(j<(jmt-1)&&i>=1)
                {
                    //d_dlu_s[k *jmt *imt+j *imt+i]=Paras->d_au0[j *imt+i] *Paras->d_atb[(k+1) *jmt *imt+j *imt+i]+Paras->d_aus[j *imt+i] *Paras->d_atb[(k+1) *jmt *imt+(j+1) *imt+i]+Paras->d_auw[j *imt+i] *Paras->d_atb[(k+1) *jmt *imt+j *imt+i-1]+Paras->d_ausw[j *imt+i] *Paras->d_atb[(k+1) *jmt *imt+(j+1) *imt+i-1];
                    //d_dlu_s[threadId_x]=Paras->d_au0[j *imt+i] *Paras->d_atb[(k+1) *jmt *imt+j *imt+i]+Paras->d_aus[j *imt+i] *Paras->d_atb[(k+1) *jmt *imt+(j+1) *imt+i]+Paras->d_auw[j *imt+i] *Paras->d_atb[(k+1) *jmt *imt+j *imt+i-1]+Paras->d_ausw[j *imt+i] *Paras->d_atb[(k+1) *jmt *imt+(j+1) *imt+i-1];

                    Paras->d_dlu[k *jmt *imt+j *imt+i]=p25 *Paras->d_atb[(k+1) *jmt *imt+j *imt+i]+
                            p25 *Paras->d_atb[(k+1) *jmt *imt+(j+1) *imt+i]+
                            p25 *Paras->d_atb[(k+1) *jmt *imt+j *imt+i-1]+
                            p25 *Paras->d_atb[(k+1) *jmt *imt+(j+1) *imt+i-1];
								
                    //d_dlv_s[k *jmt *imt+j *imt+i]=Paras->d_au0[j *imt+i] *Paras->d_atb[1 *(km+1) *jmt *imt+(k+1) *jmt *imt+j *imt+i]+Paras->d_aus[j *imt+i] *Paras->d_atb[1 *(km+1) *jmt *imt+(k+1) *jmt *imt+(j+1) *imt+i]+Paras->d_auw[j *imt+i] *Paras->d_atb[1 *(km+1) *jmt *imt+(k+1) *jmt *imt+j *imt+i-1]+Paras->d_ausw[j *imt+i] *Paras->d_atb[1 *(km+1) *jmt *imt+(k+1) *jmt *imt+(j+1) *imt+i-1];
                    //d_dlv_s[threadId_x]=Paras->d_au0[j *imt+i] *Paras->d_atb[1 *(km+1) *jmt *imt+(k+1) *jmt *imt+j *imt+i]+Paras->d_aus[j *imt+i] *Paras->d_atb[1 *(km+1) *jmt *imt+(k+1) *jmt *imt+(j+1) *imt+i]+Paras->d_auw[j *imt+i] *Paras->d_atb[1 *(km+1) *jmt *imt+(k+1) *jmt *imt+j *imt+i-1]+Paras->d_ausw[j *imt+i] *Paras->d_atb[1 *(km+1) *jmt *imt+(k+1) *jmt *imt+(j+1) *imt+i-1];

                    Paras->d_dlv[k *jmt *imt+j *imt+i]=p25 *Paras->d_atb[1 *(km+1) *jmt *imt+(k+1) *jmt *imt+j *imt+i]+
                            p25 *Paras->d_atb[1 *(km+1) *jmt *imt+(k+1) *jmt *imt+(j+1) *imt+i]+
                            p25 *Paras->d_atb[1 *(km+1) *jmt *imt+(k+1) *jmt *imt+j *imt+i-1]+
                            p25 *Paras->d_atb[1 *(km+1) *jmt *imt+(k+1) *jmt *imt+(j+1) *imt+i-1];
                }
                else{
                    //d_dlu_s[k *jmt *imt+j *imt+i]=(double)0.0e0;
                    //d_dlv_s[k *jmt *imt+j *imt+i]=(double)0.0e0;
								
//                    d_dlu_s[threadId_x]=(double)0.0e0;
//                    d_dlv_s[threadId_x]=(double)0.0e0;
                    Paras->d_dlu[k *jmt *imt+j *imt+i]=(double)0.0e0;
                    Paras->d_dlv[k *jmt *imt+j *imt+i]=(double)0.0e0;
                }
            }

        }


    }

//    ldm_free(Paras, sizeof(readyt_1_para));
}
