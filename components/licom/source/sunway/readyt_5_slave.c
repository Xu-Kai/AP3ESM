
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
__thread_local crts_rply_t l_rply = 0;
__thread_local crts_rply_t r_rply = 0;

//__userDef__ paraVarAttr(double,d_dlv,sizeof(double)*jmt*imt*km,readonly,no,no,no);	
//__userDef__ paraVarAttr(double,d_au0,sizeof(double)*ny_block*nx_block,readonly,no,no,no);	
//__userDef__ paraVarAttr(double,d_aus,sizeof(double)*ny_block*nx_block,readonly,no,no,no);	
//__userDef__ paraVarAttr(double,d_auw,sizeof(double)*ny_block*nx_block,readonly,no,no,no);	
//__userDef__ paraVarAttr(double,d_ausw,sizeof(double)*ny_block*nx_block,readonly,no,no,no);	
//__userDef__ paraVarAttr(double,d_wgp,sizeof(double)*jmt*imt,outonly,yes,no,no);	
//__userDef__ paraVarAttr(double,d_work,sizeof(double)*jmt*imt,outonly,yes,no,no);	
//__userDef__ paraVarAttr(double,d_hbx,sizeof(double)*jmt*imt,readonly,no,no,no);	
//__userDef__ paraVarAttr(double,d_hby,sizeof(double)*jmt*imt,readonly,no,no,no);	
//__userDef__ paraVarAttr(double,d_viv,sizeof(double)*jmt*imt*km,readonly,no,no,no);	
//__userDef__ paraVarAttr(double,d_whx,sizeof(double)*jmt*imt,outonly,yes,no,no);	
//__userDef__ paraVarAttr(double,d_why,sizeof(double)*jmt*imt,outonly,yes,no,no);	
//__userDef__ paraVarAttr(double,d_work1,sizeof(double)*jmt*imt,outonly,yes,no,no);	
//__userDef__ paraVarAttr(double,d_work2,sizeof(double)*jmt*imt,outonly,yes,no,no);	
//__userDef__ paraVarAttr(double,d_psa,sizeof(double)*jmt*imt,readonly,no,no,no);	
//__userDef__ paraVarAttr(double,d_dxur,sizeof(double)*ny_block*nx_block,readonly,no,no,no);	
//__userDef__ paraVarAttr(double,d_dyur,sizeof(double)*ny_block*nx_block,readonly,no,no,no);	
//__userDef__ paraVarAttr(int,d_kmu,sizeof(int)*ny_block*nx_block,readonly,no,no,no);	
//__userDef__ paraVarAttr(double,d_pax,sizeof(double)*jmt*imt,outonly,yes,no,no);	
//__userDef__ paraVarAttr(double,d_pay,sizeof(double)*jmt*imt,outonly,yes,no,no);	
//__userDef__ paraVarAttr(double,d_od0,sizeof(double),readonly,no,no,no);	
//__userDef__ paraVarAttr(int,d_jmt,sizeof(int),readonly,no,no,no);	
//__userDef__ paraVarAttr(int,d_imt,sizeof(int),readonly,no,no,no);	
typedef struct readyt_5_kernelPara{
	int cuBlock[3];
	int cuGrid[3];
	int d_jmt;
	int d_imt;
	int rank;
	int *d_kmu;
	double d_od0;
	double *d_dlv;
	double *d_au0;
	double *d_aus;
	double *d_auw;
	double *d_ausw;
	double *d_wgp;
	double *d_work;
	double *d_hbx;
	double *d_hby;
	double *d_viv;
	double *d_whx;
	double *d_why;
	double *d_work1;
	double *d_work2;
	double *d_psa;
	double *d_dxur;
	double *d_dyur;
	double *d_pax;
	double *d_pay;
}readyt_5_para;


//dynamic atomic definition
#define CP_CORE      64
#define RCV_CORE     0
#define DATAPAKAGE   1

//ATHREAD global definition
__thread_local int blockIndex;
//__thread_local Para* Para_s;

void  readyt_cu_5_(readyt_5_para *Paras){

  //ATHREAD local variable definition
  int swGrid[3];
  int swBlock[3];
  int stride;
  int blockId_x,blockId_y,blockId_z;
  int threadId_x,threadId_y,threadId_z;


  //ATHREAD basic operation
  int i,j;
  blockIndex = CRTS_tid;
  readyt_5_para* Para_s = (readyt_5_para*)ldm_malloc(sizeof(readyt_5_para));
  CRTS_dma_iget(Para_s, Paras, sizeof(readyt_5_para), &dma_rply);
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

/*
  double *d_work_s = Para_s->d_work;

  double *d_work1_s = Para_s->d_work1;
  double *d_work2_s = Para_s->d_work2;
*/
  double *d_wgp_s = Para_s->d_wgp;
  double *d_work_s = Para_s->d_work;
  //double d_work_s[jmt*imt];
  double *d_whx_s = Para_s->d_whx;
  double *d_why_s = Para_s->d_why;
  double *d_work1_s = Para_s->d_work1;
  double *d_work2_s = Para_s->d_work2;
  //double d_work1_s[jmt*imt];
  //double d_work2_s[jmt*imt];
  double *d_pax_s = Para_s->d_pax;
  double *d_pay_s = Para_s->d_pay;

  for(int block = blockStart; block < blockEnd; block++){
    blockId_z = block / (Para_s->cuGrid[0] * Para_s->cuGrid[1]);
    blockId_y = (block - blockId_z * (Para_s->cuGrid[0] * Para_s->cuGrid[1])) / Para_s->cuGrid[0];
    blockId_x = block % Para_s->cuGrid[0];
      j = blockId_y * Para_s->cuBlock[1];

    //ATHREAD serial operation replace GPU parallel mutli-thread in single thread-block
    for(i = 0; i < threadSum; i++) {

      if(j<Para_s->d_jmt&&i<Para_s->d_imt)	
      {	
        if(j<(Para_s->d_jmt-1)&&1<=i)
        {
      	  d_wgp_s[j *Para_s->d_imt+i]=Para_s->d_au0[j *Para_s->d_imt+i] *Para_s->d_dlv[j *Para_s->d_imt+i]+
                    Para_s->d_aus[j *Para_s->d_imt+i] *Para_s->d_dlv[(j+1) *Para_s->d_imt+i]+
                    Para_s->d_auw[j *Para_s->d_imt+i] *Para_s->d_dlv[j *Para_s->d_imt+i-1]+
                    Para_s->d_ausw[j *Para_s->d_imt+i] *Para_s->d_dlv[(j+1) *Para_s->d_imt+i-1];


            d_work_s[j *Para_s->d_imt+i] =
                    Para_s->d_au0[j *Para_s->d_imt+i] *Para_s->d_dlv[1 *Para_s->d_jmt *Para_s->d_imt+j *Para_s->d_imt+i]+
                    Para_s->d_aus[j *Para_s->d_imt+i] *Para_s->d_dlv[1 *Para_s->d_jmt *Para_s->d_imt+(j+1) *Para_s->d_imt+i]+
                    Para_s->d_auw[j *Para_s->d_imt+i] *Para_s->d_dlv[1 *Para_s->d_jmt *Para_s->d_imt+j *Para_s->d_imt+i-1]+
                    Para_s->d_ausw[j *Para_s->d_imt+i] *Para_s->d_dlv[1 *Para_s->d_jmt *Para_s->d_imt+(j+1) *Para_s->d_imt+i-1];

       }
       else	{
           d_wgp_s[j *Para_s->d_imt+i]=(double)0.0e0;
           d_work_s[j *Para_s->d_imt+i] = (double)0.0e0;;
       }
      	d_whx_s[j *Para_s->d_imt+i]=Para_s->d_hbx[j *Para_s->d_imt+i] *d_work_s[j *Para_s->d_imt+i] *
                  Para_s->d_viv[j *Para_s->d_imt+i];

      	d_why_s[j *Para_s->d_imt+i]=Para_s->d_hby[j *Para_s->d_imt+i] *d_work_s[j *Para_s->d_imt+i] *
                  Para_s->d_viv[j *Para_s->d_imt+i];

        d_work1_s[j *Para_s->d_imt+i]=0.0;
      	d_work2_s[j *Para_s->d_imt+i]=0.0;

      if((1<=i)&&(j<(Para_s->d_jmt-1)))	
      {	
        if(0<=Para_s->d_kmu[j *Para_s->d_imt+i]-1)
        {
      	  d_work1_s[j *Para_s->d_imt+i]=Para_s->d_dxur[j *Para_s->d_imt+i] *p5 *(Para_s->d_psa[(j+1) *Para_s->d_imt+i]-
                    Para_s->d_psa[j *Para_s->d_imt+i-1]-Para_s->d_psa[(j+1) *Para_s->d_imt+i-1]+
                    Para_s->d_psa[j *Para_s->d_imt+i]);


      	  d_work2_s[j *Para_s->d_imt+i]=Para_s->d_dyur[j *Para_s->d_imt+i] *p5 *(Para_s->d_psa[(j+1) *Para_s->d_imt+i]-
                    Para_s->d_psa[j *Para_s->d_imt+i-1]+Para_s->d_psa[(j+1) *Para_s->d_imt+i-1]-
                    Para_s->d_psa[j *Para_s->d_imt+i]);

        }
      }

      	d_pay_s[j *Para_s->d_imt+i]=-Para_s->d_od0 *d_work2_s[j *Para_s->d_imt+i];

        d_pax_s[j *Para_s->d_imt+i]=-Para_s->d_od0 *d_work1_s[j *Para_s->d_imt+i];

      }	

    }

  }

/*
    CRTS_dma_iput(Para_s->d_work, d_work_s, sizeof(double)*jmt*imt, &dma_rply);
    D_COUNT++;
    CRTS_dma_wait_value(&dma_rply, D_COUNT);

    CRTS_dma_iput(Para_s->d_work1, d_work1_s, sizeof(double)*jmt*imt, &dma_rply);
    D_COUNT++;
    CRTS_dma_wait_value(&dma_rply, D_COUNT);

    CRTS_dma_iput(Para_s->d_work2, d_work2_s, sizeof(double)*jmt*imt, &dma_rply);
    D_COUNT++;
    CRTS_dma_wait_value(&dma_rply, D_COUNT);
*/
    ldm_free(Para_s, sizeof(readyt_5_para));
}
