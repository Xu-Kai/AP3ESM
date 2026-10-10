
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

//__userDef__ paraVarAttr(double,d_ricdttms,sizeof(double)*kmm1*jmt*imt,outonly,yes,no,no);	
//__userDef__ paraVarAttr(double,d_ricdt,sizeof(double)*kmm1*jmt*imt,outonly,yes,no,no);	
//__userDef__ paraVarAttr(double,d_gg,sizeof(double)*km*jmt*imt,outonly,yes,no,no);	
//__userDef__ paraVarAttr(double,d_pdensity,sizeof(double)*(km+1)*jmt*imt,readonly,no,no,no);	
//__userDef__ paraVarAttr(double,d_vit,sizeof(double)*jmt*imt*km,readonly,no,no,no);	
//__userDef__ paraVarAttr(double,d_pp,sizeof(double)*km*jmt*imt,readonly,no,no,no);	
//__userDef__ paraVarAttr(double,d_at,sizeof(double)*ntra*km*jmt*imt,readonly,no,no,no);	
//__userDef__ paraVarAttr(double,d_po,sizeof(double)*km,readonly,no,yes,no);	
//__userDef__ paraVarAttr(double,d_dlu,sizeof(double)*jmt*imt*km,outonly,yes,no,no);	
//__userDef__ paraVarAttr(double,d_dlv,sizeof(double)*jmt*imt*km,outonly,yes,no,no);	
//__userDef__ paraVarAttr(double,d_dxur,sizeof(double)*ny_block*nx_block,readonly,no,no,no);	
//__userDef__ paraVarAttr(double,d_dyur,sizeof(double)*ny_block*nx_block,readonly,no,no,no);	
//__userDef__ paraVarAttr(int,d_kmu,sizeof(int)*ny_block*nx_block,readonly,no,no,no);	
//__userDef__ paraVarAttr(double,d_buoytur,sizeof(double)*jmt*imt,outonly,yes,no,no);	
//__userDef__ paraVarAttr(double,d_buoysol,sizeof(double)*jmt*imt,outonly,yes,no,no);	
//__userDef__ paraVarAttr(double,d_nswv,sizeof(double)*jmt*imt,readonly,no,no,no);	
//__userDef__ paraVarAttr(double,d_swv,sizeof(double)*jmt*imt,readonly,no,no,no);	
//__userDef__ paraVarAttr(double,d_alpha,sizeof(double)*km*jmt*imt,readonly,no,no,no);	
//__userDef__ paraVarAttr(double,d_beta,sizeof(double)*km*jmt*imt,readonly,no,no,no);	
//__userDef__ paraVarAttr(double,d_odzt,sizeof(double)*km,readonly,no,yes,no);	
//__userDef__ paraVarAttr(double,d_od0cp,sizeof(double),readonly,no,no,no);	
//__userDef__ paraVarAttr(double,d_od0,sizeof(double),readonly,no,no,no);	
//__userDef__ paraVarAttr(double,d_g,sizeof(double),readonly,no,no,no);	
//__userDef__ paraVarAttr(int,d_kmm1,sizeof(int),readonly,no,no,no);	
//__userDef__ paraVarAttr(int,d_km,sizeof(int),readonly,no,no,no);	
//__userDef__ paraVarAttr(int,d_jmt,sizeof(int),readonly,no,no,no);	
//__userDef__ paraVarAttr(int,d_imt,sizeof(int),readonly,no,no,no);	
typedef struct readyt_cu_3_para{
	int cuBlock[3];
	int cuGrid[3];
	int d_kmm1;
	int d_km;
	int d_jmt;
	int d_imt;
	int rank;
	int *d_kmu;
        double d_od0cp;
	double d_od0;
  	double d_g;	
	double *d_ricdttms;
	double *d_ricdt;
	double *d_gg;
	double *d_pdensity;
	double *d_vit;
	double *d_pp;
	double *d_at;
	double *d_po;
	double *d_dlu;
	double *d_dlv;
	double *d_dxur;
	double *d_dyur;
	double *d_buoytur;
	double *d_buoysol;
	double *d_nswv;
	double *d_swv;
	double *d_alpha;
	double *d_beta;
	double *d_odzt;
  double *d_akt;
#if (defined CANUTOMIXOUT)	
  	double *d_alpha_canuto;
  	double *d_beta_canuto;
#endif
}readyt_3_para;


//dynamic atomic definition
#define CP_CORE      64
#define RCV_CORE     0
#define DATAPAKAGE   1

//ATHREAD global definition 
__thread_local int warp_size = 32; 
__thread_local int blockIndex;
//__thread_local Para* Paras;

 void  readyt_cu_3_(readyt_3_para *Paras){

  //ATHREAD local variable definition
  int swGrid[3];
  int swBlock[3];
  int stride;
  int blockId_x,blockId_y,blockId_z;
  int threadId_x,threadId_y,threadId_z;


  //ATHREAD basic operation
  int i,j,k;
  blockIndex = CRTS_tid;
//  readyt_3_para* Paras = (readyt_3_para*)ldm_malloc(sizeof(readyt_3_para));
//  CRTS_dma_iget(Paras, Paras, sizeof(readyt_3_para), &dma_rply);
//  D_COUNT++;
//  CRTS_dma_wait_value(&dma_rply, D_COUNT);

  double* d_po_s = (double*)ldm_malloc(sizeof(double) *km);
  CRTS_dma_iget(d_po_s, Paras->d_po, sizeof(double) *km, &dma_rply);
  D_COUNT++;
  CRTS_dma_wait_value(&dma_rply, D_COUNT);

  double* d_odzt_s = (double*)ldm_malloc(sizeof(double) *km);
  CRTS_dma_iget(d_odzt_s, Paras->d_odzt, sizeof(double) *km, &dma_rply);
  D_COUNT++;
  CRTS_dma_wait_value(&dma_rply, D_COUNT);

//  if(Paras->rank == 0)
//      printf("[slave %d]: cuGrid:%d,%d,%d,cuBlock:%d,%d,%d,kmm1:%d,km:%d,imt:%d,jmt:%d,rank:%d\n",
//             blockIndex,Paras->cuGrid[0],Paras->cuGrid[1],Paras->cuGrid[2],
//             Paras->cuBlock[0],Paras->cuBlock[1],Paras->cuBlock[2],Paras->d_kmm1,Paras->d_km,
//             imt,jmt,Paras->rank);

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

  double *d_ricdttms_s = Paras->d_ricdttms;
  double *d_ricdt_s = Paras->d_ricdt;
  double *d_gg_s = Paras->d_gg;
  double *d_dlu_s = Paras->d_dlu;
  double *d_dlv_s = Paras->d_dlv;
  double *d_buoytur_s = Paras->d_buoytur;
  double *d_buoysol_s = Paras->d_buoysol;

  for(int block = blockStart; block < blockEnd; block++){
    blockId_z = block / (Paras->cuGrid[0] * Paras->cuGrid[1]);
    blockId_y = (block - blockId_z * (Paras->cuGrid[0] * Paras->cuGrid[1])) / Paras->cuGrid[0];
    blockId_x = block % Paras->cuGrid[0];
      j = blockId_y * Paras->cuBlock[1];
      k = blockId_z * Paras->cuBlock[2];

    //ATHREAD serial operation replace GPU parallel mutli-thread in single thread-block
    for(i= 0; i < threadSum; i++) {

      
      if(k<Paras->d_km&&j<jmt&&i<imt)	
      {	
        double epsln=1.0e-25;	
	
        if(k==0)	
        {	
      	  d_buoytur_s[j *imt+i]=Paras->d_vit[j*imt+i] *Paras->d_nswv[j*imt+i] *
      				Paras->d_g * Paras->d_alpha[j*imt+i] *Paras->d_od0cp;	
      	  d_buoysol_s[j *imt+i]=Paras->d_vit[j*imt+i] *Paras->d_swv[j *imt+i] *
      				Paras->d_g * Paras->d_alpha[j*imt+i] *Paras->d_od0cp;	
        }	
        if(k<Paras->d_kmm1&&1<=i&&i<imt-1)	
        {	
      	  d_ricdttms_s[k *jmt *imt+j *imt+i]=
      		Paras->d_vit[(k+1) *jmt *imt+j *imt+i] *Paras->d_g *
      		((Paras->d_at[k *jmt *imt+j *imt+i]-
      		Paras->d_at[(k+1) *jmt *imt+j *imt+i]) *
      		Paras->d_alpha[(k+1) *jmt *imt+j *imt+i]+
      		1000.0e0 *(Paras->d_at[1 *Paras->d_km *jmt *imt+k *jmt *imt+j *imt+i]-
      		Paras->d_at[1 *Paras->d_km *jmt *imt+(k+1) *jmt *imt+j *imt+i]) *
      		Paras->d_beta[(k+1) *jmt *imt+j *imt+i]) *d_odzt_s[k+1];

      	  d_ricdt_s[k *jmt *imt+j *imt+i]=
      		Paras->d_vit[(k+1)*jmt*imt+j*imt+i] /
      		((Paras->d_at[k *jmt *imt+j *imt+i]- Paras->d_at[(k+1) *jmt *imt+j *imt+i]+epsln) *
      		Paras->d_alpha[(k+1) *jmt *imt+j *imt+i]) *1000.0e0 *
          ((Paras->d_at[1 *Paras->d_km *jmt *imt+k *jmt *imt+j *imt+i]-
      		Paras->d_at[1 *Paras->d_km *jmt *imt+(k+1) *jmt *imt+j *imt+i]) *
      		Paras->d_beta[(k+1) *jmt *imt+j *imt+i]);
        }	
	#if ( defined CANUTOMIXOUT )	
      	//pconst_mod_mp_alpha_canuto_[iblock][k][j][i]=alpha[iblock][k][j][i];	
      	//pconst_mod_mp_beta_canuto_[iblock][k][j][i]=beta[iblock][k][j][i];	
      	Paras->d_alpha_canuto[iblock][k][j][i] = Paras->d_alpha[iblock][k][j][i];	
      	Paras->d_beta_canuto[iblock][k][j][i]   = Paras->d_beta[iblock][k][j][i];
      		
	#endif	
      	d_gg_s[k *jmt *imt+j *imt+i]=-Paras->d_od0 *Paras->d_g *
                  (Paras->d_pdensity[k *jmt *imt+j *imt+i]-d_po_s[k]-1000.0e0) *
                  Paras->d_vit[k *jmt *imt+j *imt+i];
      	d_dlu_s[k *jmt *imt+j *imt+i]=0.0;	
      	d_dlv_s[k *jmt *imt+j *imt+i]=0.0;
//        Paras->d_akt[k *jmt *imt+j *imt+i] = 0.0;
//        Paras->d_akt[km *jmt *imt+k *jmt *imt+j *imt+i] = 0.0;

        if((1<=i)&&(j<(jmt-1)))	
        {	
          if(k<=(Paras->d_kmu[j *imt+i]-1))	
          {	
      	    d_dlu_s[k *jmt *imt+j *imt+i]=Paras->d_dxur[j *imt+i] *p5 *
                  (Paras->d_pp[k *jmt *imt+(j+1) *imt+i]-
                  Paras->d_pp[k *jmt *imt+j *imt+i-1]-
                  Paras->d_pp[k *jmt *imt+(j+1) *imt+i-1]+
                  Paras->d_pp[k *jmt *imt+j *imt+i]);

      	    d_dlv_s[k *jmt *imt+j *imt+i]=Paras->d_dyur[j *imt+i] *p5 *
                  (Paras->d_pp[k *jmt *imt+(j+1) *imt+i]-
                  Paras->d_pp[k *jmt *imt+j *imt+i-1]+
                  Paras->d_pp[k *jmt *imt+(j+1) *imt+i-1]-
                  Paras->d_pp[k *jmt *imt+j *imt+i]);
          }	
        }	
      }	

    }

  }

    ldm_free(d_po_s, sizeof(double) *km);

    ldm_free(d_odzt_s, sizeof(double) *km);

 //   ldm_free(Paras, sizeof(readyt_3_para));
}
