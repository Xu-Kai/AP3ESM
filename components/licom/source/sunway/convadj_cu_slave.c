
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

//__userDef__ paraVarAttr(double,d_to,sizeof(double)*km,readonly,no,no,no);	
//__userDef__ paraVarAttr(double,d_so,sizeof(double)*km,readonly,no,no,no);	
//__userDef__ paraVarAttr(double,d_pconst_dts,sizeof(double),readonly,no,yes,no);	
//__userDef__ paraVarAttr(double,d_kmt,sizeof(int)*ny_block*nx_block,readonly,no,no,no);	
//__userDef__ paraVarAttr(double,d_vit,sizeof(double)*jmt*imt*km,readonly,no,no,no);	
//__userDef__ paraVarAttr(double,d_dzp,sizeof(double)*km,readonly,no,no,no);	
//__userDef__ paraVarAttr(double,d_at,sizeof(double)*ntra*km*jmt*imt,inout,yes,no,no);	
//__userDef__ paraVarAttr(double,d_tend,sizeof(double)*ntra*km*jmt*imt,outonly,yes,no,no);	
//__userDef__ paraVarAttr(double,d_dt_conv,sizeof(double)*ntra*km*jmt*imt,inout,yes,no,no);	
//__userDef__ paraVarAttr(double,d_atb,sizeof(double)*ntra*(km+1)*jmt*imt,inout,yes,no,no);	
//__userDef__ paraVarAttr(double,d_icmon,sizeof(double)*jmt*imt*2,outonly,yes,no,no);	
//__userDef__ paraVarAttr(double,d_c,sizeof(double)*km*9,readonly,no,no,no);	
typedef struct kernelPara{
	int cuBlock[3];
	int cuGrid[3];
	int rank;
	int pad;
	double *d_to;
	double *d_so;
	double d_pconst_dts;
	int    *d_kmt;
	double *d_vit;
	double *d_dzp;
	double *d_at;
	double *d_tend;
	double *d_dt_conv;
	double *d_atb;
	double *d_c;
#ifdef LOWRES
	float  *d_icmon;
#endif

}Para;


//dynamic atomic definition
#define CP_CORE      64
#define RCV_CORE     0
#define DATAPAKAGE   1

//ATHREAD global definition 
__thread_local int warp_size = 32; 
__thread_local int blockIndex;
__thread_local Para* Para_s;

 void  convadj_cu_(Para *Paras){	

  //ATHREAD local variable definition
     int i,j,k;
  int swGrid[3];
  int swBlock[3];
  int stride;
  int blockId_x,blockId_y,blockId_z;
  int threadId_x,threadId_y,threadId_z;


  //ATHREAD basic operation
  blockIndex = CRTS_tid;
  Para_s = (Para*)ldm_malloc(sizeof(Para));
  CRTS_dma_iget(Para_s, Paras, sizeof(Para), &dma_rply);
  D_COUNT++;
  CRTS_dma_wait_value(&dma_rply, D_COUNT);
  //printf("

  //double* d_pconst_dts_s = (double*)ldm_malloc(sizeof(double));
  //CRTS_dma_iget(d_pconst_dts_s, Para_s->d_pconst_dts, sizeof(double), &dma_rply);
  //D_COUNT++;
  //CRTS_dma_wait_value(&dma_rply, D_COUNT);

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

  double *d_at_s = Para_s->d_at;
  double *d_tend_s = Para_s->d_tend;
  double *d_dt_conv_s = Para_s->d_dt_conv;
  double *d_atb_s = Para_s->d_atb;
#if (defined LOWRES)
  float  *d_icmon_s = Para_s->d_icmon;
#endif
  for(int block = blockStart; block < blockEnd; block++){
    blockId_z = block / (Para_s->cuGrid[0] * Para_s->cuGrid[1]);
    blockId_y = (block - blockId_z * (Para_s->cuGrid[0] * Para_s->cuGrid[1])) / Para_s->cuGrid[0];
    blockId_x = block % Para_s->cuGrid[0];
      j = blockId_y * Para_s->cuBlock[1];


    //ATHREAD serial operation replace GPU parallel mutli-thread in single thread-block
    for(i = 0; i < threadSum; i++) {

      double rhoup[km] = {0.0};	
	
      double rholo[km]= {0.0};	
	
      double trasum[2];	
	
      double tup,sup,tlo,slo,dztsum,tramix;	
	
      int kcon,lctot,lcven,l1,l,lcon,lcona,lconb,lmix,n;	
	
      double c2dtts;	
	
      double tmpt1,tmpt2;	

      if(i < imt && j<jmt)	
      {	
      	kcon=Para_s->d_kmt[j *imt+i];	
      	lctot=0;	
      	lcven=0;	

        if(kcon!=0)	
        {	
      	    lcven=1;	
      	    lcon=0;	

            for(l=0;l<km-1; l++)	
	    {
      	        l1=l+1;	
      	        tup=d_at_s[l1 *jmt *imt+j *imt+i]-Para_s->d_to[l1];	
      	        sup=d_at_s[1 *km *jmt *imt+l1 *jmt *imt+j *imt+i]-Para_s->d_so[l1];	
      	        tlo=d_at_s[l *jmt *imt+j *imt+i]-Para_s->d_to[l1];	
      	        slo=d_at_s[1 *km *jmt *imt+l *jmt *imt+j *imt+i]-Para_s->d_so[l1];	
      	        rhoup[l1]=(Para_s->d_c[l1]+(Para_s->d_c[3 *km+l1]+Para_s->d_c[6 *km+l1] *sup) *sup+(Para_s->d_c[2 *km+l1]+Para_s->d_c[7 *km+l1] *sup+Para_s->d_c[5 *km+l1] *tup) *tup) *tup+(Para_s->d_c[1 *km+l1]+(Para_s->d_c[4 *km+l1]+Para_s->d_c[8 *km+l1] *sup) *sup) *sup;	
      	        rholo[l]=(Para_s->d_c[l1]+(Para_s->d_c[3 *km+l1]+Para_s->d_c[6 *km+l1] *slo) *slo+(Para_s->d_c[2 *km+l1]+Para_s->d_c[7 *km+l1] *slo+Para_s->d_c[5 *km+l1] *tlo) *tlo) *tlo+(Para_s->d_c[1 *km+l1]+(Para_s->d_c[4 *km+l1]+Para_s->d_c[8 *km+l1] *slo) *slo) *slo;	
	    }	

            for(k=kcon-1;k>=1; k--)	
	    {	
                if(rholo[k-1]> rhoup[k])	
			lcon = k;	
	    }

            if(lcon!=0)	
            {	
                for(;;)	
		{	
      	            lcona=lcon;	
      	            lconb=lcon+1;	
      	            dztsum=Para_s->d_dzp[lcona-1]+Para_s->d_dzp[lconb-1];	

                    for(n=0;n<2; n++)	
	            {	
      	                trasum[n]=d_at_s[n *km *jmt *imt+(lcona-1) *jmt *imt+j *imt+i] *Para_s->d_dzp[lcona-1]+d_at_s[n *km *jmt *imt+(lconb-1) *jmt *imt+j *imt+i] *Para_s->d_dzp[lconb-1];	
      	                tramix=trasum[n] /dztsum;	
      	                d_at_s[n *km *jmt *imt+(lcona-1) *jmt *imt+j *imt+i]=tramix;	
      	                d_at_s[n *km *jmt *imt+(lconb-1) *jmt *imt+j *imt+i]=tramix;	
	            }	

                    for(;;)	
		    {	
                        if(lconb!=kcon)	
                        {	
                        	l1=lconb+1;	
                        	tmpt1=d_at_s[(lconb-1) *jmt *imt+j *imt+i]-Para_s->d_to[l1-1];	
                        	tmpt2=d_at_s[1 *km *jmt *imt+(lconb-1) *jmt *imt+j *imt+i]-Para_s->d_so[l1-1];	
                        	rholo[lconb-1]=(Para_s->d_c[l1-1]+(Para_s->d_c[3 *km+l1-1]+Para_s->d_c[6 *km+l1-1] *tmpt2) *tmpt2+(Para_s->d_c[2 *km+l1-1]+Para_s->d_c[7 *km+l1-1] *tmpt2+Para_s->d_c[5 *km+l1-1] *tmpt1) *tmpt1) *tmpt1+(Para_s->d_c[1 *km+l1-1]+(Para_s->d_c[4 *km+l1-1]+Para_s->d_c[8 *km+l1-1] *tmpt2) *tmpt2) *tmpt2;	

                                if(rholo[lconb-1]> rhoup[l1-1])	
                                {	
                                	lconb=lconb+1;	
                                	dztsum=dztsum+Para_s->d_dzp[lconb-1];	

                                	for(n=0;n<2; n++)	
                                  	{	
                        			trasum[n]=trasum[n]+d_at_s[n *km *jmt *imt+(lconb-1) *jmt *imt+j *imt+i] *Para_s->d_dzp[lconb-1];	
                        			tramix=trasum[n] /dztsum;	

      						for(lmix=lcona-1;lmix<lconb; lmix++)	
						{	
      						    d_at_s[n *km *jmt *imt+lmix *jmt *imt+j *imt+i]=tramix;	
						}	
					}

      					goto CONV_2;	
      				}	
      			}

      			if(lcona> 1)	
                        {	
                        	l1=lcona-1;	
                        	tmpt1=d_at_s[(l1-1) *jmt *imt+j *imt+i]-Para_s->d_to[lcona-1];	
                        	tmpt2=d_at_s[1 *km *jmt *imt+(l1-1) *jmt *imt+j *imt+i]-Para_s->d_so[lcona-1];	
                        	rholo[l1-1]=(Para_s->d_c[lcona-1]+(Para_s->d_c[3 *km+lcona-1]+Para_s->d_c[6 *km+lcona-1] *tmpt2) *tmpt2+(Para_s->d_c[2 *km+lcona-1]+Para_s->d_c[7 *km+lcona-1] *tmpt2+Para_s->d_c[5 *km+lcona-1] *tmpt1) *tmpt1) *tmpt1+(Para_s->d_c[1 *km+lcona-1]+(Para_s->d_c[4 *km+lcona-1]+Para_s->d_c[8 *km+lcona-1] *tmpt2) *tmpt2) *tmpt2;	
                        	tmpt1=d_at_s[(lcona-1) *jmt *imt+j *imt+i]-Para_s->d_to[lcona-1];	
                        	tmpt2=d_at_s[1 *km *jmt *imt+(lcona-1) *jmt *imt+j *imt+i]-Para_s->d_so[lcona-1];	
                        	rhoup[lcona-1]=(Para_s->d_c[lcona-1]+(Para_s->d_c[3 *km+lcona-1]+Para_s->d_c[6 *km+lcona-1] *tmpt2) *tmpt2+(Para_s->d_c[2 *km+lcona-1]+Para_s->d_c[7 *km+lcona-1] *tmpt2+Para_s->d_c[5 *km+lcona-1] *tmpt1) *tmpt1) *tmpt1+(Para_s->d_c[1 *km+lcona-1]+(Para_s->d_c[4 *km+lcona-1]+Para_s->d_c[8 *km+lcona-1] *tmpt2) *tmpt2) *tmpt2;	

                                if(rholo[lcona-2]> rhoup[lcona-1])	
      				{	
      				    lcona=lcona-1;	
      				    dztsum=dztsum+Para_s->d_dzp[lcona-1];	

        			    for(n=0;n<2; n++)	
				    {	
      				        trasum[n]=trasum[n]+d_at_s[n *km *jmt *imt+(lcona-1) *jmt *imt+j *imt+i] *Para_s->d_dzp[lcona-1];	
      				        tramix=trasum[n] /dztsum;	

				        for(lmix=lcona;lmix<=lconb; lmix++)	
					{	
				      	    d_at_s[n *km *jmt *imt+(lmix-1) *jmt *imt+j *imt+i]=tramix;	
					}	
				    }	
      				    goto CONV_2;	
      				}	
      			}

      			goto EXIT_CONV_2;
      			CONV_2:;	
		}	
        	EXIT_CONV_2:;

      		lctot=lctot+lconb-lcona+1;	

      		if(lcona==1)	
		    lcven=lconb-lcona+1;	
	
      		if(lconb==kcon)	
      		{	
#if (defined LOWRES)
      		    d_icmon_s[j *imt+i]+=lctot;	
      		    d_icmon_s[1 *jmt *imt+j *imt+i]+=lcven;	
#endif
      		    goto III;	
      		}

      		lcon=lconb;	

      		for(;;)
		{	
      		    lcon=lcon+1;	
      		    if(lcon==kcon)	
      		    {	
#if (defined LOWRES)
      			d_icmon_s[j *imt+i]+=lctot;	
      			d_icmon_s[j *imt+i]+=lctot;	
      			d_icmon_s[1 *jmt *imt+j *imt+i]+=lcven;	
#endif
      			goto III;	
      		    }	
      		    if(rholo[lcon-1]<=rhoup[lcon])	
		        goto CONV_3;	
	
      		    goto EXIT_CONV_3;	
      		    CONV_3:;
		}

      		EXIT_CONV_3:;
	     }
III:;	
      	   }	
      }	

      c2dtts=Para_s->d_pconst_dts;	

      for(n=0;n<ntra; n++)	
      {	
          for(k=0;k<km; k++)	
          {	
          	//d_dt_conv_s[n *km *jmt *imt+k *jmt *imt+j *imt+i]=(d_at_s[n *km *jmt *imt+k *jmt *imt+j *imt+i]-d_atb_s[n *(km+1) *jmt *imt+(k+1) *jmt *imt+j *imt+i]) /c2dtts *Para_s->d_vit[k *jmt *imt+j *imt+i];	
          	//d_tend_s[n *km *jmt *imt+k *jmt *imt+j *imt+i]+=d_dt_conv_s[n *km *jmt *imt+k *jmt *imt+j *imt+i];	
          	d_atb_s[n *(km+1) *jmt *imt+(k+1) *jmt *imt+j *imt+i]=d_at_s[n *km *jmt *imt+k *jmt *imt+j *imt+i];	
          }	
      }	
    }	

   }

  }
    //todo: user need modify data size according to the practical situation

    //ldm_free(d_pconst_dts_s, sizeof(double));

    ldm_free(Para_s, sizeof(Para));
    return;
}
