//  CVS: $Id: icesnow->F90,v 1->1->1->1 2004/04/29 06:22:39 lhl Exp $
//     ==================
#include "athread.h"
#include <crts.h>
#include "slave.h"
#include <simd.h>
#include "paramsw.h"
#include "def-undef.h"
#include "math.h"
#include <stdio.h>
#include <float.h>
#include <string.h>
#include <stdbool.h>
#include "swperf.h"
#include <time.h>

//dynamic atomic definition

#define CP_CORE      64
#define RCV_CORE     0
#define DATAPAKAGE   1

//#define min(x,y) (x)<(y)?(x):(y)
#define index2D(j,i) j*imt + i
#define index3D(k,j,i) k*jmt*imt + j*imt + i

//ATHREAD DMA transfer definition 
__thread_local unsigned int D_COUNT = 0;
__thread_local crts_rply_t dma_rply = 0;
__thread_local crts_rply_t l_rply = 0;
__thread_local crts_rply_t r_rply = 0;

typedef struct
{
   int mytid;
   int tracer_n;
   int boundary_restore;
   int *start_row;
   int *start_col;
   int *offsite_row;
   int *offsite_col;
   int *start3d;
   int *offsite3d;
   int *kmt;
   double od0;
   double gamma;
   double *swv;
#ifdef SOLAR
   double *pen;
#endif
#ifdef SOLARCHLORO
   double *pen_chl;
#endif
   double *vit;
   double *tf;
   double *penetrate;
   double *odzp;
   double *odzt;
   double *wkc;
   double *atb;
   double *dz;
   double *stf;
   double *ssf;
   double *sss;
   double *ifrac;
   double *net;
   double *tarea;
   double *dt_restore;
   double *fresh;
   double *seaice;
#ifdef SSS_restore
   double *sss_restore;
#endif
}Para;

__thread_local Para* parasw;
//__thread_local double* dat;
void tracer5_cu_slave_(Para *parahost)
	
{
      int iblock,k,j,i;
      int nblocks_clinic=1;
      int startrow,endrow,startcol,endcol;
      int start3d,end3d;
      double aidif = 0.5;
#ifdef TRACER_TIME
      double t1,t2,t3,t4,t5,t6,t7;
      athread_ssync_array();
      t1=athread_stime_cycle();
#endif
      parasw = (Para*)ldm_malloc(sizeof(Para));
#ifdef TRACER_TIME
      athread_ssync_array();
      t2=athread_stime_cycle();
#endif
      athread_dma_iget(parasw,parahost,sizeof(Para),&dma_rply);
      D_COUNT++;
      athread_dma_wait_value(&dma_rply,D_COUNT);
      
      
      startrow=parasw->start_row[athread_rid];
      endrow =  startrow + parasw->offsite_row[athread_rid];
      startcol=parasw->start_col[athread_cid];
      endcol =  startcol + parasw->offsite_col[athread_cid];

      start3d=parasw->start3d[athread_tid];
      end3d =  start3d + parasw->offsite3d[athread_tid];
      
      if(athread_rid == 7) endrow-=2;
      if(athread_cid == 7) endcol-=2;
      if(athread_rid == 0) startrow+=2;
      if(athread_cid == 0) startcol+=2;

      if(parasw->tracer_n == 0){
#ifdef SOLAR
	
     for (iblock = 0; iblock < nblocks_clinic;++iblock){
       for(k = 1; k < km-1; k++){
         for(j = startrow; j < endrow ; ++j){
           for(i = startcol; i < endcol ; ++i){
	      const double wt1 = parasw->swv[j*imt+i] 
		      * parasw->pen[k-1] * parasw->vit[k*jmt*imt + j*imt + i];
	      const double wt2 = parasw->swv[j*imt+i] 
		      * parasw->pen[k  ] * parasw->vit[(k+1)*jmt*imt + j*imt + i];
	      parasw->tf[k*jmt*imt + j*imt + i] += (wt1 - wt2) * parasw->odzp[k];
	      parasw->penetrate[k*jmt*imt + j*imt + i] = (wt1 - wt2) * parasw->odzp[k];
	   }
	 }
       }
     }

     for (iblock = 0; iblock < nblocks_clinic;++iblock){
         for(j = startrow; j < endrow ; ++j){
           for(i = startcol; i < endcol ; ++i){
	      const double wt1 = parasw->swv[j*imt +i] 
		      * parasw->pen[0] * parasw->vit[1*jmt*imt + j*imt + i];
	      const double wt2 = parasw->swv[j*imt + i] 
		      * parasw->pen[km - 2] * parasw->vit[(km-1)*jmt*imt + j*imt + i];
	      
	      parasw->tf[j*imt + i] -= wt1  * parasw->odzp[0];
	      parasw->tf[(km-1)*jmt*imt + j*imt + i] += wt2  * parasw->odzp[km-1];
	      
	      parasw->penetrate[j*imt + i] = -wt1 * parasw->odzp[0];
	      parasw->penetrate[(km-1)*jmt*imt + j*imt + i] =  wt2 * parasw->odzp[km-1];
	   }
	 }
       }
     
#endif

#ifdef SOLARCHLORO
     for (iblock = 0; iblock < nblocks_clinic;++iblock){
       for(k = 1; k < km-1; k++){
         for(j = startrow; j < endrow ; ++j){
           for(i = startcol; i < endcol ; ++i){
	       double wt1 = parasw->swv[index2D(j,i)] 
		      * parasw->pen_chl[(k-1)*jmt*imt + j*imt + i] * parasw->vit[k*jmt*imt + j*imt + i];
	       double wt2 = parasw->swv[index2D(j,i)] 
		      * parasw->pen_chl[k*jmt*imt + j*imt + i] * parasw->vit[(k+1)*jmt*imt + j*imt + i];
	      parasw->tf[k*jmt*imt + j*imt + i] += (wt1 - wt2) * parasw->odzp[k];
	      parasw->penetrate[k*jmt*imt + j*imt + i] = (wt1 - wt2) * parasw->odzp[k];
	   }
	 }
       }
     }

     for (iblock = 0; iblock < nblocks_clinic;++iblock){
         for(j = startrow; j < endrow ; ++j){
           for(i = startcol; i < endcol ; ++i){
	       double wt1 = parasw->swv[index2D(j,i)] 
		      * parasw->pen_chl[j*imt + i] * parasw->vit[1*jmt*imt + j*imt + i];
	       double wt2 = parasw->swv[index2D(j,i)] 
		      * parasw->pen_chl[(km-2)*jmt*imt + j*imt + i] * parasw->vit[(km-1)*jmt*imt + j*imt + i];
	      
	      parasw->tf[j*imt + i] -= parasw->odzp[0] * wt1;
	      parasw->tf[(km-1)*jmt*imt + j*imt + i] += parasw->odzp[km-1]*wt2;
	      
	      parasw->penetrate[j*imt + i] = -wt1 * parasw->odzp[0];
	      parasw->penetrate[(km-1)*jmt*imt + j*imt + i] =  wt2 * parasw->odzp[km-1];
	   }
	 }
       }
#endif
}
     for (iblock = 0; iblock < nblocks_clinic;++iblock){
       for(k = 1; k < km-1; k++){
         for(j = startrow; j < endrow ; ++j){
           for(i = startcol; i < endcol ; ++i){
	       double wt1 = parasw->wkc[(k-1)*jmt*imt + j*imt + i] 
		      * (parasw->atb[parasw->tracer_n*(km+1)*jmt*imt + k*jmt*imt + j*imt + i]
		       - parasw->atb[parasw->tracer_n*(km+1)*jmt*imt + (k+1)*jmt*imt + j*imt + i])
		      *  parasw->odzt[k] *parasw->vit[k*jmt*imt + j*imt + i];
	       double wt2 = parasw->wkc[k*jmt*imt + j*imt + i] 
		      * (parasw->atb[parasw->tracer_n*(km+1)*jmt*imt + (k+1)*jmt*imt + j*imt + i]
		       - parasw->atb[parasw->tracer_n*(km+1)*jmt*imt + (k+2)*jmt*imt + j*imt + i])
		      *  parasw->odzt[k+1] *parasw->vit[(k+1)*jmt*imt + j*imt + i];
	      
	      
	      parasw->tf[k*jmt*imt + j*imt + i] += parasw->odzp[k] * (wt1 - wt2) *  (1.0 - aidif);
	    //  parasw->dz[parasw->tracer_n*km*jmt*imt + k*jmt*imt + j*imt + i] = parasw->odzp[k] * (wt1 - wt2) *  (1.0 -aidif);
	   }
	 }
       }
     }

     for (iblock = 0; iblock < nblocks_clinic;++iblock){
         for(j = startrow; j < endrow ; ++j){
           for(i = startcol; i < endcol ; ++i){
	      const double wt1 = parasw->wkc[j*imt + i] 
		      * (parasw->atb[parasw->tracer_n*(km+1)*jmt*imt + 1*jmt*imt + j*imt + i]
		       - parasw->atb[parasw->tracer_n*(km+1)*jmt*imt + 2*jmt*imt + j*imt + i])
		      *  parasw->odzt[1] *parasw->vit[1*jmt*imt + j*imt + i];
	      const double wt2 = parasw->wkc[(km-2)*jmt*imt + j*imt + i] 
		      * (parasw->atb[parasw->tracer_n*(km+1)*jmt*imt + (km-1)*jmt*imt + j*imt + i]
		       - parasw->atb[parasw->tracer_n*(km+1)*jmt*imt +  km*jmt*imt + j*imt + i])
		      *  parasw->odzt[km - 1] *parasw->vit[(km-1)*jmt*imt + j*imt + i];
	      
	      
	      parasw->tf[j*imt + i] -= parasw->odzp[0] * wt1 * (1.0 - aidif);
	      parasw->tf[(km-1)*jmt*imt + j*imt + i] += parasw->odzp[km-1] * wt2 * (1.0 - aidif);
	      
	    //  parasw->dz[parasw->tracer_n*km*jmt*imt + j*imt + i] = - parasw->odzp[0] * wt1 * (1.0 -aidif);
	    //  parasw->dz[parasw->tracer_n*km*jmt*imt + (km-1)*jmt*imt + j*imt + i] =  parasw->odzp[km-1] * wt2 * (1.0 -aidif);
	   }
	 }
       }
      if(athread_rid == 7) endrow+=2;
      if(athread_cid == 7) endcol+=2;
      if(athread_rid == 0) startrow-=2;
      if(athread_cid == 0) startcol-=2;

      if(parasw->tracer_n == 1){
       for (iblock = 0; iblock < nblocks_clinic;++iblock){
         for(j = startrow; j < endrow ; ++j){
           for(i = startcol; i < endcol ; ++i){
	      if(parasw->kmt[j*imt + i] > 0 ){
#ifdef COUP
//		 if(parasw->boundary_restore == 2) {
//		    parasw->stf[j*imt + i] = parasw->ssf[j*imt + i] 
//			    / parasw->odzp[0] + (parasw->gamma * 10)
//			    * (parasw->sss[j*imt + i] - parasw->atb[(km+1)*jmt*imt + jmt*imt + j*imt + i])
//			    * parasw->ifrac[j*imt + i]
//			    / parasw->odzp[0] + (parasw->gamma)*(1.0-parasw->ifrac[j*imt + i])
//			    * (parasw->sss[j*imt + i] - parasw->atb[(km+1)*jmt*imt + jmt*imt + j*imt + i])
//			    / parasw->odzp[0];
//		 }else{
		   parasw->stf[j*imt + i] = parasw->ssf[j*imt + i] / parasw->odzp[0];
//		 }
#else
#ifdef FRC_CORE
//		parasw->stf[j*imt + i] = (parasw->fresh[j*imt + i] * 34.7 * parasw->od0 * 1.0e-3 
//				+ parasw->gamma
//			            * (parasw->sss[j*imt + i] - parasw->atb[(km+1)*jmt*imt + jmt*imt + j*imt + i])
//				         * parasw->seaice[j*imt + i] / parasw->odzp[0]
//				+ parasw->gamma * 30.0 / 365.0 / 4.0 * 50.0
//			            * (parasw->sss[j*imt + i] - parasw->atb[(km+1)*jmt*imt + jmt*imt + j*imt + i])
//				         / parasw->odzp[0] * (1.0 - parasw->seaice[j*imt + i]));
      parasw->stf[j*imt + i] = parasw->fresh[j*imt + i] * 34.7 * parasw->od0 * 1.0e-3 ;
#else
		parasw->stf[j*imt + i] = parasw->gamma
			            * (parasw->sss[j*imt + i] - parasw->atb[(km+1)*jmt*imt + jmt*imt + j*imt + i])
			            / parasw->odzp[0]; 
#endif
#endif
		parasw->tf[j*imt + i] += parasw->stf[j*imt + i] * (1.0 -aidif) * parasw->odzp[0];
		parasw->net[jmt*imt + j*imt + i] = parasw->stf[j*imt + i] * parasw->odzp[0];
	      }
	   }
	 }
       }
#ifdef SSS_restore
       for (iblock = 0; iblock < nblocks_clinic;++iblock){
         for(j = startrow; j < endrow ; ++j){
          for(i = startcol; i < endcol ; ++i){
	      if(parasw->kmt[j*imt + i] > 0 ){
#ifdef COUP
          parasw->sss_restore[j*imt + i] = (parasw->gamma * 10)
			    * (parasw->sss[j*imt + i] - parasw->atb[(km+1)*jmt*imt + jmt*imt + j*imt + i])
			    * parasw->ifrac[j*imt + i]
			    / parasw->odzp[0] + (parasw->gamma)*(1.0-parasw->ifrac[j*imt + i])
			    * (parasw->sss[j*imt + i] - parasw->atb[(km+1)*jmt*imt + jmt*imt + j*imt + i])
			    / parasw->odzp[0];
#else
          parasw->sss_restore[j*imt + i] = (7.0/3.0)*parasw->gamma*(parasw->sss[j*imt + i]-parasw->atb[(km+1)*jmt*imt + jmt*imt + j*imt + i])
          *parasw->seaice[j*imt + i]+ parasw->gamma*(parasw->sss[j*imt + i]-parasw->atb[(km+1)*jmt*imt + jmt*imt + j*imt + i])
                                          *(1.0-parasw->seaice[j*imt + i]);
#endif
        }
        parasw->tf[j*imt + i] += parasw->sss_restore[j*imt + i];
          }
         }
       }
#endif
      }

#ifdef TRACER_TIME
      t3=athread_stime_cycle();
#endif
#ifdef TRACER_TIME
      t4=athread_stime_cycle();
#endif

#ifdef TRACER_TIME
      t5=athread_stime_cycle();
#endif
#ifdef TRACER_TIME
      t6=athread_stime_cycle();
#endif

ldm_free(parasw,sizeof(Para));
#ifdef TRACER_TIME
	t7=athread_stime_cycle();
#endif

#ifdef TRACER_TIME
if(parasw->mytid == 0 && athread_tid == 0){
      printf("time 1 : %f\n",t2-t1);
      printf("time 2 : %f\n",t3-t2);
      printf("time 3 : %f\n",t4-t3);
      printf("time 4 : %f\n",t5-t4);
      printf("time 5 : %f\n",t6-t5);
      printf("time 6 : %f\n",t7-t6);
}
#endif

return;
 
}



