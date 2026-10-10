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
   int nx;
   int *start_row;
   int *start_col;
   int *offsite_row;
   int *offsite_col;
   int *start3d;
   int *offsite3d;
   int *kmt;
   double err_norm2;
   double gamma;
   double od0cp;
   double *vit;
   double *tf;
   double *odzp;
   double *atb;
   double *stf;
   double *swv;
   double *sst;
   double *nswv;
   double *tsf;
   double *dqdt;
   double *net;
   double *seaice;
}Para;

__thread_local Para* parasw;
//__thread_local double* dat;
void tracer6_cu_slave_(Para *parahost)
	
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

if(parasw->tracer_n == 1){   

#ifdef SSSNORM
	
     for (iblock = 0; iblock < nblocks_clinic;++iblock){
         for(j = startrow; j < endrow ; ++j){
           for(i = startcol; i < endcol ; ++i){
	      parasw->tf[j*imt + i] += parasw->err_norm2 * parasw->vit[j*imt + i];
	   }
	 }
       }
#endif
}else{
     for (iblock = 0; iblock < nblocks_clinic;++iblock){
         for(j = startrow; j < endrow ; ++j){
           for(i = startcol; i < endcol ; ++i){
		   if(parasw->kmt[j*imt + i] > 0){
#ifdef COUP
			parasw->stf[j*imt + i] = parasw->tsf[j*imt + i];		   
#else	   
#ifdef FRC_CORE
			parasw->stf[j*imt + i] = ((parasw->swv[j*imt + i] + parasw->nswv[j*imt + i])
					* parasw->od0cp + parasw->seaice[j*imt + i] * parasw->gamma
					* (parasw->sst[j*imt + i] - parasw->atb[1*jmt*imt + j*imt + i]) / parasw->odzp[0]);		  
#else
			parasw->stf[j*imt + i] = (parasw->swv[j*imt + i] + parasw->nswv[j*imt + i]
					- parasw->dqdt[j*imt + i] * (parasw->sst[j*imt + i] - parasw->atb[1*jmt*imt + j*imt + i])) * parasw->od0cp;		  
#endif
#endif
			parasw->tf[j*imt + i] += parasw->stf[j*imt + i] * parasw->odzp[0] *(1.0 - aidif);
			parasw->net[j*imt + i] = parasw->stf[j*imt + i] * parasw->odzp[0];
			
		   }
	   }
	 }
       }
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



