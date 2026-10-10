//  CVS: $Id: icesnow->F90,v 1->1->1->1 2004/04/29 06:22:39 lhl Exp $
//     ==================
#include "athread.h"
#include <crts.h>
#include "slave.h"
#include <simd.h>
#include "paramsw.h"
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
   int pad;
   double aa; 
   double c2dtts;
   double onbc;
   double oncc;
   int *start;
   int *offsite;
   int *start3d;
   int *offsite3d;
   double *h0f;
   int *d_kmt;
   double *stf;
   double *h0l;
   double *utf;
   double *vtf;
   double *utl;
   double *vtl;
   double *wkd;
   double *wkb;
}Para;

__thread_local Para* parasw;
//__thread_local double* dat;
void tracersw_(Para *parahost)
{
      int iblock,k,j,i;
      int nblocks_clinic=1;
      int start,end;
      int start3d,end3d;
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
      
      
      start=parasw->start[athread_tid];
      end =  start + parasw->offsite[athread_tid];
      start3d=parasw->start3d[athread_tid];
      end3d =  start3d + parasw->offsite3d[athread_tid];

 //     dat = (double*)ldm_malloc(sizeof(double)*(end-start));
      for (iblock = 0; iblock < nblocks_clinic;++iblock){
          for (i = start; i < end; ++i) {
  	      parasw->h0f[i] *= parasw->onbc;
	      parasw->stf[i]  = parasw->aa * parasw->h0f[i]
		      +(1.0 - parasw->aa) * parasw->h0l[i];
	  }
      }

      for (iblock = 0; iblock < nblocks_clinic;++iblock){
          for (i = start3d; i < end3d; ++i) {
  	      parasw->utf[i] *= parasw->oncc;
	      parasw->vtf[i] *= parasw->oncc; 
	  }
      }

      for (iblock = 0; iblock < nblocks_clinic;++iblock){
          for (i = start3d; i < end3d; ++i) {
  	      parasw->wkd[i] = parasw->aa * parasw->utf[i]
		      +(1.0 - parasw->aa) * parasw->utl[i];
  	      parasw->wkb[i] = parasw->aa * parasw->vtf[i]
		      +(1.0 - parasw->aa) * parasw->vtl[i];
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



