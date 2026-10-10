//  CVS: $Id: icesnow->F90,v 1->1->1->1 2004/04/29 06:22:39 lhl Exp $
//     ==================
#include <crts.h>
#include "athread.h"
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
#undef TRACER_TIME 
//ATHREAD DMA transfer definition 
__thread_local unsigned int D_COUNT = 0;
__thread_local crts_rply_t dma_rply = 0;
__thread_local crts_rply_t l_rply = 0;
__thread_local crts_rply_t r_rply = 0;

typedef struct
{
   int mytid;
   double c2dtts;
   int *start_row;
   int *start_col;
   int *offsite_row;
   int *offsite_col;
   int *start3d;
   int *offsite3d;
   int *kmt;
   double *vit;
   double *wk;
   double *topbc;
   double *dcb;
   double *odzp;
   double *odzt;
}Para;

__thread_local Para* parasw;
//__thread_local double* dat;
void invtritsw_(Para *parahost)
{
      int iblock,k,j,i;
      int ind2d,ind3d;
      int nblocks_clinic=1;
      int startrow,endrow,startcol,endcol;
      int start3d,end3d;
      double a8[km], b8[km], c8[km], d8[km];
      double e8[km+1], f8[km+1];
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
	  
      for (iblock = 0; iblock < nblocks_clinic;++iblock){
        for(j = startrow; j < endrow ; ++j){
          for(i = startcol; i < endcol ; ++i){
	    if(parasw->kmt[j*imt + i] > 0){
	      const int kz = parasw->kmt[j*imt + i] - 1;

	      for(k = 1; k <= kz ; ++k){
	        a8[k] = parasw->dcb[(k-1)*jmt*imt + j*imt + i] * parasw->odzt[k]
			* parasw->odzp[k] * parasw->c2dtts * aidif;
	        d8[k] = parasw->wk[k*jmt*imt + j*imt + i]; 
	      }

	      for(k = 1; k <= kz-1; ++k){
	        c8[k] = parasw->dcb[k*jmt*imt + j*imt + i] * parasw->odzt[k+1]
			* parasw->odzp[k] * parasw->c2dtts * aidif;
	        b8[k] = 1.0 + a8[k] + c8[k]; 
		e8[k] = 0.0;
		f8[k] = 0.0;
	      }
                a8[0] = parasw->odzp[0] * parasw->c2dtts * aidif;
                c8[0] = parasw->dcb[j*imt + i] * parasw->odzt[1] * parasw->odzp[0] * parasw->c2dtts * aidif;
                b8[0] = 1.0 + c8[0];
                d8[0] = parasw->wk[j*imt + i];
                e8[0] = 0.0;
                f8[0] = 0.0;

                b8[kz] = 1.0 + a8[kz];
                c8[kz] = parasw->odzp[kz] * parasw->c2dtts * aidif;

                e8[kz+1] = 0.0;
                f8[kz+1] = 0.0;

              for(k = kz; k >=0; --k){
	        const double g0 = 1.0 / (b8[k] - c8[k] * e8[k+1]);
	        e8[k] = a8[k] * g0;
		f8[k] = (d8[k] + c8[k] * f8[k+1]) * g0;
	      }

	      parasw->wk[j*imt + i] = (e8[0] * parasw->topbc[j*imt + i] + f8[0])
	                          * parasw->vit[j*imt + i];
	  
	      for (k = 1; k <= kz; ++k) {
	        parasw->wk[k*jmt*imt + j*imt + i] = (e8[k] * parasw->wk[(k-1)*jmt*imt + j*imt + i]
	                          + f8[k]) * parasw->vit[k*jmt*imt + j*imt + i];

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
      athread_ssync_array();
if(parasw->mytid == 0 && athread_tid == 0){
      printf("upwell time 1 : %f\n",t2-t1);
      printf("upwell time 2 : %f\n",t3-t2);
      printf("upwell time 3 : %f\n",t4-t3);
      printf("upwell time 4 : %f\n",t5-t4);
      printf("upwell time 5 : %f\n",t6-t5);
      printf("upwell time 6 : %f\n",t7-t6);
}
#endif

return;
 
}



