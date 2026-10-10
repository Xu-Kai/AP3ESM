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
   int kvt;
   int simple_assm;
   int str_adv_tracer_flag;
   bool iau;
   double dts;
   double c2dtts;
   double d_lamda;
   int *start_row;
   int *start_col;
   int *offsite_row;
   int *offsite_col;
   int *start3d;
   int *offsite3d;
   double *tf;
   double *vit;
   double *atb;
   double *restore;
   double *restore_at;
   double *dt_restore;
//   double *ori;
   double *vtl;
   double *at;
	double *d_at_iau;
}Para;

__thread_local Para* parasw;
//__thread_local double* dat;
void tracer7_cu_slave_(Para *parahost)
	
{
      int iblock,k,j,i;
      int nblocks_clinic=1;
      int startrow,endrow,startcol,endcol;
      int start3d,end3d;
      double lamda,lamda1[km];
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

      lamda=1.0/(15.0*86400.0);
     for(k = 0; k < km && k < 45; ++k){
      lamda1[k] = 1.0 / (5.0 * 86400.0);
     }
     for(k = 45; k < km; ++k){
      lamda1[k] = 1.0 / (5.0 * 86400.0)/(k * k / 5.0 / 5.0);
     }

        if(athread_rid == 7) endrow-=2;
        if(athread_cid == 7) endcol-=2;
        if(athread_rid == 0) startrow+=2;
        if(athread_cid == 0) startcol+=2;
    
//       	if(parasw->simple_assm){
//       	  if(athread_tid == 0 && parasw->mytid == 0){
//		printf("into restoring, n = %d \n", parasw->tracer_n + 1);
//	  }
//
//
//     	  for (iblock = 0; iblock < nblocks_clinic;++iblock){
//            for(k = 1; k < km-1; k++){
//              for(j = startrow; j < endrow ; ++j){
//                for(i = startcol; i < endcol ; ++i){
//	          parasw->tf[k*jmt*imt + j*imt + i] += parasw->vit[k*jmt*imt + j*imt + i]
//			  * (parasw->restore_at[parasw->tracer_n*km*jmt*imt + k*jmt*imt + j*imt + i]
//			   - parasw->atb[parasw->tracer_n*(km+1)*jmt*imt + (k+1)*jmt*imt + j*imt + i])
//			  * lamda1[k];
//		  parasw->dt_restore[k*jmt*imt + j*imt + i] = parasw->vit[k*jmt*imt + j*imt + i]
//			  * (parasw->restore_at[parasw->tracer_n*km*jmt*imt + k*jmt*imt + j*imt + i]
//			   - parasw->atb[parasw->tracer_n*(km+1)*jmt*imt + (k+1)*jmt*imt + j*imt + i])
//			  * lamda1[k];
//		 }
//	      }
//	    }
//	  }
//     }

    if(parasw->iau){
	for (iblock = 0; iblock < nblocks_clinic;++iblock){
		for(k = 0; k < km; k++){
			for(j = startrow; j < endrow ; ++j){
				for(i = startcol; i < endcol ; ++i){
					parasw->tf[k*jmt*imt + j*imt + i] += parasw->vit[k*jmt*imt + j*imt + i]
					* parasw->d_at_iau[parasw->tracer_n*km*jmt*imt + k*jmt*imt + j*imt + i] * parasw->d_lamda;
				}
			}
		}
	}
	}

	if(parasw->boundary_restore == 1){
     	  for (iblock = 0; iblock < nblocks_clinic;++iblock){
            for(k = 1; k < km-1; k++){
              for(j = startrow; j < endrow ; ++j){
                for(i = startcol; i < endcol ; ++i){
	          parasw->tf[k*jmt*imt + j*imt + i] += parasw->vit[k*jmt*imt + j*imt + i]
			  * (parasw->restore[parasw->tracer_n*km*jmt*imt + k*jmt*imt + j*imt + i]
			   - parasw->atb[parasw->tracer_n*(km+1)*jmt*imt + (k+1)*jmt*imt + j*imt + i])
			  * lamda;
		   }
	      }
	    }
	  }
	  
     	  for (iblock = 0; iblock < nblocks_clinic;++iblock){
            for(k = 1; k < km-1; k++){
              for(j = startrow; j < endrow ; ++j){
                for(i = startcol; i < endcol ; ++i){
	          parasw->tf[k*jmt*imt + j*imt + i] += parasw->vit[k*jmt*imt + j*imt + i]
			  * (parasw->restore[parasw->tracer_n*km*jmt*imt + k*jmt*imt + j*imt + i]
			   - parasw->atb[parasw->tracer_n*(km+1)*jmt*imt + (k+1)*jmt*imt + j*imt + i])
			  * lamda;
		   }
	      }
	    }
	  }
	
	}
        if(parasw->str_adv_tracer_flag == 1){
     	  for (iblock = 0; iblock < nblocks_clinic;++iblock){
            for(k = 0; k < km; ++k){
              for(j = startrow; j < endrow ; ++j){
                for(i = startcol; i < endcol ; ++i){
	          parasw->vtl[k*jmt*imt + j*imt + i] = parasw->at[parasw->tracer_n*km*jmt*imt + k*jmt*imt + j*imt + i]
			   + parasw->dts * parasw->tf[k*jmt*imt + j*imt + i];
		   }
		}
	      }
	    }
	  }else if(parasw->str_adv_tracer_flag == 0){
     	    for (iblock = 0; iblock < nblocks_clinic;++iblock){
              for(k = 0; k < km; ++k){
                for(j = startrow; j < endrow ; ++j){
                  for(i = startcol; i < endcol ; ++i){
		     parasw->vtl[k*jmt*imt + j*imt + i] = parasw->atb[parasw->tracer_n*(km+1)*jmt*imt + k*jmt*imt + j*imt + i]
			   + parasw->c2dtts * parasw->tf[k*jmt*imt + j*imt + i];
		   }
		}
	      }
	    }
	  }
/*	  
        if(athread_rid == 7) endrow+=2;
        if(athread_cid == 7) endcol+=2;
        if(athread_rid == 0) startrow-=2;
        if(athread_cid == 0) startcol-=2;
	  
     	    for (iblock = 0; iblock < nblocks_clinic;++iblock){
              for(k = 0; k < km; ++k){
                for(j = startrow; j < endrow ; ++j){
                  for(i = startcol; i < endcol ; ++i){
	            parasw->ori[k*jmt*imt + j*imt + i] = parasw->vtl[k*jmt*imt + j*imt + i];
		   }
		}
	      }
	    }
   */     
     
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



