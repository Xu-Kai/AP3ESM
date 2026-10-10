//  CVS: $Id: icesnow->F90,v 1->1->1->1 2004/04/29 06:22:39 lhl Exp $
//     ==================
#include <crts.h>
#include "athread.h"
#include "slave.h"
#include <simd.h>
//#include "paramsw.h"
#include "math.h"
#include <stdio.h>
#include <float.h>
#include <string.h>
#include <stdbool.h>
#include "swperf.h"
#include <time.h>
#include "param_mod.h"	

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
   int pad;
   double *dzp; 
   double *vit;
   double *ohbu;
   double *ohbt;
   int *start_row;
   int *start_col;
   int *offsite_row;
   int *offsite_col;
   int *start3d;
   int *offsite3d;
   int *kmt;
   double *ws;
   double *uk;
   double *vk;
   double *wka;
   double *work;
   double *uwk;
   double *vwk;
   double *h0wk;

   double *au0;
   double *aus;
   double *auw;
   double *ausw;

   double *hts;
   double *htw;
   double *tarea_r;
}Para;

__thread_local Para* parasw;
//__thread_local double* dat;
void upwellsw_(Para *parahost)
{
      int iblock,k,j,i;
      int ind2d,ind3d;
      int nblocks_clinic=1;
      int startrow,endrow,startcol,endcol;
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
      
      
      startrow=parasw->start_row[athread_rid];
      endrow =  startrow + parasw->offsite_row[athread_rid];
      startcol=parasw->start_col[athread_cid];
      endcol =  startcol + parasw->offsite_col[athread_cid];

      start3d=parasw->start3d[athread_tid];
      end3d =  start3d + parasw->offsite3d[athread_tid];
//  if(parasw->mytid == 0)     printf("mytid : %d ,my rid: %d , my cid : %d ,rowstart: %d , rowend : %d, colstart: %d , colend : %d \n",parasw->mytid,athread_rid,athread_cid,startrow,endrow,startcol,endcol);
      	//     dat = (double*)ldm_malloc(sizeof(double)*(end-start));
    for (iblock = 0; iblock < nblocks_clinic;++iblock){
      for(j = startrow; j < endrow ; ++j){
         for(i = startcol; i < endcol ; ++i){
  	      parasw->work[index2D(j,i)] = 0.0;
	 }
      }
    }
      for (iblock = 0; iblock < nblocks_clinic;++iblock){
          for (i = start3d; i < end3d; ++i) {
	      parasw->wka[i] = 0.0;
	  }
      }
      


     
#ifdef TRACER_TIME
      t3=athread_stime_cycle();
#endif
      if(athread_rid == 7) endrow--;
      if(athread_cid == 0) startcol++;

      for (iblock = 0; iblock < nblocks_clinic; ++iblock){
      for(j = startrow; j < endrow ; ++j){
         for(i = startcol; i < endcol ; ++i){
            parasw->work[index2D(j,i)] = 
		    parasw->au0[index2D(j,i)] * parasw->h0wk[index2D(j,i)] 
		  + parasw->aus[index2D(j,i)] * parasw->h0wk[(j+1)*imt + i] 
		  + parasw->auw[index2D(j,i)] * parasw->h0wk[j * imt + i-1] 
		  + parasw->ausw[index2D(j,i)]* parasw->h0wk[(j+1)*imt + i-1] ;
	 }	          	    
      }
   }

  for (iblock = 0; iblock < nblocks_clinic;++iblock){
    for(k = 0 ; k < km ; ++k){
      for(j = startrow; j < endrow ; ++j){
         for(i = startcol; i < endcol ; ++i){
	     ind3d = index3D(k,j,i);
	     ind2d = index2D(j,i);
 	     parasw->uk[ind3d] = (1.0 + parasw->work[ind2d] * parasw->ohbu[ind2d]) * parasw->uwk[ind3d];
 	     parasw->vk[ind3d] = (1.0 + parasw->work[ind2d] * parasw->ohbu[ind2d]) * parasw->vwk[ind3d];
	 }	          	    
      }
    }
   }
      if(athread_rid == 7) endrow++;
      if(athread_cid == 0) startcol--;
#ifdef TRACER_TIME
      t4=athread_stime_cycle();
#endif

      if(athread_rid == 0) startrow++;
      if(athread_cid == 7) endcol--;

      for(i = start3d; i < end3d ; ++i){
           parasw->wka[i] = 0.0;
       }
      flush_slave_cache ();
      athread_ssync_array();
     

  for (iblock = 0; iblock < nblocks_clinic;++iblock){
    for(k = 0 ; k < km ; ++k){
       int bid = 0;
 
      for(j = startrow; j < endrow ; ++j){
         for(i = startcol; i < endcol ; ++i){
            if( k <= (parasw->kmt[bid * imt * jmt + index2D(j,i)] -1) ){
	     parasw->wka[index3D(k,j,i)] = 0.500 * (
	         (parasw->uk[k*jmt*imt + j*imt + i+1]
	        + parasw->uk[k*jmt*imt + (j-1)*imt + i+1]) * parasw->htw[bid*imt*jmt + j*imt + i+1]
	       - (parasw->uk[k*jmt*imt + j*imt + i] 
	       +  parasw->uk[k*jmt*imt + (j-1)*imt + i]) * parasw->htw[bid * imt * jmt + j*imt + i]
	       + (parasw->vk[k*jmt*imt + j * imt + i+1] 
	       +  parasw->vk[k*jmt*imt + j * imt + i]) * parasw->hts[bid * imt * jmt + j*imt + i]
	       - (parasw->vk[k*jmt*imt + (j-1)*imt + i+1]
	       +  parasw->vk[k*jmt*imt + (j-1)*imt + i]) * parasw->hts[bid * imt * jmt + (j-1)*imt+ i ])
		     * parasw->tarea_r[bid*jmt*imt + j*imt + i];
	    }
	 }
      }	 
    }
   }  
      if(athread_rid == 0) startrow--;
      if(athread_cid == 7) endcol++;
#ifdef TRACER_TIME
      t5=athread_stime_cycle();
#endif
 
      for (iblock = 0; iblock < nblocks_clinic;++iblock){
       for(j = startrow; j < endrow ; ++j){
         for(i = startcol; i < endcol ; ++i){
  	      parasw->work[index2D(j,i)] = 0;
	 }
       }
     }

      if(athread_rid == 7) endrow--;
      if(athread_cid == 7) endcol--;
      if(athread_rid == 0) startrow++;
      if(athread_cid == 0) startcol++;

     for (iblock = 0; iblock < nblocks_clinic;++iblock){
       for(k = 0; k < km; k++){
         for(j = startrow; j < endrow ; ++j){
           for(i = startcol; i < endcol ; ++i){
  	      parasw->work[index2D(j,i)] = parasw->work[index2D(j,i)] - parasw->dzp[k] * parasw->wka[index3D(k,j,i)] * parasw->vit[index3D(k,j,i)];
	   }
	 }
       }
     }

     for (iblock = 0; iblock < nblocks_clinic;++iblock){
       for(k = 1; k < km; k++){
         for(j = startrow; j < endrow ; ++j){
           for(i = startcol; i < endcol ; ++i){
  	      parasw->ws[index3D(k,j,i)] = parasw->vit[index3D(k,j,i)] * 
		   (parasw->ws[(k-1)*jmt*imt + j*imt + i] + parasw->dzp[k-1] *
		   (parasw->work[index2D(j,i)] * parasw->ohbt[index2D(j,i)] + parasw->wka[(k-1)*jmt*imt + j*imt + i])); 
	      
	   }
	 }
       }
     }

     for (iblock = 0; iblock < nblocks_clinic;++iblock){
       for(j = startrow; j < endrow ; ++j){
         for(i = startcol; i < endcol ; ++i){
  	      parasw->work[index2D(j,i)] = 1.0/(1.0 + parasw->h0wk[index2D(j,i)] * parasw->ohbt[index2D(j,i)]);
	 }
       }
     }

     for (iblock = 0; iblock < nblocks_clinic;++iblock){
       for(k = 1; k < km; k++){
         for(j = startrow; j < endrow ; ++j){
           for(i = startcol; i < endcol ; ++i){
  	      parasw->ws[index3D(k,j,i)] = parasw->ws[index3D(k,j,i)] * parasw->work[index2D(j,i)]; 
	   }
	 }
       }
     }
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



