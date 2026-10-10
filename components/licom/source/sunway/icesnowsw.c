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

//ATHREAD DMA transfer definition 
__thread_local unsigned int D_COUNT = 0;
__thread_local crts_rply_t dma_rply = 0;
__thread_local crts_rply_t l_rply = 0;
__thread_local crts_rply_t r_rply = 0;

typedef struct
{
   int mytid;
   int pad;
   double d_SHR_CONST_LATICE; 
   double d_tbice;
   double d_cp;
   int *start;
   int *offsite;
   int *start3d;
   int *offsite3d;
   double *d_at;
   int *d_kmt;
   double *d_licomqice;
   double *d_dzp;
}Para;

__thread_local Para* parasw;
__thread_local double* dat;
void icesnowsw_(Para *parahost)
{
//     ==================
//     Sea Ice Model


      double  sal_ice, sal_ocn, tdiff, heat_ice_fusion;
      int  i, k, ind3,ind2, xy;
      int start,end;
      int start3d,end3d;
      xy = imt * jmt;
      heat_ice_fusion = SHR_CONST_LATICE_;         // J/Kg
      sal_ocn = 34.7e0 ;                    //Reference salinity for ocean
      sal_ice  = 4.0e0;                     //Reference salinity for sea ice
//------------------------------------------------------------
//  if SST exceed -1->8C to restore it to -1->8C
//------------------------------------------------------------
#ifdef ICESNOW_TIME
      double t1,t2,t3,t4,t5,t6,t7;
      athread_ssync_array();
      t1=athread_stime_cycle();
#endif

      parasw = (Para*)ldm_malloc(sizeof(Para));

#ifdef ICESNOW_TIME
      athread_ssync_array();
      t2=athread_stime_cycle();
#endif
      athread_dma_iget(parasw,parahost,sizeof(Para),&dma_rply);
      D_COUNT++;
      athread_dma_wait_value(&dma_rply,D_COUNT);
      
      
      start=parasw->start[athread_tid];
      end =  start + parasw->offsite[athread_tid];

      start3d =parasw->start3d[athread_tid];
      end3d =  start3d + parasw->offsite3d[athread_tid];
    /*  
        if(parasw->mytid ==15){
     	printf("mytid : %d ,theard tid: %d , start3d: %d , end3d : %d \n",parasw->mytid,athread_tid,start3d,end3d);
     	printf("km : %d ,jmt : %d , imt: %d  \n",km,jmt,imt);
	}
	*/
   
   
#ifdef ICESNOW_TIME
      t3=athread_stime_cycle();
#endif

      tdiff = 0.0e0;
 
      for (k=0; k<1; ++k) {
             for (i = start; i < end; ++i) {
                 if ( parasw->d_at[i] < parasw->d_tbice&& parasw->d_kmt[i]>0 ) {
                      tdiff = parasw->d_tbice - parasw->d_at[i];
                      parasw->d_licomqice[i] += tdiff * parasw->d_dzp[k]/parasw->d_dzp[0] ;
                      parasw->d_at[km*xy + i] +=  tdiff*(sal_ocn - sal_ice)* parasw->d_cp/heat_ice_fusion*0.001e0 ;
                      parasw->d_at[i] = parasw->d_tbice ;
			 
                  }
	 }
      }
#ifdef ICESNOW_TIME
      t4=athread_stime_cycle();
#endif
          for (i = start; i < end; ++i) {//same as j and jmt
                  if ( parasw->d_licomqice[i]>0.0 &&
		       parasw->d_at[i] > parasw->d_tbice ) {
                        tdiff = fmin((parasw->d_at[i] - parasw->d_tbice), parasw->d_licomqice[i]);
                        parasw->d_licomqice[i] -=  tdiff ;
                        parasw->d_at[i] -=  tdiff ;
	 		parasw->d_at[km*jmt*imt +i] -=  tdiff*(sal_ocn - sal_ice)* parasw->d_cp / heat_ice_fusion * 0.001e0 ;
                  }
      }
#ifdef ICESNOW_TIME
      t5=athread_stime_cycle();
#endif
  athread_ssync_array (); 
/*	  
      for (i = start3d; i < end3d; ++i) {
	  if(parasw->d_at[i] < parasw->d_tbice){
	     parasw->d_at[i] = parasw->d_tbice;
          }
      }
      */
      doublev8 va,vb;
      vb = parasw->d_tbice;
      if(athread_tid != 63){
      dat = (double*)ldm_malloc(sizeof(double)*(end3d-start3d));
      athread_dma_get (dat, &parasw->d_at[start3d], sizeof(double)*(end3d-start3d));
       for (i = 0; i < end3d-start3d; i+=8) {
          simd_load(va,&(dat[i]));
          va=simd_smaxd(va,vb);
          simd_store(va,&(dat[i]));
       }
      athread_dma_put (&parasw->d_at[start3d],dat, sizeof(double)*(end3d-start3d));
      ldm_free(dat,sizeof(double)*(end3d-start3d));
      }else{
       for (i = start3d; i < (end3d-(end3d%8)); i+=8) {
          simd_load(va,&(parasw->d_at[i]));
          va=simd_smaxd(va,vb);
          simd_store(va,&(parasw->d_at[i]));
	//  simd_print_doublev8(va);
       }
       for(i = (end3d-start3d)%8; i > 0 ;i--)
       {
 	if(parasw->d_at[end3d-start3d-i] < parasw->d_tbice){
	   parasw->d_at[end3d-start3d-i] = parasw->d_tbice;
	}
       }
       /*
      for (k = 0; k < km ; ++k) { 
 	      int a = k *jmt * imt;
       for (i = start; i < end; ++i) {
	 	if(parasw->d_at[a+i] < parasw->d_tbice){
		   parasw->d_at[a+i] = parasw->d_tbice;
		}
      	     }
     	 }
	 */
      }

#ifdef ICESNOW_TIME
      t6=athread_stime_cycle();
#endif

ldm_free(parasw,sizeof(Para));
#ifdef ICESNOW_TIME
	t7=athread_stime_cycle();
#endif

#ifdef ICESNOW_TIME
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



