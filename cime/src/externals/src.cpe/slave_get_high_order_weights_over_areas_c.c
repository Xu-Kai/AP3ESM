#include <stdio.h>
#include <slave.h>
#include <crts.h>
#include "dma_macros.h"
#include <math.h>
#define   s_int sizeof(int)
#define   s_double sizeof(double)
//#include "get_high_order_weights_over_areas_struct.h"

#define max(x,y) ( x>y?x:y )
#define min(x,y) ( x<y?x:y )

typedef struct param_areas
{
        double *x_ptr, *dx_ptr;
        double *weights_ptr;
        double *gspts_ptr, *gsweights_ptr;
        int *num_seg_ptr;
        int num_area;
        int num_seg_max;
        int irecons;
        int ngpc;
}param_areas;

void slave_cal_(int *num_area, int *num_seg_max, int *irecons, int *ngpc,double *x,double *dx,double *gspts,double *weights,double *gsweights,int *num_seg,int *iarea);
// __ldm struct param_areas s_ghow;

void  slave_test_c_(struct param_areas *ghow)
{

    struct param_areas s_ghow;
    int num_area,num_seg_max,irecons,ngpc;
    double *x,*dx,*gspts,*weights,*gsweights;
    int *num_seg;
    
    // double *p_x,*p_dx,*p_gspts,*p_weights,*p_gsweights;
    // int *p_num_seg;
    
     int i,myid;
    // dma_syn();
    // if (!_MYID)
    //     printf("&s_ghow = %p\n", &s_ghow);
    
     dma_init();
     pe_get(ghow, &s_ghow, sizeof(param_areas));
    // __builtin_sw_slave_athread_dma(PE_MODE, DMA_GET, ghow, &s_ghow, sizeof(param_areas),
    //         (void*)(&reply_shadow), 0, 0, 0);
    // CRTS_dma_iget((&s_ghow), (ghow), (sizeof(param_areas)), (void*)&reply_shadow);
    // count_shadow ++;
    // CRTS_dma_wait_value(&reply_shadow, count_shadow);
     dma_syn();
    // if (!_MYID)
    //     printf("cpe num_area:%d\n", s_ghow.num_area);
    // return;
    // printf("cpe num_area:%d\n", s_ghow.param[0]);
    // printf("cpe num_seg_max:%d\n", s_ghow.param[1]);
    // printf("cpe irecons:%d\n", s_ghow.param[2]);
    // printf("cpe ngpc:%d\n", s_ghow.param[3]);
    
    num_area = s_ghow.num_area;
    num_seg_max = s_ghow.num_seg_max;
    irecons = s_ghow.irecons;
    ngpc = s_ghow.ngpc;
    
    x = s_ghow.x_ptr;
    dx = s_ghow.dx_ptr;
    gspts = s_ghow.gspts_ptr;
    weights = s_ghow.weights_ptr;
    gsweights = s_ghow.gsweights_ptr;
    
    num_seg = s_ghow.num_seg_ptr;
    
    // p_x       = (double*)ldm_malloc(sizeof(double)*2*num_seg_max*num_area);
    // p_dx       = (double*)ldm_malloc(sizeof(double)*2*num_seg_max*num_area);
    // p_gspts       = (double*)ldm_malloc(sizeof(double)*ngpc);
    // p_weights       = (double*)ldm_malloc(sizeof(double)*irecons*num_area);
    // p_gsweights       = (double*)ldm_malloc(sizeof(double)*ngpc);

    double p_x        [2*num_seg_max*num_area];
    double p_dx       [2*num_seg_max*num_area];
    double p_gspts    [ngpc];
    double p_weights  [irecons*num_area];
    double p_gsweights[ngpc];
    
    // p_num_seg = (int*)ldm_malloc(sizeof(int)*num_area);

     int p_num_seg[num_area];
    
     myid = athread_get_id(-1);
    
    pe_get(x, p_x,sizeof(double)*2*num_seg_max*num_area);
    pe_get(dx, p_dx,sizeof(double)*2*num_seg_max*num_area);
    pe_get(gspts, p_gspts,sizeof(double)*ngpc);
    pe_get(weights, p_weights,sizeof(double)*irecons*num_area);
    pe_get(gsweights, p_gsweights,sizeof(double)*ngpc);
    
    pe_get(num_seg, p_num_seg,sizeof(int)*num_area);
    dma_syn();
    
    for(i = myid+1; i <= num_area; i = i + 64)
    {
        slave_cal_(&num_area, &num_seg_max, &irecons, &ngpc, p_x, p_dx, p_gspts, p_weights, p_gsweights,p_num_seg, &i);
        
        //pe_put(p_weights+i*irecons,weights+i*irecons,sizeof(double)*irecons);
        pe_put(weights+(i-1)*irecons,p_weights+(i-1)*irecons,sizeof(double)*irecons);
        dma_syn();
    }
    
    // ldm_free(p_x, sizeof(double)*2*num_seg_max*num_area);
    // ldm_free(p_dx, sizeof(double)*2*num_seg_max*num_area);
    // ldm_free(p_gspts, sizeof(double)*ngpc);
    // ldm_free(p_weights, sizeof(double)*irecons*num_area);
    // ldm_free(p_gsweights, sizeof(double)*ngpc);
    
    // ldm_free(p_num_seg, sizeof(int)*num_area);

}
