//  CVS: $Id: icesnow.F90,v 1.1.1.1 2004/04/29 06:22:39 lhl Exp $
//     ==================
#include "athread.h"
#include "paramsw.h"
#include "math.h"
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

}Para_icesnow;
Para_icesnow para_icesnow;

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

}Para_tracer;
Para_tracer para_tracer;

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
}Para_tracer5;
Para_tracer5 para_tracer5;

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

}Para_tracer6;
Para_tracer6 para_tracer6;

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
}Para_tracer7;
Para_tracer7 para_tracer7;

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
}Para_invtrit;
Para_invtrit para_invtrit;

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
}Para_upwell;
Para_upwell para_upwell;

int getNextMultiple(int a,int b) {
	    int remainder = a % b;
	    int nextMultiple = a + (b - remainder);
	    return nextMultiple;
}
void computeDistribution_align(int total_line, int num_procs, int *new_myrow, int *new_rowoffsite) {
    int row_per_proc = total_line / num_procs;
    int alignment = 64 / sizeof(double);  // 64字节对齐要求，假设为int类型

    row_per_proc = getNextMultiple(row_per_proc,alignment);

//    printf("row per proc : %d \n",row_per_proc);
    int remaining_rows = total_line - ((num_procs - 1)*row_per_proc);
    if(remaining_rows <= 0){
      printf("Warning ! align distribution wrong! \n");
      printf("Warning ! align distribution wrong! \n");
      printf("Warning ! align distribution wrong! \n");
      exit(0);
    }

    int offset = 0;

    for (int i = 0; i < num_procs; i++) {
	    new_myrow[i] = i * row_per_proc;
	    new_rowoffsite[i] = row_per_proc;
/*
	    if(new_myrow[i] + new_rowoffsite[i] > total_line){
               remaining_rows = total_line - (i*row_per_proc);
	       	    
	    }
	    */
    }
    new_rowoffsite[num_procs - 1] = remaining_rows;

}


void computeDistribution(int total_line, int num_procs, int *new_myrow, int *new_rowoffsite) {
    int i;
    int row_per_proc = total_line / num_procs;
    int remaining_rows = total_line % num_procs;

    for (i = 0; i < num_procs; i++) {
       new_myrow[i] = i * row_per_proc;
       new_rowoffsite[i] = row_per_proc;
       if (i < remaining_rows) {
           new_myrow[i] += i;
           new_rowoffsite[i]++;
       } else {
           new_myrow[i] += remaining_rows;
       }
    }
}


extern void SLAVE_FUN(icesnowsw_());
extern void SLAVE_FUN(tracersw_());
extern void SLAVE_FUN(tracer5_cu_slave_());
extern void SLAVE_FUN(tracer6_cu_slave_());
extern void SLAVE_FUN(tracer7_cu_slave_());
extern void SLAVE_FUN(upwellsw_());
extern void SLAVE_FUN(invtritsw_());

void init_icesnow_para(){

    para_icesnow.mytid = __param_mod_MOD_mytid;
    para_icesnow.d_kmt = __grid_MOD_kmt;
    para_icesnow.d_cp = __pconst_mod_MOD_cp;
    para_icesnow.d_tbice = __pconst_mod_MOD_tbice;
    para_icesnow.d_dzp = __pconst_mod_MOD_dzp;
    para_icesnow.d_at = __tracer_mod_MOD_at;
    para_icesnow.d_licomqice = __tracer_mod_MOD_licomqice;
         
    computeDistribution(imt*jmt,64,new_start,new_offsite);
    computeDistribution_align(imt*jmt*km,64,new_3dstart,new_3doffsite);

    para_icesnow.start=new_start;
    para_icesnow.offsite=new_offsite;
    para_icesnow.start3d=new_3dstart;
    para_icesnow.offsite3d=new_3doffsite;
}

void init_upwell_para(){

    para_upwell.mytid = __param_mod_MOD_mytid;
    para_upwell.kmt = __grid_MOD_kmt;
    para_upwell.ws = __dyn_mod_MOD_ws;
    para_upwell.dzp = __pconst_mod_MOD_dzp;
    para_upwell.vit = __pconst_mod_MOD_vit;
    para_upwell.ohbu = __pconst_mod_MOD_ohbu;
    para_upwell.ohbt = __pconst_mod_MOD_ohbt;
    para_upwell.uk = __work_mod_MOD_uk;
    para_upwell.vk = __work_mod_MOD_vk;
    para_upwell.wka = __work_mod_MOD_wka;
    para_upwell.work = __work_mod_MOD_work;
         
    para_upwell.au0 = __grid_MOD_au0;
    para_upwell.aus = __grid_MOD_aus;
    para_upwell.auw = __grid_MOD_auw;
    para_upwell.ausw = __grid_MOD_ausw;

    para_upwell.hts = __grid_MOD_hts;
    para_upwell.htw = __grid_MOD_htw;
    para_upwell.tarea_r = __grid_MOD_tarea_r;

    computeDistribution(jmt,8,new_myrow,new_rowoffsite);
    computeDistribution(imt,8,new_mycol,new_coloffsite);

    computeDistribution_align(imt*jmt*km,64,new_3dstart,new_3doffsite);

    para_upwell.start_col=new_mycol;
    para_upwell.offsite_col=new_coloffsite;
    para_upwell.start_row=new_myrow;
    para_upwell.offsite_row=new_rowoffsite;
    para_upwell.start3d=new_3dstart;
    para_upwell.offsite3d=new_3doffsite;
}
void init_tracer_para(){

  double aa;
  double c2dtts;
  para_tracer.mytid = __param_mod_MOD_mytid;
  if (strncmp(__pconst_mod_MOD_adv_tracer,"centered",8) == 0) {
    if (__pconst_mod_MOD_ist >= 1) {
      aa = 0.5;
      c2dtts = __pconst_mod_MOD_dts * 2.0;
    } else {
      aa = 0.0;
      c2dtts = __pconst_mod_MOD_dts;
    }
    } else if ((strncmp(__pconst_mod_MOD_adv_tracer,"tspas",5) == 0 || (strncmp(__pconst_mod_MOD_adv_tracer,"prm",3)) == 0))  {
      aa = 0.5;
      c2dtts = __pconst_mod_MOD_dts;
    } else {
      if (para_tracer.mytid == 0) {
      printf("The false advection option for tracer:%s\n",__pconst_mod_MOD_adv_tracer);
      }
   exit(0);
  }
    para_tracer.aa = aa;
    para_tracer.c2dtts = c2dtts;
    para_tracer.onbc = __pconst_mod_MOD_onbc;
    para_tracer.oncc = __pconst_mod_MOD_oncc;

    para_tracer.h0f = __dyn_mod_MOD_h0f;
    para_tracer.h0l = __dyn_mod_MOD_h0l;
    para_tracer.utf = __dyn_mod_MOD_utf;
    para_tracer.vtf = __dyn_mod_MOD_vtf;
    para_tracer.utl = __dyn_mod_MOD_utl;
    para_tracer.vtl = __dyn_mod_MOD_vtl;
    para_tracer.stf = __work_mod_MOD_stf;
    para_tracer.wkd = __work_mod_MOD_wkd;
    para_tracer.wkb = __work_mod_MOD_wkb;


    computeDistribution(imt*jmt,64,new_start,new_offsite);
    computeDistribution_align(imt*jmt*km,64,new_3dstart,new_3doffsite);

    para_tracer.start=new_start;
    para_tracer.offsite=new_offsite;
    para_tracer.start3d=new_3dstart;
    para_tracer.offsite3d=new_3doffsite;
}
void init_tracer5_para(){

    para_tracer5.mytid = __param_mod_MOD_mytid;
    para_tracer5.od0 = __pconst_mod_MOD_od0;
    para_tracer5.boundary_restore = __pconst_mod_MOD_boundary_restore;
    para_tracer5.gamma = __pconst_mod_MOD_gamma;

    para_tracer5.kmt = __grid_MOD_kmt;
    para_tracer5.swv = __forc_mod_MOD_swv;
#ifdef SOLAR
    para_tracer5.pen = __pmix_mod_MOD_pen;
#endif
#ifdef SOLARCHLORO
    para_tracer5.pen_chl = __pmix_mod_MOD_pen_chl;
#endif
    para_tracer5.vit = __pconst_mod_MOD_vit;
    para_tracer5.tf = __work_mod_MOD_tf;
    para_tracer5.penetrate = __tracer_mod_MOD_penetrate;
    para_tracer5.odzp = __pconst_mod_MOD_odzp;
    para_tracer5.odzt = __pconst_mod_MOD_odzt;
    para_tracer5.atb  = __tracer_mod_MOD_atb;
    para_tracer5.dz  = __tracer_mod_MOD_dz;
    para_tracer5.ssf  = __forc_mod_MOD_ssf;
    para_tracer5.sss  = __forc_mod_MOD_sss;
    para_tracer5.ifrac  = __buf_mod_MOD_ifrac;
    para_tracer5.net  = __tracer_mod_MOD_net;
    para_tracer5.tarea  = __grid_MOD_tarea;
    para_tracer5.dt_restore  = __tracer_mod_MOD_dt_restore;
    para_tracer5.fresh  = __forc_mod_MOD_fresh;
    para_tracer5.seaice  = __forc_mod_MOD_seaice;
    para_tracer5.stf = __work_mod_MOD_stf;
    para_tracer5.wkc = __work_mod_MOD_wkc;
#ifdef SSS_restore
    para_tracer5.sss_restore = __tracer_mod_MOD_sss_restore;
#endif
    computeDistribution(jmt,8,new_myrow,new_rowoffsite);
    computeDistribution(imt,8,new_mycol,new_coloffsite);

    computeDistribution_align(imt*jmt*km,64,new_3dstart,new_3doffsite);

    para_tracer5.start_col=new_mycol;
    para_tracer5.offsite_col=new_coloffsite;
    para_tracer5.start_row=new_myrow;
    para_tracer5.offsite_row=new_rowoffsite;
    para_tracer5.start3d=new_3dstart;
    para_tracer5.offsite3d=new_3doffsite;
}
void init_tracer6_para(){

    para_tracer6.mytid = __param_mod_MOD_mytid;
    para_tracer6.nx = __buf_mod_MOD_nx;
    para_tracer6.od0cp = __pconst_mod_MOD_od0cp;
    para_tracer6.gamma = __pconst_mod_MOD_gamma;

    para_tracer6.kmt = __grid_MOD_kmt;
    para_tracer6.swv = __forc_mod_MOD_swv;
    para_tracer6.vit = __pconst_mod_MOD_vit;
    para_tracer6.tf = __work_mod_MOD_tf;
    para_tracer6.odzp = __pconst_mod_MOD_odzp;
    para_tracer6.atb  = __tracer_mod_MOD_atb;
    para_tracer6.net  = __tracer_mod_MOD_net;
    para_tracer6.seaice  = __forc_mod_MOD_seaice;
    para_tracer6.stf = __work_mod_MOD_stf;
    para_tracer6.sst = __forc_mod_MOD_sst;
    para_tracer6.nswv = __forc_mod_MOD_nswv;
    para_tracer6.tsf = __forc_mod_MOD_tsf;
    para_tracer6.dqdt = __forc_mod_MOD_dqdt;

    computeDistribution(jmt,8,new_myrow,new_rowoffsite);
    computeDistribution(imt,8,new_mycol,new_coloffsite);

    computeDistribution_align(imt*jmt*km,64,new_3dstart,new_3doffsite);

    para_tracer6.start_col=new_mycol;
    para_tracer6.offsite_col=new_coloffsite;
    para_tracer6.start_row=new_myrow;
    para_tracer6.offsite_row=new_rowoffsite;
    para_tracer6.start3d=new_3dstart;
    para_tracer6.offsite3d=new_3doffsite;
}
void init_tracer7_para(){
  double c2dtts;
  double d_iau_gamma;
  if (strncmp(__pconst_mod_MOD_adv_tracer,"centered",8) == 0) {
    para_tracer7.str_adv_tracer_flag = 0;
    if (__pconst_mod_MOD_ist >= 1) {
      c2dtts = __pconst_mod_MOD_dts * 2.0;
    } else {
      c2dtts = __pconst_mod_MOD_dts;
    }
    } else if ((strncmp(__pconst_mod_MOD_adv_tracer,"tspas",5) == 0)|| (strncmp(__pconst_mod_MOD_adv_tracer,"prm",3)) == 0)  {
      para_tracer7.str_adv_tracer_flag = 1;
      c2dtts = __pconst_mod_MOD_dts;
    } else {
      if (para_tracer.mytid == 0) {
      printf("The false advection option for tracer:%s\n",__pconst_mod_MOD_adv_tracer);
      }
   exit(0);
  }
    para_tracer7.c2dtts = c2dtts;
    para_tracer7.mytid = __param_mod_MOD_mytid;
    para_tracer7.boundary_restore = __pconst_mod_MOD_boundary_restore;
    para_tracer7.simple_assm = __pconst_mod_MOD_simple_assm;
    para_tracer7.kvt = __pconst_mod_MOD_kvt;
    para_tracer7.dts = __pconst_mod_MOD_dts;
    para_tracer7.iau = __pconst_mod_MOD_enoi_iau;

    d_iau_gamma = __pconst_mod_MOD_iau_gamma;
    para_tracer7.d_lamda = 1.0/(d_iau_gamma*86400.0);

    para_tracer7.vit = __pconst_mod_MOD_vit;
    para_tracer7.tf = __work_mod_MOD_tf;
    para_tracer7.atb  = __tracer_mod_MOD_atb;
    para_tracer7.dt_restore  = __tracer_mod_MOD_dt_restore;
    para_tracer7.restore  = __forc_mod_MOD_restore;
    para_tracer7.restore_at  = __tracer_mod_MOD_restore_at;
   // para_tracer7.ori  = __tracer_mod_MOD_dt_restore;
    para_tracer7.vtl  = __dyn_mod_MOD_vtl;
    para_tracer7.at  = __tracer_mod_MOD_at;
    para_tracer7.d_at_iau  = __tracer_mod_MOD_at_iau;

    
    computeDistribution(jmt,8,new_myrow,new_rowoffsite);
    computeDistribution(imt,8,new_mycol,new_coloffsite);

    computeDistribution_align(imt*jmt*km,64,new_3dstart,new_3doffsite);

    para_tracer7.start_col=new_mycol;
    para_tracer7.offsite_col=new_coloffsite;
    para_tracer7.start_row=new_myrow;
    para_tracer7.offsite_row=new_rowoffsite;
    para_tracer7.start3d=new_3dstart;
    para_tracer7.offsite3d=new_3doffsite;
}
void init_invtrit_para(){

    para_invtrit.mytid = __param_mod_MOD_mytid;

    para_invtrit.kmt = __grid_MOD_kmt;
    para_invtrit.vit = __pconst_mod_MOD_vit;
    para_invtrit.odzp = __pconst_mod_MOD_odzp;
    para_invtrit.odzt = __pconst_mod_MOD_odzt;
    computeDistribution(jmt,8,new_myrow,new_rowoffsite);
    computeDistribution(imt,8,new_mycol,new_coloffsite);

    computeDistribution_align(imt*jmt*km,64,new_3dstart,new_3doffsite);

    para_invtrit.start_col=new_mycol;
    para_invtrit.offsite_col=new_coloffsite;
    para_invtrit.start_row=new_myrow;
    para_invtrit.offsite_row=new_rowoffsite;
    para_invtrit.start3d=new_3dstart;
    para_invtrit.offsite3d=new_3doffsite;
}
void init_para_(){

//    init_icesnow_para();
    init_tracer_para();
    init_tracer5_para();
    init_tracer6_para();
    init_tracer7_para();
//    init_upwell_para();
//    init_invtrit_para();
}
void init_para_withp(Para_icesnow parahost){

     
    parahost.mytid = __param_mod_MOD_mytid;
	 
}

void icesnow_sw()
{
    athread_spawn(icesnowsw_,&para_icesnow);
    athread_join(); 
}

void tracer_sw()
{
    athread_spawn(tracersw_,&para_tracer);
    athread_join(); 
}
void tracer5_cu_slave(int n)
{
    para_tracer5.tracer_n = n-1;
    athread_spawn(tracer5_cu_slave_,&para_tracer5);
    athread_join(); 
}

void tracer6_cu_slave(int n,double e_norm2)
{
    para_tracer6.tracer_n = n-1;

    para_tracer6.err_norm2 = e_norm2;
    athread_spawn(tracer6_cu_slave_,&para_tracer6);
    athread_join(); 
}
void tracer7_cu_slave(int n)
{
    para_tracer7.tracer_n = n-1;
    athread_spawn(tracer7_cu_slave_,&para_tracer7);
    athread_join(); 
}
void invtrit_sw(double h_c2dtts)
{
    para_invtrit.wk = __dyn_mod_MOD_vtl;
    para_invtrit.topbc = __work_mod_MOD_stf;
    para_invtrit.dcb = __work_mod_MOD_wkc;
    para_invtrit.c2dtts = h_c2dtts;
    athread_spawn(invtritsw_,&para_invtrit);
    athread_join(); 
}
void upwell_sw()
{
    para_upwell.uwk = __work_mod_MOD_wkd; 
    para_upwell.vwk = __work_mod_MOD_wkb;  
    para_upwell.h0wk = __work_mod_MOD_stf;
    athread_spawn(upwellsw_,&para_upwell);
    athread_join(); 
}
/*
void icesnow_cpp()


  const double sal_ice = 4.0;
  const double sal_ocn = 34.7;
  const double heat_ice_fusion = SHR_CONST_LATICE_;
  int iblock = 0 ;
  double tdiff = 0.0;
    if(__param_mod_MOD_mytid == 0 )
 //icesnow 1
    for (int k = 0 ; k < 1; ++k) {   // useless, just for consistency
      for (int j = 0 ; j < jmt; ++j) {
        for (int i = 0; i < imt; ++i) {
//    if(__param_mod_MOD_mytid == 0 ){
 //         printf("at(%d,%d): %f \n",i,j,para.at[j*imt + i]);
//          printf("at(%d,%d): %f \n",i,j,__tracer_mod_MOD_at[iblock][j][i]);
//}
          if(para.at[j*imt + i] < para.tbice && para.kmt[j*imt + i] > 0) {
            tdiff = para.tbice - para.at[k*jmt*imt + j*imt + i];
	    para.licomqice[j*imt + i] += tdiff * para.dzp[k] / para.dzp[0];
            para.at[km*imt*jmt + k * imt*jmt + j*imt + i]  += tdiff * (sal_ocn - sal_ice)
                  * para.cp / heat_ice_fusion * 0.001;
            para.at[k*jmt*imt + j*imt + i] = para.tbice;
            }//end if
         }
        }
      }

    for (int j = 0 ; j < jmt; ++j) {
      for (int i = 0; i < imt; ++i) {
        if(para.licomqice[j*imt + i] > 0.0 &&
           para.at[j*imt + i] > para.tbice) {

          tdiff = fmin((para.at[j*imt + i] - para.tbice),
              para.licomqice[j*imt + i]);

          para.licomqice[j*imt + i] -= tdiff;
          para.at[j*imt + i]  -= tdiff;
          para.at[km*imt*jmt + j*imt + i]  -= tdiff * (sal_ocn - sal_ice)
               * para.cp / heat_ice_fusion * 0.001;
        }//end if
      }
    }


    for (int k = 0 ; k < km; ++k) {
       for (int j = 0 ; j < jmt; ++j) {
          for (int i = 0; i < imt; ++i) {
             if (para.at[k*imt*jmt + j*imt + i] < para.tbice) {
                para.at[k*imt*jmt + j*imt + i] = para.tbice;
              }
          }
       }
    }//end for

       return ;
}
*/

