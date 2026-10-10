#ifndef PARA_H
#define PARA_H

#include "def-undef.h"
#include <stdbool.h>
//#include "precision_mod.h"
#include "param_mod.h"

#define CHAR_LEN 256

#define SHR_CONST_LATICE_  (double)3.337e5     // latent heat of fusion      ~ J/kg
#define licom_BlockSizeX   BLCKX // size of block in first  horizontal dimension
#define licom_BlockSizeY   BLCKY // size of block in second horizontal dimension
#define nghost  2       // number of ghost cells around each block
#define nx_block (licom_BlockSizeX + 2*nghost)   //  x,y dir including ghost
#define ny_block (licom_BlockSizeY + 2*nghost)   //  cells

#define jmt ny_block

#define imt nx_block  

#define km NKM     // Number of Grid Points in Vertical Direction



#define ICESNOW_TIME
#undef TRACER_TIME
//extern struct Para para;

int new_myrow[8],new_mycol[8],new_rowoffsite[8],new_coloffsite[8];
int new_start[64],new_offsite[64];
int new_3dstart[64],new_3doffsite[64];


extern double __pconst_mod_MOD_cp;
extern double __pconst_mod_MOD_od0;
extern double __pconst_mod_MOD_akt[1][2][km][jmt][imt];
extern double __pconst_mod_MOD_ahv;
extern double __pconst_mod_MOD_dts;
extern double __pconst_mod_MOD_ist;
extern double __pconst_mod_MOD_vit[1][km][jmt][imt];
extern double __pconst_mod_MOD_od0;
extern double __pconst_mod_MOD_aft1;
extern double __pconst_mod_MOD_aft2;
extern double __pconst_mod_MOD_od0cp;
extern double __pconst_mod_MOD_odzp[km];
extern double __pconst_mod_MOD_odzt[km];
extern double __pconst_mod_MOD_onbc;
extern double __pconst_mod_MOD_oncc;
extern double __pconst_mod_MOD_ohbu[1][jmt][imt];
extern double __pconst_mod_MOD_ohbt[1][jmt][imt];
extern double __pconst_mod_MOD_gamma;
extern double __pconst_mod_MOD_dwndmix;
extern char   __pconst_mod_MOD_adv_tracer[80];
extern bool   __pconst_mod_MOD_simple_assm;
extern int    __pconst_mod_MOD_boundary_restore;
extern int    __pconst_mod_MOD_kvt;
extern double __pconst_mod_MOD_tbice;
extern double __pconst_mod_MOD_dzp[km];
extern bool   __pconst_mod_MOD_enoi_iau;
extern double __pconst_mod_MOD_iau_gamma;

extern double __tracer_mod_MOD_tend[1][2][jmt][imt];
extern double __tracer_mod_MOD_dt_diff[1][2][jmt][imt];
extern double __tracer_mod_MOD_restore_at[1][2][jmt][imt];
extern double __tracer_mod_MOD_net[1][2][jmt][imt];
extern double __tracer_mod_MOD_atb[1][2][km+1][jmt][imt];
extern double __tracer_mod_MOD_at[1][2][km][jmt][imt];
extern double __tracer_mod_MOD_dz[1][2][km][jmt][imt];
extern double __tracer_mod_MOD_dt_restore[1][2][km][jmt][imt];
extern double __tracer_mod_MOD_fw_norm2;
extern double __tracer_mod_MOD_penetrate[1][jmt][imt];
extern double __tracer_mod_MOD_licomqice[1][jmt][imt];
extern double __tracer_mod_MOD_sss_restore[1][jmt][imt];
extern double __tracer_mod_MOD_at_iau[1][2][km][jmt][imt];

extern double __work_mod_MOD_tf[1][km][jmt][imt];
extern double __work_mod_MOD_stf[1][jmt][imt];
extern double __work_mod_MOD_wkb[1][km][jmt][imt];
extern double __work_mod_MOD_wkc[1][km][jmt][imt];
extern double __work_mod_MOD_wkd[1][km][jmt][imt];
extern double __work_mod_MOD_wka[1][km][jmt][imt];
extern double __work_mod_MOD_work[1][jmt][imt];
extern double __work_mod_MOD_uk[1][km][jmt][imt];
extern double __work_mod_MOD_vk[1][km][jmt][imt];

extern double __dyn_mod_MOD_ws[1][km+1][jmt][imt];
extern double __dyn_mod_MOD_h0f[1][jmt][imt];
extern double __dyn_mod_MOD_h0l[1][jmt][imt];
extern double __dyn_mod_MOD_utf[1][km][jmt][imt];
extern double __dyn_mod_MOD_vtf[1][km][jmt][imt];
extern double __dyn_mod_MOD_utl[1][km][jmt][imt];
extern double __dyn_mod_MOD_vtl[1][km][jmt][imt];


extern double __forc_mod_MOD_ssf[1][jmt][imt];
extern double __forc_mod_MOD_sss[1][jmt][imt];
extern double __forc_mod_MOD_sst[1][jmt][imt];
extern double __forc_mod_MOD_swv[1][jmt][imt];
extern double __forc_mod_MOD_dqdt[1][jmt][imt];
extern double __forc_mod_MOD_nswv[1][jmt][imt];
extern double __forc_mod_MOD_tsf[1][jmt][imt];
extern double __forc_mod_MOD_fresh[1][jmt][imt];
extern double __forc_mod_MOD_seaice[1][jmt][imt];
extern double __forc_mod_MOD_restore[1][2][km][jmt][imt];


extern int __grid_MOD_kmt[1][ny_block][nx_block];
extern double __grid_MOD_tarea_r[1][ny_block][nx_block];
extern double __grid_MOD_tarea[1][ny_block][nx_block];
extern double __grid_MOD_hts[1][ny_block][nx_block];
extern double __grid_MOD_htw[1][ny_block][nx_block];
extern double __grid_MOD_au0[1][ny_block][nx_block];
extern double __grid_MOD_aus[1][ny_block][nx_block];
extern double __grid_MOD_auw[1][ny_block][nx_block];
extern double __grid_MOD_ausw[1][ny_block][nx_block];
extern double __grid_MOD_area_t;
extern char __grid_MOD_horiz_grid_opt[CHAR_LEN];

#ifdef SOLAR
extern double __pmix_mod_MOD_pen[km - 1];
#endif 
#ifdef SOLARCHLORO
extern double __pmix_mod_MOD_pen_chl[1][km][jmt][imt];
#endif
extern double __buf_mod_MOD_ifrac[1][jmt][imt];
extern int __buf_mod_MOD_nx;

//extern double tracer_mod_mp_at_[max_blocks_clinic][ntra][km][jmt][imt];
//extern double tracer_mod_mp_licomqice_[max_blocks_clinic][jmt][imt];
extern int __domain_MOD_nblocks_clinic;

#endif /* PARA_H */
