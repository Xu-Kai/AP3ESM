#include <crts.h>
#include <slave.h>
#include <simd.h>
#include <math.h>
#include <stdio.h>
#include <float.h>
#include <string.h>
#include <stdbool.h>
#include "swperf.h" 
#include "param_mod.h"
#include "canuto_mod.h" 
#include "pconst_mod.h"
#include "constant_mod.h"
#include "def-undef.h"

typedef struct turb_kernelPara{
	int cuGrid[3];
	int cuBlock[3];
    int d_km;
	int d_imt;
	int d_jmt;
	int rank;
    int *d_irimax;
    int *d_kmt;
	double *d_akmt;
	double *d_akt;
	double *d_ak_tide;
	double *d_amld;
	double *d_amtaun2a1;
	double *d_and2on2a1;
	double *d_at;
	double *d_back_ra_r;
	double *d_buoysol;
	double *d_buoytur;
	double *d_dzp;
	double *d_fcort;
	double *d_fztidal;
	double *d_fz_tide;
	double *d_pdensity;
	double *d_rib;
    double *d_ridb;
	double *d_ricdttms;
	double *d_richardson;
	double *d_rict;
	double *d_rit;
	double *d_s2t;
	double *d_sha1;
	double *d_shb;
	double *d_sh_r1;
	double *d_slq2b;
	double *d_slq2_r1;
	double *d_sma1;
	double *d_smb;
	double *d_sm_r1;
	double *d_ssa1;
	double *d_ssb;
	double *d_ss_r1;
	double *d_turb_param;
	double *d_vit;
	double *d_wave_dis;
	double *d_wp3_tidal;
	double *d_zkp;
#ifdef BCKMEX
	double *d_diff_back;
#endif
	//double *d_wk1;
	//double *d_wk2;
	//double *d_wk3;
}turb_Para;

//dynamic atomic definition
#define CP_CORE      64
#define RCV_CORE     0
#define DATAPAKAGE   1

//ATHREAD global definition 
__thread_local int warp_size = 32; 
__thread_local int blockIndex;


__inline__  double  wavelat(double xf, double yn){	
     // return xf *acosh(yn /xf);
     double X = yn/xf;
    return xf * (log(X + sqrt(X*X-1.0e0)));

}
 double  eplatidepend_fun(double f, double an){	
      double an0,f_30;	
	
      	an0=5.24e-3;	
      double anum,den;	
	
      double omega_tmp;

      	omega_tmp=4.0e0 *atan(1.0e0) /43082.0e0;
      	f_30=omega_tmp;	
      	den=wavelat(f_30,an0);	
      	anum=wavelat(f,an);	
      return anum /den;	

}
__inline__  int  sign_int(int x, int y){	
      return (y) >= 0 ? abs(x) : -abs(x);	

}
__inline__  double  sign(double x, double y){	
      return (y) >= 0.0e0 ? fabs(x) : -fabs(x);	

}
 void  formld(int n, double *z, double *t, double *amld){	
      int k;	
	
      double tm;	
	
      for(	k=0;k<n; k++)	
	{	
      if(fabs(t[k]-t[0])> 0.1)	
      {	
      	tm=t[0]-sign(0.1e0,t[0]-t[k]);	
      	*amld=z[k]+(z[k-1]-z[k]) *(tm-t[k]) /(t[k-1]-t[k]+1.0e-20);	
      return ;	
      }	
	}	
      	*amld=z[n-1];	

}
 void  interp1d_expabs(double x, double *x_1, double *slq2_1, double *sm_1, double *sh_1, 
		 double *ss_1, double *slq2, double *sm, double *sh, double *ss, double delta, double rat){	
      double deltax,deltaxta,dsh_x,dslq2_x,dsm_x,dss_x,tabindx;	
	
      int lx0,lx1;	
	
      if(x> x_1[2 *mt])	
      	x=x_1[2 *mt];	
      else if(x<x_1[0])	
      	x=x_1[0];	

      if(fabs(x) < x_1[mt + mt0]){
          lx1=(int)(x /delta)+(int)sign((double)1,x);
      }
      else if(fabs(x)>=x_1[2 *mt]){
          lx0=(int)sign((double)mt,x);
          lx1=lx0;
      }
      else{
          tabindx=sign((double)(mt0)+((log(fabs(x))-log(x_1[mt+mt0])) /log(rat)),x);
          lx1=(int)(tabindx)+(int)sign((double)1,x);
      }

      if(fabs(x_1[mt+lx1])<(fabs(x)))
          lx1+=sign_int(1,lx1);
      else if(fabs(x_1[mt+lx1-sign_int(1,lx1)])> fabs(x))
          lx1+=-sign_int(1,lx1);
     lx0=lx1-(int)sign((double)1,x);
     if(fabs(x)<=0.0e0)
         lx1=1;

/*     if(fabs(x)>=x_1[mt+mt0]&&(fabs(x)>=x_1[2 *mt]))
     {
      	lx0=(int)sign((double)mt,x);	
      	lx1=lx0;	
     }
      else{
	    if(fabs(x)<x_1[mt+mt0])
      	    lx1=(int)(x /delta)+(int)sign((double)1,x);	
        
	else{
	    tabindx=sign((double)(mt0)+((log(fabs(x))-log(x_1[mt+mt0])) /log(rat)),x);	
      	    lx1=(int)(tabindx)+(int)sign((double)1,x);	
        }	
      
	if(fabs(x_1[mt+lx1])<(fabs(x)))	
            lx1+=sign_int(1,lx1);	
      
	else if(fabs(x_1[mt+lx1-sign_int(1,lx1)])> fabs(x))	
      	    lx1+=-sign_int(1,lx1);	
	
      	lx0=lx1-(int)sign((double)1,x);	
      
	if(fabs(x)<=0.0e0)	
      	    lx1=1;	
      }	*/
      
      if((x> 0.0e0&&(x<x_1[mt+lx0]||x> x_1[mt+lx1]))||(x<0.0e0&&(x> x_1[mt+lx0]||x<x_1[mt+lx1])))	
          return ;	
      	
      deltaxta=x_1[mt+lx1]-x_1[mt+lx0];	
      deltax=x-x_1[mt+lx0];	
      
      if(lx1==lx0)	
      {	
      	dslq2_x=0.0e0;	
      	dsm_x=0.0e0;	
      	dsh_x=0.0e0;	
      	dss_x=0.0e0;	
      }	
      else{
	    dslq2_x=(slq2_1[mt+lx1]-slq2_1[mt+lx0]) /deltaxta;
      	dsm_x=(sm_1[mt+lx1]-sm_1[mt+lx0]) /deltaxta;	
      	dsh_x=(sh_1[mt+lx1]-sh_1[mt+lx0]) /deltaxta;	
      	dss_x=(ss_1[mt+lx1]-ss_1[mt+lx0]) /deltaxta;	
      }	
      	*slq2=slq2_1[mt+lx0]+dslq2_x *deltax;	
      	*sm=sm_1[mt+lx0]+dsm_x *deltax;	
      	*sh=sh_1[mt+lx0]+dsh_x *deltax;	
      	*ss=ss_1[mt+lx0]+dss_x *deltax;	
      return ;	

}
 void  interp2d_expabs(double *ri, double *rid, double *slq2, double *sm, double *sh, double *ss, 
		 double *d_rib, double *d_ridb, double *d_shb, double *d_slq2b, double *d_smb, double *d_ssb,
         int *d_irimax, double d_dri, double d_rri){
      double deltari,deltarid,deltaridta,deltarita,dsh_ri,dsh_rid,dslq2_rid,dsm_ri,dsm_rid,dss_ri,dss_rid,tabindrid,dslq2_ri,tabindri;	
	
      int lri0,lri1,lrid0,lrid1;	
	
      if(*ri> d_rib[2 *mt])
      {	
        if(fabs(*rid)<=*ri)	
        {	
      	  *rid=d_rib[2 *mt] *(*rid / *ri);
      	  *ri=d_rib[2 *mt];	
        }	
        else if(*rid> *ri)	
        {	
      	  *ri=d_ridb[2 *mt] *(*ri / *rid);
      	  *rid=d_ridb[2 *mt];
        }	
        else if(*rid<-*ri)	
        {	
      	  *ri=d_ridb[0] *(*ri / *rid);
      	  *rid=d_ridb[0];
        }

      }	
      else if(*ri<d_rib[0])
      {	
        if(fabs(*rid)<=-*ri)	
        {	
      	  *rid=d_rib[0] *(*rid / *ri);	
      	  *ri=d_rib[0];
        }	
        else if(*rid> -*ri)	
        {	
      	  *ri=d_ridb[2 *mt] *(*ri / *rid);
      	  *rid=d_ridb[2 *mt];
        }	
        else if(*rid<*ri)	
        {	
      	  *ri=d_ridb[0] *(*ri / *rid);
      	  *rid=d_ridb[0];
        }

      }	
      else if(*rid> d_ridb[2 *mt])
      {	
      	*ri=d_ridb[2 *mt] *(*ri / *rid);
      	*rid=d_ridb[2 *mt];
      }
      else if(*rid<d_ridb[0])
      {	
      	*ri=d_ridb[0] *(*ri / *rid);
      	*rid=d_ridb[0];
      }
	
      double ri_tmp=*ri;	
	
      double rid_tmp=*rid;	
	
      double tmp=d_ridb[mt+mt0];
	
      double tmp1=d_ridb[2 *mt];

     if(fabs(rid_tmp)<d_ridb[mt+mt0])
         lrid1=(int)(rid_tmp /d_dri)+(int)sign((double)1,rid_tmp);
     else if(fabs(rid_tmp)>=d_ridb[2 *mt]){
         lrid0=(int)sign((double)(mt),rid_tmp);
         lrid1=lrid0;
     }
     else{
         tabindrid=sign((double)(mt0)+((log(fabs(rid_tmp))-log(d_ridb[mt+mt0])) /log(d_rri)),rid_tmp);
         lrid1=(int)(tabindrid)+(int)sign((double)1,rid_tmp);
     }

     if((fabs(d_ridb[mt+lrid1]))<(fabs(rid_tmp)))
         lrid1+=sign_int(1,lrid1);
     else if(fabs(d_ridb[mt+lrid1-sign_int(1,lrid1)])> fabs(rid_tmp))
         lrid1+=-sign_int(1,lrid1);

     lrid0=lrid1-(int)sign((double)1,rid_tmp);

     if(rid_tmp==0.0e0)
         lrid1=1;
      
      if((rid_tmp> 0.0e0&&(rid_tmp<d_ridb[mt+lrid0]||rid_tmp> d_ridb[mt+lrid1]))||
		      (rid_tmp<0.0e0&&(rid_tmp> d_ridb[mt+lrid0]||rid_tmp<d_ridb[mt+lrid1])))
        return ;	
      
      if(ri_tmp> fmin(d_rib[mt+d_irimax[mt+lrid0]],d_rib[mt+d_irimax[mt+lrid1]]))
      {	
      	*slq2=0.0e0;	
      	*sm=0.0e0;	
      	*sh=0.0e0;	
      	*ss=0.0e0;	
        return ;	
      }

      if(fabs(ri_tmp)<d_rib[mt+mt0])
          lri1=(int)(ri_tmp /d_dri)+(int)sign((double)1,ri_tmp);
      else if(fabs(ri_tmp)>=d_rib[mt+mt]){
          lri0=(int)sign((double)(mt),ri_tmp);
          lri1=lri0;
      }
      else{
          tabindri=sign((double)(mt0)+((log(fabs(ri_tmp))-log(d_rib[mt+mt0])) /log(d_rri)),ri_tmp);
          lri1=(int)(tabindri)+(int)sign((double)1,ri_tmp);
      }

     if((fabs(d_rib[mt+lri1]))<(fabs(ri_tmp)))
         lri1+=sign_int(1,lri1);
     else if(fabs(d_rib[mt+lri1-sign_int(1,lri1)])> fabs(ri_tmp))
         lri1-=sign_int(1,lri1);

     lri0=lri1-(int)sign((double)1,ri_tmp);

     if(ri_tmp==0.0e0)
         lri1=1;

      if((ri_tmp> 0.0e0&&(ri_tmp<d_rib[mt+lri0]||ri_tmp> d_rib[mt+lri1]))||
		      (ri_tmp<0.0e0&&(ri_tmp> d_rib[mt+lri0]||ri_tmp<d_rib[mt+lri1])))	
        return ;	
      	deltaridta=d_ridb[mt+lrid1]-d_ridb[mt+lrid0];
      	deltarita=d_rib[mt+lri1]-d_rib[mt+lri0];	
      	deltarid=rid_tmp-d_ridb[mt+lrid0];
      	deltari=ri_tmp-d_rib[mt+lri0];

      if(lrid1==lrid0){
      	dslq2_rid=0.0e0;
      	dsm_rid=0.0e0;
      	dsh_rid=0.e0;
      	dss_rid=0.e0;
      }	
      else{
        dslq2_rid=(d_slq2b[(2 *mt+1) *(lrid1+mt)+lri0+mt]-d_slq2b[(2 *mt+1) *(lrid0+mt)+lri0+mt]) /deltaridta;
      	dsm_rid=(d_smb[(2 *mt+1) *(lrid1+mt)+lri0+mt]-d_smb[(2 *mt+1) *(lrid0+mt)+lri0+mt]) /deltaridta;
      	dsh_rid=(d_shb[(2 *mt+1) *(lrid1+mt)+lri0+mt]-d_shb[(2 *mt+1) *(lrid0+mt)+lri0+mt]) /deltaridta;
      	dss_rid=(d_ssb[(2 *mt+1) *(lrid1+mt)+lri0+mt]-d_ssb[(2 *mt+1) *(lrid0+mt)+lri0+mt]) /deltaridta;
      }
      if(lri1 == lri0){
          dslq2_ri=0.0e0;
          dsm_ri=0.0e0;
          dsh_ri=0.e0;
          dss_ri=0.e0;
      }
      else{
          dslq2_ri=(d_slq2b[(2 *mt+1) *(lrid0+mt)+lri1+mt]-d_slq2b[(2 *mt+1) *(lrid0+mt)+lri0+mt]) /deltarita;
          dsm_ri=(d_smb[(2 *mt+1) *(lrid0+mt)+lri1+mt]-d_smb[(2 *mt+1) *(lrid0+mt)+lri0+mt]) /deltarita;
          dsh_ri=(d_shb[(2 *mt+1) *(lrid0+mt)+lri1+mt]-d_shb[(2 *mt+1) *(lrid0+mt)+lri0+mt]) /deltarita;
          dss_ri=(d_ssb[(2 *mt+1) *(lrid0+mt)+lri1+mt]-d_ssb[(2 *mt+1) *(lrid0+mt)+lri0+mt]) /deltarita;
      }
      	*slq2=d_slq2b[(2 *mt+1) *(lrid0+mt)+lri0+mt]+dslq2_ri *deltari+dslq2_rid *deltarid;
      	*sm=d_smb[(2 *mt+1) *(lrid0+mt)+lri0+mt]+dsm_ri *deltari+dsm_rid *deltarid;
      	*sh=d_shb[(2 *mt+1) *(lrid0+mt)+lri0+mt]+dsh_ri *deltari+dsh_rid *deltarid;
      	*ss=d_ssb[(2 *mt+1) *(lrid0+mt)+lri0+mt]+dss_ri *deltari+dss_rid *deltarid;
      return ;	

}
 void  turb_device(double *z, double *t, double *rh, double *ri, double *rid, double *s2, double fricmx, 
		 double wndmix, double *an2, double buoytur, double buoysol, double coriol, double *amld, 
		 double *akm, double *akh, double *aks, int n, int na, int nmax, double *d_amtaun2a1, 
		 double *d_and2on2a1, double *d_back_ra_r, double *d_rib, double *d_ridb, double *d_sha1, double *d_shb,
		 double *d_sh_r1, double *d_slq2b, double *d_slq2_r1, double *d_sma1, double *d_smb, 
		 double *d_sm_r1, double *d_ssa1, double *d_ssb, double *d_ss_r1, double *d_turb_param, int *d_irimax,int ii, int jj){	
      int k;	
	
      double v_back=0.0e0,t_back=0.0e0,s_back=0.0e0;	
	
      double visc_cbu_limit1=0.0e0,diff_cbt_limit1=0.0e0;	
	
      bool lifupper,lifepsy;
	
      double al2=0.0e0,tmp=0.0e0,slq2=0.0e0,sm=0.0e0,ss=0.0e0,sh=0.0e0;	
	
      double aldeep=0.0e0;	
	
      int ilmldpoint=0,ilmldpointneg=0,ipoint=0,icall=0;	
	
      double amtaun2=0.0e0;	
	
      double deltheta_r_1,d_rri=d_turb_param[0],d_rnd2on2=d_turb_param[1],d_dri=d_turb_param[2],d_deltheta_r=d_turb_param[3],d_b1=d_turb_param[4],d_theta_rcrp=d_turb_param[5],d_theta_rcrn=d_turb_param[6],d_dand2on2=d_turb_param[10];	
	
      double anlq2 = 0.0e0,back_ra_r1 = 0.0e0,back_ri1=0.0e0,back_ric1=0.0e0,back_rid1=0.0e0,back_rit1=0.0e0;
	
      double dback_ra_r_o_dtheta,delback_ra_r,delsh_back,delslq2_back,delsm_back,delss_back,dsh_back_o_dtheta,dslq2_back_o_dtheta,dsm_back_o_dtheta,dss_back_o_dtheta,eplatidepend,epson2,epson2_;	
	
      double epsy = 0.0e0;
	
      int ifpureshear=0,ifrafglt=0,inegproblem=0,iproblem=0,itheta_r0=0,itheta_r1=0,jtheta_r0=0;
	
      double ra_r = 0.0e0,ric = 0.0e0,rit = 0.0e0,sh_back=0.0e0;
      double slq2_back=0.0e0,sm_back=0.0e0,ss_back=0.0e0,theta_r = 0.0e0;
      double theta_r0 = 0.0e0,theta_r1 = 0.0e0,tmp_back = 0.0e0;
	
      int nb=n > 0 ? n : 0;	
	
      visc_cbu_limit1=fricmx;	
      diff_cbt_limit1=fricmx;	
      double buoytot=buoytur+buoysol;	
	
      if(n> 0)	
      {	
      	formld(n,z,t,amld);	
      }	
      else if(n==0)
      	*amld=z[0];	

      else
	    *amld=0.0e0;

	
      double al0=0.17 **amld;	
	
      int ifbelow=0;
	
      double ifnofsmall,an;	
	
      double ri1 = 0.0e0,rid1 = 0.0e0;
	
      for(k=0;k<n; k++)	
      {	
      	ifnofsmall=0;	
        
	    if(an2[k]>=0.0e0) an=sqrt(an2[k]);
	
        if(an<fabs(coriol)) ifnofsmall=1;

      	ri1=ri[k];	
        double and2,and2on2;	
	
      	rid1=rid[k];	
      	and2=rid[k] /(ri[k]+1.0e-25) *an2[k];
      	and2on2=and2 /(an2[k]+1.0e-25);

	    if(an2[k]<d_rib[0] *s2[k])
      	  ifpureshear=1;
        else
	      ifpureshear=0;

        if(ifpureshear==1)	
        {	
      	  interp1d_expabs(and2on2,d_and2on2a1,d_amtaun2a1,d_sma1,
                            d_sha1,d_ssa1,&amtaun2,&sm,&sh,&ss,d_dand2on2,d_rnd2on2);
      	  slq2=(-amtaun2) /((d_b1 *d_b1) *ri[k]+1.0e-25);
        }	
        else{
	        interp2d_expabs(&ri1,&rid1,&slq2,&sm,&sh,&ss,d_rib,d_ridb,
                            d_shb,d_slq2b,d_smb,d_ssb,d_irimax,d_dri,d_rri);
          if(slq2<0.0e0)
      	    return ;	

          if(slq2==0.0e0)
	        ifbelow=1;
        }	
      	
	    lifupper=((ifbelow==0||ifnofsmall==1||((ri1<0.0e0)&&(k<=1)))&&(slq2> 0.0e0));

	    if(lifupper)
        {	
      	  ++ilmldpoint;	
          if(ri1<0.0e0)	
	        ++ilmldpointneg;
        }
      	++ipoint;	
        if(k==0)	
	      ++icall;
	
      	epsy=0.0e0;	
      	lifepsy=false;	
        if(n> 0)	
        {	
          if(slq2==0.0e0)
      	    epsy=0.0e0;
          else
	        epsy=-buoytot /((1.0e0 /((d_b1 *d_b1) *(slq2)))-0.5e0 *sm);

      	  lifepsy=((epsy>=0.0e0)&&lifupper);
          if((epsy<0.0e0)&&lifupper)
          {
            ++iproblem;
            if(ri1<0.0e0)
		        ++inegproblem;
          }
        }
        double akz=0.4 *z[k];
	
        double al=akz *al0 /(al0+akz);
	
      	al2=al *al;	
        if(ifbelow!=1&&!lifepsy)
        {
          if(ri1> 0.0e0)
          {
      	    anlq2=slq2 *ri1;
            if(anlq2> 0.281e0)
            {
      	      al2=0.281e0 /anlq2 * al2;
      	      slq2=0.281e0 /(ri1+1.0e-20);
            }
          }
        }
        if(an2[k]<0.0e0)
      	  epson2_=epson2__;

      else{
	    if(ifnofsmall==1)
      	  eplatidepend=0.0e0;	

        else
	      	eplatidepend=eplatidepend_fun(fabs(coriol),an);	

      	eplatidepend=fmax(eplatidepend,eplatidependmin);	
      	epson2_=epson2__ *eplatidepend;	
      }

      epson2=epson2_;	
      if(ri[k]> 0.0e0)	
      {	
      	rit=(ri[k]+rid[k]) /2.0e0;	
      	ric=(ri[k]-rid[k]) /2.0e0;	
      	ra_r=sqrt((rit *rit)+(ric *ric));
        if(fabs(rit)==0.0e0)
        {	
          if(fabs(ric)==0.0e0)
          {	
      	    theta_r=atan(1.0e0);	
          }	
          else{
	        //theta_r=2.0e0 *atan(1.0e0);
	        theta_r=pi/2.0e0 ;
          }	
        }	
        else{
	        theta_r=atan(ric /rit);
        }	
        //if(fabs(theta_r)> 2.0e0 *atan(1.0e0))	
        if(fabs(theta_r)> (pi/2.0e0))
          return ;	

        //if(theta_r<-atan(1.0e0))	
        if(theta_r<(-pi/4.0e0))
	    //theta_r+=4.0e0 *atan(1.0e0);
	    theta_r+= pi;
	
      	//jtheta_r0=(int)((theta_r+atan(1.0e0)) /d_deltheta_r);	
      	jtheta_r0=(int)((theta_r+(pi/4.0e0)) /d_deltheta_r);
      	itheta_r0=jtheta_r0-n_theta_r_oct;
      	itheta_r1=itheta_r0+1;
      	theta_r0=itheta_r0 *d_deltheta_r;
      	theta_r1=itheta_r1 *d_deltheta_r;
        if((theta_r0<=d_theta_rcrp)&&(theta_r> d_theta_rcrp))
        {	
      	  theta_r=theta_r1;	
      	  theta_r0=theta_r1;	
      	  itheta_r0=itheta_r1;	
      	  ++itheta_r1;	
      	  theta_r1+=d_deltheta_r;
        }	
        else if((theta_r1>=d_theta_rcrn)&&(theta_r<d_theta_rcrn)){
      	  theta_r=theta_r0;	
      	  theta_r1=theta_r0;	
      	  itheta_r1=itheta_r0;	
          --itheta_r0;	
      	  theta_r0+=-d_deltheta_r;
        }	
	
        if((itheta_r1> 3 *n_theta_r_oct)||(itheta_r0<-n_theta_r_oct))
          return ;	

      	deltheta_r_1=theta_r-theta_r0;	
      	delback_ra_r=d_back_ra_r[n_theta_r_oct+itheta_r1]-d_back_ra_r[n_theta_r_oct+itheta_r0];	
      	dback_ra_r_o_dtheta=delback_ra_r /d_deltheta_r;	
      	back_ra_r1=d_back_ra_r[n_theta_r_oct+itheta_r0]+deltheta_r_1 *dback_ra_r_o_dtheta;	
      	ifrafglt=0;

	    if(ifrafgmax==1)
        {	
          if((theta_r<=d_theta_rcrp)||(theta_r>=d_theta_rcrn))	
          {	
      	    if(back_ra_r1> ra_r)	
            {	
      	      ifrafglt=1;	
      	      back_ra_r1=ra_r;	
            }	
          }
        }	
        if(back_ra_r1<0.0e0)
          return ;	

      	back_rit1=cos(theta_r) *back_ra_r1;	
      	back_ric1=sin(theta_r) *back_ra_r1;	
      	back_ri1=back_rit1+back_ric1;	
      	back_rid1=back_rit1-back_ric1;

	    if(ifrafglt==1)
        {	
      	    interp2d_expabs(&back_ri1,&back_rid1,&slq2_back,&sm_back,&sh_back,&ss_back,
                              d_rib,d_ridb,d_shb,d_slq2b,d_smb,d_ssb,d_irimax,d_dri,d_rri);
        }	
        else{
	        deltheta_r_1=theta_r-itheta_r0 *d_deltheta_r;
      	    delsm_back=d_sm_r1[n_theta_r_oct+itheta_r1]-d_sm_r1[n_theta_r_oct+itheta_r0];	
      	    dsm_back_o_dtheta=delsm_back /d_deltheta_r;	
      	    sm_back=d_sm_r1[n_theta_r_oct+itheta_r0]+deltheta_r_1 *dsm_back_o_dtheta;	
      	    delsh_back=d_sh_r1[n_theta_r_oct+itheta_r1]-d_sh_r1[n_theta_r_oct+itheta_r0];	
      	    dsh_back_o_dtheta=delsh_back /d_deltheta_r;	
      	    sh_back=d_sh_r1[n_theta_r_oct+itheta_r0]+deltheta_r_1 *dsh_back_o_dtheta;	
      	    delss_back=d_ss_r1[n_theta_r_oct+itheta_r1]-d_ss_r1[n_theta_r_oct+itheta_r0];	
      	    dss_back_o_dtheta=delss_back /d_deltheta_r;	
      	    ss_back=d_ss_r1[n_theta_r_oct+itheta_r0]+deltheta_r_1 *dss_back_o_dtheta;	
      	    delslq2_back=d_slq2_r1[n_theta_r_oct+itheta_r1]-d_slq2_r1[n_theta_r_oct+itheta_r0];	
      	    dslq2_back_o_dtheta=delslq2_back /d_deltheta_r;	
      	    slq2_back=d_slq2_r1[n_theta_r_oct+itheta_r0]+deltheta_r_1 *dslq2_back_o_dtheta;

        }	
        if(slq2_back<0.0e0)
            return ;

      }	
      if(ri[k]<=0.0e0)	
      {	
      	back_ra_r1=0.0e0;	
      	back_rit1=0.0e0;	
      	back_ric1=0.0e0;	
      	back_ri1=0.0e0;	
      	back_rid1=0.0e0;	
      }	
      if(ri1<=0.0e0)	
      {	
      	sm_back=0.0e0;	
      	sh_back=0.0e0;	
      	ss_back=0.0e0;	
      }	
      if((sm_back<0.0e0)||(sh_back<0.0e0)||(ss_back<0.0e0))
        return ;
      
      tmp_back=0.5e0 *(d_b1 *d_b1) *back_ri1 *slq2_back *epson2;
      v_back=tmp_back *sm_back;	
      t_back=tmp_back *sh_back;	
      s_back=tmp_back *ss_back;
      if((v_back<0.0e0)||(t_back<0.0e0)||(s_back<0.0e0))
        return ;	

      if((ri[k]> 0.0e0)&&((fabs(v_back)<=0.0e0)||(fabs(t_back)<=0.0e0)||(fabs(s_back)<0.0e0)))
        return;

      aldeep=0.0e0;
      double delz,delrh,del2rh;	
	
      if(ifbelow==1)	
      {	
        if(ri1>=0.0e0)	
        {	
      	  tmp=0.5e0 *(d_b1 *d_b1) *ri1 *slq2 *epson2;
        }	
        else if(k>=2)	
        {	
          if(k==n-1)	
          {	
      	    delz=z[k]-z[k-1];	
      	    delrh=rh[k]-rh[k-1];	
      	    del2rh=rh[k]-2.0e0 *rh[k-1]+rh[k-2];
          }	
          else{
	        delz=z[k+1]-z[k-1];
      	    delrh=rh[k+1]-rh[k-1];	
      	    del2rh=rh[k+1]-2.0e0 *rh[k]+rh[k-1];
          }	
          double dzrh=delrh /delz;	
	
          double d2zrh=4.0e0 *del2rh /(delz *delz);
	
          double rdzlndzrh=dzrh /d2zrh;	
	
          double al0deep=0.17e0 *fabs(rdzlndzrh);	
	
          double akz=0.4e0 *z[k];	
	
      	  aldeep=akz *al0deep /(al0deep+akz);	
      	  al2=aldeep *aldeep;
        
	      if(ifpureshear!=1)
          {	
            if(ifshearmin)	
            {	
      	      s2[k]=fmax(s2[k],s2min);	
            }

	        tmp=0.5e0 *d_b1 *al2 *sqrt(s2[k] /(slq2+1.0e-40));
          }	
        }	
        else{
	      if(ifpureshear!=1)
          {
      	    s2[k]=fmax(s2[k],s2min);	
      	    if(lifepsy)
      	      tmp=0.5e0 *epsy /(s2[k]+1.0e-40);
      	    else
              tmp=0.5e0 *d_b1 *al2 *sqrt(s2[k] /(slq2+1.0e-40));
          }	
        }	
	
      }	
      else{
	    if(ifpureshear!=1)
        {	
      	  s2[k]=fmax(s2[k],s2min);	
          if(lifepsy)
      	    tmp=0.5e0 *epsy /(s2[k]+1.0e-40);

          else
	        tmp=0.5e0 *d_b1 *al2 *sqrt(s2[k] /(slq2+1.0e-40));
        }	
      }	
      if(ifpureshear==1)	
	    tmp=0.5e0 *(d_b1 *d_b1) *al2 *sqrt(-an2[k] /amtaun2);

      akm[k]=fmin(tmp *sm+v_back,visc_cbu_limit1);	
      akh[k]=fmin(tmp *sh+t_back,diff_cbt_limit1);	
      aks[k]=fmin(tmp *ss+s_back,diff_cbt_limit1);

//      if(Para_s->rank == 0)
//        printf("[slave %d]:i:%d,j:%d,k:%d,ifbelow:%d,ifpureshear:%d,lifepsy:%d,ifrafglt:%d,tmp:%.16lf,tmp_back:%.16lf,sm_back:%.16lf,v_back:%.16lf,akm:%.16lf,epsy:%.16lf,s2:%.16lf,sm:%.16lf\n",
//			blockIndex,ii,jj,k,ifbelow,ifpureshear,lifepsy,ifrafglt,tmp,tmp_back,sm_back,v_back,akm[k],epsy,s2[k],sm);
      }
    for(k=nb;k<nmax; k++)	
    {

      akm[k]=0.0e0;	
      akh[k]=0.0e0;	
      aks[k]=0.0e0;	
    }	
    
    if(n> 0)	
    {	
      if(akm[0]<wndmix)	
	akm[0]=wndmix;	
	
      if(akh[0]<wndmix)	
	akh[0]=wndmix;	
	
      if(aks[0]<wndmix)	
	aks[0]=wndmix;	
	
    }	
    for(k=0;k<n; k++)	
    {	
      if((akm[k]<0.0e0)||(akh[k]<0.0e0)||(aks[k]<0.0e0))	
      {	
        return ;	
      }	
    }	
      return ;	

}


 void  turb_cu_(turb_Para *Paras){

  //ATHREAD local variable definition
  int swGrid[3];
  int swBlock[3];
  int stride;
  int blockId_x,blockId_y,blockId_z;
  int threadId_x,threadId_y,threadId_z;

  //ATHREAD basic operation
     int i,j,k;
     blockIndex = CRTS_tid;
     turb_Para *Para_s = (turb_Para*)ldm_malloc(sizeof(turb_Para));
  CRTS_dma_get(Para_s, Paras, sizeof(turb_Para));
  //D_COUNT++;
  //CRTS_dma_wait_value(&dma_rply, D_COUNT);

  //cuda grid and block affine translation
  int blockSum   = Para_s->cuGrid[0] * Para_s->cuGrid[1] * Para_s->cuGrid[2];
  int threadSum  = Para_s->cuBlock[0] * Para_s->cuBlock[1] * Para_s->cuBlock[2];
  int quotient   = blockSum / 64; 
  int remainder  = blockSum % 64; 
  int blockStart = quotient * blockIndex;
  if(remainder > 0)
    blockStart += (remainder > blockIndex ? blockIndex : remainder);
  int blockEnd = blockStart + quotient;
  if(remainder > 0)
    blockEnd   += (remainder > blockIndex ? 1 : 0); 

  //double *d_wp12_s = Para_s->d_wp12;
  //double *d_wp13_s = Para_s->d_wp13;
  double wp1[km],wp3[km],wp4[km],wp5[km],wp6[km],wp7[km],wp8[km],wk1[km],wk2[km],wk3[km];
    
  double mldtmp,wp10,wp11;
    
  int d_ncc=Para_s->d_turb_param[7];
    
  double dfricmx=Para_s->d_turb_param[8],dwndmix=Para_s->d_turb_param[9];   //to check

  double ak_tide_s[km];
  for(int block = blockStart; block < blockEnd; block++){
    blockId_z = block / (Para_s->cuGrid[0] * Para_s->cuGrid[1]);
    blockId_y = (block - blockId_z * (Para_s->cuGrid[0] * Para_s->cuGrid[1])) / Para_s->cuGrid[0];
    blockId_x = block % Para_s->cuGrid[0];
      j = blockId_y * Para_s->cuBlock[1];

    //ATHREAD serial operation replace GPU parallel mutli-thread in single thread-block
    for(i = 0; i < threadSum; i++) {


//      if(Para_s->rank == 0)
//          printf("[slave %d]:turb_cu_ i:%d,j:%d, thread, x,y:%d,%d,%d,block,x,y:%d,%d,%d\n",
//                 blockIndex,i,j,threadId,threadId_x,threadId_y,block,blockId_x,blockId_y);

      if(i> 0&&i<imt-1&&j> 0&&j<jmt-1)	
      {	
        unsigned int tid=j *imt+i;	
	
        if(Para_s->d_vit[tid]> 0.5)	
        {	
	        int iwk=Para_s->d_kmt[tid]-1;
	        int iwk1=(Para_s->d_kmt[tid] > 2) ? iwk : 1;
	  
	        for(k=0;k<km; k++)
	        {
                wp1[k]=0.0e0;
	            wp3[k]=0.0e0;
	            wp7[k]=0.0e0;
	            wp8[k]=0.0e0;
	            wp4[k]=0.0e0;
	            wp5[k]=0.0e0;
	            wp6[k]=0.0e0;
	            wk1[k]=0.0e0;
                wk2[k]=0.0e0;
	            wk3[k]=0.0e0;
		        ak_tide_s[k] = 0.0e0;
            }
	  
	        for(k=0;k<iwk; k++)
	        {
	            wp8[k]=-Para_s->d_vit[(k+1) *imt *jmt+tid] *Para_s->d_zkp[k+1] *1.0e+2;
	            wp1[k]=Para_s->d_vit[(k+1) *imt *jmt+tid] *(Para_s->d_at[k *imt *jmt+tid]-
		            (Para_s->d_at[k *imt *jmt+tid]-Para_s->d_at[(k+1) *imt *jmt+tid]) *
		            Para_s->d_dzp[k] /(Para_s->d_dzp[k]+Para_s->d_dzp[k+1]));

		        wp4[k]=Para_s->d_vit[(k+1) *imt *jmt+tid] *Para_s->d_rit[k *imt *jmt+tid];
		        wp5[k]=Para_s->d_vit[(k+1) *imt *jmt+tid] *Para_s->d_ricdttms[k *imt *jmt+tid];
		        wp6[k]=Para_s->d_vit[(k+1) *imt *jmt+tid] *Para_s->d_s2t[k *imt *jmt+tid];
		        wp7[k]=Para_s->d_vit[(k+1) *imt *jmt+tid] *Para_s->d_rict[k *imt *jmt+tid];
		        wp3[k]=Para_s->d_vit[(k+1) *imt *jmt+tid] *(Para_s->d_pdensity[k *imt *jmt+tid]+
				    (Para_s->d_pdensity[(k+1) *imt *jmt+tid]-Para_s->d_pdensity[k *imt *jmt+tid]) *
				    Para_s->d_dzp[k] /(Para_s->d_dzp[k]+Para_s->d_dzp[k+1])) *1.0e-3;
            }
      	    
		    wp10=Para_s->d_vit[tid] *Para_s->d_buoytur[tid] *1.0e+4;
            wp11=Para_s->d_vit[tid] *Para_s->d_buoysol[tid] *1.0e+4;

		    //inout:wk1 wk2 wk3
            turb_device(&wp8[0],&wp1[0],&wp3[0],&wp4[0],&wp5[0],&wp6[0],dfricmx *1.0e+4,dwndmix *1.0e+4,&wp7[0],
			    wp10,wp11,Para_s->d_fcort[tid],&mldtmp,&wk1[0],&wk2[0],&wk3[0],iwk,iwk,km,
			    Para_s->d_amtaun2a1,Para_s->d_and2on2a1,Para_s->d_back_ra_r,Para_s->d_rib,Para_s->d_ridb,
                Para_s->d_sha1,Para_s->d_shb,Para_s->d_sh_r1,Para_s->d_slq2b,Para_s->d_slq2_r1,Para_s->d_sma1,
			    Para_s->d_smb,Para_s->d_sm_r1,Para_s->d_ssa1,Para_s->d_ssb,Para_s->d_ss_r1,
                	    Para_s->d_turb_param,Para_s->d_irimax,i,j);

//	        Para_s->d_amld[tid]=mldtmp *1.0e-2;

              for(k=0;k<iwk; k++)	
              {
                ak_tide_s[k] = back_tidalmixing+mixing_ef *local_mixing_fraction *
	          //Para_s->d_ak_tide[k *imt *jmt+tid]=back_tidalmixing+mixing_ef *local_mixing_fraction *
				Para_s->d_wave_dis[tid] *Para_s->d_fz_tide[k *imt *jmt+tid] /
				(fmax(Para_s->d_rict[k *imt *jmt+tid],1.0e-8) *(wp3[k] *1000.0));

    	        //if(Para_s->d_ak_tide[k *imt *jmt+tid]> max_tidalmixing)
                if(ak_tide_s[k] > max_tidalmixing) {
                  Para_s->d_ak_tide[k * imt * jmt + tid] = max_tidalmixing;
                  ak_tide_s[k] = max_tidalmixing;
                }
                else
                  Para_s->d_ak_tide[k *imt *jmt+tid] = ak_tide_s[k];
	
	          //  Para_s->d_richardson[k *imt *jmt+tid]=Para_s->d_rict[k *imt *jmt+tid];
	           // Para_s->d_fztidal[k *imt *jmt+tid]=Para_s->d_fz_tide[k *imt *jmt+tid];
	           // Para_s->d_wp3_tidal[k *imt *jmt+tid]=wp3[k];
              }

	      for(k=iwk-2;k>=0; k--)
	      {
	         //Para_s->d_ak_tide[k *imt *jmt+tid]=fmin(Para_s->d_ak_tide[k *imt *jmt+tid],Para_s->d_ak_tide[(k+1) *imt *jmt+tid]);
             ak_tide_s[k] = fmin(ak_tide_s[k],ak_tide_s[k+1]);
             Para_s->d_ak_tide[k *imt *jmt+tid] = ak_tide_s[k];
          }

          for(k=0;k<km-1; k++)
	      {
	         //Para_s->d_akmt[k *imt *jmt+tid]=wk1[k] *1.0e-4+Para_s->d_ak_tide[k *imt *jmt+tid] *5.0;
                 Para_s->d_akmt[k *imt *jmt+tid]=wk1[k] *1.0e-4+ak_tide_s[k] *5.0;
                 //Para_s->d_akt[k *imt *jmt+tid]=(wk2[k] *1.0e-4+Para_s->d_ak_tide[k *imt *jmt+tid]) /(double)d_ncc;
                 Para_s->d_akt[k *imt *jmt+tid]+=(wk2[k] *1.0e-4+ak_tide_s[k]) /(double)d_ncc;
                 //Para_s->d_akt[km *imt *jmt+k *imt *jmt+tid]=(wk3[k] *1.0e-4+Para_s->d_ak_tide[k *imt *jmt+tid]) /(double)d_ncc;
                 Para_s->d_akt[km *imt *jmt+k *imt *jmt+tid]+=(wk3[k] *1.0e-4+ak_tide_s[k]) /(double)d_ncc;
//              Para_s->d_wk1[k *imt *jmt+tid] = wk1[k];
//              Para_s->d_wk2[k *imt *jmt+tid] = wk2[k];
//              Para_s->d_wk3[k *imt *jmt+tid] = wk3[k];

	      }
        }
      }	

    }

  }
    ldm_free(Para_s, sizeof(turb_Para));
}

