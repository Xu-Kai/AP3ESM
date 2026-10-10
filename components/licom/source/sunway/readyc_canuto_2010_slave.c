
//ATHREAD include file
#include <crts.h>
#include <slave.h>
//#include <simd.h>
#include <math.h>
#include <stdio.h>
#include <float.h>
#include <string.h>
#include <stdbool.h>
#include "swperf.h" 

#include "param_mod.h"
#include "constant_mod.h"
#include "pconst_mod.h"

//ATHREAD DMA transfer definition 
__thread_local unsigned int D_COUNT = 0;
__thread_local crts_rply_t dma_rply = 0;


typedef struct readyc_canuto2010_kernelPara{
	int cuGrid[3];
  int cuBlock[3];
  int *d_mldlev;
  int *d_kmt;
  double d_pi0_1;
  double *d_ak_tide_mixing;
  double *d_dzp;
  double *d_zkp;
  double *d_ulat;
  double *d_amld;
  double *d_rit;
  double *d_rict;
  double *d_ricdt;
  double *d_s2t;
  double *d_vit;
  double *d_akt;
  double *d_at;
  double *d_to;
  double *d_so;
  double *d_c;
  double *d_wave_dis;
  double *d_fz_tide;
  double *d_akmt;
  //double *d_wk1;
  //double *d_wk2;
  //double *d_wk3;
  //double *d_wp3;
  int rank;
}canuto2010_Para;


int myrank;
//ATHREAD global definition
__thread_local int blockIndex;
//__thread_local Para* Para_s;

__inline__ double  sign(double x, double y){
    return (y) >= 0.0e0 ? fabs(x) : -fabs(x);

}
void mixing_efficiency (double* mix_eff_var, double struct_var, double Ri, double Gm) {
    *mix_eff_var = 0.5e0 * Ri * Gm * struct_var;
    return ;
}

 double dens_k0 ( double tq, double sq,double *d_c_s)  {

     double dens;
     dens = (d_c_s[0]  + (d_c_s[3]      + d_c_s[6] * sq) * sq +
        (d_c_s[2]  +  d_c_s[7] * sq + d_c_s[5] * tq) * tq) * tq +
       (d_c_s[1]  + (d_c_s[4]      + d_c_s[8] * sq) * sq) * sq;
     return dens;
}

void prepare_pi (double Ri, double Rrho, double* out_pi,double pi0_1) {
    double a     = 10.0e0;
    double sig_t = 0.72e0;
    double zero  = 1.0e-10;
    double pi0_2 = 1.0e0/3.0e0;
    double pi0_3 = sig_t;
    double pi0_4 = pi0_1;
    double pi0_5 = sig_t;
    double Ko =   1.66e0;

   // double  pi0_1_s   =   1.0e0/(sqrt(27.0e0/5.0e0*(Ko * Ko * Ko))*(1+1.0e0/sig_t));

    if (Ri > zero && Rrho > zero) {
        out_pi[0] = pi0_1 / (1.0e0 + Ri * Rrho / (a + Rrho));
        // out_pi[1] = PI0_2_ / (1.0 + Ri) * (1.0 + 2.0 * Ri * Rrho
        //     / (1.0 + Rrho * Rrho));
        out_pi[1] = pi0_2 /(1.0e0 + Ri)  * (1.0e0 + 2.0e0 * Ri * Rrho
                                         / (1.0e0 + Rrho * Rrho));
        out_pi[2] = pi0_3;
        out_pi[3] = pi0_4 / (1.0e0 + Ri / (1.0e0 + a * Rrho));
        out_pi[4] = pi0_5;
    }

    if (Ri > zero && Rrho <= zero) {
        // out_pi[0] = PI0_1_ / (1.0 + Ri);
        out_pi[0] = pi0_1 / (1.0e0 + Ri);
        out_pi[1] = pi0_2;
        out_pi[2] = pi0_3;
        out_pi[3] = pi0_4 / (1.0e0 + Ri);
        //out_pi[3] = pi0_4 * tmp;
        out_pi[4] = pi0_5;
    }

    if (Ri <= zero) {
        out_pi[0] = pi0_1;
        out_pi[1] = pi0_2;
        out_pi[2] = pi0_3;
        out_pi[3] = pi0_4;
        out_pi[4] = pi0_5;
    }
    return ;
}
void Rf_calc (double* Rf_out, double Ri_in, double Rrho_in,
              double struct_h, double struct_m, double r_in) {
    *Rf_out = Ri_in * struct_h / struct_m * (1.0 - Rrho_in / r_in);
    return ;
}

void struct_function (double Ri, double Rrho, double *pi_n,
                      double Gm, double* struct_m, double* struct_h, double* struct_s,
                      double* r_out, double* r_out_new) {

    // hwy_need_check
    double x = Ri * Gm * pow(1.0 - Rrho, -1);
    // const double p = pi_n[3] * pi_n[4] - pi_n[3] * pi_n[1] * (1.0 + Rrho);
    // const double q = pi_n[0] * pi_n[1] * (1.0 + Rrho) - pi_n[0] * pi_n[2] * Rrho;
    //double pi_n1_times_1addRrho = pi_n[1] * (1.0 + Rrho);
    double p = pi_n[3] * pi_n[4] - pi_n[3] * pi_n[1] * (1.0e0 + Rrho);
    double q = pi_n[0] * pi_n[1] * (1.0e0 + Rrho) - pi_n[0] * pi_n[2] * Rrho;

    // hwy_need_check
    // const double r     = (pi_n[3] / pi_n[0]) / Rrho * (1.0 + q * x) / (1.0 + p * x);
    // hwy_need_check
    // const double r_new = (pi_n[3] / pi_n[0])        * (1.0 + q * x) / (1.0 + p * x);
    double r =(pi_n[3] / pi_n[0])/Rrho * (1.0e0 + q * x) / (1.0e0 + p * x);
    double r_new = (pi_n[3] / pi_n[0]) * (1.0e0 + q * x) / (1.0e0 + p * x);
    //double r     = r_new / Rrho;
    *r_out     = r;
    *r_out_new = r_new;

    //double tmp = 1.0 - 1.0 / r;
    // const double Ah = pi_n[3] / (1.0 + p * x + pi_n[3] * pi_n[1] * x * (1.0 - 1.0 / r));
    double Ah = pi_n[3] / (1.0e0 + p * x + pi_n[3] * pi_n[1] * x * (1.0e0 - 1.0e0 / r));
    double As = Ah / r_new;

    // const double Am1 = 4.0 / 5.0 - (pi_n[3] - pi_n[0] + (pi_n[0] - 1.0 / 150.0)
    //     * (1.0 - 1.0 / r)) * x * Ah;
    // const double Am2 = 10.0 + (pi_n[3] - pi_n[0] * Rrho) * x + 1.0 / 50.0 * Gm;
    double Am1 = 4.0e0 / 5.0e0 -  (pi_n[3] - pi_n[0] + (pi_n[0] - 1.0e0 / 150.0e0)
                                             * (1.0e0 - 1.0e0 / r)) * x * Ah;
    double Am2 = 10.0e0 + (pi_n[3] - pi_n[0] * Rrho) * x + 1.0e0/50.0e0 * Gm;
    double Am  = Am1 / Am2;

    // const double LX = (1.0 - 1.0 / r) * x * Ah;
    double LX = (1.0e0 - 1.0e0 / r) * x * Ah;

    // const double ratio = 2.0 / 3.0 / (1.0 + 2.0 / 15.0 * LX + 1.0 / 10.0 * Am * Gm);
    double ratio = 2.0e0 / 3.0e0 / (1.0e0 + 2.0e0 / 15.0e0 * LX + 1.0e0 / 10.0e0 * Am * Gm);

    // struct_m = ratio * Am1 / Am2;
    *struct_m = ratio * Am1 / Am2;
    *struct_h = ratio * Ah;
    *struct_s = ratio * As;
    return ;
}

void shengjin_calc (double* array, double *real_x) {
    const double epsln1 = 1.0e-6;

    const double sa = array[0];
    const double sb = array[1];
    const double sc = array[2];
    const double sd = array[3];

    *real_x = 1.0;

    if (fabs(sa) <= epsln1 && fabs(sb) <= epsln1) {
        if ((-sd / sc) < 0.0) {
            // printf ("No real positive root at position 1\n");
            return;
        } else {
            // printf ("Only one real positive root exists 1!\n");
            *real_x = -sd / sc;
            return;
        }
    }

    double delta;
    double x[3];
    if (fabs(sa) <= epsln1 && fabs(sb) >= epsln1) {
        delta = sc * sc - 4.0 * sb * sd;
        if (delta < 0.0) {
            // printf ("No real roots in this secondary order equation 2!\n");
            return;
        } else if (delta >= 0.0) {
            // const double tmp = std::sqrt(delta);
            // x[0] = (-sc + tmp) / (2.0 * sb);
            // x[1] = (-sc - tmp) / (2.0 * sb);
            //const double tmpD = sqrt(delta);
            //const double tmp2 = 1.0 / (2.0 * sb);
            x[0] = (-sc + sqrt(delta)) / (2.0e0 * sb);
            x[1] = (-sc - sqrt(delta)) / (2.0e0 * sb);
            if (x[0] < 0.0 && x[1] < 0.0) {
                // std::cout << "No real positive root at position 2" << std::endl;
                return;
            } else {
                // std::cout << "Two real with one positive at least 2!" << std::endl;
                if (x[0] > 0.0 && x[1] > 0.0) {
                    *real_x = fmin(x[0], x[1]);
                } else if (x[0] >  0.0 && x[1] <= 0.0) {
                    *real_x = x[0];
                } else if (x[0] <= 0.0 && x[1]  > 0.0) {
                    *real_x = x[1];
                }
                return;
            }
        }
    }

    const double A = sb * sb - 3.0e0 * sa * sc;
    const double B = sb * sc - 9.0e0 * sa * sd;
    const double C = sc * sc - 3.0e0 * sb * sd;

    if (fabs(A) <= epsln1 && fabs(B) <= epsln1) {
        if ((- sc / sb) < 0.0) {
            // std::cout << "No real positive root at position 3" << std::endl;
            return;
        } else {
            // std::cout << "Only one real positive root exists 3!" << std::endl;
            *real_x = - sc / sb;
            return;
        }
    }

    delta = B * B - 4.0e0 * A * C;

    double y[2];
//    const double three_times_sa = 3.0 * sa;
//    const double one_three_times_sa = 1.0 / three_times_sa;
    if (delta > 0.0) {
        // double tmp = std::sqrt(delta);
        // y[0] = A * sb + 3.0 * sa * (-B + tmp) / 2.0;
        // y[1] = A * sb + 3.0 * sa * (-B - tmp) / 2.0;
        // const double tmp = sqrt(delta);
        y[0] = A * sb + 3.0e0 * sa * (-B + sqrt(delta)) / 2.0e0;
        y[1] = A * sb + 3.0e0 * sa * (-B - sqrt(delta)) / 2.0e0;
        // real_x = (-sb - sign(pow(fabs(y[0]), 1.0 / 3.0), y[0])
        //               - sign(pow(fabs(y[1]), 1.0 / 3.0), y[1])) / (3.0 * sa);
        *real_x = (-sb - sign(pow(fabs(y[0]), 1.0e0 / 3.0e0), y[0])
                   - sign(pow(fabs(y[1]), 1.0e0 / 3.0e0), y[1]))
                  / (3.0e0 * sa);
        if (*real_x < 0.0) {
            *real_x = 1.0;
            return;
        } else {
            // std::cout << "Only one real positive root exists 4!" << std::endl;
            // real_x = real_x;
            return;
        }
    }
    // const double T = (2.0 * A * sb - 3.0 * sa * B) / (2.0 * pow(A, 1.5));
    const double T = (2.0e0 * A * sb - 3.0e0 * sa * B) / (2.0e0 * pow(A, 1.5e0));
    const double the  = acos(T);
    //const double tmp1 = sqrt(A);
    //const double tmp2 = cos(the / 3.0);
    // const double tmp3 = std::sqrt(3.0) * std::sin(the / 3.0);
    //const double tmp3 = SQRT_3 * sin(the / 3.0);
    //const double tmp3 = sqrt(3.0) * sin(the / 3.0);
    // x[0] = (-sb - 2.0  *  tmp1 *  tmp2        ) / 3.0 / sa;
    // x[1] = (-sb +         tmp1 * (tmp2 + tmp3)) / 3.0 / sa;
    // x[2] = (-sb +         tmp1 * (tmp2 - tmp3)) / 3.0 / sa;
    x[0] = (-sb - 2.0e0  *  sqrt(A) *  cos(the / 3.0e0)        ) / 3.0e0 / sa ;
    x[1] = (-sb +         sqrt(A) * (cos(the / 3.0e0) + sqrt(3.0e0) * sin(the / 3.0e0))) / 3.0e0 / sa;
    x[2] = (-sb +         sqrt(A) * (cos(the / 3.0e0) - sqrt(3.0e0) * sin(the / 3.0e0))) / 3.0e0 / sa;

    *real_x = 500.0e0;
    for (int k = 0; k < 3; ++k) {
        if (x[k] >= 0.0 && x[k] <= 1000.0) {
            if (x[k] < *real_x) {
                *real_x = x[k];
            }
        }
    }

    // real_x = std::numeric_limits<double>::max();
    *real_x = DBL_MAX;
    for (int i = 0; i < 3; ++i) {
        if (x[i] > 0.0) {
            *real_x = fmin (*real_x, x[i]);
        }
    }

    *real_x = fmin (*real_x, 2.0e3);
    *real_x = fmax (*real_x, 1.0);

    return ;
}

void dyn_time_scale_calc (double Ri, double Rrho,
                          double* pi_n, double* Gm, double* cubic_coef) {

    double A1, A2, A3, A4, A5, A6;
    double temp1, temp2, temp3, temp4, temp5;
    double temp6, temp7, temp8, temp9, temp10;
    double temp11;

    temp1 = pi_n[0] * pi_n[3] * (pi_n[3] - pi_n[0] * Rrho);
    temp2 = pi_n[1] * (15.0e0 * pi_n[2] + 7.0e0) * (Rrho * Rrho + 1.0e0);
    temp3 = (14.0e0 * (pi_n[1] - pi_n[2]) - 15.0e0 * (pi_n[2] * pi_n[2])) * Rrho;
    // hwy_need_check
    double one_minus_Rrho = 1.0e0 - Rrho;
    temp4 = 150.0e0 * (one_minus_Rrho * one_minus_Rrho * one_minus_Rrho);

    A1 = temp1 * (temp2 + temp3) / temp4;

    temp1 = pi_n[0] * pi_n[3];
    temp2 = pi_n[1] * (210.0e0 * pi_n[0] - 150.0e0 * pi_n[2] + 7.0e0)
            * (Rrho * Rrho + 1.0e0);
    temp3 = 14.0e0 * (pi_n[1] - pi_n[2]) * (1.0e0 + 15.0e0 * pi_n[0] + 15.0e0 * pi_n[3]);
    temp4 = 150.0e0 * (pi_n[2] * pi_n[2]);
    temp5 = (temp3 + temp4) * Rrho;
    temp6 = 210.0e0 * pi_n[1] * (pi_n[3] - pi_n[0]);

    // hwy_need_check
    temp7 = 9000.0e0 * (one_minus_Rrho * one_minus_Rrho);

    A2 = temp1 * (temp2 + temp5 + temp6) / temp7;

    temp1 = pi_n[0];
    temp2 = 5.0e0 * pi_n[1] * pi_n[3] * (30.0e0 * pi_n[2] + 17.0e0);
    temp3 = pi_n[0] * (15.0e0 * pi_n[2] + 7.0e0);
    temp4 = Rrho * Rrho + 1.0e0;
    temp5 = temp1 * (temp2 + temp3) * temp4;
    temp6 = -(15.0e0 * pi_n[2] + 7.0) * (pi_n[0] * pi_n[0] - pi_n[3] * pi_n[3]);
    temp7 = 10.0e0 * pi_n[0] * pi_n[2] * pi_n[3] * (15.0e0 * pi_n[2] + 17.0);
    temp8 = 15.0e0 * pi_n[1]          * (pi_n[0] * pi_n[0] + pi_n[3] * pi_n[3]);
    temp9 = 14.0e0 * pi_n[0] * pi_n[3] * (1.0e0 - 10.0e0 * pi_n[1]);
    temp10 = -(temp7 + temp8 + temp9) * Rrho;
    temp11 = 150.0e0 * (one_minus_Rrho * one_minus_Rrho);  // hwy_need_check

    A3 = (temp5 + temp6 + temp10) / temp11;

    temp1 = 150.0e0 * (pi_n[0] * pi_n[2] + pi_n[1] * pi_n[3]);
    temp2 = -7.0e0 * pi_n[0] * (1.0e0 + 30.0e0 * pi_n[0]);
    temp3 = (temp1 + temp2) * Rrho;
    temp4 = -150.0e0 * (pi_n[0] * pi_n[1] + pi_n[2] * pi_n[3]);
    temp5 = 7.0e0 * pi_n[3] * (1.0e0 + 30.0e0 * pi_n[3]);
    temp6 = 9000.0e0 * one_minus_Rrho;  // hwy_need_check

    A4 = (temp3 + temp4 + temp5) / temp6;

    temp1 = -30.0e0 * (pi_n[0] * pi_n[2] + pi_n[1] * pi_n[3]);
    temp2 = -17.0e0 * pi_n[0];
    temp3 = (temp1 + temp2) * Rrho;
    temp4 = 30.0e0 * (pi_n[0] * pi_n[1] + pi_n[2] * pi_n[3]);
    temp5 = 17.0e0 * pi_n[3];
    temp6 = 30.0e0 * one_minus_Rrho;  // hwy_need_check

    A5 = (temp3 + temp4 + temp5) / temp6;

    A6 = - 1.0e0 / 60.0e0;

    cubic_coef[0] = A1 * (Ri * Ri * Ri) + A2 * (Ri * Ri);
    cubic_coef[1] = A3 * (Ri * Ri) + A4 * Ri;
    cubic_coef[2] = A5 * Ri + A6;
    cubic_coef[3] = 1.0e0;
    shengjin_calc(cubic_coef, Gm);
    return ;
}
void find_mix_layer_depth(double *mld_out,int *mld_lev_out,double *den,double delta,double *zlev_in,int n){
    double den_m;
    double zlev[km];
    int k;
    double VERY_SMALL = 1.0e-15;
    //zlev(:) = abs(zlev_in(:))
    for(k=0; k<km; k++) {
        //if(k<(n-1))
            zlev[k] = fabs(zlev_in[k]);
        //else
        //    zlev[k] = 0.0;
    }

    for(k=0; k<n; ++k) {
        if (fabs(den[k] - den[0]) > delta) {
            den_m = den[0] - sign(delta, den[0] - den[k]);
            *mld_out = zlev[k] + (zlev[k - 1] - zlev[k]) *
                    (den_m - den[k]) / (den[k - 1] - den[k] + VERY_SMALL);
            *mld_lev_out = k+1;
            break;
        }
        *mld_out = zlev[n-1];
        *mld_lev_out = n;
    }

//        *mld_out = zlev[n-1];
//        *mld_lev_out = n;
}

void thermocline_mixing_coeff_calc (double* Km_out, double* Kh_out, double* Ks_out,
                                    double mix_eff_m, double mix_eff_h, double mix_eff_s,
                                    double n2_in,     double lat_in) {
    // const double Tke_n2 = 0.288e-4;
    double N0 = 5.24e-3;
    double f30 = 2.0e0*omega*sin(30.0e0*DEGtoRAD);
    double Tke_n2 = 0.288e-4;
    double f_lat = fabs (2.0e0 * omega * sin(lat_in * DEGtoRAD) + 1.0e-15);

    double cond1 = fmax(sqrt(n2_in) / f_lat, 1.0e0);

    // const double L_lat = (f_lat * acosh(cond1)) / L_LAT_CONST_ * Tke_n2;
    //double L_lat = (f_lat * acosh(cond1)) * L_LAT_CONST;
    double L_lat = (f_lat*acosh(cond1))/(f30*acosh(N0/(f30)));

    *Km_out = L_lat * Tke_n2 * mix_eff_m;
    *Kh_out = L_lat * Tke_n2 * mix_eff_h;
    *Ks_out = L_lat * Tke_n2 * mix_eff_s;
    return ;
}
void mixed_layer_TKE_calc (double* TKE_out, double mld_in,
                           double Gm_in, double s2_in, double lev_in,
                           double Rf_in, double Rf_inf_in) {
    double l0 = 0.17e0 * mld_in;
    double KARMAN = 0.4e0;
    double lB = KARMAN * fabs(lev_in * l0) / (l0 + KARMAN * fabs(lev_in));
//    double tmp = pow(pow(1.0 - Rf_in / Rf_inf_in, 4.0), 1.0 / 3.0);
//    double l = lB * (pow(pow(1.0e0 - Rf_in / Rf_inf_in, 4.0), 1.0e0 / 3.0e0));
    double tmp=1.0e0 - Rf_in / Rf_inf_in;
    double l = lB * (pow(tmp*tmp*tmp*tmp, 1.0e0 / 3.0e0));
    const double B1 = 21.6e0;
    // TKE_out = B1 * B1 * pow(Gm_in, - 3.0 / 2.0)
    //     * l * l * pow(s2_in, 3.0/2.0);
    *TKE_out = (B1 * B1) * pow(Gm_in, (-3.0e0/2.0e0)) * (l * l) * pow(s2_in, (3.0e0/2.0e0));
    return ;
}

void canuto_2010_interface (
        int ii,
        int jj,
        double *Km_out,  // wk1
        double *Kh_out,  // wk2
        double *Ks_out,  // wk3
        double  *mld_out, // amld
        double  *rho_in,   // wp3
        double  *ri_in,   // wp4
        double  *rrho_in, // wp5
        double  *n2_in,   // wp7
        double  *s2_in,   // wp6
        double  lat_in,  // ulat[iblock][j][i] / DegToRad
        double  *lev_in,  // wp8
        int     num_lev, // kmt[iblock][j][i]
        double  *tmp_dzp,// 0.008 * dzp[k] * dzp[k]
        double pi0_1,
        int *d_mldlev)
{

    double Ri[km];
    double Rrho[km];
    double Ri_low  = - 1.0e+10;
    double Ri_high =   1.0e+10;
    double Rrho_bound = 1.0e-3;
    //int     *mld_lev = 0;
    double delta_s = 0.03e0;

    find_mix_layer_depth(mld_out,d_mldlev,rho_in,delta_s,lev_in,num_lev);
    for(int k=0; k < km; k++) {
        Ri[k]   = ri_in[k];
        Rrho[k] = rrho_in[k];
        if(Ri[k] > Ri_high)
            Ri[k] = Ri_high;
        if(Ri[k] < Ri_low)
            Ri[k] = Ri_low;

    }
    for(int k=0; k < num_lev - 1; k++) {

        if (fabs(rrho_in[k] - 1.0) < Rrho_bound) {
            if (rrho_in[k] >= 1.0) {
                Rrho[k] = 1.001;
            } else if (rrho_in[k] <= 1.0) {
                Rrho[k] = 0.999;
            }
        }
    }

    double out_pi[5];
    double cube[4];
    double struct_m_inf = 0.0;
    double struct_h_inf = 0.0;
    double struct_s_inf = 0.0;
    double R_inf = 0.0;
    double Rnew_inf = 0.0;
    double Rf_inf = 0.0;
    double struct_m = 0.0;
    double struct_h = 0.0;
    double struct_s = 0.0;
    double R = 0.0;
    double Rnew = 0.0;
    double mix_eff_m = 0.0;
    double mix_eff_h = 0.0;
    double mix_eff_s = 0.0;
    double Rf = 0.0;
    double TKE_mld = 0.0;
    double Gm[km];

    for(int k=0; k < num_lev - 1; k++) {
        prepare_pi (1.0e10, Rrho[k], out_pi,pi0_1);

        dyn_time_scale_calc (1.0e10, Rrho[k], out_pi, &Gm[k], cube);

        struct_function(1.0e10, Rrho[k], out_pi, Gm[k],
                        &struct_m_inf, &struct_h_inf, &struct_s_inf, &R_inf, &Rnew_inf);

        Rf_calc(&Rf_inf, 1.0e10, Rrho[k], struct_h_inf, struct_m_inf, Rnew_inf);

        prepare_pi(Ri[k], Rrho[k], out_pi,pi0_1);

        dyn_time_scale_calc(Ri[k], Rrho[k], out_pi, &Gm[k], cube);

        struct_function(Ri[k], Rrho[k], out_pi, Gm[k],
                        &struct_m, &struct_h, &struct_s, &R, &Rnew);

        mixing_efficiency(&mix_eff_m, struct_m, Ri[k], Gm[k]);
        mixing_efficiency(&mix_eff_h, struct_h, Ri[k], Gm[k]);
        mixing_efficiency(&mix_eff_s, struct_s, Ri[k], Gm[k]);

        Rf_calc(&Rf, Ri[k], Rrho[k], struct_h, struct_m, Rnew);

        double VERY_SMALL = 1.0e-15;

        if (k <= (*d_mldlev - 1)) {
            mixed_layer_TKE_calc (&TKE_mld, *mld_out, Gm[k],
                                  s2_in[k], lev_in[k], Rf, Rf_inf);

            //double tmp = TKE_mld / (n2_in[k] + VERY_SMALL);
            Km_out[k] = mix_eff_m * TKE_mld / (n2_in[k] + VERY_SMALL);
            Kh_out[k] = mix_eff_h * TKE_mld / (n2_in[k] + VERY_SMALL);
            Ks_out[k] = mix_eff_s * TKE_mld / (n2_in[k] + VERY_SMALL);
        }

        if (k > (*d_mldlev - 1)) {
            thermocline_mixing_coeff_calc(
                    &Km_out[k], &Kh_out[k], &Ks_out[k],
                    mix_eff_m, mix_eff_h, mix_eff_s, n2_in[k], lat_in);
        }

        if (k < 3) {
            Km_out[k] = fmax(1.0e-3, Km_out[k]);
            Kh_out[k] = fmax(1.0e-3, Kh_out[k]);
            Ks_out[k] = fmax(1.0e-3, Ks_out[k]);
        } else {
            Km_out[k] = fmax(1.0e-4, Km_out[k]);
            Kh_out[k] = fmax(1.0e-5, Kh_out[k]);
            Ks_out[k] = fmax(1.0e-5, Ks_out[k]);
        }

        Km_out[k]= fmin(1.2e-1, Km_out[k]);
        Kh_out[k] = fmin(1.2e-1, Kh_out[k]);
        Ks_out[k] = fmin(1.2e-1, Ks_out[k]);

        if (n2_in[k] < VERY_SMALL) {
            Kh_out[k] = fmin(tmp_dzp[k], 8.0);
            Ks_out[k] = fmin(tmp_dzp[k], 8.0);
            Km_out[k] = fmin(tmp_dzp[k], 8.0);
        }
    }
    return;
}

 void  readyc_canuto2010_(canuto2010_Para *Paras){

  int blockId_x,blockId_y,blockId_z;
  int threadId_x,threadId_y,threadId_z;
  double wk1[km];
  double wk2[km];
  double wk3[km];
  double wp3[km];
  double wp4[km];
  double wp5[km];
  double wp6[km];
  double wp7[km];
  double wp8[km];
  double d_c[9];
  double tmp_dzp[km];
  double tq =0.0;
  double sq = 0.0;
  double amld = 0.0;
  double ak_tide_mixing[km];
  //ATHREAD basic operation
  int i,j,k;
  blockIndex = CRTS_tid;
  canuto2010_Para* Para_s = (canuto2010_Para*)ldm_malloc(sizeof(canuto2010_Para));
  CRTS_dma_iget(Para_s, Paras, sizeof(canuto2010_Para), &dma_rply);
  D_COUNT++;
  CRTS_dma_wait_value(&dma_rply, D_COUNT);
  double to = Para_s->d_to[0];
  double so = Para_s->d_so[0];
//  CRTS_dma_iget(d_c, Paras->d_c, sizeof(double) * 9, &dma_rply);
//  D_COUNT++;
//  CRTS_dma_wait_value(&dma_rply, D_COUNT);
  d_c[0] = Para_s->d_c[0];
  d_c[1] = Para_s->d_c[km];
  d_c[2] = Para_s->d_c[km*2];
  d_c[3] = Para_s->d_c[km*3];
  d_c[4] = Para_s->d_c[km*4];
  d_c[5] = Para_s->d_c[km*5];
  d_c[6] = Para_s->d_c[km*6];
  d_c[7] = Para_s->d_c[km*7];
  d_c[8] = Para_s->d_c[km*8];

     for(k=0; k<km; k++)
         tmp_dzp[k] = 0.008 * Paras->d_dzp[k] * Paras->d_dzp[k];

     int d_mldlev_s = 0;
  double d_pi0_1 = Para_s->d_pi0_1;
     myrank = Para_s->rank;

    int tile  = imt ;
//    int avaiLdm = 98304; //available ldm: 96kB
    int avaiLdm = get_allocatable_size(); //available ldm: 64kB
    int varNum = 3; // variable number of k direction inconsistent access
    int tileBlk = 1; //blocks number in i direction
    int tileLimit = avaiLdm/(varNum*km*sizeof(double));
    int remTile = 0;
    if(tile > tileLimit) {
        remTile = tile % tileLimit;
        if(remTile == 0)
            tileBlk = tile / tileLimit;
        else
            tileBlk = tile / tileLimit + 1;
        tile = tileLimit;
    }
    int start_blk;
    int tile_size  = tile * sizeof(double);
    int tileSum   = tile * km * sizeof(double);
    int size_stride = (imt * jmt - tile) * sizeof(double);
    int tileV;
    int tileSumV;
    int size_strideV;
    // target variable for k direction inconsistent access

    double* ldm_akmt       = (double*) ldm_malloc (tileSum);
    double* ldm_akt1       = (double*) ldm_malloc (tileSum);
    double* ldm_akt2       = (double*) ldm_malloc (tileSum);
    // each cpe is responsible for i*k data
    for(int j = blockIndex; j < jmt; j += 64) {

        // blocks in i direction
        for(int idx=0; idx < tileBlk; idx ++) {

            start_blk = j * imt + idx * tile;
            //distinguish tile and remainder for stride dma
            if((remTile != 0) && (idx == tileBlk -1)) {
                tileV = remTile;
                tileSumV   = remTile * km * sizeof(double);
                tile_size  = remTile * sizeof(double);
                size_stride = (imt * jmt - remTile) * sizeof(double);
           }
           else {
                tileV = tile;
                tileSumV   = tile * km * sizeof(double);
                tile_size  = tile * sizeof(double);
                size_stride = (imt * jmt - tile) * sizeof(double);
           }
//           if(Para_s->rank == 71){
//               printf("[slave %d]: j:%d,idx:%d,tileV:%d,remTile:%d,tileBlk:%d\n",
//                    blockIndex,j,idx,tileV,remTile,tileBlk);
//           }
            athread_dma_iget_stride(ldm_akmt, &Para_s->d_akmt[start_blk], tileSumV, tile_size, size_stride, &dma_rply);
            D_COUNT++;
            athread_dma_wait_value(&dma_rply, D_COUNT);

            athread_dma_iget_stride(ldm_akt1, &Para_s->d_akt[start_blk], tileSumV, tile_size, size_stride, &dma_rply);
            D_COUNT++;
            athread_dma_wait_value(&dma_rply, D_COUNT);

            athread_dma_iget_stride(ldm_akt2, &Para_s->d_akt[km*imt*jmt+start_blk],tileSumV, tile_size, size_stride, &dma_rply);
            D_COUNT++;
            athread_dma_wait_value(&dma_rply, D_COUNT);

            for (int i = 0; i < tileV; ++i) {
                int index = start_blk + i;
                int ii = idx * tile + i;
                //if(i> 0&&i<imt-1&&j> 0&&j<jmt-1)
                //if(((idx == 0 && i>0) || (idx !=0)) && ((idx == tileBlk -1 && i<tileV - 1) || (idx != tileBlk -1)) && j> 0&&j<jmt-1){
                if(ii > 0 && ii < imt -1 && j> 0&&j<jmt-1){
                    if(Para_s->d_vit[index]> 0.5){
                        int iwk=Para_s->d_kmt[index];

                        for(k=0; k<km; k++) {
                            wk1[k] = 0.0;
                            wk2[k] = 0.0;
                            wk3[k] = 0.0;
                            wp3[k] = 0.0;
                            wp4[k] = 0.0;
                            wp5[k] = 0.0;
                            wp6[k] = 0.0;
                            wp7[k] = 0.0;
                            if(k<km-1)
                                wp8[k] = -Para_s->d_vit[(k + 1) * imt * jmt + index] * Paras->d_zkp[k + 1];
                               // wp8[k] = -ldm_vit[(k + 1) * tileV + i] * Paras->d_zkp[k + 1];
                            else
                                wp8[k] = 0.0;
                            ak_tide_mixing[k] = 0.0;
                        }

                        for(k=0; k<iwk; k++){
                            //wp4(k)=vit(i,j,k+1,iblock)*RIT(i,j,k,iblock)
                            wp4[k] = Para_s->d_vit[(k + 1) * imt * jmt + index] * Para_s->d_rit[k * imt * jmt + index];
                            //wp5(k)=vit(i,j,k+1,iblock)*RICDT(i,j,k,iblock)
                            wp5[k] = Para_s->d_vit[(k + 1) * imt * jmt + index] * Para_s->d_ricdt[k * imt * jmt + index];
                            //wp6(k)=vit(i,j,k+1,iblock)*S2T(i,j,k,iblock)
                            wp6[k] = Para_s->d_vit[(k + 1) * imt * jmt + index] * Para_s->d_s2t[k * imt * jmt + index];
                            //wp7(k)=vit(i,j,k+1,iblock)*RICT(i,j,k,iblock)
                            wp7[k] = Para_s->d_vit[(k + 1) * imt * jmt + index] * Para_s->d_rict[k * imt * jmt + index];
                            //tq = at(i,j,k,1,iblock) - to(1)
                            tq = Para_s->d_at[k * imt * jmt + index] - to;
                            //sq = at(i,j,k,2,iblock) - so(1)
                            sq = Para_s->d_at[km * imt * jmt + k * imt * jmt + index] - so;
                            wp3[k] = dens_k0(tq,sq,d_c)+1000.0;
                            if(wp7[k] < 0.0) wp7[k] = 0.0;
                        }
                        amld = 0.0;
                        d_mldlev_s = 0;
                        canuto_2010_interface (i,j,wk1, wk2, wk3, &amld,wp3,
                                       wp4, wp5, wp7, wp6, Para_s->d_ulat[index] / DEGtoRAD,
                                       wp8, iwk, tmp_dzp,d_pi0_1,&d_mldlev_s);
                        Para_s->d_amld[index] = amld;
//                        Para_s->d_mldlev[index] = d_mldlev_s;
                        for(k=0; k<iwk-1; k++){
                            ak_tide_mixing[k]=back_tidalmixing+mixing_ef*local_mixing_fraction *
                                    Para_s->d_wave_dis[index]* Para_s->d_fz_tide[k * imt * jmt + index]/
                                    (fmax(Para_s->d_rict[k * imt * jmt + index],1.0e-8)*wp3[k]);

                            ak_tide_mixing[k]=fmin(ak_tide_mixing[k],max_tidalmixing);

                        }

                        //DO K = INT(KMT(I,J,IBLOCK))-2,1,-1     !yuzp-2016/11/23
                        for(k=iwk-3; k>=0; k--)
                            ak_tide_mixing[k]=fmin(ak_tide_mixing[k],ak_tide_mixing[k+1]);

                        for(k=0; k<km; k++) {
//                            Para_s->d_wk1[k*imt*jmt + index] = wk1[k];
//                            Para_s->d_wk2[k*imt*jmt + index] = wk2[k];
//                            Para_s->d_wk3[k*imt*jmt + index] = wk3[k];
//                            Para_s->d_wp3[k*imt*jmt + index] = wp3[k];
                            //Para_s->d_akmt[k*imt*jmt + index] = wk1[k] + ak_tide_mixing[k] * 5.0;
                            ldm_akmt[k * tileV + i] = wk1[k] + ak_tide_mixing[k] * 5.0;
                            //Para_s->d_akt[k*imt*jmt+index] = Para_s->d_akt[k*imt*jmt+index] + wk2[k] + ak_tide_mixing[k];
                            ldm_akt1[k * tileV + i] = ldm_akt1[k * tileV + i] + wk2[k] + ak_tide_mixing[k];
                            //Para_s->d_akt[km*imt*jmt+k*imt*jmt+index] =
                            //Para_s->d_akt[km*imt*jmt+k*imt*jmt+index] +wk3[k] + ak_tide_mixing[k];
                            ldm_akt2[k * tileV + i] =ldm_akt2[k * tileV + i] + wk3[k] + ak_tide_mixing[k];
//                    Para_s->d_akmt[k*imt*jmt + tid] = Para_s->d_akmt[k*imt*jmt + tid] + ak_tide_mixing[k] * 5.0;
//                    Para_s->d_akt[k*imt*jmt+tid] = Para_s->d_akt[k*imt*jmt+tid] + ak_tide_mixing[k];
//                    Para_s->d_akt[km*imt*jmt+k*imt*jmt+tid]  =
//                            Para_s->d_akt[km*imt*jmt+k*imt*jmt+tid]  + ak_tide_mixing[k];

                           // Para_s->d_akmt[k*imt*jmt + index] = fmin(Para_s->d_akmt[k*imt*jmt + index], 8.0e-3 *
                           //         Paras->d_dzp[k] * Paras->d_dzp[k]);
                            ldm_akmt[k * tileV + i] = fmin(ldm_akmt[k * tileV + i],tmp_dzp[k]);

                           // Para_s->d_akt[k*imt*jmt+index] = fmin(Para_s->d_akt[k*imt*jmt+index], 8.0e-3 *
                           //                 Para_s->d_dzp[k] * Para_s->d_dzp[k]);
                            ldm_akt1[k * tileV + i] = fmin(ldm_akt1[k * tileV + i],tmp_dzp[k]);

                           // Para_s->d_akt[km*imt*jmt+k*imt*jmt+index] = fmin(Para_s->d_akt[km*imt*jmt+k*imt*jmt+index],
                           //                 8.0e-3 *  Para_s->d_dzp[k] * Para_s->d_dzp[k]);
                            ldm_akt2[k * tileV + i] = fmin(ldm_akt2[k * tileV + i],tmp_dzp[k]);
                        }

                    }
                }
            }

            athread_dma_iput_stride (&Para_s->d_akmt[start_blk],
                                     ldm_akmt, tileSumV, tile_size, size_stride,&dma_rply);
            D_COUNT++;
            athread_dma_wait_value(&dma_rply, D_COUNT);

            athread_dma_iput_stride (&Para_s->d_akt[start_blk],
                                     ldm_akt1, tileSumV, tile_size, size_stride,&dma_rply);
            D_COUNT++;
            athread_dma_wait_value(&dma_rply, D_COUNT);

            athread_dma_iput_stride (&Para_s->d_akt[km*imt*jmt+start_blk],
                                     ldm_akt2, tileSumV, tile_size, size_stride,&dma_rply);
            D_COUNT++;
            athread_dma_wait_value(&dma_rply, D_COUNT);

        }

    }
    ldm_free(ldm_akmt, tileSum);
    ldm_free(ldm_akt1, tileSum);
    ldm_free(ldm_akt2, tileSum);
    ldm_free(Para_s, sizeof(canuto2010_Para));
}

void  readyc_canuto2010_old_(canuto2010_Para *Paras){

    int blockId_x,blockId_y,blockId_z;
    int threadId_x,threadId_y,threadId_z;
    double wk1[km];
    double wk2[km];
    double wk3[km];
    double wp3[km];
    double wp4[km];
    double wp5[km];
    double wp6[km];
    double wp7[km];
    double wp8[km];
    double d_c[9];
    double tmp_dzp[km];
    double tq =0.0;
    double sq = 0.0;
    double amld = 0.0;
    double ak_tide_mixing[km];
    //ATHREAD basic operation
    int i,j,k;
    blockIndex = CRTS_tid;
    canuto2010_Para* Para_s = (canuto2010_Para*)ldm_malloc(sizeof(canuto2010_Para));
    CRTS_dma_iget(Para_s, Paras, sizeof(canuto2010_Para), &dma_rply);
    D_COUNT++;
    CRTS_dma_wait_value(&dma_rply, D_COUNT);
    double to = Para_s->d_to[0];
    double so = Para_s->d_so[0];
//  CRTS_dma_iget(d_c, Paras->d_c, sizeof(double) * 9, &dma_rply);
//  D_COUNT++;
//  CRTS_dma_wait_value(&dma_rply, D_COUNT);
    d_c[0] = Para_s->d_c[0];
    d_c[1] = Para_s->d_c[km];
    d_c[2] = Para_s->d_c[km*2];
    d_c[3] = Para_s->d_c[km*3];
    d_c[4] = Para_s->d_c[km*4];
    d_c[5] = Para_s->d_c[km*5];
    d_c[6] = Para_s->d_c[km*6];
    d_c[7] = Para_s->d_c[km*7];
    d_c[8] = Para_s->d_c[km*8];


    int d_mldlev_s = 0;
    double d_pi0_1 = Para_s->d_pi0_1;
    myrank = Para_s->rank;
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


    for(int block = blockStart; block < blockEnd; block++){
        blockId_z = block / (Para_s->cuGrid[0] * Para_s->cuGrid[1]);
        blockId_y = (block - blockId_z * (Para_s->cuGrid[0] * Para_s->cuGrid[1])) / Para_s->cuGrid[0];
        blockId_x = block % Para_s->cuGrid[0];
        j = blockId_y * Para_s->cuBlock[1];



        for(i = 0; i < threadSum; i++) {
            int tid=j *imt+i;
            if(i> 0&&i<imt-1&&j> 0&&j<jmt-1){
                if(Para_s->d_vit[tid]> 0.5){

                    int iwk=Para_s->d_kmt[tid];

                    for(k=0; k<km; k++) {
                        wk1[k] = 0.0;
                        wk2[k] = 0.0;
                        wk3[k] = 0.0;
                        wp3[k] = 0.0;
                        wp4[k] = 0.0;
                        wp5[k] = 0.0;
                        wp6[k] = 0.0;
                        wp7[k] = 0.0;
                        if(k<km-1)
                            wp8[k] = -Para_s->d_vit[(k + 1) * imt * jmt + tid] * Paras->d_zkp[k + 1];
                        else
                            wp8[k] = 0.0;
                        ak_tide_mixing[k] = 0.0;
                    }
                    for(k=0; k<iwk; k++){
                        if(k < iwk -1) {
                            //wp4(k)=vit(i,j,k+1,iblock)*RIT(i,j,k,iblock)
                            wp4[k] = Para_s->d_vit[(k + 1) * imt * jmt + tid] * Para_s->d_rit[k * imt * jmt + tid];
                            //wp5(k)=vit(i,j,k+1,iblock)*RICDT(i,j,k,iblock)
                            wp5[k] = Para_s->d_vit[(k + 1) * imt * jmt + tid] * Para_s->d_ricdt[k * imt * jmt + tid];
                            //wp6(k)=vit(i,j,k+1,iblock)*S2T(i,j,k,iblock)
                            wp6[k] = Para_s->d_vit[(k + 1) * imt * jmt + tid] * Para_s->d_s2t[k * imt * jmt + tid];
                            //wp7(k)=vit(i,j,k+1,iblock)*RICT(i,j,k,iblock)
                            wp7[k] = Para_s->d_vit[(k + 1) * imt * jmt + tid] * Para_s->d_rict[k * imt * jmt + tid];
                        }
                        tq = Para_s->d_at[k*imt*jmt + tid] - to;
                        sq = Para_s->d_at[km*imt*jmt + k*imt*jmt + tid] - so;
                        wp3[k] = dens_k0(tq,sq,d_c)+1000.0;
                        if(wp7[k] < 0.0) wp7[k] = 0.0;
                        tmp_dzp[k] = 0.008 * Paras->d_dzp[k] * Paras->d_dzp[k];
                    }
                    amld = 0.0;
                    d_mldlev_s = 0;
                    canuto_2010_interface (i,j,wk1, wk2, wk3, &amld,wp3,
                                           wp4, wp5, wp7, wp6, Para_s->d_ulat[tid] / DEGtoRAD,
                                           wp8, iwk, tmp_dzp,d_pi0_1,&d_mldlev_s);
                  Para_s->d_amld[tid] = amld;
//                Para_s->d_mldlev[tid] = d_mldlev_s;
                    for(k=0; k<iwk-1; k++){
                        ak_tide_mixing[k]=back_tidalmixing+mixing_ef*local_mixing_fraction *
                                                           Para_s->d_wave_dis[tid]* Para_s->d_fz_tide[k*imt*jmt+tid]/
                                                           (fmax(Para_s->d_rict[k*imt*jmt + tid],1.0e-8)*wp3[k]);
                        ak_tide_mixing[k]=fmin(ak_tide_mixing[k],max_tidalmixing);

                    }

                    //DO K = INT(KMT(I,J,IBLOCK))-2,1,-1     !yuzp-2016/11/23
                    for(k=iwk-3; k>=0; k--)
                        ak_tide_mixing[k]=fmin(ak_tide_mixing[k],ak_tide_mixing[k+1]);

                    for(k=0; k<km; k++) {
//                    Para_s->d_wk1[k*imt*jmt + tid] = wk1[k];
//                    Para_s->d_wk2[k*imt*jmt + tid] = wk2[k];
//                    Para_s->d_wk3[k*imt*jmt + tid] = wk3[k];
//                    Para_s->d_wp3[k*imt*jmt + tid] = wp3[k];
                        Para_s->d_akmt[k*imt*jmt + tid] = wk1[k] + ak_tide_mixing[k] * 5.0;
                        Para_s->d_akt[k*imt*jmt+tid] = Para_s->d_akt[k*imt*jmt+tid] + wk2[k] + ak_tide_mixing[k];
                        Para_s->d_akt[km*imt*jmt+k*imt*jmt+tid] =
                                Para_s->d_akt[km*imt*jmt+k*imt*jmt+tid] +wk3[k] + ak_tide_mixing[k];

//                    Para_s->d_akmt[k*imt*jmt + tid] = Para_s->d_akmt[k*imt*jmt + tid] + ak_tide_mixing[k] * 5.0;
//                    Para_s->d_akt[k*imt*jmt+tid] = Para_s->d_akt[k*imt*jmt+tid] + ak_tide_mixing[k];
//                    Para_s->d_akt[km*imt*jmt+k*imt*jmt+tid]  =
//                            Para_s->d_akt[km*imt*jmt+k*imt*jmt+tid]  + ak_tide_mixing[k];

                        Para_s->d_akmt[k*imt*jmt + tid] = fmin(Para_s->d_akmt[k*imt*jmt + tid], 8.0e-3 *
                                Paras->d_dzp[k] * Paras->d_dzp[k]);
                        Para_s->d_akt[k*imt*jmt+tid] = fmin(Para_s->d_akt[k*imt*jmt+tid], 8.0e-3 *
                                                                                          Para_s->d_dzp[k] * Para_s->d_dzp[k]);
                        Para_s->d_akt[km*imt*jmt+k*imt*jmt+tid] = fmin(Para_s->d_akt[km*imt*jmt+k*imt*jmt+tid],
                                                                       8.0e-3 *  Para_s->d_dzp[k] * Para_s->d_dzp[k]);
                    }

                }
            }

        }

    }
    ldm_free(Para_s, sizeof(canuto2010_Para));
}
