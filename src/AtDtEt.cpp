#include <stdio.h>
#include <string>
#include <vector>
#include <math.h>
#include <iostream>

//#include <Rcpp.h>
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace std;

/*============GLOBAL VARIABLES===========*/
/*=====END OF GLOBAL VARIABLES===========*/

/*============FUNCTION DEFINITIONS*/

RcppExport SEXP loglik_AtDtEt_esp_c(SEXP b_a, SEXP b_d, SEXP b_e, SEXP pheno_m, SEXP pheno_d, SEXP B_a_m, SEXP B_a_d, SEXP B_d_m, SEXP B_d_d, SEXP B_e_m, SEXP B_e_d) {
  
  arma::vec ba = as<arma::vec>(b_a); // Number of columns
  arma::vec bd = as<arma::vec>(b_d);
  arma::vec be = as<arma::vec>(b_e);
  
  arma::vec p_m = as<arma::vec>(pheno_m);
  arma::vec p_d = as<arma::vec>(pheno_d);
  arma::mat b_a_m = as<arma::mat>(B_a_m);
  arma::mat b_a_d = as<arma::mat>(B_a_d);
  arma::mat b_d_m = as<arma::mat>(B_d_m);
  arma::mat b_d_d = as<arma::mat>(B_d_d);
  arma::mat b_e_m = as<arma::mat>(B_e_m);
  arma::mat b_e_d = as<arma::mat>(B_e_d);
  
  int numt_m = p_m.n_elem;
  int numt_d = p_d.n_elem;
  int nump_m = numt_m/2;
  int nump_d = numt_d/2;
  int num_a = ba.n_elem;
  int num_d = bd.n_elem;
  int num_e = be.n_elem;
  arma::mat k_m(2,2);
  arma::mat k_d(2,2);
  arma::mat k_d_a(2,2);
  // k_d_a << 1 << 0.5 << arma::endr << 0.5 << 1 << arma::endr;
  k_d_a.fill(0.5);
  k_d_a(0,0) = 1;
  k_d_a(1,1) = 1;
  arma::vec  v = arma::ones<arma::vec>(2);
  arma::mat diag = arma::diagmat(v);
  
  k_m.fill(1);
  k_d.fill(0.25);
  k_d(0,0) = 1;
  k_d(1,1) = 1;
  
  double YSY_m = 0;
  double YSY_d = 0;
  double D_m = 0;
  double D_d = 0;
  
  arma::mat V_m(2,2);
  arma::mat inv_V_m(2,2);
  
  for(int i = 0; i<nump_m; i++)
  {
    int start = 2*i;
    int end = start + 1;
    double temp_com_a = exp(as_scalar(b_a_m.row(start)*ba));
    if(num_a==1)
    {
      temp_com_a = exp(ba(0));
    }
    double temp_com_d = exp(as_scalar(b_d_m.row(start)*bd));
    if(num_d==1)
    {
      temp_com_d = exp(bd(0));
    }
    double temp_com_e = exp(as_scalar(b_e_m.row(start)*be));
    if(num_e==1)
    {
      temp_com_e = exp(be(0));
    }
    V_m = temp_com_e*diag+temp_com_d*k_m+temp_com_a*k_m;
    double det_m = V_m(0,0)*V_m(0,0)-V_m(0,1)*V_m(0,1);
    //inv_V_m = inv(V_m);
    inv_V_m = V_m/det_m;
    inv_V_m(0,1) = (-1)*inv_V_m(0,1);
    inv_V_m(1,0) = inv_V_m(0,1);
    arma::vec p_m_v({p_m[start], p_m[end]});
    arma::rowvec p_m_r({p_m[start], p_m[end]});
    arma::mat temp = p_m_r*inv_V_m*p_m_v;
    YSY_m = YSY_m + as_scalar(temp);
    D_m = D_m + log(det_m);
  }
  
  for(int i = 0; i<nump_d; i++)
  {
    int start = 2*i;
    int end = start + 1;
    double temp_com_a = exp(as_scalar(b_a_d.row(start)*ba));
    if(num_a==1)
    {
      temp_com_a = exp(ba(0));
    }
    double temp_com_d = exp(as_scalar(b_d_d.row(start)*bd));
    if(num_d==1)
    {
      temp_com_d = exp(bd(0));
    }
    double temp_com_e = exp(as_scalar(b_e_d.row(start)*be));
    if(num_e==1)
    {
      temp_com_e = exp(be(0));
    }
    arma::mat V_d = temp_com_e*diag+temp_com_d*k_d+temp_com_a*k_d_a;
    //arma::mat inv_V_d = inv(V_d);
    arma::mat inv_V_d(2,2);
    double det_d = V_d(0,0)*V_d(0,0)-V_d(0,1)*V_d(0,1);
    inv_V_d = V_d/det_d;
    inv_V_d(0,1) = (-1)*inv_V_d(0,1);
    inv_V_d(1,0) = inv_V_d(0,1);
    arma::vec p_d_v({p_d[start], p_d[end]});
    arma::rowvec p_d_r({p_d[start], p_d[end]});
    arma::mat temp = p_d_r*inv_V_d*p_d_v;
    YSY_d = YSY_d + as_scalar(temp);
    D_d = D_d + log(det_d);
  }
  
  double res = (D_m + YSY_m + D_d + YSY_d)/2;
  
  return(Rcpp::wrap(res));
}


RcppExport SEXP gr_AtDtEt_esp_c(SEXP b_a, SEXP b_d, SEXP b_e, SEXP pheno_m, SEXP pheno_d, SEXP B_a_m, SEXP B_a_d, SEXP B_d_m, SEXP B_d_d, SEXP B_e_m, SEXP B_e_d) {
  
  arma::vec ba = as<arma::vec>(b_a); // Number of columns
  arma::vec bd = as<arma::vec>(b_d);
  arma::vec be = as<arma::vec>(b_e);
  
  arma::vec p_m = as<arma::vec>(pheno_m);
  arma::vec p_d = as<arma::vec>(pheno_d);
  arma::mat b_a_m = as<arma::mat>(B_a_m);
  arma::mat b_a_d = as<arma::mat>(B_a_d);
  arma::mat b_d_m = as<arma::mat>(B_d_m);
  arma::mat b_d_d = as<arma::mat>(B_d_d);
  arma::mat b_e_m = as<arma::mat>(B_e_m);
  arma::mat b_e_d = as<arma::mat>(B_e_d);
  
  int numt_m = p_m.n_elem;
  int numt_d = p_d.n_elem;
  int nump_m = numt_m/2;
  int nump_d = numt_d/2;
  int num_a = ba.n_elem;
  int num_d = bd.n_elem;
  int num_e = be.n_elem;
  arma::mat k_m(2,2);
  arma::mat k_d(2,2);
  arma::mat k_d_a(2,2);
  k_d_a.fill(0.5);
  k_d_a(0,0) = 1;
  k_d_a(1,1) = 1;
  arma::vec  v = arma::ones<arma::vec>(2);
  arma::mat diag = arma::diagmat(v);
  
  k_m.fill(1);
  k_d.fill(0.25);
  k_d(0,0) = 1;
  k_d(1,1) = 1;
  
  arma::vec var_bd_m = arma::zeros<arma::vec>(num_d);
  arma::vec var_ba_m = arma::zeros<arma::vec>(num_a);
  arma::vec var_be_m = arma::zeros<arma::vec>(num_e);
  
  arma::mat V_m(2,2);
  arma::mat inv_V_m(2,2);
  
  for(int i = 0; i<nump_m; i++)
  {
    int start = 2*i;
    int end = start + 1;
    double temp_com_a = exp(as_scalar(b_a_m.row(start)*ba));
    if(num_a==1)
    {
      temp_com_a = exp(ba(0));
    }
    double temp_com_d = exp(as_scalar(b_d_m.row(start)*bd));
    if(num_d==1)
    {
      temp_com_d = exp(bd(0));
    }
    double temp_com_e = exp(as_scalar(b_e_m.row(start)*be));
    if(num_e==1)
    {
      temp_com_e = exp(be(0));
    }
    V_m = temp_com_e*diag+temp_com_d*k_m+temp_com_a*k_m;
    inv_V_m = inv(V_m);
    arma::vec p_m_v({p_m[start], p_m[end]});
    arma::rowvec p_m_r({p_m[start], p_m[end]});
    arma::mat inv_ph_m = inv_V_m*p_m_v;
    arma::mat ph_inv_m = p_m_r*inv_V_m;
    
    arma::mat temp_prod = inv_V_m*k_m;
    //double temp_c = exp(as_scalar(b_c_m.row(start)*bc))*(arma::sum(temp_prod.diag())-as_scalar(ph_inv_m*k_m*inv_ph_m));
    double temp_d = arma::sum(temp_prod.diag())-as_scalar(ph_inv_m*k_m*inv_ph_m);
    double temp_a = temp_d;
    //if(num_c>1)
    //{
    temp_d = temp_com_d*temp_d;
    //}
    for(int j = 0; j < num_d; j++)
    {
      var_bd_m(j) = var_bd_m(j) + b_d_m(start,j)*temp_d;
    }
    // double temp_a = arma::sum(temp_prod.diag())-as_scalar(ph_inv_m*k_m*inv_ph_m);
    //if(num_a>1)
    //{
    temp_a = temp_com_a*temp_a;
    //}
    for(int j = 0; j < num_a; j++)
    {
      var_ba_m(j) = var_ba_m(j) + b_a_m(start,j)*temp_a;
    }
    double temp_e = arma::sum(inv_V_m.diag())-as_scalar(ph_inv_m*inv_ph_m);
    //if(num_e>1)
    //{
    temp_e = temp_com_e*temp_e;
    //}
    for(int j = 0; j < num_e; j++)
    {
      var_be_m(j) = var_be_m(j) + b_e_m(start,j)*temp_e;
    }
    
  }
  
  arma::vec var_bd_d = arma::zeros<arma::vec>(num_d);
  arma::vec var_ba_d = arma::zeros<arma::vec>(num_a);
  arma::vec var_be_d = arma::zeros<arma::vec>(num_e);
  
  for(int i = 0; i<nump_d; i++)
  {
    int start = 2*i;
    int end = start + 1;
    double temp_com_a = exp(as_scalar(b_a_d.row(start)*ba));
    if(num_a==1)
    {
      temp_com_a = exp(ba(0));
    }
    double temp_com_d = exp(as_scalar(b_d_d.row(start)*bd));
    if(num_d==1)
    {
      temp_com_d = exp(bd(0));
    }
    double temp_com_e = exp(as_scalar(b_e_d.row(start)*be));
    if(num_e==1)
    {
      temp_com_e = exp(be(0));
    }
    arma::mat V_d = temp_com_e*diag+temp_com_d*k_d+temp_com_a*k_d_a;
    arma::mat inv_V_d = inv(V_d);
    arma::vec p_d_v({p_d[start], p_d[end]});
    arma::rowvec p_d_r({p_d[start], p_d[end]});
    arma::mat inv_ph_d = inv_V_d*p_d_v;
    arma::mat ph_inv_d = p_d_r*inv_V_d;
    
    arma::mat temp_prod = inv_V_d*k_d;
    //double temp_c = exp(as_scalar(b_c_d.row(start)*bc))*(arma::sum(temp_prod.diag())-as_scalar(ph_inv_d*k_d*inv_ph_d));
    double temp_d = arma::sum(temp_prod.diag())-as_scalar(ph_inv_d*k_d*inv_ph_d);
    //if(num_c>1)
    //{
    temp_d = temp_com_d*temp_d;
    //}
    for(int j = 0; j < num_d; j++)
    {
      var_bd_d(j) = var_bd_d(j) + b_d_d(start,j)*temp_d;
    }
    temp_prod = inv_V_d*k_d_a;
    double temp_a = arma::sum(temp_prod.diag())-as_scalar(ph_inv_d*k_d_a*inv_ph_d);
    //if(num_a>1)
    //{
    temp_a = temp_com_a*temp_a;
    //}
    for(int j = 0; j < num_a; j++)
    {
      var_ba_d(j) = var_ba_d(j) + b_a_d(start,j)*temp_a;
    }
    double temp_e = arma::sum(inv_V_d.diag())-as_scalar(ph_inv_d*inv_ph_d);
    //if(num_e>1)
    //{
    temp_e = temp_com_e*temp_e;
    //}
    for(int j = 0; j < num_e; j++)
    {
      var_be_d(j) = var_be_d(j) + b_e_d(start,j)*temp_e;
    }
  }
  
  arma::vec d_var_d = (var_bd_m + var_bd_d)/2;
  arma::vec d_var_a = (var_ba_m + var_ba_d)/2;
  arma::vec d_var_e = (var_be_m + var_be_d)/2;
  
  arma::vec res(num_e+num_d+num_a);
  for(int i = 0; i < num_a; i++)
  {
    res(i) = d_var_a(i);
  }
  for(int i = 0; i < num_d; i++)
  {
    res(num_a+i) = d_var_d(i);
  }
  for(int i = 0; i < num_e; i++)
  {
    res(num_a+num_d+i) = d_var_e(i);
  }
  
  return(Rcpp::wrap(res));
}


RcppExport SEXP hessian_AtDtEt_esp_c(SEXP b_a, SEXP b_d, SEXP b_e, SEXP pheno_m, SEXP pheno_d, SEXP B_a_m, SEXP B_a_d, SEXP B_d_m, SEXP B_d_d, SEXP B_e_m, SEXP B_e_d) 
{
  
  arma::vec p_m = as<arma::vec>(pheno_m);
  arma::vec p_d = as<arma::vec>(pheno_d);
  arma::mat b_a_m = as<arma::mat>(B_a_m);
  arma::mat b_a_d = as<arma::mat>(B_a_d);
  arma::mat b_d_m = as<arma::mat>(B_d_m);
  arma::mat b_d_d = as<arma::mat>(B_d_d);
  arma::mat b_e_m = as<arma::mat>(B_e_m);
  arma::mat b_e_d = as<arma::mat>(B_e_d);
  
  arma::vec ba = as<arma::vec>(b_a);
  arma::vec bd = as<arma::vec>(b_d);
  arma::vec be = as<arma::vec>(b_e);
  
  int numt_m = p_m.n_elem;
  int numt_d = p_d.n_elem;
  int nump_m = numt_m/2;
  int nump_d = numt_d/2;
  int num_a = b_a_m.n_cols;
  int num_d = b_d_m.n_cols;
  int num_e = b_e_m.n_cols;
  arma::mat k_m(2,2);
  arma::mat k_1d(2,2);
  arma::mat k_2d(2,2);
  arma::mat k_3d(2,2);
  arma::mat k_4d(2,2);
  // k_1d << 0.5 << 0 << arma::endr << 0 << 0.5 << arma::endr;
  k_1d.zeros();
  k_1d(0,0) = 0.5;
  k_1d(1,1) = 0.5;
  arma::vec  v = arma::ones<arma::vec>(2);
  arma::mat diag = arma::diagmat(v);
  
  k_m.fill(1);
  k_2d.fill(0.5);
  k_3d.zeros();
  k_3d(0,0) = 0.75;
  k_3d(1,1) = 0.75;
  k_4d.fill(0.25);
  
  arma::rowvec  gsd_a2m2a = arma::zeros<arma::rowvec>(nump_m);
  arma::rowvec  gsd_d2m2a = arma::zeros<arma::rowvec>(nump_m);
  arma::rowvec  gsd_d2m2d = arma::zeros<arma::rowvec>(nump_m);
  arma::rowvec  gsd_dd = arma::zeros<arma::rowvec>(numt_d);
  arma::rowvec  gsd_2dd2 = arma::zeros<arma::rowvec>(nump_d);
  arma::rowvec  gsd_3dd3 = arma::zeros<arma::rowvec>(nump_d);
  arma::rowvec  gsd_4dd4 = arma::zeros<arma::rowvec>(nump_d);
  arma::rowvec  gsd_2d2 = arma::zeros<arma::rowvec>(nump_d);
  arma::rowvec  gsd_3d2 = arma::zeros<arma::rowvec>(nump_d);
  arma::rowvec  gsd_3d3 = arma::zeros<arma::rowvec>(nump_d);
  
  arma::rowvec  gsd_e2mm2a = arma::zeros<arma::rowvec>(nump_m);
  arma::rowvec  gsd_d2mm2e = arma::zeros<arma::rowvec>(nump_m);
  arma::rowvec  gsd_emme = arma::zeros<arma::rowvec>(numt_m);
  arma::rowvec  gsd_adda = arma::zeros<arma::rowvec>(numt_d);
  arma::rowvec  gsd_dddd = arma::zeros<arma::rowvec>(numt_d);
  arma::rowvec  gsd_addd = arma::zeros<arma::rowvec>(numt_d);
  arma::rowvec  gsd_adde = arma::zeros<arma::rowvec>(numt_d);
  arma::rowvec  gsd_edde = arma::zeros<arma::rowvec>(numt_d);
  arma::rowvec  gsd_ddde = arma::zeros<arma::rowvec>(numt_d);
  arma::rowvec  gsd_a2dd2e = arma::zeros<arma::rowvec>(nump_d);
  arma::rowvec  gsd_d3dd3e = arma::zeros<arma::rowvec>(nump_d);
  
  for(int i = 0; i<nump_m; i++)
  {
    int start = 2*i;
    int end = start + 1;
    double ba_m = exp(as_scalar(b_a_m.row(start)*ba));
    double bd_m = exp(as_scalar(b_d_m.row(start)*bd));
    double be_m = exp(as_scalar(b_e_m.row(start)*be));
    arma::mat V_m = be_m*diag + ba_m*k_m + bd_m*k_m;
    arma::mat inv_V_m = inv(V_m);
    arma::vec p_m_v({p_m[start], p_m[end]});
    arma::rowvec p_m_r({p_m[start], p_m[end]});
    double accu_temp = arma::accu(inv_V_m*k_m*inv_V_m);
    gsd_a2m2a(i) = ba_m*ba_m*accu_temp;
    gsd_d2m2a(i) = ba_m*accu_temp*bd_m;
    gsd_d2m2d(i) = bd_m*bd_m*accu_temp;
    arma::mat temp_mm = inv_V_m*inv_V_m;
    gsd_emme(start) = be_m*be_m*temp_mm(0,0);
    gsd_emme(end) = be_m*be_m*temp_mm(1,1);
    accu_temp = arma::accu(temp_mm);
    gsd_e2mm2a(i) = ba_m*be_m*accu_temp;
    gsd_d2mm2e(i) = bd_m*be_m*accu_temp;
  }
  
  for(int i = 0; i<nump_d; i++)
  {
    int start = 2*i;
    int end = start + 1;
    double ba_d = exp(as_scalar(b_a_d.row(start)*ba));
    double bd_d = exp(as_scalar(b_d_d.row(start)*bd));
    double be_d = exp(as_scalar(b_e_d.row(start)*be));
    double ba2_d = ba_d*ba_d;
    double bad_d = ba_d*bd_d;
    double bd2_d = bd_d*bd_d;
    double bae_d = ba_d*be_d;
    double be2_d = be_d*be_d;
    double bde_d = be_d*bd_d;
    arma::mat V_d = be_d*diag + ba_d*(k_1d+k_2d) + bd_d*(k_3d+k_4d);
    arma::mat inv_V_d = inv(V_d);
    arma::vec p_d_v({p_d[start], p_d[end]});
    arma::rowvec p_d_r({p_d[start], p_d[end]});
    arma::mat temp_dd = inv_V_d*inv_V_d;
    gsd_dd(start) = ba2_d*temp_dd(0,0);
    gsd_dd(end) = ba2_d*temp_dd(1,1);
    gsd_adde(start) = bae_d*temp_dd(0,0);
    gsd_adde(end) = bae_d*temp_dd(1,1);
    gsd_dddd(start) = bd2_d*temp_dd(0,0);
    gsd_dddd(end) = bd2_d*temp_dd(1,1);
    gsd_addd(start) = bad_d*temp_dd(0,0);
    gsd_addd(end) = bad_d*temp_dd(1,1);
    gsd_edde(start) = be2_d*temp_dd(0,0);
    gsd_edde(end) = be2_d*temp_dd(1,1);
    gsd_ddde(start) = bde_d*temp_dd(0,0);
    gsd_ddde(end) = bde_d*temp_dd(1,1);
    double accu_temp = arma::accu(temp_dd);
    gsd_2dd2(i) = ba2_d*0.5*accu_temp;
    gsd_3dd3(i) = bad_d*0.5*accu_temp;
    gsd_4dd4(i) = bd2_d*0.375*accu_temp;
    gsd_a2dd2e(i) = bae_d*0.5*accu_temp;
    gsd_d3dd3e(i) = bde_d*0.25*accu_temp;
    double temp_si2i = arma::accu(inv_V_d*k_2d*inv_V_d);
    gsd_2d2(i) = ba2_d*0.5*temp_si2i;
    gsd_3d2(i) = bad_d*0.25*temp_si2i;
    gsd_3d3(i) = bd2_d*0.125*temp_si2i;
  }
  
  arma::mat gsd_max(num_a+num_d+num_e,num_a+num_d+num_e);
  arma::mat b_a_m_h(nump_m,num_a);
  arma::mat b_d_m_h(nump_m,num_d);
  arma::mat b_e_m_h(nump_m,num_e);
  arma::mat b_a_d_h(nump_d,num_a);
  arma::mat b_d_d_h(nump_d,num_d);
  arma::mat b_e_d_h(nump_d,num_e);
  for(int i = 0; i < nump_m; i++)
  {
    b_a_m_h.row(i) = b_a_m.row(2*i);
    b_d_m_h.row(i) = b_d_m.row(2*i);
    b_e_m_h.row(i) = b_e_m.row(2*i);
  }
  for(int i = 0; i < nump_d; i++)
  {
    b_a_d_h.row(i) = b_a_d.row(2*i);
    b_d_d_h.row(i) = b_d_d.row(2*i);
    b_e_d_h.row(i) = b_e_d.row(2*i);
  }
  //arma::mat b_a_m_t = b_a_m.t();
  //arma::mat b_c_m_t = b_c_m.t();
  arma::mat b_e_m_t = b_e_m.t();
  arma::mat b_a_d_t = b_a_d.t();
  arma::mat b_e_d_t = b_e_d.t();
  arma::mat b_d_d_t = b_d_d.t();
  arma::mat b_a_m_ht = b_a_m_h.t();
  arma::mat b_d_m_ht = b_d_m_h.t();
  arma::mat b_e_m_ht = b_e_m_h.t();
  arma::mat b_a_d_ht = b_a_d_h.t();
  arma::mat b_d_d_ht = b_d_d_h.t();
  arma::mat b_e_d_ht = b_e_d_h.t();
  arma::mat temp1 = b_a_m_ht;
  arma::mat temp4 = b_d_m_ht;
  arma::mat temp6 = b_d_m_ht;
  arma::mat temp8 = b_e_m_ht;
  arma::mat temp13 = b_d_m_ht;
  for(int i = 0; i < nump_m; i++)
  {
    temp1.col(i) = gsd_a2m2a(i)*b_a_m_ht.col(i);
    temp4.col(i) = gsd_d2m2a(i)*b_d_m_ht.col(i);
    temp6.col(i) = gsd_d2m2d(i)*b_d_m_ht.col(i);
    temp8.col(i) = gsd_e2mm2a(i)*b_e_m_ht.col(i);
    temp13.col(i) = gsd_d2mm2e(i)*b_d_m_ht.col(i);
  }
  arma::mat temp3 = b_a_d_ht;
  arma::mat temp5 = b_d_d_ht;
  arma::mat temp7 = b_d_d_ht;
  arma::mat temp10 = b_e_d_ht;
  arma::mat temp14 = b_d_d_ht;
  for(int i = 0; i < nump_d; i++)
  {
    temp3.col(i) = (gsd_2dd2(i)+gsd_2d2(i))*b_a_d_ht.col(i);
    temp5.col(i) = (gsd_3dd3(i)+gsd_3d2(i))*b_d_d_ht.col(i);
    temp7.col(i) = (gsd_3d3(i)+gsd_4dd4(i))*b_d_d_ht.col(i);
    temp10.col(i) = gsd_a2dd2e(i)*b_e_d_ht.col(i);
    temp14.col(i) = gsd_d3dd3e(i)*b_d_d_ht.col(i);
  }
  arma::mat temp2 = b_a_d_t;
  arma::mat temp9 = b_e_d_t;
  arma::mat temp12 = b_e_d_t;
  arma::mat temp15 = b_d_d_t;
  arma::mat temp16 = b_d_d_t;
  arma::mat temp17 = b_d_d_t;
  for(int i = 0; i < numt_d; i++)
  {
    temp2.col(i) = gsd_dd(i)*b_a_d_t.col(i);
    temp15.col(i) = gsd_dddd(i)*b_d_d_t.col(i);
    temp16.col(i) = gsd_addd(i)*b_d_d_t.col(i);
    temp17.col(i) = gsd_ddde(i)*b_d_d_t.col(i);
    temp9.col(i) = gsd_adde(i)*b_e_d_t.col(i);
    temp12.col(i) = gsd_edde(i)*b_e_d_t.col(i);
  }
  arma::mat temp11 = b_e_m_t;
  for(int i = 0; i < numt_m; i++)
  {
    temp11.col(i) = gsd_emme(i)*b_e_m_t.col(i);
  }
  gsd_max.submat(0,0,num_a-1,num_a-1) = temp1*b_a_m_h + 0.25*temp2*b_a_d + temp3*b_a_d_h;
  gsd_max.submat(num_a,0,num_a+num_d-1,num_a-1) = temp4*b_a_m_h + 0.375*temp16*b_a_d + temp5*b_a_d_h;
  gsd_max.submat(0,num_a,num_a-1,num_a+num_d-1) = trans(gsd_max.submat(num_a,0,num_a+num_d-1,num_a-1));
  gsd_max.submat(num_a,num_a,num_a+num_d-1,num_a+num_d-1) = temp6*b_d_m_h + 0.5625*temp15*b_d_d + temp7*b_d_d_h;
  gsd_max.submat(num_a+num_d,0,num_a+num_d+num_e-1,num_a-1) = temp8*b_a_m_h + 0.5*temp9*b_a_d + temp10*b_a_d_h;
  gsd_max.submat(0,num_a+num_d,num_a-1,num_a+num_d+num_e-1) = trans(gsd_max.submat(num_a+num_d,0,num_a+num_d+num_e-1,num_a-1));
  gsd_max.submat(num_a+num_d,num_a+num_d,num_a+num_d+num_e-1,num_a+num_d+num_e-1) = temp11*b_e_m + temp12*b_e_d;
  gsd_max.submat(num_a,num_a+num_d,num_a+num_d-1,num_a+num_d+num_e-1) = temp13*b_e_m_h + 0.75*temp17*b_e_d + temp14*b_e_d_h;
  gsd_max.submat(num_a+num_d, num_a, num_a+num_d+num_e-1,num_a+num_d-1) = trans(gsd_max.submat(num_a,num_a+num_d,num_a+num_d-1,num_a+num_d+num_e-1));
  
  gsd_max = 0.5*gsd_max;
  int num_t = num_a+num_d+num_e;
  arma::vec res((num_t+1)*num_t/2);
  int k = 0;
  for(int i=0; i<num_t; i++)
    for(int j=i; j<num_t; j++)
    {
      res[k] = gsd_max(i,j);
      k++;
    }
    
    return(Rcpp::wrap(res));
}
