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

RcppExport SEXP loglik_AtCtE_esp_c(SEXP var, SEXP b_c, SEXP b_a, SEXP pheno_m, SEXP pheno_d, SEXP B_a_m, SEXP B_a_d, SEXP B_c_m, SEXP B_c_d) {
    
    double e = as<double>(var); // Number of rows
    arma::vec bc = as<arma::vec>(b_c); // Number of columns
    arma::vec ba = as<arma::vec>(b_a);
	
    arma::vec p_m = as<arma::vec>(pheno_m);
    arma::vec p_d = as<arma::vec>(pheno_d);
	arma::mat b_a_m = as<arma::mat>(B_a_m);
    arma::mat b_a_d = as<arma::mat>(B_a_d);
	arma::mat b_c_m = as<arma::mat>(B_c_m);
    arma::mat b_c_d = as<arma::mat>(B_c_d);
	
	int num_m = p_m.n_elem;
	int num_d = p_d.n_elem;
	arma::mat k_m(2,2);
	arma::mat k_d(2,2);
	arma::mat k_d_a(2,2);
	k_d_a << 1 << 0.5 << arma::endr << 0.5 << 1 << arma::endr;
	arma::vec  v = arma::ones<arma::vec>(2);
	arma::mat diag = arma::diagmat(v);
	
	k_m.fill(1);
	k_d.fill(1);
	
	double YSY_m = 0;
	double YSY_d = 0;
	double D_m = 0;
	double D_d = 0;

	arma::mat V_m(2,2);
	arma::mat inv_V_m(2,2);

	for(int i = 0; i<(0.5*num_m); i++)
	{
		int start = 2*i;
		int end = start + 1;
		V_m = e*diag+exp(as_scalar(b_c_m.row(start)*bc))*k_m+exp(as_scalar(b_a_m.row(start)*ba))*k_m;
		double det_m = V_m(0,0)*V_m(0,0)-V_m(0,1)*V_m(0,1);
		//inv_V_m = inv(V_m);
		inv_V_m = V_m/det_m;
		inv_V_m(0,1) = (-1)*inv_V_m(0,1);
		inv_V_m(1,0) = inv_V_m(0,1);
		arma::vec p_m_v(2);
		arma::rowvec p_m_r(2);
		p_m_v << p_m[start] << p_m[end];
		p_m_r << p_m[start] << p_m[end];
		arma::mat temp = p_m_r*inv_V_m*p_m_v;
		YSY_m = YSY_m + as_scalar(temp);
		D_m = D_m + log(det_m);
	}

	for(int i = 0; i<(0.5*num_d); i++)
	{
		int start = 2*i;
		int end = start + 1;
		arma::mat V_d = e*diag+exp(as_scalar(b_c_d.row(start)*bc))*k_d+exp(as_scalar(b_a_d.row(start)*ba))*k_d_a;
		//arma::mat inv_V_d = inv(V_d);
		arma::mat inv_V_d(2,2);
		double det_d = V_d(0,0)*V_d(0,0)-V_d(0,1)*V_d(0,1);
		inv_V_d = V_d/det_d;
		inv_V_d(0,1) = (-1)*inv_V_d(0,1);
		inv_V_d(1,0) = inv_V_d(0,1);
		arma::vec p_d_v(2);
		arma::rowvec p_d_r(2);
		p_d_v << p_d[start] << p_d[end];
		p_d_r << p_d[start] << p_d[end];
		arma::mat temp = p_d_r*inv_V_d*p_d_v;
		YSY_d = YSY_d + as_scalar(temp);
		D_d = D_d + log(det_d);
	}
	
	double res = (D_m + YSY_m + D_d + YSY_d)/2;
   
    return(Rcpp::wrap(res));
}


RcppExport SEXP gr_AtCtE_esp_c(SEXP var, SEXP b_c, SEXP b_a, SEXP pheno_m, SEXP pheno_d, SEXP B_a_m, SEXP B_a_d, SEXP B_c_m, SEXP B_c_d) {
    
    double e = as<double>(var); // Number of rows
    arma::vec bc = as<arma::vec>(b_c); // Number of columns
    arma::vec ba = as<arma::vec>(b_a);
	
    arma::vec p_m = as<arma::vec>(pheno_m);
    arma::vec p_d = as<arma::vec>(pheno_d);
	arma::mat b_a_m = as<arma::mat>(B_a_m);
    arma::mat b_a_d = as<arma::mat>(B_a_d);
	arma::mat b_c_m = as<arma::mat>(B_c_m);
    arma::mat b_c_d = as<arma::mat>(B_c_d);
	
	int num_m = p_m.n_elem;
	int num_d = p_d.n_elem;
	int num_a = ba.n_elem;
	int num_c = bc.n_elem;
	arma::mat k_m(2,2);
	arma::mat k_d(2,2);
	arma::mat k_d_a(2,2);
	k_d_a << 1 << 0.5 << arma::endr << 0.5 << 1 << arma::endr;
	arma::vec  v = arma::ones<arma::vec>(2);
	arma::mat diag = arma::diagmat(v);
	
	k_m.fill(1);
	k_d.fill(1);
	
	double var_YSY_m = 0;
	double var_YSY_d = 0;
	double var_D_m = 0;
	double var_D_d = 0;
	arma::vec var_bc_m = arma::zeros<arma::vec>(num_c);
	arma::vec var_ba_m = arma::zeros<arma::vec>(num_a);

	arma::mat V_m(2,2);
	arma::mat inv_V_m(2,2);

	for(int i = 0; i<(0.5*num_m); i++)
	{
		int start = 2*i;
		int end = start + 1;
		V_m = e*diag+exp(as_scalar(b_c_m.row(start)*bc))*k_m+exp(as_scalar(b_a_m.row(start)*ba))*k_m;
		inv_V_m = inv(V_m);
		arma::vec p_m_v(2);
		arma::rowvec p_m_r(2);
		p_m_v << p_m[start] << p_m[end];
		p_m_r << p_m[start] << p_m[end];
		arma::mat inv_ph_m = inv_V_m*p_m_v;
		arma::mat ph_inv_m = p_m_r*inv_V_m;
		var_YSY_m = var_YSY_m + as_scalar(ph_inv_m*inv_ph_m);
		var_D_m = var_D_m + arma::sum(inv_V_m.diag());
		arma::mat temp_prod = inv_V_m*k_m;
		double temp_c = exp(as_scalar(b_c_m.row(start)*bc))*(arma::sum(temp_prod.diag())-as_scalar(ph_inv_m*k_m*inv_ph_m));
		for(int j = 0; j < num_c; j++)
		{
			var_bc_m(j) = var_bc_m(j) + b_c_m(start,j)*temp_c;
		}
		double temp_a = exp(as_scalar(b_a_m.row(start)*ba))*(arma::sum(temp_prod.diag())-as_scalar(ph_inv_m*k_m*inv_ph_m));
		for(int j = 0; j < num_a; j++)
		{
			var_ba_m(j) = var_ba_m(j) + b_a_m(start,j)*temp_a;
		}
		
	}

	arma::vec var_bc_d = arma::zeros<arma::vec>(num_c);
	arma::vec var_ba_d = arma::zeros<arma::vec>(num_a);

	for(int i = 0; i<(0.5*num_d); i++)
	{
		int start = 2*i;
		int end = start + 1;
		arma::mat V_d = e*diag+exp(as_scalar(b_c_d.row(start)*bc))*k_d+exp(as_scalar(b_a_d.row(start)*ba))*k_d_a;
		arma::mat inv_V_d = inv(V_d);
		arma::vec p_d_v(2);
		arma::rowvec p_d_r(2);
		p_d_v << p_d[start] << p_d[end];
		p_d_r << p_d[start] << p_d[end];
		arma::mat inv_ph_d = inv_V_d*p_d_v;
		arma::mat ph_inv_d = p_d_r*inv_V_d;
		var_YSY_d = var_YSY_d + as_scalar(ph_inv_d*inv_ph_d);
		var_D_d = var_D_d + arma::sum(inv_V_d.diag());
		arma::mat temp_prod = inv_V_d*k_d;
		double temp_c = exp(as_scalar(b_c_d.row(start)*bc))*(arma::sum(temp_prod.diag())-as_scalar(ph_inv_d*k_d*inv_ph_d));
		for(int j = 0; j < num_c; j++)
		{
			var_bc_d(j) = var_bc_d(j) + b_c_d(start,j)*temp_c;
		}
		temp_prod = inv_V_d*k_d_a;
		double temp_a = exp(as_scalar(b_a_d.row(start)*ba))*(arma::sum(temp_prod.diag())-as_scalar(ph_inv_d*k_d_a*inv_ph_d));
		for(int j = 0; j < num_a; j++)
		{
			var_ba_d(j) = var_ba_d(j) + b_a_d(start,j)*temp_a;
		}
	}

	double d_var = (var_D_m - var_YSY_m + var_D_d - var_YSY_d)/2;
	arma::vec d_var_c = (var_bc_m + var_bc_d)/2;
	arma::vec d_var_a = (var_ba_m + var_ba_d)/2;
	
	arma::vec res(1+num_c+num_a);
	res(0) = d_var;
	for(int i = 0; i < num_c; i++)
	{
		res(1+i) = d_var_c(i);
	}
	for(int i = 0; i < num_a; i++)
	{
		res(1+num_c+i) = d_var_a(i);
	}
   
    return(Rcpp::wrap(res));
}
