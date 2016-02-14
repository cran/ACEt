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

RcppExport SEXP loglik_AtCtE_epsp_c(SEXP var, SEXP v_b_a, SEXP v_b_c, SEXP pheno_m, SEXP pheno_d, SEXP B_a_m, SEXP B_a_d, SEXP b_a, SEXP D_a, SEXP B_c_m, SEXP B_c_d, SEXP b_c, SEXP D_c) {
    
    double e = as<double>(var); // Number of rows
    double va = as<double>(v_b_a); // Number of columns
    double vc = as<double>(v_b_c);
	
    arma::vec p_m = as<arma::vec>(pheno_m);
    arma::vec p_d = as<arma::vec>(pheno_d);
	arma::mat b_a_m = as<arma::mat>(B_a_m);
    arma::mat b_a_d = as<arma::mat>(B_a_d);
	arma::mat b_c_m = as<arma::mat>(B_c_m);
    arma::mat b_c_d = as<arma::mat>(B_c_d);
	arma::mat d_a = as<arma::mat>(D_a);
	arma::mat d_c = as<arma::mat>(D_c);
	arma::vec ba = as<arma::vec>(b_a);
	arma::vec bc = as<arma::vec>(b_c);
	
	int num_m = p_m.n_elem;
	int num_d = p_d.n_elem;
	int num_a = b_a_m.n_cols;
	int num_c = b_c_m.n_cols;
	arma::mat k_m(2,2);
	arma::mat k_1d(2,2);
	arma::mat k_2d(2,2);
	arma::mat k_3d(2,2);
	k_1d << 0.5 << 0 << arma::endr << 0 << 0.5 << arma::endr;
	arma::vec  v = arma::ones<arma::vec>(2);
	arma::mat diag = arma::diagmat(v);
	
	k_m.fill(1);
	k_2d.fill(0.5);
	k_3d.fill(1);
	
	double YSY_m = 0;
	double YSY_d = 0;
	double D_m = 0;
	double D_d = 0;
	arma::rowvec  gsd_a2m2a = arma::zeros<arma::rowvec>(num_m/2);
	arma::rowvec  gsd_c2m2a = arma::zeros<arma::rowvec>(num_m/2);
	arma::rowvec  gsd_c2m2c = arma::zeros<arma::rowvec>(num_m/2);
	arma::rowvec  gsd_dd = arma::zeros<arma::rowvec>(num_d);
	arma::rowvec  gsd_2dd2 = arma::zeros<arma::rowvec>(num_d/2);
	arma::rowvec  gsd_3dd3 = arma::zeros<arma::rowvec>(num_d/2);
	arma::rowvec  gsd_2d2 = arma::zeros<arma::rowvec>(num_d/2);
	arma::rowvec  gsd_3d2 = arma::zeros<arma::rowvec>(num_d/2);
	arma::rowvec  gsd_3d3 = arma::zeros<arma::rowvec>(num_d/2);

	for(int i = 0; i<(0.5*num_m); i++)
	{
		int start = 2*i;
		int end = start + 1;
		double ba_m = exp(as_scalar(b_a_m.row(start)*ba));
		double bc_m = exp(as_scalar(b_c_m.row(start)*bc));
		arma::mat V_m = e*diag + ba_m*k_m + bc_m*k_m;
		arma::mat inv_V_m = inv(V_m);
		arma::vec p_m_v(2);
		arma::rowvec p_m_r(2);
		p_m_v << p_m[start] << p_m[end];
		p_m_r << p_m[start] << p_m[end];
		arma::mat temp = p_m_r*inv_V_m*p_m_v;
		YSY_m = YSY_m + as_scalar(temp);
		D_m = D_m + log(V_m(0,0)*V_m(0,0)-V_m(0,1)*V_m(0,1));
		double accu_temp = arma::accu(inv_V_m*k_m*inv_V_m);
		gsd_a2m2a(i) = ba_m*ba_m*accu_temp;
		gsd_c2m2a(i) = ba_m*accu_temp*bc_m;
		gsd_c2m2c(i) = bc_m*bc_m*accu_temp;
	}

	for(int i = 0; i<(0.5*num_d); i++)
	{
		int start = 2*i;
		int end = start + 1;
		double ba_d = exp(as_scalar(b_a_d.row(start)*ba));
		double bc_d = exp(as_scalar(b_c_d.row(start)*bc));
		double ba2_d = ba_d*ba_d;
		double bac_d = ba_d*bc_d;
		double bc2_d = bc_d*bc_d;
		arma::mat V_d = e*diag + ba_d*(k_1d+k_2d) + bc_d*k_3d;
		arma::mat inv_V_d = inv(V_d);
		arma::vec p_d_v(2);
		arma::rowvec p_d_r(2);
		p_d_v << p_d[start] << p_d[end];
		p_d_r << p_d[start] << p_d[end];
		arma::mat temp = p_d_r*inv_V_d*p_d_v;
		YSY_d = YSY_d + as_scalar(temp);
		D_d = D_d + log(V_d(0,0)*V_d(0,0)-V_d(0,1)*V_d(0,1));
		arma::mat temp_dd = inv_V_d*inv_V_d;
		gsd_dd(start) = ba2_d*temp_dd(0,0);
		gsd_dd(end) = ba2_d*temp_dd(1,1);
		gsd_2dd2(i) = ba2_d*0.5*arma::accu(temp_dd);
		gsd_3dd3(i) = bac_d*arma::accu(temp_dd);
		double temp_si2i = arma::accu(inv_V_d*k_2d*inv_V_d);
		gsd_2d2(i) = ba2_d*0.5*temp_si2i;
		gsd_3d2(i) = bac_d*temp_si2i;
		gsd_3d3(i) = bc2_d*2*temp_si2i;
	}

	arma::mat gsd_max(num_a+num_c,num_a+num_c);
	arma::mat b_a_m_h(num_m/2,num_a);
	arma::mat b_c_m_h(num_m/2,num_c);
	arma::mat b_a_d_h(num_d/2,num_a);
	arma::mat b_c_d_h(num_d/2,num_c);
	for(int i = 0; i < (0.5*num_m); i++)
	{
		b_a_m_h.row(i) = b_a_m.row(2*i);
		b_c_m_h.row(i) = b_c_m.row(2*i);
	}
	for(int i = 0; i < (0.5*num_d); i++)
	{
		b_a_d_h.row(i) = b_a_d.row(2*i);
		b_c_d_h.row(i) = b_c_d.row(2*i);
	}
	//arma::mat b_a_m_t = b_a_m.t();
	//arma::mat b_c_m_t = b_c_m.t();
	arma::mat b_a_d_t = b_a_d.t();
	//arma::mat b_c_d_t = b_c_d.t();
	arma::mat b_a_m_ht = b_a_m_h.t();
	arma::mat b_c_m_ht = b_c_m_h.t();
	arma::mat b_a_d_ht = b_a_d_h.t();
	arma::mat b_c_d_ht = b_c_d_h.t();
	arma::mat temp1 = b_a_m_ht;
	arma::mat temp4 = b_c_m_ht;
	arma::mat temp6 = b_c_m_ht;
	for(int i = 0; i < (0.5*num_m); i++)
	{
		temp1.col(i) = gsd_a2m2a(i)*b_a_m_ht.col(i);
		temp4.col(i) = gsd_c2m2a(i)*b_c_m_ht.col(i);
		temp6.col(i) = gsd_c2m2c(i)*b_c_m_ht.col(i);
	}
	arma::mat temp3 = b_a_d_ht;
	arma::rowvec tempv1 = gsd_2dd2+gsd_2d2;
	arma::mat temp5 = b_c_d_ht;
	arma::rowvec tempv2 = 0.5*gsd_3dd3+gsd_3d2;
	arma::mat temp7 = b_c_d_ht;
	for(int i = 0; i < (0.5*num_d); i++)
	{
		temp3.col(i) = tempv1(i)*b_a_d_ht.col(i);
		temp5.col(i) = tempv2(i)*b_c_d_ht.col(i);
		temp7.col(i) = gsd_3d3(i)*b_c_d_ht.col(i);
	}
	arma::mat temp2 = b_a_d_t;
	for(int i = 0; i < num_d; i++)
	{
		temp2.col(i) = gsd_dd(i)*b_a_d_t.col(i);
	}
	gsd_max.submat(0,0,num_a-1,num_a-1) = arma::diagmat(2*d_a.diag())/va + temp1*b_a_m_h + 0.25*temp2*b_a_d + temp3*b_a_d_h;
	gsd_max.submat(num_a,0,num_a+num_c-1,num_a-1) = temp4*b_a_m_h + temp5*b_a_d_h;
	gsd_max.submat(0,num_a,num_a-1,num_a+num_c-1) = trans(gsd_max.submat(num_a,0,num_a+num_c-1,num_a-1));
	gsd_max.submat(num_a,num_a,num_a+num_c-1,num_a+num_c-1) = arma::diagmat(2*d_c.diag())/vc + temp6*b_c_m_h + temp7*b_c_d_h;

	double gsd = log(fabs(det(0.5*gsd_max)));

	double res = D_m + YSY_m + D_d + YSY_d + (num_a-2)*log(va) + (num_c-2)*log(vc) + as_scalar(trans(ba)*d_a*ba)/va + as_scalar(trans(bc)*d_c*bc)/vc + gsd;

	return(Rcpp::wrap(res));
}

RcppExport SEXP loglik_AtCtE_epsp_g_c(SEXP b_a, SEXP b_c, SEXP pheno_m, SEXP pheno_d, SEXP B_a_m, SEXP B_a_d, SEXP B_c_m, SEXP B_c_d, SEXP var, SEXP v_b_a, SEXP v_b_c, SEXP D_a, SEXP D_c) {

    arma::vec bc = as<arma::vec>(b_c); // Number of columns
    arma::vec ba = as<arma::vec>(b_a);
	
    arma::vec p_m = as<arma::vec>(pheno_m);
    arma::vec p_d = as<arma::vec>(pheno_d);
	arma::mat b_a_m = as<arma::mat>(B_a_m);
    arma::mat b_a_d = as<arma::mat>(B_a_d);
	arma::mat b_c_m = as<arma::mat>(B_c_m);
    arma::mat b_c_d = as<arma::mat>(B_c_d);

	double e = as<double>(var);
	double va = as<double>(v_b_a);
	double vc = as<double>(v_b_c);
	arma::mat d_a = as<arma::mat>(D_a);
	arma::mat d_c = as<arma::mat>(D_c);
	
	int num_m = p_m.n_elem;
	int num_d = p_d.n_elem;
	int num_a = b_a_m.n_cols;
	int num_c = b_c_m.n_cols;
	arma::mat k_m(2,2);
	arma::mat k_1d(2,2);
	arma::mat k_2d(2,2);
	arma::mat k_3d(2,2);
	k_1d << 0.5 << 0 << arma::endr << 0 << 0.5 << arma::endr;
	arma::vec  v = arma::ones<arma::vec>(2);
	arma::mat diag = arma::diagmat(v);
	
	k_m.fill(1);
	k_2d.fill(0.5);
	k_3d.fill(1);
	
	double YSY_m = 0;
	double YSY_d = 0;
	double D_m = 0;
	double D_d = 0;

	for(int i = 0; i < (0.5*num_m); i++)
	{
		int start = 2*i;
		int end = start + 1;
		arma::mat V_m = e*diag+exp(as_scalar(b_c_m.row(start)*bc))*k_m+exp(as_scalar(b_a_m.row(start)*ba))*k_m;
		double det_m = V_m(0,0)*V_m(0,0)-V_m(0,1)*V_m(0,1);
		//inv_V_m = inv(V_m);
		arma::mat inv_V_m = V_m/det_m;
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

	for(int i = 0; i < (0.5*num_d); i++)
	{
		int start = 2*i;
		int end = start + 1;
		arma::mat V_d = e*diag+exp(as_scalar(b_c_d.row(start)*bc))*k_3d+exp(as_scalar(b_a_d.row(start)*ba))*(k_1d+k_2d);
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

	double res = D_m + YSY_m + D_d + YSY_d + as_scalar(trans(ba)*d_a*ba)/va + as_scalar(trans(bc)*d_c*bc)/vc;

	return(Rcpp::wrap(res));
}

RcppExport SEXP gr_AtCtE_epsp_g_c(SEXP b_a, SEXP b_c, SEXP pheno_m, SEXP pheno_d, SEXP B_a_m, SEXP B_a_d, SEXP B_c_m, SEXP B_c_d, SEXP var, SEXP v_b_a, SEXP v_b_c, SEXP D_a, SEXP D_c)
{
	arma::vec bc = as<arma::vec>(b_c); // Number of columns
    arma::vec ba = as<arma::vec>(b_a);
	
    arma::vec p_m = as<arma::vec>(pheno_m);
    arma::vec p_d = as<arma::vec>(pheno_d);
	arma::mat b_a_m = as<arma::mat>(B_a_m);
    arma::mat b_a_d = as<arma::mat>(B_a_d);
	arma::mat b_c_m = as<arma::mat>(B_c_m);
    arma::mat b_c_d = as<arma::mat>(B_c_d);

	double e = as<double>(var);
	double va = as<double>(v_b_a);
	double vc = as<double>(v_b_c);
	arma::mat d_a = as<arma::mat>(D_a);
	arma::mat d_c = as<arma::mat>(D_c);
	
	int num_m = p_m.n_elem;
	int num_d = p_d.n_elem;
	int num_a = b_a_m.n_cols;
	int num_c = b_c_m.n_cols;
	arma::mat k_m(2,2);
	arma::mat k_1d(2,2);
	arma::mat k_2d(2,2);
	arma::mat k_3d(2,2);
	k_1d << 0.5 << 0 << arma::endr << 0 << 0.5 << arma::endr;
	arma::vec  v = arma::ones<arma::vec>(2);
	arma::mat diag = arma::diagmat(v);
	
	k_m.fill(1);
	k_2d.fill(0.5);
	k_3d.fill(1);
	
	arma::vec v_b_a_m = arma::ones<arma::vec>(num_a);
	arma::vec v_b_c_m = arma::ones<arma::vec>(num_c);
	
	for(int i = 0; i<(0.5*num_m); i++)
	{
		int start = 2*i;
		int end = start + 1;
		double ba_m = exp(as_scalar(b_a_m.row(start)*ba));
		double bc_m = exp(as_scalar(b_c_m.row(start)*bc));
		arma::mat V_m = e*diag + ba_m*k_m + bc_m*k_m;
		arma::mat inv_V_m = inv(V_m);
		arma::vec p_m_v(2);
		arma::rowvec p_m_r(2);
		p_m_v << p_m[start] << p_m[end];
		p_m_r << p_m[start] << p_m[end];
		arma::mat inv_ph_m = inv_V_m*p_m_v;
		arma::mat ph_inv_m = p_m_r*inv_V_m;
		arma::mat temp_prod = inv_V_m*k_m;
		double temp_a = ba_m*(arma::sum(temp_prod.diag())-as_scalar(ph_inv_m*k_m*inv_ph_m));
		for(int j = 0; j < num_a; j++)
		{
			v_b_a_m(j) = v_b_a_m(j) + b_a_m(start,j)*temp_a;
		}
		double temp_c = bc_m*(arma::sum(temp_prod.diag())-as_scalar(ph_inv_m*k_m*inv_ph_m));
		for(int j = 0; j < num_c; j++)
		{
			v_b_c_m(j) = v_b_c_m(j) + b_c_m(start,j)*temp_c;
		}
	}

	arma::vec v_b_a_d = arma::ones<arma::vec>(num_a);
	arma::vec v_b_c_d = arma::ones<arma::vec>(num_c);
	
	for(int i = 0; i<(0.5*num_d); i++)
	{
		int start = 2*i;
		int end = start + 1;
		double ba_d = exp(as_scalar(b_a_d.row(start)*ba));
		double bc_d = exp(as_scalar(b_c_d.row(start)*bc));
		arma::mat V_d = e*diag + ba_d*(k_1d+k_2d) + bc_d*k_3d;
		arma::mat inv_V_d = inv(V_d);
		arma::vec p_d_v(2);
		arma::rowvec p_d_r(2);
		p_d_v << p_d[start] << p_d[end];
		p_d_r << p_d[start] << p_d[end];
		arma::mat inv_ph_d = inv_V_d*p_d_v;
		arma::mat ph_inv_d = p_d_r*inv_V_d;
		arma::mat temp_prod = inv_V_d*(k_1d+k_2d);
		double temp_a = ba_d*(arma::sum(temp_prod.diag())-as_scalar(ph_inv_d*(k_1d+k_2d)*inv_ph_d));
		for(int j = 0; j < num_a; j++)
		{
			v_b_a_d(j) = v_b_a_d(j) + b_a_d(start,j)*temp_a;
		}
		temp_prod = inv_V_d*k_3d;
		double temp_c = bc_d*(arma::sum(temp_prod.diag())-as_scalar(ph_inv_d*k_3d*inv_ph_d));
		for(int j = 0; j < num_c; j++)
		{
			v_b_c_d(j) = v_b_c_d(j) + b_c_d(start,j)*temp_c;
		}
	}

	arma::vec res(num_a+num_c);
	arma::vec res_a = v_b_a_m + v_b_a_d + 2*d_a*ba/va;
	arma::vec res_c = v_b_c_m + v_b_c_d + 2*d_c*bc/vc;
	for(int i = 0; i < num_a; i++)
	{
		res(i) = res_a(i);
	}
	for(int i = 0; i < num_c; i++)
	{
		res(i+num_a) = res_c(i);
	}
	return(Rcpp::wrap(res));
}

RcppExport SEXP gr_AtCtE_epsp_c(SEXP var, SEXP v_b_a, SEXP v_b_c, SEXP pheno_m, SEXP pheno_d, SEXP B_a_m, SEXP B_a_d, SEXP b_a, SEXP D_a, SEXP B_c_m, SEXP B_c_d, SEXP b_c, SEXP D_c)
{
	 double e = as<double>(var); // Number of rows
    double va = as<double>(v_b_a); // Number of columns
    double vc = as<double>(v_b_c);
	
    arma::vec p_m = as<arma::vec>(pheno_m);
    arma::vec p_d = as<arma::vec>(pheno_d);
	arma::mat b_a_m = as<arma::mat>(B_a_m);
    arma::mat b_a_d = as<arma::mat>(B_a_d);
	arma::mat b_c_m = as<arma::mat>(B_c_m);
    arma::mat b_c_d = as<arma::mat>(B_c_d);
	arma::mat d_a = as<arma::mat>(D_a);
	arma::mat d_c = as<arma::mat>(D_c);
	arma::vec ba = as<arma::vec>(b_a);
	arma::vec bc = as<arma::vec>(b_c);
	
	int num_m = p_m.n_elem;
	int num_d = p_d.n_elem;
	int num_a = b_a_m.n_cols;
	int num_c = b_c_m.n_cols;
	arma::mat k_m(2,2);
	arma::mat k_1d(2,2);
	arma::mat k_2d(2,2);
	arma::mat k_3d(2,2);
	k_1d << 0.5 << 0 << arma::endr << 0 << 0.5 << arma::endr;
	arma::vec  v = arma::ones<arma::vec>(2);
	arma::mat diag = arma::diagmat(v);
	
	k_m.fill(1);
	k_2d.fill(0.5);
	k_3d.fill(1);
	
	double YSY_m = 0;
	double YSY_d = 0;
	double D_m = 0;
	double D_d = 0;
	arma::rowvec  gsd_a2m2a = arma::zeros<arma::rowvec>(num_m/2);
	arma::rowvec  gsd_c2m2a = arma::zeros<arma::rowvec>(num_m/2);
	arma::rowvec  gsd_c2m2c = arma::zeros<arma::rowvec>(num_m/2);
	arma::rowvec  gsd_amm2ma = arma::zeros<arma::rowvec>(num_m/2);
	arma::rowvec  gsd_am2mma = arma::zeros<arma::rowvec>(num_m/2);
	arma::rowvec  gsd_cmm2ma = arma::zeros<arma::rowvec>(num_m/2);
	arma::rowvec  gsd_cm2mma = arma::zeros<arma::rowvec>(num_m/2);
	arma::rowvec  gsd_cmm2mc = arma::zeros<arma::rowvec>(num_m/2);
	arma::rowvec  gsd_cm2mmc = arma::zeros<arma::rowvec>(num_m/2);
	arma::rowvec  gsd_dd = arma::zeros<arma::rowvec>(num_d);
	arma::rowvec  gsd_2dd2 = arma::zeros<arma::rowvec>(num_d/2);
	arma::rowvec  gsd_3dd3 = arma::zeros<arma::rowvec>(num_d/2);
	arma::rowvec  gsd_2d2 = arma::zeros<arma::rowvec>(num_d/2);
	arma::rowvec  gsd_3d2 = arma::zeros<arma::rowvec>(num_d/2);
	arma::rowvec  gsd_3d3 = arma::zeros<arma::rowvec>(num_d/2);
	arma::rowvec  gsd_ddd = arma::zeros<arma::rowvec>(num_d);
	arma::rowvec  gsd_2ddd2 = arma::zeros<arma::rowvec>(num_d/2);
	arma::rowvec  gsd_2dd2d2 = arma::zeros<arma::rowvec>(num_d/2);
	arma::rowvec  gsd_2d2dd2 = arma::zeros<arma::rowvec>(num_d/2);
	arma::rowvec  gsd_3ddd3 = arma::zeros<arma::rowvec>(num_d/2);
	arma::rowvec  gsd_3dd2d3 = arma::zeros<arma::rowvec>(num_d/2);
	arma::rowvec  gsd_3d2dd3 = arma::zeros<arma::rowvec>(num_d/2);
	arma::rowvec  gsd_3dd3d3 = arma::zeros<arma::rowvec>(num_d/2);
	arma::rowvec  gsd_3d3dd3 = arma::zeros<arma::rowvec>(num_d/2);

	for(int i = 0; i<(0.5*num_m); i++)
	{
		int start = 2*i;
		int end = start + 1;
		double ba_m = exp(as_scalar(b_a_m.row(start)*ba));
		double bc_m = exp(as_scalar(b_c_m.row(start)*bc));
		double ba2_m = ba_m*ba_m;
		double bc2_m = bc_m*bc_m;
		double bac_m = ba_m*bc_m;
		arma::mat V_m = e*diag + ba_m*k_m + bc_m*k_m;
		arma::mat inv_V_m = inv(V_m);
		arma::vec p_m_v(2);
		arma::rowvec p_m_r(2);
		p_m_v << p_m[start] << p_m[end];
		p_m_r << p_m[start] << p_m[end];
		arma::mat inv_ph_m = inv_V_m*p_m_v;
		arma::mat ph_inv_m = p_m_r*inv_V_m;
		YSY_m = YSY_m + as_scalar(ph_inv_m*inv_ph_m);
		D_m = D_m + arma::sum(inv_V_m.diag());
		arma::mat temp_prod = inv_V_m*k_m*inv_V_m;
		gsd_a2m2a(i) = ba2_m*arma::accu(temp_prod);
		gsd_c2m2a(i) = bac_m*arma::accu(temp_prod);
		gsd_c2m2c(i) = bc2_m*arma::accu(temp_prod);
		double temp_mm2m = arma::accu(inv_V_m*temp_prod);
		double temp_m2mm = arma::accu(temp_prod*inv_V_m);
		gsd_amm2ma(i) = ba2_m*temp_mm2m;
		gsd_am2mma(i) = ba2_m*temp_m2mm;
		gsd_cmm2ma(i) = bac_m*temp_mm2m;
		gsd_cm2mma(i) = bac_m*temp_m2mm;
		gsd_cmm2mc(i) = bc2_m*temp_mm2m;
		gsd_cm2mmc(i) = bc2_m*temp_m2mm;
	}

	for(int i = 0; i < (0.5*num_d); i++)
	{
		int start = 2*i;
		int end = start + 1;
		double ba_d = exp(as_scalar(b_a_d.row(start)*ba));
		double bc_d = exp(as_scalar(b_c_d.row(start)*bc));
		double ba2_d = ba_d*ba_d;
		double bc2_d = bc_d*bc_d;
		double bac_d = ba_d*bc_d;
		arma::mat V_d = e*diag + ba_d*(k_1d+k_2d) + bc_d*k_3d;
		arma::mat inv_V_d = inv(V_d);
		arma::vec p_d_v(2);
		arma::rowvec p_d_r(2);
		p_d_v << p_d[start] << p_d[end];
		p_d_r << p_d[start] << p_d[end];
		arma::mat inv_ph_d = inv_V_d*p_d_v;
		arma::mat ph_inv_d = p_d_r*inv_V_d;
		YSY_d = YSY_d + as_scalar(ph_inv_d*inv_ph_d);
		D_d = D_d + arma::sum(inv_V_d.diag());
		arma::mat temp_dd = inv_V_d*inv_V_d;
		gsd_dd(start) = ba2_d*temp_dd(0,0);
		gsd_dd(end) = ba2_d*temp_dd(1,1);
		gsd_2dd2(i) = ba2_d*0.5*arma::accu(temp_dd);
		gsd_3dd3(i) = bac_d*arma::accu(temp_dd);
		arma::mat temp_i2i = inv_V_d*k_2d*inv_V_d;
		gsd_2d2(i) = ba2_d*0.5*arma::accu(temp_i2i);
		gsd_3d2(i) = bac_d*arma::accu(temp_i2i);
		arma::mat temp_i3i = 2*temp_i2i;
		gsd_3d3(i) = bc2_d*arma::accu(temp_i3i);
		arma::mat temp_ddd = temp_dd*inv_V_d;
		gsd_ddd(start) = ba2_d*temp_ddd(0,0);
		gsd_ddd(end) = ba2_d*temp_ddd(1,1);
		gsd_2ddd2(i) = ba2_d*0.5*arma::accu(temp_ddd);
		gsd_3ddd3(i) = bac_d*arma::accu(temp_ddd);
		arma::mat temp_dd2d = inv_V_d*temp_i2i;
		arma::mat temp_d2dd = temp_i2i*inv_V_d;
		arma::mat temp_dd3d = inv_V_d*temp_i3i;
		arma::mat temp_d3dd = temp_i3i*inv_V_d;
		gsd_2dd2d2(i) = ba2_d*0.5*arma::accu(temp_dd2d);
		gsd_2d2dd2(i) = ba2_d*0.5*arma::accu(temp_d2dd);
		gsd_3dd2d3(i) = bac_d*arma::accu(temp_dd2d);
		gsd_3d2dd3(i) = bac_d*arma::accu(temp_d2dd);
		gsd_3dd3d3(i) = bc2_d*arma::accu(temp_dd3d);
		gsd_3d3dd3(i) = bc2_d*arma::accu(temp_d3dd);
	}

	arma::mat gsd_max(num_a+num_c,num_a+num_c);
	arma::mat gsd_b_a(num_a+num_c,num_a+num_c);
	arma::mat gsd_b_c(num_a+num_c,num_a+num_c);
	arma::mat gsd_e(num_a+num_c,num_a+num_c);
	gsd_b_a.fill(0);
	gsd_b_c.fill(0);
	gsd_e.fill(0);

	arma::mat b_a_m_h(num_m/2,num_a);
	arma::mat b_c_m_h(num_m/2,num_c);
	arma::mat b_a_d_h(num_d/2,num_a);
	arma::mat b_c_d_h(num_d/2,num_c);
	for(int i = 0; i < (0.5*num_m); i++)
	{
		b_a_m_h.row(i) = b_a_m.row(2*i);
		b_c_m_h.row(i) = b_c_m.row(2*i);
	}
	for(int i = 0; i < (0.5*num_d); i++)
	{
		b_a_d_h.row(i) = b_a_d.row(2*i);
		b_c_d_h.row(i) = b_c_d.row(2*i);
	}
	//arma::mat b_a_m_t = b_a_m.t();
	//arma::mat b_c_m_t = b_c_m.t();
	arma::mat b_a_d_t = b_a_d.t();
	//arma::mat b_c_d_t = b_c_d.t();
	arma::mat b_a_m_ht = b_a_m_h.t();
	arma::mat b_c_m_ht = b_c_m_h.t();
	arma::mat b_a_d_ht = b_a_d_h.t();
	arma::mat b_c_d_ht = b_c_d_h.t();
	arma::mat temp1 = b_a_m_ht;
	arma::mat temp4 = b_c_m_ht;
	arma::mat temp6 = b_c_m_ht;
	for(int i = 0; i < (0.5*num_m); i++)
	{
		temp1.col(i) = gsd_a2m2a(i)*b_a_m_ht.col(i);
		temp4.col(i) = gsd_c2m2a(i)*b_c_m_ht.col(i);
		temp6.col(i) = gsd_c2m2c(i)*b_c_m_ht.col(i);
	}
	arma::mat temp3 = b_a_d_ht;
	arma::rowvec tempv1 = gsd_2dd2+gsd_2d2;
	arma::mat temp5 = b_c_d_ht;
	arma::rowvec tempv2 = 0.5*gsd_3dd3+gsd_3d2;
	arma::mat temp7 = b_c_d_ht;
	for(int i = 0; i < (0.5*num_d); i++)
	{
		temp3.col(i) = tempv1(i)*b_a_d_ht.col(i);
		temp5.col(i) = tempv2(i)*b_c_d_ht.col(i);
		temp7.col(i) = gsd_3d3(i)*b_c_d_ht.col(i);
	}
	arma::mat temp2 = b_a_d_t;
	for(int i = 0; i < num_d; i++)
	{
		temp2.col(i) = gsd_dd(i)*b_a_d_t.col(i);
	}

	gsd_max.submat(0,0,num_a-1,num_a-1) = arma::diagmat(2*d_a.diag())/va + temp1*b_a_m_h + 0.25*temp2*b_a_d + temp3*b_a_d_h;
	gsd_max.submat(num_a,0,num_a+num_c-1,num_a-1) = temp4*b_a_m_h + temp5*b_a_d_h;
	gsd_max.submat(0,num_a,num_a-1,num_a+num_c-1) = trans(gsd_max.submat(num_a,0,num_a+num_c-1,num_a-1));
	gsd_max.submat(num_a,num_a,num_a+num_c-1,num_a+num_c-1) = arma::diagmat(2*d_c.diag())/vc + temp6*b_c_m_h + temp7*b_c_d_h;
	arma::mat inv_gsd_max = arma::inv(gsd_max);

	gsd_b_a.submat(0,0,num_a-1,num_a-1) = arma::diagmat(2*d_a.diag());
	gsd_b_c.submat(num_a,num_a,num_a+num_c-1,num_a+num_c-1) = arma::diagmat(2*d_c.diag());

	temp1 = b_a_m_ht;
	temp2 = b_a_d_t;
	temp3 = b_a_d_ht;
	temp4 = b_c_m_ht;
	temp5 = b_c_d_ht;
	temp6 = b_c_m_ht;
	temp7 = b_c_d_ht;
	tempv1 = gsd_amm2ma+gsd_am2mma;
	arma::rowvec tempv3 = 2*gsd_2ddd2+gsd_2dd2d2+gsd_2d2dd2;
	arma::rowvec tempv4 = gsd_cmm2ma+gsd_cm2mma;
	arma::rowvec tempv5 = gsd_3ddd3+gsd_3dd2d3+gsd_3d2dd3;
	arma::rowvec tempv6 = gsd_cmm2mc+gsd_cm2mmc;
	arma::rowvec tempv7 = gsd_3dd3d3+gsd_3d3dd3;

	for(int i = 0; i < (0.5*num_m); i++)
	{
		temp1.col(i) = tempv1(i)*b_a_m_ht.col(i);
		temp4.col(i) = tempv4(i)*b_c_m_ht.col(i);
		temp6.col(i) = tempv6(i)*b_c_m_ht.col(i);
	}
	for(int i = 0; i < (0.5*num_d); i++)
	{
		temp3.col(i) = tempv3(i)*b_a_d_ht.col(i);
		temp5.col(i) = tempv5(i)*b_c_d_ht.col(i);
		temp7.col(i) = tempv7(i)*b_c_d_ht.col(i);
	}
	for(int i = 0; i < num_d; i++)
	{
		temp2.col(i) = gsd_ddd(i)*b_a_d_t.col(i);
	}

	gsd_e.submat(0,0,num_a-1,num_a-1) = temp1*b_a_m_h + 0.5*temp2*b_a_d + temp3*b_a_d_h;
	gsd_e.submat(num_a,0,num_a+num_c-1,num_a-1) = temp4*b_a_m_h + temp5*b_a_d_h;
	gsd_e.submat(0,num_a,num_a-1,num_a+num_c-1) = trans(gsd_e.submat(num_a,0,num_a+num_c-1,num_a-1));
	gsd_e.submat(num_a,num_a,num_a+num_c-1,num_a+num_c-1) = temp6*b_c_m_h + temp7*b_c_d_h;

	gsd_e = (-0.5)*gsd_e;

	arma::mat temp_m = inv_gsd_max*gsd_b_a;
	double d_v_b_a = ((-1)*as_scalar(trans(ba)*d_a*ba)-0.5*arma::sum(temp_m.diag()))/(va*va) + (num_a-2)/va;
	temp_m = inv_gsd_max*gsd_b_c;
	double d_v_b_c = ((-1)*as_scalar(trans(bc)*d_c*bc)-0.5*arma::sum(temp_m.diag()))/(vc*vc) + (num_c-2)/vc;
	temp_m = inv_gsd_max*gsd_e;
	double d_v = D_m - YSY_m + D_d - YSY_d +arma::sum(temp_m.diag());

	arma::vec res(3);
	res(0) = d_v;
	res(1) = d_v_b_a;
	res(2) = d_v_b_c;
   
    return(Rcpp::wrap(res));

}
