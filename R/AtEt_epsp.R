loglik_AtEt_epsp <-
function(param, pheno_m, pheno_d, B_des_a_m, B_des_a_d, beta_a, D_a, B_des_e_m, B_des_e_d, beta_e, D_e)
{
var_b_a <- param[1]
var_b_e <- param[2]

nll <- .Call('loglik_AtEt_epsp_c', var_b_a, var_b_e, pheno_m, pheno_d, B_des_a_m, B_des_a_d, beta_a, D_a, B_des_e_m, B_des_e_d, beta_e, D_e)
return(nll)
}

loglik_AtEt_epsp_g <-
function(param, pheno_m, pheno_d, B_des_a_m, B_des_a_d, B_des_e_m, B_des_e_d, var_b_a, var_b_e, D_a, D_e)
{
beta_a <- param[1:ncol(D_a)]
beta_e <- param[(ncol(D_a)+1):(ncol(D_a)+ncol(D_e))]

nll <- .Call('loglik_AtEt_epsp_g_c', beta_a, beta_e, pheno_m, pheno_d, B_des_a_m, B_des_a_d, B_des_e_m, B_des_e_d, var_b_a, var_b_e, D_a, D_e)
return(nll)
}

gr_AtEt_epsp_g <- 
function(param, pheno_m, pheno_d, B_des_a_m, B_des_a_d, B_des_e_m, B_des_e_d, var_b_a, var_b_e, D_a, D_e)
{
	beta_a <- param[1:ncol(D_a)]
	beta_e <- param[(ncol(D_a)+1):(ncol(D_a)+ncol(D_e))]
	
d <- .Call('gr_AtEt_epsp_g_c', beta_a, beta_e, pheno_m, pheno_d, B_des_a_m, B_des_a_d, B_des_e_m, B_des_e_d, var_b_a, var_b_e, D_a, D_e)
		
	return(d)
}

gr_AtEt_epsp <- function(param, pheno_m, pheno_d, B_des_a_m, B_des_a_d, beta_a, D_a, B_des_e_m, B_des_e_d, beta_e, D_e)
{

	var_b_a <- param[1]
	var_b_e <- param[2]
	d <- .Call('gr_AtEt_epsp_c', var_b_a, var_b_e, pheno_m, pheno_d, B_des_a_m, B_des_a_d, beta_a, D_a, B_des_e_m, B_des_e_d, beta_e, D_e)
	return(d)
}

