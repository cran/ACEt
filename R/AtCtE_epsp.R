loglik_AtCtE_epsp <-
function(param, pheno_m, pheno_d, B_des_a_m, B_des_a_d, beta_a, D_a, B_des_c_m, B_des_c_d, beta_c, D_c)
{
var <- param[1]
var_b_a <- param[2]
var_b_c <- param[3]

nll <- .Call('loglik_AtCtE_epsp_c', var, var_b_a, var_b_c, pheno_m, pheno_d, B_des_a_m, B_des_a_d, beta_a, D_a, B_des_c_m, B_des_c_d, beta_c, D_c)

return(nll)
}

loglik_AtCtE_epsp_g <-
function(param, pheno_m, pheno_d, B_des_a_m, B_des_a_d, B_des_c_m, B_des_c_d, var, var_b_a, var_b_c, D_a, D_c)
{
beta_a <- param[1:ncol(D_a)]
beta_c <- param[(ncol(D_a)+1):(ncol(D_a)+ncol(D_c))]

nll <- .Call('loglik_AtCtE_epsp_g_c', beta_a, beta_c, pheno_m, pheno_d, B_des_a_m, B_des_a_d, B_des_c_m, B_des_c_d, var, var_b_a, var_b_c, D_a, D_c)

return(nll)
}

gr_AtCtE_epsp <-
function(param, pheno_m, pheno_d, B_des_a_m, B_des_a_d, beta_a, D_a, B_des_c_m, B_des_c_d, beta_c, D_c)
{
var <- param[1]
var_b_a <- param[2]
var_b_c <- param[3]

d <- .Call('gr_AtCtE_epsp_c', var, var_b_a, var_b_c, pheno_m, pheno_d, B_des_a_m, B_des_a_d, beta_a, D_a, B_des_c_m, B_des_c_d, beta_c, D_c)

return(d)
}

gr_AtCtE_epsp_g <-
function(param, pheno_m, pheno_d, B_des_a_m, B_des_a_d, B_des_c_m, B_des_c_d, var, var_b_a, var_b_c, D_a, D_c)
{
beta_a <- param[1:ncol(D_a)]
beta_c <- param[(ncol(D_a)+1):(ncol(D_a)+ncol(D_c))]

d <- .Call('gr_AtCtE_epsp_g_c', beta_a, beta_c, pheno_m, pheno_d, B_des_a_m, B_des_a_d, B_des_c_m, B_des_c_d, var, var_b_a, var_b_c, D_a, D_c)

return(d)
}
