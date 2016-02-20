AtCtEt <-
function(data_m, data_d, model = c('d','d','d'), knot_a=5, knot_c=5, knot_e=5, boot=FALSE, num_b = 100)
{

pheno_m <- c(t(data_m[,1:2]))
pheno_d <- c(t(data_d[,1:2]))
T_m <- rep(data_m[,3], each=2)
T_d <- rep(data_d[,3], each=2)

if((is.vector(model)==FALSE) | (length(model)!=3) )
{stop('The model parameter must be a vector of length 3.')}

if(!(model[1] %in% c('d','c','n')))
{stop('The \'model\' parameter for the A component must be \'d\'(dynamic), \'c\'(constant) or \'n\'(NA).')}

if(!(model[2] %in% c('d','c','n')))
{stop('The \'model\' parameter for the C component must be \'d\'(dynamic), \'c\'(constant) or \'n\'(NA).')}

if(!(model[3] %in% c('d','c')))
{stop('The \'model\' parameter for the E component must be \'d\'(dynamic), \'c\'(constant).')}

order <- 3
if(model[1]=='d')
{
order <- 3

if(knot_a < 3)
{stop('The number of interior knots must be no less than 3.')}

}else
{
order <- 1
}
#knot <- 8
if(model[1]=='d')
{
knots_a <- seq(from=min(T_m, T_d), to=max(T_m, T_d), length.out=knot_a)
interval_a <- knots_a[2] - knots_a[1]
knots_a <- c(c(min(T_m, T_d)-interval_a*2,min(T_m, T_d)-interval_a), knots_a)
knots_a <- c(knots_a, c(max(T_m, T_d)+interval_a,max(T_m, T_d)+interval_a*2))
B_des_a_m <- splineDesign(knots_a, x=T_m, ord=order)
B_des_a_d <- splineDesign(knots_a, x=T_d, ord=order)
}else{
knots_a <- c(min(T_m, T_d),max(T_m, T_d))
B_des_a_m <- splineDesign(knots_a, x=T_m, ord=order)
B_des_a_d <- splineDesign(knots_a, x=T_d, ord=order)
}

if(model[2]=='d')
{
order <- 3

if(knot_c < 3)
{stop('The number of interior knots must be no less than 3.')}

}else
{
order <- 1
}
if(model[2]=='d')
{
knots_c <- seq(from=min(T_m, T_d), to=max(T_m, T_d), length.out=knot_c)
interval_c <- knots_c[2] - knots_c[1]
knots_c <- c(c(min(T_m, T_d)-interval_c*2,min(T_m, T_d)-interval_c), knots_c)
knots_c <- c(knots_c, c(max(T_m, T_d)+interval_c,max(T_m, T_d)+interval_c*2))
B_des_c_m <- splineDesign(knots_c, x=T_m, ord=order)
B_des_c_d <- splineDesign(knots_c, x=T_d, ord=order)
}else{
knots_c <- c(min(T_m, T_d),max(T_m, T_d))
B_des_c_m <- splineDesign(knots_c, x=T_m, ord=order)
B_des_c_d <- splineDesign(knots_c, x=T_d, ord=order)
}

if(model[3]=='d')
{
order <- 3

if(knot_e < 3)
{stop('The number of interior knots must be no less than 3.')}

}else
{
order <- 1
}
if(model[3]=='d')
{
knots_e <- seq(from=min(T_m, T_d), to=max(T_m, T_d), length.out=knot_e)
interval_e <- knots_e[2] - knots_e[1]
knots_e <- c(c(min(T_m, T_d)-interval_e*2,min(T_m, T_d)-interval_e), knots_e)
knots_e <- c(knots_e, c(max(T_m, T_d)+interval_e,max(T_m, T_d)+interval_e*2))
B_des_e_m <- splineDesign(knots_e, x=T_m, ord=order)
B_des_e_d <- splineDesign(knots_e, x=T_d, ord=order)
}else{
knots_e <- c(min(T_m, T_d),max(T_m, T_d))
B_des_e_m <- splineDesign(knots_e, x=T_m, ord=order)
#knots_e <- c(min(T_d),max(T_d))
B_des_e_d <- splineDesign(knots_e, x=T_d, ord=order)

}

n_c <- ncol(B_des_c_m)
n_a <- ncol(B_des_a_m)
n_e <- ncol(B_des_e_m)

init_a <- 0
init_c <- 0
init_e <- 0

up_a <- up_c <- up_e <- 5
lo_a <- lo_c <- -Inf
lo_e <- -20

if(model[1]=='n')
{
up_a <- lo_a <- 0
init_a <- 0
}

if(model[1]=='c')
{
up_a <- 1000
lo_a <- 0
init_a <- 1
}

if(model[2]=='n')
{
up_c <- lo_c <- 0
init_c <- 0
}

if(model[2]=='c')
{
up_c <- 1000
lo_c <- 0
init_c <- 1
}

if(model[3]=='c')
{
up_e <- 1000
lo_e <- 0.000001
init_e <- 1
}

result <- optim(c(rep(init_a, n_a),rep(init_c, n_c),rep(init_e, n_e)), loglik_AtCtEt_esp, gr_AtCtEt_esp, pheno_m = matrix(pheno_m), pheno_d = matrix(pheno_d), B_des_a_m = B_des_a_m, B_des_a_d = B_des_a_d, B_des_c_m = B_des_c_m, B_des_c_d = B_des_c_d, B_des_e_m = B_des_e_m, B_des_e_d = B_des_e_d,lower = c(rep(lo_a, n_a),rep(lo_c, n_c),rep(lo_e, n_e)), upper = c(rep(up_a, n_a),rep(up_c, n_c),rep(up_e, n_e)), method = "L-BFGS-B", hessian = TRUE, control=list(maxit = 3000))

res_a <- result$par[1:n_a]
res_c <- result$par[(1+n_a):(n_c+n_a)]
res_e <- result$par[(1+n_a+n_c):(n_e+n_c+n_a)]

AtCtEt_model <- list(n_beta_a=n_a, n_beta_c=n_c, n_beta_e=n_e, beta_a=res_a, beta_c=res_c, beta_e=res_e, hessian=result$hessian, con=result$convergence, lik=result$value, knots_a =knots_a, knots_c = knots_c, knots_e = knots_e, min_t = min(T_m, T_d), max_t = max(T_m, T_d), boot = NULL )

if(AtCtEt_model$con!=0)
{
warning('The optimization algorithm might not converge.')
}

class(AtCtEt_model) <- 'AtCtEt_model'

if(boot==TRUE)
{
	boot_res <- AtCtEt_boot(res = AtCtEt_model, model, data_m, data_d, knot_a, knot_c, knot_e, B=num_b,alpha=0.05,m=500)
	AtCtEt_model$boot <- boot_res
}

return(invisible(AtCtEt_model))

}
