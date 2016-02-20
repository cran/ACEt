AtEtp <-
function(data_m, data_d, knot_a=8, knot_e=8, eps = 0.1)
{
num_m <- nrow(data_m)*2
num_d <- nrow(data_d)*2
pheno_m <- matrix(NA, num_m, 1)
pheno_d <- matrix(NA, num_d, 1)
pheno_m[seq(from=1, to=num_m, by=2),1] <- data_m[,1]
pheno_m[seq(from=2, to=num_m, by=2),1] <- data_m[,2]
pheno_d[seq(from=1, to=num_d, by=2),1] <- data_d[,1]
pheno_d[seq(from=2, to=num_d, by=2),1] <- data_d[,2]
T_m <- rep(data_m[,3],each=2)
T_d <- rep(data_d[,3],each=2)

order <- 3
#knot_a <- 12
#knot_c <- 12
penal <- 2

delta_a <- matrix(0, knot_a+order-2-penal, knot_a+order-2)
for(i in 1:nrow(delta_a))
{
	delta_a[i, i:(i+2)] <- c(1,-2,1)
}
D_a <- t(delta_a)%*%delta_a

knots_a <- seq(from=min(T_m, T_d), to=max(T_m, T_d), length.out=knot_a)
interval_a <- knots_a[2] - knots_a[1]
knots_a <- c(c(min(T_m, T_d)-interval_a*2,min(T_m, T_d)-interval_a), knots_a)
knots_a <- c(knots_a, c(max(T_m, T_d)+interval_a,max(T_m, T_d)+interval_a*2))
B_des_a_m <- splineDesign(knots_a, x=T_m, ord=order)	
B_des_a_d <- splineDesign(knots_a, x=T_d, ord=order)

delta_e <- matrix(0, knot_e+order-2-penal, knot_e+order-2)
for(i in 1:nrow(delta_e))
{
	delta_e[i, i:(i+2)] <- c(1,-2,1)
}
D_e <- t(delta_e)%*%delta_e

knots_e <- seq(from=min(T_m, T_d), to=max(T_m, T_d), length.out=knot_e)
interval_e <- knots_e[2] - knots_e[1]
knots_e <- c(c(min(T_m, T_d)-interval_e*2,min(T_m, T_d)-interval_e), knots_e)
knots_e <- c(knots_e, c(max(T_m, T_d)+interval_e,max(T_m, T_d)+interval_e*2))
B_des_e_m <- splineDesign(knots_e, x=T_m, ord=order)	
B_des_e_d <- splineDesign(knots_e, x=T_d, ord=order)	

n_a <- ncol(D_a)
n_e <- ncol(D_e)

var_b_a <- 1
var_b_e <- 1
beta_a <- rep(0, n_a)
beta_e <- rep(0, n_e)

lik <- 100000
lik_pre <- 200000

liks <- c()
betas <- matrix(0,0,n_a+n_e)
vars <- matrix(0,0,2)

while(abs(lik-lik_pre)>eps)
{
lik_pre <- lik
result <- optim(c(beta_a,beta_e), loglik_AtEt_epsp_g, gr_AtEt_epsp_g, pheno_m = pheno_m, pheno_d = pheno_d, B_des_a_m=B_des_a_m, B_des_a_d=B_des_a_d, B_des_e_m=B_des_e_m, B_des_e_d=B_des_e_d, var_b_a=var_b_a, var_b_e=var_b_e, D_a=D_a, D_e=D_e, lower = rep(-5,n_a+n_e), upper = rep(5,n_a+n_e), method = "L-BFGS-B", control=list(maxit = 3000))
betas <- rbind(betas,result$par)
beta_a <- result$par[1:n_a]
beta_e <- result$par[(1+n_a):(n_a+n_e)]
result <- optim(c(var_b_a,var_b_e), loglik_AtEt_epsp, gr_AtEt_epsp, pheno_m = pheno_m, pheno_d = pheno_d, B_des_a_m=B_des_a_m, B_des_a_d=B_des_a_d, B_des_e_m=B_des_e_m, B_des_e_d=B_des_e_d, beta_a=beta_a, beta_e=beta_e, D_a=D_a, D_e=D_e, lower = c(0.0001, 0.0001), upper = c(100, 100), method = "L-BFGS-B", control=list(maxit = 3000))
vars <- rbind(vars, result$par)
var_b_a <- result$par[1]
var_b_e <- result$par[2]
lik <- result$value
liks <- c(liks, result$value)

}

min_i <- match(min(liks), liks)

result <- optim(betas[min_i,], loglik_AtEt_epsp_g, gr_AtEt_epsp_g, pheno_m = pheno_m, pheno_d = pheno_d, B_des_a_m=B_des_a_m, B_des_a_d=B_des_a_d, B_des_e_m=B_des_e_m, B_des_e_d=B_des_e_d, var_b_a=vars[min_i,1], var_b_e=vars[min_i,2], D_a=D_a, D_e=D_e, lower = rep(-10,n_e+n_a), upper = rep(10,n_e+n_a), method = "L-BFGS-B", control=list(maxit = 3000))

AtEtp_model <- list(D_a = D_a, D_e = D_e, pheno_m = pheno_m, pheno_d = pheno_d, T_m = T_m, T_d = T_d, knot_a=knots_a, knot_e=knots_e, beta_a=result$par[1:n_a], beta_e=result$par[(1+n_a):(n_e+n_a)], con=result$convergence, lik=result$value, iter=liks, var_b_a=vars[min_i,1], var_b_e=vars[min_i,2])

class(AtEtp_model) <- 'AtEtp_model'

return(AtEtp_model)

}
