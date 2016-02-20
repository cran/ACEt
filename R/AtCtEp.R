AtCtEp <-
function(data_m, data_d, knot_a=8, knot_c=8, eps = 0.1)
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

delta_c <- matrix(0, knot_c+order-2-penal, knot_c+order-2)
for(i in 1:nrow(delta_c))
{
delta_c[i, i:(i+2)] <- c(1,-2,1)
}
D_c <- t(delta_c)%*%delta_c

knots_c <- seq(from=min(T_m, T_d), to=max(T_m, T_d), length.out=knot_c)
interval_c <- knots_c[2] - knots_c[1]
knots_c <- c(c(min(T_m, T_d)-interval_c*2,min(T_m, T_d)-interval_c), knots_c)
knots_c <- c(knots_c, c(max(T_m, T_d)+interval_c,max(T_m, T_d)+interval_c*2))
B_des_c_m <- splineDesign(knots_c, x=T_m, ord=order)
B_des_c_d <- splineDesign(knots_c, x=T_d, ord=order)

n_a <- ncol(B_des_a_m)
n_c <- ncol(B_des_c_m)

var <- 1
var_b_a <- 1
var_b_c <- 1
beta_a <- rep(0, n_a)
beta_c <- rep(0, n_c)

lik <- 100000
lik_pre <- 200000

liks <- c()
betas <- matrix(0,0,n_a+n_c)
vars <- matrix(0,0,3)

while(abs(lik-lik_pre)>eps)
{
lik_pre <- lik
result <- optim(c(beta_a,beta_c), loglik_AtCtE_epsp_g, gr_AtCtE_epsp_g, pheno_m = pheno_m, pheno_d = pheno_d, B_des_a_m=B_des_a_m, B_des_a_d=B_des_a_d, B_des_c_m=B_des_c_m, B_des_c_d=B_des_c_d, var=var, var_b_a=var_b_a, var_b_c=var_b_c, D_a=D_a, D_c=D_c, lower = rep(-Inf,n_a+n_c), upper = rep(10,n_a+n_c), method = "L-BFGS-B", control=list(maxit = 3000))
betas <- rbind(betas, result$par)
beta_a <- result$par[1:n_a]
beta_c <- result$par[(1+n_a):(n_a+n_c)]
result <- optim(c(var,var_b_a,var_b_c), loglik_AtCtE_epsp, gr_AtCtE_epsp, pheno_m = pheno_m, pheno_d = pheno_d, B_des_a_m=B_des_a_m, B_des_a_d=B_des_a_d, B_des_c_m=B_des_c_m, B_des_c_d=B_des_c_d, beta_a=beta_a, beta_c=beta_c, D_a=D_a, D_c=D_c, lower = rep(0.00001,3), upper = rep(10,3), method = "L-BFGS-B", control=list(maxit = 3000))
vars <- rbind(vars, result$par)
var <- result$par[1]
var_b_a <- result$par[2]
var_b_c <- result$par[3]
lik <- result$value
liks <- c(liks, result$value)
#print(lik)
}

min_i <- match(min(liks), liks)

result <- optim(betas[min_i,], loglik_AtCtE_epsp_g, gr_AtCtE_epsp_g, pheno_m = matrix(pheno_m), pheno_d = matrix(pheno_d), B_des_a_m=B_des_a_m, B_des_a_d=B_des_a_d, B_des_c_m=B_des_c_m, B_des_c_d=B_des_c_d, var=vars[min_i,1], var_b_a=vars[min_i,2], var_b_c=vars[min_i,3], D_a=D_a, D_c=D_c, lower = rep(-Inf,n_c+n_a), upper = rep(10,n_c+n_a), method = "L-BFGS-B", control=list(maxit = 3000))

AtCtEp_model <- list(D_a = D_a, D_c = D_c, pheno_m = pheno_m, pheno_d = pheno_d, T_m = T_m, T_d = T_d, knot_a=knots_a, knot_c=knots_c, beta_a=result$par[1:n_a], beta_c=result$par[(1+n_a):(n_c+n_a)], con=result$convergence, lik=(result$value)/2, iter=(liks)/2, var=vars[min_i,1], var_b_a=vars[min_i,2], var_b_c=vars[min_i,3])

class(AtCtEp_model) <- 'AtCtEp_model'

print('Variance of the E component:')
print(var)
print('Estimates of beta_a:')
print(beta_a)
print('Estimates of beta_c')
print(beta_c)

return(invisible(AtCtEp_model))

}