## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----eval=FALSE---------------------------------------------------------------
#  install.packages("ACEt")

## ----eval=FALSE---------------------------------------------------------------
#  install.packages("devtools")
#  library(devtools)
#  install_github("lhe17/ACEt")

## -----------------------------------------------------------------------------
library(ACEt)
data(data_ace)


## -----------------------------------------------------------------------------
attributes(data_ace)
head(data_ace$mz)
head(data_ace$dz)

## -----------------------------------------------------------------------------
# fitting the ACE(t) model
re <- AtCtEt(data_ace$mz, data_ace$dz, mod = c('d','d','c'), knot_a = 6, knot_c = 4)
summary(re)

## -----------------------------------------------------------------------------
# part of the expected information matrix
re$hessian[1:8,1:8]
# part the observed information matrix approximated by the L-BFGS algorithm
re$hessian_ap[1:8,1:8]

## -----------------------------------------------------------------------------
re_cc <- AtCtEt(data_ace$mz, data_ace$dz, mod = c('d','c','c'), knot_a = 6, knot_c = 4)
p1 <- pchisq(2*(re_cc$lik-re$lik), 4, lower.tail=FALSE)
p1
re_ac <- AtCtEt(data_ace$mz, data_ace$dz, mod = c('c','d','c'), knot_a = 6, knot_c = 4)
p2 <- pchisq(2*(re_ac$lik-re$lik), 6, lower.tail=FALSE)
p2
re_cn <- AtCtEt(data_ace$mz, data_ace$dz, mod = c('d','n','c'), knot_a = 6, knot_c = 4)
p3 <- 0.5*pchisq(2*(re_cn$lik-re_cc$lik), 1, lower.tail=FALSE)
p3

## -----------------------------------------------------------------------------
plot_acet(re, ylab='Var', xlab='Age (1-50)')

## -----------------------------------------------------------------------------
## fitting an ACE(t) model with the CIs esitmated by the bootstrap method 
re_b <- AtCtEt(data_ace$mz, data_ace$dz, mod = c('d','d','c'), knot_a = 6, knot_c = 4, boot = TRUE, 
               num_b = 60)
plot_acet(re_b, boot = TRUE)

## -----------------------------------------------------------------------------
## plot dynamic heritability with the CIs using the delta method 
plot_acet(re_b, heri=TRUE, boot = FALSE)
## plot dynamic heritability with the CIs using the bootstrap method 
plot_acet(re_b, heri=TRUE, boot = TRUE)

## ----eval=FALSE---------------------------------------------------------------
#  ## fitting an ADE(t) model with the CIs esitmated by the bootstrap method
#  re_b <- AtDtEt(data_ace$mz, data_ace$dz, mod = c('d','d','c'), boot = TRUE, num_b = 60)
#  plot_acet(re_b, boot = TRUE)

## -----------------------------------------------------------------------------
## fitting an ACE(t)-p model
re <- AtCtEtp(data_ace$mz, data_ace$dz, knot_a = 8, knot_c = 8, mod=c('d','d','l'))
summary(re)

## -----------------------------------------------------------------------------
re_mcmc <- acetp_mcmc(re, iter_num = 5000, burnin = 500)
summary(re_mcmc)

## -----------------------------------------------------------------------------
plot_acet(re_mcmc)
plot_acet(re_mcmc, heri=TRUE)

## ----knot_10, echo=FALSE, fig.cap="Plots of variance curves of the example data set fitted by the ACE(t) and ACE(t)-p model with 10 interior knots for each component. Left: the ACE(t) model. Right: the ACE(t)-p model."----
knitr::include_graphics("knot_10.jpg")

## -----------------------------------------------------------------------------
test <- test_acetp(re, comp = 'e')
test$p


## ----eval=FALSE---------------------------------------------------------------
#  test <- test_acetp(re, comp = 'c', sim = 100, robust = 0)
#  test$p

