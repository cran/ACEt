plot_AtCtEt <- function(AtCtEt)
{
	if(class(AtCtEt)!='AtCtEt_model')
	{
		stop('The first parameter must be an object obtained from the AtCtEt function.')
	}
	
	model_cur <- AtCtEt

	#pheno_m <- c(t(data_m[,1:2]))
	#pheno_d <- c(t(data_d[,1:2]))
	#T_m <- rep(data_m[,3], each=2)
	#T_d <- rep(data_d[,3], each=2)

	order <- 3
	x <- seq(from=model_cur$min_t, to=model_cur$max_t, length.out=500)
	if(model_cur$n_beta_c>1)
	{
		bb_c <- splineDesign(model_cur$knots_c, x = x, ord=order, outer.ok = TRUE)
	}else{
		bb_c <- splineDesign(model_cur$knots_c, x = x, ord=1, outer.ok = TRUE)
	}
	points_c <- exp(bb_c%*%model_cur$beta_c)
	if(model_cur$n_beta_a>1)
	{
		bb_a <- splineDesign(model_cur$knots_a, x = x, ord=order, outer.ok = TRUE)
	}else{
		bb_a <- splineDesign(model_cur$knots_a, x = x, ord=1, outer.ok = TRUE)
	}
	points_a <- exp(bb_a%*%model_cur$beta_a)
	if(model_cur$n_beta_e>1)
	{
		bb_e <- splineDesign(model_cur$knots_e, x = x, ord=order, outer.ok = TRUE)
	}else{
		bb_e <- splineDesign(model_cur$knots_e, x = x, ord=1, outer.ok = TRUE)
	}
	points_e <- exp(bb_e%*%model_cur$beta_e)

	l_a <- model_cur$n_beta_a
	l_c <- model_cur$n_beta_c
	l_e <- model_cur$n_beta_e
	
	h_v <- c()
	if((l_a==1)&(model_cur$beta_a[1]==-Inf))
	{
		h_v <- c(h_v,FALSE)
	}else{
		h_v <- c(h_v, rep(TRUE, l_a))
	}
	if((l_c==1)&(model_cur$beta_c[1]==-Inf))
	{
		h_v <- c(h_v,FALSE)
	}else{
		h_v <- c(h_v, rep(TRUE, l_c))
	}
	h_v <- c(h_v, rep(TRUE, l_e))

	#fisher <- solve(model_cur$hessian[2:(1+l_a+l_c),2:(1+l_a+l_c)])
	fisher <- solve(model_cur$hessian[h_v,h_v])

	plot(range(x), c(0,max(points_c, points_a, points_e)+1), type = "n", xlab = "Age", ylab = "Variance",main =  "Variance curves of the A, C, and E components")
	
	index <- 1
	# bb <- splineDesign(model_cur$knots_a, x = x, ord=order, outer.ok = TRUE)
	lines(x, points_a, col = "red", lwd = 2)
	if((l_a>1)|(model_cur$beta_a[1]!=-Inf))
	{
	
	fisher_a <- fisher[index:l_a,index:l_a]
	index <- index + l_a
	lower <- rep(NA, length(x))
	upper <- rep(NA, length(x))
	sd <- rep(NA, length(x))

	for(i in 1:length(x))
	{
		sd[i] <- sqrt(exp(2*sum(bb_a[i,]*model_cur$beta_a))*(t(bb_a[i,])%*%fisher_a%*%bb_a[i,]))
		lower[i] <- exp(sum(bb_a[i,]*model_cur$beta_a)) - 1.96*sd[i]
		upper[i] <- exp(sum(bb_a[i,]*model_cur$beta_a)) + 1.96*sd[i]
	}

	lines(x, lower, col = "orange" ,lty = 2 , lwd = 0.6)
	lines(x, upper, col = "orange" ,lty = 2 , lwd = 0.6)
	polygon(c(x, rev(x)),c(upper, rev(lower)),col='grey',border = NA, lty=3, density=20)
	}
	
	lines(x, points_c, col = "blue", lwd = 2)
	if((l_c>1)|(model_cur$beta_c[1]!=-Inf))
	{
	fisher_c <- fisher[index:(index+l_c-1),index:(index+l_c-1)]
	index <- index + l_c
	lower <- rep(NA, length(x))
	upper <- rep(NA, length(x))
	sd <- rep(NA, length(x))

	for(i in 1:length(x))
	{
		sd[i] <- sqrt(exp(2*sum(bb_c[i,]*model_cur$beta_c))*(t(bb_c[i,])%*%fisher_c%*%bb_c[i,]))
		lower[i] <- exp(sum(bb_c[i,]*model_cur$beta_c)) - 1.96*sd[i]
		upper[i] <- exp(sum(bb_c[i,]*model_cur$beta_c)) + 1.96*sd[i]
	}

	lines(x, lower, col = "green" ,lty = 2 , lwd = 0.6)
	lines(x, upper, col = "green" ,lty = 2 , lwd = 0.6)
	polygon(c(x, rev(x)),c(upper, rev(lower)),col='grey',border = NA, lty=3, density=20)

	}
	
	lines(x, points_e, col = "pink", lwd = 2)
	
	fisher_e <- fisher[index:(index+l_e-1),index:(index+l_e-1)]
	lower <- rep(NA, length(x))
	upper <- rep(NA, length(x))
	sd <- rep(NA, length(x))

	for(i in 1:length(x))
	{
		sd[i] <- sqrt(exp(2*sum(bb_e[i,]*model_cur$beta_e))*(t(bb_e[i,])%*%fisher_e%*%bb_e[i,]))
		lower[i] <- exp(sum(bb_e[i,]*model_cur$beta_e)) - 1.96*sd[i]
		upper[i] <- exp(sum(bb_e[i,]*model_cur$beta_e)) + 1.96*sd[i]
	}
	
	lines(x, lower, col = "yellow" ,lty = 2 , lwd = 0.6)
	lines(x, upper, col = "yellow" ,lty = 2 , lwd = 0.6)
	polygon(c(x, rev(x)),c(upper, rev(lower)),col='grey',border = NA, lty=3, density=20)

	legend(x[1], max(points_c, points_a, points_e)+1, c('Additive genetic component','Common environmental component', 'Unique environmental component'), col = c('red','blue','pink'), lty=c(1,1,1), lwd=c(2,2,2))

}