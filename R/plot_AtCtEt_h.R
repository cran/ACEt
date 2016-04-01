plot_AtCtEt_h <- function(AtCtEt, boot=FALSE)
{
	if((class(AtCtEt)!='AtCtEt_model')&(class(AtCtEt)!='AtCtEtp_mc_model'))
	{
		stop('The first parameter must be an object obtained from the AtCtEt or acetp_mcmc function.')
	}

	if((boot==TRUE)&(class(AtCtEt)=='AtCtEt_model'))
	{
		if(is.null(AtCtEt$boot)==TRUE)
		{
			stop('Please first run the AtCtEt model with bootstrapping.')
		}
	}
	
	if(class(AtCtEt)=='AtCtEt_model')
	{
		model_cur <- AtCtEt

		l_a <- model_cur$n_beta_a
		l_c <- model_cur$n_beta_c
		l_e <- model_cur$n_beta_e
		if((l_a==1)&(model_cur$beta_a[1]==-Inf))
		{stop('The current model has no additive genetic component.')}

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
	
		if(model_cur$n_beta_a>1)
		{
			bb_a <- splineDesign(model_cur$knots_a, x = x, ord=order, outer.ok = TRUE)
		}else{
			bb_a <- splineDesign(model_cur$knots_a, x = x, ord=1, outer.ok = TRUE)
		}
	
		if(model_cur$n_beta_e>1)
		{
			bb_e <- splineDesign(model_cur$knots_e, x = x, ord=order, outer.ok = TRUE)
		}else{
			bb_e <- splineDesign(model_cur$knots_e, x = x, ord=1, outer.ok = TRUE)
		}

		points_a <- exp(bb_a%*%model_cur$beta_a)
		points_c <- exp(bb_c%*%model_cur$beta_c)
		points_e <- exp(bb_e%*%model_cur$beta_e)
		points_h <- points_a/(points_a+points_c+points_e)	

		#fisher <- solve(model_cur$hessian[2:(1+l_a+l_c),2:(1+l_a+l_c)])
		fisher <- solve(model_cur$hessian)
	
		max_v <- max(points_h)*1.2
		plot(range(x), c(0,max_v), type = "n", xlab = "Age", ylab = "Heritability",main =  "Dynamic heritability")
	
		index <- 1
		# bb <- splineDesign(model_cur$knots_a, x = x, ord=order, outer.ok = TRUE)
		lines(x, points_h, col = "black", lwd = 2)
		#lines(x, points_c, col = "blue", lwd = 2)
		#lines(x, points_e, col = "pink", lwd = 2)

		if(boot == TRUE)
		{
			lines(AtCtEt$boot$x, AtCtEt$boot$lower.ci_h, col = "grey" ,lty = 2 , lwd = 0.6)
			lines(AtCtEt$boot$x, AtCtEt$boot$upper.ci_h, col = "grey" ,lty = 2 , lwd = 0.6)
			polygon(c(AtCtEt$boot$x, rev(AtCtEt$boot$x)),c(AtCtEt$boot$upper.ci_h, rev(AtCtEt$boot$lower.ci_h)),col='grey',border = NA, lty=3, density=20)
		
		}else{
	
			e_a <- exp(bb_c%*%model_cur$beta_c-bb_a%*%model_cur$beta_a)
			e_b <- exp(bb_e%*%model_cur$beta_e-bb_a%*%model_cur$beta_a)
			lower <- rep(NA, length(x))
			upper <- rep(NA, length(x))
			sd <- rep(NA, length(x))
			flag <- 0
			for(i in 1:length(x))
			{
				P <- matrix(NA,2,l_a+l_c+l_e)
				P[1,] <- c((-1)*bb_a[i,],bb_c[i,],rep(0,l_e))
				P[2,] <- c((-1)*bb_a[i,],rep(0,l_c),bb_e[i,])

				Sigma <- P%*%fisher%*%t(P)
	
				delta <- t(c(e_a[i],e_b[i]))%*%Sigma%*%c(e_a[i],e_b[i])
				delta <- delta*((1+e_a[i]+e_b[i])^(-4))
				if(delta>=0)
				{
				sd[i] <- sqrt(delta)
				esti <- 1/(1+e_a[i]+e_b[i])
				lower[i] <- esti - 1.96*sd[i]
				upper[i] <- esti + 1.96*sd[i]
				}else{flag <- 1}
			}
	
			lower <- ifelse(lower<0, 0, lower)
			upper <- ifelse(upper>max_v, max_v, upper)	
	
			if(flag == 0)
			{
			lines(x, lower, col = "grey" ,lty = 2 , lwd = 0.6)
			lines(x, upper, col = "grey" ,lty = 2 , lwd = 0.6)
			polygon(c(x, rev(x)),c(upper, rev(lower)),col='grey',border = NA, lty=3, density=20)
			}else{warning('Please try the bootstrap method for the confidence interval or use a different model.')}
	
		} # boot == FALSE
	
	# legend(x[1], max_v, c('Additive genetic component','Common environmental component', 'Unique environmental component'), col = c('red','blue','pink'), lty=c(1,1,1), lwd=c(2,2,2))
	}else{
		model_cur <- AtCtEt

		l_a <- length(model_cur$beta_a_mc)
		l_c <- length(model_cur$beta_c_mc)
		l_e <- length(model_cur$beta_e_mc)
		if((l_a==1)&(model_cur$beta_a_mc[1]==-Inf))
		{stop('The current model has no additive genetic component.')}


		order <- 3
		p_n <- 500
		x <- seq(from=model_cur$min_t, to=model_cur$max_t, length.out=p_n)
		t_int <- model_cur$max_t-model_cur$min_t
		l_m_1 <- (model_cur$max_t-x)/t_int
		l_m_2 <- (x-model_cur$min_t)/t_int

		if(length(model_cur$beta_a_mc)>2)
		{
			bb_a <- splineDesign(model_cur$knots_a, x = x, ord=order, outer.ok = TRUE)
		}else{
			if(length(model_cur$beta_a_mc)==2)
			{
				bb_a <- matrix(NA, p_n, 2)
				bb_a[,1] <- l_m_1
				bb_a[,2] <- l_m_2
			}else{
				bb_a <- matrix(1, p_n, 1)
			}
		}
		points_a <- exp(bb_a%*%model_cur$beta_a_mc)

		if(length(model_cur$beta_c_mc)>2)
		{
			bb_c <- splineDesign(model_cur$knots_c, x = x, ord=order, outer.ok = TRUE)
		}else{
			if(length(model_cur$beta_c_mc)==2)
			{
				bb_c <- matrix(NA, p_n, 2)
				bb_c[,1] <- l_m_1
				bb_c[,2] <- l_m_2
			}else{
				bb_c <- matrix(1, p_n, 1)
			}
		}
		points_c <- exp(bb_c%*%model_cur$beta_c_mc)

		if(length(model_cur$beta_e_mc)>2)
		{
			bb_e <- splineDesign(model_cur$knots_e, x = x, ord=order, outer.ok = TRUE)
		}else{
			if(length(model_cur$beta_e_mc)==2)
			{
				bb_e <- matrix(NA, p_n, 2)
				bb_e[,1] <- l_m_1
				bb_e[,2] <- l_m_2
			}else{
				bb_e <- matrix(1, p_n, 1)
			}
		}
		points_e <- exp(bb_e%*%model_cur$beta_e_mc)

		points_h <- points_a/(points_a+points_c+points_e)	

		
	}
}