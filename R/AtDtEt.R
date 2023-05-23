AtDtEt <-
  function(data_m, data_d, mod = c('d','d','d'), knot_a=5, knot_d=5, knot_e=5, loc = c('e','e','e'), boot=FALSE, num_b = 100, init = rep(0,3), robust = 0)
  {
    
    pheno_m <- c(t(data_m[,1:2]))
    pheno_d <- c(t(data_d[,1:2]))
    T_m <- rep(data_m[,3], each=2)
    T_d <- rep(data_d[,3], each=2)
    
    mag <- var(pheno_m)
    init_max <- log(mag)
    init_min <- log(mag) - abs(log(mag))*1.3
    
    if((is.vector(mod)==FALSE) | (length(mod)!=3) )
    {stop('The \'mod\' argument must be a vector of length 3.')}
    
    if(!(mod[1] %in% c('d','c','n')))
    {stop('The \'mod\' argument for the A component must be \'d\'(dynamic), \'c\'(constant) or \'n\'(NA).')}
    
    if(!(mod[2] %in% c('d','c','n')))
    {stop('The \'mod\' argument for the D component must be \'d\'(dynamic), \'c\'(constant) or \'n\'(NA).')}
    
    if(!(mod[3] %in% c('d','c')))
    {stop('The \'mod\' argument for the E component must be \'d\'(dynamic), \'c\'(constant).')}
    
    if((is.vector(loc)==FALSE) | (length(loc)!=3) )
    {stop('The \'loc\' argument must be a vector of length 3.')}
    
    order <- 3
    if(mod[1]=='d')
    {
      order <- 3
      
      if(knot_a < 3)
      {stop('The number of interior knots must be no less than 3.')}
      
    }else
    {
      order <- 1
    }
    #knot <- 8
    min_T <- min(T_m, T_d)
    max_T <- max(T_m, T_d)
    
    if(mod[1]=='d')
    {
      if(loc[1]=='e')
      {
        knots_a <- seq(from=min_T, to=max_T, length.out=knot_a)
        interval_a <- knots_a[2] - knots_a[1]
        knots_a <- c(c(min_T-interval_a*2,min_T-interval_a), knots_a)
        knots_a <- c(knots_a, c(max_T+interval_a,max_T+interval_a*2))
      }else{
        knots_a <- quantile(unique(T_m,T_d), probs = seq(from=0,to=1,length.out=knot_a))
        knots_a <- c(knots_a[1], knots_a[1], knots_a)
        knots_a <- c(knots_a, knots_a[knot_a+2], knots_a[knot_a+2])
      }
      
      B_des_a_m <- splineDesign(knots_a, x=T_m, ord=order)
      B_des_a_d <- splineDesign(knots_a, x=T_d, ord=order)
    }else{
      knots_a <- c(min_T,max_T)
      B_des_a_m <- splineDesign(knots_a, x=T_m, ord=order)
      B_des_a_d <- splineDesign(knots_a, x=T_d, ord=order)
    }
    
    if(mod[2]=='d')
    {
      order <- 3
      
      if(knot_d < 3)
      {stop('The number of interior knots must be no less than 3.')}
      
    }else
    {
      order <- 1
    }
    if(mod[2]=='d')
    {
      if(loc[2]=='e')
      {
        knots_d <- seq(from=min_T, to=max_T, length.out=knot_d)
        interval_c <- knots_d[2] - knots_d[1]
        knots_d <- c(c(min_T-interval_c*2,min_T-interval_c), knots_d)
        knots_d <- c(knots_d, c(max_T+interval_c,max_T+interval_c*2))
      }else{
        knots_d <- quantile(unique(T_m,T_d), probs = seq(from=0,to=1,length.out=knot_d))
        knots_d <- c(knots_d[1], knots_d[1], knots_d)
        knots_d <- c(knots_d, knots_d[knot_d+2], knots_d[knot_d+2])
      }
      
      B_des_d_m <- splineDesign(knots_d, x=T_m, ord=order)
      B_des_d_d <- splineDesign(knots_d, x=T_d, ord=order)
    }else{
      knots_d <- c(min(T_m, T_d),max(T_m, T_d))
      B_des_d_m <- splineDesign(knots_d, x=T_m, ord=order)
      B_des_d_d <- splineDesign(knots_d, x=T_d, ord=order)
    }
    
    if(mod[3]=='d')
    {
      order <- 3
      
      if(knot_e < 3)
      {stop('The number of interior knots must be no less than 3.')}
      
    }else
    {
      order <- 1
    }
    if(mod[3]=='d')
    {
      if(loc[3]=='e')
      {
        knots_e <- seq(from=min_T, to=max_T, length.out=knot_e)
        interval_e <- knots_e[2] - knots_e[1]
        knots_e <- c(c(min_T-interval_e*2,min_T-interval_e), knots_e)
        knots_e <- c(knots_e, c(max_T+interval_e,max_T+interval_e*2))
      }else{
        knots_e <- quantile(unique(T_m,T_d), probs = seq(from=0,to=1,length.out=knot_e))
        knots_e <- c(knots_e[1], knots_e[1], knots_e)
        knots_e <- c(knots_e, knots_e[knot_e+2], knots_e[knot_e+2])
      }
      
      B_des_e_m <- splineDesign(knots_e, x=T_m, ord=order)
      B_des_e_d <- splineDesign(knots_e, x=T_d, ord=order)
    }else{
      knots_e <- c(min(T_m, T_d),max(T_m, T_d))
      B_des_e_m <- splineDesign(knots_e, x=T_m, ord=order)
      B_des_e_d <- splineDesign(knots_e, x=T_d, ord=order)
      
    }
    
    n_d <- ncol(B_des_d_m)
    n_a <- ncol(B_des_a_m)
    n_e <- ncol(B_des_e_m)
    
    init_a <- rep(init[1],n_a)
    init_d <- rep(init[2],n_d)
    init_e <- rep(init[3],n_e)
    
    up_a <- up_d <- up_e <- 10
    lo_a <- lo_d <- -50
    lo_e <- -15
    
    if(mod[1]=='n')
    {
      up_a <- lo_a <- -50
      init_a <- -50
    }
    
    if(mod[1]=='c')
    {
      up_a <- 20
      lo_a <- -50
    }
    
    if(mod[2]=='n')
    {
      up_d <- lo_d <- -50
      init_d <- -50
    }
    
    if(mod[2]=='c')
    {
      up_d <- 20
      lo_d <- -50
    }
    
    if(mod[3]=='c')
    {
      up_e <- 20
      lo_e <- -10
    }
    

    result <- optim(c(init_a,init_d,init_e), loglik_AtDtEt_esp, gr_AtDtEt_esp, pheno_m = matrix(pheno_m), pheno_d = matrix(pheno_d), B_des_a_m = B_des_a_m, B_des_a_d = B_des_a_d, B_des_d_m = B_des_d_m, B_des_d_d = B_des_d_d, B_des_e_m = B_des_e_m, B_des_e_d = B_des_e_d,lower = c(rep(lo_a, n_a),rep(lo_d, n_d),rep(lo_e, n_e)), upper = c(rep(up_a, n_a),rep(up_d, n_d),rep(up_e, n_e)), method = "L-BFGS-B", hessian = TRUE, control=list(maxit = 3000))

    if(robust>0)
    {
      for(i in 1:ceiling(robust))
      {
        init <- runif(n_a+n_d+n_e,min=init_min,max=init_max)
        if(mod[1]!='n')
        {init_a <- init[1:n_a]}
        if(mod[2]!='n')
        {init_d <- init[(n_a+1):(n_a+n_d)]}
        init_e <- init[(n_a+n_d+1):(n_a+n_d+n_e)]
        
        result_r <- optim(c(init_a,init_d,init_e), loglik_AtDtEt_esp, gr_AtDtEt_esp, pheno_m = matrix(pheno_m), pheno_d = matrix(pheno_d), B_des_a_m = B_des_a_m, B_des_a_d = B_des_a_d, B_des_d_m = B_des_d_m, B_des_d_d = B_des_d_d, B_des_e_m = B_des_e_m, B_des_e_d = B_des_e_d,lower = c(rep(lo_a, n_a),rep(lo_d, n_d),rep(lo_e, n_e)), upper = c(rep(up_a, n_a),rep(up_d, n_d),rep(up_e, n_e)), method = "L-BFGS-B", hessian = TRUE, control=list(maxit = 3000))
      
        if(result_r$value < result$value)
        {
          result <- result_r
        }
      }
    }
    
    res_a <- result$par[1:n_a]
    res_d <- result$par[(1+n_a):(n_d+n_a)]
    res_e <- result$par[(1+n_a+n_d):(n_e+n_d+n_a)]
    
    hes <- .Call('hessian_AtDtEt_esp_c', res_a, res_d, res_e, matrix(pheno_m), matrix(pheno_d), B_des_a_m, B_des_a_d, B_des_d_m, B_des_d_d, B_des_e_m, B_des_e_d)
      
    n_t <- n_a+n_d+n_e
    hes_m <- matrix(0, n_t, n_t)
    k <- 1
    for(i in 1:n_t)
    {
      for(j in i:n_t)
      {
        hes_m[i,j] <- hes[k]
        k <- k + 1
      }
    }
    
    hes_m_t <- t(hes_m)
    diag(hes_m_t) <- 0
    hes_m <- hes_m_t + hes_m
    
    if(mod[1]=='n')
    {res_a <- -Inf}
    if(mod[2]=='n')
    {res_d <- -Inf}
    
    AtCtEt_model <- list(n_beta_a=n_a, n_beta_d=n_d, n_beta_e=n_e, beta_a=res_a, beta_d=res_d, beta_e=res_e, hessian_ap=result$hessian, hessian=hes_m, con=result$convergence, lik=result$value, knots_a =knots_a, knots_d = knots_d, knots_e = knots_e, min_t = min_T, max_t = max_T, boot = NULL )
    class(AtCtEt_model) <- 'AtDtEt_model'
    
    if(boot==TRUE)
    {
      boot_res <- AtDtEt_boot(res = AtCtEt_model, mod, data_m, data_d, knot_a, knot_d, knot_e, loc, B=num_b,alpha=0.05,m=500)
      AtCtEt_model$boot <- boot_res
    }
    
    return(invisible(AtCtEt_model))
    
  }
