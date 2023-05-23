plot_acet <- function(acet, boot = FALSE, heri = FALSE, xlab, ylab, main, col, legend = TRUE)
{
	if(!(class(acet) %in% c('AtCtEt_model', 'AtDtEt_model', 'AtCtEp_mc_model','AtEtp_mc_model','AtCtEtp_mc_model')))
	{
		stop('The first argument must be an AtCtEt_model, AtDtEt_model, or AtCtEtp_mc_model object.')
	}
  
  if(missing(xlab))
  {
    xlab_t <- 'Age'
  }else{
    xlab_t <- xlab
  }
  
	if(missing(ylab))
	{
	  if(heri == FALSE)
    {
      ylab_t <- 'Variance'     
     }else{
      ylab_t <- 'Heritability'   
    }
	}else{
	  ylab_t <- ylab
	}
  
	if(missing(main))
	{
	  if(heri == FALSE)
	  {  
	    if(is(acet,'AtDtEt_model')==FALSE)
	    {main_t <- "Variance curves of the A, C, and E components"}else{
	      main_t <- "Variance curves of the A, D, and E components"
	    }
	  }else{   
	    main_t <- "Dynamic heritability"
	  }
	}else{
	  main_t <- main
	}
  
  if(missing(col))
  {
    if(heri == FALSE)
    {  
      col <- c("red","blue","pink")
    }else{   
      col <- "black"
    }
  }else
  {
    if(heri == FALSE)
    {
      if(length(col)<3)
      {
        stop('The \'col\' argument should have three elements.')
      }
    }
  }

	# if(class(acet)=='AtCtEt_model')
  if(is(acet,'AtCtEt_model'))
	{
		if(heri == FALSE)
		{
			plot_AtCtEt(acet, boot, xlab=xlab_t, ylab=ylab_t, main=main_t, col=col, legend =legend)
		}else{
			plot_AtCtEt_h(acet, boot, xlab=xlab_t, ylab=ylab_t, main=main_t, col=col)
		}
  }
  
  if(is(acet,'AtDtEt_model'))
  {
    if(heri == FALSE)
    {
      plot_AtDtEt(acet, boot, xlab=xlab_t, ylab=ylab_t, main=main_t, col=col, legend =legend)
    }else{
      plot_AtDtEt_h(acet, boot, xlab=xlab_t, ylab=ylab_t, main=main_t, col=col)
    }
  }
	
	# if(class(acet)=='AtCtEtp_mc_model')
  if(is(acet,'AtCtEtp_mc_model'))
	{
		if(heri==FALSE)
		{
			plot_AtCtEtp(acet, xlab=xlab_t, ylab=ylab_t, main=main_t, col=col, legend =legend)
		}else{
			plot_AtCtEt_h(acet, boot, xlab=xlab_t, ylab=ylab_t, main=main_t, col=col)
		}
	}


}