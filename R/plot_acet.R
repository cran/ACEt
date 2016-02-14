plot_acet <- function(acet)
{
	if(!(class(acet) %in% c('AtCtEt_model', 'AtCtEp_mc_model','AtEtp_mc_model','AtCtEtp_mc_model')))
	{
		stop('The first parameter must be an acet object.')
	}

	if(class(acet)=='AtCtEt_model')
	{
		plot_AtCtEt(acet)
	}
	
	
	if(class(acet)=='AtCtEp_mc_model')
	{
		plot_AtCtEp(acet)
	}

	if(class(acet)=='AtCtEtp_mc_model')
	{
		plot_AtCtEtp(acet)
	}

	if(class(acet)=='AtEtp_mc_model')
	{
		plot_AtEtp(acet)
	}

}