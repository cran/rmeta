
meta.MH <- function ( ntrt, nctrl, ptrt, pctrl, conf.level = .95,
		      names = NULL, data = NULL, 
		      subset = NULL, na.action = na.fail ) {
    if (conf.level>1 & conf.level<100)
        conf.level<-conf.level/100
    if ( is.null( data ) ) 
        data <- sys.frame( sys.parent() )
    mf <- match.call()
    mf$data <- NULL 
    mf$subset <- NULL
    mf[[1]] <- as.name( "data.frame" )
    mf <- eval( mf,data )
    if ( !is.null( subset ) ) 
        mf <- mf[subset,]
    mf <- na.action( mf )
    A <- mf$ptrt
    B <- mf$pctrl
    C <- mf$ntrt - mf$ptrt
    D <- mf$nctrl - mf$pctrl
    logORs <- log( ( A * D ) / ( B * C ) )
    vars <- 1/A + 1/B + 1/C + 1/D
    Ti <- A + B + C + D
    P <- ( A + D ) / Ti
    Q <- ( B + C ) / Ti
    R <- A * D / Ti
    S <- B * C / Ti
    logMH <- log( sum( R ) / sum( S ) )
    varMH <- sum( P * R ) / ( 2 * sum( R )^2 ) + 
    	     sum( P * S + Q * R ) / ( 2 * sum( R ) * sum( S ) ) + 
    	     sum( Q * S ) / ( 2 * sum( S ) ^ 2 )
    MHchisq <- ( sum( A - ( A + B ) * ( A + C ) / Ti )^2 ) /
    		 sum( ( A + B ) * ( C + D ) * ( A + C ) * ( B + D ) /
    		      ( Ti * Ti * ( Ti - 1 ) 
    		    ) 
    	       )
    ok <- is.finite( logORs )
    heterog <- sum( ( ( logORs[ok] - logMH ) ^ 2 ) / vars[ok] )
    df <- sum( ok ) - 1
    rval <- list( logOR = logORs, selogOR = sqrt( vars ), 
    		  logMH = logMH, selogMH = sqrt( varMH ), 
    		  MHtest = c( MHchisq, 1 - pchisq( MHchisq, 1 ) ), 
    		  het = c( heterog, df, 1 - pchisq( heterog, df ) ),
    		  call = match.call(), names = as.character(mf$names),
    		  conf.level = conf.level )
    class( rval ) <- "meta.MH"
    rval
}

# Print out the output of meta.MH function.
print.meta.MH <- function( x, ... ) {
    cat( "Fixed effects ( Mantel-Haenszel ) Meta-Analysis\n" )
    cat( "Call: " )
    print( x$call )
    conf.level <- x$conf.level
    ci.value <- -qnorm( ( 1 - conf.level ) / 2 )
    #lower <- mn - ci.value * se
    #upper <- mn + ci.value * se
    ci <- exp( x$logMH + c( -ci.value, 0, ci.value ) * x$selogMH )
    cat( paste( "Mantel-Haensel OR = ", round( ci[2],2 ), "    ",
    	        conf.level * 100, "% CI ( ", round( ci[1],2 ), 
    	        ", ", round( ci[3],2 ), " )\n", sep = "" 
    	      ) 
       )
    cat( paste( "Test for heterogeneity: X^2( ",x$het[2]," ) = ",
    	        round( x$het[1],2 ), " ( p-value ", 
    	        round( x$het[3],4 ), " )\n", sep="" ) )
} 

# Summarise meta.MH function, same with print.meta.MH, but with odds
# ratios for each variable.
summary.meta.MH <- function( object ,conf.level=NULL, ...) {
    if (is.null(conf.level))
        conf.level <- object$conf.level
    if (conf.level>1 & conf.level<100)
        conf.level<-conf.level/100
    ci.value <- -qnorm( ( 1 - conf.level ) / 2 )
    m <- exp( outer( object$selogOR, c( 0, -ci.value, ci.value ), "*" ) +
              cbind( object$logOR, object$logOR, object$logOR ) 
            )
    dimnames( m ) <- list( as.character( object$names ),
                          c( "OR", "(lower ", paste(100*conf.level,"% upper)",sep=""))  )
    rval <- list( ors = m, call = object$call, 
    		  MHci = exp( object$logMH + c( -ci.value, 0, ci.value ) *
    		              object$selogMH 
    			    ),
    		  het = object$het,
                  conf.level=conf.level
    	        )
    class( rval ) <- "summary.meta.MH"
    rval
}

# Print out the summarised meta.MH function.
print.summary.meta.MH <- function( x, ... ) {
    cat( "Fixed effects ( Mantel-Haenszel ) meta-analysis\n" )
    cat( "Call: " )
    print( x$call )
    cat( "------------------------------------\n" )
    print( round( x$ors,2 ) )
    cat( "------------------------------------\n" )
    conf.level <- x$conf.level
    cat( paste( "Mantel-Haensel OR = ", round( x$MHci[2], 2 ), " ",
    	        conf.level*100, "% CI ( ", round( x$MHci[1],2 ), ",", 
    	        round( x$MHci[3],2 ), " )\n", sep="" 
    	      ) 
       )
    cat( paste( "Test for heterogeneity: X^2( ",x$het[2]," ) = ",
       	        round( x$het[1],2 ), " ( p-value ",
       	        round( x$het[3],4 ), " )\n", sep = "" 
       	      ) 
       )
}

meta.colors<-function(all.elements,box="black",lines="gray",summary="black",zero="lightgray",mirror="lightblue",text="black"){
    if (!missing(all.elements)){
        if (is.null(all.elements))
            all.elements<-par("fg")
        list(box=all.elements,lines=all.elements,summary=all.elements,zero=all.elements,mirror=all.elements,text=all.elements)
    }
    else
        list(box=box,lines=lines,summary=summary,zero=zero,mirror=mirror,text=text)    
}


# Plot the Odds Ratio of meta.MH
plot.meta.MH <- function( x, summary = TRUE, summlabel = "Summary", conf.level=NULL,colors=meta.colors(),... ){
    if (is.null(conf.level))
        conf.level <- x$conf.level
    if (conf.level>1 & conf.level<100)
        conf.level<-conf.level/100
 
    if ( summary )
        metaplot( x$logOR, x$selogOR, labels=x$names,
        	  summn = x$logMH, sumse = x$selogMH, 
        	  conf.level = conf.level, sumnn = x$selogMH ^ ( -2 ),
        	  summlabel = summlabel, logeffect = TRUE,colors=colors,... )
    else 
       metaplot( x$logOR, x$selogOR, labels=x$names, logeffect = TRUE,colors=colors,... )
 }
 
# Meta Analysis using DerSimonian-Laird Method.
meta.DSL <- function( ntrt, nctrl, ptrt, pctrl, conf.level = .95,
		      names = NULL, data = NULL,
		      subset = NULL, na.action = na.fail ) {
    if (conf.level>1 & conf.level<100)
        conf.level<-conf.level/100
    if ( is.null( data ) ) 
        data <- sys.frame( sys.parent() )
    mf <- match.call()
    mf$data <- NULL 
    mf$subset <- NULL
    mf[[1]] <- as.name( "data.frame" )
    mf <- eval( mf,data )
    if ( !is.null( subset ) ) 
        mf <- mf[subset,]
    mf <- na.action( mf )
    A <- mf$ptrt
    B <- mf$pctrl
    C <- mf$ntrt - mf$ptrt
    D <- mf$nctrl - mf$pctrl
    logORs <- log( ( A * D ) / ( B * C ) )
    ok <- is.finite( logORs )
    if ( any( !ok ) )
        warning( "Studies with 0/Inf odds ratio omitted" )
    vars <- 1/A + 1/B + 1/C + 1/D
    tau2 <- max( 0, var( logORs[ok] ) - sum( vars[ok] ) / ( length( ok )-1 ) )
    wts <- 1 / ( vars + tau2 )
    logDSL <- sum( wts[ok] * logORs[ok] ) / sum( wts[ok] )
    varDSL <- 1 / sum( wts[ok] )
    summary.test <- logDSL / sqrt( varDSL )
    heterog <- sum( ( ( logORs[ok] - logDSL )^2 ) / vars[ok] )
    df <- sum( ok ) - 1
    rval <- list( logOR = logORs, selogOR = sqrt( vars ), 
    		  logDSL = logDSL, selogDSL = sqrt( varDSL ), 
    		  test = c( summary.test,
    		   	    1 - pchisq( summary.test ^ 2, 1 ) ), 
    		  het = c( heterog, df, 1 - pchisq( heterog, df ) ),
                  call = match.call(),
                  names = as.character(mf$names), conf.level = conf.level,
                  omitted = !ok, tau2=tau2 )
    class( rval ) <- "meta.DSL"
    rval
}

# Print out the output of meta.DSL function.
print.meta.DSL <- function( x, ... ){
    conf.level <- x$conf.level
    ci.value <- -qnorm( ( 1 - conf.level ) / 2 )
    cat( "Random effects ( DerSimonian-Laird ) meta-analysis\n" )
    cat( "Call: " )
    print( x$call )
    ci <- exp( x$logDSL + c( -ci.value, 0, ci.value ) * x$selogDSL )
    cat( paste( "Summary OR = ", round( ci[2],2 ), "    ",
    	        conf.level * 100, "% CI ( ",round( ci[1],2 ), ", ",
    	        round( ci[3],2 )," )\n", sep="" ) )
    cat( paste( "Estimated random effects variance:",round( x$tau2,2 ),"\n" ) )
} 

# Summarise meta.DSL function, same with print.meta.DSL, but with odds
# ratios for each variable.
summary.meta.DSL <- function( object,conf.level=NULL, ... ){
    if (is.null(conf.level))
        conf.level <- object$conf.level
    if (conf.level>1 & conf.level<100)
        conf.level<-conf.level/100
    ci.value <- -qnorm( ( 1 - conf.level ) / 2 )
    m <- exp( outer( object$selogOR, c( 0,-ci.value,ci.value ), "*" ) +
    		     cbind( object$logOR, object$logOR, object$logOR ) )
    dimnames( m ) <- list( as.character( object$names ),
                          c( "OR", "(lower ", paste(100*conf.level,"% upper)",sep=""))  )
    rval <- list( ors = m, call = object$call, 
                  ci = exp( object$logDSL + c( -ci.value, 0, ci.value ) *
                	    object$selogDSL ),
                  het = object$het, 
                  omitted = as.character( object$names )[object$omitted],
                  conf.level=conf.level,
                  tau2 = object$tau2 )
    class( rval ) <- "summary.meta.DSL"
    rval
}

# Print out the summarised meta.MH function.
print.summary.meta.DSL <- function( x, ... ) {
    conf.level <- x$conf.level
    cat( "Random effects ( DerSimonian-Laird ) meta-analysis\n" )
    cat( "Call: " )
    print( x$call )
    cat( "------------------------------------\n" )
    print( round( x$ors,2 ) )
    cat( "------------------------------------\n" )
    cat( paste( "Summary OR = ",round( x$ci[2],2 ), "  ",
    		conf.level * 100, "% CI ( ",round( x$ci[1],2 ), ",",
    		round( x$ci[3],2 )," )\n",sep="" ) )
    cat( paste( "Test for heterogeneity: X^2( ",x$het[2]," ) = ",
    		round( x$het[1],2 ),
    		" ( p-value ", round( x$het[3],4 ), " )\n", sep="" ) )
    cat( paste( "Estimated random effects variance:",
    		round( x$tau2,2 ), "\n" ) )
    if ( length( x$omitted )>0 ){
        cat( paste( "( ", length( x$omitted ),
    		    "studies with zero or infinite odds ratio omitted )\n" ) )
    }
}

# Plot the Odds Ratio of meta.DSL
plot.meta.DSL <- function( x,summary=TRUE,summlabel="Summary",conf.level=NULL,colors=meta.colors(),... ){
    if (is.null(conf.level))
        conf.level <- x$conf.level
    if (conf.level>1 & conf.level<100)
        conf.level<-conf.level/100
    if ( summary )
        metaplot( x$logOR, x$selogOR, conf.level = conf.level,
                 labels = x$names, summn = x$logDSL,
            	  sumse = x$selogDSL, sumnn = x$selogDSL^( -2 ),
            	  summlabel = summlabel, logeffect = TRUE,colors=colors,... )
    else 
        metaplot( x$logOR, x$selogOR, labels = x$names,
        	  logeffect = TRUE,colors=colors,conf.level=conf.level,... )
 }
 

meta.summaries<-function(d,se,method=c("fixed","random"),
			weights=NULL,logscale=FALSE,names=NULL,data=NULL,
			conf.level=.95,subset=NULL,na.action=na.fail)
{
      if (conf.level>1 & conf.level<100)
        conf.level<-conf.level/100
    if ( is.null( data ) ) 
        data <- sys.frame( sys.parent() )
    mf <- match.call()
    mf$data <- NULL 
    mf$method<-NULL  
    mf$logscale<-NULL
    mf$subset <- NULL
    mf[[1]] <- as.name( "data.frame" )
    mf <- eval( mf,data )
    if ( !is.null( subset ) ) 
        mf <- mf[subset,]
    mf <- na.action( mf )
  
    if (is.null(mf$names)){
	if (is.null(mf$d))
	   mf$names<-seq(along=mf$d)
	else
	  mf$names<-names(mf$d)
    }
   mf$names<-as.character(mf$names)
   method<-match.arg(method)
   vars<-mf$se^2
   tau2<-max(0,var(mf$d)-sum(vars)/NROW(mf))
   if(is.null(mf$weights)){
	if (method=="fixed"){
   		wt<-1/vars
	} else {
		wt<-1/(vars+tau2)
	}
   } else
	wt<-mf$weights

  summ<-sum(wt*mf$d)/sum(wt)
  if (method=="fixed")
	varsum<-sum(wt*wt*vars)/(sum(wt)^2)
  else
	varsum<-sum(wt*wt*(vars+tau2))/(sum(wt)^2)
  
  summtest<-summ/sqrt(varsum)
  heterog<-sum(((mf$d-summ)^2)/vars)
  df<-length(vars)-1
  rval<-list(effects=mf$d, stderrs=mf$se, summary=summ,se.summary=sqrt(varsum),
	     test=c(summtest,1-pchisq(summtest^2,1)),
	     het=c(heterog,df,1-pchisq(heterog,df)),
	     call=match.call(), names=mf$names,tau2=tau2,
	     variance.method=method, weights=wt, 
	     weight.method=if(is.null(mf$weights)) method else "user",
	     conf.level=conf.level,logscale=logscale)
  class(rval)<-"meta.summaries"
  rval
}

summary.meta.summaries <- function( object ,conf.level=NULL, ...) {
    if (is.null(conf.level))
        conf.level <- object$conf.level
    if (conf.level>1 & conf.level<100)
        conf.level<-conf.level/100
    ci.value <- -qnorm( ( 1 - conf.level ) / 2 )
    m <- outer( object$stderrs, c( 0, -ci.value, ci.value ), "*" ) +
              cbind( object$effects, object$effects, object$effects )
    label<-"Effect"
    if (object$logscale){
        m<-exp(m)
        label<-"exp(Effect)"
    }
    wt<-object$weights
    m<-cbind(m, round(wt*length(wt)/sum(wt),1))
    dimnames( m ) <- list( as.character( object$names ),
                          c( label, "(lower ", paste(100*conf.level,"% upper)",sep=""),"weights")  )
    summci<-object$summary + c( -ci.value, 0, ci.value ) *object$se.summary
    if (object$logscale) summci<-exp(summci)
    rval <- list( studies = m, call = object$call, 
    		  summci = summci,
    		  het = object$het,tau2=object$tau2,
                  conf.level=conf.level,logscale=object$logscale,
                 weight.method=object$weight.method,
                 variance.method=object$variance.method
    	        )
    class( rval ) <- "summary.meta.summaries"
    rval
}

# Print out the summarised meta.summaries function.
print.summary.meta.summaries <- function( x, ... ) {
    if (x$weight.method=="user")
        cat(paste("Weighted meta-analysis with ",x$variance.method,"-effects standard errors\n",sep=""))
    else
        cat(paste(switch(x$variance.method,fixed="Fixed",random="Random"),"-effects meta-analysis\n",sep=""))

    cat( "Call: " )
    print( x$call )
    cat( "----------------------------------------------------\n" )
    print( round( x$studies,2 ) )
    cat( "----------------------------------------------------\n" )
    conf.level <- x$conf.level
    if (x$logscale)
        cat("Summary exp(effect): ")
    else
        cat("Summary effect: ")
    cat( paste(  round( x$summci[2], 2 ), " ",
    	        conf.level*100, "% CI ( ", round( x$summci[1],2 ), ",", 
    	        round( x$summci[3],2 ), " )\n", sep="" 
    	      ) 
       )
    if (x$variance.method=="fixed")
        cat( paste( "Test for heterogeneity: X^2( ",x$het[2]," ) = ",
                   round( x$het[1],2 ), " ( p-value ",
                   round( x$het[3],4 ), " )\n", sep = "" 
                   ) 
            )
    else
        cat("Estimated heterogeneity variance:",signif(x$tau2,2)," p=",round(x$het[3],3),"\n")
}

print.meta.summaries<-function(x, ...){
   conf.level<-x$conf.level
   ci.value<- -qnorm((1-conf.level)/2)
   if (x$variance.method=="fixed")
	cat("Fixed-effects meta-analysis")
   else	
	cat("Random-effects meta-analysis")
   if(x$weight.method!="user")
	cat("\n")
   else
	cat(" with user-supplied weights\n")
   cat("Call: ")
   print(x$call)
   ci<-x$summary+c(-ci.value,0,ci.value)*x$se.summary
   cat(paste("Summary effect=",signif(ci[2],3), "   ",
	conf.level*100,"% CI (",signif(ci[1],3),", ",signif(ci[3],3),
	")\n",sep=""))
   cat("Estimated heterogeneity variance:",signif(x$tau2,2)," p=",round(x$het[3],3),"\n")
}

plot.meta.summaries<-function(x,summary=TRUE,summlabel="Summary",
	conf.level=NULL,colors=meta.colors(),xlab=NULL,logscale=NULL,... )
{
   if (is.null(conf.level))
        conf.level <- x$conf.level
   if (conf.level>1 & conf.level<100)
        conf.level<-conf.level/100
   if(is.null(logscale)) logscale<-x$logscale
   if (is.null(xlab)){
       if (logscale) xlab<-"exp(Effect)" else xlab<-"Effect"
   }
   if ( summary )
        metaplot( x$effects, x$stderrs, conf.level = conf.level,
                 labels = x$names, summn = x$summary,
            	  sumse = x$se.summary, sumnn = sum(x$weights),
            	  summlabel = summlabel, logeffect = logscale,
		colors=colors,nn=x$weights,xlab,... )
    else 
        metaplot( x$effects, x$stderrs, labels = x$names,
        	  logeffect = logscale,colors=colors,
		nn=sum(x$weights),conf.level=conf.level,xlab,... )

}




metaplot <- function( mn, se, nn=NULL, labels=NULL, conf.level = .95,
		      xlab = "Odds ratio", ylab = "Study Reference",
		       xlim = NULL, summn = NULL,
		      sumse = NULL, sumnn = NULL, 
		      summlabel = "Summary", logeffect = FALSE,
		      lwd = 2, boxsize = 1, 
		      zero = if ( logeffect ) 
		                 1 
		             else 
		                 0,
                      colors=meta.colors(),
		      ... ) {
    ci.value <- -qnorm( ( 1 - conf.level ) / 2 )
    ok <- is.finite( mn + se )
    if ( is.null( xlim ) ) 
        xlim <- c( min( mn[ok] - ci.value * se[ok], na.rm = TRUE ),
      		     max( mn[ok] + ci.value * se[ok], na.rm = TRUE ) )
    ##par( pty="s" )
    n <- length( mn )
    if ( logeffect ) {
        xlog <- "x"
        nxlim <- exp( xlim )
    }
    else {
        xlog <- ""
        nxlim <- xlim
    }
    if ( !is.null( labels ) ) {
        if ( logeffect ) 
            nxlim[1] <- nxlim[1] / sqrt( nxlim[2] / nxlim[1] )
        else
          nxlim[1] <- nxlim[1] - 0.5 * ( nxlim[2] - nxlim[1] )
        
    }
    par( xaxt = "n",yaxt = "n" )
    plot( nxlim,c( 1,-n-2-3 * !is.null( summn ) ),
          type = "n", bty = "n", xaxt = "n", yaxt = "n",
          log = xlog, xlab=xlab,ylab=ylab,... )
    par( xaxt = "s" )
    if (logeffect)
    	axis( 1,at = round( 10 ^ pretty( log( exp( xlim ),10 ), 6 ), 2 ) )
    else
    	axis( 1,at = pretty(xlim, 6 ) )
	
    if ( !is.null( zero ) )
        abline( v = zero, lty = 2, lwd = 2 ,col=colors$zero)
    ci.value <- -qnorm( ( 1 - conf.level ) / 2 )
    lower <- mn - ci.value * se
    upper <- mn + ci.value * se
    if ( logeffect ){
        lower <- exp( lower )
        upper <- exp( upper )
    }
    for ( i in 1:n ){
        if ( is.na( lower[i]+upper[i] ) ) 
            next
        lines( c( lower[i], upper[i] ), c( -i, -i ), lwd = lwd, col=colors$lines,... )
    }
    if ( !is.null( labels ) )
        text( rep( nxlim[1], n ), -( 1:n ), labels,..., col=colors$text,adj=0 )
    if ( is.null( nn ) ) 
        nn <- se ^ -2
    yscale <- 0.3 * boxsize / max( sqrt( nn ), na.rm = TRUE )
    if ( logeffect ) { 
        scale <- ( nxlim[2] / nxlim[1] ) ^ ( yscale / ( 4 + n ) )
        xl <- exp( mn ) * ( scale ^ -sqrt( nn ) )
        xr <- exp( mn ) * ( scale ^ sqrt( nn ) )
    }
    else {
        scale <- yscale * ( nxlim[2] - nxlim[1] ) / ( 4 + n )
        xl <- mn - scale * sqrt( nn )
        xr <- mn + scale * sqrt( nn )
    }
    yb <- ( 1:n ) - yscale * sqrt( nn )
    yt <- ( 1:n ) + yscale * sqrt( nn )
    for ( i in 1:n ) {
        if ( !is.finite( mn[i] ) ) 
            next  
     rect( xl[i], -yb[i], xr[i], -yt[i], col = colors$box,border=  colors$box)
    }
    if ( !is.null( summn ) ) {
        if ( logeffect ) {
            x0 <- exp( summn )
            xl <- exp( summn - ci.value * sumse )
            xr <- exp( summn + ci.value * sumse )
        }
        else{
            x0 <- summn
            xl <- summn - ci.value * sumse
            xr <- summn + ci.value * sumse
        }
        y0 <- n + 3
        yb <- n + 3 - sqrt( sumnn ) * yscale
        yt <- n + 3 + sqrt( sumnn ) * yscale
        polygon( c( xl, x0, xr, x0 ), -c( y0, yt, y0, yb ),
    	         col = colors$summary, border = colors$summary )
        text( nxlim[1], -y0, labels = summlabel, adj = 0,col=colors$text )
    }
}

funnelplot<-function(x,...)
    UseMethod("funnelplot")

funnelplot.meta.MH<-function(x,...){
    funnelplot.default(x$logOR,x$selogOR,summ=x$logMH,...)
}
funnelplot.meta.DSL<-function(x,...){
    funnelplot.default(x$logOR,x$selogOR,summ=x$logDSL,...)
}

funnelplot.meta.summaries<-function(x,...){
    funnelplot.default(x$effects,x$stderrs,summ=x$summary,...)
}


funnelplot.default<-function(x,se,size=1/se,summ=NULL,xlab="Effect",ylab="Size",
		colors=meta.colors(),
		conf.level=0.95,plot.conf=FALSE,zero=NULL,mirror=FALSE,...)
{
   finite<-function(x) x[is.finite(x)]
    
   if (mirror && is.null(summ))
	stop("Can't do a mirror plot without a summary value")

   if (plot.conf){
	ci<--qnorm((1-conf.level)/2)
	xlim<-range(finite(c(zero,x-ci*se,x+ci*se)))
	if (mirror)
	  xlim<-range(finite(c(xlim,2*summ-x-ci*se,2*summ-x+ci*se)))
   }else{
      xlim<-range(finite(c(zero,x)))
      if (mirror)
	xlim<-range(finite(c(xlim,2*summ-x,2*summ+x)))
   }
   plot(x,size,ylim=c(0,max(size)*1.1),xlim=xlim,xlab=xlab,
	ylab=ylab,col=if(is.null(colors$points)) par("fg") else colors$points)

   if (plot.conf)
       segments(x-ci*se,size,x+ci*se,size,col=if(is.null(colors$conf)) par("fg") else colors$conf,lwd=2)
   if (!is.null(summ))
	   abline(v=summ,col=if(is.null(colors$summary)) par("fg") else colors$summary,lty=2,lwd=2)

   if(!is.null(zero))
	abline(v=zero,col=if(is.null(colors$zero)) par("fg") else colors$zero,lwd=2)

   if(mirror){
	points(2*summ-x,size,col=if(is.null(colors$mirror)) par("fg") else colors$mirror)
       if (plot.conf)
	segments(2*summ-x-ci*se,size,2*summ-x+ci*se,size,col=if(is.null(colors$mirror)) par("fg") else colors$mirror,lwd=2)
  }
}




