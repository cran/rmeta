
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
    		  call = match.call(), names = mf$names,
    		  conf.level = conf.level )
    class( rval ) <- "meta.MH"
    rval
}

# Print out the output of meta.MH function.
print.meta.MH <- function( obj ) {
    cat( "Fixed effects ( Mantel-Haenszel ) Meta-Analysis\n" )
    cat( "Call: " )
    print( obj$call )
    conf.level <- obj$conf.level
    ci.value <- -qnorm( ( 1 - conf.level ) / 2 )
    #lower <- mn - ci.value * se
    #upper <- mn + ci.value * se
    ci <- exp( obj$logMH + c( -ci.value, 0, ci.value ) * obj$selogMH )
    cat( paste( "Mantel-Haensel OR = ", round( ci[2],2 ), "    ",
    	        conf.level * 100, "% CI ( ", round( ci[1],2 ), 
    	        ", ", round( ci[3],2 ), " )\n", sep = "" 
    	      ) 
       )
    cat( paste( "Test for heterogeneity: X^2( ",obj$het[2]," ) = ",
    	        round( obj$het[1],2 ), " ( p-value ", 
    	        round( obj$het[3],4 ), " )\n", sep="" ) )
} 

# Summarise meta.MH function, same with print.meta.MH, but with odds
# ratios for each variable.
summary.meta.MH <- function( obj ,conf.level=NULL) {
    if (is.null(conf.level))
        conf.level <- obj$conf.level
    if (conf.level>1 & conf.level<100)
        conf.level<-conf.level/100
    ci.value <- -qnorm( ( 1 - conf.level ) / 2 )
    m <- exp( outer( obj$selogOR, c( 0, -ci.value, ci.value ), "*" ) +
              cbind( obj$logOR, obj$logOR, obj$logOR ) 
            )
    dimnames( m ) <- list( as.character( obj$names ),
                          c( "OR", "(lower ", paste(100*conf.level,"% upper)",sep=""))  )
    rval <- list( ors = m, call = obj$call, 
    		  MHci = exp( obj$logMH + c( -ci.value, 0, ci.value ) *
    		              obj$selogMH 
    			    ),
    		  het = obj$het,
                  conf.level=conf.level
    	        )
    class( rval ) <- "summary.meta.MH"
    rval
}

# Print out the summarised meta.MH function.
print.summary.meta.MH <- function( obj ) {
    cat( "Fixed effects ( Mantel-Haenszel ) meta-analysis\n" )
    cat( "Call: " )
    print( obj$call )
    cat( "------------------------------------\n" )
    print( round( obj$ors,2 ) )
    cat( "------------------------------------\n" )
    conf.level <- obj$conf.level
    cat( paste( "Mantel-Haensel OR = ", round( obj$MHci[2], 2 ), " ",
    	        conf.level*100, "% CI ( ", round( obj$MHci[1],2 ), ",", 
    	        round( obj$MHci[3],2 ), " )\n", sep="" 
    	      ) 
       )
    cat( paste( "Test for heterogeneity: X^2( ",obj$het[2]," ) = ",
       	        round( obj$het[1],2 ), " ( p-value ",
       	        round( obj$het[3],4 ), " )\n", sep = "" 
       	      ) 
       )
}

# Plot the Odds Ratio of meta.MH
plot.meta.MH <- function( obj, summary = T, summlabel = "Summary", conf.level=NULL,colors=list(box="black",lines="gray",summary="black",zero="lightgray"),... ){
    if (is.null(conf.level))
        conf.level <- obj$conf.level
    if (conf.level>1 & conf.level<100)
        conf.level<-conf.level/100
    if (is.null(colors))
        colors<-par("fg")
    if (!is.list(colors))
        colors<-list(box=colors,lines=colors,summary=colors,zero=colors)

    if ( summary )
        metaplot( obj$logOR, obj$selogOR, labels=obj$names,
        	  summn = obj$logMH, sumse = obj$selogMH, 
        	  conf.level = conf.level, sumnn = obj$selogMH ^ ( -2 ),
        	  summlabel = summlabel, logeffect = T,colors=colors,... )
    else 
       metaplot( obj$logOR, obj$selogOR, labels=obj$names, logeffect = T,colors=colors,... )
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
                  names = mf$names, conf.level = conf.level,
                  omitted = !ok, tau2=tau2 )
    class( rval ) <- "meta.DSL"
    rval
}

# Print out the output of meta.DSL function.
print.meta.DSL <- function( obj ){
    conf.level <- obj$conf.level
    ci.value <- -qnorm( ( 1 - conf.level ) / 2 )
    cat( "Random effects ( DerSimonian-Laird ) meta-analysis\n" )
    cat( "Call: " )
    print( obj$call )
    ci <- exp( obj$logDSL + c( -ci.value, 0, ci.value ) * obj$selogDSL )
    cat( paste( "Summary OR = ", round( ci[2],2 ), "    ",
    	        conf.level * 100, "% CI ( ",round( ci[1],2 ), ", ",
    	        round( ci[3],2 )," )\n", sep="" ) )
    cat( paste( "Estimated random effects variance:",round( obj$tau2,2 ),"\n" ) )
} 

# Summarise meta.DSL function, same with print.meta.DSL, but with odds
# ratios for each variable.
summary.meta.DSL <- function( obj,conf.level=NULL ){
    if (is.null(conf.level))
        conf.level <- obj$conf.level
    if (conf.level>1 & conf.level<100)
        conf.level<-conf.level/100
    ci.value <- -qnorm( ( 1 - conf.level ) / 2 )
    m <- exp( outer( obj$selogOR, c( 0,-ci.value,ci.value ), "*" ) +
    		     cbind( obj$logOR, obj$logOR, obj$logOR ) )
    dimnames( m ) <- list( as.character( obj$names ),
                          c( "OR", "(lower ", paste(100*conf.level,"% upper)",sep=""))  )
    rval <- list( ors = m, call = obj$call, 
                  ci = exp( obj$logDSL + c( -ci.value, 0, ci.value ) *
                	    obj$selogDSL ),
                  het = obj$het, 
                  omitted = as.character( obj$names )[obj$omitted],
                  conf.level=conf.level,
                  tau2 = obj$tau2 )
    class( rval ) <- "summary.meta.DSL"
    rval
}

# Print out the summarised meta.MH function.
print.summary.meta.DSL <- function( obj ) {
    conf.level <- obj$conf.level
    cat( "Random effects ( DerSimonian-Laird ) meta-analysis\n" )
    cat( "Call: " )
    print( obj$call )
    cat( "------------------------------------\n" )
    print( round( obj$ors,2 ) )
    cat( "------------------------------------\n" )
    cat( paste( "Summary OR = ",round( obj$ci[2],2 ), "  ",
    		conf.level * 100, "% CI ( ",round( obj$ci[1],2 ), ",",
    		round( obj$ci[3],2 )," )\n",sep="" ) )
    cat( paste( "Test for heterogeneity: X^2( ",obj$het[2]," ) = ",
    		round( obj$het[1],2 ),
    		" ( p-value ", round( obj$het[3],4 ), " )\n", sep="" ) )
    cat( paste( "Estimated random effects variance:",
    		round( obj$tau2,2 ), "\n" ) )
    if ( length( obj$omitted )>0 ){
        cat( paste( "( ", length( obj$omitted ),
    		    "studies with zero or infinite odds ratio omitted )\n" ) )
    }
}

# Plot the Odds Ratio of meta.DSL
plot.meta.DSL <- function( obj,summary=T,summlabel="Summary",conf.level=NULL,colors=list(box="black",lines="gray",summary="black",zero="lightgray"),... ){
    if (is.null(conf.level))
        conf.level <- obj$conf.level
    if (conf.level>1 & conf.level<100)
        conf.level<-conf.level/100
    if (is.null(colors))
        colors<-par("fg")
    if (!is.list(colors))
        colors<-list(box=colors,lines=colors,summary=colors,zero=colors)
    if ( summary )
        metaplot( obj$logOR, obj$selogOR, conf.level = conf.level,
                 labels = obj$names, summn = obj$logDSL,
            	  sumse = obj$selogDSL, sumnn = obj$selogDSL^( -2 ),
            	  summlabel = summlabel, logeffect = T,colors=colors,... )
    else 
        metaplot( obj$logOR, obj$selogOR, labels = obj$names,
        	  logeffect = T,colors=colors,... )
 }
 
 
metaplot <- function( mn, se, nn=NULL, labels=NULL, conf.level = .95,
		      xlab = "Odds ratio", ylab = "Study Reference",
		      boxwidth = NULL, xlim = NULL, summn = NULL,
		      sumse = NULL, sumnn = NULL, 
		      summlabel = "Summary", logeffect = F,
		      lwd = 2, boxsize = 1, 
		      zero = if ( logeffect ) 
		                 1 
		             else 
		                 0,
                      colors=NULL,
		      ... ) {
    ci.value <- -qnorm( ( 1 - conf.level ) / 2 )
    ok <- is.finite( mn + se )
    if ( is.null( xlim ) ) 
        xlim <- c( min( mn[ok] - ci.value * se[ok], na.rm = T ),
      		     max( mn[ok] + ci.value * se[ok], na.rm = T ) )
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
    axis( 1,at = round( 10 ^ pretty( log( exp( xlim ),10 ), 6 ), 2 ) )
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
        text( rep( nxlim[1], n ), -( 1:n ), labels,..., adj=0 )
    if ( is.null( nn ) ) 
        nn <- se ^ -2
    yscale <- 0.3 * boxsize / max( sqrt( nn ), na.rm = T )
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
        text( nxlim[1], -y0, labels = summlabel, adj = 0 )
    }
}
