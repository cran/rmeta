meta.DSL<-function(ntrt, nctrl, ptrt, pctrl,names=NULL,data=NULL,subset=NULL,na.action=na.fail)
 {
        if (is.null(data)) data<-sys.frame(sys.parent())
        mf<-match.call()
        mf$data<-NULL 
        mf$subset<-NULL
        mf[[1]]<-as.name("data.frame")
        mf<-eval(mf,data)
        if (!is.null(subset)) mf<-mf[subset,]
        mf<-na.action(mf)
        A <- mf$ptrt
        B <- mf$pctrl
        C <- mf$ntrt -mf$ptrt
        D <- mf$nctrl - mf$pctrl
        logORs <- log((A * D)/(B * C))
        ok <- is.finite(logORs)
        if (any(!ok))
          warning("Studies with 0/Inf odds ratio omitted")
        vars <- 1/A + 1/B + 1/C + 1/D
        tau2<-max(0,var(logORs[ok])-sum(vars[ok])/(length(ok)-1))
        wts<-1/(vars+tau2)
        logDSL <- sum(wts[ok]*logORs[ok])/sum(wts[ok])
        varDSL<-1/sum(wts[ok])
        summary.test<-logDSL/sqrt(varDSL)
        heterog <- sum(((logORs[ok] - logDSL)^2)/vars[ok])
        df <- sum(ok) - 1
        rval<-list(logOR = logORs, selogOR = sqrt(vars), logDSL = logDSL, 
                selogDSL = sqrt(varDSL), test = c(summary.test, 1 - 
                        pchisq(summary.test^2, 1)), het = c(heterog, df,
                        1 - pchisq(heterog, df)),call=match.call(),names=mf$names,omitted=!ok,tau2=tau2)
        class(rval)<-"meta.DSL"
        rval
}
   
  

"meta.MH" <-function (ntrt, nctrl, ptrt, pctrl,names=NULL,data=NULL,subset=NULL,na.action=na.fail) 
{
        if (is.null(data)) data<-sys.frame(sys.parent())
        mf<-match.call()
        mf$data<-NULL 
        mf$subset<-NULL
        mf[[1]]<-as.name("data.frame")
        mf<-eval(mf,data)
        if (!is.null(subset)) mf<-mf[subset,]
        mf<-na.action(mf)
        A <- mf$ptrt
        B <- mf$pctrl
        C <- mf$ntrt -mf$ptrt
        D <- mf$nctrl - mf$pctrl
        logORs <- log((A * D)/(B * C))
        vars <- 1/A + 1/B + 1/C + 1/D
        Ti <- A + B + C + D
        P <- (A + D)/Ti
        Q <- (B + C)/Ti
        R <- A * D/Ti
        S <- B * C/Ti
        logMH <- log(sum(R)/sum(S))
        varMH <- sum(P * R)/(2 * sum(R)^2) + sum(P * S + Q * 
                R)/(2 * sum(R) * sum(S)) + sum(Q * S)/(2 * sum(S)^2)
        MHchisq <- (sum(A - (A + B) * (A + C)/Ti)^2)/sum((A + 
                B) * (C + D) * (A + C) * (B + D)/(Ti * Ti * (Ti - 
                1)))
        ok <- is.finite(logORs)
        heterog <- sum(((logORs[ok] - logMH)^2)/vars[ok])
        df <- sum(ok) - 1
        rval<-list(logOR = logORs, selogOR = sqrt(vars), logMH = logMH, 
                selogMH = sqrt(varMH), MHtest = c(MHchisq, 1 - 
                        pchisq(MHchisq, 1)), het = c(heterog, df,
                        1 - pchisq(heterog, df)),call=match.call(),names=mf$names)
        class(rval)<-"meta.MH"
        rval
}
  
print.meta.MH<-function(obj){
  cat("Fixed effects (Mantel-Haenszel) meta-analysis\n")
  cat("Call: ")
  print(obj$call)
  ci<-exp(obj$logMH+c(-1.96,0,1.96)*obj$selogMH)
  cat(paste("Mantel-Haensel OR = ",round(ci[2],2),"  95% CI (",round(ci[1],2),",",round(ci[3],2),")\n",sep=""))
  cat(paste("Test for heterogeneity: X^2(",obj$het[2],") = ",round(obj$het[1],2)," (p-value ",round(obj$het[3],4),")\n",sep=""))
} 

summary.meta.MH<-function(obj){
  m<-exp(outer(obj$selogOR,c(0,-1.96,1.96),"*")+cbind(obj$logOR,obj$logOR,obj$logOR))
  dimnames(m)<-list(as.character(obj$names),c("OR","(lower,","upper)"))
  rval<-list(ors=m,call=obj$call,MHci=exp(obj$logMH+c(-1.96,0,1.96)*obj$selogMH),het=obj$het)
  class(rval)<-"summary.meta.MH"
  rval
}

print.summary.meta.MH<-function(obj){
  cat("Fixed effects (Mantel-Haenszel) meta-analysis\n")
  cat("Call: ")
  print(obj$call)
  cat("------------------------------------\n")
  print(round(obj$ors,2))
  cat("------------------------------------\n")
  cat(paste("Mantel-Haensel OR = ",round(obj$MHci[2],2),"  95% CI (",round(obj$MHci[1],2),",",round(obj$MHci[3],2),")\n",sep=""))
  cat(paste("Test for heterogeneity: X^2(",obj$het[2],") = ",round(obj$het[1],2)," (p-value ",round(obj$het[3],4),")\n",sep=""))
}

plot.meta.MH<-function(obj,summary=T,summlabel="Summary",...){
  if (summary)
    metaplot(obj$logOR,obj$selogOR,labels=obj$names,summn=obj$logMH,sumse=obj$selogMH,sumnn=obj$selogMH^(-2),summlabel=summlabel,logeffect=T,...)
 else 
   metaplot(obj$logOR,obj$selogOR,labels=obj$names,logeffect=T,...)
 }

print.meta.DSL<-function(obj){
  cat("Random effects (DerSimonian-Laird) meta-analysis\n")
  cat("Call: ")
  print(obj$call)
  ci<-exp(obj$logDSL+c(-1.96,0,1.96)*obj$selogDSL)
  cat(paste("Summary OR = ",round(ci[2],2),"  95% CI (",round(ci[1],2),",",round(ci[3],2),")\n",sep=""))
  cat(paste("Estimated random effects variance:",round(obj$tau2,2),"\n"))

} 

summary.meta.DSL<-function(obj){
  m<-exp(outer(obj$selogOR,c(0,-1.96,1.96),"*")+cbind(obj$logOR,obj$logOR,obj$logOR))
  dimnames(m)<-list(as.character(obj$names),c("OR","(lower,","upper)"))
  rval<-list(ors=m,call=obj$call,ci=exp(obj$logDSL+c(-1.96,0,1.96)*obj$selogDSL),het=obj$het,omitted=as.character(obj$names)[obj$omitted],tau2=obj$tau2)
  class(rval)<-"summary.meta.DSL"
  rval
}

print.summary.meta.DSL<-function(obj){
  cat("Random effects (DerSimonian-Laird) meta-analysis\n")
  cat("Call: ")
  print(obj$call)
  cat("------------------------------------\n")
  print(round(obj$ors,2))
  cat("------------------------------------\n")
  cat(paste("Summary OR = ",round(obj$ci[2],2),"  95% CI (",round(obj$ci[1],2),",",round(obj$ci[3],2),")\n",sep=""))
  cat(paste("Test for heterogeneity: X^2(",obj$het[2],") = ",round(obj$het[1],2)," (p-value ",round(obj$het[3],4),")\n",sep=""))
  cat(paste("Estimated random effects variance:",round(obj$tau2,2),"\n"))
  if (length(obj$omitted)>0){
    cat(paste("(",length(obj$omitted),"studies with zero or infinite odds ratio omitted)\n"))
  }}

plot.meta.DSL<-function(obj,summary=T,summlabel="Summary",...){
  if (summary)
    metaplot(obj$logOR,obj$selogOR,labels=obj$names,summn=obj$logDSL,sumse=obj$selogDSL,sumnn=obj$selogDSL^(-2),summlabel=summlabel,logeffect=T,...)
 else 
   metaplot(obj$logOR,obj$selogOR,labels=obj$names,logeffect=T,...)
 }

metaplot<-function(mn,se,nn=NULL,labels=NULL,xlabel="Odds ratio",ylabel="",boxwidth=NULL,xlim=NULL,summn=NULL,sumse=NULL,sumnn=NULL,summlabel="Summary",logeffect=F,lwd=3,boxsize=1,zero=if (logeffect) 1 else 0,...){
  ok<-is.finite(mn+se)
  if (is.null(xlim)) xlim<-c(min(mn[ok]-2*se[ok],na.rm=T),max(mn[ok]+2*se[ok],na.rm=T))
  ##par(pty="s")
  n<-length(mn)
  if (logeffect) {
    xlog<-"x"
    nxlim<-exp(xlim)
  }
  else {
    xlog<-""
    nxlim<-xlim
  }
  if (!is.null(labels)) {
    if (logeffect){
      nxlim[1]<-nxlim[1]/sqrt(nxlim[2]/nxlim[1])
    }
    else{
      nxlim[1]<-nxlim[1]-0.5*(nxlim[2]-nxlim[1])
    }
  }
  par(xaxt="n",yaxt="n")
  plot(nxlim,c(1,-n-2-3*!is.null(summn)),type="n",bty="n",xaxt="n",yaxt="n",log=xlog,xlab=xlabel,ylab=ylabel,...)
  par(xaxt="s")
  axis(1,at=round(10^pretty(log(exp(xlim),10),6),2))
  if (!is.null(zero))
    abline(v=zero,lty=2,lwd=2)
  lower<-mn-1.96*se
  upper<-mn+1.96*se
  if (logeffect){
    lower<-exp(lower)
    upper<-exp(upper)
  }
  for (i in 1:n){
    if (is.na(lower[i]+upper[i])) next
    lines(c(lower[i],upper[i]),c(-i,-i),lwd=lwd,...)
  }
  if (!is.null(labels))
    text(rep(nxlim[1],n),-(1:n),labels,...,adj=0)
  if (is.null(nn)) nn<-se^-2
  yscale<-0.3*boxsize/max(sqrt(nn),na.rm=T)
  if (logeffect){
    scale<-(nxlim[2]/nxlim[1])^(yscale/(4+n))
    xl<-exp(mn)*(scale^-sqrt(nn))
    xr<-exp(mn)*(scale^sqrt(nn))
   
  }
  else{
    scale<-yscale*(nxlim[2]-nxlim[1])/(4+n)
    xl<-mn-scale*sqrt(nn)
    xr<-mn+scale*sqrt(nn)
  }
  yb<-(1:n)-yscale*sqrt(nn)
  yt<-(1:n)+yscale*sqrt(nn)
  for (i in 1:n){
    if (!is.finite(mn[i])) next
   rect(xl[i],-yb[i],xr[i],-yt[i],col=par("fg"))
  }
  if (!is.null(summn)){
    if (logeffect){
      x0<-exp(summn)
      xl<-exp(summn-1.96*sumse)
      xr<-exp(summn+1.96*sumse)
      
    }
    else{
      x0<-summn
      xl<-summn-1.96*sumse
      xr<-summn+1.96*sumse
    }
    y0<-n+3
    yb<-n+3-sqrt(sumnn)*yscale
    yt<-n+3+sqrt(sumnn)*yscale
    polygon(c(xl,x0,xr,x0),-c(y0,yt,y0,yb),col=par("fg"),border=par("fg"))
    text(nxlim[1],-y0,labels=summlabel,adj=0)
  }
  
}


